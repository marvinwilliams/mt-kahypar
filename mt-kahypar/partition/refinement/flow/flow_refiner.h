
#pragma once

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_do.h"

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flow/quotient_graph_block_scheduler.h"
#include "mt-kahypar/partition/refinement/flow/flow_region_build_policy.h"
#include "mt-kahypar/partition/refinement/flow/maximum_flow.h"

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_array.h"

#include "mt-kahypar/datastructures/flow_network.h"

namespace mt_kahypar {

template< typename TypeTraits,
          typename ExecutionPolicy,
          template<typename> class GainPolicy >
class FlowRefinerT final : public IRefiner{
    private:
        using HyperGraph = typename TypeTraits::HyperGraph;
        using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
        using TBB = typename TypeTraits::TBB;
        using HwTopology = typename TypeTraits::HwTopology;
        using GainCalculator = GainPolicy<HyperGraph>;
        using EdgeList = std::vector<std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>>;
        using Edge = std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>;

        using FlowNetwork = ds::FlowNetwork<TypeTraits>;
        using MaximumFlow = IBFS<TypeTraits, FlowNetwork>;
        using ThreadLocalFlowNetwork = tbb::enumerable_thread_specific<FlowNetwork>;
        using ThreadLocalMaximumFlow = tbb::enumerable_thread_specific<MaximumFlow>;

    public:
        explicit FlowRefinerT(HyperGraph& hypergraph, const Context& context):
            _hg(hypergraph),
            _context(context),
            _flow_network(_hg.initialNumNodes(), _hg.initialNumEdges(), _hg.initialNumNodes() + 2 * _hg.initialNumEdges()),
            _maximum_flow(_hg.initialNumNodes() + 2 * _hg.initialNumEdges(), _hg.initialNumNodes()),
            _visited(_hg.initialNumNodes() + _hg.initialNumEdges()),
            _current_num_nodes(0),
            _current_level(0),
            _execution_policy(context.refinement.flow.execution_policy_alpha),
            _num_improvements(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
            _round_delta(0) {
                initialize();
        }

        FlowRefinerT(const FlowRefinerT&) = delete;
        FlowRefinerT(FlowRefinerT&&) = delete;

        FlowRefinerT& operator= (const FlowRefinerT&) = delete;
        FlowRefinerT& operator= (FlowRefinerT&&) = delete;

    private:
        void initialize() {
            for ( const HypernodeID& hn : _hg.nodes() ) {
                unused(hn);
                ++_current_num_nodes;
            }
            _execution_policy.initialize(_hg, _current_num_nodes);
        }

        bool refineImpl(const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                        kahypar::Metrics& best_metrics) override final {
            // flow refinement is not executed on all levels of the n-level hierarchy.
            // If flow should be executed on the current level is determined by the execution policy.
            ASSERT(refinement_nodes.size() % 2 == 0);
            _current_level += refinement_nodes.size() / 2;
            _current_num_nodes += refinement_nodes.size() / 2;
            if ( !_execution_policy.execute(_current_level) ) {
                return false;
            }

            utils::Timer::instance().start_timer("flow", "Flow");

            // Initialize Quotient Graph
            // 1.) Contains edges between each adjacent block of the partition
            // 2.) Contains for each edge all hyperedges, which are cut in the
            //     corresponding bipartition
            // NOTE(heuer): If anything goes wrong in integration experiments,
            //              this should be moved inside while loop.
            utils::Timer::instance().start_timer("build_quotient_graph", "Build Quotient Graph");
            QuotientGraphBlockScheduler<TypeTraits> scheduler(_hg, _context);
            scheduler.buildQuotientGraph();
            utils::Timer::instance().stop_timer("build_quotient_graph");

            // Active Block Scheduling
            bool improvement = false;
            bool active_block_exist = true;
            size_t current_round = 1;

            utils::Timer::instance().start_timer("flow_refinement", "Flow Refinement ");
            while (active_block_exist) {
                scheduler.randomShuffleQoutientEdges();
                auto edges = scheduler.getInitialParallelEdges();
                active_block_exist = false;
                _round_delta = 0;
                //parallel here
                // TODO(reister): this looks like the right parallel primitive to parallelize the flow
                // computations. However, we should think of a threshold to abort the parallel flow computations
                // when the feeder does not contain as many edges to fully utilize the cores.
                // Furthermore, we have to think of numa-awareness. One possible solution could be to schedule
                // a pairwise flow-based refinement on the numa node which contains the most nodes of the two
                // blocks.
                tbb::parallel_do(edges,
                    [&](Edge e,
                        tbb::parallel_do_feeder<Edge>& feeder){
                        const PartitionID block_0 = e.first;
                        const PartitionID block_1 = e.second;

                        _hg.updateLocalPartInfos(block_0, block_1);

                        //LOG<< V(sched_getcpu()) << "computing:" V(block_0) << "and" << V(block_1);

                        // Heuristic: If a flow refinement never improved a bipartition,
                        //            we ignore the refinement for these block in the
                        //            second iteration of active block scheduling
                        if (_context.refinement.flow.use_improvement_history &&
                            current_round > 1 && _num_improvements[block_0][block_1] == 0) {
                            scheduler.scheduleNextBlock(feeder, block_0, block_1);
                            return;
                        }

                        const bool improved = executeAdaptiveFlow(block_0, block_1, scheduler);
                        if (improved) {

                            improvement = true;
                            active_block_exist = true;
                            scheduler.setActiveBlock(block_0, true);
                            scheduler.setActiveBlock(block_1, true);
                            _num_improvements[block_0][block_1]++;
                        }

                        scheduler.scheduleNextBlock(feeder, block_0, block_1);
                    }
                );
                //LOG << "ROUND done_______________________________________________________";


                _hg.updateGlobalPartInfos();

                // TODO(reister): I agree with you that you cannot verify the metric inside the adaptive
                // flow iterations, but what you can do is to sum up delta inside the flow iterations
                // (cut_flow_network_after - cut_flow_network_before and use an atomic to sum up the deltas
                // globally) and compare here the metric before minus the delta with
                // metrics::objective(_hg, _context.partition.objective) in an assertion.

                HyperedgeWeight current_metric = metrics::objective(_hg, _context.partition.objective);
                double current_imbalance = metrics::imbalance(_hg, _context);

                //check if the metric improved as exspected
                ASSERT(best_metrics.getMetric(_context.partition.mode, _context.partition.objective) - _round_delta
                        == current_metric,
                        "Sum of deltas is not the global improvement!"
                        << V(_context.partition.objective)
                        << V(best_metrics.getMetric(_context.partition.mode, _context.partition.objective))
                        << V(_round_delta)
                        << V(current_metric));

                //check if imbalance is feasable
                ASSERT(current_imbalance <= _context.partition.epsilon,
                        "Imbalance got to bad!"
                        << V(current_imbalance));

                //Update bestmetrics
                best_metrics.updateMetric(current_metric, _context.partition.mode, _context.partition.objective);
                best_metrics.imbalance = current_imbalance;

                current_round++;
            }
            utils::Timer::instance().stop_timer("flow_refinement");
            //LOG << "REFINEMENT done_______________________________________________________";
            utils::Timer::instance().stop_timer("flow");
            return improvement;
        }


        bool executeAdaptiveFlow(PartitionID block_0, PartitionID block_1,
         QuotientGraphBlockScheduler<TypeTraits> & quotientGraph){

            bool improvement = false;
            double alpha = _context.refinement.flow.alpha * 2.0;
            HyperedgeWeight thread_local_delta = 0;
            FlowNetwork& flow_network = _flow_network.local();
            MaximumFlow& maximum_flow = _maximum_flow.local();
            kahypar::ds::FastResetFlagArray<>& visited = _visited.local();


            do {
                alpha /= 2.0;
                flow_network.reset(block_0, block_1);
                const double old_imbalance = metrics::localBlockImbalance(_hg, _context, block_0, block_1);
                //check if imbalance is feasable
                ASSERT(old_imbalance <= _context.partition.epsilon,
                        "Old local_imbalance got to bad!"
                        << V(old_imbalance)
                        << V(metrics::imbalance(_hg, _context))
                        <<V(printPartWeights(block_0, block_1)));

                // Initialize set of cut hyperedges for blocks 'block_0' and 'block_1'
                std::vector<HyperedgeID> cut_hes;
                HyperedgeWeight cut_weight = 0;
                for (const HyperedgeID& he : quotientGraph.blockPairCutHyperedges(block_0, block_1)) {
                    cut_weight += _hg.edgeWeight(he);
                    cut_hes.push_back(he);
                }

                // Heurist 1: Don't execute 2way flow refinement for adjacent blocks
                //            in the quotient graph with a small cut
                //
                // always use heuristic
                if (cut_weight <= 10 && !isRefinementOnLastLevel()) {
                    break;
                }

                // If cut is 0 no improvement is possible
                if (cut_hes.size() == 0) {
                    break;
                }

                utils::Randomize::instance().shuffleVector(cut_hes);

                // Build Flow Problem
                CutBuildPolicy<TypeTraits>::buildFlowNetwork(_hg, _context, flow_network,
                                                cut_hes, alpha, block_0, block_1,
                                                visited);
                const HyperedgeWeight cut_flow_network_before = flow_network.build(_hg, _context, block_0, block_1);

                // Find minimum (S,T)-bipartition
                const HyperedgeWeight cut_flow_network_after = maximum_flow.minimumSTCut(_hg, flow_network, _context, block_0, block_1);

                // Maximum Flow algorithm returns infinity, if all
                // hypernodes contained in the flow problem are either
                // sources or sinks
                if (cut_flow_network_after == FlowNetwork::kInfty) {
                    break;
                }

                const HyperedgeWeight delta = cut_flow_network_before - cut_flow_network_after;
                ASSERT(cut_flow_network_before >= cut_flow_network_after,
                        "Flow calculation should not increase cut!"
                        << V(cut_flow_network_before) << V(cut_flow_network_after));

                const double current_imbalance = metrics::localBlockImbalance(_hg, _context, block_0, block_1);
                const bool equal_metric = delta == 0;
                const bool improved_metric = delta > 0;
                const bool improved_imbalance = current_imbalance < old_imbalance;
                const bool is_feasible_partition = current_imbalance <= _context.partition.epsilon;

                bool current_improvement = false;
                if ((improved_metric && (is_feasible_partition || improved_imbalance)) ||
                    (equal_metric && improved_imbalance)) {
                    //check if imbalance is feasable
                    ASSERT(current_imbalance <= _context.partition.epsilon,
                        "new local_imbalance got to bad!"
                        << V(current_imbalance));
                    improvement = true;
                    current_improvement = true;
                    thread_local_delta += delta;
                    alpha *= (alpha == _context.refinement.flow.alpha ? 2.0 : 4.0);
                }

                maximum_flow.rollback(_hg, flow_network, current_improvement);

                // Perform moves in quotient graph in order to update
                // cut hyperedges between adjacent blocks.
                if (current_improvement) {
                    for (const HypernodeID& hn : flow_network.hypernodes()) {
                        const PartitionID from = _hg.partID(hn);
                        const PartitionID to = maximum_flow.getOriginalPartition(hn);
                        if (from != to) {
                            quotientGraph.changeNodePart(hn, from, to);
                        }
                    }
                }

                // Heuristic 2: If no improvement was found, but the cut before and
                //              after is equal, we assume that the partition is close
                //              to the optimum and break the adaptive flow iterations.
                //
                // always use
                if (!improvement && cut_flow_network_before == cut_flow_network_after) {
                    break;
                }

                // TODO(reister): I added this line to update the local part weights. However, it seems
                // that the flow refiner produces imbalanced partitions. This should be a point we have to
                // address.
                //_hg.updateLocalPartInfos();
                // would update the localPartInfos with the current state of other threads. The Threads can be in the middle of
                // a flow_calculation and the Graph can be in an infeasable imbalance. I instead update only the two relevant
                // blocks in the local part infos at the beginning of the parallel flow refinement.

            } while (alpha > 1.0);

            //atomic-add improvements to _round_delta
            if(thread_local_delta > 0){
                _round_delta.fetch_and_add(thread_local_delta);
            }

            return improvement;
        }

        bool isRefinementOnLastLevel() {
          return _current_num_nodes == _hg.initialNumNodes();
        }

        //debug only
        int printPartWeights(PartitionID block_0, PartitionID block_1){
                int localWeight = 0,gloablweight = 0;
                for(int k = 0; k < _context.partition.k; k++){
                    localWeight += _hg.localPartWeight(k);
                    gloablweight += _hg.partWeight(k);
                }
                LOG << V(localWeight) << V(gloablweight);

                LOG << V(block_0) << V(_hg.localPartSize(block_0)) << V(_hg.localPartWeight(block_0));
                LOG << V(block_0) << V(_hg.partSize(block_0)) << V(_hg.partWeight(block_0));
                _hg.printBlockInfos(block_0);

                LOG << V(block_1) << V(_hg.localPartSize(block_1)) << V(_hg.localPartWeight(block_1));
                LOG << V(block_1) << V(_hg.partSize(block_1)) << V(_hg.partWeight(block_1));
                _hg.printBlockInfos(block_1);
            return 0;
        }

    HyperGraph& _hg;
    const Context& _context;
    ThreadLocalFlowNetwork _flow_network;
    ThreadLocalMaximumFlow _maximum_flow;
    ThreadLocalFastResetFlagArray _visited;
    HypernodeID _current_num_nodes;
    size_t _current_level;
    ExecutionPolicy _execution_policy;
    std::vector<std::vector<size_t> > _num_improvements;
    tbb::atomic<HyperedgeWeight> _round_delta;
};



template< typename ExecutionPolicy = Mandatory >
using FlowRefiner = FlowRefinerT<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;

} //namespace mt_kahypar