
#pragma once

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
        using Network = ds::FlowNetwork<TypeTraits>;
        using EdgeList = std::vector<std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>>;
        using Edge = std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>;

    public:
        explicit FlowRefinerT(HyperGraph& hypergraph, const Context& context):
            _hg(hypergraph),
            _context(context),
            _current_level(0),
            _execution_policy(context.refinement.flow.execution_policy_alpha),
            _num_improvements(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
            _iteration(0) {
                initialize();
        }

        FlowRefinerT(const FlowRefinerT&) = delete;
        FlowRefinerT(FlowRefinerT&&) = delete;

        FlowRefinerT& operator= (const FlowRefinerT&) = delete;
        FlowRefinerT& operator= (FlowRefinerT&&) = delete;

    private:
        void initialize() {
            _execution_policy.initialize(_hg, _hg.currentNumNodes());
        }

        bool refineImpl(const std::vector<HypernodeID>&, kahypar::Metrics& best_metrics) override final {
            // flow refinement is not executed on all levels of the n-level hierarchy.
            // If flow should be executed on the current level is determined by the execution policy.
            ++_current_level;
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
            ++_iteration;
            bool improvement = false;
            bool active_block_exist = true;
            size_t current_round = 1;

            utils::Timer::instance().start_timer("flow_refinement_" + std::to_string(_iteration),
                                                 "Flow Refinement " + std::to_string(_iteration));
            while (active_block_exist) {
                scheduler.randomShuffleQoutientEdges();
                auto edges = scheduler.getInitialParallelEdges();
                active_block_exist = false;
                //parallel here
                tbb::parallel_do(edges,
                    [&](Edge e,
                        tbb::parallel_do_feeder<Edge>& feeder){
                        const PartitionID block_0 = e.first;
                        const PartitionID block_1 = e.second;

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

                //Update bestmetrics
                _hg.updateGlobalPartInfos();
                HyperedgeWeight current_metric = metrics::objective(_hg, _context.partition.objective);
                HyperedgeWeight current_imbalance = metrics::imbalance(_hg, _context);

                best_metrics.updateMetric(current_metric, _context.partition.mode, _context.partition.objective);
                best_metrics.imbalance = current_imbalance;

                current_round++;
            }
            utils::Timer::instance().stop_timer("flow_refinement_" + std::to_string(_iteration));
            //LOG << "REFINEMENT done_______________________________________________________";

            utils::Timer::instance().stop_timer("flow");
            return improvement;
        }


        bool executeAdaptiveFlow(PartitionID block_0, PartitionID block_1,
         QuotientGraphBlockScheduler<TypeTraits> & quotientGraph){

            bool improvement = false;
            double alpha = _context.refinement.flow.alpha * 2.0;
            Network flow_network = FlowNetwork<TypeTraits>(_hg, _context, static_cast<size_t>(_hg.initialNumNodes()) + 2 * _hg.initialNumEdges());
            IBFS<TypeTraits,Network> maximum_flow = IBFS<TypeTraits,Network>(_hg, _context, flow_network);
            kahypar::ds::FastResetFlagArray<> visited = kahypar::ds::FastResetFlagArray<>(static_cast<size_t>(_hg.initialNumNodes() + _hg.initialNumEdges()));


            do {
                alpha /= 2.0;
                flow_network.reset(block_0, block_1);
                const double old_imbalance = metrics::localImbalance(_hg, _context);

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
                    return improvement;
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
                const HyperedgeWeight cut_flow_network_before = flow_network.build(block_0, block_1);

                // Find minimum (S,T)-bipartition
                const HyperedgeWeight cut_flow_network_after = maximum_flow.minimumSTCut(block_0, block_1);

                // Maximum Flow algorithm returns infinity, if all
                // hypernodes contained in the flow problem are either
                // sources or sinks
                if (cut_flow_network_after == Network::kInfty) {
                    break;
                }

                const HyperedgeWeight delta = cut_flow_network_before - cut_flow_network_after;
                ASSERT(cut_flow_network_before >= cut_flow_network_after,
                        "Flow calculation should not increase cut!"
                        << V(cut_flow_network_before) << V(cut_flow_network_after));
                /*
                ASSERT(best_metrics.getMetric(_context.partition.mode, _context.partition.objective) - delta
                        == metrics::objective(_hg, _context.partition.objective),
                        "Maximum Flow is not the minimum cut!"
                        << V(_context.partition.objective)
                        << V(best_metrics.getMetric(_context.partition.mode, _context.partition.objective))
                        << V(delta)
                        << V(metrics::objective(_hg, _context.partition.objective)));*/

                const double current_imbalance = metrics::localImbalance(_hg, _context);
                //const HyperedgeWeight old_metric = best_metrics.getMetric(_context.partition.mode, _context.partition.objective);
                //const HyperedgeWeight current_metric = old_metric - delta;


                //const bool equal_metric = current_metric == best_metrics.getMetric(_context.partition.mode, _context.partition.objective);
                const bool equal_metric = delta == 0;
                //const bool improved_metric = current_metric < best_metrics.getMetric(_context.partition.mode, _context.partition.objective);
                const bool improved_metric = delta > 0;
                //const bool improved_imbalance = current_imbalance < best_metrics.imbalance;
                const bool improved_imbalance = current_imbalance < old_imbalance;
                const bool is_feasible_partition = current_imbalance <= _context.partition.epsilon;

                bool current_improvement = false;
                if ((improved_metric && (is_feasible_partition || improved_imbalance)) ||
                    (equal_metric && improved_imbalance)) {
                    //best_metrics.updateMetric(current_metric, _context.partition.mode, _context.partition.objective);
                    //best_metrics.imbalance = current_imbalance;
                    improvement = true;
                    current_improvement = true;

                    alpha *= (alpha == _context.refinement.flow.alpha ? 2.0 : 4.0);
                }

                maximum_flow.rollback(current_improvement);

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
            } while (alpha > 1.0);

            return improvement;
        }

        bool isRefinementOnLastLevel() {
            return _hg.currentNumNodes() == _hg.initialNumNodes();
        }

        /*
        void printMetric(bool newline = false, bool endline = false) {
            if (newline) {
            DBG << "";
            }
            DBG << V(metrics::imbalance(_hg, _context))
                << V(_context.partition.objective)
                << V(metrics::objective(_hg, _context.partition.objective));
            if (endline) {
            DBG << "-------------------------------------------------------------";
            }
        }*/

    HyperGraph& _hg;
    const Context& _context;
    size_t _current_level;
    ExecutionPolicy _execution_policy;
    std::vector<std::vector<size_t> > _num_improvements;
    size_t _iteration;
};



template< typename ExecutionPolicy = Mandatory >
using FlowRefiner = FlowRefinerT<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;

} //namespace mt_kahypar