
#pragma once

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_do.h"
#include "tbb/task_group.h"

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flow/quotient_graph_block_scheduler.h"
#include "mt-kahypar/partition/refinement/flow/flow_region_build_policy.h"
#include "mt-kahypar/partition/refinement/flow/maximum_flow.h"

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_array.h"

#include "mt-kahypar/parallel/tbb_numa_arena.h"

#include "mt-kahypar/datastructures/flow_network.h"

namespace mt_kahypar {

template< typename TypeTraits, typename Scheduler>
class FlowRefinerT final : public IRefiner{
    private:
        using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
        using TBB = typename TypeTraits::TBB;
        using HwTopology = typename TypeTraits::HwTopology;
        using EdgeList = std::vector<std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>>;
        using Edge = std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>;

        using FlowNetwork = ds::FlowNetwork<TypeTraits>;
        using MaximumFlow = IBFS<TypeTraits, FlowNetwork>;
        using ThreadLocalFlowNetwork = tbb::enumerable_thread_specific<FlowNetwork>;
        using ThreadLocalMaximumFlow = tbb::enumerable_thread_specific<MaximumFlow>;

        struct FlowConfig {
            explicit FlowConfig(HyperGraph& hg) :
                hypergraph(hg),
                flow_network(
                   hypergraph.initialNumNodes(), hypergraph.initialNumEdges(),
                   hypergraph.initialNumNodes() + 2 * hypergraph.initialNumEdges()),
                maximum_flow(
                   hypergraph.initialNumNodes() + 2 * hypergraph.initialNumEdges(),
                   hypergraph.initialNumNodes()),
                visited(hypergraph.initialNumNodes() + hypergraph.initialNumEdges()) { }

            HyperGraph& hypergraph;
            ThreadLocalFlowNetwork flow_network;
            ThreadLocalMaximumFlow maximum_flow;
            ThreadLocalFastResetFlagArray visited;
        };

    public:
        explicit FlowRefinerT(HyperGraph&,
                              const Context& context,
                              const TaskGroupID task_group_id) :
            _context(context),
            _task_group_id(task_group_id),
            _num_improvements(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
            _round_delta(0),
            _start_new_parallel_do(TBB::instance().num_used_numa_nodes()) { }

        FlowRefinerT(const FlowRefinerT&) = delete;
        FlowRefinerT(FlowRefinerT&&) = delete;

        FlowRefinerT& operator= (const FlowRefinerT&) = delete;
        FlowRefinerT& operator= (FlowRefinerT&&) = delete;

    private:
        bool refineImpl(HyperGraph& hypergraph,
                        const parallel::scalable_vector<HypernodeID>&,
                        kahypar::Metrics& best_metrics) override final {

            utils::Timer::instance().start_timer("flow", "Flow");

            FlowConfig config(hypergraph);

            // Initialize Quotient Graph
            // 1.) Contains edges between each adjacent block of the partition
            // 2.) Contains for each edge all hyperedges, which are cut in the
            //     corresponding bipartition
            // NOTE(heuer): If anything goes wrong in integration experiments,
            //              this should be moved inside while loop.
            utils::Timer::instance().start_timer("build_quotient_graph", "Build Quotient Graph");
            Scheduler scheduler(hypergraph, _context);
            scheduler.buildQuotientGraph();
            utils::Timer::instance().stop_timer("build_quotient_graph");

            // Active Block Scheduling
            bool improvement = false;
            size_t active_blocks = _context.partition.k;
            size_t current_round = 1;

            utils::Timer::instance().start_timer("flow_refinement", "Flow Refinement ");
            while (active_blocks >= 2) {

                scheduler.randomShuffleQoutientEdges();
                auto scheduling_edges = scheduler.getInitialParallelEdges();
                scheduler.init_block_weights();
                _round_delta = 0;

                //parallel here
                TBB::instance().execute_parallel_on_all_numa_nodes(TBBNumaArena::GLOBAL_TASK_GROUP, [&](const int node) {
                    tbb::parallel_do(scheduling_edges[node],
                    [&](Edge e,
                        tbb::parallel_do_feeder<Edge>& feeder){
                            parallelFlowCalculation(
                                config, e, node, current_round,
                                improvement, scheduler, feeder);
                        });
                });

                while(scheduler.getNumberOfActiveTasks() > 0){
                    TBB::instance().wait(0);
                }

                //LOG << "ROUND done_______________________________________________________";

                HyperedgeWeight current_metric = metrics::objective(hypergraph, _context.partition.objective);
                //best_metrics.getMetric(_context.partition.mode, _context.partition.objective) - _round_delta;
                double current_imbalance = metrics::imbalance(hypergraph, _context);

                //check if the metric improved as exspected
                /*ASSERT(current_metric == metrics::objective(hypergraph, _context.partition.objective),
                        "Sum of deltas is not the global improvement!"
                        << V(_context.partition.objective)
                        << V(metrics::objective(hypergraph, _context.partition.objective))
                        << V(_round_delta)
                        << V(current_metric));*/

                //Update bestmetrics
                best_metrics.updateMetric(current_metric, _context.partition.mode, _context.partition.objective);
                best_metrics.imbalance = current_imbalance;

                //update number of active blocks
                active_blocks = scheduler.getNumberOfActiveBlocks();

                current_round++;
            }
            utils::Timer::instance().stop_timer("flow_refinement");
            //LOG << "REFINEMENT done_______________________________________________________";
            utils::Timer::instance().stop_timer("flow");
            return improvement;
        }


        void initializeImpl(HyperGraph&) override final { 
            
        }

        void parallelFlowCalculation(FlowConfig& config,
                                     const Edge& edge,
                                     const int node,
                                     const size_t& current_round,
                                     bool& improvement,
                                     Scheduler& scheduler,
                                     tbb::parallel_do_feeder<Edge>& feeder) {
            const PartitionID block_0 = edge.first;
            const PartitionID block_1 = edge.second;

            //LOG<< V(sched_getcpu()) << "computing:" V(block_0) << "and" << V(block_1);

            // Heuristic: If a flow refinement never improved a bipartition,
            //            we ignore the refinement for these block in the
            //            second iteration of active block scheduling
            if (_context.refinement.flow.use_improvement_history &&
                current_round > 1 && _num_improvements[block_0][block_1] == 0) {
                scheduleNextBlocks(
                    config, edge, node, current_round,
                    improvement,scheduler, feeder);
                //LOG << "Done with job on Numa Node" << sched_edge.first << " , Blocks:" << sched_edge.second.first << " " << sched_edge.second.second;
                return;
            }

            const bool improved = executeAdaptiveFlow(config, block_0, block_1, scheduler);
            if (improved) {
                improvement = true;
                scheduler.setActiveBlock(block_0, true);
                scheduler.setActiveBlock(block_1, true);
                _num_improvements[block_0][block_1]++;
            }

            scheduleNextBlocks(
                config, edge, node, current_round,
                improvement,scheduler, feeder);
            //LOG << "Done with job on Numa Node" << sched_edge.first << " , Blocks:" << sched_edge.second.first << " " << sched_edge.second.second;
        }

        void scheduleNextBlocks(FlowConfig& config,
                                Edge old_edge,
                                int node,
                                const size_t & current_round,
                                bool & improvement,
                                Scheduler & scheduler,
                                tbb::parallel_do_feeder<Edge>& feeder) {
            utils::Timer::instance().start_timer("schedule", "Scheduling Next Blocks ", true);
            //start new tasks on this numa node
            auto sched_edges = scheduler.scheduleNextBlocks(old_edge, node, feeder);
            //start new parallel do on stalled numa nodes
            for(const auto sched_edge:sched_edges){
                _start_new_parallel_do[sched_edge.first] = std::vector<Edge>{sched_edge.second};
                int numa_to_start = sched_edge.first;
                TBB::instance().numa_task_arena(numa_to_start).execute([&] {
                    TBB::instance().numa_task_group(TBBNumaArena::GLOBAL_TASK_GROUP, numa_to_start).run([&, numa_to_start] {
                        tbb::parallel_do(_start_new_parallel_do[numa_to_start],
                            [&](Edge e,
                                tbb::parallel_do_feeder<Edge>& newfeeder){
                                    parallelFlowCalculation(
                                        config, e, numa_to_start, current_round,
                                        improvement, scheduler, newfeeder);
                        });
                    });
                });
            }
            utils::Timer::instance().stop_timer("schedule");
        }

        bool executeAdaptiveFlow(FlowConfig& config,
                                const PartitionID block_0,
                                const PartitionID block_1,
                                Scheduler & scheduler ) {
            bool improvement = false;
            double alpha = _context.refinement.flow.alpha * 2.0;
            HyperedgeWeight thread_local_delta = 0;
            HyperGraph& hypergraph = config.hypergraph;
            FlowNetwork& flow_network = config.flow_network.local();
            MaximumFlow& maximum_flow = config.maximum_flow.local();
            kahypar::ds::FastResetFlagArray<>& visited = config.visited.local();


            do {
                alpha /= 2.0;
                flow_network.reset(block_0, block_1);
                //const double old_imbalance = metrics::localBlockImbalance(
                //hypergraph, _context, block_0, block_1);

                // Initialize set of cut hyperedges for blocks 'block_0' and 'block_1'
                std::vector<HyperedgeID> cut_hes;
                HyperedgeWeight cut_weight = 0;
                for (const HyperedgeID& he : scheduler.blockPairCutHyperedges(block_0, block_1)) {
                    cut_weight += hypergraph.edgeWeight(he);
                    cut_hes.push_back(he);
                }

                // Heurist 1: Don't execute 2way flow refinement for adjacent blocks
                //            in the quotient graph with a small cut
                //
                // always use heuristic
                if (cut_weight <= 10 ) {
                    break;
                }

                // If cut is 0 no improvement is possible
                if (cut_hes.size() == 0) {
                    break;
                }

                utils::Randomize::instance().shuffleVector(cut_hes);

                // Build Flow Problem
                CutBuildPolicy<TypeTraits>::buildFlowNetwork(
                    hypergraph, _context, flow_network,
                    cut_hes, alpha, block_0, block_1,
                    visited, scheduler);

                const HyperedgeWeight cut_flow_network_before =
                    flow_network.build(
                        hypergraph, _context, block_0, block_1);

                // Find minimum (S,T)-bipartition
                const HyperedgeWeight cut_flow_network_after =
                    maximum_flow.minimumSTCut(
                        hypergraph, flow_network, _context, block_0, block_1, scheduler);

                // Maximum Flow algorithm returns infinity, if all
                // hypernodes contained in the flow problem are either
                // sources or sinks
                if (cut_flow_network_after == FlowNetwork::kInfty) {
                    flow_network.releaseHyperNodes(hypergraph, scheduler, block_0, block_1);
                    break;
                }

                const HyperedgeWeight delta = cut_flow_network_before - cut_flow_network_after;
                ASSERT(cut_flow_network_before >= cut_flow_network_after,
                        "Flow calculation should not increase cut!"
                        << V(cut_flow_network_before) << V(cut_flow_network_after));

                //const double current_imbalance = metrics::localBlockImbalance(hypergraph, _context, block_0, block_1);
                std::vector<HypernodeWeight> aquired_part_weight = flow_network.get_aquired_part_weight(hypergraph, block_0, block_1);
                double current_imbalance = metrics::localBlockImbalanceParallel(_context, block_0, block_1, scheduler,
                    aquired_part_weight[0], aquired_part_weight[1]);

                //const bool equal_metric = delta == 0;
                const bool improved_metric = delta > 0;
                //const bool improved_imbalance = current_imbalance < old_imbalance;
                const bool is_feasible_partition = current_imbalance <=  _context.partition.epsilon;

                bool current_improvement = false;
                if (improved_metric && is_feasible_partition) {
                    improvement = true;
                    current_improvement = true;
                    thread_local_delta += delta;
                    alpha *= (alpha == _context.refinement.flow.alpha ? 2.0 : 4.0);
                }

                maximum_flow.rollback(hypergraph, flow_network, current_improvement);

                // Perform moves in quotient graph in order to update
                // cut hyperedges between adjacent blocks.
                if (current_improvement) {
                    for (const HypernodeID& ogHn : flow_network.hypernodes()) {
                        const HypernodeID& hn = hypergraph.globalNodeID(ogHn);
                        const PartitionID from = hypergraph.partID(hn);
                        const PartitionID to = maximum_flow.getOriginalPartition(ogHn);
                        if (from != to) {
                            scheduler.changeNodePart(hn, from, to);
                        }
                    }
                }

                // Heuristic 2: If no improvement was found, but the cut before and
                //              after is equal, we assume that the partition is close
                //              to the optimum and break the adaptive flow iterations.
                //
                // always use
                if (!improvement && cut_flow_network_before == cut_flow_network_after) {
                    flow_network.releaseHyperNodes(hypergraph, scheduler, block_0, block_1);
                    break;
                }

                flow_network.releaseHyperNodes(hypergraph, scheduler, block_0, block_1);
                
            } while (alpha > 1.0);

            //atomic-add improvements to _round_delta
            if(thread_local_delta > 0){
                _round_delta.fetch_and_add(thread_local_delta);
            }

            return improvement;
        }

    const Context& _context;
    const TaskGroupID _task_group_id;

    std::vector<std::vector<size_t> > _num_improvements;
    tbb::atomic<HyperedgeWeight> _round_delta;
    std::vector<std::vector<Edge>> _start_new_parallel_do;
};

using FlowRefinerMatch = FlowRefinerT<GlobalTypeTraits, MatchingScheduler<GlobalTypeTraits>>;
using FlowRefinerOpt = FlowRefinerT<GlobalTypeTraits, OptScheduler<GlobalTypeTraits>>;

} //namespace mt_kahypar