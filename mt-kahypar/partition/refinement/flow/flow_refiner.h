
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

template<typename FlowTypeTraits>
class FlowRefiner final : public IRefiner<>{
    private:
        using HyperGraph = PartitionedHypergraph<>;
        using Scheduler = typename FlowTypeTraits::Scheduler;
        using RegionBuildPolicy =  typename FlowTypeTraits::RegionBuildPolicy;
        using FlowNetwork = typename FlowTypeTraits::FlowNetwork;

        using EdgeList = parallel::scalable_vector<std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>>;
        using Edge = std::pair<mt_kahypar::PartitionID, mt_kahypar::PartitionID>;

        using MaximumFlow = IBFS<FlowTypeTraits>;
        using GainPolicy = Km1Policy<HyperGraph>;
        using ThreadLocalFlowNetwork = tbb::enumerable_thread_specific<FlowNetwork>;
        using ThreadLocalMaximumFlow = tbb::enumerable_thread_specific<MaximumFlow>;
        using ThreadLocalGain = tbb::enumerable_thread_specific<GainPolicy>;

        struct FlowConfig {
            explicit FlowConfig(PartitionedHypergraph<>& hg, const Context& context) :
                hypergraph(hg),
                flow_network(
                   hypergraph.initialNumNodes(), hypergraph.initialNumEdges(),
                   hypergraph.initialNumNodes() + 2 * hypergraph.initialNumEdges()),
                maximum_flow(
                   hypergraph.initialNumNodes() + 2 * hypergraph.initialNumEdges(),
                   hypergraph.initialNumNodes()),
                visited(hypergraph.initialNumNodes() + hypergraph.initialNumEdges()),
                gain(context) { }

            PartitionedHypergraph<>& hypergraph;
            ThreadLocalFlowNetwork flow_network;
            ThreadLocalMaximumFlow maximum_flow;
            ThreadLocalFastResetFlagArray visited;
            ThreadLocalGain gain;
        };

    public:
        explicit FlowRefiner(PartitionedHypergraph<>&,
                              const Context& context,
                              const TaskGroupID task_group_id) :
            _context(context),
            _task_group_id(task_group_id),
            _num_improvements(context.partition.k, parallel::scalable_vector<size_t>(context.partition.k, 0)),
            _round_delta(0),
            _start_new_parallel_do(TBBNumaArena::instance().num_used_numa_nodes()),
            _times_per_block(context.partition.k, parallel::scalable_vector<double>(context.partition.k, 0)),
            _improved_per_block(context.partition.k, parallel::scalable_vector<size_t>(context.partition.k, 0)),
            _rounds_per_block(context.partition.k, parallel::scalable_vector<size_t>(context.partition.k, 0)),
            _iterations_per_block(context.partition.k, parallel::scalable_vector<size_t>(context.partition.k, 0)) {
                utils::Stats::instance().add_stat("blocks_refined", 0);
                utils::Stats::instance().add_stat("time_per_block", double(0));
                utils::Stats::instance().add_stat("iterations_per_block", 0);
            }

        FlowRefiner(const FlowRefiner&) = delete;
        FlowRefiner(FlowRefiner&&) = delete;

        FlowRefiner& operator= (const FlowRefiner&) = delete;
        FlowRefiner& operator= (FlowRefiner&&) = delete;

    private:
        bool refineImpl(PartitionedHypergraph<>& hypergraph,
                        kahypar::Metrics& best_metrics) override final {

            FlowConfig config(hypergraph, _context);
            resetStats();
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
            scheduler.init_block_weights();

            utils::Timer::instance().start_timer("flow_refinement", "Flow Refinement ");
            _round_delta = 0;
            do{       
                scheduler.randomShuffleQoutientEdges();
                auto scheduling_edges = scheduler.getInitialParallelEdges();                                

                //parallel here
                tbb::parallel_do(scheduling_edges,
                [&](Edge e,
                    tbb::parallel_do_feeder<Edge>& feeder){
                        parallelFlowCalculation(
                            config, e,
                            improvement, scheduler, feeder);
                    });
                //LOG << "ROUND done_______________________________________________________";
            }while (scheduler.hasNextRound());

            HyperedgeWeight current_metric = best_metrics.getMetric(_context.partition.mode, _context.partition.objective) - _round_delta;        

            double current_imbalance = metrics::imbalance(hypergraph, _context);

            //check if the metric improved as exspected
            ASSERT(current_metric == metrics::objective(hypergraph, _context.partition.objective),
                    "Sum of deltas is not the global improvement!"
                    << V(_context.partition.objective)
                    << V(metrics::objective(hypergraph, _context.partition.objective))
                    << V(_round_delta)
                    << V(current_metric));

            //Update bestmetrics
            best_metrics.updateMetric(current_metric, _context.partition.mode, _context.partition.objective);
            best_metrics.imbalance = current_imbalance;

            utils::Timer::instance().stop_timer("flow_refinement");
            //printStats();
            //LOG << "REFINEMENT done_______________________________________________________";
            updateStats();
            return improvement;
        }

        void initializeImpl(PartitionedHypergraph<>&) override final {

        }

        void parallelFlowCalculation(FlowConfig& config,
                                     const Edge& edge,
                                     bool& improvement,
                                     Scheduler& scheduler,
                                     tbb::parallel_do_feeder<Edge>& feeder) {
            const PartitionID block_0 = edge.first;
            const PartitionID block_1 = edge.second;
            _rounds_per_block[block_0][block_1] ++;

            auto start = std::chrono::high_resolution_clock::now();
            //LOG<< V(sched_getcpu()) << "computing:" V(block_0) << "and" << V(block_1);

            const bool improved = executeAdaptiveFlow(config, block_0, block_1, scheduler);
            if (improved) {
                improvement = true;
                scheduler.setImprovement(block_0, block_1);
                scheduler.setBlocksActive(block_0, block_1, feeder);   
            }
            
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            _times_per_block[block_0][block_1] += elapsed.count();
            scheduler.scheduleNextBlocks(edge, feeder);
            //LOG << "Done with job on Numa Node" << sched_edge.first << " , Blocks:" << sched_edge.second.first << " " << sched_edge.second.second;
        }

        bool executeAdaptiveFlow(FlowConfig& config,
                                const PartitionID block_0,
                                const PartitionID block_1,
                                Scheduler & scheduler ) {
            bool improvement = false;
            double alpha = _context.refinement.flow.alpha * 2.0;
            HyperedgeWeight thread_local_delta = 0;
            size_t times_no_real_improvement = 0;
            PartitionedHypergraph<>& hypergraph = config.hypergraph;
            FlowNetwork& flow_network = config.flow_network.local();
            MaximumFlow& maximum_flow = config.maximum_flow.local();
            kahypar::ds::FastResetFlagArray<>& visited = config.visited.local();
            GainPolicy gain = config.gain.local();


            do {
                utils::Timer::instance().start_timer("resetNetwork", "Reset FlowNetwork ", true);
                _iterations_per_block[block_0][block_1] ++;
                alpha /= 2.0;
                flow_network.reset(block_0, block_1);
                utils::Timer::instance().stop_timer("resetNetwork");

                // Initialize set of cut hyperedges for blocks 'block_0' and 'block_1'
                utils::Timer::instance().start_timer("getCutHe", "Get Cut He's ", true);
                parallel::scalable_vector<HyperedgeID> cut_hes;
                HyperedgeWeight cut_weight = 0;
                for (const HyperedgeID& he : scheduler.blockPairCutHyperedges(block_0, block_1)) {
                    cut_weight += hypergraph.edgeWeight(he);
                    cut_hes.push_back(he);
                }
                utils::Timer::instance().stop_timer("getCutHe");


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
                utils::Timer::instance().start_timer("aquHn", "Building Region ", true);
                RegionBuildPolicy::buildFlowNetwork(
                    hypergraph, _context, flow_network,
                    cut_hes, alpha, block_0, block_1,
                    visited, scheduler);
                utils::Timer::instance().stop_timer("aquHn");

                utils::Timer::instance().start_timer("buildFlow", "Building Flow Network ", true);
                const HyperedgeWeight cut_flow_network_before =
                    flow_network.build(
                        hypergraph, _context, block_0, block_1, scheduler);
                utils::Timer::instance().stop_timer("buildFlow");

                // Find minimum (S,T)-bipartition
                const HyperedgeWeight cut_flow_network_after =
                    maximum_flow.minimumSTCut(
                        hypergraph, flow_network, _context, block_0, block_1, scheduler, cut_flow_network_before);

                utils::Timer::instance().start_timer("apply", "Applying Improvement ", true);
                // Maximum Flow algorithm returns infinity, if all
                // hypernodes contained in the flow problem are either
                // sources or sinks
                if (cut_flow_network_after == kInfty) {
                    flow_network.release(hypergraph, block_0, block_1, scheduler);
                    break;
                }

                const HyperedgeWeight flownetwork_delta = cut_flow_network_before - cut_flow_network_after;
                ASSERT(cut_flow_network_before >= cut_flow_network_after,
                        "Flow calculation should not increase cut!"
                        << V(cut_flow_network_before) << V(cut_flow_network_after));

                double new_imbalance = maximum_flow.get_new_imbalance();

                //const bool equal_metric = flownetwork_delta == 0;
                const bool flownetwork_improved_metric = flownetwork_delta > 0;
                //const bool improved_imbalance = current_imbalance < old_imbalance;
                const bool is_feasible_partition = new_imbalance <=  _context.partition.epsilon;

                // This function is passed as lambda to the changeNodePart function and used
                // to calculate the "real" delta of a move (in terms of the used objective function).
                auto objective_delta = [&](const HyperedgeID he,
                                        const HyperedgeWeight edge_weight,
                                        const HypernodeID edge_size,
                                        const HypernodeID pin_count_in_from_part_after,
                                        const HypernodeID pin_count_in_to_part_after) {
                                        gain.computeDeltaForHyperedge(he, edge_weight, edge_size,
                                                                        pin_count_in_from_part_after, pin_count_in_to_part_after);
                                    };

                HyperedgeWeight real_delta = 0;
                parallel::scalable_vector<std::pair<HypernodeID, std::pair<PartitionID, PartitionID>>> moves;
                // Perform moves in quotient graph in order to update
                // cut hyperedges between adjacent blocks.
                if (flownetwork_improved_metric && is_feasible_partition) {
                    auto& assignment = maximum_flow.get_assignment();
                    for (const HypernodeID& hn : flow_network.hypernodes()) {
                        const PartitionID from = hypergraph.partID(hn);
                        const PartitionID to = assignment[hn]? block_1: block_0;
                        if (from != to) {
                            HyperedgeWeight delta_before = gain.localDelta();
                            scheduler.changeNodePart(hn, from, to, objective_delta);
                            real_delta += delta_before - gain.localDelta();
                            moves.push_back(std::make_pair(hn, std::make_pair(from, to)));
                        }
                    }
                    
                    // Heuristic: Abort Round if there are to many flownetwork improvements without a
                    //            real improvement to prevent a busy deadlock.
                    if(real_delta <= 0){
                        if(++times_no_real_improvement > 3){
                            flow_network.release(hypergraph, block_0, block_1, scheduler);
                            utils::Timer::instance().stop_timer("apply");
                            break;
                        }
                    }

                    if(!_context.refinement.flow.only_real){
                        improvement = true;
                        alpha *= (alpha == _context.refinement.flow.alpha ? 2.0 : 4.0);
                    }
                }
                
                if (real_delta > 0) {
                    // update local delta
                    thread_local_delta += real_delta;
                    _improved_per_block[block_0][block_1] += real_delta;
                    if(_context.refinement.flow.only_real){
                        improvement = true;
                        alpha *= (alpha == _context.refinement.flow.alpha ? 2.0 : 4.0);
                    }
                } else if(real_delta < 0){
                    //reverse changes
                    for(auto move:moves){
                        scheduler.changeNodePart(move.first, move.second.second, move.second.first);
                    }
                }
                utils::Timer::instance().stop_timer("apply");

                // Heuristic 2: If no improvement was found, but the cut before and
                //              after is equal, we assume that the partition is close
                //              to the optimum and break the adaptive flow iterations.
                //
                // always use
                if (!improvement && real_delta <= 0) {
                    flow_network.release(hypergraph, block_0, block_1, scheduler);
                    break;
                }

                flow_network.release(hypergraph, block_0, block_1, scheduler);

            } while (alpha > 1.0);

            //atomic-add improvements to _round_delta
            if(thread_local_delta > 0){
                _round_delta.fetch_and_add(thread_local_delta);
            }

            return improvement;
        }
        void printStats(){
            for (int i = 0; i < _context.partition.k; i++){
                for (int j = 0; j < _context.partition.k; j++){
                    if(i < j){
                        LOG << "[" << i << "," << j << "]:" << _times_per_block[i][j] << " improved:" << _improved_per_block[i][j]
                        << " rounds:" << _rounds_per_block[i][j];
                    }   
                }
            }
        }

        void updateStats(){
            int blocks_refined = 0;
            double time_per_block = 0;
            int iterations_per_block = 0;
            for (int i = 0; i < _context.partition.k; i++){
                for (int j = 0; j < _context.partition.k; j++){
                    if(i < j){
                        time_per_block += _times_per_block[i][j];
                        blocks_refined += _rounds_per_block[i][j];
                        iterations_per_block += _iterations_per_block[i][j];

                    }   
                }
            }
            time_per_block = time_per_block / (double) blocks_refined;
            utils::Stats::instance().update_stat("blocks_refined", blocks_refined);
            utils::Stats::instance().update_stat("time_per_block", time_per_block);
            utils::Stats::instance().update_stat("iterations_per_block", iterations_per_block);
        }

        void resetStats(){
            for (int i = 0; i < _context.partition.k; i++){
                for (int j = 0; j < _context.partition.k; j++){
                    if(i < j){
                        _times_per_block[i][j] = 0;
                        _rounds_per_block[i][j] = 0;
                        _improved_per_block[i][j] = 0;
                        _iterations_per_block[i][j] = 0;
                    }   
                }
            }
        }

    const Context& _context;
    const TaskGroupID _task_group_id;

    parallel::scalable_vector<parallel::scalable_vector<size_t> > _num_improvements;
    tbb::atomic<HyperedgeWeight> _round_delta;
    parallel::scalable_vector<parallel::scalable_vector<Edge>> _start_new_parallel_do;
    parallel::scalable_vector<parallel::scalable_vector<double>> _times_per_block;
    parallel::scalable_vector<parallel::scalable_vector<size_t>> _improved_per_block;
    parallel::scalable_vector<parallel::scalable_vector<size_t>> _rounds_per_block;
    parallel::scalable_vector<parallel::scalable_vector<size_t>> _iterations_per_block;
};

struct FlowMatchingTypeTraits{
    using Scheduler = MatchingScheduler;
    using RegionBuildPolicy = MatchingFlowRegionBuildPolicy;
    using FlowNetwork = ds::MatchingFlowNetwork<FlowMatchingTypeTraits>;
    using MostBalancedMinimumCut = MatchingMostBalancedMinimumCut<FlowMatchingTypeTraits>;
};

struct FlowOptTypeTraits{
    using Scheduler = OptScheduler;
    using RegionBuildPolicy = OptFlowRegionBuildPolicy;
    using FlowNetwork = ds::OptFlowNetwork<FlowOptTypeTraits>;
    using MostBalancedMinimumCut = OptMostBalancedMinimumCut<FlowOptTypeTraits>;
};

struct FlowOneRoundTypeTraits{
    using Scheduler = OneRoundScheduler;
    using RegionBuildPolicy = OptFlowRegionBuildPolicy;
    using FlowNetwork = ds::OptFlowNetwork<FlowOneRoundTypeTraits>;
    using MostBalancedMinimumCut = OptMostBalancedMinimumCut<FlowOneRoundTypeTraits>;
};

using FlowRefinerMatch = FlowRefiner<FlowMatchingTypeTraits>;
using FlowRefinerOpt = FlowRefiner<FlowOptTypeTraits>;
using FlowRefinerOneRound = FlowRefiner<FlowOneRoundTypeTraits>;

} //namespace mt_kahypar