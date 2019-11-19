
#pragma once

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flow/quotient_graph_block_scheduler.h"
#include "mt-kahypar/partition/refinement/flow/flow_region_build_policy.h"
#include "mt-kahypar/partition/refinement/flow/maximum_flow.h"


#include "external_tools/kahypar/kahypar/datastructure/flow_network.h"
#include "external_tools/kahypar/kahypar/datastructure/fast_reset_array.h"

#include "mt-kahypar/datastructures/flow_network.h"

namespace mt_kahypar {

template< typename TypeTraits,
          typename ExecutionPolicy,
          template<typename> class GainPolicy >
class FlowRefiner final : public IRefiner{
    private:
        using HyperGraph = typename TypeTraits::HyperGraph;
        using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
        using TBB = typename TypeTraits::TBB;
        using HwTopology = typename TypeTraits::HwTopology;
        using GainCalculator = GainPolicy<HyperGraph>;
        using Network = ds::FlowNetwork<TypeTraits>;
    
    public:
        explicit FlowRefiner(Hypergraph& hypergraph, const Context& context):
            _hg(hypergraph),
            _context(context),
            _flow_network(_hg, _context, static_cast<size_t>(hypergraph.initialNumNodes()) + 2 * hypergraph.initialNumEdges()),
            _current_level(0),
            _execution_policy(context.refinement.execution_policy_alpha),
            _num_improvements(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
            _maximum_flow(IBFS(hypergraph, context, _flow_network)),
            _visited(static_cast<size_t>(hypergraph.initialNumNodes() + hypergraph.initialNumEdges())) {
                initialize();
        }

        FlowRefiner(const FlowRefiner&) = delete;
        FlowRefiner(FlowRefiner&&) = delete;

        FlowRefiner& operator= (const FlowRefiner&) = delete;
        FlowRefiner& operator= (FlowRefiner&&) = delete;

    private:
        void initialize() {
            HypernodeID current_num_nodes = 0;
            for ( const HypernodeID& hn : _hg.nodes() ) {
            ++current_num_nodes;
            }
            _execution_policy.initialize(_hg, current_num_nodes);
        }

        bool refineImpl(const std::vector<HypernodeID>& refinement_nodes, kahypar::Metrics& best_metrics) override final {
            // flow refinement is not executed on all levels of the n-level hierarchy.
            // If flow should be executed on the current level is determined by the execution policy.
            ++_current_level;
            if ( !_execution_policy.execute(_current_level) ) {
                return false;
            }
            // Initialize Quotient Graph
            // 1.) Contains edges between each adjacent block of the partition
            // 2.) Contains for each edge all hyperedges, which are cut in the
            //     corresponding bipartition
            // NOTE(heuer): If anything goes wrong in integration experiments,
            //              this should be moved inside while loop.
            QuotientGraphBlockScheduler scheduler(_hg, _context);
            scheduler.buildQuotientGraph();

            // Active Block Scheduling
            bool improvement = false;
            bool active_block_exist = true;
            std::vector<bool> active_blocks(_context.partition.k, true);
            size_t current_round = 1;
            while (active_block_exist) {
                scheduler.randomShuffleQoutientEdges();
                std::vector<bool> tmp_active_blocks(_context.partition.k, false);
                active_block_exist = false;
                for (const auto& e : scheduler.qoutientGraphEdges()) {
                    const PartitionID block_0 = e.first;
                    const PartitionID block_1 = e.second;

                    // Heuristic: If a flow refinement never improved a bipartition,
                    //            we ignore the refinement for these block in the
                    //            second iteration of active block scheduling
                    if (_context.refinement.flow.use_improvement_history &&
                        current_round > 1 && _num_improvements[block_0][block_1] == 0) {
                    continue;
                    }

                    if (active_blocks[block_0] || active_blocks[block_1]) {
                        /*_twoway_flow_refiner.updateConfiguration(block_0, block_1, &scheduler, true);
                        const bool improved = _twoway_flow_refiner.refine(refinement_nodes,
                                                                            max_allowed_part_weights,
                                                                            changes,
                                                                            best_metrics);*/
                        const bool improved = executeAdaptiveFlow(block_0, block_1, scheduler, best_metrics);
                        if (improved) {
                            /*DBG << "Improvement found beetween blocks " << block_0 << " and "
                                << block_1 << " in round #"
                                << current_round;
                            printMetric();*/
                            improvement = true;
                            active_block_exist = true;
                            tmp_active_blocks[block_0] = true;
                            tmp_active_blocks[block_1] = true;
                            _num_improvements[block_0][block_1]++;
                        }
                    }
                }
                current_round++;
                std::swap(active_blocks, tmp_active_blocks);
            }

            //printMetric(true, true);

            return improvement;
        }

        
        bool executeAdaptiveFlow(PartitionID block_0, PartitionID block_1,
         QuotientGraphBlockScheduler & quotientGraph, kahypar::Metrics& best_metrics){
            
            bool improvement = false;
            double alpha = _context.refinement.flow.alpha * 2.0;
            
            do {
            alpha /= 2.0;
            _flow_network.reset(block_0, block_1);

            //DBG << "";
            //DBG << V(alpha);
            
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
                //DBG << "Cut is zero";
                break;
            }
            
            
            utils::Randomize::instance().shuffleVector(cut_hes, cut_hes.size(), 1);
            //std::shuffle(cut_hes.begin(), cut_hes.end(),Randomize::instance().getGenerator());

            // Build Flow Problem
            CutBuildPolicy::buildFlowNetwork(_hg, _context, _flow_network,
                                            cut_hes, alpha, block_0, block_1,
                                            _visited);
            const HyperedgeWeight cut_flow_network_before = _flow_network.build(block_0, block_1);
            
            //DBG << V(_flow_network.numNodes()) << V(_flow_network.numEdges());

            //printMetric();

            // Find minimum (S,T)-bipartition
            const HyperedgeWeight cut_flow_network_after = _maximum_flow->minimumSTCut(block_0, block_1);

            // Maximum Flow algorithm returns infinity, if all
            // hypernodes contained in the flow problem are either
            // sources or sinks
            if (cut_flow_network_after == Network::kInfty) {
                //DBG << "Trivial Cut";
                break;
            }

            const HyperedgeWeight delta = cut_flow_network_before - cut_flow_network_after;
            ASSERT(cut_flow_network_before >= cut_flow_network_after,
                    "Flow calculation should not increase cut!"
                    << V(cut_flow_network_before) << V(cut_flow_network_after));
            ASSERT(best_metrics.getMetric(_context.partition.mode, _context.partition.objective) - delta
                    == metrics::objective(_hg, _context.partition.objective),
                    "Maximum Flow is not the minimum cut!"
                    << V(_context.partition.objective)
                    << V(best_metrics.getMetric(_context.partition.mode, _context.partition.objective))
                    << V(delta)
                    << V(metrics::objective(_hg, _context.partition.objective)));

            const double current_imbalance = metrics::imbalance(_hg, _context);
            const HyperedgeWeight old_metric = best_metrics.getMetric(_context.partition.mode, _context.partition.objective);
            const HyperedgeWeight current_metric = old_metric - delta;

            /*DBG << V(cut_flow_network_before)
                << V(cut_flow_network_after)
                << V(delta)
                << V(old_metric)
                << V(current_metric);

            printMetric();*/

            const bool equal_metric = current_metric == best_metrics.getMetric(_context.partition.mode, _context.partition.objective);
            const bool improved_metric = current_metric < best_metrics.getMetric(_context.partition.mode, _context.partition.objective);
            const bool improved_imbalance = current_imbalance < best_metrics.imbalance;
            const bool is_feasible_partition = current_imbalance <= _context.partition.epsilon;

            bool current_improvement = false;
            if ((improved_metric && (is_feasible_partition || improved_imbalance)) ||
                (equal_metric && improved_imbalance)) {
                best_metrics.updateMetric(current_metric, _context.partition.mode, _context.partition.objective);
                best_metrics.imbalance = current_imbalance;
                improvement = true;
                current_improvement = true;
                
                alpha *= (alpha == _context.refinement.flow.alpha ? 2.0 : 4.0);
            }

            _maximum_flow->rollback(current_improvement);

            // Perform moves in quotient graph in order to update
            // cut hyperedges between adjacent blocks.
            if (current_improvement) {
                for (const HypernodeID& hn : _flow_network.hypernodes()) {
                const PartitionID from = _hg.partID(hn);
                const PartitionID to = _maximum_flow->getOriginalPartition(hn);
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

            //printMetric(true, true);

            // Delete quotient graph
            /*if (delete_quotientgraph_after_flow) {
            delete _quotient_graph;
            _quotient_graph = nullptr;
            }*/

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
    Network _flow_network;
    size_t _current_level;
    ExecutionPolicy _execution_policy;
    std::vector<std::vector<size_t> > _num_improvements;
    std::unique_ptr<MaximumFlow<Network> > _maximum_flow;
    kahypar::ds::FastResetFlagArray<> _visited;
};



template< typename ExecutionPolicy = Mandatory >
using flowKm1Refiner = FlowRefiner<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;

} //namespace mt_kahypar