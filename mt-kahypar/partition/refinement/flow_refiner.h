
#pragma once

#include "mt-kahypar/partition/refinement/i_refiner.h"

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
    
    public:
        explicit FlowRefiner(Hypergraph& hypergraph, const Context& context):
            _hg(hypergraph),
            _context(context){

        }

        FlowRefiner(const FlowRefiner&) = delete;
        FlowRefiner(FlowRefiner&&) = delete;

        FlowRefiner& operator= (const FlowRefiner&) = delete;
        FlowRefiner& operator= (FlowRefiner&&) = delete;

    private:
        bool refineImpl(const std::vector<HypernodeID>& refinement_nodes, kahypar::Metrics& best_metrics) override final {
            // flow refinement is not executed on all levels of the n-level hierarchy.
            // If flow should be executed on the current level is determined by the execution policy.
            ++_current_level;
            if ( !_execution_policy.execute(_current_level) ) {
                return false;
            }

            return false;
        }

    HyperGraph& _hg;
    const Context& _context;
};



template< typename ExecutionPolicy = Mandatory >
using flowKm1Refiner = FlowRefiner<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;

} //namespace mt_kahypar