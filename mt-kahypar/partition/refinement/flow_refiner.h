
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
        bool refineImpl(const std::vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics) override final {
        }
};

} //namespace mt_kahypar