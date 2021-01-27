/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include <string>

#include "tbb/concurrent_queue.h"
#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/multilevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/multilevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
class MultilevelCoarsener : public ICoarsener,
                            private MultilevelCoarsenerBase {
 private:

  using Base = MultilevelCoarsenerBase;
  using Rater = MultilevelVertexPairRater<ScorePolicy,
                                          HeavyNodePenaltyPolicy,
                                          AcceptancePolicy>;
  using Rating = typename Rater::Rating;

  struct CornerCaseCounter {
      int edgeTooBig;
      int vertexTooBig;
      int notInCommunity;
      int degreeZero;
  };

  enum class MatchingState : uint8_t {
    UNMATCHED = 0,
    MATCHING_IN_PROGRESS = 1,
    MATCHED = 2
  };

  #define STATE(X) static_cast<uint8_t>(X)
  using AtomicMatchingState = parallel::IntegralAtomicWrapper<uint8_t>;
  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  MultilevelCoarsener(Hypergraph& hypergraph,
                      const Context& context,
                      const TaskGroupID task_group_id,
                      const bool top_level) :
    Base(hypergraph, context, task_group_id, top_level),
    _rater(hypergraph, context),
    _current_vertices(),
    _matching_state(),
    _cluster_weight(),
    _matching_partner(),
    _max_allowed_node_weight(context.coarsening.max_allowed_node_weight),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true),
    _use_two_hop_matching(false),
    _use_large_edge_matching(false),
    _remove_degree_zero_hns(false),
    _use_matching_without_community_detection(false),
    _min_hash_seeds(std::vector<int>(128)) {
    _progress_bar += hypergraph.numRemovedHypernodes();

    // Initialize internal data structures parallel
    tbb::parallel_invoke([&] {
      _current_vertices.resize(hypergraph.initialNumNodes());
    }, [&] {
      _matching_state.resize(hypergraph.initialNumNodes());
    }, [&] {
      _cluster_weight.resize(hypergraph.initialNumNodes());
    }, [&] {
      _matching_partner.resize(hypergraph.initialNumNodes());
    });

    if ( _context.coarsening.use_adaptive_max_allowed_node_weight &&
          hypergraph.totalWeight() !=
          static_cast<HypernodeWeight>(hypergraph.initialNumNodes()) ) {
      // If we have a weighted instance and adaptive maximum node weight is
      // enabled we adapt the maximum allowed node such that it is greater
      // than the heaviest node of the hypergraph.
      const HypernodeWeight max_vertex_weight = tbb::parallel_reduce(
      tbb::blocked_range<HypernodeID>(ID(0), hypergraph.initialNumNodes()), 0,
      [&](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
        HypernodeWeight weight = init;
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
          if ( hypergraph.nodeIsEnabled(hn) ) {
            weight = std::max(weight, hypergraph.nodeWeight(hn));
          }
        }
        return weight;
      }, [](const HypernodeWeight lhs, const HypernodeWeight rhs) {
        return std::max(lhs, rhs);
      });
      double node_weight_multiplier = std::pow(2.0, std::ceil(std::log2(
        static_cast<double>(max_vertex_weight) / static_cast<double>(_max_allowed_node_weight))));
      if ( node_weight_multiplier > 1.0 ) {
        _max_allowed_node_weight = increaseMaximumAllowedNodeWeight(node_weight_multiplier);
      }
    }
    // By default all fallback strategies are disabled
    if ( _context.coarsening.enabled_fallback_stratgies.size() == 4 ) {
      _use_two_hop_matching = _context.coarsening.enabled_fallback_stratgies[0];
      _use_large_edge_matching = _context.coarsening.enabled_fallback_stratgies[1];
      _remove_degree_zero_hns = _context.coarsening.enabled_fallback_stratgies[2];
      _use_matching_without_community_detection = _context.coarsening.enabled_fallback_stratgies[3];
    }
    if (_use_large_edge_matching) {
      std::mt19937 rng(_context.partition.seed);
      for (int i = 0; i < 128; i++) {
        _min_hash_seeds[i] = rng();
      }
    }
  }

  MultilevelCoarsener(const MultilevelCoarsener&) = delete;
  MultilevelCoarsener(MultilevelCoarsener&&) = delete;
  MultilevelCoarsener & operator= (const MultilevelCoarsener &) = delete;
  MultilevelCoarsener & operator= (MultilevelCoarsener &&) = delete;

  ~MultilevelCoarsener() {
    parallel::parallel_free(
      _current_vertices, _matching_state,
      _cluster_weight, _matching_partner);
  };

  void disableRandomization() {
    _enable_randomization = false;
  }

 private:
  void coarsenImpl() override {
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar ) {
      _progress_bar.enable();
    }

    int pass_nr = 0;
    const HypernodeID initial_num_nodes = Base::currentNumNodes();
    while ( Base::currentNumNodes() > _context.coarsening.contraction_limit ) {
      HighResClockTimepoint round_start = std::chrono::high_resolution_clock::now();
      Hypergraph& current_hg = Base::currentHypergraph();
      DBG << V(pass_nr)
          << V(current_hg.initialNumNodes())
          << V(current_hg.initialNumEdges())
          << V(current_hg.initialNumPins());

      // Random shuffle vertices of current hypergraph
      utils::Timer::instance().start_timer("shuffle_vertices", "Shuffle Vertices");
      _current_vertices.resize(current_hg.initialNumNodes());
      vec<HypernodeID> cluster_ids(current_hg.initialNumNodes());
      tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID hn) {
        ASSERT(hn < _current_vertices.size());
        // Reset clustering
        _current_vertices[hn] = hn;
        _matching_state[hn] = STATE(MatchingState::UNMATCHED);
        _matching_partner[hn] = hn;
        cluster_ids[hn] = hn;
        if ( current_hg.nodeIsEnabled(hn) ) {
          _cluster_weight[hn] = current_hg.nodeWeight(hn);
        }
      });

      if ( _enable_randomization ) {
        utils::Randomize::instance().parallelShuffleVector(
          _current_vertices, 0UL, _current_vertices.size());
      }
      utils::Timer::instance().stop_timer("shuffle_vertices");

      // We iterate in parallel over all vertices of the hypergraph and compute its contraction partner.
      // Matched vertices are linked in a concurrent union find data structure, that also aggregates
      // weights of the resulting clusters and keep track of the number of nodes left, if we would
      // contract all matched vertices.
      utils::Timer::instance().start_timer("parallel_clustering", "Parallel Clustering");
      if ( _context.partition.show_detailed_clustering_timings ) {
        utils::Timer::instance().start_timer("clustering_level_" + std::to_string(pass_nr), "Level " + std::to_string(pass_nr));
      }
      _rater.resetMatches();
      _rater.setCurrentNumberOfNodes(current_hg.initialNumNodes());
      const HypernodeID num_hns_before_pass = current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
      const HypernodeID num_pins_before_pass = current_hg.initialNumPins();
      const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
      DBG << V(current_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
      HypernodeID current_num_nodes = num_hns_before_pass;
      tbb::enumerable_thread_specific<HypernodeID> contracted_nodes(0);
      tbb::enumerable_thread_specific<HypernodeID> num_nodes_update_threshold(0);
      tbb::enumerable_thread_specific<CornerCaseCounter> counter;
      tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID id) {
        ASSERT(id < _current_vertices.size());
        const HypernodeID hn = _current_vertices[id];
        if ( current_hg.nodeIsEnabled(hn) ) {
          // We perform rating if ...
          //  1.) The contraction limit of the current level is not reached
          //  2.) Vertex hn is not matched before
          const HypernodeID u = hn;
          if ( _matching_state[u] == STATE(MatchingState::UNMATCHED) ) {
            if ( current_num_nodes > hierarchy_contraction_limit ) {
              ASSERT(current_hg.nodeIsEnabled(hn));
              const Rating rating = _rater.rate(current_hg, hn,
                cluster_ids, _cluster_weight, _max_allowed_node_weight);
              if ( rating.target != kInvalidHypernode ) {
                const HypernodeID v = rating.target;
                HypernodeID& local_contracted_nodes = contracted_nodes.local();
                matchVertices(current_hg, u, v, cluster_ids, local_contracted_nodes);

                // To maintain the current number of nodes of the hypergraph each PE sums up
                // its number of contracted nodes locally. To compute the current number of
                // nodes, we have to sum up the number of contracted nodes of each PE. This
                // operation becomes more expensive the more PEs are participating in coarsening.
                // In order to prevent expensive updates of the current number of nodes, we
                // define a threshold which the local number of contracted nodes have to exceed
                // before the current PE updates the current number of nodes. This threshold is defined
                // by the distance to the current contraction limit divided by the number of PEs.
                // Once one PE exceeds this bound the first time it is not possible that the
                // contraction limit is reached, because otherwise an other PE would update
                // the global current number of nodes before. After update the threshold is
                // increased by the new difference (in number of nodes) to the contraction limit
                // divided by the number of PEs.
                if (  local_contracted_nodes >= num_nodes_update_threshold.local() ) {
                  current_num_nodes = num_hns_before_pass -
                    contracted_nodes.combine(std::plus<HypernodeID>());
                  num_nodes_update_threshold.local() +=
                    (current_num_nodes - hierarchy_contraction_limit) /
                    _context.shared_memory.num_threads;
                }
              } else {
                CornerCaseCounter& c = counter.local();
                //std::string out = "";
                switch (rating.state) {
                  case 11:
                    c.vertexTooBig++;
                    break;
                  case 12:
                    c.edgeTooBig++;
                    break;
                  case 13:
                    c.degreeZero++;
                    break;
                  case 15:
                    c.notInCommunity++;
                    break;
                  default:
                    break;
                }
                _matching_state[u] = rating.state;
                if ( rating.state == STATE(Rater::RatingState::VERTEX_TOO_BIG) ) {
                  // Store preferred cluster for all nodes that don't get matched, because their target is too big,
                  // so they can be used for two hop matching
                  _matching_partner[u] = rating.opt_target;
                } else if ( rating.state == STATE(Rater::RatingState::DIFFERENT_COMMUNITY) ) {
                  // Store preferred cluster for all nodes that don't get matched, because their target is in a different
                  // community, so they can be matched ignoring the community detection
                  _matching_partner[u] = rating.community_target;
                }
              }
            }
          }
        }
      });
      if ( _context.partition.show_detailed_clustering_timings ) {
        utils::Timer::instance().stop_timer("clustering_level_" + std::to_string(pass_nr));
      }
      utils::Timer::instance().stop_timer("parallel_clustering");

      current_num_nodes = num_hns_before_pass -
        contracted_nodes.combine(std::plus<HypernodeID>());
      DBG << V(current_num_nodes);

      HEAVY_COARSENING_ASSERT([&] {
        parallel::scalable_vector<HypernodeWeight> expected_weights(current_hg.initialNumNodes());
        // Verify that clustering is correct
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          const HypernodeID u = hn;
          const HypernodeID root_u = cluster_ids[u];
          if ( root_u != cluster_ids[root_u] ) {
            LOG << "Hypernode" << u << "is part of cluster" << root_u << ", but cluster"
                << root_u << "is also part of cluster" << cluster_ids[root_u];
            return false;
          }
          expected_weights[root_u] += current_hg.nodeWeight(hn);
        }

        // Verify that cluster weights are aggregated correct
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          const HypernodeID u = hn;
          const HypernodeID root_u = cluster_ids[u];
          if ( root_u == u && expected_weights[u] != _cluster_weight[u] ) {
            LOG << "The expected weight of cluster" << u << "is" << expected_weights[u]
                << ", but currently it is" << _cluster_weight[u];
            return false;
          }
        }
        return true;
      }(), "Parallel clustering computed invalid cluster ids and weights");


      parallel::scalable_vector<HypernodeID> degree_zero_hns = {};
      const double reduction_vertices_percentage =
        static_cast<double>(num_hns_before_pass) /
        static_cast<double>(current_num_nodes);
      if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
        CornerCaseCounter corner_cases = counter.combine(
          [&](CornerCaseCounter a, CornerCaseCounter b) {
            CornerCaseCounter sum;
            sum.vertexTooBig = a.vertexTooBig + b.vertexTooBig;
            sum.edgeTooBig = a.edgeTooBig + b.edgeTooBig;
            sum.notInCommunity = a.notInCommunity + b.notInCommunity;
            sum.degreeZero = a.degreeZero + b.degreeZero;
            return sum;
          });
        current_num_nodes -= executeFallbackStrategies(corner_cases, degree_zero_hns, cluster_ids, current_hg);

        const double reduction_vertices_percentage =
            static_cast<double>(num_hns_before_pass) /
            static_cast<double>(current_num_nodes);
        if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
          //TODO remove debug output
          std::cout << "Contraction limit not reached" << std::endl;
          std::cout << "num_hns_before_pass: " << num_hns_before_pass << std::endl;
          std::cout << "current_num_nodes: " << current_num_nodes << std::endl;
          std::cout << "contraction_limit: " << _context.coarsening.contraction_limit << std::endl;
          std::cout << "current_num_nodes/context.coarsening.contraction_limit: " << ((1.0*(current_num_nodes))/(_context.coarsening.contraction_limit)) << std::endl;
          break;
        }
      }
      _progress_bar += (num_hns_before_pass - current_num_nodes);

      utils::Timer::instance().start_timer("parallel_multilevel_contraction", "Parallel Multilevel Contraction");
      // Perform parallel contraction
      Base::performMultilevelContraction(std::move(cluster_ids), std::move(degree_zero_hns), round_start);
      utils::Timer::instance().stop_timer("parallel_multilevel_contraction");

      if ( _context.coarsening.use_adaptive_max_allowed_node_weight ) {
        // If the reduction ratio of the number of vertices or pins is below
        // a certain threshold, we increase the maximum allowed node weight by
        // a factor of two. Idea behind this is that if we are not able to reduce
        // the number of nodes or pins by a significant ratio, then some vertices
        // reach their maximum allowed node weight and are not able to contract
        // with other nodes, which prevents some high score contractions.
        const double reduction_pins_percentage =
          static_cast<double>(num_pins_before_pass) /
          static_cast<double>(Base::currentHypergraph().initialNumPins());
        const bool reduction_vertices_below_threshold = reduction_vertices_percentage <
          _context.coarsening.adaptive_node_weight_shrink_factor_threshold;
        const bool reduction_pins_below_threshold = reduction_pins_percentage <
          _context.coarsening.adaptive_node_weight_shrink_factor_threshold;
        if ( ( reduction_vertices_below_threshold && reduction_pins_below_threshold ) ||
             ( !reduction_vertices_below_threshold && reduction_pins_below_threshold ) ) {
          _max_allowed_node_weight = increaseMaximumAllowedNodeWeight(2.0);
        }
        DBG << V(reduction_vertices_percentage)
            << V(reduction_pins_percentage)
            << V(_max_allowed_node_weight);
      }
      ++pass_nr;
    }

    _progress_bar += (initial_num_nodes - _progress_bar.count());
    _progress_bar.disable();
    Base::finalize();
  }

  /*!
   * We maintain the invariant during clustering that each cluster has a unique
   * representative and all vertices also part of that cluster point to that
   * representative. Let v be the representative of a cluster C_v, then for
   * all nodes u \in C_v follows that cluster_ids[u] = v.
   * If we perform sequential clustering, we can simply set
   * cluster_ids[u] = cluster_ids[v] to maintain our invariant. However,
   * things become more complicated if we perform parallel clustering.
   * Especially, if two neighbors u and v are concurrently matched, we have
   * to guarantee that our clustering fullfils our invariant. There are mainly
   * two different cases, which needs special attention:
   *   1.) u is matched with v and v is matched with u concurrently
   *   2.) u is matched with v and v is matched an other vertex w concurrently
   * The following functions guarantees that our invariant is fullfilled, if
   * vertices are matched concurrently.
   */
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool matchVertices(const Hypergraph& hypergraph,
                                                        const HypernodeID u,
                                                        const HypernodeID v,
                                                        parallel::scalable_vector<HypernodeID>& cluster_ids,
                                                        HypernodeID& contracted_nodes) {
    ASSERT(u < hypergraph.initialNumNodes());
    ASSERT(v < hypergraph.initialNumNodes());
    uint8_t unmatched = STATE(MatchingState::UNMATCHED);
    uint8_t match_in_progress = STATE(MatchingState::MATCHING_IN_PROGRESS);

    // Indicates that u wants to join the cluster of v.
    // Will be important later for conflict resolution.
    bool success = false;
    const HypernodeWeight weight_u = hypergraph.nodeWeight(u);
    HypernodeWeight weight_v = _cluster_weight[v];
    if ( weight_u + weight_v <= _max_allowed_node_weight ) {

      if ( _matching_state[u].compare_exchange_strong(unmatched, match_in_progress) ) {
        _matching_partner[u] = v;
        // Current thread gets "ownership" for vertex u. Only threads with "ownership"
        // can change the cluster id of a vertex.

        uint8_t matching_state_v = _matching_state[v].load();
        // If vertex v hasn't found a matching target, but we have a vertex u that wants to match with v,
        // then we don't need the corner case handling for vertex v.
        if (matching_state_v >= 10) {
          _matching_state[v].store(unmatched);
        }
        if ( matching_state_v == STATE(MatchingState::MATCHED) ) {
          // Vertex v is already matched and will not change it cluster id any more.
          // In that case, it is safe to set the cluster id of u to the cluster id of v.
          if ( v == cluster_ids[v] ) {
            // In case v is also the representative of the cluster,
            // we change the cluster id of u to v, ...
            cluster_ids[u] = v;
            _cluster_weight[v] += weight_u;
            ++contracted_nodes;
            success = true;
          } else {
            // ... otherwise, we try again to match u with the
            // representative of the cluster.
            const HypernodeID cluster_v = cluster_ids[v];
            weight_v = _cluster_weight[cluster_v];
            if ( weight_u + weight_v <= _max_allowed_node_weight ) {
              ASSERT(_matching_state[cluster_v] == STATE(MatchingState::MATCHED));
              cluster_ids[u] = cluster_v;
              _cluster_weight[cluster_v] += weight_u;
              ++contracted_nodes;
              success = true;
            }
          }
        } else if ( _matching_state[v].compare_exchange_strong(unmatched, match_in_progress) ) {
          // Current thread has the "ownership" for u and v and can change the cluster id
          // of both vertices thread-safe.
          cluster_ids[u] = v;
          _cluster_weight[v] += weight_u;
          ++contracted_nodes;
          _matching_state[v] = STATE(MatchingState::MATCHED);
          success = true;
        } else {
          // State of v must be either MATCHING_IN_PROGRESS or an other thread changed the state
          // in the meantime to MATCHED. We have to wait until the state of v changed to
          // MATCHED or resolve the conflict if u is matched within a cyclic matching dependency

          // Conflict Resolution
          while ( _matching_state[v] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {

            // Check if current vertex is in a cyclic matching dependency
            HypernodeID cur_u = u;
            HypernodeID smallest_node_id_in_cycle = cur_u;
            while ( _matching_partner[cur_u] != u && _matching_partner[cur_u] != cur_u ) {
              cur_u = _matching_partner[cur_u];
              smallest_node_id_in_cycle = std::min(smallest_node_id_in_cycle, cur_u);
            }

            // Resolve cyclic matching dependency
            // Vertex with smallest id starts to resolve conflict
            const bool is_in_cyclic_dependency = _matching_partner[cur_u] == u;
            if ( is_in_cyclic_dependency && u == smallest_node_id_in_cycle) {
              cluster_ids[u] = v;
              _cluster_weight[v] += weight_u;
              ++contracted_nodes;
              _matching_state[v] = STATE(MatchingState::MATCHED);
              _matching_state[u] = STATE(MatchingState::MATCHED);
              success = true;
            }
          }

          // If u is still in state MATCHING_IN_PROGRESS its matching partner v
          // must be matched in the meantime with an other vertex. Therefore,
          // we try to match u with the representative v's cluster.
          if ( _matching_state[u] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {
            ASSERT( _matching_state[v] == STATE(MatchingState::MATCHED) );
            const HypernodeID cluster_v = cluster_ids[v];
            const HypernodeWeight weight_v = _cluster_weight[cluster_v];
            if ( weight_u + weight_v <= _max_allowed_node_weight ){
              cluster_ids[u] = cluster_v;
              _cluster_weight[cluster_v] += weight_u;
              ++contracted_nodes;
              success = true;
            }
          }
        }
        _rater.markAsMatched(u);
        _rater.markAsMatched(v);
        _matching_partner[u] = u;
        _matching_state[u] = STATE(MatchingState::MATCHED);
      }
    }

    return success;
  }

  PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm) override {
    return Base::doUncoarsen(label_propagation, fm);
  }

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
  }

  HypernodeID hierarchyContractionLimit(const Hypergraph& hypergraph) const {
    return std::max( static_cast<HypernodeID>( static_cast<double>(hypergraph.initialNumNodes() -
      hypergraph.numRemovedHypernodes()) / _context.coarsening.maximum_shrink_factor ),
      _context.coarsening.contraction_limit );
  }

  HypernodeWeight increaseMaximumAllowedNodeWeight(const double multiplier) {
    HypernodeWeight max_part_weight = 0;
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      max_part_weight = std::max(max_part_weight,
        _context.partition.max_part_weights[block]);
    }
    return std::min( multiplier * static_cast<double>(_max_allowed_node_weight),
      std::max( max_part_weight / _context.coarsening.max_allowed_weight_fraction,
        static_cast<double>(_context.coarsening.max_allowed_node_weight ) ) );
  }

  HypernodeID executeFallbackStrategies(
        CornerCaseCounter& corner_cases,
        parallel::scalable_vector<HypernodeID>& degree_zero_hns,
        parallel::scalable_vector<HypernodeID>& cluster_ids,
        Hypergraph& current_hg) {

    //TODO remove debug output
    HypernodeID contracted_hns = 0;

    tbb::enumerable_thread_specific < std::vector < std::pair < uint32_t, HypernodeID>>> thread_edge_targets;
    tbb::enumerable_thread_specific < std::vector < std::pair < HypernodeID, HypernodeID>>> thread_opt_targets;
    tbb::enumerable_thread_specific < std::vector < std::pair < HypernodeID, HypernodeID>>> thread_community_targets;
    tbb::enumerable_thread_specific <std::vector<HypernodeID>> thread_degree_zero_hns;
    tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID hn) {
      if (_use_two_hop_matching && _matching_state[hn] == STATE(Rater::RatingState::VERTEX_TOO_BIG)) {
        // Store preferred cluster for all nodes that don't get matched, because their target is too big,
        // so they can be used for two hop matching
        thread_opt_targets.local().push_back(std::pair<HypernodeID, HypernodeID>(_matching_partner[hn], hn));
        _matching_partner[hn] = hn;
      } else if (_use_large_edge_matching &&_matching_state[hn] == STATE(Rater::RatingState::EDGE_TOO_BIG)) {
        // Store nodes that didn't find a contraction target, because they had an edge that was too big,
        // and their minhash fingerprint
        thread_edge_targets.local().push_back(
          std::pair<uint32_t, HypernodeID>(minhash(current_hg.incidentEdges(hn)), hn));
      } else if (_remove_degree_zero_hns &&_matching_state[hn] == STATE(Rater::RatingState::NO_NEIGHBOURS)) {
        // Store nodes that have a degree of zero and thus can't be contracted anymore
        thread_degree_zero_hns.local().push_back(hn);
      } else if (_use_matching_without_community_detection &&_matching_state[hn] == STATE(Rater::RatingState::DIFFERENT_COMMUNITY)) {
        // Store preferred cluster for all nodes that don't get matched, because their target is in a different
        // community, so they can be matched ignoring the community detection
        thread_community_targets.local().push_back(std::pair<HypernodeID, HypernodeID>(_matching_partner[hn], hn));
        _matching_partner[hn] = hn;
      }
    });
    std::vector<int> sizes;
    int size_total;
    // Contract vertices which want to join a cluster that is to big with vertices that want to join the same cluster
    if (_use_two_hop_matching) {
      sizes = {0};
      size_total = 0;
      for (auto &item : thread_opt_targets) {
        size_total += item.size();
        sizes.push_back(size_total);
      }
      std::vector<std::pair<HypernodeID, HypernodeID>> opt_targets;
      opt_targets.resize(corner_cases.vertexTooBig);
      tbb::parallel_for(0, (int) thread_opt_targets.size(), [&](int i) {
        auto iter = thread_opt_targets.begin() + i;
        int offset = sizes[i];
        for (size_t j = 0; j < iter->size(); j++) {
          opt_targets[offset + j] = (*iter)[j];
        }
      });

      tbb::parallel_sort(opt_targets.begin(), opt_targets.end());
      std::vector<int> bounds = {0};
      for (int i = 1; i < (int) opt_targets.size(); i++) {
        if (opt_targets[i].first != opt_targets[i - 1].first) {
          bounds.push_back(i);
        }
      }
      bounds.push_back(opt_targets.size());
      tbb::enumerable_thread_specific<HypernodeID> two_hop_contracted_nodes(0);
      tbb::parallel_for(0, (int) bounds.size() - 1, [&](const int index) {
        int bucket_begin = bounds[index];
        int bucket_size = bounds[index + 1] - bucket_begin;

        std::sort(opt_targets.begin() + bucket_begin, opt_targets.begin() + bucket_begin + bucket_size,
                  [&](const std::pair<HypernodeID, HypernodeID> &a, const std::pair<HypernodeID, HypernodeID> &b) {
                    return current_hg.nodeWeight(a.second) < current_hg.nodeWeight(b.second);
                  });
        int l = 0;
        int r = bucket_size - 1;
        while (l < r) {
          HypernodeID u = opt_targets[bucket_begin + r].second;
          HypernodeID v = opt_targets[bucket_begin + l].second;
          const HypernodeWeight weight_u = current_hg.nodeWeight(u);
          HypernodeWeight weight_v = current_hg.nodeWeight(v);
          if (weight_u + weight_v >= _max_allowed_node_weight /*/ 2*/) {
            r--;
          } else {
            _matching_partner[v] = u;
            cluster_ids[v] = u;
            _cluster_weight[u] += weight_v;
            _matching_state[u].store(STATE(MatchingState::MATCHED));
            _matching_state[v].store(STATE(MatchingState::MATCHED));
            ++(two_hop_contracted_nodes.local());
            l++;
            r--;
          }
        }
      });
      contracted_hns += two_hop_contracted_nodes.combine(std::plus<HypernodeID>());
      //debug
      std::string two_hop =
        "Two hop contractions: " + std::to_string(two_hop_contracted_nodes.combine(std::plus<HypernodeID>())) + "\n";
      std::cout << two_hop;
    }

    // Contract vertices that don't have a contraction target, but are both pins of the same edge that is too big
    // too be considered during rating
    if (_use_large_edge_matching) {
      sizes = {0};
      size_total = 0;
      for (auto &item : thread_edge_targets) {
        size_total += item.size();
        sizes.push_back(size_total);
      }
      std::vector<std::pair<uint32_t, HypernodeID>> edge_targets;
      edge_targets.resize(corner_cases.edgeTooBig);
      tbb::parallel_for(0, (int) thread_edge_targets.size(), [&](int i) {
        auto iter = thread_edge_targets.begin() + i;
        int offset = sizes[i];
        for (size_t j = 0; j < iter->size(); j++) {
          edge_targets[offset + j] = (*iter)[j];
        }
      });

      tbb::parallel_sort(edge_targets.begin(), edge_targets.end());
      tbb::enumerable_thread_specific<HypernodeID> edge_contracted_nodes(0);
      tbb::parallel_for(0, (int) edge_targets.size() - 1, 2, [&](const int index) {
        HypernodeID u = edge_targets[index].second;
        HypernodeID v = edge_targets[index + 1].second;
        if (edge_targets[index].first == edge_targets[index + 1].first) {
          const HypernodeWeight weight_u = current_hg.nodeWeight(u);
          HypernodeWeight weight_v = current_hg.nodeWeight(v);
          if (weight_u + weight_v <= _max_allowed_node_weight) {
            _matching_partner[v] = u;
            cluster_ids[v] = u;
            _cluster_weight[u] += weight_v;
            _matching_state[u].store(STATE(MatchingState::MATCHED));
            _matching_state[v].store(STATE(MatchingState::MATCHED));
            ++(edge_contracted_nodes.local());
          }
        }
      });
      contracted_hns += edge_contracted_nodes.combine(std::plus<HypernodeID>());
      //debug
      std::string edge =
        "Edge contractions: " + std::to_string(edge_contracted_nodes.combine(std::plus<HypernodeID>())) + "\n";
      std::cout << edge;
    }

    //Prepare to remove degree zero vertices
    if (_remove_degree_zero_hns) {
      sizes = {0};
      size_total = 0;
      for (auto &item : thread_degree_zero_hns) {
        size_total += item.size();
        sizes.push_back(size_total);
      }
      degree_zero_hns.resize(corner_cases.degreeZero);
      tbb::parallel_for(0, (int) thread_degree_zero_hns.size(), [&](int i) {
        auto iter = thread_degree_zero_hns.begin() + i;
        int offset = sizes[i];
        for (size_t j = 0; j < iter->size(); j++) {
          degree_zero_hns[offset + j] = (*iter)[j];
        }
      });
      contracted_hns += degree_zero_hns.size();
      //debug
      std::string degree_zero = "Degree zero hns: " + std::to_string(degree_zero_hns.size()) + "\n";
      std::cout << degree_zero;
    }

    // Contract vertices which want to join a cluster that is in a different community
    if (_use_matching_without_community_detection) {
      sizes = {0};
      size_total = 0;
      for (auto &item : thread_community_targets) {
        size_total += item.size();
        sizes.push_back(size_total);
      }
      std::vector<std::pair<HypernodeID, HypernodeID>> community_targets;
      community_targets.resize(corner_cases.notInCommunity);
      tbb::parallel_for(0, (int) thread_community_targets.size(), [&](int i) {
        auto iter = thread_community_targets.begin() + i;
        int offset = sizes[i];
        for (size_t j = 0; j < iter->size(); j++) {
          community_targets[offset + j] = (*iter)[j];
        }
      });
      tbb::enumerable_thread_specific<HypernodeID> community_contracted_nodes(0);
      tbb::parallel_for(0, (int) community_targets.size(), [&](const int index) {
        const HypernodeID u = community_targets[index].second;
        const HypernodeID v = community_targets[index].first;
        if (_matching_state[u].load() != STATE(MatchingState::MATCHED)
            && _matching_state[v].load() != STATE(MatchingState::MATCHED)) {
          _matching_state[u].store(STATE(MatchingState::UNMATCHED));
          if (v != kInvalidHypernode) {
            HypernodeID &local_contracted_nodes = community_contracted_nodes.local();
            matchVertices(current_hg, u, v, cluster_ids, local_contracted_nodes);
          }
        }
      });
      contracted_hns += community_contracted_nodes.combine(std::plus<HypernodeID>());
      //debug
      std::string community = "Contractions with different Communities: " +
                              std::to_string(community_contracted_nodes.combine(std::plus<HypernodeID>())) + "\n";
      std::cout << community;
    }
    return contracted_hns;
  }

  template<class T>
  uint32_t minhash(T edges) {
    uint32_t min_hash = 0;
    for (const auto& seed : _min_hash_seeds) {
      uint32_t min_value = std::numeric_limits<int32_t>::max();
      for (const auto& edge : edges) {
        uint32_t hash_value = hash(edge ^ seed);
        min_value = std::min(min_value, hash_value);
      }
      min_hash = combine(min_hash, min_value);
    }
    return min_hash;
  }

  // from parlay
  inline uint32_t hash(uint32_t a) {
    uint32_t z = a + 0x9e3779b9;
    z ^= z >> 15;
    z *= 0x85ebca6b;
    z ^= z >> 13;
    z *= 0xc2b2ae3d;  // 0xc2b2ae35 for murmur3
    return z ^= z >> 16;
  }

  // from boost::hash_combine
  inline uint32_t combine(uint32_t left, uint32_t hashed_right) {
    return left ^ (hashed_right + 0x9e3779b9 + (left << 6) + (left >> 2));
  }

  using Base::_context;
  using Base::_task_group_id;
  Rater _rater;
  parallel::scalable_vector<HypernodeID> _current_vertices;
  parallel::scalable_vector<AtomicMatchingState> _matching_state;
  parallel::scalable_vector<AtomicWeight> _cluster_weight;
  parallel::scalable_vector<HypernodeID> _matching_partner;
  HypernodeWeight _max_allowed_node_weight;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
  bool _use_two_hop_matching;
  bool _use_large_edge_matching;
  bool _remove_degree_zero_hns;
  bool _use_matching_without_community_detection;
  std::vector<int> _min_hash_seeds;
};

}  // namespace mt_kahypar
