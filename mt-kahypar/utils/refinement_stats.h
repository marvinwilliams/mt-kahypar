/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <mutex>
#include <string>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace utils {

namespace {
  double improvement_share(const HyperedgeWeight improvement, const HyperedgeWeight total) {
    return total > 0 ? static_cast<double>(improvement) / total : 0;
  }
}

struct SearchStats {
  explicit SearchStats(const std::string& algo) :
    algorithm(algo),
    moved_nodes(0),
    num_searches(0),
    expected_improvements(0),
    conflicts(0),
    balance_violations(0),
    zero_gain_improvements(0),
    positive_gain_improvements(0),
    correct_gains(0) { }

  void operator+=(const SearchStats& other) {
    moved_nodes += other.moved_nodes;
    num_searches += other.num_searches;
    expected_improvements += other.expected_improvements;
    conflicts += other.conflicts;
    balance_violations += other.balance_violations;
    zero_gain_improvements += other.zero_gain_improvements;
    positive_gain_improvements += other.positive_gain_improvements;
    correct_gains += other.correct_gains;
  }

  std::string algorithm;
  size_t moved_nodes;
  size_t num_searches;
  size_t expected_improvements;
  size_t conflicts;
  size_t balance_violations;
  size_t zero_gain_improvements;
  size_t positive_gain_improvements;
  size_t correct_gains;
};

inline std::ostream & operator<< (std::ostream& str, const SearchStats& stats) {
  str << stats.algorithm << "_num_searches=" << stats.num_searches
      << " " << stats.algorithm << "_expected_improvements=" << stats.expected_improvements
      << " " << stats.algorithm << "_moved_nodes=" << stats.moved_nodes
      << " " << stats.algorithm << "_conflicts=" << stats.conflicts
      << " " << stats.algorithm << "_balance_violations=" << stats.balance_violations
      << " " << stats.algorithm << "_zero_gain_improvements=" << stats.zero_gain_improvements
      << " " << stats.algorithm << "_positive_gain_improvements=" << stats.positive_gain_improvements
      << " " << stats.algorithm << "_correct_gain=" << stats.correct_gains;
  return str;
}

struct LabelPropagationRoundStats {
  explicit LabelPropagationRoundStats(const int i) :
    cur_round(i),
    active_nodes(0),
    moved_nodes(0),
    improvement(0) { }

  void operator+=(const LabelPropagationRoundStats& other) {
    active_nodes += other.active_nodes;
    moved_nodes += other.moved_nodes;
    improvement += other.improvement;
  }

  int cur_round;
  size_t active_nodes;
  size_t moved_nodes;
  HyperedgeWeight improvement;
};

inline std::ostream & operator<< (std::ostream& str, const LabelPropagationRoundStats& lp) {
  str << "lp_round_" << lp.cur_round << "_active_nodes=" << lp.active_nodes
      << " lp_round_" << lp.cur_round << "_moved_nodes=" << lp.moved_nodes
      << " lp_round_" << lp.cur_round << "_improvement=" << lp.improvement;
  return str;
}

struct LabelPropagationStats {
  explicit LabelPropagationStats() :
    stats("lp"),
    improvement(0),
    lp_rounds() { }

  void operator+=(const LabelPropagationStats& other) {
    stats += other.stats;
    improvement += other.improvement;
  }

  LabelPropagationRoundStats& currentRound() {
    return lp_rounds.back();
  }

  SearchStats stats;
  HyperedgeWeight improvement;
  vec<LabelPropagationRoundStats> lp_rounds;
};

inline std::ostream & operator<< (std::ostream& str, const LabelPropagationStats& lp) {
  str << lp.stats
      << " lp_improvement=" << lp.improvement;
  return str;
}

struct FMStats {
  explicit FMStats() :
    stats("fm"),
    global_rollback_triggered(false),
    improvement(0) { }

  void operator+=(const FMStats& other) {
    stats += other.stats;
    improvement += other.improvement;
  }

  SearchStats stats;
  bool global_rollback_triggered;
  HyperedgeWeight improvement;
};

inline std::ostream & operator<< (std::ostream& str, const FMStats& fm) {
  str << fm.stats
      << " fm_global_rollback_triggered=" << fm.global_rollback_triggered
      << " fm_improvement=" << fm.improvement;
  return str;
}

struct FlowStats {
  explicit FlowStats() :
    stats("flows"),
    num_time_limits(0),
    improvement(0) { }

  void operator+=(const FlowStats& other) {
    stats += other.stats;
    num_time_limits += other.num_time_limits;
    improvement += other.improvement;
  }

  SearchStats stats;
  size_t num_time_limits;
  HyperedgeWeight improvement;
};

inline std::ostream & operator<< (std::ostream& str, const FlowStats& flows) {
  str << flows.stats
      << " flows_time_limit=" << flows.num_time_limits
      << " flows_improvement=" << flows.improvement;
  return str;
}

struct Search {
  explicit Search(const HypernodeID nodes,
                  const HyperedgeID edges,
                  const HypernodeID pins) :
    num_nodes(nodes),
    num_edges(edges),
    num_pins(pins),
    max_lp_rounds(0),
    num_lp_rounds(0),
    lp_stats(std::thread::hardware_concurrency()),
    fm_stats(std::thread::hardware_concurrency()),
    flow_stats(std::thread::hardware_concurrency()) { }

  void setMaxLPRounds(const size_t max_rounds) {
    max_lp_rounds = max_rounds;
  }

  void newLPRound(const bool is_main) {
    if ( is_main ) {
      for ( LabelPropagationStats& lp : lp_stats ) {
        lp.lp_rounds.emplace_back(lp.lp_rounds.size() + 1);
      }
      ++num_lp_rounds;
    } else if ( lp_stats[0].lp_rounds.size() == 0 ) {
      for ( LabelPropagationStats& lp : lp_stats ) {
        lp.lp_rounds.emplace_back(1);
      }
    }
  }

  LabelPropagationRoundStats label_propagation_round_stats(const int round) const {
    LabelPropagationRoundStats lp_round(round + 1);
    for ( const LabelPropagationStats& lp : lp_stats ) {
      lp_round += lp.lp_rounds[round];
    }
    return lp_round;
  }

  LabelPropagationStats label_propagation_stats() const {
    LabelPropagationStats stats;
    for ( const LabelPropagationStats& lp : lp_stats ) {
      stats += lp;
    }
    return stats;
  }

  FMStats fm_search_stats() const {
    FMStats stats;
    for ( const FMStats& fm : fm_stats ) {
      stats += fm;
    }
    return stats;
  }

  FlowStats flow_search_stats() const {
    FlowStats stats;
    for ( const FlowStats& flow : flow_stats ) {
      stats += flow;
    }
    return stats;
  }

  HypernodeID num_nodes;
  HyperedgeID num_edges;
  HypernodeID num_pins;
  size_t max_lp_rounds;
  size_t num_lp_rounds;
  vec<LabelPropagationStats> lp_stats;
  vec<FMStats> fm_stats;
  vec<FlowStats> flow_stats;
};

inline std::ostream & operator<< (std::ostream& str, const Search& search) {
  LabelPropagationStats lp = search.label_propagation_stats();
  FMStats fm = search.fm_search_stats();
  FlowStats flows = search.flow_search_stats();
  HyperedgeWeight total_improvement = lp.improvement + fm.improvement + flows.improvement;
  str << "num_nodes=" << search.num_nodes
      << " num_edges=" << search.num_edges
      << " num_pins=" << search.num_pins
      << " total_improvement=" << total_improvement
      << " lp_share_improvement=" << improvement_share(lp.improvement, total_improvement)
      << " fm_share_improvement=" << improvement_share(fm.improvement, total_improvement)
      << " flows_share_improvement=" << improvement_share(flows.improvement, total_improvement)
      << " num_lp_rounds=" << search.max_lp_rounds
      << " " << lp;
  vec<LabelPropagationRoundStats> lp_rounds;
  for ( size_t i = 0; i < search.num_lp_rounds; ++i ) {
    lp_rounds.emplace_back(search.label_propagation_round_stats(i));
  }
  for ( size_t i = search.num_lp_rounds; i < search.max_lp_rounds; ++i ) {
    lp_rounds.emplace_back(LabelPropagationRoundStats(i + 1));
  }
  for ( const LabelPropagationRoundStats& lp_round : lp_rounds ) {
    str << " " << lp_round
        << " lp_round_" << lp_round.cur_round << "_share_improvement="
        << improvement_share(lp_round.improvement, lp.improvement);
  }
  str << " " << fm
      << " " << flows;
  return str;
}

class RefinementStats {
  static constexpr bool debug = false;

 public:
  RefinementStats(const RefinementStats&) = delete;
  RefinementStats & operator= (const RefinementStats &) = delete;

  RefinementStats(RefinementStats&&) = delete;
  RefinementStats & operator= (RefinementStats &&) = delete;

  static RefinementStats& instance() {
    static RefinementStats instance;
    return instance;
  }

  void newSearch(const HypernodeID num_nodes,
                 const HyperedgeID num_edges,
                 const HypernodeID num_pins,
                 const bool is_main) {
    if ( is_main ) {
      _searches.emplace_back(num_nodes,num_edges,num_pins);
    }
  }

  Search& currentSearch() {
    return _searches.back();
  }

  void clear() {
    _searches.clear();
  }

  friend std::ostream & operator<< (std::ostream& str, const RefinementStats& stats);

 public:
  std::string _graph;
  PartitionID _k;
  double _epsilon;
  int _num_threads;
  int _seed;
  HyperedgeWeight _initial_km1;

 private:
  explicit RefinementStats() :
    _graph(),
    _k(kInvalidPartition),
    _epsilon(0),
    _num_threads(0),
    _seed(-1),
    _initial_km1(std::numeric_limits<HyperedgeWeight>::max()),
    _searches() {
    _searches.emplace_back(0,0,0);
  }


  vec<Search> _searches;
};


inline std::ostream & operator<< (std::ostream& str, const RefinementStats& stats) {
  HyperedgeWeight lp_improvement = 0;
  HyperedgeWeight fm_improvement = 0;
  HyperedgeWeight flows_improvement = 0;
  for ( const Search& search : stats._searches ) {
    str << "SEARCH_STATS"
        << " graph=" << stats._graph
        << " k=" << stats._k
        << " epsilon=" << stats._epsilon
        << " num_threads=" << stats._num_threads
        << " seed=" << stats._seed
        << " " << search << std::endl;
    lp_improvement += search.label_propagation_stats().improvement;
    fm_improvement += search.fm_search_stats().improvement;
    flows_improvement += search.flow_search_stats().improvement;
  }
  HyperedgeWeight total_improvement = lp_improvement + fm_improvement + flows_improvement;
  str << "REFINEMENT_STATS"
      << " graph=" << stats._graph
      << " k=" << stats._k
      << " epsilon=" << stats._epsilon
      << " num_threads=" << stats._num_threads
      << " seed=" << stats._seed
      << " total_improvement=" << total_improvement
      << " lp_improvement=" << lp_improvement
      << " fm_improvement=" << fm_improvement
      << " flows_improvement=" << flows_improvement
      << " lp_share_improvement=" << improvement_share(lp_improvement, total_improvement)
      << " fm_share_improvement=" << improvement_share(fm_improvement, total_improvement)
      << " flows_share_improvement=" << improvement_share(flows_improvement, total_improvement)
      << " relative_lp_improvement=" << improvement_share(stats._initial_km1, stats._initial_km1 - lp_improvement)
      << " relative_fm_improvement=" << improvement_share(stats._initial_km1, stats._initial_km1 - fm_improvement)
      << " relative_flows_improvement=" << improvement_share(stats._initial_km1, stats._initial_km1 - flows_improvement);
  return str;
}


}  // namespace utils
}  // namespace mt_kahypar
