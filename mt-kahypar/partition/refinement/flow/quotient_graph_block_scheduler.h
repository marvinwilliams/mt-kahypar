/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <list>
#include <tbb/spin_mutex.h>

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h"
#include "external_tools/kahypar/kahypar/datastructure/sparse_set.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/refinement/flow/flow_refiner.h"

namespace mt_kahypar {

template<typename TypeTraits>
class QuotientGraphBlockScheduler {
  typedef std::pair<PartitionID, PartitionID> edge;
  typedef std::pair<int, edge> scheduling_edge;
  using ConstIncidenceIterator = std::vector<edge>::const_iterator;
  using ConstCutHyperedgeIterator = std::vector<HyperedgeID>::const_iterator;
  using TBB = typename TypeTraits::TBB;

 private:
  using HyperGraph = typename TypeTraits::HyperGraph;

 public:
  QuotientGraphBlockScheduler(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _quotient_graph(),
    _round_edges(),
    _active_blocks(_context.partition.k, true),
    _locked_blocks(_context.partition.k, false),

    _block_pair_cut_he(context.partition.k,
                       std::vector<std::vector<HyperedgeID> >(context.partition.k,
                                                              std::vector<HyperedgeID>())),
    _schedule_mutex(),
    //_task_group(),
    _tasks_on_numa(TBB::instance().num_used_numa_nodes(), 0) { }

  QuotientGraphBlockScheduler(const QuotientGraphBlockScheduler&) = delete;
  QuotientGraphBlockScheduler(QuotientGraphBlockScheduler&&) = delete;
  QuotientGraphBlockScheduler& operator= (const QuotientGraphBlockScheduler&) = delete;
  QuotientGraphBlockScheduler& operator= (QuotientGraphBlockScheduler&&) = delete;

  void buildQuotientGraph() {
    std::set<edge> edge_list;
    for (const HyperedgeID& he : _hg.edges()) {
      if (_hg.connectivity(he) > 1) {
        for (const PartitionID& block0 : _hg.connectivitySet(he)) {
          for (const PartitionID& block1 : _hg.connectivitySet(he)) {
            if (block0 < block1) {
              edge_list.insert(std::make_pair(block0, block1));
              _block_pair_cut_he[block0][block1].push_back(he);
            }
          }
        }
      }
    }
    for (const edge& e : edge_list) {
      _quotient_graph.push_back(e);
    }
  }

  std::vector<scheduling_edge> getInitialParallelEdges(){
    std::vector<scheduling_edge> initialEdges;
    for(auto const edge:_quotient_graph){
      if(_active_blocks[edge.first] && _active_blocks[edge.second]){
        scheduling_edge sched_edge = std::make_pair(_hg.get_numa_node_of_blockpair(edge.first, edge.second),edge);
        _round_edges.push_back(sched_edge);
      }
    }

    size_t N = _round_edges.size();
    for (size_t i = 0; i < N; ++i) {
      const edge e = _round_edges[i].second;
      if (!_locked_blocks[e.first] && !_locked_blocks[e.second]){
        initialEdges.push_back(_round_edges[i]);
        _tasks_on_numa[_round_edges[i].first] ++;
        _locked_blocks[e.first] = true;
        _locked_blocks[e.second] = true;

        std::swap(_round_edges[i], _round_edges[N - 1]);
        _round_edges.pop_back();
        --i;
        --N;
      }
    }
    //reset active-array before each round
    //blocks are set active, if improvement was found
    _active_blocks.assign(_context.partition.k, false);

    return initialEdges;
  }

  std::vector<scheduling_edge> getNextBlockToSchedule(const scheduling_edge old_sched_edge){
    tbb::spin_mutex::scoped_lock lock{_schedule_mutex};

    //unlock the blocks
    _locked_blocks[old_sched_edge.second.first] = false;
    _locked_blocks[old_sched_edge.second.second] = false;

    //start new flow-calclation for blocks
    std::vector<scheduling_edge> new_sched_edges;
    size_t N = _round_edges.size();
    for (size_t i = 0; i < N; ++i) {
      auto sched_edge = _round_edges[i];
      const edge e = sched_edge.second;
      if (!_locked_blocks[e.first] && !_locked_blocks[e.second]){
        _locked_blocks[e.first] = true;
        _locked_blocks[e.second] = true;
        _tasks_on_numa[sched_edge.first] ++;

        new_sched_edges.push_back(sched_edge);

        std::swap(_round_edges[i], _round_edges[N - 1]);
        _round_edges.pop_back();
        --i;
        --N;
      }
    }
    _tasks_on_numa[old_sched_edge.first] --;
    return new_sched_edges;
  }

  size_t getNumberOfActiveTasks(){
    size_t tasks = 0;
    for(auto numa:_tasks_on_numa){
      tasks += numa;
    }
    return tasks;
  }

  /* tbb::task_group& get_task_group(){
    return _task_group;
  }*/

  void randomShuffleQoutientEdges() {
    utils::Randomize::instance().shuffleVector(_quotient_graph);
  }

  std::pair<ConstIncidenceIterator, ConstIncidenceIterator> qoutientGraphEdges() const {
    return std::make_pair(_quotient_graph.cbegin(), _quotient_graph.cend());
  }

  std::pair<ConstCutHyperedgeIterator, ConstCutHyperedgeIterator> blockPairCutHyperedges(const PartitionID block0, const PartitionID block1) {
    ASSERT(block0 < block1, V(block0) << " < " << V(block1));
    updateBlockPairCutHyperedges(block0, block1);

    ASSERT([&]() {
        std::set<HyperedgeID> cut_hyperedges;
        for (const HyperedgeID& he : _block_pair_cut_he[block0][block1]) {
          if (cut_hyperedges.find(he) != cut_hyperedges.end()) {
            LOG << "Hyperedge " << he << " is contained more than once!";
            return false;
          }
          cut_hyperedges.insert(he);
        }    
        for (const HyperedgeID& he : _hg.edges()) {
          if (_hg.pinCountInPart(he, block0) > 0 &&
              _hg.pinCountInPart(he, block1) > 0) {
            /*if (cut_hyperedges.find(he) == cut_hyperedges.end()) {
              LOG << V(_hg.pinCountInPart(he, block0));
              LOG << V(_hg.pinCountInPart(he, block1));
              LOG << V(_hg.initialNumNodes()) << V(_hg.currentNumNodes());
              LOG << V(he) << "should be inside the incidence set of"
                  << V(block0) << "and" << V(block1);
              return false;
            }*/
          } else {
            if (cut_hyperedges.find(he) != cut_hyperedges.end()) {
              LOG << V(_hg.pinCountInPart(he, block0));
              LOG << V(_hg.pinCountInPart(he, block1));
              LOG << V(he) << "shouldn't be inside the incidence set of"
                  << V(block0) << "and" << V(block1);
              return false;
            }
          }
        }
        return true;
      } (), "Cut hyperedge set between " << V(block0) << " and " << V(block1) << " is wrong!");

    return std::make_pair(_block_pair_cut_he[block0][block1].cbegin(),
                          _block_pair_cut_he[block0][block1].cend());
  }

  void changeNodePart(const HypernodeID hn, const PartitionID from, const PartitionID to) {
    if (from != to) {
      bool success = _hg.changeNodePart(hn, from, to);
      unused(success);
      ASSERT(success);

      for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
        if (_hg.pinCountInPart(he, to) == 1) {
          // This is not thread-safe
          // Can cause that a Hyperedge is missing in _block_pair_cut_he
          // Decided to accept the downside this causes as it happens very rarely and does not break the Algorithm
          for (const PartitionID& part : _hg.connectivitySet(he)) {
            if (to < part) {
              _block_pair_cut_he[to][part].push_back(he);
            } else if (to > part) {
              _block_pair_cut_he[part][to].push_back(he);
            }
          }
        }
      }
    }
  }

  void setActiveBlock(size_t blockId, bool active){
    _active_blocks[blockId] = active;
  }

  size_t getNumberOfActiveBlocks(){
    size_t active_blocks = 0;
    for(auto active:_active_blocks){
      if(active){
        active_blocks++;
      }
    }
    return active_blocks;
  }

 private:
  static constexpr bool debug = false;

  void updateBlockPairCutHyperedges(const PartitionID block0, const PartitionID block1) {
    kahypar::ds::FastResetFlagArray<> visited = kahypar::ds::FastResetFlagArray<>(_hg.initialNumEdges());
    size_t N = _block_pair_cut_he[block0][block1].size();
    for (size_t i = 0; i < N; ++i) {
      const HyperedgeID he = _block_pair_cut_he[block0][block1][i];
      if (_hg.pinCountInPart(he, block0) == 0 ||
          _hg.pinCountInPart(he, block1) == 0 ||
          visited[_hg.originalEdgeID(he)]) {
        std::swap(_block_pair_cut_he[block0][block1][i],
                  _block_pair_cut_he[block0][block1][N - 1]);
        _block_pair_cut_he[block0][block1].pop_back();
        --i;
        --N;
      }
      visited.set(_hg.originalEdgeID(he), true);
    }
  }

  HyperGraph& _hg;
  const Context& _context;
  std::vector<edge> _quotient_graph;
  //holds all eges, that are executed in that round (both blocks are active)
  std::vector<scheduling_edge> _round_edges;

  std::vector<bool> _active_blocks;
  std::vector<bool> _locked_blocks;
  

  // Contains the cut hyperedges for each pair of blocks.
  std::vector<std::vector<std::vector<HyperedgeID> > > _block_pair_cut_he;

  tbb::spin_mutex _schedule_mutex;
  std::vector<size_t> _tasks_on_numa;
  //tbb::task_group _task_group;
};
}  // namespace mt-kahypar
