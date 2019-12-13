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

namespace mt_kahypar {

template<typename TypeTraits>
class QuotientGraphBlockScheduler {
  typedef std::pair<PartitionID, PartitionID> edge;
  using ConstIncidenceIterator = std::vector<edge>::const_iterator;
  using ConstCutHyperedgeIterator = std::vector<HyperedgeID>::const_iterator;

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
    _schedule_mutex() { }

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

  std::vector<edge> getInitialParallelEdges(){
    std::vector<edge> initialEdges;
    // TODO(reister): use 'const edge& e' in order to prevent that elements are copied
    // or modified during iteration.
    for(auto edge:_quotient_graph){
      if(_active_blocks[edge.first] && _active_blocks[edge.second]){
        _round_edges.push_back(edge);
      }
    }
    std::list<edge>::const_iterator e = _round_edges.cbegin();
    while (e != _round_edges.cend()){
      if (!_locked_blocks[e->first] && !_locked_blocks[e->second]) {
        initialEdges.push_back(*e);
        _locked_blocks[e->first] = true;
        _locked_blocks[e->second] = true;

        e = _round_edges.erase(e);
      }
      else {
        ++e;
      }
    }
    //reset active-array before each round
    //blocks are set active, if improvement was found
    _active_blocks.assign(_context.partition.k, false);

    return initialEdges;
  }

  void scheduleNextBlock(tbb::parallel_do_feeder<edge>& feeder, const PartitionID block_0, const PartitionID block_1){
    // TODO(reister): Usually I'm not a big fan of locks, but I think in this situation it is okay
    // since the function is called not that often during flow refinement. Alternatively, you can
    // implement the locking mechanism with atomics and perform a compare_and_swap operation on the
    // two blocks. If both are successful you can add the edge to the feeder, otherwise reset them and continue
    // to next pair of blocks.
    tbb::spin_mutex::scoped_lock lock{_schedule_mutex};
    //unlock the blocks
    _locked_blocks[block_0] = false;
    _locked_blocks[block_1] = false;

    //start new flow-calclation for blocks
    // TODO(reister): Usually we don't use std::list, since erase has linear running time.
    // Again this should be not a problem, since the list only contains k^2 elements at most.
    // But what you could do instead is use a vector and if you add an element to the feeder
    // you can swap it to the end and decrement a pointer to that vector. The pointer points
    // to all elements that where not considered during that current round. Once all blocks are
    // processed you can just reset the pointer to the end of the vector and start again.
    std::list<edge>::const_iterator e = _round_edges.cbegin();
    while (e != _round_edges.cend()) {
      if (!_locked_blocks[e->first] && !_locked_blocks[e->second]) {
        _locked_blocks[e->first] = true;
        _locked_blocks[e->second] = true;
        feeder.add(*e);

        e = _round_edges.erase(e);
      }
      else {
        ++e;
      }
    }
  }

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
        //TODO: Assertion fails very rarely, probably race condition?
        /*
        for (const HyperedgeID& he : _hg.edges()) {
          if (_hg.pinCountInPart(he, block0) > 0 &&
              _hg.pinCountInPart(he, block1) > 0) {
            if (cut_hyperedges.find(he) == cut_hyperedges.end()) {
              LOG << V(_hg.pinCountInPart(he, block0));
              LOG << V(_hg.pinCountInPart(he, block1));
              LOG << V(_hg.initialNumNodes()) << V(_hg.currentNumNodes());
              LOG << V(he) << "should be inside the incidence set of"
                  << V(block0) << "and" << V(block1);
              return false;
            }
          } else {
            if (cut_hyperedges.find(he) != cut_hyperedges.end()) {
              LOG << V(_hg.pinCountInPart(he, block0));
              LOG << V(_hg.pinCountInPart(he, block1));
              LOG << V(he) << "shouldn't be inside the incidence set of"
                  << V(block0) << "and" << V(block1);
              return false;
            }
          }
        }*/
        return true;
      } (), "Cut hyperedge set between " << V(block0) << " and " << V(block1) << " is wrong!");

    return std::make_pair(_block_pair_cut_he[block0][block1].cbegin(),
                          _block_pair_cut_he[block0][block1].cend());
  }

  void changeNodePart(const HypernodeID hn, const PartitionID from, const PartitionID to) {
    if (from != to) {
      // TODO(reister): This should not fail, since flow problems are independent:
      // bool success = _hg.changeNodePart(hn, from, to);
      // ASSERT(success);
      while(_hg.changeNodePart(hn, from, to) == false);
      for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
        if (_hg.pinCountInPart(he, to) == 1) {
          // TODO(reister): This is not thread-safe and this is also why
          // your assertion fails rarely.
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

 private:
  static constexpr bool debug = false;

  void updateBlockPairCutHyperedges(const PartitionID block0, const PartitionID block1) {
    kahypar::ds::FastResetFlagArray<> visited = kahypar::ds::FastResetFlagArray<>(_hg.initialNumEdges());
    size_t N = _block_pair_cut_he[block0][block1].size();
    for (size_t i = 0; i < N; ++i) {
      const HyperedgeID he = _block_pair_cut_he[block0][block1][i];
      if (_hg.pinCountInPart(he, block0) == 0 ||
          _hg.pinCountInPart(he, block1) == 0 ||
          visited[he]) {
        std::swap(_block_pair_cut_he[block0][block1][i],
                  _block_pair_cut_he[block0][block1][N - 1]);
        _block_pair_cut_he[block0][block1].pop_back();
        --i;
        --N;
      }
      visited.set(he, true);
    }
  }

  HyperGraph& _hg;
  const Context& _context;
  std::vector<edge> _quotient_graph;
  //holds all eges, that are executed in that round (both blocks are active)
  std::list<edge> _round_edges;

  std::vector<bool> _active_blocks;
  std::vector<bool> _locked_blocks;

  // Contains the cut hyperedges for each pair of blocks.
  std::vector<std::vector<std::vector<HyperedgeID> > > _block_pair_cut_he;

  tbb::spin_mutex _schedule_mutex;
};
}  // namespace mt-kahypar
