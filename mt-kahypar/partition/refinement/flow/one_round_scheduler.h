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
#include <tbb/spin_rw_mutex.h>

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h"
#include "external_tools/kahypar/kahypar/datastructure/sparse_set.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"



namespace mt_kahypar {
using edge = std::pair<PartitionID, PartitionID>;
using scheduling_edge = std::pair<int, edge>;
using ConstIncidenceIterator = std::vector<edge>::const_iterator;
using ConstCutHyperedgeIterator = std::vector<HyperedgeID>::const_iterator;

class OneRoundScheduler {

 public:
  OneRoundScheduler(PartitionedHypergraph<>& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _quotient_graph(),
    _round_edges(),
    _finished(false),
    _current_rounds(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
    _active_blocks(1 , std::vector<tbb::atomic<bool>>(context.partition.k, false)),
    _active_blocks_mutex(),
    _max_round(0),
    _running_edges(context.partition.k, std::vector<bool>(context.partition.k, false)),
    _block_pair_cut_he(context.partition.k,
                       std::vector<std::vector<HyperedgeID> >(context.partition.k,
                                                              std::vector<HyperedgeID>())),
    _schedule_mutex(),
    _block_weights(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
    _rw_locks(context.partition.k),
    _tasks_on_block(context.partition.k, 0),
    _node_lock(hypergraph.initialNumNodes(), false) { }

  OneRoundScheduler(const OneRoundScheduler&) = delete;
  OneRoundScheduler(OneRoundScheduler&&) = delete;
  OneRoundScheduler& operator= (const OneRoundScheduler&) = delete;
  OneRoundScheduler& operator= (OneRoundScheduler&&) = delete;

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
    //push edges from active blocks in _round_edges
    for(auto const edge : _quotient_graph) {
        _round_edges.push_back(edge);
    }
    std::vector<edge> initial_edges;
    const size_t num_threads = TBBNumaArena::instance().total_number_of_threads();
    for (size_t i = 0; i < num_threads; i++){
      edge e = getMostIndependentEdge();
      if(e.first != kInvalidPartition && e.second != kInvalidPartition) {
        initial_edges.push_back(e);
      }
    }
    return initial_edges;
  }



  void scheduleNextBlocks(const edge old_edge, tbb::parallel_do_feeder<edge>& feeder){
    utils::Timer::instance().start_timer("schedNext", "Scheduling next block ", true);

    {tbb::spin_mutex::scoped_lock lock{_schedule_mutex};

    _tasks_on_block[old_edge.first] --;
    _tasks_on_block[old_edge.second] --;
    _running_edges[old_edge.first][old_edge.second] = false;

    size_t round = _current_rounds[old_edge.first][old_edge.second];
    if(_active_blocks[round][old_edge.first] && _active_blocks[round][old_edge.second]){
        _round_edges.push_back(old_edge);
    }

    if(getNumberOfActiveTasks() == 0 && _round_edges.empty())
        _finished = true;
    }//unlock here

    while(!_finished){
        tbb::spin_mutex::scoped_lock lock{_schedule_mutex};
        edge e = getMostIndependentEdge();
        if(e.first != kInvalidPartition && e.second != kInvalidPartition) {
        feeder.add(e);
        utils::Timer::instance().stop_timer("maxFlow");
        return;
        }
    }
    utils::Timer::instance().stop_timer("maxFlow");
  }

  /**
  * try to aquire a Hypernode. Return true on success, false if already aquired.
  * Not aquired Nodes have value 0.
  * Aquired Nodes have the value (block_0 * k) + block_1 to store the blocks of the flow calculation holding the Node
  **/
  bool tryAquireNode(HypernodeID node, const int blocks_idx){
    bool already_aquired = _node_lock[node].compare_and_swap(blocks_idx, 0);
    return !already_aquired;
  }

  bool isAquired(HypernodeID node){
    return _node_lock[node];
  }

  void releaseNode(HypernodeID node){
    ASSERT(_node_lock[node] > 0, "Tryed to release Node that is not aquired!");
    _node_lock[node] = 0;
  }

  void randomShuffleQoutientEdges() {
    utils::Randomize::instance().shuffleVector(_quotient_graph);
  }

  std::pair<ConstCutHyperedgeIterator, ConstCutHyperedgeIterator> blockPairCutHyperedges(const PartitionID block0, const PartitionID block1) {
    ASSERT(block0 < block1, V(block0) << " < " << V(block1));
    updateBlockPairCutHyperedges(block0, block1);

    return std::make_pair(_block_pair_cut_he[block0][block1].cbegin(),
                          _block_pair_cut_he[block0][block1].cend());
  }

  template<typename F>
  void changeNodePart(const HypernodeID hn, const PartitionID from, const PartitionID to, const F& objective_delta) {
    if (from != to) {
      bool success = _hg.changeNodePart(hn, from, to, objective_delta);
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

  void setBlocksActive(PartitionID block_0, PartitionID block_1){
    size_t round = _current_rounds[block_0][block_1];
    bool old_block_0 = _active_blocks[round][block_0].fetch_and_store(true);
    bool old_block_1 = _active_blocks[round][block_1].fetch_and_store(true);

    std::vector<edge> edges_to_schedule;
    
    //if block_0 was not active before, add the eges for the next round
    if(!old_block_0){
        tbb::spin_mutex::scoped_lock lock{_active_blocks_mutex};
        for (int i = 0; i < _context.partition.k; i++){
            if(_active_blocks[round][i] && i != block_0 && i != block_1){
                if(i < block_0)
                    edges_to_schedule.push_back(std::make_pair(i, block_0));
                else
                    edges_to_schedule.push_back(std::make_pair(block_0, i));
            }
        } 
    }
    //if block_1 was not active before, add the eges for the next round
    if(!old_block_1){
        tbb::spin_mutex::scoped_lock lock{_active_blocks_mutex};
        for (int i = 0; i < _context.partition.k; i++){
            if(_active_blocks[round][i] && i != block_0 && i != block_1){
                if(i < block_1)
                    edges_to_schedule.push_back(std::make_pair(i, block_1));
                else
                    edges_to_schedule.push_back(std::make_pair(block_1, i));
            }
        } 
    }
    if(!edges_to_schedule.empty()){
        tbb::spin_mutex::scoped_lock lock{_schedule_mutex};
        for(auto edge:edges_to_schedule){
            for(auto r_edge:_round_edges){
                if(r_edge == edge)
                    goto edge_already_contained;
            }
            if(_running_edges[edge.first][edge.second])
                goto edge_already_contained;

            _round_edges.push_back(edge);
            edge_already_contained:;
        }
    }
  }

  //return zero to stop while loop in flow refiner after one iteration
  size_t getNumberOfActiveBlocks(){
    return 0;
  }

  /**
   * Initialise _blockweights to deal with imbalance in parallel environment.
   * Flow calculations aquire the weight of the Hypernodes they hold from both blocks and save these weights to make it possible for other blocks to calculate and optimize the imbalance.
   * After a calculation, the modified weight gets written back to the block to make them available again. The operations are saved using a r/w lock.
   * This method is not safe to keep a balanced Hypergraph. When two calculations try to correct an imbalance by increasing a blockweight of the same block concurrently, the imbalance can
   * exceed epsilon. In practice this was never observed. The imbalance would be corrected in a following label propagation step.
   **/
  void init_block_weights(){
    for (int i = 0; i < _context.partition.k; i++){
      for (int j = 0; j < _context.partition.k; j++){
        if(j==i){
          _block_weights[i][j] = _hg.partWeight(i);
        }else{
          _block_weights[i][j] = 0;
        }
      }
    }
  }

  void aquire_block_weight(size_t block_to_aquire, size_t other_block, size_t amount){
    tbb::spin_rw_mutex::scoped_lock lock{_rw_locks[block_to_aquire], true};

    _block_weights[block_to_aquire][other_block] = amount;
    _block_weights[block_to_aquire][block_to_aquire] -= amount;
  }

  void release_block_weight(size_t block_to_release, size_t other_block, size_t amount){
    tbb::spin_rw_mutex::scoped_lock lock{_rw_locks[block_to_release], true};

    _block_weights[block_to_release][other_block] = 0;
    _block_weights[block_to_release][block_to_release] += amount;
  }

  size_t get_not_aquired_weight(PartitionID block, PartitionID other_block){
    tbb::spin_rw_mutex::scoped_lock lock{_rw_locks[block], false};

    size_t weight = 0;
    for (int i = 0; i < _context.partition.k; i++){
      if(i != other_block){
        weight += _block_weights[block][i];
      }
    }
    return weight;
  }

  std::pair<HypernodeWeight, HypernodeWeight> get_aquired_part_weight(PartitionID  block_0, PartitionID block_1){
    return std::make_pair(_block_weights[block_0][block_1], _block_weights[block_1][block_0]);
  }

  bool is_block_overlap(HypernodeID node,const PartitionID block_0, const PartitionID block_1){
    int blocks_idx = _node_lock[node];
    PartitionID other_block_0 = blocks_idx / _context.partition.k;
    PartitionID other_block_1 = blocks_idx % _context.partition.k;
    return (block_0 == other_block_0 || block_0 == other_block_1 || block_1 == other_block_0 || block_1 == other_block_1);
  }

 protected:
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

  edge getMostIndependentEdge(){
    bestEdge best_edge;
    size_t independence = std::numeric_limits<size_t>::max();

    for(size_t i = 0; i < _round_edges.size(); i++){
        edge e = _round_edges[i];
        size_t temp_independence = std::max(_tasks_on_block[e.first], _tasks_on_block[e.second]);
        if(temp_independence < independence){
        independence = temp_independence;
        best_edge = bestEdge{e, i};
        }
    }

    if(independence < std::numeric_limits<size_t>::max()) {
        _tasks_on_block[best_edge.e.first] ++;
        _tasks_on_block[best_edge.e.second] ++;
        removeElement(best_edge.index, _round_edges);
        size_t round = ++_current_rounds[best_edge.e.first][best_edge.e.second];
        if(round > _max_round){
            _active_blocks.push_back(std::vector<tbb::atomic<bool>>(_context.partition.k, false));
        }
        _running_edges[best_edge.e.first][best_edge.e.second] = true;
        return best_edge.e;
    }

    return std::make_pair(kInvalidPartition, kInvalidPartition);
  }

  size_t getNumberOfActiveTasks(){
      size_t tasks = 0;
      for (int i = 0; i < _context.partition.k; i++){
          tasks += _tasks_on_block[i];
      }
      return tasks;
  }

  struct bestEdge{
    edge e;
    size_t index;
  };

  template<typename T>
  void removeElement(size_t index, std::vector<T> &vector){
    size_t N = vector.size();
    std::swap(vector[index], vector[N - 1]);
    vector.pop_back();
  }

  PartitionedHypergraph<>& _hg;
  const Context& _context;
  std::vector<edge> _quotient_graph;
  //holds all eges, that are executed in that round (both blocks are active)
  std::vector<edge> _round_edges;
  tbb::atomic<bool> _finished;
  std::vector<std::vector<size_t>> _current_rounds;
  std::vector<std::vector<tbb::atomic<bool>>> _active_blocks;
  tbb::spin_mutex _active_blocks_mutex;
  size_t _max_round;
  std::vector<std::vector<bool>> _running_edges;

  // Contains the cut hyperedges for each pair of blocks.
  std::vector<std::vector<std::vector<HyperedgeID> > > _block_pair_cut_he;

  tbb::spin_mutex _schedule_mutex;
  std::vector<std::vector<size_t>> _block_weights;
  std::vector<tbb::spin_rw_mutex> _rw_locks;
  std::vector<tbb::atomic<int>> _tasks_on_block;
  std::vector<tbb::atomic<int>> _node_lock;
};

}  // namespace mt-kahypar
