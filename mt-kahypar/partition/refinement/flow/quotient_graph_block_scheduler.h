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
#include "mt-kahypar/partition/refinement/flow/flow_refiner.h"



namespace mt_kahypar {
using edge = std::pair<PartitionID, PartitionID>;
using scheduling_edge = std::pair<int, edge>;
using ConstIncidenceIterator = std::vector<edge>::const_iterator;
using ConstCutHyperedgeIterator = std::vector<HyperedgeID>::const_iterator;

template <typename TypeTraits = Mandatory, 
          typename Derived = Mandatory>
class SchedulerBase {

 private:
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;

 public:
  SchedulerBase(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _num_numa_nodes(TypeTraits::TBB::instance().num_used_numa_nodes()),
    _quotient_graph(_num_numa_nodes),
    _round_edges(_num_numa_nodes),
    _active_blocks(_context.partition.k, true),
    _locked_blocks(_context.partition.k, false),

    _block_pair_cut_he(context.partition.k,
                       std::vector<std::vector<HyperedgeID> >(context.partition.k,
                                                              std::vector<HyperedgeID>())),
    _schedule_mutex(),
    _tasks_on_numa(_num_numa_nodes, 0),
    _empty_numas(),
    _block_weights(context.partition.k, std::vector<size_t>(context.partition.k, 0)),
    _rw_locks(context.partition.k) { }

  SchedulerBase(const SchedulerBase&) = delete;
  SchedulerBase(SchedulerBase&&) = delete;
  SchedulerBase& operator= (const SchedulerBase&) = delete;
  SchedulerBase& operator= (SchedulerBase&&) = delete;

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

   std::vector<std::vector<int>> nodes_per_numa(_context.partition.k,std::vector<int>(_num_numa_nodes, 0));
    for(const HypernodeID& hn:_hg.nodes()){
      nodes_per_numa[_hg.partID(hn)][common::get_numa_node_of_vertex(hn)] ++;
    }


    for (const edge& e : edge_list) {
      std::vector<size_t> block_nodes_per_numa(_num_numa_nodes, 0);
      for(size_t i = 0; i < _num_numa_nodes; i++){
        block_nodes_per_numa[i] = nodes_per_numa[e.first][i] + nodes_per_numa[e.second][i];
      }
      int numa_node = distance(std::begin(block_nodes_per_numa), max_element(std::begin(block_nodes_per_numa), std::end(block_nodes_per_numa)));
      _quotient_graph[numa_node].push_back(e);
    }
  }

  std::vector<std::vector<edge>> getInitialParallelEdges(){
    //push edges from active blocks in _round_edges
    for (size_t i = 0; i < this->_num_numa_nodes; i++){
      for(auto const edge:this->_quotient_graph[i]){
        if(this->_active_blocks[edge.first] && this->_active_blocks[edge.second]){
          this->_round_edges[i].push_back(edge);
        }
      }
    }
    return static_cast<Derived*>(this)->getInitialParallelEdgesImpl();
  }

  std::vector<scheduling_edge> scheduleNextBlocks(const edge old_edge, const int node, tbb::parallel_do_feeder<edge>& feeder){
    return static_cast<Derived*>(this)->scheduleNextBlocksImpl(old_edge, node, feeder);
  }

  size_t getNumberOfActiveTasks(){
    return static_cast<Derived*>(this)->getNumberOfActiveTasksImpl();
  }

  bool tryAquireNode(HypernodeID node){
    return static_cast<Derived*>(this)->tryAquireNodeImpl(node);
  }

  void releaseNode(HypernodeID node){
    static_cast<Derived*>(this)->releaseNodeImpl(node);
  }

  void randomShuffleQoutientEdges() {
    for (size_t numa = 0; numa < _num_numa_nodes; numa++){
      utils::Randomize::instance().shuffleVector(_quotient_graph[numa]);
    }
  }

  

  /*std::pair<ConstIncidenceIterator, ConstIncidenceIterator> qoutientGraphEdges() const {
    return std::make_pair(_quotient_graph.cbegin(), _quotient_graph.cend());
  }*/

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
            // other threads can move pins after the update
            /*
            if (cut_hyperedges.find(he) != cut_hyperedges.end()) {
              LOG << V(_hg.pinCountInPart(he, block0));
              LOG << V(_hg.pinCountInPart(he, block1));
              LOG << V(he) << "shouldn't be inside the incidence set of"
                  << V(block0) << "and" << V(block1);
              return false;
            }*/
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

  //TODO: add read write lock or make sth atomic
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

 protected:
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

  template<typename T>
  void removeElement(size_t index, std::vector<T> &vector){
    size_t N = vector.size();
    std::swap(vector[index], vector[N - 1]);
    vector.pop_back();
  }

  HyperGraph& _hg;
  const Context& _context;
  const size_t _num_numa_nodes;
  std::vector<std::vector<edge>> _quotient_graph;
  //holds all eges, that are executed in that round (both blocks are active)
  std::vector<std::vector<edge>> _round_edges;

  std::vector<bool> _active_blocks;
  std::vector<bool> _locked_blocks;


  // Contains the cut hyperedges for each pair of blocks.
  std::vector<std::vector<std::vector<HyperedgeID> > > _block_pair_cut_he;

  tbb::spin_mutex _schedule_mutex;
  std::vector<size_t> _tasks_on_numa;
  std::vector<int> _empty_numas;

  std::vector<std::vector<size_t>> _block_weights;
  std::vector<tbb::spin_rw_mutex> _rw_locks;
};

template <typename TypeTraits>
class MatchingScheduler : public SchedulerBase<TypeTraits, MatchingScheduler<TypeTraits>> {
  using Base = SchedulerBase<TypeTraits, MatchingScheduler<TypeTraits>>;
  public:
  //using base constructor
  using Base::Base;
  
  std::vector<std::vector<edge>> getInitialParallelEdgesImpl(){
    std::vector<std::vector<edge>> initialEdges(this->_num_numa_nodes);

    //split work round-robin style among numa nodes

    //vector for right access after Edges get scheduled
    std::vector<int> scheduled_on_numa(this->_num_numa_nodes, 0);
    int i = 0;
    bool numa_has_edges = true;
    while(numa_has_edges){
      numa_has_edges = false;
      for (size_t numa = 0; numa < this->_num_numa_nodes; numa++){
        size_t numa_i = i - scheduled_on_numa[numa];
        if(numa_i < this->_round_edges[numa].size()){
          numa_has_edges = true;
          const edge e = this->_round_edges[numa][numa_i];
          if (!this->_locked_blocks[e.first] && !this->_locked_blocks[e.second]){
            initialEdges[numa].push_back(this->_round_edges[numa][numa_i]);
            this->_tasks_on_numa[numa] ++;
            scheduled_on_numa[numa]++;
            this->_locked_blocks[e.first] = true;
            this->_locked_blocks[e.second] = true;

            //delete edge in round_edges
            size_t N = this->_round_edges[numa].size();
            std::swap(this->_round_edges[numa][numa_i], this->_round_edges[numa][N - 1]);
            this->_round_edges[numa].pop_back();
          }
        }
      }
      i++;
    }

    //test if there is a NUMA node with no initial work but with work left in this round.
    //add to _empty_numas list to start new parallel do on them in schedulenextblock if possible
    for (size_t numa = 0; numa < this->_num_numa_nodes; numa++){
      if(initialEdges[numa].size() == 0 && this->_round_edges[numa].size() > 0){
        this->_empty_numas.push_back(numa);
      }
    }

    //reset active-array before each round
    //blocks are set active, if improvement was found
    this->_active_blocks.assign(this->_context.partition.k, false);

    return initialEdges;
  }

  std::vector<scheduling_edge> scheduleNextBlocksImpl(const edge old_edge, const int node, tbb::parallel_do_feeder<edge>& feeder){
    tbb::spin_mutex::scoped_lock lock{this->_schedule_mutex};
    //unlock the blocks
    this->_locked_blocks[old_edge.first] = false;
    this->_locked_blocks[old_edge.second] = false;

    //start new flow-calclation for blocks on the same numa node
    size_t N = this->_round_edges[node].size();
    for (size_t i = 0; i < N; ++i) {
      auto e = this->_round_edges[node][i];
      if (!this->_locked_blocks[e.first] && !this->_locked_blocks[e.second]){
        this->_locked_blocks[e.first] = true;
        this->_locked_blocks[e.second] = true;
        this->_tasks_on_numa[node] ++;

        //start new
        feeder.add(e);

        std::swap(this->_round_edges[node][i], this->_round_edges[node][N - 1]);
        this->_round_edges[node].pop_back();
        --i;
        --N;
      }
    }

    std::vector<scheduling_edge> sched_edges;
    //build vector to start new parallel do's on empty numa nodes with work left
    size_t numa_N = this->_empty_numas.size();
    for (size_t j = 0; j < numa_N; j++){
      int numa = this->_empty_numas[j];
      size_t N = this->_round_edges[numa].size();
      for (size_t i = 0; i < N; ++i) {
        auto e = this->_round_edges[numa][i];
        if (!this->_locked_blocks[e.first] && !this->_locked_blocks[e.second]){
          this->_locked_blocks[e.first] = true;
          this->_locked_blocks[e.second] = true;
          this->_tasks_on_numa[numa] ++;

          //add new sched_edge to start
          sched_edges.push_back(std::make_pair(numa, e));

          //remove round_edge
          std::swap(this->_round_edges[numa][i], this->_round_edges[numa][N - 1]);
          this->_round_edges[numa].pop_back();
          i = N;
          //remove numa from empty numas
          std::swap(this->_empty_numas[j], this->_empty_numas[numa_N - 1]);
          this->_empty_numas.pop_back();
          --j;
          --numa_N;
        }
      }
    }

    // check if this is the last running task on the numa node and if there are still tasks
    // left to do on the numa node. Try to steal a task from another numa-node to keep running
    // if not possible add this numa node to _empty_numas
    if(this->_tasks_on_numa[node] == 1 && this->_round_edges[node].size() > 0){
      for (size_t numa = 0; numa < this->_num_numa_nodes; numa++){
      size_t N = this->_round_edges[numa].size();
        for (size_t i = 0; i < N; ++i) {
          auto e = this->_round_edges[numa][i];
          if (!this->_locked_blocks[e.first] && !this->_locked_blocks[e.second]){
            this->_locked_blocks[e.first] = true;
            this->_locked_blocks[e.second] = true;
            this->_tasks_on_numa[node] ++;

            //start new
            feeder.add(e);

            std::swap(this->_round_edges[numa][i], this->_round_edges[numa][N - 1]);
            this->_round_edges[numa].pop_back();
            --i;
            --N;
            // return if a task could be stolen
            this->_tasks_on_numa[node] --;
            return sched_edges;
          }
        }
      }
      LOG << "numa got empty" << V(node);
      this->_empty_numas.push_back(node);
    }

    this->_tasks_on_numa[node] --;
    return sched_edges;
  }

  size_t getNumberOfActiveTasksImpl(){
    size_t tasks = 0;
    for(auto numa:this->_tasks_on_numa){
      tasks += numa;
    }
    return tasks;
  }

  bool tryAquireNodeImpl(HypernodeID node){
    unused(node);
    return true;
  }

  void releaseNodeImpl(HypernodeID node){
    unused(node);
  }

};

template <typename TypeTraits>
class OptScheduler : public SchedulerBase<TypeTraits, OptScheduler<TypeTraits>> {
  using Base = SchedulerBase<TypeTraits, OptScheduler<TypeTraits>>;
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using HwTopology = typename TypeTraits::HwTopology;
  public:
  OptScheduler(HyperGraph& hypergraph, const Context& context) :
    Base(hypergraph, context),
    _tasks_on_block(context.partition.k, 0),
    _node_lock(hypergraph.initialNumNodes(), false){ }

  std::vector<std::vector<edge>> getInitialParallelEdgesImpl(){
    std::vector<std::vector<edge>> initial_edges (this->_num_numa_nodes);
    for (size_t numa = 0; numa < this->_num_numa_nodes; numa++){
      //TODO: also use Hyperthreads?
      size_t num_cpus = HwTopology::instance().num_cores_on_numa_node(numa);
      for (size_t i = 0; i < num_cpus; i++){
        edge e = getMostIndependentEdge(numa);
        if(e.first >= 0){
          initial_edges[numa].push_back(e);
        }
      }
    }
    //reset active-array before each round
    //blocks are set active, if improvement was found
    this->_active_blocks.assign(this->_context.partition.k, false);
    return initial_edges;
  }

  std::vector<scheduling_edge> scheduleNextBlocksImpl(const edge old_edge, const int node, tbb::parallel_do_feeder<edge>& feeder){
    tbb::spin_mutex::scoped_lock lock{this->_schedule_mutex};

    _tasks_on_block[old_edge.first] --;
    _tasks_on_block[old_edge.second] --;

    edge e = getMostIndependentEdge(node);
    if(e.first >= 0)
      feeder.add(e);
    return std::vector<scheduling_edge>(0);
  }

  size_t getNumberOfActiveTasksImpl(){
    size_t tasks = 0;
    for(auto block:_tasks_on_block){
      tasks += block;
    }
    return tasks;
  }

  bool tryAquireNodeImpl(HypernodeID node){
    bool already_aquired = _node_lock[node].compare_and_swap(true, false);
    return !already_aquired;
  }

  void releaseNodeImpl(HypernodeID node){
    ASSERT(_node_lock[node] == true, "Tryed to release Node that is not aquired!");
    _node_lock[node] = false;
  }

  private:
  struct bestEdge{
    edge e;
    size_t numa;
    size_t index;
  };

  edge getMostIndependentEdge(size_t own_numa){
    bestEdge best_edge;
    size_t independence = std::numeric_limits<size_t>::max();

    //try to find edge from the same numa
    for(size_t i = 0; i < this->_round_edges[own_numa].size(); i++){
      edge e = this->_round_edges[own_numa][i];
      size_t temp_independence = std::max(_tasks_on_block[e.first], _tasks_on_block[e.second]);
      if(temp_independence < independence){
        independence = temp_independence;
        best_edge = bestEdge{e, own_numa, i};
      } 
    }
    if(independence < std::numeric_limits<size_t>::max()){
      _tasks_on_block[best_edge.e.first] ++;
      _tasks_on_block[best_edge.e.second] ++;
      this->removeElement(best_edge.index, this->_round_edges[best_edge.numa]);
      return best_edge.e;
    }

    //search all numas
    for (size_t numa = 0; numa < this->_num_numa_nodes; numa++){
      for(size_t i = 0; i < this->_round_edges[numa].size(); i++){
        edge e = this->_round_edges[numa][i];
        size_t temp_independence = std::max(_tasks_on_block[e.first], _tasks_on_block[e.second]);
        if(temp_independence < independence){
          independence = temp_independence;
          best_edge = bestEdge{e, numa, i};
        } 
      }
    }
    if(independence < std::numeric_limits<size_t>::max()){
      _tasks_on_block[best_edge.e.first] ++;
      _tasks_on_block[best_edge.e.second] ++;
      this->removeElement(best_edge.index, this->_round_edges[best_edge.numa]);
      return best_edge.e;
    } else {
      return std::make_pair(-1, -1);
    }
  }

  std::vector<size_t> _tasks_on_block;
  std::vector<tbb::atomic<bool>> _node_lock;
};

}  // namespace mt-kahypar
