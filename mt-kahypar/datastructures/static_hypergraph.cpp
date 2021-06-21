/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "static_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/memory_tree.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

namespace mt_kahypar::ds {

  StaticHypergraph StaticHypergraph::contract(vec<HypernodeID>& clusters) {
    auto& timer = utils::Timer::instance();
    timer.start_timer("hypergraph_contraction","Contraction");


    if (!_tmp_contraction_buffer) {
      allocateTmpContractionBuffer();
    }

    timer.start_timer("compactify","compactify");

    auto& mapping = _tmp_contraction_buffer->mapping;
    tbb::parallel_for(0U, initialNumNodes(), [&](HypernodeID u) { mapping[u] = 0;});
    tbb::parallel_for(0U, initialNumNodes(), [&](HypernodeID u) { mapping[clusters[u]] = 1; });
    parallel_prefix_sum(mapping.begin(), mapping.begin() + initialNumNodes(), mapping.begin(), std::plus<>(), 0);
    HypernodeID num_coarse_nodes = mapping[initialNumNodes() - 1];
    _tmp_contraction_buffer->num_coarse_nodes = num_coarse_nodes;
    // apply mapping to cluster IDs. subtract one because prefix sum is inclusive
    tbb::parallel_for(0U, initialNumNodes(), [&](HypernodeID u) {
      clusters[u] = nodeIsEnabled(u) ? mapping[clusters[u]] - 1 : kInvalidHypernode;
    });

    timer.stop_timer("compactify");
    timer.start_timer("generate pinlists","generate pinlists");

    auto get_cluster = [&](HypernodeID u) { assert(u < clusters.size()); return clusters[u]; };
    auto cs2 = [](const size_t x) { return x * x; };

    auto& coarse_pin_lists = _tmp_contraction_buffer->coarse_pin_lists;
    auto& permutation = _tmp_contraction_buffer->permutation;
    tbb::enumerable_thread_specific<boost::dynamic_bitset<>>& local_maps = _tmp_contraction_buffer->local_maps;

    // map coarse pin lists and insert into hash map for work distribution
    tbb::parallel_for(0U, initialNumEdges(), [&](HyperedgeID he) {
      vec<HypernodeID>& pin_list = coarse_pin_lists[he];
      pin_list.clear();
      if (!edgeIsEnabled(he)) {
        permutation[he] = { std::numeric_limits<size_t>::max(), 0, he,false };
        return;
      }
      boost::dynamic_bitset<>& contained = local_maps.local();
      pin_list.reserve(edgeSize(he) / 2);
      for (HypernodeID v : pins(he)) {
        const HypernodeID cv = get_cluster(v);
        if (cv != kInvalidHypernode && !contained[cv]) {
          contained.set(cv);
          pin_list.push_back(cv);
        }
        // pin_list.push_back(cv);
      }
      for (HypernodeID v : pin_list) {
        contained.reset(v);
      }

      // std::sort(pin_list.begin(), pin_list.end());
      // pin_list.erase(std::unique(pin_list.begin(), pin_list.end()), pin_list.end());
      //if (!pin_list.empty() && pin_list.back() == kInvalidHypernode)
      //  pin_list.pop_back();
      if (pin_list.size() > 1) {
        size_t edge_hash = 420; for (const HypernodeID v : pin_list) { edge_hash += cs2(v); }
        permutation[he] = { edge_hash, pin_list.size(), he, true };
        // net_map.insert(edge_hash, ContractedHyperedgeInformation{ he, edge_hash, pin_list.size(), true });
      } else {
        pin_list.clear();   // globally mark net as removed
        permutation[he] = { std::numeric_limits<size_t>::max(), 0, he, false };
      }
    });

    timer.stop_timer("generate pinlists");
    timer.start_timer("identical net detection","identical net detection");

    timer.start_timer("sort", "sort");

    tbb::parallel_sort(permutation.begin(), permutation.begin() + initialNumEdges());

    timer.stop_timer("sort");
    timer.start_timer("detect deduplicates", "detect deduplicates");

    auto& coarse_edge_weights = _tmp_contraction_buffer->coarse_edge_weights;
    HypernodeID num_coarse_nets = 0;
    size_t num_coarse_pins = 0;

    // identical net detection
    doParallelForAllEdges([&](HyperedgeID pos) {
      if ((pos == 0 || permutation[pos].hash != permutation[pos - 1].hash) && permutation[pos].valid ) {
        size_t num_local_nets = 0, num_local_pins = 0;
        size_t hash = permutation[pos].hash;

        for ( ; pos < permutation.size() && hash == permutation[pos].hash; ++pos) {
          const auto& rep = permutation[pos];
          HyperedgeWeight rep_weight = edgeWeight(rep.he);
          if (rep.valid) {
            auto& contained = local_maps.local();
            for (HypernodeID v : coarse_pin_lists[rep.he]) {
              contained.set(v);
            }

            for (size_t j = pos + 1; j < permutation.size() && hash == permutation[j].hash; ++j) {
              auto& cand = permutation[j];
              const auto& cand_pins = coarse_pin_lists[cand.he];
              if (cand.valid && cand_pins.size() == rep.size &&
                  std::all_of(cand_pins.begin(), cand_pins.end(), [&](HypernodeID v) { return contained[v]; })) {
                cand.valid = false;
                rep_weight += edgeWeight(cand.he);
                coarse_pin_lists[cand.he].clear();    // globally mark net as removed
              }
            }

            coarse_edge_weights[rep.he] = rep_weight;
            num_local_nets++;
            num_local_pins += coarse_pin_lists[rep.he].size();

            for (HypernodeID v : coarse_pin_lists[rep.he]) {
              contained.reset(v);
            }
          }

        }

        __atomic_fetch_add(&num_coarse_nets, num_local_nets, __ATOMIC_RELAXED);
        __atomic_fetch_add(&num_coarse_pins, num_local_pins, __ATOMIC_RELAXED);
      }
    });

    timer.stop_timer("detect duplicates");

    /*
    // identical net detection
    tbb::parallel_for(0UL, net_map.numBuckets(), [&](const size_t bucket_id) {
      size_t num_local_nets = 0, num_local_pins = 0;
      auto& bucket = net_map.getBucket(bucket_id);
      std::sort(bucket.begin(), bucket.end());
      for (size_t i = 0; i < bucket.size(); ++i) {
        const auto& rep = bucket[i];
        HyperedgeWeight rep_weight = edgeWeight(rep.he);
        if (rep.valid) {
          auto& contained = local_maps.local();
          for (HypernodeID v : coarse_pin_lists[rep.he]) { contained.set(v); }

          for (size_t j = i+1; j < bucket.size(); ++j) {
            auto& cand = bucket[j];
            if (cand.hash != rep.hash) { break; }
            const auto& cand_pins = coarse_pin_lists[cand.he];
            if (cand.valid && coarse_pin_lists[rep.he].size() == cand_pins.size()
                && std::all_of(cand_pins.begin(), cand_pins.end(), [&](HypernodeID v) { return contained[v];})) {
              cand.valid = false;
              rep_weight += edgeWeight(cand.he);
              coarse_pin_lists[cand.he].clear();    // globally mark net as removed
            }
          }
          coarse_edge_weights[rep.he] = rep_weight;
          num_local_nets++;
          num_local_pins += coarse_pin_lists[rep.he].size();
          for (HypernodeID v : coarse_pin_lists[rep.he]) { contained.reset(v); }
        }
      }
      net_map.free(bucket_id);
      __atomic_fetch_add(&num_coarse_nets, num_local_nets, __ATOMIC_RELAXED);
      __atomic_fetch_add(&num_coarse_pins, num_local_pins, __ATOMIC_RELAXED);
    });

    */

    timer.stop_timer("identical net detection");
    timer.start_timer("allocs","allocs");

    StaticHypergraph chg;
    chg._num_hypernodes = num_coarse_nodes;
    chg._num_hyperedges = num_coarse_nets;
    chg._num_pins = num_coarse_pins;
    chg._total_degree = num_coarse_pins;
    chg._total_weight = _total_weight;   // didn't lose any vertices. or at least we don't want to let the imbalance calculation know about it...

    tbb::parallel_invoke([&] {
      chg._incident_nets.resize(num_coarse_pins);
    }, [&] {
      chg._incidence_array.resize(num_coarse_pins);
    }, [&]{
      chg._community_ids.resize(num_coarse_nodes);
    }, [&] {
      chg._hyperedges.resize(num_coarse_nets);
    }, [&] {
      chg._hypernodes.resize(num_coarse_nodes);
    });

    timer.stop_timer("allocs");
    timer.start_timer("write pin lists", "write pin lists and count degrees");

    auto& offsets_for_fine_nets = _tmp_contraction_buffer->offsets_for_fine_nets;

    auto net_size_prefix_sum = [&](const tbb::blocked_range<HyperedgeID>& r,
            std::pair<size_t, HyperedgeID> sums, bool is_final_scan) -> std::pair<size_t,HyperedgeID> {
      size_t net_size_sum = sums.first;
      HyperedgeID coarse_net_id = sums.second;
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {
        if (!coarse_pin_lists[he].empty()) {
          if (is_final_scan) {
            chg._hyperedges[coarse_net_id].enable();
            chg._hyperedges[coarse_net_id].setSize(coarse_pin_lists[he].size());
            chg._hyperedges[coarse_net_id].setFirstEntry(net_size_sum);
            chg._hyperedges[coarse_net_id].setWeight(coarse_edge_weights[he]);
            offsets_for_fine_nets[he] = net_size_sum;
          }
          net_size_sum += coarse_pin_lists[he].size();
          coarse_net_id++;
        }
      }
      return std::make_pair(net_size_sum, coarse_net_id);
    };
    auto sum_pair = [](std::pair<size_t, HyperedgeID> l, std::pair<size_t, HyperedgeID> r) {
      return std::make_pair(l.first + r.first, l.second + r.second);
    };
    tbb::parallel_scan(tbb::blocked_range<HyperedgeID>(0U, initialNumEdges()), std::make_pair(0UL,0U),
            net_size_prefix_sum, sum_pair);

    doParallelForAllEdges([&](HyperedgeID he) {
      // removed nets are marked via empty pin list
      if (!coarse_pin_lists[he].empty()) {
        size_t pos = offsets_for_fine_nets[he];
        for (HypernodeID v : coarse_pin_lists[he]) {
          chg._incidence_array[pos++] = v;                                        // copy pin list
          __atomic_fetch_add(&chg._hypernodes[v]._size, 1, __ATOMIC_RELAXED);     // increment pin's degree
        }
      }
    });

    timer.stop_timer("write pin lists");
    timer.start_timer("write incident nets", "write incident nets");

    auto degree_prefix_sum = [&](const tbb::blocked_range<HypernodeID>& r, size_t sum, bool is_final_scan) -> size_t {
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if (is_final_scan) {
          chg._hypernodes[u].enable();
          chg._hypernodes[u].setFirstEntry(sum);
        }
        sum += chg._hypernodes[u]._size;
      }
      return sum;
    };
    tbb::parallel_scan(tbb::blocked_range<HypernodeID>(0U, num_coarse_nodes), 0UL, degree_prefix_sum, std::plus<>());

    tbb::parallel_for(0U, num_coarse_nets, [&](HyperedgeID he) {
      // pin lists fully constructed --> safe to use
      for (HypernodeID v : chg.pins(he)) {
        const size_t pos = __atomic_fetch_add(&chg._hypernodes[v]._begin, 1, __ATOMIC_RELAXED);
        chg._incident_nets[pos] = he;
      }
    });

    // reset begin pointers of nodes that we destroyed when writing the incident nets, and make the order deterministic
    tbb::parallel_for(0U, num_coarse_nodes, [&](HypernodeID u) {
      chg._hypernodes[u]._weight = 0;
      chg._hypernodes[u]._begin -= chg._hypernodes[u].size();
      if (chg._hypernodes[u].size() > 100000) {
        tbb::parallel_sort(chg._incident_nets.begin() + chg._hypernodes[u].firstEntry(),
                           chg._incident_nets.begin() + chg._hypernodes[u].firstInvalidEntry());
      } else {
        std::sort(chg._incident_nets.begin() + chg._hypernodes[u].firstEntry(),
                  chg._incident_nets.begin() + chg._hypernodes[u].firstInvalidEntry());
      }
    });

    timer.stop_timer("write incident nets");
    timer.start_timer("find max edge size", "find max edge size");

    auto find_max_net_size = [&](const tbb::blocked_range<HyperedgeID>& r, HypernodeID max_net_size) {
      for (HyperedgeID e = r.begin(); e < r.end(); ++e) { max_net_size = std::max(max_net_size, chg.edgeSize(e)); }
      return max_net_size;
    };
    auto get_max = [&](HypernodeID lhs, HypernodeID rhs) { return std::max(lhs, rhs); };
    chg._max_edge_size = tbb::parallel_reduce(tbb::blocked_range<HyperedgeID>(0U, num_coarse_nets),
                                      0,find_max_net_size, get_max);

    timer.stop_timer("find max edge size");
    timer.start_timer("aggregate node weights", "aggregate node weights");

    doParallelForAllNodes([&](HypernodeID u) {
      __atomic_fetch_add(&chg._hypernodes[get_cluster(u)]._weight, nodeWeight(u), __ATOMIC_RELAXED);
      chg.setCommunityID(get_cluster(u), communityID(u));
    });

    timer.stop_timer("aggregate node weights");
    timer.stop_timer("hypergraph_contraction");

    chg._tmp_contraction_buffer = _tmp_contraction_buffer;
    _tmp_contraction_buffer = nullptr;
    return chg;
  }



  // ! Copy static hypergraph in parallel
  StaticHypergraph StaticHypergraph::copy(parallel_tag_t) {
    StaticHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    tbb::parallel_invoke([&] {
      hypergraph._hypernodes.resize(_hypernodes.size());
      memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
             sizeof(Hypernode) * _hypernodes.size());
    }, [&] {
      hypergraph._incident_nets.resize(_incident_nets.size());
      memcpy(hypergraph._incident_nets.data(), _incident_nets.data(),
             sizeof(HyperedgeID) * _incident_nets.size());
    }, [&] {
      hypergraph._hyperedges.resize(_hyperedges.size());
      memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
             sizeof(Hyperedge) * _hyperedges.size());
    }, [&] {
      hypergraph._incidence_array.resize(_incidence_array.size());
      memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
             sizeof(HypernodeID) * _incidence_array.size());
    }, [&] {
      hypergraph._community_ids = _community_ids;
    });
    return hypergraph;
  }

  // ! Copy static hypergraph sequential
  StaticHypergraph StaticHypergraph::copy() {
    StaticHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    hypergraph._hypernodes.resize(_hypernodes.size());
    memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
           sizeof(Hypernode) * _hypernodes.size());
    hypergraph._incident_nets.resize(_incident_nets.size());
    memcpy(hypergraph._incident_nets.data(), _incident_nets.data(),
           sizeof(HyperedgeID) * _incident_nets.size());

    hypergraph._hyperedges.resize(_hyperedges.size());
    memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
           sizeof(Hyperedge) * _hyperedges.size());
    hypergraph._incidence_array.resize(_incidence_array.size());
    memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
           sizeof(HypernodeID) * _incidence_array.size());

    hypergraph._community_ids = _community_ids;

    return hypergraph;
  }




  void StaticHypergraph::memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
    parent->addChild("Incident Nets", sizeof(HyperedgeID) * _incident_nets.size());
    parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
    parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
    parent->addChild("Communities", sizeof(PartitionID) * _community_ids.capacity());
  }

  // ! Computes the total node weight of the hypergraph
  void StaticHypergraph::computeAndSetTotalNodeWeight(parallel_tag_t) {
    _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_hypernodes), 0,
                                         [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
                                           HypernodeWeight weight = init;
                                           for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
                                             if (nodeIsEnabled(hn)) {
                                               weight += this->_hypernodes[hn].weight();
                                             }
                                           }
                                           return weight;
                                         }, std::plus<>());
  }

} // namespace
