/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/datastructures/incident_edge_array.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

namespace mt_kahypar {
namespace ds {

void print_header(const IncidentEdgeArray::Header& header) {
  LOG << V(header.prev) << ", " << V(header.next) << ", " << V(header.it_prev) << ", " << V(header.it_next);
  LOG << V(header.tail) << ", " << V(header.first_active) << ", " << V(header.first_inactive) << V(header.degree);
  LOG << V(header.current_version) << V(header.is_head);
}

IncidentEdgeIterator::IncidentEdgeIterator(const HypernodeID u,
                                           const IncidentEdgeArray* incident_edge_array,
                                           const size_t pos,
                                           const bool end):
    _u(u),
    _current_u(u),
    _current_size(incident_edge_array->header(u)->size()),
    _current_pos(pos),
    _incident_edge_array(incident_edge_array),
    _end(end) {
  if ( end ) {
    _current_pos = _current_size;
  }

  ASSERT(pos <= incident_edge_array->nodeDegree(u));
  traverse_headers();
}

HyperedgeID IncidentEdgeIterator::operator* () const {
  return _incident_edge_array->firstActiveEdge(_current_u) + _current_pos;
}

IncidentEdgeIterator & IncidentEdgeIterator::operator++ () {
  ASSERT(!_end);
  ++_current_pos;
  traverse_headers();
  return *this;
}

bool IncidentEdgeIterator::operator!= (const IncidentEdgeIterator& rhs) {
  return !(*this == rhs);
}

bool IncidentEdgeIterator::operator== (const IncidentEdgeIterator& rhs) {
  return _u == rhs._u && _end == rhs._end;
}

void IncidentEdgeIterator::traverse_headers() {
  while ( _current_pos >= _current_size ) {
    const HypernodeID last_u = _current_u;
    _current_u = _incident_edge_array->header(last_u)->it_next;
    _current_pos -= _current_size;
    _current_size = _incident_edge_array->header(_current_u)->size();
    // It can happen that due to a contraction the current vertex
    // we iterate over becomes empty or the head of the current vertex
    // changes. Therefore, we set the end flag if we reach the current
    // head of the list or it_next is equal with the current vertex (means
    // that list becomes empty due to a contraction)
    if ( _incident_edge_array->header(_current_u)->is_head ||
         last_u == _current_u ) {
      _end = true;
      break;
    }
  }
}

void IncidentEdgeArray::construct(const EdgeVector& edge_vector, const HyperedgeWeight* edge_weight) {
  // Accumulate degree of each vertex thread local
  const HyperedgeID num_hyperedges = edge_vector.size();
  ThreadLocalCounter local_incident_nets_per_vertex(_num_nodes + 1, 0);
  AtomicCounter current_incident_net_pos;
  tbb::parallel_invoke([&] {
    tbb::parallel_for(ID(0), num_hyperedges, [&](const size_t pos) {
      parallel::scalable_vector<size_t>& num_incident_nets_per_vertex =
        local_incident_nets_per_vertex.local();
        ++num_incident_nets_per_vertex[edge_vector[pos].first + 1];
        ++num_incident_nets_per_vertex[edge_vector[pos].second + 1];
    });
  }, [&] {
    _index_array.assign(_num_nodes + 1, sizeof(Header) / sizeof(Edge));
    current_incident_net_pos.assign(
      _num_nodes, parallel::IntegralAtomicWrapper<size_t>(0));
  });

  // We sum up the number of incident nets per vertex only thread local.
  // To obtain the global number of incident nets per vertex, we iterate
  // over each thread local counter and sum it up.
  for ( const parallel::scalable_vector<size_t>& c : local_incident_nets_per_vertex ) {
    tbb::parallel_for(ID(0), _num_nodes + 1, [&](const size_t pos) {
      _index_array[pos] += c[pos];
    });
  }

  // Compute start positon of the incident nets of each vertex via a parallel prefix sum
  parallel::TBBPrefixSum<HyperedgeID, Array> incident_net_prefix_sum(_index_array);
  tbb::parallel_scan(tbb::blocked_range<size_t>(
          0UL, UI64(_num_nodes + 1)), incident_net_prefix_sum);
  _size_in_bytes = incident_net_prefix_sum.total_sum() * sizeof(Edge);
  _incident_edge_array = parallel::make_unique<Edge>(_size_in_bytes / sizeof(Edge));

  // Insert incident nets into incidence array
  tbb::parallel_for(ID(0), num_hyperedges, [&](const HyperedgeID he) {
    HypernodeID source = edge_vector[he].first;
    HypernodeID target = edge_vector[he].second;
    const HyperedgeWeight weight = edge_weight == nullptr ? 1 : edge_weight[he];
    Edge& e1 = edge(firstEdge(source) + current_incident_net_pos[source].fetch_add(1));
    e1.source = source;
    e1.target = target;
    e1.weight = weight;
    e1.version = 0;
    Edge& e2 = edge(firstEdge(target) + current_incident_net_pos[target].fetch_add(1));
    e2.source = target;
    e2.target = source;
    e2.weight = weight;
    e2.version = 0;
  });

  // Setup Header of each vertex
  tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
    Header* head = header(u);
    head->prev = u;
    head->next = u;
    head->it_prev = u;
    head->it_next = u;
    head->degree = current_incident_net_pos[u].load(std::memory_order_relaxed);
    head->first_active = 0;
    head->first_inactive = head->degree;
    head->current_version = 0;
    head->is_head = true;
  });
}

void IncidentEdgeArray::contract(const HypernodeID u,
                                 const HypernodeID v,
                                 kahypar::ds::FastResetFlagArray<>& node_bitset,
                                 const AcquireLockFunc& acquire_lock,
                                 const ReleaseLockFunc& release_lock) {
  node_bitset.reset();

  // iterate over edges of u, add neighbors to the bitset and remove the contracted edge (if present)
  Header* head_u = header(u);
  for (HypernodeID current_u: headers(u)) {
    Header* head = header(current_u);
    const HypernodeID new_version = ++head->current_version;
    for ( HyperedgeID current_edge = firstActiveEdge(current_u); current_edge < firstInactiveEdge(current_u); ) {
      Edge& e = edge(current_edge);
      if ( e.target == v ) {
        swap_to_back(current_u, current_edge);
        --head_u->degree;
      } else {
        e.version = new_version;
        // add neighbor to bitset
        ASSERT(!node_bitset[e.target], V(u) << "contains duplicate edge to node" << e.target);
        node_bitset.set(e.target, true);
        ++current_edge;
      }
    }

    if ( head->size() == 0 && current_u != u ) {
      removeEmptyIncidentEdgeList(current_u);
    }
  }
  // ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");

  // iterate over edges of v and remove duplicate edges
  Header* head_v = header(v);
  for (HypernodeID current_v: headers(v)) {
    Header* head = header(current_v);
    const HypernodeID new_version = ++head->current_version;
    for ( HyperedgeID current_edge = firstActiveEdge(current_v); current_edge < firstInactiveEdge(current_v); ) {
      Edge& e = edge(current_edge);
      if ( e.target == u ) {
        // contracted edge
        swap_to_back(current_v, current_edge);
        --head_v->degree;
      } else if (node_bitset[e.target]) {
        // duplicate edge
        // TODO(maas): locking?!
        HyperedgeID backwardsEdge = findOutgoingEdge(e.target, v);
        swap_to_front(edge(backwardsEdge).source, backwardsEdge);
        --header(e.target)->degree;

        swap_to_back(current_v, current_edge);
        --head_v->degree;
      } else {
        e.version = new_version;
        // TODO(maas): locking?!
        HyperedgeID backwardsEdge = findOutgoingEdge(e.target, v);
        edge(backwardsEdge).target = u;
        ++current_edge;
      }
    }

    if ( head->size() == 0 && current_v != v ) {
      removeEmptyIncidentEdgeList(current_v);
    }
  }

  acquire_lock(u);
  // Concatenate double-linked list of u and v
  append(u, v);
  header(u)->degree += header(v)->degree;
  // ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");
  release_lock(u);
}

void IncidentEdgeArray::uncontract(const HypernodeID u,
                  const HypernodeID v,
                  const AcquireLockFunc& acquire_lock,
                  const ReleaseLockFunc& release_lock) {
  ASSERT(header(v)->prev != v);
  Header* head_u = header(u);
  Header* head_v = header(v);
  acquire_lock(u);
  // Restores the incident list of v to the time before it was appended
  // to the double-linked list of u.
  splice(u, v);
  header(u)->degree -= head_v->degree;
  // ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");
  release_lock(u);

  // iterate over linked list of u to restore the contracted edge (if present)
  HypernodeID last_non_empty_u = u;
  for (HypernodeID current_u: headers(u)) {
    Header* head = header(current_u);
    const HypernodeID new_version = --head->current_version;
    for (HyperedgeID current_edge = firstActiveEdge(current_u); current_edge < current_u; ++current_edge) {
      Edge& e = edge(current_edge);
      e.version = new_version;
    }
    Edge& e = edge(firstInactiveEdge(current_u));
    if (e.version == new_version) {
      ASSERT(e.version == head->current_version);
      e.target = v;
      ++head->first_inactive;
      ++head_u->degree;
    }
    if (head->size() > 0) {
      restoreItLink(u, last_non_empty_u, current_u);
      last_non_empty_u = current_u;
    }
  }
  // ASSERT(verifyIteratorPointers(u), "Iterator pointers of vertex" << u << "are corrupted");

  // iterate over edges of v, update backwards edges and restore removed edges
  HypernodeID last_non_empty_v = v;
  for (HypernodeID current_v: headers(v)) {
    Header* head = header(current_v);
    ASSERT(head->current_version > 0);
    const HypernodeID new_version = --head->current_version;
    const HyperedgeID first_inactive = firstInactiveEdge(current_v);
    for (HyperedgeID current_edge = firstActiveEdge(current_v); current_edge < first_inactive; ++current_edge) {
      Edge& e = edge(current_edge);
      e.version = new_version;
      // TODO(maas): locking?!
      HyperedgeID backwardsEdge = findOutgoingEdge(e.target, u);
      edge(backwardsEdge).target = v;
    }

    const HyperedgeID last_entry = _index_array[current_v + 1];
    for (HyperedgeID current_edge = first_inactive; current_edge < last_entry; ++current_edge) {
      Edge& e = edge(current_edge);
      if (e.version != new_version) {
        break;
      }
      ++head->first_inactive;
      ++head_v->degree;

      if (e.target != u) {
        const HyperedgeID backwardsEdge = findRemovedOutgoingEdge(e.target, v);
        const HypernodeID neighbor = edge(backwardsEdge).source;
        swap(edge(backwardsEdge), edge(firstActiveEdge(neighbor) - 1));
        --header(neighbor)->first_active;
        ++header(e.target)->degree;
      }
    }

    if (head->size() > 0) {
      restoreItLink(v, last_non_empty_v, current_v);
      last_non_empty_v = current_v;
    }
  }
}

HyperedgeID IncidentEdgeArray::findOutgoingEdge(const HypernodeID u, const HypernodeID target) const {
  for (HyperedgeID e: incidentEdges(u)) {
    if (edge(e).target == target) {
      return e;
    }
  }
  ASSERT(false, "Hypernode" << u << "has no outgoing edge with target" << target);
}

HyperedgeID IncidentEdgeArray::findRemovedOutgoingEdge(const HypernodeID u, const HypernodeID target) const {
  for (HypernodeID current_u: headers(u)) {
    const HyperedgeID num_edges = firstActiveEdge(current_u) - firstEdge(current_u);
    // we consider the edges in reverse order because
    // the searched edge is more likely at the end of the range
    for (HyperedgeID offset = 1; offset <= num_edges; ++offset) {
      HyperedgeID e = firstActiveEdge(current_u) - offset;
      if (edge(e).target == target) {
        return e;
      }
    }
  }
  ASSERT(false, "Hypernode" << u << "has no outgoing edge with target" << target);
}

void IncidentEdgeArray::append(const HypernodeID u, const HypernodeID v) {
  const HypernodeID tail_u = header(u)->prev;
  const HypernodeID tail_v = header(v)->prev;
  header(tail_u)->next = v;
  header(v)->prev = tail_u;
  header(tail_v)->next = u;
  header(u)->prev = tail_v;

  const HypernodeID it_tail_u = header(u)->it_prev;
  const HypernodeID it_tail_v = header(v)->it_prev;
  header(it_tail_u)->it_next = v;
  header(v)->it_prev = it_tail_u;
  header(it_tail_v)->it_next = u;
  header(u)->it_prev = it_tail_v;

  header(v)->tail = tail_v;
  header(v)->is_head = false;

  if ( header(v)->size() == 0 ) {
    removeEmptyIncidentEdgeList(v);
  }
}

void IncidentEdgeArray::splice(const HypernodeID u, const HypernodeID v) {
  // Restore the iterator double-linked list of u such that it does not contain
  // any incident net list of v
  const HypernodeID tail = header(v)->tail;
  HypernodeID non_empty_entry_prev_v = v;
  HypernodeID non_empty_entry_next_tail = tail;
  while ( ( non_empty_entry_prev_v == v ||
          header(non_empty_entry_prev_v)->size() == 0 ) &&
          non_empty_entry_prev_v != u ) {
    non_empty_entry_prev_v = header(non_empty_entry_prev_v)->prev;
  }
  while ( ( non_empty_entry_next_tail == tail ||
          header(non_empty_entry_next_tail)->size() == 0 ) &&
          non_empty_entry_next_tail != u ) {
    non_empty_entry_next_tail = header(non_empty_entry_next_tail)->next;
  }
  header(non_empty_entry_prev_v)->it_next = non_empty_entry_next_tail;
  header(non_empty_entry_next_tail)->it_prev = non_empty_entry_prev_v;

  // Cut out incident list of v
  const HypernodeID prev_v = header(v)->prev;
  const HypernodeID next_tail = header(tail)->next;
  header(v)->prev = tail;
  header(tail)->next = v;
  header(next_tail)->prev = prev_v;
  header(prev_v)->next = next_tail;
  header(v)->is_head = true;
}

void IncidentEdgeArray::removeEmptyIncidentEdgeList(const HypernodeID u) {
  ASSERT(!header(u)->is_head);
  ASSERT(header(u)->size() == 0, V(u) << V(header(u)->size()));
  Header* head = header(u);
  header(head->it_prev)->it_next = head->it_next;
  header(head->it_next)->it_prev = head->it_prev;
  head->it_next = u;
  head->it_prev = u;
}

void IncidentEdgeArray::restoreItLink(const HypernodeID u, const HypernodeID prev, const HypernodeID current) {
  header(prev)->it_next = current;
  header(current)->it_prev = prev;
  header(current)->it_next = u;
  header(u)->it_prev = current;
}

}  // namespace ds
}  // namespace mt_kahypar
