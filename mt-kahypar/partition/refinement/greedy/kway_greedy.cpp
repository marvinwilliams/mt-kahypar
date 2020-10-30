/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
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

#include "mt-kahypar/partition/refinement/greedy/kway_greedy.h"

namespace mt_kahypar {

bool KWayGreedy::findMoves(PartitionedHypergraph &phg, size_t taskID,
                           size_t numSeeds) {
  localMoves.clear();
  thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

  HypernodeID seedNode;
  while (runStats.pushes < numSeeds &&
         sharedData.refinementNodes.try_pop(seedNode, taskID)) {
    SearchID previousSearchOfSeedNode =
        sharedData.nodeTracker.searchOfNode[seedNode].load(
            std::memory_order_relaxed);
    if (sharedData.nodeTracker.tryAcquireNode(seedNode, thisSearch)) {
      fm_strategy.insertIntoPQ(phg, seedNode, previousSearchOfSeedNode);
    }
  }
  fm_strategy.updatePQs(phg);

  if (runStats.pushes > 0) {
    internalFindMoves(phg);
    return true;
  } else {
    return false;
  }
}

template <typename Partition>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::pair<PartitionID, HypernodeWeight>
heaviestPartAndWeight(const Partition &partition) {
  PartitionID p = kInvalidPartition;
  HypernodeWeight w = std::numeric_limits<HypernodeWeight>::min();
  for (PartitionID i = 0; i < partition.k(); ++i) {
    if (partition.partWeight(i) > w) {
      w = partition.partWeight(i);
      p = i;
    }
  }
  return std::make_pair(p, w);
}

template <typename PHG>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void
KWayGreedy::updateNeighbors(PHG &phg, const Move &move) {
  // Note: In theory we should acquire/update all neighbors. It just turned out
  // that this works fine Actually: only vertices incident to edges with gain
  // changes can become new boundary vertices. Vertices that already were
  // boundary vertices, can still be considered later since they are in the task
  // queue
  // --> actually not that bad
  for (HyperedgeID e : edgesWithGainChanges) {
    if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
      for (HypernodeID v : phg.pins(e)) {
        if (neighborDeduplicator[v] != deduplicationTime) {
          SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(
              std::memory_order_acq_rel);
          if (searchOfV == thisSearch) {
            fm_strategy.updateGain(phg, v, move);
          }
          /* TODO: dont aquire neighbor, something else to do? Inform other
           * thread about update? <26-10-20, @noahares> */
          //            else if (sharedData.nodeTracker.tryAcquireNode(v,
          //            thisSearch)) {
          //              fm_strategy.insertIntoPQ(phg, v, searchOfV);
          //            }
          neighborDeduplicator[v] = deduplicationTime;
        }
      }
    }
  }
  edgesWithGainChanges.clear();

  if (++deduplicationTime == 0) {
    neighborDeduplicator.assign(neighborDeduplicator.size(), 0);
    deduplicationTime = 1;
  }
}

void KWayGreedy::internalFindMoves(PartitionedHypergraph &phg) {
  Move move;

  auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight edge_weight,
                        const HypernodeID,
                        const HypernodeID pin_count_in_from_part_after,
                        const HypernodeID pin_count_in_to_part_after) {
    // Gains of the pins of a hyperedge can only change in the following
    // situations.
    if (pin_count_in_from_part_after == 0 ||
        pin_count_in_from_part_after == 1 || pin_count_in_to_part_after == 1 ||
        pin_count_in_to_part_after == 2) {
      edgesWithGainChanges.push_back(he);
    }

    fm_strategy.deltaGainUpdates(phg, he, edge_weight, move.from,
                                 pin_count_in_from_part_after, move.to,
                                 pin_count_in_to_part_after);
  };

  // we can almost make this function take a generic partitioned hypergraph
  // we would have to add the success func to the interface of DeltaPhg (and
  // then ignore it there...) and do the local rollback outside this function

  size_t bestImprovementIndex = 0;
  Gain estimatedImprovement = 0;
  Gain bestImprovement = 0;
  Gain lastImprovement = std::numeric_limits<Gain>::max();

  HypernodeWeight heaviestPartWeight = 0;
  HypernodeWeight fromWeight = 0, toWeight = 0;

  /* TODO: is this enough? Any more checks to prevent the negative move?
   * <27-10-20, @noahares> */
  while (lastImprovement > 0 &&
         sharedData.finishedTasks.load(std::memory_order_relaxed) <
             sharedData.finishedTasksLimit) {

    if (!fm_strategy.findNextMoveNoRetry(phg, move))
      break;

    sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
    MoveID move_id = std::numeric_limits<MoveID>::max();
    bool moved = false;
    if (move.to != kInvalidPartition) {
      heaviestPartWeight = heaviestPartAndWeight(phg).second;
      fromWeight = phg.partWeight(move.from);
      toWeight = phg.partWeight(move.to);
      moved = phg.changeNodePart(
          move.node, move.from, move.to,
          context.partition.max_part_weights[move.to],
          [&] { move_id = sharedData.moveTracker.insertMove(move); },
          delta_func);
    }

    if (moved) {
      runStats.moves++;
      estimatedImprovement += move.gain;
      localMoves.emplace_back(move, move_id);
      lastImprovement = move.gain;
      const bool improved_km1 = estimatedImprovement > bestImprovement;
      const bool improved_balance_less_equal_km1 =
          estimatedImprovement >= bestImprovement &&
          fromWeight == heaviestPartWeight &&
          toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

      if (improved_km1 || improved_balance_less_equal_km1) {
        //          stopRule.reset();
        bestImprovement = estimatedImprovement;
        bestImprovementIndex = localMoves.size();
      }

      updateNeighbors(phg, move);
    }

    fm_strategy.updatePQs(phg);
  }

  runStats.estimated_improvement = bestImprovement;
  fm_strategy.clearPQs(bestImprovementIndex);
  runStats.merge(stats);
}

void KWayGreedy::memoryConsumption(utils::MemoryTreeNode *parent) const {
  ASSERT(parent);

  utils::MemoryTreeNode *localized_fm_node = parent->addChild("k-Way Greedy");

  utils::MemoryTreeNode *deduplicator_node =
      localized_fm_node->addChild("Deduplicator");
  deduplicator_node->updateSize(neighborDeduplicator.capacity() *
                                sizeof(HypernodeID));
  utils::MemoryTreeNode *edges_to_activate_node =
      localized_fm_node->addChild("edgesWithGainChanges");
  edges_to_activate_node->updateSize(edgesWithGainChanges.capacity() *
                                     sizeof(HyperedgeID));

  utils::MemoryTreeNode *local_moves_node = parent->addChild("Local FM Moves");
  local_moves_node->updateSize(localMoves.capacity() *
                               sizeof(std::pair<Move, MoveID>));

  fm_strategy.memoryConsumption(localized_fm_node);
  // TODO fm_strategy.memoryConsumptiom(..)
  /*
  utils::MemoryTreeNode* block_pq_node = localized_fm_node->addChild("Block
  PQ"); block_pq_node->updateSize(blockPQ.size_in_bytes());
  utils::MemoryTreeNode* vertex_pq_node = localized_fm_node->addChild("Vertex
  PQ"); for ( const VertexPriorityQueue& pq : vertexPQs ) {
    vertex_pq_node->updateSize(pq.size_in_bytes());
  }
   */
}

} // namespace mt_kahypar
