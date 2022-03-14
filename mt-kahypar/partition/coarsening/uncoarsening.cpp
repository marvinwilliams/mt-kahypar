#include <mt-kahypar/partition/refinement/do_nothing_refiner.h>
#include <mt-kahypar/partition/refinement/label_propagation/async_lp_refiner.h>
#include "mt-kahypar/partition/coarsening/multilevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener_base.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/partitioning_output.h"

#include "mt-kahypar/datastructures/async/array_lock_manager.h"
#include "mt-kahypar/datastructures/async/group_pool.h"

namespace mt_kahypar {

  void MultilevelCoarsenerBase::finalize() {
    utils::Timer::instance().start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
    // Free memory of temporary contraction buffer and
    // release coarsening memory in memory pool
    currentHypergraph().freeTmpContractionBuffer();
    if (_top_level) {
      parallel::MemoryPool::instance().release_mem_group("Coarsening");
    }

    // Construct top level partitioned hypergraph (memory is taken from memory pool)
    _partitioned_hg = PartitionedHypergraph(
            _context.partition.k, _hg, parallel_tag_t());

    // Construct partitioned hypergraphs parallel
    tbb::task_group group;
    // Construct partitioned hypergraph for each coarsened hypergraph in the hierarchy
    for (size_t i = 0; i < _hierarchy.size(); ++i) {
      group.run([&, i] {
        _hierarchy[i].contractedPartitionedHypergraph() = PartitionedHypergraph(
                _context.partition.k, _hierarchy[i].contractedHypergraph(), parallel_tag_t());
      });
    }
    group.wait();

    // Set the representative partitioned hypergraph for each hypergraph
    // in the hierarchy
    if (_hierarchy.size() > 0) {
      _hierarchy[0].setRepresentativeHypergraph(&_partitioned_hg);
      for (size_t i = 1; i < _hierarchy.size(); ++i) {
        _hierarchy[i].setRepresentativeHypergraph(&_hierarchy[i - 1].contractedPartitionedHypergraph());
      }
    }
    _is_finalized = true;
    utils::Timer::instance().stop_timer("finalize_multilevel_hierarchy");
  }

  void NLevelCoarsenerBase::finalize() {
    // Create compactified hypergraph containing only enabled vertices and hyperedges
    // with consecutive IDs => Less complexity in initial partitioning.
    utils::Timer::instance().start_timer("compactify_hypergraph", "Compactify Hypergraph");
    auto compactification = HypergraphFactory::compactify(_hg);
    _compactified_hg = std::move(compactification.first);
    _compactified_hn_mapping = std::move(compactification.second);
    _compactified_phg = PartitionedHypergraph(_context.partition.k, _compactified_hg, parallel_tag_t());
    utils::Timer::instance().stop_timer("compactify_hypergraph");

    if (_context.uncoarsening.use_asynchronous_uncoarsening) {
        utils::Timer::instance().start_timer("create_uncontraction_pools", "Create Uncontraction Group Pools");
        _group_pools_for_versions = _hg.createUncontractionGroupPoolsForVersions(_context);
//        HEAVY_COARSENING_ASSERT(_hg.verifyIncidenceArraySortedness(_group_pools_for_versions));
        auto num_nodes = _hg.initialNumNodes();
        _lock_manager_for_async = std::make_unique<ds::ArrayLockManager<HypernodeID, ds::ContractionGroupID>>(num_nodes, ds::invalidGroupID);
        utils::Timer::instance().stop_timer("create_uncontraction_pools");
    } else {
        // Create n-level batch uncontraction hierarchy
        utils::Timer::instance().start_timer("create_batch_uncontraction_hierarchy", "Create n-Level Hierarchy");
        _hierarchy = _hg.createBatchUncontractionHierarchy(_context.refinement.max_batch_size);
        ASSERT(_removed_hyperedges_batches.size() == _hierarchy.size() - 1);
        utils::Timer::instance().stop_timer("create_batch_uncontraction_hierarchy");
    }

    _is_finalized = true;
  }

  namespace {
    double refinementTimeLimit(const Context& context, const double time) {
      if ( context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
        const double time_limit_factor = std::max(1.0,  context.refinement.fm.time_limit_factor * context.partition.k);
        return std::max(5.0, time_limit_factor * time);
      } else {
        return std::numeric_limits<double>::max();
      }
    }
  } // namespace

  void MultilevelCoarsenerBase::performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities,
          const HighResClockTimepoint& round_start) {
    ASSERT(!_is_finalized);
    Hypergraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities);
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), elapsed_time);
  }

  PartitionedHypergraph&& MultilevelCoarsenerBase::doUncoarsen(
          std::unique_ptr<IRefiner>& label_propagation,
          std::unique_ptr<IRefiner>& fm) {
    PartitionedHypergraph& coarsest_hg = currentPartitionedHypergraph();
    kahypar::Metrics current_metrics = initialize(coarsest_hg);

    if (_top_level) {
      _context.initial_km1 = current_metrics.km1;
    }

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                              _context.partition.objective == kahypar::Objective::km1
                                              ? current_metrics.km1 : current_metrics.cut,
                                              _context.partition.verbose_output &&
                                              _context.partition.enable_progress_bar && !debug);
    uncontraction_progress += coarsest_hg.initialNumNodes();

    // Refine Coarsest Partitioned Hypergraph
    double time_limit = refinementTimeLimit(_context, _hierarchy.back().coarseningTime());
    refine(coarsest_hg, label_propagation, fm, current_metrics, time_limit);

    for (int i = _hierarchy.size() - 1; i >= 0; --i) {
      // Project partition to next level finer hypergraph
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
      PartitionedHypergraph& representative_hg = _hierarchy[i].representativeHypergraph();
      PartitionedHypergraph& contracted_hg = _hierarchy[i].contractedPartitionedHypergraph();
      representative_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = _hierarchy[i].mapToContractedHypergraph(hn);
        const PartitionID block = contracted_hg.partID(coarse_hn);
        ASSERT(block != kInvalidPartition && block < representative_hg.k());
        representative_hg.setOnlyNodePart(hn, block);
      });
      representative_hg.initializePartition();

      ASSERT(metrics::objective(representative_hg, _context.partition.objective) ==
             metrics::objective(contracted_hg, _context.partition.objective),
             V(metrics::objective(representative_hg, _context.partition.objective)) <<
                                                                                    V(metrics::objective(
                                                                                            contracted_hg,
                                                                                            _context.partition.objective)));
      ASSERT(metrics::imbalance(representative_hg, _context) ==
             metrics::imbalance(contracted_hg, _context),
             V(metrics::imbalance(representative_hg, _context)) <<
                                                                V(metrics::imbalance(contracted_hg, _context)));
      utils::Timer::instance().stop_timer("projecting_partition");

      // Refinement
      time_limit = refinementTimeLimit(_context, _hierarchy[i].coarseningTime());
      refine(representative_hg, label_propagation, fm, current_metrics, time_limit);

      // Update Progress Bar
      uncontraction_progress.setObjective(
              current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective));
      uncontraction_progress += representative_hg.initialNumNodes() - contracted_hg.initialNumNodes();
    }

    // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
    if (_top_level && !metrics::isBalanced(_partitioned_hg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
              kahypar::Mode::direct_kway, _context.partition.objective);
      if (_context.partition.verbose_output) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
            << metrics::imbalance(_partitioned_hg, _context) << ")" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(_partitioned_hg, _context);
      }

      if (_context.partition.deterministic) {
        if (_context.partition.verbose_output) {
          LOG << RED << "Skip rebalancing since deterministic mode is activated" << END;
        }
      } else {
        if (_context.partition.verbose_output) {
          LOG << RED << "Start rebalancing!" << END;
        }
        utils::Timer::instance().start_timer("rebalance", "Rebalance");
        if (_context.partition.objective == kahypar::Objective::km1) {
          Km1Rebalancer rebalancer(_partitioned_hg, _context);
          rebalancer.rebalance(current_metrics);
        } else if (_context.partition.objective == kahypar::Objective::cut) {
          CutRebalancer rebalancer(_partitioned_hg, _context);
          rebalancer.rebalance(current_metrics);
        }
        utils::Timer::instance().stop_timer("rebalance");

        const HyperedgeWeight quality_after = current_metrics.getMetric(
                kahypar::Mode::direct_kway, _context.partition.objective);
        if (_context.partition.verbose_output) {
          const HyperedgeWeight quality_delta = quality_after - quality_before;
          if (quality_delta > 0) {
            LOG << RED << "Rebalancer decreased solution quality by" << quality_delta
                << "(Current Imbalance:" << metrics::imbalance(_partitioned_hg, _context) << ")" << END;
          } else {
            LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
                << "(Current Imbalance:" << metrics::imbalance(_partitioned_hg, _context) << ")" << END;
          }
        }
      }

      ASSERT(metrics::objective(_partitioned_hg, _context.partition.objective) ==
             current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
             V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective))
                     << V(metrics::objective(_partitioned_hg, _context.partition.objective)));
    }
    return std::move(_partitioned_hg);
  }

  PartitionedHypergraph&& NLevelCoarsenerBase::doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                                           std::unique_ptr<IRefiner>& fm) {

      // Switch to asynchronous uncoarsening if the option is set
      if (_context.uncoarsening.use_asynchronous_uncoarsening) {
            // Global Label Propagation not needed in asynch case, clear the memory
//            label_propagation.reset();
            return doAsynchronousUncoarsen(fm);
      }


    ASSERT(_is_finalized);
    kahypar::Metrics current_metrics = initialize(_compactified_phg);
    if (_top_level) {
      _context.initial_km1 = current_metrics.km1;
    }

    // Project partition from compactified hypergraph to original hypergraph
    utils::Timer::instance().start_timer("initialize_partition", "Initialize Partition");
    _phg = PartitionedHypergraph(_context.partition.k, _hg, parallel_tag_t());
    _phg.doParallelForAllNodes([&](const HypernodeID hn) {
      ASSERT(static_cast<size_t>(hn) < _compactified_hn_mapping.size());
      const HypernodeID compactified_hn = _compactified_hn_mapping[hn];
      const PartitionID block_id = _compactified_phg.partID(compactified_hn);
      ASSERT(block_id != kInvalidPartition && block_id < _context.partition.k);
      _phg.setOnlyNodePart(hn, block_id);
    });
    _phg.initializePartition();

    if ( _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache ) {
      _phg.initializeGainCache();
    }

    ASSERT(metrics::objective(_compactified_phg, _context.partition.objective) ==
            metrics::objective(_phg, _context.partition.objective),
            V(metrics::objective(_compactified_phg, _context.partition.objective)) <<
            V(metrics::objective(_phg, _context.partition.objective)));
    ASSERT(metrics::imbalance(_compactified_phg, _context) ==
            metrics::imbalance(_phg, _context),
            V(metrics::imbalance(_compactified_phg, _context)) <<
            V(metrics::imbalance(_phg, _context)));
    utils::Timer::instance().stop_timer("initialize_partition");

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
      _context.partition.objective == kahypar::Objective::km1 ? current_metrics.km1 : current_metrics.cut,
      _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug);
    uncontraction_progress += _compactified_hg.initialNumNodes();

    // Initialize Refiner
    if ( label_propagation ) {
      label_propagation->initialize(_phg);
    }
    if ( fm ) {
      fm->initialize(_phg);
    }

    // Perform batch uncontractions
    bool is_timer_disabled = false;
    bool force_measure_timings = _context.partition.measure_detailed_uncontraction_timings && _top_level;
    if ( utils::Timer::instance().isEnabled() ) {
      utils::Timer::instance().disable();
      is_timer_disabled = true;
    }

    ASSERT(_round_coarsening_times.size() == _removed_hyperedges_batches.size());
    _round_coarsening_times.push_back(_round_coarsening_times.size() > 0 ?
      _round_coarsening_times.back() : std::numeric_limits<double>::max()); // Sentinel

    size_t num_batches = 0;
    size_t total_batches_size = 0;
    const size_t minimum_required_number_of_border_vertices = std::max(_context.refinement.max_batch_size,
      _context.shared_memory.num_threads * _context.refinement.min_border_vertices_per_thread);
    ds::StreamingVector<HypernodeID> tmp_refinement_nodes;
    kahypar::ds::FastResetFlagArray<> border_vertices_of_batch(_phg.initialNumNodes());

    // todo mlaupichler remove debug stats
    FMStats localized_fm_stats;
    size_t total_localized_refine_calls = 0;
    size_t total_localized_refine_iterations = 0;

    auto do_localized_refinement = [&]() {
      parallel::scalable_vector<HypernodeID> refinement_nodes = tmp_refinement_nodes.copy_parallel();
      tmp_refinement_nodes.clear_parallel();
      border_vertices_of_batch.reset();
      auto total_fm_stats_before_localized_refine = fm->getTotalFMStats();
      ++total_localized_refine_calls;
        localizedRefine(_phg, refinement_nodes, label_propagation,
                        fm, current_metrics, force_measure_timings, total_localized_refine_iterations);
      auto stats_diff_from_localized_refine = fm->getTotalFMStats();
      stats_diff_from_localized_refine.subtract(total_fm_stats_before_localized_refine);
      stats_diff_from_localized_refine.merge(localized_fm_stats);
    };

    while ( !_hierarchy.empty() ) {
      BatchVector& batches = _hierarchy.back();

      // Uncontract all batches of a specific version of the hypergraph
      while ( !batches.empty() ) {
        const Batch& batch = batches.back();
        if ( batch.size() > 0 ) {
          HEAVY_REFINEMENT_ASSERT(metrics::objective(_phg, _context.partition.objective) ==
                current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                V(metrics::objective(_phg, _context.partition.objective)));
          utils::Timer::instance().start_timer("batch_uncontractions", "Batch Uncontractions", false, force_measure_timings);
          _phg.uncontract(batch);
          utils::Timer::instance().stop_timer("batch_uncontractions", force_measure_timings);
          HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
          HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());
          HEAVY_REFINEMENT_ASSERT(metrics::objective(_phg, _context.partition.objective) ==
                current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                V(metrics::objective(_phg, _context.partition.objective)));

          utils::Timer::instance().start_timer("collect_border_vertices", "Collect Border Vertices", false, force_measure_timings);
          tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
            const Memento& memento = batch[i];
            if ( !border_vertices_of_batch[memento.u] && _phg.isBorderNode(memento.u) ) {
              border_vertices_of_batch.set(memento.u, true);
              tmp_refinement_nodes.stream(memento.u);
            }
            if ( !border_vertices_of_batch[memento.v] && _phg.isBorderNode(memento.v) ) {
              border_vertices_of_batch.set(memento.v, true);
              tmp_refinement_nodes.stream(memento.v);
            }
          });
          utils::Timer::instance().stop_timer("collect_border_vertices", force_measure_timings);

          if ( tmp_refinement_nodes.size() >= minimum_required_number_of_border_vertices ) {
            // Perform localized refinement if we uncontract more
            // than the minimum required number of border vertices
            do_localized_refinement();
          }

          ++num_batches;
          total_batches_size += batch.size();
          // Update Progress Bar
          uncontraction_progress.setObjective(current_metrics.getMetric(
            _context.partition.mode, _context.partition.objective));
          uncontraction_progress += batch.size();
        }
        batches.pop_back();
      }

      if ( tmp_refinement_nodes.size() > 0 ) {
        // Perform localized refinement on remaining border vertices
        do_localized_refinement();
      }

      // Restore single-pin and parallel nets to continue with the next vector of batches
      if ( !_removed_hyperedges_batches.empty() ) {
        utils::Timer::instance().start_timer("restore_single_pin_and_parallel_nets", "Restore Single Pin and Parallel Nets", false, force_measure_timings);
        _phg.restoreSinglePinAndParallelNets(_removed_hyperedges_batches.back());
        _removed_hyperedges_batches.pop_back();
        utils::Timer::instance().stop_timer("restore_single_pin_and_parallel_nets", force_measure_timings);
        HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
        HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());

        // Perform refinement on all vertices
        const double time_limit = refinementTimeLimit(_context, _round_coarsening_times.back());
        globalRefine(_phg, fm, current_metrics, time_limit);
        uncontraction_progress.setObjective(current_metrics.getMetric(
          _context.partition.mode, _context.partition.objective));
        _round_coarsening_times.pop_back();
      }
      _hierarchy.pop_back();
    }

    utils::Stats::instance().update_stat("localized_fm_moves", static_cast<int64_t>(localized_fm_stats.moves));
    utils::Stats::instance().update_stat("localized_fm_local_reverts", static_cast<int64_t>(localized_fm_stats.local_reverts));
    utils::Stats::instance().update_stat("localized_fm_total_find_moves_calls", static_cast<int64_t>(localized_fm_stats.find_moves_calls));
    utils::Stats::instance().update_stat("localized_fm_find_moves_calls_with_good_prefix", static_cast<int64_t>(localized_fm_stats.find_moves_calls_with_good_prefix));
    utils::Stats::instance().update_stat("localized_fm_find_move_retries", static_cast<int64_t>(localized_fm_stats.retries));
    utils::Stats::instance().update_stat("localized_fm_total_pushes", static_cast<int64_t>(localized_fm_stats.pushes_with_pos_gain + localized_fm_stats.pushes_with_non_pos_gain));
    utils::Stats::instance().update_stat("localized_fm_pushes_with_pos_gain", static_cast<int64_t>(localized_fm_stats.pushes_with_pos_gain));
    utils::Stats::instance().update_stat("localized_fm_pushes_with_non_pos_gain", static_cast<int64_t>(localized_fm_stats.pushes_with_non_pos_gain));
    utils::Stats::instance().update_stat("localized_fm_pins_touched_by_delta_gain_cache_updates", static_cast<int64_t>(localized_fm_stats.pins_touched_by_delta_gain_cache_updates));

    size_t sum_gc_cases = localized_fm_stats.num_case_from_zero_gc_updates + localized_fm_stats.num_case_from_one_gc_updates + localized_fm_stats.num_case_to_one_gc_updates + localized_fm_stats.num_case_to_two_gc_updates;
    utils::Stats::instance().update_stat("localized_fm_total_delta_gain_cache_updates", static_cast<int64_t>(sum_gc_cases));
    utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_from_zero", static_cast<int64_t>(localized_fm_stats.num_case_from_zero_gc_updates));
    utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_from_one", static_cast<int64_t>(localized_fm_stats.num_case_from_one_gc_updates));
    utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_to_one", static_cast<int64_t>(localized_fm_stats.num_case_to_one_gc_updates));
    utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_to_two", static_cast<int64_t>(localized_fm_stats.num_case_to_two_gc_updates));

//    //todo mlaupichler remove debug
    if (_context.type == kahypar::ContextType::main) {
      LOG << localized_fm_stats.serialize();
      double frac_local_reverted = (double) localized_fm_stats.local_reverts / (double) localized_fm_stats.moves;
      double frac_pos_gain_pushes = (double) localized_fm_stats.pushes_with_pos_gain / (double) localized_fm_stats.pushes;
      LOG << "Fraction of local reverts: " << frac_local_reverted;
      LOG << "Fraction of pos gain pushes: " << frac_pos_gain_pushes;
      LOG << std::setprecision(5) << std::fixed
                << "FM calls: " << localized_fm_stats.find_moves_calls
                << ", With Good Prefix: " << localized_fm_stats.find_moves_calls_with_good_prefix
                << ", Fraction: " << ((double) localized_fm_stats.find_moves_calls_with_good_prefix / (double) localized_fm_stats.find_moves_calls);
      LOG << "Pins touched in DeltaPHG Gain Cache updates: " << localized_fm_stats.pins_touched_by_delta_gain_cache_updates;
      LOG << std::setprecision(5) << std::fixed << "Number of DeltaPHG Gain Cache updates triggered: "
          << "\n\tTotal: " << sum_gc_cases
          << "\n\tCase: Pins in from=0: " << localized_fm_stats.num_case_from_zero_gc_updates << ", Fraction: " << ((double) localized_fm_stats.num_case_from_zero_gc_updates / (double) sum_gc_cases)
          << "\n\tCase: Pins in from=1: " << localized_fm_stats.num_case_from_one_gc_updates << ", Fraction: " << ((double) localized_fm_stats.num_case_from_one_gc_updates / (double) sum_gc_cases)
          << "\n\tCase: Pins in to=1:   " << localized_fm_stats.num_case_to_one_gc_updates << ", Fraction: " << ((double) localized_fm_stats.num_case_to_one_gc_updates / (double) sum_gc_cases)
          << "\n\tCase: Pins in to=2:   " << localized_fm_stats.num_case_to_two_gc_updates << ", Fraction: " << ((double) localized_fm_stats.num_case_to_two_gc_updates / (double) sum_gc_cases);
      LOG << std::setprecision(5) << std::fixed
          << "localizedRefine() calls: " << total_localized_refine_calls
          << ", lR() Loop Iterations: " << total_localized_refine_iterations
          << ", Avg. Iterations per Call: " << ((double) total_localized_refine_iterations / (double) total_localized_refine_calls);
    }

    // Top-Level Refinement on all vertices
    const HyperedgeWeight objective_before = current_metrics.getMetric(
      _context.partition.mode, _context.partition.objective);
    const double time_limit = refinementTimeLimit(_context, _round_coarsening_times.back());
    globalRefine(_phg, fm, current_metrics, time_limit);
    _round_coarsening_times.pop_back();
    ASSERT(_round_coarsening_times.size() == 0);
    const HyperedgeWeight objective_after = current_metrics.getMetric(
      _context.partition.mode, _context.partition.objective);
    if ( _context.partition.verbose_output && objective_after < objective_before ) {
      LOG << GREEN << "Top-Level Refinment improved objective from"
          << objective_before << "to" << objective_after << END;
    }

    if ( is_timer_disabled ) {
      utils::Timer::instance().enable();
    }

    // If we finish batch uncontractions and partition is imbalanced, we try to rebalance it
    if ( _top_level && !metrics::isBalanced(_phg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
            << metrics::imbalance(_phg, _context) << ") ->"
            << "Rebalancer is activated" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(_phg, _context);
      }

      utils::Timer::instance().start_timer("rebalance", "Rebalance");
      if ( _context.partition.objective == kahypar::Objective::km1 ) {
        Km1Rebalancer rebalancer(_phg, _context);
        rebalancer.rebalance(current_metrics);
      } else if ( _context.partition.objective == kahypar::Objective::cut ) {
        CutRebalancer rebalancer(_phg, _context);
        rebalancer.rebalance(current_metrics);
      }
      utils::Timer::instance().stop_timer("rebalance");

      const HyperedgeWeight quality_after = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        const HyperedgeWeight quality_delta = quality_after - quality_before;
        if ( quality_delta > 0 ) {
          LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
              << "(Current Imbalance:" << metrics::imbalance(_phg, _context) << ")" << END;
        } else {
          LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
              << "(Current Imbalance:" << metrics::imbalance(_phg, _context) << ")" << END;
        }
      }
    }

    double avg_batch_size = static_cast<double>(total_batches_size) / num_batches;
    utils::Stats::instance().add_stat("num_batches", static_cast<int64_t>(num_batches));
    utils::Stats::instance().add_stat("avg_batch_size", avg_batch_size);
    DBG << V(num_batches) << V(avg_batch_size);

    ASSERT(metrics::objective(_phg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
           V(metrics::objective(_phg, _context.partition.objective)));

    return std::move(_phg);
  }


  void
  NLevelCoarsenerBase::uncontractGroupAsyncSubtask(const ds::ContractionGroup &group,
                                                   const ds::ContractionGroupID groupID) {

      _phg.uncontract(group, groupID, _context.uncoarsening.use_bitcopy_snapshots);

      auto repr_part_id = _phg.partID(group.getRepresentative());
      unused(repr_part_id);
      ASSERT(repr_part_id != kInvalidPartition);
      ASSERT(std::all_of(ds::GroupNodeIDIterator::getAtBegin(group),
                         ds::GroupNodeIDIterator::getAtEnd(group),
                         [&](const HypernodeID &hn) {
                             return _hg.nodeIsEnabled(hn) && _phg.partID(hn) != kInvalidPartition;
                         }),
             "After uncontracting a group, either the representative or any of the contracted nodes is not enabled or not assigned a partition!");
  }

  void NLevelCoarsenerBase::uncoarsenAsyncTask(TreeGroupPool *pool, pq_type& pq, metrics::ThreadSafeMetrics &current_metrics,
                                               IAsyncRefiner *async_lp_refiner, IAsyncRefiner *async_fm_refiner,
                                               HypernodeID &uncontraction_counter,
                                               utils::ProgressBar &uncontraction_progress,
                                               AsyncNodeTracker &async_node_tracker,
                                               RegionComparator &node_region_comparator,
                                               SeedDeduplicator &seed_deduplicator, const size_t task_id,
                                               const bool alwaysInsertIntoPQ, size_t &local_calls_to_localized_refine,
                                               size_t &local_iterations_in_localized_refine) {

      if (pool->taskFinished()) return;

      ds::ContractionGroupID groupID = ds::invalidGroupID;
      bool pick_new_group = true;

      const bool use_old_seeds = _context.uncoarsening.use_old_refinement_seeds && _context.type == kahypar::ContextType::main;

      const size_t seed_nodes_step_length = _context.refinement.fm.num_seed_nodes;
      const size_t min_seed_nodes = _context.partition.k * seed_nodes_step_length;
      auto local_refinement_nodes = ds::NoDownsizeIntegralTypeVector<HypernodeID>(2 * min_seed_nodes, kInvalidHypernode);
      auto num_edges_activated_per_refinement_node = ds::NoDownsizeIntegralTypeVector<HyperedgeID>(2 * min_seed_nodes, kInvalidHyperedge);
      bool reached_min_seeds_once = false;

      auto handle = pq.get_handle();
      if (task_id == 0) {
        for (const auto& root : pool->getPtrToHierarchyForQueries()->roots()) {
          pool->activate(root, handle);
        }
      }
      while (true) {

        // Attempt to pick a group from the pool. Stop if pool is completed.
        while (pick_new_group && !pool->tryToPickActiveID(groupID, task_id, handle)) {
          if (pool->taskFinished()) return;
        }

          ASSERT(groupID != ds::invalidGroupID);
          if (groupID == ds::invalidGroupID || groupID >= pool->getNumTotal()) {
            ERROR("groupID is invalid!" << V(groupID));
          }

        const ds::ContractionGroup &group = pool->group(groupID);

          // Attempt to acquire lock for representative of the group. If the lock cannot be
          // acquired, revert to previous state and attempt to pick an id again
          bool acquired = _lock_manager_for_async->tryToAcquireLock(group.getRepresentative(), groupID);
          if (!acquired) {

            pool->activate(groupID, handle);
            pick_new_group = true;
            continue;
          }
          ASSERT(acquired);

          pool->markAccepted(groupID);
          const bool repIsBorderNodeBeforeUncontraction = _phg.isBorderNode(group.getRepresentative());
          HyperedgeID num_edges_activated_in_task = 0;
          if (repIsBorderNodeBeforeUncontraction) {
            node_region_comparator.markActive(_phg.incidentEdges(group.getRepresentative()), task_id, num_edges_activated_in_task);
          }

          uncontractGroupAsyncSubtask(group, groupID);

          // Release lock (Locks will be reacquired for moves during refinement)
          _lock_manager_for_async->strongReleaseLock(group.getRepresentative(), groupID);

          if (!async_lp_refiner || !async_fm_refiner) {
            ERROR("Local Async Refiners are nullptr!");
          }

          // Extract refinement seeds
          size_t num_extracted_seeds = 0;
          auto begin = ds::GroupNodeIDIterator::getAtBegin(group);
          auto end = ds::GroupNodeIDIterator::getAtEnd(group);
          for (auto it = begin; it != end; ++it) {
            HypernodeID hn = *it;
            ASSERT(_phg.nodeIsEnabled(hn));
            if (_phg.isBorderNode(hn)) {
              // Deduplicate
              if (!seed_deduplicator[hn]) {
                ++num_extracted_seeds;
                seed_deduplicator.set(hn);
                local_refinement_nodes.push_back(hn);
                num_edges_activated_per_refinement_node.push_back(0);
              }
            }
          }

          ASSERT(local_refinement_nodes.size() == num_edges_activated_per_refinement_node.size());

          // Mark edges incident to extracted seeds active for task
          if (num_extracted_seeds > 0) {
            ASSERT(!num_edges_activated_per_refinement_node.empty());
            ASSERT(num_edges_activated_per_refinement_node.back() == 0);

            for (size_t i = 0; i < num_extracted_seeds; ++i) {
              const HypernodeID seed = *(local_refinement_nodes.end() - i - 1);
              // Skip representative if its incident hyperedges have already been activated before uncontraction
              if (seed == group.getRepresentative() && repIsBorderNodeBeforeUncontraction) continue;
              node_region_comparator.markActive(_phg.incidentEdges(seed), task_id, num_edges_activated_in_task);
            }

            // Give last extracted seed (if any) the entry for number of edges activated
            num_edges_activated_per_refinement_node.back() = num_edges_activated_in_task;
          } else if (repIsBorderNodeBeforeUncontraction) {
            // Representative was a border node before uncontraction but no refinement seeds due to intermediate changes => deactivate hyperedges incident to representative again as it will not be part of refinement
            node_region_comparator.markAllEdgesForTaskInactive(task_id);
          }

          // Refine only once enough seeds are available
          if ((_context.type != kahypar::ContextType::main && !local_refinement_nodes.empty()) // case for Initial Partitioning => refine right away
            || (!use_old_seeds && local_refinement_nodes.size() >= seed_nodes_step_length) // case for not using old seed nodes => refine every step_length seeds
            || (use_old_seeds && local_refinement_nodes.size() >= min_seed_nodes) // case for using old seeds when min_seed_nodes are reached
            || (use_old_seeds && !reached_min_seeds_once && !local_refinement_nodes.empty() && local_refinement_nodes.size() % seed_nodes_step_length == 0))  // case for using old seed nodes for before min_seeds have been reached
          {
            ++local_calls_to_localized_refine;
            reached_min_seeds_once |= (use_old_seeds && local_refinement_nodes.size() >= min_seed_nodes);

            auto cur_seeds = IteratorRange(local_refinement_nodes.cbegin(), local_refinement_nodes.cend());

            localizedRefineForAsync(_phg, cur_seeds, async_lp_refiner, async_fm_refiner, groupID,
                                    current_metrics, local_iterations_in_localized_refine);

            if (use_old_seeds && local_refinement_nodes.size() >= min_seed_nodes) {
              // Deactivate nodes that expired, i.e. the oldest seeds
              size_t num_old_edges_to_deactivate = 0;
              for (size_t i = 0; i < seed_nodes_step_length; ++i) {
                const HypernodeID expired_seed = local_refinement_nodes.pop_front();
                num_old_edges_to_deactivate += num_edges_activated_per_refinement_node.pop_front();
                seed_deduplicator.set(expired_seed, false);
              }
              node_region_comparator.markFirstActivatedEdgesForTaskInactive(task_id, num_old_edges_to_deactivate);
            } else if (!use_old_seeds) {
              // When not using old seed nodes clear all refinement seeds
              local_refinement_nodes.clear();
              num_edges_activated_per_refinement_node.clear();
              seed_deduplicator.reset();
              node_region_comparator.markAllEdgesForTaskInactive(task_id);
            }
            // else if {} : When using old seeds but fewer than min_seed_nodes available, none expire yet
          }

          async_node_tracker.incrementTime(group.size());

          // If number of uncontractions in this thread surpasses threshold, update the uncontraction progress bar
          // and reset counter
          uncontraction_counter += group.size();
          if (uncontraction_counter >= ASYNC_UPDATE_PROGRESS_BAR_THRESHOLD) {
            uncontraction_progress.setObjective(current_metrics.getMetric(
                _context.partition.mode, _context.partition.objective));
            uncontraction_progress += uncontraction_counter;
            uncontraction_counter = 0;
          }

          // If the group has successors, activate the successors (i.e. put them in the queue) and potentially work on one of them.
          // If the group has no successors, simply pick a new group for the next iteration.
          if (pool->numSuccessors(groupID) > 0) {
              auto successors = pool->successors(groupID);
              auto suc_begin = successors.begin();
              auto suc_end = successors.end();

              if (alwaysInsertIntoPQ) {
                HypernodeID num_succ = 0;
                // Insert all successors into the PQ and pick a new group for the next iteration.
                for (auto it = suc_begin; it < suc_end; ++it) {
                  num_succ++;
                  pool->activate(*it, handle);
                }
                pick_new_group = true;
              } else {
                // Insert all but one successor into the PQ and work on the remaining successor. Do not pick a new group next iteration.
                groupID = *suc_begin;
                for (auto it = suc_begin + 1; it < suc_end; ++it) {
                  pool->activate(*it, handle);
                }
                pick_new_group = false;
              }
          } else {
            pick_new_group = true;
          }
      }
  }

  PartitionedHypergraph&& NLevelCoarsenerBase::doAsynchronousUncoarsen(std::unique_ptr<IRefiner>& global_fm) {

      ASSERT(_is_finalized);
      metrics::ThreadSafeMetrics current_metrics = initializeForAsync(_compactified_phg);
      // Used where the non-thread-safe metrics type is needed. Convert ThreadSafeMetrics using unsafeLoadMetrics() and
      // back using unsafeStoreMetrics() only in single-threaded environments!
      kahypar::Metrics tmp_unsafe_metrics;

      // Project partition from compactified hypergraph to original hypergraph
      utils::Timer::instance().start_timer("initialize_partition", "Initialize Partition");
      _phg = PartitionedHypergraph(_context.partition.k, _hg, parallel_tag_t());
      _phg.doParallelForAllNodes([&](const HypernodeID hn) {
          ASSERT(static_cast<size_t>(hn) < _compactified_hn_mapping.size());
          const HypernodeID compactified_hn = _compactified_hn_mapping[hn];
          const PartitionID block_id = _compactified_phg.partID(compactified_hn);
          ASSERT(block_id != kInvalidPartition && block_id < _context.partition.k);
          _phg.setOnlyNodePart(hn, block_id);
      });
      _phg.initializePartition();
      _phg.setSnapshotEdgeSizeThreshold(_context.uncoarsening.snapshot_edge_size_threshold);

      if ( _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache) {
          _phg.initializeGainCache();
      }

      ASSERT(metrics::objective(_compactified_phg, _context.partition.objective) ==
             metrics::objective(_phg, _context.partition.objective),
             V(metrics::objective(_compactified_phg, _context.partition.objective)) <<
                                                                                    V(metrics::objective(_phg, _context.partition.objective)));
      ASSERT(metrics::imbalance(_compactified_phg, _context) ==
             metrics::imbalance(_phg, _context),
             V(metrics::imbalance(_compactified_phg, _context)) <<
                                                                V(metrics::imbalance(_phg, _context)));
      utils::Timer::instance().stop_timer("initialize_partition");

      utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                                _context.partition.objective == kahypar::Objective::km1 ? current_metrics.loadKm1() : current_metrics.loadCut(),
                                                _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug);
      uncontraction_progress += _compactified_hg.initialNumNodes();

      // Initialize Refiner for global refinement
      if ( global_fm ) {
        global_fm->initialize(_phg);
      }

      // Perform uncontractions version by version
      bool is_timer_disabled = false;
      bool force_measure_timings = _context.partition.measure_detailed_uncontraction_timings && _top_level;
      if ( utils::Timer::instance().isEnabled() ) {
          utils::Timer::instance().disable();
          is_timer_disabled = true;
      }

      ASSERT(_round_coarsening_times.size() == _removed_hyperedges_batches.size());
      _round_coarsening_times.push_back(_round_coarsening_times.size() > 0 ?
                                        _round_coarsening_times.back() : std::numeric_limits<double>::max()); // Sentinel

    auto uncoarsen_tg = tbb::task_group();
    size_t total_uncontractions = 0;

    // Thread specific asynchronous refiners.
    // Indirectly managed through unique_ptr's so the type of gain policy can be abstracted still.
    auto fm_shared_data = std::make_unique<AsyncFMSharedData>(_phg.initialNumNodes(), _context);
    auto node_anti_duplicator = std::make_unique<ds::ThreadSafeFlagArray<HypernodeID>>(_phg.initialNumNodes());
    auto edge_anti_duplicator = std::make_unique<ds::ThreadSafeFlagArray<HyperedgeID>>(_phg.initialNumEdges());

    const size_t num_threads = _context.shared_memory.num_threads;

    auto async_lp_refiners = std::vector<std::unique_ptr<IAsyncRefiner>>(num_threads);
    auto async_fm_refiners = std::vector<std::unique_ptr<IAsyncRefiner>>(num_threads);
    auto seed_deduplicator_arrays = std::vector<std::unique_ptr<SeedDeduplicator>>(num_threads);
    tbb::parallel_for(size_t(0), num_threads, [&](const size_t t) {
      async_lp_refiners[t] = AsyncLPRefinerFactory::getInstance().createObject(
              _context.refinement.label_propagation.algorithm,
              _phg.hypergraph(),
              _context,
              _lock_manager_for_async.get(),
              node_anti_duplicator.get(),
              edge_anti_duplicator.get(),
              &fm_shared_data->nodeTracker);
      async_fm_refiners[t] = AsyncFMRefinerFactory::getInstance().createObject(
              _context.refinement.fm.algorithm,
              _phg.hypergraph(),
              _context,
              _lock_manager_for_async.get(),
              *fm_shared_data);
      seed_deduplicator_arrays[t] = std::make_unique<kahypar::ds::FastResetFlagArray<HypernodeID>>(_phg.initialNumNodes());
    });

    auto async_uncontraction_counters = ds::Array<HypernodeID>(num_threads, 0);

    auto node_region_comparator = RegionComparator(_phg.hypergraph(), _context.uncoarsening.node_region_similarity_threshold, num_threads);
    size_t total_calls_to_pick = 0;
    size_t calls_to_pick_that_reached_max_retries = 0;
    size_t calls_to_pick_with_empty_pq = 0;

    auto calls_to_localized_refine_per_task = std::vector<size_t>(num_threads, 0);
    auto iterations_of_localized_refine_per_task = std::vector<size_t>(num_threads, 0);
      for (size_t inv_version = 0; inv_version < _group_pools_for_versions.size(); ++inv_version) {

        size_t version = _group_pools_for_versions.size() - inv_version - 1;

          utils::Timer::instance().start_timer("uncontract_and_refine_version",
                                               "Uncontracting and Refining Version", false,
                                               force_measure_timings);

          ASSERT(_phg.version() == _removed_hyperedges_batches.size());
          TreeGroupPool* pool = _group_pools_for_versions[version].get();
          ASSERT(_phg.version() == pool->getVersion());
          ASSERT(_phg.version() == _hg.version());
          ASSERT(_phg.version() == _phg.hypergraph().version());

          _phg.hypergraph().sortStableActivePinsToBeginning();

          pool->setNodeRegionComparator(&node_region_comparator);

          auto num_uncontractions = static_cast<size_t>(pool->getTotalNumUncontractions());

          for (auto& seed_deduplicator : seed_deduplicator_arrays) {
            seed_deduplicator.get()->reset(); // NOLINT(readability-redundant-smartptr-get)
          }
#if defined PQ_MQ_RANDOM || defined PQ_MQ_STICKY || defined PQ_MQ_SWAPPING
          auto params = pq_type::param_type{};
          params.c = MQ_C;
#ifndef PQ_MQ_RANDOM
          params.stickiness = MQ_STICKINESS;
#endif
          auto pq = pq_type(1'000'000, _context.shared_memory.num_threads, params);
#else
          auto pq = pq_type(1'000'000, _context.shared_memory.num_threads);
#endif

          auto uncoarsen_task = [&](const size_t task_id){
              // Setting alwaysInsertIntoPQ to true for Initial Partitioning as well cryptically does not work (spits out Seg faults only in release mode). I assume it's due to TBB internals.
              uncoarsenAsyncTask(pool, pq,
                                 current_metrics, async_lp_refiners[task_id].get(), async_fm_refiners[task_id].get(),
                                 async_uncontraction_counters[task_id],
                                 uncontraction_progress, fm_shared_data->nodeTracker, node_region_comparator,
                                 *seed_deduplicator_arrays[task_id], task_id,
                                 _context.uncoarsening.always_insert_groups_into_pq && (_context.type ==
                                                                                        kahypar::ContextType::main) /* Do not use this option in initial partitioning*/
                      , calls_to_localized_refine_per_task[task_id], iterations_of_localized_refine_per_task[task_id]);
          };

          for (size_t tid = 0; tid < num_threads; ++tid) {
              uncoarsen_tg.run([uncoarsen_task, tid] { return uncoarsen_task(tid); });
          }
          uncoarsen_tg.wait();

          // clean up RegionComparator
          for (size_t tid = 0; tid < num_threads; ++tid) {
            node_region_comparator.markAllEdgesForTaskInactive(tid);
          }

          total_calls_to_pick += pool->getTotalCallsToPick();
          calls_to_pick_that_reached_max_retries += pool->getCallsToPickThatReachedMaxRetries();
          calls_to_pick_with_empty_pq += pool->getCallsToPickWithEmptyPQ();


          // Update Progress Bar With Rest of Uncontractions
          HypernodeID uncontractions_not_in_progress_bar = 0;
          for (HypernodeID& counter : async_uncontraction_counters) {
            uncontractions_not_in_progress_bar += counter;
            counter = 0;
          }
          uncontraction_progress.setObjective(current_metrics.getMetric(
              _context.partition.mode, _context.partition.objective));
          uncontraction_progress += uncontractions_not_in_progress_bar;

          total_uncontractions += num_uncontractions;
//          ASSERT(total_uncontractions == uncontraction_progress.count() - _compactified_hg.initialNumNodes());

//          HEAVY_REFINEMENT_ASSERT(node_anti_duplicator->checkAllFalse());
//          HEAVY_REFINEMENT_ASSERT(edge_anti_duplicator->checkAllFalse());
          HEAVY_REFINEMENT_ASSERT(pool->allAccepted(), V(pool->getNumAccepted()) << V(pool->getTotalNumUncontractions()));
          HEAVY_REFINEMENT_ASSERT(node_region_comparator.checkNoneActiveParallel());
          HEAVY_REFINEMENT_ASSERT(fm_shared_data->nodeTracker.checkNoneLocked());
          HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
          HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());
          HEAVY_REFINEMENT_ASSERT(metrics::objective(_phg, _context.partition.objective) ==
                                  current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                                  V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective))
                                  << V(metrics::objective(_phg, _context.partition.objective)));

          utils::Timer::instance().stop_timer("uncontract_and_refine_version", force_measure_timings);

          // Restore single-pin and parallel nets to continue with the next version
          if ( !_removed_hyperedges_batches.empty() ) {
              utils::Timer::instance().start_timer("restore_single_pin_and_parallel_nets",
                                                   "Restore Single Pin and Parallel Nets", false,
                                                   force_measure_timings);
              _phg.restoreSinglePinAndParallelNets(_removed_hyperedges_batches.back());
              _removed_hyperedges_batches.pop_back();
              utils::Timer::instance().stop_timer("restore_single_pin_and_parallel_nets", force_measure_timings);
              HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
              HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());

              // Perform refinement on all vertices
              const double time_limit = refinementTimeLimit(_context, _round_coarsening_times.back());
              tmp_unsafe_metrics = current_metrics.unsafeLoadMetrics();
              globalRefine(_phg, global_fm, tmp_unsafe_metrics, time_limit);
              current_metrics.unsafeStoreMetrics(tmp_unsafe_metrics);
              uncontraction_progress.setObjective(current_metrics.getMetric(
                      _context.partition.mode, _context.partition.objective));
              _round_coarsening_times.pop_back();
          }


//          _group_pools_for_versions.pop_back();
      }

      node_anti_duplicator.reset();
      edge_anti_duplicator.reset();

      int64_t total_lp_attempted_moves = 0;
      int64_t total_lp_moved_nodes = 0;
      // Calculate total LP attempted moves and actual moves and update stats
      for (const auto & async_lp_refiner : async_lp_refiners) {
          total_lp_attempted_moves += async_lp_refiner->getNumTotalAttemptedMoves();
          total_lp_moved_nodes += async_lp_refiner->getNumTotalMovedNodes();
      }

      size_t total_localized_refine_calls = 0;
      size_t total_iterations_of_localized_refine = 0;
      for (size_t i = 0; i < num_threads; ++i) {
        total_localized_refine_calls += calls_to_localized_refine_per_task[i];
        total_iterations_of_localized_refine += iterations_of_localized_refine_per_task[i];
      }

      utils::Stats::instance().update_stat("lp_attempted_moves", static_cast<int64_t>(total_lp_attempted_moves));
      utils::Stats::instance().update_stat("lp_moved_nodes", static_cast<int64_t>(total_lp_moved_nodes));
//      std::cout << "Total LP attempted moves: " << total_lp_attempted_moves << std::endl;
//      std::cout << "Total LP moves: " << total_lp_moved_nodes << std::endl;

      utils::Stats::instance().update_stat("localized_fm_moves", static_cast<int64_t>(fm_shared_data->total_moves));
      utils::Stats::instance().update_stat("localized_fm_local_reverts", static_cast<int64_t>(fm_shared_data->total_reverts));
      utils::Stats::instance().update_stat("localized_fm_total_find_moves_calls", static_cast<int64_t>(fm_shared_data->total_find_moves_calls));
      utils::Stats::instance().update_stat("localized_fm_find_moves_calls_with_good_prefix", static_cast<int64_t>(fm_shared_data->find_moves_calls_with_good_prefix));
      utils::Stats::instance().update_stat("localized_fm_find_move_retries", static_cast<int64_t>(fm_shared_data->find_move_retries));
      utils::Stats::instance().update_stat("localized_fm_total_pushes", static_cast<int64_t>(fm_shared_data->total_pushes_pos_gain + fm_shared_data->total_pushes_non_pos_gain));
      utils::Stats::instance().update_stat("localized_fm_pushes_with_pos_gain", static_cast<int64_t>(fm_shared_data->total_pushes_pos_gain));
      utils::Stats::instance().update_stat("localized_fm_pushes_with_non_pos_gain", static_cast<int64_t>(fm_shared_data->total_pushes_non_pos_gain));
      utils::Stats::instance().update_stat("localized_fm_pins_touched_by_delta_gain_cache_updates", static_cast<int64_t>(fm_shared_data->num_pins_touched_by_delta_gain_cache_updates));

      size_t sum_gc_cases = fm_shared_data->num_case_from_zero_gc_updates + fm_shared_data->num_case_from_one_gc_updates + fm_shared_data->num_case_to_one_gc_updates + fm_shared_data->num_case_to_two_gc_updates;
      utils::Stats::instance().update_stat("localized_fm_total_delta_gain_cache_updates", static_cast<int64_t>(sum_gc_cases));
      utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_from_zero", static_cast<int64_t>(fm_shared_data->num_case_from_zero_gc_updates));
      utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_from_one", static_cast<int64_t>(fm_shared_data->num_case_from_one_gc_updates));
      utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_to_one", static_cast<int64_t>(fm_shared_data->num_case_to_one_gc_updates));
      utils::Stats::instance().update_stat("localized_fm_delta_gain_cache_updates_case_to_two", static_cast<int64_t>(fm_shared_data->num_case_to_two_gc_updates));

      utils::Stats::instance().update_stat("stable_pins_seen", static_cast<int64_t>(_phg.getNumStablePinsSeen()));
      utils::Stats::instance().update_stat("volatile_pins_seen", static_cast<int64_t>(_phg.getNumVolatilePinsSeen()));
      utils::Stats::instance().update_stat("total_calls_to_pick_next_group", static_cast<int64_t>(total_calls_to_pick));
      utils::Stats::instance().update_stat("calls_to_pick_next_group_with_max_retries", static_cast<int64_t>(calls_to_pick_that_reached_max_retries));
      utils::Stats::instance().update_stat("calls_to_pick_next_group_with_empty_pq", static_cast<int64_t>(calls_to_pick_with_empty_pq));


      if ( _context.partition.verbose_output && _context.type == kahypar::ContextType::main) {
        LOG << std::setprecision(5) << std::fixed
                  << "FM moves: " << fm_shared_data->total_moves
                  << ", Reverts: " << fm_shared_data->total_reverts
                  << ", Fraction: " << fm_shared_data->getFractionOfRevertedMoves();
        LOG << std::setprecision(5) << std::fixed
                  << "FM calls: " << fm_shared_data->total_find_moves_calls
                  << ", With Good Prefix: " << fm_shared_data->find_moves_calls_with_good_prefix
                  << ", Fraction: " << fm_shared_data->getFractionOfFMCallsWithGoodPrefix();
        LOG << "FM find move retries: " << fm_shared_data->find_move_retries;
        LOG << std::setprecision(5) << std::fixed
                  << "FM pushes: " << (fm_shared_data->total_pushes_pos_gain + fm_shared_data->total_pushes_non_pos_gain)
                  << ", With pos gain: " << fm_shared_data->total_pushes_pos_gain
                  << ", Fraction: " << fm_shared_data->getFractionOfPosGainPushes();

        size_t num_stable_pins = _phg.getNumStablePinsSeen();
        size_t num_volatile_pins = _phg.getNumVolatilePinsSeen();
        double volatile_rel_to_stable_pins = (double) num_volatile_pins / (double) num_stable_pins;
        LOG << std::setprecision(5) << std::fixed
                  << "Stable pins seen: " << num_stable_pins
                  << ", Volatile pins seen: " << num_volatile_pins
                  << ", Volatile rel to Stable: " << volatile_rel_to_stable_pins;

        double fraction_calls_to_pick_that_reached_max_retries = (double) calls_to_pick_that_reached_max_retries / (double) total_calls_to_pick;
        double fraction_calls_to_pick_with_empty_pq = (double) calls_to_pick_with_empty_pq / (double) total_calls_to_pick;
        LOG << std::setprecision(5) << std::fixed
            << "Total calls to pickActiveID(): " << total_calls_to_pick
            << "\n\t Reached Max Retries: " << calls_to_pick_that_reached_max_retries << ", Fraction: " << fraction_calls_to_pick_that_reached_max_retries
            << "\n\t With Empty PQ: " << calls_to_pick_with_empty_pq << ", Fraction: " << fraction_calls_to_pick_with_empty_pq;
        LOG << "Pins touched in DeltaPHG Gain Cache updates: " << fm_shared_data->num_pins_touched_by_delta_gain_cache_updates;
        LOG << std::setprecision(5) << std::fixed << "Number of DeltaPHG Gain Cache updates triggered: "
        << "\n\tTotal: " << sum_gc_cases
        << "\n\tCase: Pins in from=0: " << fm_shared_data->num_case_from_zero_gc_updates << ", Fraction: " << ((double) fm_shared_data->num_case_from_zero_gc_updates / (double) sum_gc_cases)
        << "\n\tCase: Pins in from=1: " << fm_shared_data->num_case_from_one_gc_updates << ", Fraction: " << ((double) fm_shared_data->num_case_from_one_gc_updates / (double) sum_gc_cases)
        << "\n\tCase: Pins in to=1:   " << fm_shared_data->num_case_to_one_gc_updates << ", Fraction: " << ((double) fm_shared_data->num_case_to_one_gc_updates / (double) sum_gc_cases)
        << "\n\tCase: Pins in to=2:   " << fm_shared_data->num_case_to_two_gc_updates << ", Fraction: " << ((double) fm_shared_data->num_case_to_two_gc_updates / (double) sum_gc_cases);
        LOG << std::setprecision(5) << std::fixed
            << "localizedRefine() calls: " << total_localized_refine_calls
            << ", lR() Loop Iterations: " << total_iterations_of_localized_refine
            << ", Avg. Iterations per Call: " << ((double) total_iterations_of_localized_refine / (double) total_localized_refine_calls);
      }

      size_t total_num_nodes = _hg.initialNumNodes() - _hg.numRemovedHypernodes();
      size_t num_nodes_after_coarsening = _compactified_hg.initialNumNodes();
      unused(total_num_nodes);
      unused(num_nodes_after_coarsening);
      ASSERT(total_num_nodes == total_uncontractions + num_nodes_after_coarsening, V(total_num_nodes) << ", " << V(total_uncontractions) << ", " << V(num_nodes_after_coarsening));
      auto checkAllAssigned = [&](){
          CAtomic<uint8_t> all_assigned(uint8_t(true));
          _hg.doParallelForAllNodes([&](const HypernodeID& i) {
              if(_phg.partID(i) == kInvalidPartition) {
                std::cout << "Node " << i << " has an invalid partition." << std::endl;
                all_assigned.fetch_and(uint8_t(false), std::memory_order_relaxed);
              }
          });
          return (bool) all_assigned.load(std::memory_order_relaxed);
      };
      unused(checkAllAssigned);
      HEAVY_REFINEMENT_ASSERT(checkAllAssigned());
      auto checkNoneLocked = [&](){
        for (HypernodeID i = 0; i < _hg.initialNumNodes(); ++i) {
          if(_lock_manager_for_async->isLocked(i)) return false;
        }
        return true;
      };
      unused(checkNoneLocked);
      HEAVY_REFINEMENT_ASSERT(checkNoneLocked());

      // Top-Level Refinement on all vertices
      const HyperedgeWeight objective_before = current_metrics.getMetric(
              _context.partition.mode, _context.partition.objective);
      const double time_limit = refinementTimeLimit(_context, _round_coarsening_times.back());
      tmp_unsafe_metrics = current_metrics.unsafeLoadMetrics();
      globalRefine(_phg, global_fm, tmp_unsafe_metrics, time_limit);
      current_metrics.unsafeStoreMetrics(tmp_unsafe_metrics);
      _round_coarsening_times.pop_back();
      ASSERT(_round_coarsening_times.size() == 0);
      const HyperedgeWeight objective_after = current_metrics.getMetric(
              _context.partition.mode, _context.partition.objective);
      if ( _context.partition.verbose_output && objective_after < objective_before ) {
          LOG << GREEN << "Top-Level Refinment improved objective from"
              << objective_before << "to" << objective_after << END;
      }

      if ( is_timer_disabled ) {
          utils::Timer::instance().enable();
      }

      // If we finish uncontractions and partition is imbalanced, we try to rebalance it
      if ( _top_level && !metrics::isBalanced(_phg, _context)) {
          const HyperedgeWeight quality_before = current_metrics.getMetric(
                  kahypar::Mode::direct_kway, _context.partition.objective);
          if ( _context.partition.verbose_output ) {
              LOG << RED << "Partition is imbalanced (Current Imbalance:"
                  << metrics::imbalance(_phg, _context) << ") ->"
                  << "Rebalancer is activated" << END;

              LOG << "Part weights: (violations in red)";
              io::printPartWeightsAndSizes(_phg, _context);
          }

          utils::Timer::instance().start_timer("rebalance", "Rebalance");
          tmp_unsafe_metrics = current_metrics.unsafeLoadMetrics();
          if ( _context.partition.objective == kahypar::Objective::km1 ) {
              Km1Rebalancer rebalancer(_phg, _context);
              rebalancer.rebalance(tmp_unsafe_metrics);
          } else if ( _context.partition.objective == kahypar::Objective::cut ) {
              CutRebalancer rebalancer(_phg, _context);
              rebalancer.rebalance(tmp_unsafe_metrics);
          }
          current_metrics.unsafeStoreMetrics(tmp_unsafe_metrics);
          utils::Timer::instance().stop_timer("rebalance");

          const HyperedgeWeight quality_after = current_metrics.getMetric(
                  kahypar::Mode::direct_kway, _context.partition.objective);
          if ( _context.partition.verbose_output ) {
              const HyperedgeWeight quality_delta = quality_after - quality_before;
              if ( quality_delta > 0 ) {
                  LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
                      << "(Current Imbalance:" << metrics::imbalance(_phg, _context) << ")" << END;
              } else {
                  LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
                      << "(Current Imbalance:" << metrics::imbalance(_phg, _context) << ")" << END;
              }
          }
      }

      ASSERT(metrics::objective(_phg, _context.partition.objective) ==
             current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
             V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                                                                                                    V(metrics::objective(_phg, _context.partition.objective)));

      return std::move(_phg);

  }


  void MultilevelCoarsenerBase::refine(
          PartitionedHypergraph& partitioned_hypergraph,
          std::unique_ptr<IRefiner>& label_propagation,
          std::unique_ptr<IRefiner>& fm,
          kahypar::Metrics& current_metrics,
          const double time_limit) {

    if ( debug && _top_level ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
          << ", imbalance = " << current_metrics.imbalance;
    }

    parallel::scalable_vector<HypernodeID> dummy;
    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;

      if ( label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        label_propagation->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_lp_refiner");

        utils::Timer::instance().start_timer("label_propagation", "Label Propagation");
        improvement_found |= label_propagation->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        utils::Timer::instance().stop_timer("label_propagation");
      }

      if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        fm->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_fm_refiner");

        utils::Timer::instance().start_timer("fm", "FM");
        improvement_found |= fm->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        utils::Timer::instance().stop_timer("fm");
      }

      if ( _top_level ) {
        ASSERT(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)
               == metrics::objective(partitioned_hypergraph, _context.partition.objective),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
                               << "does not match the metric updated by the refiners" << V(current_metrics.km1));
      }

      if ( !_context.refinement.refine_until_no_improvement ) {
        break;
      }
    }

    if ( _top_level) {
      DBG << "--------------------------------------------------\n";
    }
  }

  void NLevelCoarsenerBase::localizedRefine(PartitionedHypergraph &partitioned_hypergraph,
                                            const parallel::scalable_vector <HypernodeID> &refinement_nodes,
                                            std::unique_ptr<IRefiner> &label_propagation, std::unique_ptr<IRefiner> &fm,
                                            kahypar::Metrics &current_metrics, const bool force_measure_timings,
                                            size_t &num_refinement_iterations) {
    if ( debug && _top_level ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
          << ", imbalance = " << current_metrics.imbalance;
    }

    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;
      ++num_refinement_iterations;

      if ( label_propagation &&
           _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("label_propagation", "Label Propagation", false, force_measure_timings);
        improvement_found |= label_propagation->refine(partitioned_hypergraph,
          refinement_nodes, current_metrics, std::numeric_limits<double>::max());
        utils::Timer::instance().stop_timer("label_propagation", force_measure_timings);
      }

      if ( fm &&
           _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("fm", "FM", false, force_measure_timings);
        improvement_found |= fm->refine(partitioned_hypergraph,
          refinement_nodes, current_metrics, std::numeric_limits<double>::max());
        utils::Timer::instance().stop_timer("fm", force_measure_timings);
      }

      if ( _top_level ) {
        ASSERT(current_metrics.km1 == metrics::km1(partitioned_hypergraph),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
                               << "does not match the metric updated by the refiners" << V(current_metrics.km1));
      }

      if ( !_context.refinement.refine_until_no_improvement ) {
        break;
      }
    }

    if ( _top_level) {
      DBG << "--------------------------------------------------\n";
    }
  }

  void NLevelCoarsenerBase::localizedRefineForAsync(PartitionedHypergraph &partitioned_hypergraph,
                                                    const IteratorRange<IAsyncRefiner::NodeIteratorT> &refinement_nodes,
                                                    IAsyncRefiner *async_lp, IAsyncRefiner *async_fm,
                                                    ds::ContractionGroupID group_id,
                                                    metrics::ThreadSafeMetrics &current_metrics,
                                                    size_t &num_iterations) {

      bool improvement_found = true;
      while( improvement_found ) {
          improvement_found = false;
          ++num_iterations;

          if (async_lp &&
               _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
              improvement_found |= async_lp->refine(partitioned_hypergraph,
                                                    refinement_nodes, current_metrics, std::numeric_limits<double>::max(), group_id);
          }

          if (async_fm &&
            _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
            improvement_found |= async_fm->refine(partitioned_hypergraph,
                                                refinement_nodes, current_metrics, std::numeric_limits<double>::max(), group_id);
          }

          if ( !_context.refinement.refine_until_no_improvement ) {
              break;
          }
      }

  }


  namespace {
    NLevelGlobalFMParameters applyGlobalFMParameters(const FMParameters& fm,
                                                     const NLevelGlobalFMParameters global_fm) {
      NLevelGlobalFMParameters tmp_global_fm;
      tmp_global_fm.num_seed_nodes = fm.num_seed_nodes;
      tmp_global_fm.obey_minimal_parallelism = fm.obey_minimal_parallelism;
      fm.num_seed_nodes = global_fm.num_seed_nodes;
      fm.obey_minimal_parallelism = global_fm.obey_minimal_parallelism;
      return tmp_global_fm;
    }
  } // namespace

  void NLevelCoarsenerBase::globalRefine(PartitionedHypergraph& partitioned_hypergraph,
                                         std::unique_ptr<IRefiner>& fm,
                                         kahypar::Metrics& current_metrics,
                                         const double time_limit) {
    if ( _context.refinement.global_fm.use_global_fm ) {
      if ( debug && _top_level ) {
        io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
        DBG << "Start Refinement - km1 = " << current_metrics.km1
            << ", imbalance = " << current_metrics.imbalance;
      }

      // Apply global FM parameters to FM context and temporary store old fm context
      NLevelGlobalFMParameters tmp_global_fm = applyGlobalFMParameters(
        _context.refinement.fm, _context.refinement.global_fm);
      bool improvement_found = true;
      while( improvement_found ) {
        improvement_found = false;

        if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
          utils::Timer::instance().start_timer("global_fm", "Global FM", false, _top_level);
          improvement_found |= fm->refine(partitioned_hypergraph, {}, current_metrics, time_limit);
          utils::Timer::instance().stop_timer("global_fm", _top_level);
        }

        if ( _top_level ) {
          ASSERT(current_metrics.km1 == metrics::km1(partitioned_hypergraph),
                "Actual metric" << V(metrics::km1(partitioned_hypergraph))
                                << "does not match the metric updated by the refiners" << V(current_metrics.km1));
        }

        if ( !_context.refinement.global_fm.refine_until_no_improvement ) {
          break;
        }
      }
      // Reset FM context
      applyGlobalFMParameters(_context.refinement.fm, tmp_global_fm);

      if ( _top_level) {
        DBG << "--------------------------------------------------\n";
      }
    }
  }

  kahypar::Metrics MultilevelCoarsenerBase::initialize(PartitionedHypergraph& phg) {
    kahypar::Metrics m = { 0, 0, 0.0 };
    tbb::parallel_invoke([&] {
      m.cut = metrics::hyperedgeCut(phg);
    }, [&] {
      m.km1 = metrics::km1(phg);
    });
    m.imbalance = metrics::imbalance(phg, _context);

    int64_t num_nodes = phg.initialNumNodes();
    int64_t num_edges = phg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", m.cut);
    utils::Stats::instance().add_stat("initial_km1", m.km1);
    utils::Stats::instance().add_stat("initial_imbalance", m.imbalance);
    return m;
  }

  kahypar::Metrics NLevelCoarsenerBase::initialize(PartitionedHypergraph& current_hg) {
    kahypar::Metrics current_metrics = computeMetrics(current_hg);
    int64_t num_nodes = current_hg.initialNumNodes();
    int64_t num_edges = current_hg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", current_metrics.cut);
    utils::Stats::instance().add_stat("initial_km1", current_metrics.km1);
    utils::Stats::instance().add_stat("initial_imbalance", current_metrics.imbalance);
    return current_metrics;
  }

  metrics::ThreadSafeMetrics NLevelCoarsenerBase::initializeForAsync(PartitionedHypergraph& current_hg) {
      metrics::ThreadSafeMetrics current_metrics = computeMetricsForAsync(current_hg);
      int64_t num_nodes = current_hg.initialNumNodes();
      int64_t num_edges = current_hg.initialNumEdges();
      utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
      utils::Stats::instance().add_stat("initial_num_edges", num_edges);
      utils::Stats::instance().add_stat("initial_cut", current_metrics.loadCut());
      utils::Stats::instance().add_stat("initial_km1", current_metrics.loadKm1());
      utils::Stats::instance().add_stat("initial_imbalance", current_metrics.loadImbalance());
      return current_metrics;
  }

}
