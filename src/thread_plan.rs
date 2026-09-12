const PRESSURE_STREAK_TARGET: i8 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ThreadPlan {
    pub(crate) requested_threads: usize,
    pub(crate) parser_threads: usize,
    pub(crate) shared_slots: usize,
    pub(crate) initial_input_threads: usize,
    pub(crate) initial_worker_threads: usize,
    pub(crate) max_input_threads: usize,
    pub(crate) worker_headroom: usize,
    pub(crate) parallel_gzip: bool,
    pub(crate) adaptive: bool,
    pub(crate) budget_overcommit: usize,
    pub(crate) cram_decode_threads: usize,
    pub(crate) cram_order_threads: usize,
    pub(crate) output_writer_threads: usize,
}

impl ThreadPlan {
    pub(crate) fn new(
        requested_threads: usize,
        paired: bool,
        has_gzip: bool,
    ) -> Result<Self, String> {
        Self::with_parser_threads(
            requested_threads,
            if paired { 2 } else { 1 },
            paired,
            has_gzip,
        )
    }

    pub(crate) fn with_parser_threads(
        requested_threads: usize,
        parser_threads: usize,
        paired: bool,
        has_gzip: bool,
    ) -> Result<Self, String> {
        if requested_threads == 0 {
            return Err("Thread budget must be greater than 0".into());
        }
        if parser_threads == 0 {
            return Err("Parser thread count must be greater than 0".into());
        }

        let shared_slots = requested_threads.saturating_sub(parser_threads).max(1);
        let budget_overcommit = parser_threads
            .checked_add(shared_slots)
            .ok_or("Thread accounting overflow")?
            .saturating_sub(requested_threads);

        let parallel_gzip = has_gzip && shared_slots >= 4;

        let (initial_input_threads, initial_worker_threads, max_input_threads, adaptive) =
            if parallel_gzip {
                let desired_input = if paired {
                    (shared_slots / 3).max(2)
                } else {
                    shared_slots.div_ceil(4).max(2)
                };

                // Routing and output compression are separate live stages.
                // Reserve one shared slot for each before allowing adaptive
                // input decompression to grow.
                let max_input_threads = shared_slots - 2;
                let initial_input_threads = desired_input.min(max_input_threads);
                let initial_worker_threads = shared_slots - initial_input_threads;

                (
                    initial_input_threads,
                    initial_worker_threads,
                    max_input_threads,
                    max_input_threads > 1,
                )
            } else {
                (0, shared_slots, 0, false)
            };

        let worker_headroom = if adaptive {
            // Input decompression never drops below one active worker, so the
            // remaining shared slots are the maximum routing/compression
            // allocation that can become active.
            shared_slots - 1
        } else {
            initial_worker_threads
        };

        Ok(Self {
            requested_threads,
            parser_threads,
            shared_slots,
            initial_input_threads,
            initial_worker_threads,
            max_input_threads,
            worker_headroom,
            parallel_gzip,
            adaptive,
            budget_overcommit,
            cram_decode_threads: 0,
            cram_order_threads: 0,
            output_writer_threads: 0,
        })
    }

    pub(crate) fn for_cram(requested_threads: usize, output_is_cram: bool) -> Result<Self, String> {
        Self::for_cram_with_writer_threads(requested_threads, output_is_cram, None)
    }

    pub(crate) fn for_cram_with_writer_threads(
        requested_threads: usize,
        output_is_cram: bool,
        writer_override: Option<usize>,
    ) -> Result<Self, String> {
        if requested_threads == 0 {
            return Err("Thread budget must be greater than 0".into());
        }
        if writer_override == Some(0) {
            return Err("CRAM writer thread count must be greater than zero".into());
        }

        if requested_threads == 1 {
            // CRAM has a dedicated inline path at this budget: the caller's
            // thread reads, routes, and writes each item without stage threads.
            return Ok(Self {
                requested_threads,
                parser_threads: 1,
                shared_slots: 0,
                initial_input_threads: 0,
                initial_worker_threads: 0,
                max_input_threads: 0,
                worker_headroom: 0,
                parallel_gzip: false,
                adaptive: false,
                budget_overcommit: 0,
                cram_decode_threads: 0,
                cram_order_threads: 0,
                output_writer_threads: 0,
            });
        }

        let output_writer_threads = if output_is_cram {
            writer_override.unwrap_or(if requested_threads >= 8 { 2 } else { 1 })
        } else {
            1
        };
        let cram_order_threads = usize::from(output_is_cram && output_writer_threads > 1);
        if output_is_cram && output_writer_threads > 1 {
            let minimum = 1usize
                .checked_add(cram_order_threads)
                .and_then(|value| value.checked_add(output_writer_threads))
                .and_then(|value| value.checked_add(1))
                .ok_or("CRAM thread accounting overflow")?;
            if requested_threads < minimum {
                return Err(format!(
                    "CRAM output requires at least {minimum} total threads for {output_writer_threads} writers"
                ));
            }
        }
        let available_after_fixed = requested_threads
            .saturating_sub(1)
            .saturating_sub(cram_order_threads)
            .saturating_sub(output_writer_threads);
        let requested_decode_threads = if requested_threads < 4 {
            0
        } else if output_is_cram {
            // Controlled CRAM-output measurements found that an HTSlib decode
            // worker reduced whole-pipeline throughput. Keep decoding inline
            // so the global budget is available to routing and output.
            0
        } else {
            // Whole-pipeline measurements favored roughly one quarter of the
            // global budget for CRAM decoding, with no benefit beyond eight
            // workers on a 24-thread qualification host.
            (requested_threads / 4).min(8)
        };
        let cram_decode_threads =
            requested_decode_threads.min(available_after_fixed.saturating_sub(1));
        let initial_worker_threads = available_after_fixed
            .saturating_sub(cram_decode_threads)
            .max(1);
        let accounted = 1usize
            .checked_add(cram_decode_threads)
            .and_then(|value| value.checked_add(initial_worker_threads))
            .and_then(|value| value.checked_add(cram_order_threads))
            .and_then(|value| value.checked_add(output_writer_threads))
            .ok_or("CRAM thread accounting overflow")?;

        Ok(Self {
            requested_threads,
            parser_threads: 1,
            shared_slots: initial_worker_threads,
            initial_input_threads: 0,
            initial_worker_threads,
            max_input_threads: 0,
            worker_headroom: initial_worker_threads,
            parallel_gzip: false,
            adaptive: false,
            budget_overcommit: accounted.saturating_sub(requested_threads),
            cram_decode_threads,
            cram_order_threads,
            output_writer_threads,
        })
    }

    pub(crate) fn adaptive_allocation(self) -> Option<AdaptiveAllocation> {
        self.adaptive.then_some(AdaptiveAllocation {
            shared_slots: self.shared_slots,
            input_threads: self.initial_input_threads,
            max_input_threads: self.max_input_threads,
            pressure_streak: 0,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AllocationSnapshot {
    pub(crate) input_threads: usize,
    pub(crate) worker_threads: usize,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct AdaptiveAllocation {
    shared_slots: usize,
    input_threads: usize,
    max_input_threads: usize,
    pressure_streak: i8,
}

impl AdaptiveAllocation {
    pub(crate) fn current(self) -> AllocationSnapshot {
        AllocationSnapshot {
            input_threads: self.input_threads,
            worker_threads: self.shared_slots - self.input_threads,
        }
    }

    pub(crate) fn observe(
        &mut self,
        queued_batches: usize,
        queue_capacity: usize,
    ) -> Option<AllocationSnapshot> {
        if queue_capacity == 0 {
            return None;
        }

        let low_watermark = queue_capacity / 4;
        let high_watermark = queue_capacity.saturating_sub(low_watermark);

        if queued_batches >= high_watermark {
            self.pressure_streak = self.pressure_streak.saturating_add(1);
        } else if queued_batches <= low_watermark {
            self.pressure_streak = self.pressure_streak.saturating_sub(1);
        } else {
            self.pressure_streak = 0;
            return None;
        }

        if self.pressure_streak >= PRESSURE_STREAK_TARGET {
            self.pressure_streak = 0;

            if self.input_threads > 1 {
                self.input_threads -= 1;
                return Some(self.current());
            }

            return None;
        }

        if self.pressure_streak <= -PRESSURE_STREAK_TARGET {
            self.pressure_streak = 0;

            if self.input_threads < self.max_input_threads
                && self.shared_slots - self.input_threads > 1
            {
                self.input_threads += 1;
                return Some(self.current());
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eight_thread_single_end_budget_counts_parser_thread() {
        let plan = ThreadPlan::new(8, false, true).unwrap();
        assert_eq!(
            (
                plan.parser_threads,
                plan.initial_input_threads,
                plan.initial_worker_threads,
                plan.worker_headroom,
            ),
            (1, 2, 5, 6)
        );
        assert!(plan.parallel_gzip);
        assert!(plan.adaptive);
    }

    #[test]
    fn eight_thread_paired_budget_counts_both_parsers() {
        let plan = ThreadPlan::new(8, true, true).unwrap();
        assert_eq!(
            (
                plan.parser_threads,
                plan.initial_input_threads,
                plan.initial_worker_threads,
                plan.worker_headroom,
            ),
            (2, 2, 4, 5)
        );
        assert!(plan.parallel_gzip);
        assert!(plan.adaptive);
    }

    #[test]
    fn sixty_four_thread_server_plans_scale() {
        let se = ThreadPlan::new(64, false, true).unwrap();
        let pe = ThreadPlan::new(64, true, true).unwrap();

        assert_eq!(
            (
                se.parser_threads,
                se.initial_input_threads,
                se.initial_worker_threads,
            ),
            (1, 16, 47)
        );

        assert_eq!(
            (
                pe.parser_threads,
                pe.initial_input_threads,
                pe.initial_worker_threads,
            ),
            (2, 20, 42)
        );
    }

    #[test]
    fn small_budgets_skip_rapidgzip() {
        let se = ThreadPlan::new(4, false, true).unwrap();
        let pe = ThreadPlan::new(4, true, true).unwrap();

        assert_eq!(
            (
                se.parser_threads,
                se.initial_input_threads,
                se.initial_worker_threads,
            ),
            (1, 0, 3)
        );
        assert_eq!(
            (
                pe.parser_threads,
                pe.initial_input_threads,
                pe.initial_worker_threads,
            ),
            (2, 0, 2)
        );
        assert!(!se.parallel_gzip);
        assert!(!pe.parallel_gzip);
    }

    #[test]
    fn cram_fastq_thread_plan_uses_bounded_quarter_for_decode() {
        let expected = [
            (1, 0, 0, 0),
            (2, 0, 1, 1),
            (4, 1, 1, 0),
            (8, 2, 4, 0),
            (12, 3, 7, 0),
            (16, 4, 10, 0),
            (24, 6, 16, 0),
            (32, 8, 22, 0),
            (64, 8, 54, 0),
        ];
        for (budget, decode, workers, overcommit) in expected {
            let plan = ThreadPlan::for_cram(budget, false).unwrap();
            assert_eq!(
                (
                    plan.cram_decode_threads,
                    plan.initial_worker_threads,
                    plan.output_writer_threads,
                    plan.budget_overcommit,
                ),
                (decode, workers, usize::from(budget > 1), overcommit),
                "budget {budget}"
            );
        }
    }

    #[test]
    fn cram_output_thread_plan_keeps_decode_inline() {
        let expected = [
            (1, 0, 0, 0, 0, 0),
            (2, 0, 1, 1, 0, 1),
            (4, 0, 2, 1, 0, 0),
            (8, 0, 4, 2, 1, 0),
            (12, 0, 8, 2, 1, 0),
            (16, 0, 12, 2, 1, 0),
            (24, 0, 20, 2, 1, 0),
            (32, 0, 28, 2, 1, 0),
            (64, 0, 60, 2, 1, 0),
        ];
        for (budget, decode, workers, writers, order, overcommit) in expected {
            let plan = ThreadPlan::for_cram(budget, true).unwrap();
            assert_eq!(
                (
                    plan.cram_decode_threads,
                    plan.initial_worker_threads,
                    plan.output_writer_threads,
                    plan.cram_order_threads,
                    plan.budget_overcommit,
                ),
                (decode, workers, writers, order, overcommit),
                "budget {budget}"
            );
        }
    }

    #[test]
    fn cram_output_writer_experiment_stays_within_global_budget() {
        let expected = [(2, 20), (4, 18), (6, 16), (8, 14)];
        for (writers, routing_workers) in expected {
            let plan = ThreadPlan::for_cram_with_writer_threads(24, true, Some(writers)).unwrap();
            assert_eq!(plan.cram_decode_threads, 0);
            assert_eq!(plan.cram_order_threads, 1);
            assert_eq!(plan.output_writer_threads, writers);
            assert_eq!(plan.initial_worker_threads, routing_workers);
            assert_eq!(plan.budget_overcommit, 0);
        }
    }

    #[test]
    fn cram_output_writer_experiment_rejects_impossible_budget() {
        assert!(ThreadPlan::for_cram_with_writer_threads(4, true, Some(2)).is_err());
        assert!(ThreadPlan::for_cram_with_writer_threads(8, true, Some(0)).is_err());
        assert_eq!(
            ThreadPlan::for_cram_with_writer_threads(1, true, Some(2))
                .unwrap()
                .output_writer_threads,
            0
        );
    }

    #[test]
    fn high_queue_pressure_moves_slot_to_workers() {
        let plan = ThreadPlan::new(8, true, true).unwrap();
        let mut allocation = plan.adaptive_allocation().unwrap();

        assert!(allocation.observe(16, 16).is_none());
        assert!(allocation.observe(16, 16).is_none());
        assert_eq!(
            allocation.observe(16, 16),
            Some(AllocationSnapshot {
                input_threads: 1,
                worker_threads: 5,
            })
        );
    }

    #[test]
    fn low_queue_pressure_moves_slot_to_input() {
        let plan = ThreadPlan::new(8, false, true).unwrap();
        let mut allocation = plan.adaptive_allocation().unwrap();

        assert!(allocation.observe(0, 16).is_none());
        assert!(allocation.observe(0, 16).is_none());
        assert_eq!(
            allocation.observe(0, 16),
            Some(AllocationSnapshot {
                input_threads: 3,
                worker_threads: 4,
            })
        );
    }

    #[test]
    fn adaptive_input_growth_preserves_two_output_pipeline_workers() {
        let plan = ThreadPlan::new(8, false, true).unwrap();
        let mut allocation = plan.adaptive_allocation().unwrap();

        for _ in 0..12 {
            allocation.observe(0, 16);
        }

        let current = allocation.current();
        assert_eq!(current.input_threads, 5);
        assert_eq!(current.worker_threads, 2);
    }
}
