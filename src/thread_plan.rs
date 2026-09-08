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
}

impl ThreadPlan {
    pub(crate) fn new(
        requested_threads: usize,
        paired: bool,
        has_gzip: bool,
    ) -> Result<Self, String> {
        if requested_threads == 0 {
            return Err("Thread budget must be greater than 0".into());
        }

        let parser_threads = if paired { 2 } else { 1 };

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

                let max_input_threads = shared_slots - 1;
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
            requested_threads
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
            (1, 2, 5, 8)
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
            (2, 2, 4, 8)
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
}
