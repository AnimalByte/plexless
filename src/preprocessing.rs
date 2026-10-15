//! Shared preprocessing for demux and trim. Coordinates always refer to input.
use crate::cli::PreprocessingArgs;
use crate::overlap::{self, OverlapResult, OverlapScratch};
use crate::pass_through::{BiologicalBounds, PassThroughPlans, detect_pass_through};
use crate::stats::FastqStats;
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Default)]
pub(crate) struct Preprocessor {
    pub adapters: PassThroughPlans,
    pub options: PreprocessingArgs,
    pub qc: bool,
    index_lists: [Vec<Vec<u8>>; 2],
}

#[derive(Debug, Default)]
pub(crate) struct PreprocessingMetrics {
    pub processed: u64,
    pub retained: u64,
    pub rejected: u64,
    pub assigned_rejected: u64,
    reasons: [u64; 8],
    pending_index: bool,
    pending_duplicate: bool,
    pending_merged: bool,
    last_failed_mates: [bool; 2],
    pub corrected_bases: u64,
    pub merged_fragments: u64,
    pub duplicate_fragments: u64,
    merged_stats: FastqStats,
    detailed: [crate::distribution::Distribution; 2],
    overrep: [Option<crate::overrepresentation::Overrepresentation>; 2],
    merged_detailed: crate::distribution::Distribution,
    merged_overrep: Option<crate::overrepresentation::Overrepresentation>,
    pub scanned: u64,
    pub overlap_confirmed: u64,
    pub overlap_only: [u64; 2],
    pub overlap_only_bases: [u64; 2],
    pub gapped_pairs: u64,
    pub dimers: [u64; 2],
    pub overlap_rejected: u64,
    pub overlap_ambiguous: u64,
    overlap_scratch: OverlapScratch,
    cached_overlap: Option<([BiologicalBounds; 2], OverlapResult)>,
    pub full: [u64; 2],
    pub partial: [u64; 2],
    pub full_bases: [u64; 2],
    pub partial_bases: [u64; 2],
    proximal: [u64; 2],
    quality_removed: [u64; 2],
    poly_g_removed: [u64; 2],
    poly_x_removed: [u64; 2],
    max_length_removed: [u64; 2],
    after: [FastqStats; 2],
    // Reused across mates/records, grown at most once for the longest signature
    // in a worker batch. No per-read adapter construction.
    scratch: Vec<u64>,
}

impl Preprocessor {
    pub(crate) fn new(
        adapters: PassThroughPlans,
        options: PreprocessingArgs,
        qc: bool,
    ) -> Result<Self, String> {
        let mut index_lists = [Vec::new(), Vec::new()];
        for (mate, path) in [
            options.fastp.filter_by_index1.as_deref(),
            options.fastp.filter_by_index2.as_deref(),
        ]
        .into_iter()
        .enumerate()
        {
            if let Some(path) = path {
                use std::io::BufRead;
                for line in std::io::BufReader::new(
                    std::fs::File::open(path)
                        .map_err(|e| format!("Could not open index blacklist: {e}"))?,
                )
                .lines()
                {
                    let line = line.map_err(|e| e.to_string())?;
                    let sequence = crate::technical::dna(&line)?;
                    if index_lists[mate].contains(&sequence) {
                        return Err("Duplicate index blacklist sequence".into());
                    }
                    index_lists[mate].push(sequence);
                }
                if index_lists[mate].is_empty() {
                    return Err("Index blacklist is empty".into());
                }
            }
        }
        Ok(Self {
            adapters,
            options,
            qc,
            index_lists,
        })
    }
    pub(crate) fn check_indexes(
        &self,
        ids: [&[u8]; 2],
        paired: bool,
        metrics: &mut PreprocessingMetrics,
    ) -> Result<(), String> {
        for mate in 0..2 {
            if self.index_lists[mate].is_empty() {
                continue;
            }
            let index = crate::fastp_runtime::validated_index(
                ids[if paired { mate } else { 0 }],
                mate == 1,
            )?;
            if self.index_lists[mate].iter().any(|expected| {
                index.len() == expected.len()
                    && index
                        .iter()
                        .zip(expected)
                        .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
                        .count()
                        <= self.options.fastp.filter_by_index_threshold
            }) {
                metrics.pending_index = true;
            }
        }
        Ok(())
    }
    fn prepare(
        &self,
        mate: usize,
        seq: &[u8],
        qual: &[u8],
        start: usize,
        metrics: &mut PreprocessingMetrics,
    ) -> (BiologicalBounds, usize) {
        let start = start.min(seq.len());
        let start = start
            .saturating_add(crate::fastp_runtime::umi_bases(&self.options.fastp, mate))
            .min(seq.len());
        metrics.proximal[mate] += start as u64;
        let initial = BiologicalBounds {
            start,
            end: seq.len(),
        };
        let mut bounds =
            crate::quality_trim::global_and_quality(seq, qual, initial, mate, &self.options.fastp);
        metrics.quality_removed[mate] += (initial.len() - bounds.len()) as u64;
        if self.options.fastp.trim_poly_g {
            let end = bounds.end;
            bounds = crate::quality_trim::poly_g(seq, bounds, self.options.fastp.poly_g_min_len);
            metrics.poly_g_removed[mate] += (end - bounds.end) as u64;
        }
        (bounds, bounds.start - start)
    }
    fn finish(
        &self,
        seq: &[u8],
        mut bounds: BiologicalBounds,
        mate: usize,
        metrics: &mut PreprocessingMetrics,
    ) -> BiologicalBounds {
        if self.options.fastp.trim_poly_x {
            let end = bounds.end;
            bounds = crate::quality_trim::poly_x(seq, bounds, self.options.fastp.poly_x_min_len);
            metrics.poly_x_removed[mate] += (end - bounds.end) as u64;
        }
        let maximum = if mate == 0 {
            self.options.fastp.max_len1
        } else {
            self.options.fastp.max_len2.or(self.options.fastp.max_len1)
        };
        if let Some(maximum) = maximum {
            let end = bounds.end;
            bounds.end = bounds.start + bounds.len().min(maximum);
            metrics.max_length_removed[mate] += (end - bounds.end) as u64;
        }
        bounds
    }
    pub(crate) fn bounds(
        &self,
        sample: u32,
        mate: usize,
        seq: &[u8],
        qual: &[u8],
        start: usize,
        metrics: &mut PreprocessingMetrics,
    ) -> BiologicalBounds {
        let (mut bounds, _) = self.prepare(mate, seq, qual, start, metrics);
        if self.adapters.signatures(sample, mate).next().is_some() {
            metrics.scanned += 1;
            let biological = &seq[bounds.range()];
            let full = self
                .adapters
                .signatures(sample, mate)
                .filter_map(|sig| detect_pass_through(biological, sig))
                .min();
            if let Some(relative) = full {
                let end = bounds.start + relative;
                metrics.full[mate] += 1;
                metrics.full_bases[mate] += (bounds.end - end) as u64;
                bounds.end = end;
                metrics.dimers[mate] += u64::from(bounds.len() <= self.options.fastp.dimer_max_len);
            } else {
                let partial = self
                    .adapters
                    .signatures(sample, mate)
                    .filter_map(|sig| {
                        sig.detect_partial(biological, &mut metrics.scratch)
                            .map(|start| (sig.partial_evidence(biological.len(), start), start))
                    })
                    .max()
                    .map(|(_, start)| start);
                if let Some(relative) = partial {
                    let end = bounds.start + relative;
                    metrics.partial[mate] += 1;
                    metrics.partial_bases[mate] += (bounds.end - end) as u64;
                    bounds.end = end;
                    metrics.dimers[mate] +=
                        u64::from(bounds.len() <= self.options.fastp.dimer_max_len);
                }
            }
        }
        self.finish(seq, bounds, mate, metrics)
    }

    /// Mate confirmation is applied once per pair, before counters and filtering.
    pub(crate) fn pair_bounds(
        &self,
        sample: u32,
        sequences: [&[u8]; 2],
        qualities: [&[u8]; 2],
        starts: [usize; 2],
        metrics: &mut PreprocessingMetrics,
    ) -> [BiologicalBounds; 2] {
        metrics.cached_overlap = None;
        if self.options.no_overlap {
            return [
                self.bounds(sample, 0, sequences[0], qualities[0], starts[0], metrics),
                self.bounds(sample, 1, sequences[1], qualities[1], starts[1], metrics),
            ];
        }
        let (a, front1) = self.prepare(0, sequences[0], qualities[0], starts[0], metrics);
        let (b, front2) = self.prepare(1, sequences[1], qualities[1], starts[1], metrics);
        let mut bounds = [a, b];
        let biological = [&sequences[0][a.range()], &sequences[1][b.range()]];
        let mut full = [None, None];
        let mut partial = [None, None];
        for mate in 0..2 {
            if self.adapters.signatures(sample, mate).next().is_some() {
                metrics.scanned += 1;
                full[mate] = self
                    .adapters
                    .signatures(sample, mate)
                    .filter_map(|sig| detect_pass_through(biological[mate], sig))
                    .min();
                if full[mate].is_none() {
                    partial[mate] = self
                        .adapters
                        .signatures(sample, mate)
                        .filter_map(|sig| {
                            sig.detect_partial(biological[mate], &mut metrics.scratch)
                                .map(|start| {
                                    (sig.partial_evidence(biological[mate].len(), start), start)
                                })
                        })
                        .max()
                        .map(|(_, start)| start);
                }
            }
        }
        let candidate = full.iter().chain(&partial).any(Option::is_some);
        let overlap_only = self.adapters.enabled && !self.options.fastp.require_adapter_for_overlap;
        if candidate || overlap_only {
            let overlap = overlap::analyze(
                biological[0],
                biological[1],
                self.options.overlap_len_require,
                self.options.overlap_diff_limit,
                self.options.overlap_diff_percent_limit,
                &mut metrics.overlap_scratch,
            );
            metrics.cached_overlap = Some(([a, b], overlap));
            let mut confirmed = false;
            let mut gapped = false;
            let ends = if let OverlapResult::Unique(overlap) = overlap {
                overlap
                    .insert_length
                    .checked_add(front2)
                    .zip(overlap.insert_length.checked_add(front1))
                    .map(|(a, b)| [a, b])
            } else if overlap == OverlapResult::None
                && self.options.fastp.allow_gap_overlap_trimming
            {
                match overlap::analyze_gap(
                    biological,
                    [front1, front2],
                    self.options.overlap_len_require,
                    self.options.overlap_diff_limit,
                    self.options.overlap_diff_percent_limit,
                    &mut metrics.overlap_scratch,
                ) {
                    Ok(ends) => {
                        gapped = ends.is_some();
                        ends
                    }
                    Err(()) => {
                        metrics.overlap_ambiguous += 1;
                        None
                    }
                }
            } else {
                None
            };
            if let Some(ends) = ends {
                for mate in 0..2 {
                    let start = ends[mate];
                    let full_match = full[mate].is_some()
                        && self
                            .adapters
                            .signatures(sample, mate)
                            .any(|signature| signature.full_at(biological[mate], start));
                    let partial_match = full[mate].is_none()
                        && self
                            .adapters
                            .signatures(sample, mate)
                            .any(|signature| signature.partial_at(biological[mate], start));
                    let inferred = full[mate].is_none()
                        && partial[mate].is_none()
                        && overlap_only
                        && start < biological[mate].len();
                    if full_match || partial_match || inferred {
                        let end = bounds[mate].start + start;
                        let removed = (bounds[mate].end - end) as u64;
                        bounds[mate].end = end;
                        metrics.dimers[mate] +=
                            u64::from(bounds[mate].len() <= self.options.fastp.dimer_max_len);
                        if full_match {
                            metrics.full[mate] += 1;
                            metrics.full_bases[mate] += removed;
                        } else if partial_match {
                            metrics.partial[mate] += 1;
                            metrics.partial_bases[mate] += removed;
                        } else {
                            metrics.overlap_only[mate] += 1;
                            metrics.overlap_only_bases[mate] += removed;
                        }
                        confirmed = true;
                    }
                }
            }
            if confirmed {
                metrics.overlap_confirmed += 1;
                metrics.gapped_pairs += u64::from(gapped);
            } else {
                metrics.overlap_rejected += 1;
            }
            if overlap == OverlapResult::Ambiguous {
                metrics.overlap_ambiguous += 1;
            }
        }
        [
            self.finish(sequences[0], bounds[0], 0, metrics),
            self.finish(sequences[1], bounds[1], 1, metrics),
        ]
    }

    pub(crate) fn feature_overlap(
        &self,
        sequences: [&[u8]; 2],
        bounds: [BiologicalBounds; 2],
        metrics: &mut PreprocessingMetrics,
    ) -> OverlapResult {
        if let Some((cached_bounds, result)) = metrics.cached_overlap
            && cached_bounds == bounds
        {
            return result;
        }
        overlap::analyze(
            &sequences[0][bounds[0].range()],
            &sequences[1][bounds[1].range()],
            self.options.overlap_len_require,
            self.options.overlap_diff_limit,
            self.options.overlap_diff_percent_limit,
            &mut metrics.overlap_scratch,
        )
    }
    /// Pair-wide decision: collect all failed reasons once per fragment, and
    /// update post-filter QC only when every mate passes.
    pub(crate) fn retain(
        &self,
        mates: &[(&[u8], &[u8], BiologicalBounds)],
        assigned: bool,
        metrics: &mut PreprocessingMetrics,
    ) -> Result<bool, String> {
        metrics.processed += 1;
        if self.options.fastp.verbose && metrics.processed.is_multiple_of(1_000_000) {
            eprintln!("Fragments processed: {}", metrics.processed);
        }
        let o = &self.options;
        let scan = self.qc
            || o.max_n.is_some()
            || o.max_low_quality_percent.is_some()
            || o.min_average_quality.is_some();
        if mates.iter().any(|(seq, qual, _)| seq.len() != qual.len()) {
            return Err("Sequence and quality lengths do not match".into());
        }
        if !scan
            && o.min_length.is_none()
            && !o.fastp.low_complexity_filter
            && o.fastp.length_limit.is_none()
            && !metrics.pending_duplicate
            && !metrics.pending_index
        {
            metrics.retained += 1;
            return Ok(true);
        }
        metrics.last_failed_mates = [false; 2];
        let mut failed = [false; 8];
        failed[7] = metrics.pending_index;
        failed[6] = metrics.pending_duplicate;
        let mut observations = [FastqStats::default(), FastqStats::default()];
        for (mate, &(seq, qual, bounds)) in mates.iter().enumerate() {
            let seq = &seq[bounds.range()];
            let qual = &qual[bounds.range()];
            let stats = &mut observations[mate];
            let low = if scan {
                stats.update_with_low(seq, qual, o.qualified_quality.unwrap_or(15))?
            } else {
                0
            };
            let before = failed;
            let mut mate_failed = [false; 8];
            mate_failed[7] = metrics.pending_index;
            mate_failed[0] = o.min_length.is_some_and(|min| seq.len() < min);
            mate_failed[1] = o.max_n.is_some_and(|max| stats.n_bases > max as u64);
            mate_failed[2] = o
                .max_low_quality_percent
                .is_some_and(|max| !seq.is_empty() && low as f64 * 100.0 > max * seq.len() as f64);
            mate_failed[3] = o
                .min_average_quality
                .is_some_and(|min| stats.mean_quality() < min);
            mate_failed[4] = o.fastp.length_limit.is_some_and(|max| seq.len() > max);
            mate_failed[5] = o.fastp.low_complexity_filter
                && (seq.len() < 2
                    || seq
                        .windows(2)
                        .filter(|w| !w[0].eq_ignore_ascii_case(&w[1]))
                        .count() as f64
                        * 100.0
                        < o.fastp.complexity_threshold * (seq.len() - 1) as f64);
            mate_failed[6] = metrics.pending_duplicate;
            metrics.last_failed_mates[mate] = mate_failed.iter().any(|&failed| failed);
            for i in 0..8 {
                failed[i] = before[i] || mate_failed[i];
            }
        }
        let rejected = failed.iter().any(|&f| f);
        if rejected {
            metrics.rejected += 1;
            metrics.assigned_rejected += u64::from(assigned);
            for (count, failed) in metrics.reasons.iter_mut().zip(failed) {
                *count += u64::from(failed);
            }
        } else {
            metrics.retained += 1;
            if self.qc && metrics.pending_merged {
                metrics.merged_stats.merge(&observations[0]);
                if o.fastp.json.is_some() || o.fastp.html.is_some() {
                    metrics.merged_detailed.update(
                        &mates[0].0[mates[0].2.range()],
                        &mates[0].1[mates[0].2.range()],
                    );
                }
            } else if self.qc {
                if o.fastp.json.is_some() || o.fastp.html.is_some() {
                    for (mate, &(seq, qual, bounds)) in mates.iter().enumerate() {
                        metrics.detailed[mate].update(&seq[bounds.range()], &qual[bounds.range()]);
                    }
                }
                for (total, observation) in metrics
                    .after
                    .iter_mut()
                    .zip(&observations)
                    .take(mates.len())
                {
                    total.merge(observation);
                }
            }
        }
        if !rejected && o.fastp.overrepresentation_analysis {
            for (mate, &(seq, _, bounds)) in mates.iter().enumerate() {
                let overrep = if metrics.pending_merged {
                    &mut metrics.merged_overrep
                } else {
                    &mut metrics.overrep[mate]
                };
                overrep
                    .get_or_insert_with(|| {
                        crate::overrepresentation::Overrepresentation::new(
                            o.fastp.overrepresentation_sampling,
                        )
                    })
                    .update(&seq[bounds.range()]);
            }
        }
        Ok(!rejected)
    }
}

impl PreprocessingMetrics {
    pub(crate) fn set_duplicate(&mut self, duplicate: bool, dedup: bool) {
        self.duplicate_fragments += u64::from(duplicate);
        self.pending_duplicate = duplicate && dedup;
        self.pending_merged = false;
        self.cached_overlap = None;
        self.pending_index = false;
        self.last_failed_mates = [false; 2];
    }
    pub(crate) fn mark_merged(&mut self) {
        self.pending_merged = true;
    }
    pub(crate) fn failed_mates(&self) -> [bool; 2] {
        self.last_failed_mates
    }
    pub(crate) fn merge(&mut self, other: &Self) {
        self.corrected_bases += other.corrected_bases;
        self.merged_fragments += other.merged_fragments;
        self.duplicate_fragments += other.duplicate_fragments;
        self.merged_stats.merge(&other.merged_stats);
        self.merged_detailed.merge(&other.merged_detailed);
        if let Some(source) = other.merged_overrep.as_ref() {
            self.merged_overrep
                .get_or_insert_with(|| crate::overrepresentation::Overrepresentation::new(1))
                .merge(source);
        }
        for mate in 0..2 {
            self.detailed[mate].merge(&other.detailed[mate]);
        }
        self.scanned += other.scanned;
        self.overlap_confirmed += other.overlap_confirmed;
        self.gapped_pairs += other.gapped_pairs;
        self.overlap_rejected += other.overlap_rejected;
        self.overlap_ambiguous += other.overlap_ambiguous;
        self.processed += other.processed;
        self.retained += other.retained;
        self.rejected += other.rejected;
        self.assigned_rejected += other.assigned_rejected;
        for (a, b) in self.reasons.iter_mut().zip(other.reasons) {
            *a += b;
        }
        for mate in 0..2 {
            self.overlap_only[mate] += other.overlap_only[mate];
            self.overlap_only_bases[mate] += other.overlap_only_bases[mate];
            self.dimers[mate] += other.dimers[mate];
            self.full[mate] += other.full[mate];
            self.partial[mate] += other.partial[mate];
            self.full_bases[mate] += other.full_bases[mate];
            self.partial_bases[mate] += other.partial_bases[mate];
            self.proximal[mate] += other.proximal[mate];
            self.quality_removed[mate] += other.quality_removed[mate];
            self.poly_g_removed[mate] += other.poly_g_removed[mate];
            self.poly_x_removed[mate] += other.poly_x_removed[mate];
            self.max_length_removed[mate] += other.max_length_removed[mate];
            self.after[mate].merge(&other.after[mate]);
            if let Some(source) = other.overrep[mate].as_ref() {
                self.overrep[mate]
                    .get_or_insert_with(|| crate::overrepresentation::Overrepresentation::new(1))
                    .merge(source);
            }
        }
    }
    pub(crate) fn report_json(&self) -> serde_json::Value {
        use serde_json::json;
        let mut report = json!({"processing":{"fragments_processed":self.processed,"fragments_retained":self.retained,"fragments_rejected":self.rejected,"merged_fragments":self.merged_fragments,"corrected_bases":self.corrected_bases},
            "duplication":{"duplicate_fragments":self.duplicate_fragments,"rate":if self.processed==0{0.0}else{self.duplicate_fragments as f64/self.processed as f64},"method":"exact original-read sequence catalog with explicit memory ceiling"},
            "trimming":{"R1_full_adapter_reads":self.full[0],"R2_full_adapter_reads":self.full[1],"R1_partial_adapter_reads":self.partial[0],"R2_partial_adapter_reads":self.partial[1],"R1_full_adapter_bases":self.full_bases[0],"R2_full_adapter_bases":self.full_bases[1],"overlap_confirmed_pairs":self.overlap_confirmed,"overlap_unconfirmed_pairs":self.overlap_rejected,"overlap_ambiguous_pairs":self.overlap_ambiguous},
            "filtering_result":{"too_short":self.reasons[0],"too_many_N":self.reasons[1],"low_quality_percent":self.reasons[2],"low_mean_quality":self.reasons[3],"too_long":self.reasons[4],"low_complexity":self.reasons[5],"duplicates":self.reasons[6],"index_blacklist":self.reasons[7]}});
        for mate in 0..2 {
            let mut row = self.detailed[mate].json();
            row["summary"] = crate::report::summary(&self.after[mate]);
            report[format!("read{}_after_filtering", mate + 1)] = row;
        }
        if self.merged_stats.reads > 0 {
            let mut row = self.merged_detailed.json();
            row["summary"] = crate::report::summary(&self.merged_stats);
            report["merged_after_filtering"] = row;
        }
        report["trimming"]["partial_adapter_bases_removed"] = json!(self.partial_bases);
        report["trimming"]["proximal_bases_removed"] = json!(self.proximal);
        report["trimming"]["quality_and_fixed_bases_removed"] = json!(self.quality_removed);
        report["trimming"]["poly_g_bases_removed"] = json!(self.poly_g_removed);
        report["trimming"]["poly_x_bases_removed"] = json!(self.poly_x_removed);
        report["trimming"]["maximum_length_bases_removed"] = json!(self.max_length_removed);
        report["overrepresented_sequences"] = json!({"R1_after":self.overrep[0].as_ref().map(|v|v.json()),"R2_after":self.overrep[1].as_ref().map(|v|v.json()),"Merged_after":self.merged_overrep.as_ref().map(|v|v.json())});
        report
    }
    pub(crate) fn write(&self, dir: &Path, paired: bool, qc: bool) -> Result<(), String> {
        let mut file = File::create(dir.join("preprocessing.tsv")).map_err(|e| e.to_string())?;
        writeln!(file, "Metric\tCount").map_err(|e| e.to_string())?;
        for (key, value) in [
            ("duplicate_fragments", self.duplicate_fragments),
            ("merged_fragments", self.merged_fragments),
            ("corrected_bases", self.corrected_bases),
            ("overlap_confirmed_pairs", self.overlap_confirmed),
            ("gapped_overlap_pairs", self.gapped_pairs),
            ("overlap_unconfirmed_pairs", self.overlap_rejected),
            ("overlap_ambiguous_pairs", self.overlap_ambiguous),
            ("fragments_processed", self.processed),
            ("fragments_retained", self.retained),
            ("fragments_rejected", self.rejected),
            ("rejected_duplicates", self.reasons[6]),
            ("rejected_index", self.reasons[7]),
            ("rejected_min_length", self.reasons[0]),
            ("rejected_max_length", self.reasons[4]),
            ("rejected_low_complexity", self.reasons[5]),
            ("rejected_max_n", self.reasons[1]),
            ("rejected_low_quality_percent", self.reasons[2]),
            ("rejected_mean_quality", self.reasons[3]),
        ] {
            writeln!(file, "{key}\t{value}").map_err(|e| e.to_string())?;
        }
        for mate in 0..if paired { 2 } else { 1 } {
            for (key, value) in [
                ("overlap_only_reads", self.overlap_only[mate]),
                ("overlap_only_bases_removed", self.overlap_only_bases[mate]),
                ("adapter_dimers", self.dimers[mate]),
                ("full_adapter_reads", self.full[mate]),
                ("partial_adapter_reads", self.partial[mate]),
                ("full_adapter_bases_removed", self.full_bases[mate]),
                ("partial_adapter_bases_removed", self.partial_bases[mate]),
                ("proximal_bases_removed", self.proximal[mate]),
                ("quality_global_bases_removed", self.quality_removed[mate]),
                ("poly_g_bases_removed", self.poly_g_removed[mate]),
                ("poly_x_bases_removed", self.poly_x_removed[mate]),
                ("max_length_bases_removed", self.max_length_removed[mate]),
            ] {
                writeln!(file, "R{}_{key}\t{value}", mate + 1).map_err(|e| e.to_string())?;
            }
        }
        if qc {
            let mut file =
                File::create(dir.join("retained_fastq_stats.tsv")).map_err(|e| e.to_string())?;
            writeln!(file,"Mate\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent").map_err(|e| e.to_string())?;
            if self.merged_stats.reads > 0 {
                crate::demux::write_stats_row(&mut file, "Merged", &self.merged_stats)?;
            }
            for mate in 0..if paired { 2 } else { 1 } {
                crate::demux::write_stats_row(
                    &mut file,
                    &format!("R{}", mate + 1),
                    &self.after[mate],
                )?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn retained(
        options: PreprocessingArgs,
        seq: &[u8],
        qual: &[u8],
    ) -> (bool, PreprocessingMetrics) {
        let preprocessor = Preprocessor {
            options,
            ..Default::default()
        };
        let mut metrics = PreprocessingMetrics::default();
        let result = preprocessor
            .retain(
                &[(seq, qual, BiologicalBounds::full(seq.len()))],
                true,
                &mut metrics,
            )
            .unwrap();
        (result, metrics)
    }
    #[test]
    fn zero_length_and_disabled_filter_semantics() {
        assert!(retained(PreprocessingArgs::default(), b"", b"").0);
        assert!(
            retained(
                PreprocessingArgs {
                    max_n: Some(0),
                    max_low_quality_percent: Some(0.0),
                    min_length: Some(0),
                    min_average_quality: Some(0.0),
                    ..Default::default()
                },
                b"",
                b""
            )
            .0
        );
        let (result, metrics) = retained(
            PreprocessingArgs {
                min_length: Some(1),
                min_average_quality: Some(1.0),
                ..Default::default()
            },
            b"",
            b"",
        );
        assert!(!result);
        assert_eq!(metrics.rejected, 1);
        assert_eq!(metrics.assigned_rejected, 1);
        assert_eq!(metrics.reasons, [1, 0, 0, 1, 0, 0, 0, 0]);
        assert!(retained(PreprocessingArgs::default(), b"NNNN", b"!!!!").0);
    }
    #[test]
    fn paired_reasons_count_once_and_no_post_qc_for_rejected_pair() {
        let preprocessor = Preprocessor {
            options: PreprocessingArgs {
                max_n: Some(0),
                min_length: Some(4),
                ..Default::default()
            },
            qc: true,
            ..Default::default()
        };
        let mut metrics = PreprocessingMetrics::default();
        assert!(
            !preprocessor
                .retain(
                    &[
                        (b"NN", b"II", BiologicalBounds::full(2)),
                        (b"NN", b"II", BiologicalBounds::full(2))
                    ],
                    true,
                    &mut metrics
                )
                .unwrap()
        );
        assert_eq!(metrics.processed, 1);
        assert_eq!(metrics.rejected, 1);
        assert_eq!(metrics.reasons, [1, 1, 0, 0, 0, 0, 0, 0]);
        assert_eq!(metrics.after[0].reads, 0);
        assert_eq!(metrics.after[1].reads, 0);
    }
}
