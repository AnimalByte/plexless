use std::cmp::Ordering;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::output::OutputMate;
use crate::samples::SampleSheet;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct SampleMetrics {
    pub(crate) assigned_fragments: u64,
    pub(crate) single_records: u64,
    pub(crate) r1_records: u64,
    pub(crate) r2_records: u64,
    pub(crate) assigned_bases_single: u64,
    pub(crate) assigned_bases_r1: u64,
    pub(crate) assigned_bases_r2: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SampleStatus {
    Ok,
    LowRepresentation,
    Missing,
}

impl SampleStatus {
    fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "OK",
            Self::LowRepresentation => "LOW_REPRESENTATION",
            Self::Missing => "MISSING",
        }
    }
}

#[derive(Debug)]
pub(crate) struct SampleQc {
    paired: bool,
    metrics: Vec<SampleMetrics>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SampleQcSummary {
    pub(crate) expected: usize,
    pub(crate) populated: usize,
    pub(crate) missing: usize,
    pub(crate) low: usize,
    pub(crate) median_nonzero: f64,
}

impl SampleQc {
    pub(crate) fn new(sample_count: usize, paired: bool) -> Self {
        Self {
            paired,
            metrics: vec![SampleMetrics::default(); sample_count],
        }
    }

    pub(crate) fn observe_output(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        records: u64,
        bases: u64,
    ) -> Result<(), String> {
        let index = usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?;
        let metrics = self
            .metrics
            .get_mut(index)
            .ok_or_else(|| format!("Unknown sample ID {sample_id}"))?;
        match (self.paired, mate) {
            (false, OutputMate::Single) => {
                metrics.assigned_fragments = metrics
                    .assigned_fragments
                    .checked_add(records)
                    .ok_or("Sample fragment count overflow")?;
                metrics.single_records = metrics
                    .single_records
                    .checked_add(records)
                    .ok_or("Sample record count overflow")?;
                metrics.assigned_bases_single = metrics
                    .assigned_bases_single
                    .checked_add(bases)
                    .ok_or("Sample base count overflow")?;
            }
            (true, OutputMate::R1) => {
                metrics.assigned_fragments = metrics
                    .assigned_fragments
                    .checked_add(records)
                    .ok_or("Sample fragment count overflow")?;
                metrics.r1_records = metrics
                    .r1_records
                    .checked_add(records)
                    .ok_or("Sample R1 count overflow")?;
                metrics.assigned_bases_r1 = metrics
                    .assigned_bases_r1
                    .checked_add(bases)
                    .ok_or("Sample R1 base count overflow")?;
            }
            (true, OutputMate::R2) => {
                metrics.r2_records = metrics
                    .r2_records
                    .checked_add(records)
                    .ok_or("Sample R2 count overflow")?;
                metrics.assigned_bases_r2 = metrics
                    .assigned_bases_r2
                    .checked_add(bases)
                    .ok_or("Sample R2 base count overflow")?;
            }
            (false, _) => return Err("Single-end metrics require OutputMate::Single".into()),
            (true, OutputMate::Single) => {
                return Err("Paired-end metrics require R1 or R2".into());
            }
        }
        Ok(())
    }

    pub(crate) fn merge(&mut self, other: Self) -> Result<(), String> {
        if self.paired != other.paired || self.metrics.len() != other.metrics.len() {
            return Err("Cannot merge incompatible sample QC accumulators".into());
        }
        for (target, source) in self.metrics.iter_mut().zip(other.metrics) {
            target.assigned_fragments = target
                .assigned_fragments
                .checked_add(source.assigned_fragments)
                .ok_or("Sample fragment count overflow")?;
            target.single_records = target
                .single_records
                .checked_add(source.single_records)
                .ok_or("Sample record count overflow")?;
            target.r1_records = target
                .r1_records
                .checked_add(source.r1_records)
                .ok_or("Sample R1 count overflow")?;
            target.r2_records = target
                .r2_records
                .checked_add(source.r2_records)
                .ok_or("Sample R2 count overflow")?;
            target.assigned_bases_single = target
                .assigned_bases_single
                .checked_add(source.assigned_bases_single)
                .ok_or("Sample base count overflow")?;
            target.assigned_bases_r1 = target
                .assigned_bases_r1
                .checked_add(source.assigned_bases_r1)
                .ok_or("Sample R1 base count overflow")?;
            target.assigned_bases_r2 = target
                .assigned_bases_r2
                .checked_add(source.assigned_bases_r2)
                .ok_or("Sample R2 base count overflow")?;
        }
        Ok(())
    }

    pub(crate) fn verify(&self, global_assigned: u64) -> Result<(), String> {
        let sample_total = self.metrics.iter().try_fold(0u64, |total, sample| {
            total
                .checked_add(sample.assigned_fragments)
                .ok_or("Sample fragment total overflow")
        })?;
        if sample_total != global_assigned {
            return Err(format!(
                "Assigned-count reconciliation failed: global={global_assigned}, samples={sample_total}"
            ));
        }
        if self.paired {
            for (index, sample) in self.metrics.iter().enumerate() {
                if sample.assigned_fragments != sample.r1_records
                    || sample.assigned_fragments != sample.r2_records
                {
                    return Err(format!(
                        "Paired sample {index} has inconsistent fragment/R1/R2 counts"
                    ));
                }
            }
        }
        Ok(())
    }

    pub(crate) fn write_report(
        &self,
        output_dir: &Path,
        samples: &SampleSheet,
        low_fraction: f64,
    ) -> Result<SampleQcSummary, String> {
        if self.metrics.len() != samples.samples.len() {
            return Err("Sample metrics do not match the worklist".into());
        }
        let median = median_nonzero(&self.metrics);
        let assigned_total = self.metrics.iter().try_fold(0u64, |total, item| {
            total
                .checked_add(item.assigned_fragments)
                .ok_or("Assigned sample total overflow")
        })?;
        let statuses: Vec<_> = self
            .metrics
            .iter()
            .map(|item| classify(item.assigned_fragments, median, low_fraction))
            .collect();
        let populated = statuses
            .iter()
            .filter(|&&status| status != SampleStatus::Missing)
            .count();
        let missing = statuses
            .iter()
            .filter(|&&status| status == SampleStatus::Missing)
            .count();
        let low = statuses
            .iter()
            .filter(|&&status| status == SampleStatus::LowRepresentation)
            .count();

        let path = output_dir.join("sample_metrics.tsv");
        let mut file = File::create(&path)
            .map_err(|error| format!("Could not create '{}': {error}", path.display()))?;
        writeln!(
            file,
            "sample\tassigned_fragments\tsingle_records\tr1_records\tr2_records\tassigned_bases_single\tassigned_bases_r1\tassigned_bases_r2\tfraction_of_assigned_fragments\tfraction_of_median\tstatus"
        )
        .map_err(|error| format!("Could not write sample metrics: {error}"))?;

        for ((sample, metrics), status) in samples.samples.iter().zip(&self.metrics).zip(&statuses)
        {
            let assigned_fraction = ratio(metrics.assigned_fragments, assigned_total as f64);
            let median_fraction = ratio(metrics.assigned_fragments, median);
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{assigned_fraction:.8}\t{median_fraction:.8}\t{}",
                sample.name,
                metrics.assigned_fragments,
                metrics.single_records,
                metrics.r1_records,
                metrics.r2_records,
                metrics.assigned_bases_single,
                metrics.assigned_bases_r1,
                metrics.assigned_bases_r2,
                status.as_str(),
            )
            .map_err(|error| format!("Could not write sample metrics: {error}"))?;
        }

        file.flush()
            .map_err(|error| format!("Could not flush sample metrics: {error}"))?;

        if missing > 0 {
            eprintln!("WARNING: {missing} expected sample(s) received zero assigned reads:");
            for (sample, status) in samples.samples.iter().zip(&statuses) {
                if *status == SampleStatus::Missing {
                    eprintln!("  {}", sample.name);
                }
            }
        }
        if low > 0 {
            eprintln!(
                "LOW REPRESENTATION: {low} sample(s) are below {:.2}% of the nonzero median ({median:.0} fragments):",
                low_fraction * 100.0
            );
            for ((sample, metrics), status) in
                samples.samples.iter().zip(&self.metrics).zip(&statuses)
            {
                if *status == SampleStatus::LowRepresentation {
                    eprintln!(
                        "  {}: {} fragments ({:.2}% of median)",
                        sample.name,
                        metrics.assigned_fragments,
                        ratio(metrics.assigned_fragments, median) * 100.0
                    );
                }
            }
        }

        Ok(SampleQcSummary {
            expected: self.metrics.len(),
            populated,
            missing,
            low,
            median_nonzero: median,
        })
    }
}

fn ratio(numerator: u64, denominator: f64) -> f64 {
    if denominator == 0.0 {
        0.0
    } else {
        numerator as f64 / denominator
    }
}

fn median_nonzero(metrics: &[SampleMetrics]) -> f64 {
    let mut counts: Vec<u64> = metrics
        .iter()
        .map(|item| item.assigned_fragments)
        .filter(|&count| count > 0)
        .collect();
    counts.sort_unstable();
    match counts.len() {
        0 => 0.0,
        length if length % 2 == 1 => counts[length / 2] as f64,
        length => {
            let left = counts[length / 2 - 1] as f64;
            let right = counts[length / 2] as f64;
            (left + right) / 2.0
        }
    }
}

fn classify(count: u64, median: f64, low_fraction: f64) -> SampleStatus {
    match count.cmp(&0) {
        Ordering::Equal => SampleStatus::Missing,
        Ordering::Greater if (count as f64) < median * low_fraction => {
            SampleStatus::LowRepresentation
        }
        _ => SampleStatus::Ok,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metrics(counts: &[u64]) -> SampleQc {
        let mut qc = SampleQc::new(counts.len(), false);
        for (index, &count) in counts.iter().enumerate() {
            qc.observe_output(index as u32, OutputMate::Single, count, count * 10)
                .unwrap();
        }
        qc
    }

    #[test]
    fn completeness_handles_all_present_missing_all_zero_and_single_sample() {
        assert_eq!(median_nonzero(&metrics(&[10, 20, 30]).metrics), 20.0);
        let several_missing = metrics(&[0, 20, 0]);
        let median = median_nonzero(&several_missing.metrics);
        assert_eq!(median, 20.0);
        assert_eq!(
            several_missing
                .metrics
                .iter()
                .map(|item| classify(item.assigned_fragments, median, 0.05))
                .collect::<Vec<_>>(),
            vec![
                SampleStatus::Missing,
                SampleStatus::Ok,
                SampleStatus::Missing
            ]
        );

        let all_zero = metrics(&[0, 0]);
        let median = median_nonzero(&all_zero.metrics);
        assert_eq!(median, 0.0);
        assert!(all_zero.metrics.iter().all(|item| classify(
            item.assigned_fragments,
            median,
            0.05
        ) == SampleStatus::Missing));

        let single = metrics(&[7]);
        let median = median_nonzero(&single.metrics);
        assert_eq!(median, 7.0);
        assert_eq!(classify(7, median, 0.05), SampleStatus::Ok);
        assert_eq!(classify(0, 0.0, 0.05), SampleStatus::Missing);
    }

    #[test]
    fn low_representation_uses_nonzero_median_and_threshold() {
        let qc = metrics(&[4, 100, 100, 110]);
        let median = median_nonzero(&qc.metrics);
        assert_eq!(median, 100.0);
        assert_eq!(classify(4, median, 0.05), SampleStatus::LowRepresentation);
        assert_eq!(classify(5, median, 0.05), SampleStatus::Ok);
    }

    #[test]
    fn sample_counts_reconcile_with_global_assigned_count() {
        let qc = metrics(&[3, 4, 0]);
        assert!(qc.verify(7).is_ok());
        assert!(qc.verify(8).is_err());
    }

    #[test]
    fn paired_metrics_count_fragments_once_and_require_both_mates() {
        let mut qc = SampleQc::new(1, true);
        qc.observe_output(0, OutputMate::R1, 3, 30).unwrap();
        qc.observe_output(0, OutputMate::R2, 3, 36).unwrap();

        assert!(qc.verify(3).is_ok());
        assert_eq!(qc.metrics[0].assigned_fragments, 3);
        assert_eq!(qc.metrics[0].r1_records, 3);
        assert_eq!(qc.metrics[0].r2_records, 3);
        assert_eq!(qc.metrics[0].assigned_bases_r1, 30);
        assert_eq!(qc.metrics[0].assigned_bases_r2, 36);

        qc.observe_output(0, OutputMate::R1, 1, 10).unwrap();
        assert!(qc.verify(4).is_err());
    }
}
