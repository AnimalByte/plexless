//! Standalone machine-readable and interactive HTML QC reports.
use crate::fastp_options::FastpOptions;
use crate::preprocessing::PreprocessingMetrics;
use crate::stats::{FastqStats, MateQcStats};
use serde_json::{Value, json};
use std::io::Write;
pub(crate) fn summary(s: &FastqStats) -> Value {
    json!({"total_reads":s.reads,"total_bases":s.bases,"q20_bases":s.q20_bases,"q30_bases":s.q30_bases,"q20_rate":s.q20_percent()/100.0,"q30_rate":s.q30_percent()/100.0,"gc_content":s.gc_percent()/100.0,"n_content":s.n_percent()/100.0,"mean_length":s.mean_length(),"min_length":s.min_length.unwrap_or(0),"max_length":s.max_length,"mean_quality":s.mean_quality()})
}
fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
fn write_file(path: &std::path::Path, bytes: &[u8]) -> Result<(), String> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| format!("Could not create report '{}': {e}", path.display()))?;
    file.write_all(bytes)
        .and_then(|_| file.flush())
        .map_err(|e| e.to_string())
}
pub(crate) fn write(
    o: &FastpOptions,
    before: [Option<&MateQcStats>; 2],
    metrics: &PreprocessingMetrics,
) -> Result<(), String> {
    if o.json.is_none() && o.html.is_none() {
        return Ok(());
    }
    let mut document = metrics.report_json();
    document["duplication"]["method"] = json!(if o.dedup {
        "exact original-read sequences, bounded memory with temporary SQLite spill"
    } else {
        "bounded four-hash Bloom estimate; false positives affect reporting only"
    });
    document["duplication"]["enabled"] =
        json!(o.dedup || o.evaluate_duplication || (o.fastp_defaults && !o.dont_eval_duplication));
    document["trimming"]["overlap_only_reads"] = json!(metrics.overlap_only);
    document["trimming"]["overlap_only_bases_removed"] = json!(metrics.overlap_only_bases);
    document["trimming"]["gapped_pairs"] = json!(metrics.gapped_pairs);
    document["trimming"]["adapter_dimers"] = json!(metrics.dimers);
    document["report_limits"] = json!({"overrepresentation":"bounded 4096-candidate sampled Space-Saving occurrence counts, at most 100 reported sequences per mate; full reads capped at 512 bases","alignment":"ungapped unique overlap for merge and correction; optional single gap for short-insert adapter clipping only","defaults":"Plexless retention defaults are preserved unless --fastp-defaults is selected"});
    document["schema_version"] = json!(1);
    document["title"] = json!(o.report_title);
    document["tool"] = json!({"name":"plexless","version":env!("CARGO_PKG_VERSION"),"fastp_reference":"ee125a52fddae9d8189e13cd8a33ad04988a1351"});
    document["input_coordinates"] =
        json!("biological suffix after declared proximal structure, before preprocessing");
    for (mate, stats) in before.into_iter().enumerate() {
        if let Some(stats) = stats {
            let mut row = stats.distribution().json();
            row["summary"] = summary(stats.biological());
            document["overrepresented_sequences"][format!("R{}_before", mate + 1)] =
                stats.overrepresentation();
            document[format!("read{}_before_filtering", mate + 1)] = row;
        }
    }
    if let Some(path) = &o.json {
        let bytes = serde_json::to_vec_pretty(&document).map_err(|e| e.to_string())?;
        write_file(path, &bytes)?;
    }
    if let Some(path) = &o.html {
        let data = serde_json::to_string(&document)
            .map_err(|e| e.to_string())?
            .replace('<', "\\u003c")
            .replace('>', "\\u003e")
            .replace('&', "\\u0026");
        let template = include_str!("report_template.html");
        let page = template
            .replace("{{TITLE}}", &escape(&o.report_title))
            .replace("{{DATA}}", &data);
        write_file(path, page.as_bytes())?;
    }
    Ok(())
}
