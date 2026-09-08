use std::collections::BTreeMap;

use crate::barcodes::BarcodeCatalog;
use crate::decoder::{DecodeResult, Decoder};
use crate::encoder::{MAX_BARCODE_LENGTH, encode};
use crate::samples::SampleSheet;
use crate::structure::{BarcodeExtractionPlans, ReadLayout};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RouteResult {
    Assigned { sample_id: u32 },
    Unmatched,
    Ambiguous,
    Unrouted,
}

#[derive(Debug, Clone, Copy)]
enum RouteTarget {
    Node(u32),
    Sample(u32),
    Unrouted,
}

/// One routing decision point. Decoder-local IDs index `targets` directly, so
/// successful traversal requires no path vector or route hash lookup.
#[derive(Debug)]
struct CompiledNode {
    symbol: u8,
    plan_index: usize,
    decoder: Decoder,
    targets: Vec<RouteTarget>,
}

/// Immutable startup-compiled hierarchical router shared by all workers.
#[derive(Debug)]
pub struct RoutingTree {
    nodes: Vec<CompiledNode>,
    root: u32,
    extraction: BarcodeExtractionPlans,
}

#[derive(Debug, Default)]
struct BuildNode {
    children: BTreeMap<u32, BuildTarget>,
}

#[derive(Debug)]
enum BuildTarget {
    Node(Box<BuildNode>),
    Sample(u32),
}

impl RoutingTree {
    pub fn new(
        layout: &ReadLayout,
        catalog: &BarcodeCatalog,
        samples: &SampleSheet,
        max_mismatches: u8,
    ) -> Result<Self, String> {
        let symbols = layout.barcode_symbols();
        if symbols.is_empty() {
            return Err("Read structure contains no sample barcode segments".into());
        }

        let extraction = layout.compile_extraction_plans()?;
        let mut build_root = BuildNode::default();

        for (index, sample) in samples.samples.iter().enumerate() {
            let sample_id = u32::try_from(index).map_err(|_| "Too many samples")?;
            if sample.barcode_ids.len() != symbols.len() {
                return Err(format!(
                    "Sample '{}' has {} routing levels, but the read structure requires {}",
                    sample.name,
                    sample.barcode_ids.len(),
                    symbols.len()
                ));
            }
            insert_sample_path(
                &mut build_root,
                &sample.barcode_ids,
                sample_id,
                &sample.name,
            )?;
        }

        let mut nodes = Vec::new();
        let mut parent_path = Vec::new();
        let root = compile_node(
            &build_root,
            0,
            &symbols,
            &extraction,
            catalog,
            max_mismatches,
            &mut parent_path,
            &mut nodes,
        )?;

        Ok(Self {
            nodes,
            root,
            extraction,
        })
    }

    /// Traverses one route lazily. Each level is assembled and decoded only
    /// after its parent has resolved to a child node.
    ///
    /// `None` means a read is shorter than its compiled structured prefix.
    pub fn route_read(
        &self,
        r1_seq: &[u8],
        r2_seq: Option<&[u8]>,
    ) -> Result<Option<RouteResult>, String> {
        if !self.extraction.reads_are_long_enough(r1_seq, r2_seq)? {
            return Ok(None);
        }

        let mut buffer = [0u8; MAX_BARCODE_LENGTH];
        let mut node_index = self.root;

        loop {
            let node = self
                .nodes
                .get(usize::try_from(node_index).expect("node index fits usize"))
                .ok_or("Internal error: invalid routing node")?;
            let Some(length) =
                self.extraction
                    .assemble(node.plan_index, r1_seq, r2_seq, &mut buffer)?
            else {
                return Ok(None);
            };

            let observed = match encode(&buffer[..length]) {
                Ok(observed) => observed,
                Err(_) => return Ok(Some(RouteResult::Unmatched)),
            };

            let local_barcode_id = match node.decoder.decode(observed) {
                DecodeResult::Exact { barcode_id } | DecodeResult::Corrected { barcode_id, .. } => {
                    barcode_id
                }
                DecodeResult::Ambiguous => return Ok(Some(RouteResult::Ambiguous)),
                DecodeResult::Unmatched => return Ok(Some(RouteResult::Unmatched)),
            };

            let target = node
                .targets
                .get(
                    usize::try_from(local_barcode_id)
                        .expect("decoder-local barcode index fits usize"),
                )
                .ok_or("Internal error: decoder target is missing")?;

            match *target {
                RouteTarget::Node(next) => node_index = next,
                RouteTarget::Sample(sample_id) => {
                    return Ok(Some(RouteResult::Assigned { sample_id }));
                }
                RouteTarget::Unrouted => return Ok(Some(RouteResult::Unrouted)),
            }
        }
    }

    pub fn r1_prefix_len(&self) -> usize {
        self.extraction.r1_prefix_len()
    }

    pub fn r2_prefix_len(&self) -> usize {
        self.extraction.r2_prefix_len()
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub fn root_symbol(&self) -> u8 {
        self.nodes[usize::try_from(self.root).expect("root index fits usize")].symbol
    }
}

fn insert_sample_path(
    root: &mut BuildNode,
    barcode_ids: &[u32],
    sample_id: u32,
    sample_name: &str,
) -> Result<(), String> {
    let mut node = root;

    for (level, &barcode_id) in barcode_ids.iter().enumerate() {
        let is_leaf = level + 1 == barcode_ids.len();

        if is_leaf {
            if node
                .children
                .insert(barcode_id, BuildTarget::Sample(sample_id))
                .is_some()
            {
                return Err(format!(
                    "Complete barcode path for sample '{sample_name}' is already assigned"
                ));
            }
            return Ok(());
        }

        let target = node
            .children
            .entry(barcode_id)
            .or_insert_with(|| BuildTarget::Node(Box::default()));
        match target {
            BuildTarget::Node(child) => node = child,
            BuildTarget::Sample(_) => {
                return Err(format!(
                    "Sample '{sample_name}' has a routing level below an existing complete path"
                ));
            }
        }
    }

    Err(format!("Sample '{sample_name}' has no barcode path"))
}

#[allow(clippy::too_many_arguments)]
fn compile_node(
    build: &BuildNode,
    level: usize,
    symbols: &[u8],
    extraction: &BarcodeExtractionPlans,
    catalog: &BarcodeCatalog,
    max_mismatches: u8,
    parent_path: &mut Vec<(u8, u32)>,
    nodes: &mut Vec<CompiledNode>,
) -> Result<u32, String> {
    let symbol = *symbols
        .get(level)
        .ok_or("Internal error: routing level exceeds sample columns")?;
    let set = catalog.set(symbol).ok_or_else(|| {
        format!(
            "Required barcode symbol {} has no whitelist entries",
            symbol as char
        )
    })?;
    let plan_index = extraction.plan_index(symbol).ok_or_else(|| {
        format!(
            "Required barcode symbol {} has no extraction plan",
            symbol as char
        )
    })?;

    // At the root every whitelist barcode is a valid first-level barcode;
    // unused roots terminate as Unrouted. Below the root, candidates are
    // exactly the children present beneath that parent path.
    let candidate_ids: Vec<u32> = if level == 0 {
        (0..set.barcodes.len())
            .map(|index| u32::try_from(index).map_err(|_| "Too many barcodes"))
            .collect::<Result<_, _>>()?
    } else {
        build.children.keys().copied().collect()
    };

    if candidate_ids.is_empty() {
        return Err(format!(
            "Required barcode symbol {} has no whitelist entries under {}",
            symbol as char,
            format_parent_path(parent_path, catalog)?
        ));
    }

    let mut targets = Vec::with_capacity(candidate_ids.len());
    for &barcode_id in &candidate_ids {
        let target = match build.children.get(&barcode_id) {
            Some(BuildTarget::Sample(sample_id)) => RouteTarget::Sample(*sample_id),
            Some(BuildTarget::Node(child)) => {
                parent_path.push((symbol, barcode_id));
                let child_index = compile_node(
                    child,
                    level + 1,
                    symbols,
                    extraction,
                    catalog,
                    max_mismatches,
                    parent_path,
                    nodes,
                )?;
                parent_path.pop();
                RouteTarget::Node(child_index)
            }
            None => RouteTarget::Unrouted,
        };
        targets.push(target);
    }

    let location = format_parent_path(parent_path, catalog)?;
    let decoder = Decoder::new_subset(set, &candidate_ids, max_mismatches, &location)?;
    let node_index = u32::try_from(nodes.len()).map_err(|_| "Too many routing nodes")?;
    nodes.push(CompiledNode {
        symbol,
        plan_index,
        decoder,
        targets,
    });
    Ok(node_index)
}

fn format_parent_path(
    parent_path: &[(u8, u32)],
    catalog: &BarcodeCatalog,
) -> Result<String, String> {
    if parent_path.is_empty() {
        return Ok("routing root".into());
    }

    let mut result = String::from("path ");
    for (index, &(symbol, barcode_id)) in parent_path.iter().enumerate() {
        if index != 0 {
            result.push_str(" / ");
        }
        let barcode = catalog
            .set(symbol)
            .and_then(|set| set.barcode(barcode_id))
            .ok_or("Internal error: invalid barcode in parent path")?;
        result.push(symbol as char);
        result.push('=');
        result.push_str(&barcode.id);
    }
    Ok(result)
}
