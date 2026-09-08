use std::collections::HashMap;

use crate::decoder::DecodeResult;
use crate::samples::SampleSheet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RouteResult {
    Assigned { sample_id: u32 },
    Unmatched,
    Ambiguous,
    Unrouted,
}

#[derive(Debug)]
pub struct RoutingTable {
    routes: HashMap<Vec<u32>, u32>,
}

impl RoutingTable {
    pub fn new(samples: &SampleSheet) -> Result<Self, String> {
        let mut routes = HashMap::with_capacity(samples.samples.len());

        for (index, sample) in samples.samples.iter().enumerate() {
            let sample_id = u32::try_from(index).map_err(|_| "Too many samples")?;

            if routes
                .insert(sample.barcode_ids.clone(), sample_id)
                .is_some()
            {
                return Err(format!(
                    "Duplicate barcode combination for sample '{}'",
                    sample.name
                ));
            }
        }

        Ok(Self { routes })
    }

    pub fn route(&self, calls: &[DecodeResult]) -> RouteResult {
        let mut barcode_ids = Vec::with_capacity(calls.len());

        for call in calls {
            match call {
                DecodeResult::Exact { barcode_id } | DecodeResult::Corrected { barcode_id, .. } => {
                    barcode_ids.push(*barcode_id);
                }

                DecodeResult::Ambiguous => {
                    return RouteResult::Ambiguous;
                }

                DecodeResult::Unmatched => {
                    return RouteResult::Unmatched;
                }
            }
        }

        match self.routes.get(&barcode_ids) {
            Some(&sample_id) => RouteResult::Assigned { sample_id },

            None => RouteResult::Unrouted,
        }
    }
}
