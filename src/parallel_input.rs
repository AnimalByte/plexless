use std::path::Path;

use rapidgzip_core::{Decoder, DecoderPool};

use crate::fastq::InputReader;
use crate::input::is_gzip_file;

#[derive(Clone)]
pub(crate) struct ParallelInput {
    decoder: Option<Decoder>,
    pool: Option<DecoderPool>,
    requested_workers: usize,
}

impl ParallelInput {
    pub(crate) fn new(
        max_input_threads: usize,
        initial_input_threads: usize,
        parallel_gzip: bool,
    ) -> Result<Self, String> {
        if !parallel_gzip {
            return Ok(Self {
                decoder: None,
                pool: None,
                requested_workers: 1,
            });
        }

        if max_input_threads < 2 {
            return Err("Parallel gzip requires at least 2 decoder slots".into());
        }

        if initial_input_threads == 0 || initial_input_threads > max_input_threads {
            return Err(format!(
                "Invalid initial parallel gzip allocation: {initial_input_threads} of \
                 {max_input_threads} slots"
            ));
        }

        let pool = DecoderPool::builder()
            .workers(max_input_threads)
            .initial_worker_limit(initial_input_threads)
            .build()
            .map_err(|e| format!("Could not create parallel gzip decoder pool: {e}"))?;

        let decoder = Decoder::builder()
            .decoder_threads(max_input_threads)
            .decoder_pool(pool.clone())
            .build()
            .map_err(|e| format!("Could not configure parallel gzip decoder: {e}"))?;

        Ok(Self {
            decoder: Some(decoder),
            pool: Some(pool),
            requested_workers: max_input_threads,
        })
    }

    pub(crate) fn set_worker_limit(&self, workers: usize) -> Result<(), String> {
        let Some(pool) = self.pool.as_ref() else {
            return Ok(());
        };

        pool.set_worker_limit(workers)
            .map_err(|e| format!("Could not update parallel gzip worker limit: {e}"))
    }

    pub(crate) fn open(&self, path: &Path) -> Result<InputReader, String> {
        let Some(decoder) = self.decoder.as_ref() else {
            return InputReader::try_open(path);
        };

        if !is_gzip_file(path)? {
            return InputReader::try_open(path);
        }

        let reader = decoder
            .open(path)
            .map_err(|e| format!("Could not open gzip input '{}': {e}", path.display()))?;

        reader
            .request_workers(self.requested_workers)
            .map_err(|e| {
                format!(
                    "Could not request parallel gzip workers for '{}': {e}",
                    path.display()
                )
            })?;

        InputReader::from_reader(reader, path)
    }
}
