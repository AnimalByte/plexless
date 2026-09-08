# Routing architecture

## Imported simplex data flow

The imported implementation already had a sound streaming pipeline:

```text
FASTQ input
  -> Needletail parsing / parallel gzip decoding
  -> bounded record batches
  -> demultiplexing workers
  -> per-batch FASTQ buffers and gzip compression
  -> ordered batch emission
  -> bounded output-writer cache
```

Paired input validates normalized record IDs and uses bounded lookahead to
resynchronize, emitting records present in only one input as orphans. Assigned
reads have each mate's declared structured prefix trimmed; unassigned and
orphan reads remain untrimmed. Those parts of the pipeline are unchanged.

The imported matching path was flat:

```text
scan structures and assemble every symbol
  -> decode each symbol against its complete global whitelist
  -> collect decoder IDs in a fixed [DecodeResult; 3]
  -> allocate Vec<u32>
  -> HashMap<Vec<u32>, sample_id>
```

That made unrelated child namespaces participate in correction validation and
matching.

## Compiled hierarchical model

`ReadLayout::compile_extraction_plans` converts parsed read structures into one
`BarcodeExtractionPlan` per logical symbol. Each plan contains its logical
length and an ordered list of `(mate, start, end, orientation)` pieces. The
piece order is all R1 segments in structure order followed by all R2 segments
in structure order.

`RoutingTree::new` compiles the ordinary tabular sample sheet into immutable
contiguous `CompiledNode` values. Each node stores:

- its generic barcode symbol;
- the index of its extraction plan;
- a `Decoder` built only from barcode candidates reachable at that parent
  path;
- a target array indexed directly by the decoder-local barcode ID.

A target is another node, a sample ID, or an unrouted terminal. The root uses
the complete root whitelist so a known but unused root barcode remains
unrouted. Every deeper node contains only sample-sheet children beneath its
specific parent path.

## Read hot path

```text
check R1/R2 structured-prefix lengths
  -> assemble root symbol into a 32-byte stack buffer
  -> 2-bit encode and decode in the root node
  -> index the resolved target
  -> if it is a node, lazily assemble that node's symbol
  -> repeat until sample, unmatched, ambiguous, or unrouted
```

Traversal allocates no strings, vectors, or maps per read. An exact/corrected
decoder result is a node-local integer used directly to index the target
array. Invalid observed bases classify the read as unmatched. Observed `N`
bases retain the existing mismatch-budget behavior, but the N-aware scan sees
only the current node's candidates.

## Node-local correction safety

Each node independently validates minimum sibling Hamming distance
`2 * max_mismatches + 1` and independently precomputes exact and correction
lookups. An error identifies the parent path plus both conflicting IDs.

Consequently, close or identical child sequences are allowed in independent
parent namespaces. They conflict only when both are candidates in the same
node. This is the defining difference from the imported global-decoder model.

## Extension point

`Decoder` still exposes exact, corrected, ambiguous, and unmatched outcomes,
while `RoutingTree` owns traversal and node targets. A future path-aware rescue
implementation can maintain multiple `(node, score)` states over the same
compiled representation instead of replacing the extraction or tree model.
