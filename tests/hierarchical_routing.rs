mod common;

use plexless::barcodes::BarcodeCatalog;
use plexless::routing::{RouteResult, RoutingTree};
use plexless::samples::SampleSheet;
use plexless::structure::ReadLayout;

use common::TestDir;

fn compile_tree(
    test: &TestDir,
    structure: &str,
    barcodes: &str,
    samples: &str,
    max_mismatches: u8,
) -> Result<RoutingTree, String> {
    let layout = ReadLayout::single(Some(structure))?;
    let barcode_path = test.write("barcodes.tsv", barcodes);
    let sample_path = test.write("samples.tsv", samples);
    let catalog = BarcodeCatalog::load(&barcode_path, &layout)?;
    let sheet = SampleSheet::load(&sample_path, &layout, &catalog)?;
    RoutingTree::new(&layout, &catalog, &sheet, max_mismatches)
}

fn route(tree: &RoutingTree, sequence: &str) -> RouteResult {
    tree.route_read(sequence.as_bytes(), None)
        .expect("Routing should not fail")
        .expect("Read should be long enough")
}

#[test]
fn basic_hierarchical_exact_routing() {
    let test = TestDir::new("hierarchical-exact");
    let tree = compile_tree(
        &test,
        "R1_4A4B",
        "Set\tID\tSequence\n\
         A\tA1\tAAAA\n\
         A\tA2\tCCCC\n\
         B\tB1\tGGGG\n\
         B\tB2\tTTTT\n",
        "Sample\tA\tB\n\
         sample1\tA1\tB1\n\
         sample2\tA1\tB2\n\
         sample3\tA2\tB1\n\
         sample4\tA2\tB2\n",
        0,
    )
    .expect("Hierarchy should compile");

    assert_eq!(
        route(&tree, "AAAAGGGG"),
        RouteResult::Assigned { sample_id: 0 }
    );
    assert_eq!(
        route(&tree, "AAAATTTT"),
        RouteResult::Assigned { sample_id: 1 }
    );
    assert_eq!(
        route(&tree, "CCCCGGGG"),
        RouteResult::Assigned { sample_id: 2 }
    );
    assert_eq!(
        route(&tree, "CCCCTTTT"),
        RouteResult::Assigned { sample_id: 3 }
    );
    assert_eq!(tree.node_count(), 3);
}

#[test]
fn child_sequence_reuse_is_allowed_across_parent_namespaces() {
    let test = TestDir::new("child-sequence-reuse");
    let tree = compile_tree(
        &test,
        "R1_4A4B",
        "Set\tID\tSequence\n\
         A\tA1\tAAAA\n\
         A\tA2\tCCCC\n\
         B\tB1\tACTG\n\
         B\tB7\tACTG\n",
        "Sample\tA\tB\n\
         sample1\tA1\tB1\n\
         sample2\tA2\tB7\n",
        0,
    )
    .expect("Reused child sequence should compile in independent namespaces");

    assert_eq!(
        route(&tree, "AAAAACTG"),
        RouteResult::Assigned { sample_id: 0 }
    );
    assert_eq!(
        route(&tree, "CCCCACTG"),
        RouteResult::Assigned { sample_id: 1 }
    );
}

#[test]
fn correction_safety_is_validated_per_parent_node() {
    let test = TestDir::new("parent-local-correction");
    let tree = compile_tree(
        &test,
        "R1_6A6B",
        "Set\tID\tSequence\n\
         A\tA1\tAAAAAA\n\
         A\tA2\tCCCCCC\n\
         B\tB1\tAAAAAA\n\
         B\tB2\tCCCCCC\n\
         B\tB3\tAAAAAT\n\
         B\tB4\tGGGGGG\n",
        "Sample\tA\tB\n\
         sample1\tA1\tB1\n\
         sample2\tA1\tB2\n\
         sample3\tA2\tB3\n\
         sample4\tA2\tB4\n",
        1,
    )
    .expect("Globally close B barcodes in separate parents should be safe");

    assert_eq!(
        route(&tree, "CCCCCCAAAAAG"),
        RouteResult::Assigned { sample_id: 2 }
    );
}

#[test]
fn correction_collision_between_true_siblings_is_rejected_with_path() {
    let test = TestDir::new("sibling-correction-collision");
    let error = compile_tree(
        &test,
        "R1_6A6B",
        "Set\tID\tSequence\n\
         A\tA1\tAAAAAA\n\
         B\tB1\tAAAAAA\n\
         B\tB2\tAAAAAT\n",
        "Sample\tA\tB\n\
         sample1\tA1\tB1\n\
         sample2\tA1\tB2\n",
        1,
    )
    .expect_err("Close siblings must be rejected");

    assert!(error.contains("path A=A1"), "unexpected error: {error}");
    assert!(error.contains("B1") && error.contains("B2"));
    assert!(error.contains("minimum Hamming distance 3"));
}

#[test]
fn three_level_hierarchy_routes_to_leaf_sample() {
    let test = TestDir::new("three-level");
    let tree = compile_tree(
        &test,
        "R1_2A2B2C",
        "Set\tID\tSequence\n\
         A\tA1\tAA\n\
         A\tA2\tCC\n\
         B\tB1\tGG\n\
         B\tB2\tTT\n\
         C\tC1\tAC\n\
         C\tC2\tGT\n\
         C\tC3\tCA\n",
        "Sample\tA\tB\tC\n\
         sample1\tA1\tB1\tC1\n\
         sample2\tA1\tB1\tC2\n\
         sample3\tA1\tB2\tC1\n\
         sample4\tA2\tB1\tC3\n",
        0,
    )
    .expect("Three-level hierarchy should compile");

    assert_eq!(
        route(&tree, "AAGGGT"),
        RouteResult::Assigned { sample_id: 1 }
    );
    assert_eq!(route(&tree, "CCTTAC"), RouteResult::Unmatched);
}

#[test]
fn more_than_three_generic_barcode_levels_are_supported() {
    let test = TestDir::new("generic-levels");
    let tree = compile_tree(
        &test,
        "R1_1A1B1C1D1E1T",
        "Set\tID\tSequence\n\
         A\tA1\tA\n\
         B\tB1\tC\n\
         C\tC1\tG\n\
         D\tD1\tT\n\
         E\tE1\tA\n",
        "Sample\tA\tB\tC\tD\tE\n\
         sample1\tA1\tB1\tC1\tD1\tE1\n",
        0,
    )
    .expect("Five generic levels should compile");

    assert_eq!(
        route(&tree, "ACGTAA"),
        RouteResult::Assigned { sample_id: 0 }
    );
}

#[test]
fn n_aware_matching_is_restricted_to_the_current_node() {
    let test = TestDir::new("hierarchical-n");
    let tree = compile_tree(
        &test,
        "R1_4A4B",
        "Set\tID\tSequence\n\
         A\tA1\tAAAA\n\
         A\tA2\tCCCC\n\
         B\tB1\tGGGG\n\
         B\tB2\tGGGT\n",
        "Sample\tA\tB\n\
         sample1\tA1\tB1\n\
         sample2\tA2\tB2\n",
        1,
    )
    .expect("Close B sequences in independent nodes should compile");

    assert_eq!(
        route(&tree, "AAAAGGGG"),
        RouteResult::Assigned { sample_id: 0 }
    );
    assert_eq!(
        route(&tree, "AAANGGGG"),
        RouteResult::Assigned { sample_id: 0 }
    );
    assert_eq!(route(&tree, "AANNGGGG"), RouteResult::Unmatched);
}

#[test]
fn unused_root_barcode_is_unrouted() {
    let test = TestDir::new("unrouted-root");
    let tree = compile_tree(
        &test,
        "R1_4A4B",
        "Set\tID\tSequence\n\
         A\tA1\tAAAA\n\
         A\tA2\tCCCC\n\
         B\tB1\tGGGG\n",
        "Sample\tA\tB\n\
         sample1\tA1\tB1\n",
        0,
    )
    .expect("Hierarchy should compile");

    assert_eq!(route(&tree, "CCCCGGGG"), RouteResult::Unrouted);
}
