use plexless::encoder::MAX_BARCODE_LENGTH;
use plexless::structure::{Orientation, ReadLayout, ReadMate, ReadStructure};

#[test]
fn repeated_symbols_are_planned_in_r1_then_r2_order() {
    let layout = ReadLayout::paired(Some("R1_2A2B"), Some("R2_2A2B")).expect("Layout should parse");
    let plans = layout
        .compile_extraction_plans()
        .expect("Plans should compile");
    let mut buffer = [0; MAX_BARCODE_LENGTH];

    let a = plans.plan_index(b'A').expect("A plan should exist");
    let a_len = plans
        .assemble(a, b"ACGT", Some(b"TGCA"), &mut buffer)
        .expect("Assembly should succeed")
        .expect("Reads should be long enough");
    assert_eq!(&buffer[..a_len], b"ACTG");

    let b = plans.plan_index(b'B').expect("B plan should exist");
    let b_len = plans
        .assemble(b, b"ACGT", Some(b"TGCA"), &mut buffer)
        .expect("Assembly should succeed")
        .expect("Reads should be long enough");
    assert_eq!(&buffer[..b_len], b"GTCA");
}

#[test]
fn reverse_complement_applies_only_to_the_marked_piece() {
    let layout = ReadLayout::paired(Some("R1_4A"), Some("R2_4A(rc)")).expect("Layout should parse");
    let plans = layout
        .compile_extraction_plans()
        .expect("Plans should compile");
    let mut buffer = [0; MAX_BARCODE_LENGTH];
    let a = plans.plan_index(b'A').expect("A plan should exist");
    let length = plans
        .assemble(a, b"ACGT", Some(b"CCAT"), &mut buffer)
        .expect("Assembly should succeed")
        .expect("Reads should be long enough");

    assert_eq!(&buffer[..length], b"ACGTATGG");
    assert_eq!(
        plans.plan(a).expect("plan").pieces[0].orientation,
        Orientation::Forward
    );
    assert_eq!(
        plans.plan(a).expect("plan").pieces[1].orientation,
        Orientation::ReverseComplement
    );
}

#[test]
fn asymmetric_layout_has_correct_piece_mates_and_trim_prefixes() {
    let layout =
        ReadLayout::paired(Some("R1_10A11B"), Some("R2_10A(rc)")).expect("Layout should parse");
    let plans = layout
        .compile_extraction_plans()
        .expect("Plans should compile");
    let a = plans
        .plan(plans.plan_index(b'A').expect("A plan"))
        .expect("A plan");
    let b = plans
        .plan(plans.plan_index(b'B').expect("B plan"))
        .expect("B plan");

    assert_eq!(a.logical_length, 20);
    assert_eq!(a.pieces.len(), 2);
    assert_eq!(a.pieces[0].mate, ReadMate::R1);
    assert_eq!(a.pieces[1].mate, ReadMate::R2);
    assert_eq!(b.logical_length, 11);
    assert_eq!(b.pieces.len(), 1);
    assert_eq!(plans.r1_prefix_len(), 21);
    assert_eq!(plans.r2_prefix_len(), 10);
}

#[test]
fn malformed_orientation_syntax_is_rejected() {
    for structure in ["R1_4A(rc", "R1_4A(RC)", "R1_4A(foo)", "R1_4A(rc)(rc)"] {
        let error = ReadStructure::parse(structure).expect_err("Malformed modifier must fail");
        assert!(error.contains("orientation"), "unexpected error: {error}");
    }

    let error = ReadStructure::parse("R1_4T(rc)").expect_err("T cannot be oriented");
    assert!(error.contains("only valid for barcode"));
}

#[test]
fn generic_symbols_except_t_are_accepted() {
    let layout = ReadLayout::single(Some("R1_1A1S1T1U1Z")).expect("Generic symbols should parse");
    assert_eq!(layout.barcode_symbols(), vec![b'A', b'S', b'U', b'Z']);
}
