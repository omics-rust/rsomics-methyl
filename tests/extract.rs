use rsomics_methyl::extract::{ExtractEvent, ExtractOptions, SiteMetric, extract, extract_events};

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/extract")
        .join(name)
}

fn data_line(metric: &SiteMetric) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        metric.chromosome(),
        metric.start(),
        metric.end(),
        metric.percentage(),
        metric.methylated(),
        metric.unmethylated()
    )
}

#[test]
fn default_cpg_calls_match_methyldackel_golden() {
    let mut observed = Vec::new();
    let stats = extract(
        &fixture("synthetic.bam"),
        &fixture("synthetic.fa"),
        ExtractOptions::default(),
        |metric| {
            observed.push(data_line(&metric));
            Ok(())
        },
    )
    .unwrap();
    let expected = std::fs::read_to_string(fixture("expected.bedGraph"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with("track"))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(observed, expected);
    assert_eq!(stats.input_records, 8);
    assert_eq!(stats.filtered_records, 3);
    assert_eq!(stats.emitted_sites, 20);
}

#[test]
fn stricter_mapq_filters_all_records() {
    let mut observed = Vec::new();
    let stats = extract(
        &fixture("synthetic.bam"),
        &fixture("synthetic.fa"),
        ExtractOptions {
            minimum_mapping_quality: 61,
            ..ExtractOptions::default()
        },
        |metric| {
            observed.push(metric);
            Ok(())
        },
    )
    .unwrap();
    assert!(observed.is_empty());
    assert_eq!(stats.input_records, 8);
    assert_eq!(stats.input_records, stats.filtered_records);
}

#[test]
fn variant_events_report_the_usable_opposite_strand_evidence() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/variant");
    let mut events = Vec::new();
    let stats = extract_events(
        &fixture.join("input.bam"),
        &fixture.join("reference.fa"),
        ExtractOptions {
            minimum_opposite_depth: 4,
            maximum_variant_fraction: 0.35,
            ..ExtractOptions::default()
        },
        |event| {
            events.push(event);
            Ok(())
        },
    )
    .unwrap();

    assert_eq!(stats.excluded_variant_sites, 1);
    assert_eq!(stats.emitted_sites, 1);
    let ExtractEvent::ExcludedVariant(excluded) = &events[0] else {
        panic!("first event should exclude the forward cytosine");
    };
    assert_eq!(excluded.start(), 0);
    assert_eq!(excluded.opposite_depth(), 4);
    assert_eq!(excluded.variant_bases(), 2);
    let ExtractEvent::Site(metric) = &events[1] else {
        panic!("second event should retain the reverse cytosine");
    };
    assert_eq!(metric.start(), 1);
    assert_eq!(metric.methylated(), 5);
}
