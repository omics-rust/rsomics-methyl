use rsomics_methyl::extract::{ExtractOptions, SiteMetric, extract};

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
