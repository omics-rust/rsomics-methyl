use std::process::Command;

fn extract_fixture() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/extract")
}

fn cytosine_fixture() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cytosine-report")
}

#[test]
fn extract_cli_matches_methyldackel_golden() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("result");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let observed = std::fs::read_to_string(directory.path().join("result_CpG.bedGraph"))
        .unwrap()
        .lines()
        .skip(1)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let expected = std::fs::read_to_string(fixture.join("expected.bedGraph"))
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(observed, expected);
}

#[test]
fn extract_region_bounds_the_reported_sites() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("region");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--region",
            "chrSynthetic:5-10",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let lines = std::fs::read_to_string(directory.path().join("region_CpG.bedGraph"))
        .unwrap()
        .lines()
        .skip(1)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let expected = std::fs::read_to_string(fixture.join("expected.bedGraph"))
        .unwrap()
        .lines()
        .filter(|line| {
            let start = line.split('\t').nth(1).unwrap().parse::<u64>().unwrap();
            (4..10).contains(&start)
        })
        .map(str::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(lines, expected);
}

#[test]
fn inclusion_bounds_are_one_based_inclusive_and_accept_zero_sentinels() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("inclusion");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--OT",
            "5,0,0,0",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let observed =
        std::fs::read_to_string(directory.path().join("inclusion_CpG.bedGraph")).unwrap();
    assert!(!observed.contains("chrSynthetic\t1\t2\t"));
    assert!(observed.contains("chrSynthetic\t4\t5\t33\t1\t2\n"));
}

#[test]
fn fixed_end_trimming_matches_live_methyldackel_goldens() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();

    let extract_prefix = directory.path().join("extract");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            extract_prefix.to_str().unwrap(),
            "--nOT",
            "5,1,1,1",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let observed = std::fs::read_to_string(directory.path().join("extract_CpG.bedGraph"))
        .unwrap()
        .lines()
        .skip(1)
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        format!("{observed}\n"),
        std::fs::read_to_string(fixture.join("expected.trim-fixed.bedGraph")).unwrap()
    );

    let mbias_prefix = directory.path().join("mbias");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "mbias",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            mbias_prefix.to_str().unwrap(),
            "--nOT",
            "5,1,1,1",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        std::fs::read(directory.path().join("mbias_mbias.tsv")).unwrap(),
        std::fs::read(fixture.join("expected.trim-fixed.mbias.tsv")).unwrap()
    );
}

#[test]
fn cytosine_reports_match_live_methyldackel_including_zero_coverage() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let reference = fixture.join("synthetic.fa");
    let input = fixture.join("synthetic.bam");

    for (name, extra, expected) in [
        (
            "all",
            vec!["--cytosine_report", "--chg", "--chh"],
            "expected.cytosine-report.tsv",
        ),
        (
            "zero",
            vec!["--format", "cytosine-report", "--chg", "--chh", "-q", "61"],
            "expected.cytosine-report.zero.tsv",
        ),
    ] {
        let prefix = directory.path().join(name);
        let mut arguments = vec![
            "extract",
            reference.to_str().unwrap(),
            input.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ];
        arguments.extend(extra);
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args(arguments)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            std::fs::read(format!("{}.cytosine_report.txt", prefix.display())).unwrap(),
            std::fs::read(fixture.join(expected)).unwrap(),
            "{name}"
        );
    }
}

#[test]
fn cytosine_report_covers_all_contexts_and_reference_boundaries() {
    let fixture = cytosine_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("contexts");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("reference.fa").to_str().unwrap(),
            fixture.join("empty.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--cytosine-report",
            "--chg",
            "--chh",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        std::fs::read(format!("{}.cytosine_report.txt", prefix.display())).unwrap(),
        std::fs::read(fixture.join("expected.tsv")).unwrap()
    );
}

#[test]
fn cytosine_report_region_is_exhaustive_only_inside_the_interval() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("region");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--format",
            "cytosine-report",
            "--chg",
            "--chh",
            "--region",
            "chrSynthetic:5-10",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let expected = std::fs::read_to_string(fixture.join("expected.cytosine-report.tsv"))
        .unwrap()
        .lines()
        .filter(|line| {
            let position = line.split('\t').nth(1).unwrap().parse::<u64>().unwrap();
            (5..=10).contains(&position)
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        std::fs::read_to_string(format!("{}.cytosine_report.txt", prefix.display())).unwrap(),
        format!("{expected}\n")
    );
}

#[test]
fn failed_cytosine_report_preserves_existing_output() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("report");
    let output = directory.path().join("report.cytosine_report.txt");
    std::fs::write(&output, b"keep\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            directory.path().join("missing.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--cytosine-report",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"keep\n");
}

#[test]
fn alternative_extract_formats_match_methyldackel_goldens() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    for (format, suffix, expected_name) in [
        (
            "fraction",
            "_CpG.meth.bedGraph",
            "expected.fraction.bedGraph",
        ),
        ("counts", "_CpG.counts.bedGraph", "expected.counts.bedGraph"),
        ("logit", "_CpG.logit.bedGraph", "expected.logit.bedGraph"),
        ("methyl-kit", "_CpG.methylKit", "expected.methylKit"),
    ] {
        let prefix = directory.path().join(format);
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args([
                "extract",
                fixture.join("synthetic.fa").to_str().unwrap(),
                fixture.join("synthetic.bam").to_str().unwrap(),
                "--output-prefix",
                prefix.to_str().unwrap(),
                "--format",
                format,
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{format}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let observed = std::fs::read_to_string(format!("{}{suffix}", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .map(str::to_owned)
            .collect::<Vec<_>>();
        let expected = std::fs::read_to_string(fixture.join(expected_name))
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        assert_eq!(observed, expected, "{format}");
    }
}

#[test]
fn merged_extract_formats_match_methyldackel_after_combined_depth_filtering() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    for (format, suffix, expected_name) in [
        ("standard", "_CpG.bedGraph", "expected.merged.bedGraph"),
        (
            "fraction",
            "_CpG.meth.bedGraph",
            "expected.merged.fraction.bedGraph",
        ),
        (
            "counts",
            "_CpG.counts.bedGraph",
            "expected.merged.counts.bedGraph",
        ),
        (
            "logit",
            "_CpG.logit.bedGraph",
            "expected.merged.logit.bedGraph",
        ),
    ] {
        let prefix = directory.path().join(format!("merged-{format}"));
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args([
                "extract",
                fixture.join("synthetic.fa").to_str().unwrap(),
                fixture.join("synthetic.bam").to_str().unwrap(),
                "--output-prefix",
                prefix.to_str().unwrap(),
                "--format",
                format,
                "--minimum-depth",
                "4",
                "--merge-context",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{format}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let observed = std::fs::read_to_string(format!("{}{suffix}", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .map(str::to_owned)
            .collect::<Vec<_>>();
        let expected = std::fs::read_to_string(fixture.join(expected_name))
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        assert_eq!(observed, expected, "{format}");
    }
}

#[test]
fn methylkit_rejects_context_merging() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            directory.path().join("result").to_str().unwrap(),
            "--format",
            "methyl-kit",
            "--merge-context",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("cannot merge complementary contexts")
    );
}

#[test]
fn cytosine_report_rejects_context_merging() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            directory.path().join("result").to_str().unwrap(),
            "--format",
            "cytosine-report",
            "--merge-context",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("cannot merge complementary contexts")
    );
}

#[test]
fn per_read_matches_the_documented_methyldackel_contract() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("per-read.tsv");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "per-read",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        std::fs::read(output).unwrap(),
        std::fs::read(fixture.join("expected.per-read.tsv")).unwrap()
    );

    let output = directory.path().join("per-read-all.tsv");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "per-read",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--ignore-nh",
        ])
        .output()
        .unwrap();
    assert!(result.status.success());
    assert!(
        std::fs::read_to_string(output)
            .unwrap()
            .contains("multimapper\tchrSynthetic\t0\t100.000000\t10\n")
    );
}

#[test]
fn mbias_cli_matches_metrics_and_suggestions_from_methyldackel() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("result");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "mbias",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        std::fs::read(directory.path().join("result_mbias.tsv")).unwrap(),
        std::fs::read(fixture.join("expected.mbias.tsv")).unwrap()
    );
    for (strand, suggestion) in [("OT", "--OT 0,0,0,0"), ("OB", "--OB 0,0,0,0")] {
        let svg =
            std::fs::read_to_string(directory.path().join(format!("result_{strand}.svg"))).unwrap();
        assert!(svg.contains(&format!("<title>{strand} M-bias</title>")));
        assert!(svg.contains("id=\"read-1\""));
        assert!(svg.contains(suggestion));
        assert!(svg.ends_with("</svg>\n"));
    }
}

#[test]
fn failed_mbias_preserves_every_existing_output() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("result");
    let table = directory.path().join("result_mbias.tsv");
    let svg = directory.path().join("result_OT.svg");
    std::fs::write(&table, b"keep\n").unwrap();
    std::fs::create_dir(&svg).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "mbias",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(table).unwrap(), b"keep\n");
    assert!(svg.is_dir());
}

#[test]
fn per_read_region_requires_the_alignment_start_inside_the_interval() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let included = directory.path().join("included.tsv");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "per-read",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output",
            included.to_str().unwrap(),
            "--region",
            "chrSynthetic:1-1",
        ])
        .output()
        .unwrap();
    assert!(result.status.success());
    assert_eq!(
        std::fs::read(included).unwrap(),
        std::fs::read(fixture.join("expected.per-read.tsv")).unwrap()
    );

    let excluded = directory.path().join("excluded.tsv");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "per-read",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output",
            excluded.to_str().unwrap(),
            "--region",
            "chrSynthetic:2-10",
        ])
        .output()
        .unwrap();
    assert!(result.status.success());
    assert!(std::fs::read(excluded).unwrap().is_empty());
}

#[test]
fn region_rejects_unknown_and_outside_references() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    for region in ["missing", "chrSynthetic:31"] {
        let output = directory
            .path()
            .join(format!("{}.tsv", region.replace(':', "-")));
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args([
                "per-read",
                fixture.join("synthetic.fa").to_str().unwrap(),
                fixture.join("synthetic.bam").to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--region",
                region,
            ])
            .output()
            .unwrap();
        assert!(!result.status.success(), "{region}");
    }
}

#[test]
fn failed_per_read_preserves_existing_output() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("per-read.tsv");
    std::fs::write(&output, b"keep\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "per-read",
            directory.path().join("missing.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"keep\n");
}

#[test]
fn failed_extract_preserves_every_existing_context_output() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("result");
    let cpg = directory.path().join("result_CpG.bedGraph");
    let chg = directory.path().join("result_CHG.bedGraph");
    std::fs::write(&cpg, b"keep-cpg\n").unwrap();
    std::fs::write(&chg, b"keep-chg\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            directory.path().join("missing.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--chg",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(cpg).unwrap(), b"keep-cpg\n");
    assert_eq!(std::fs::read(chg).unwrap(), b"keep-chg\n");
}

#[test]
fn merge_context_cli_matches_golden() {
    let fixture =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/merge-context");
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("merged.bedGraph");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "merge-context",
            fixture.join("reference.fa").to_str().unwrap(),
            fixture.join("input.bedGraph").to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        std::fs::read(output).unwrap(),
        std::fs::read(fixture.join("expected.bedGraph")).unwrap()
    );
}

#[test]
fn failed_merge_preserves_existing_output() {
    let fixture =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/merge-context");
    let directory = tempfile::tempdir().unwrap();
    let invalid = directory.path().join("invalid.bedGraph");
    let output = directory.path().join("merged.bedGraph");
    std::fs::write(&invalid, b"chr1\t2\t3\t50\t1\t1\nchr1\t1\t2\t50\t1\t1\n").unwrap();
    std::fs::write(&output, b"keep\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "merge-context",
            fixture.join("reference.fa").to_str().unwrap(),
            invalid.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"keep\n");
}

#[test]
fn bed_selection_matches_live_extract_and_mbias_goldens() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    for (label, bed, keep_strand) in [
        ("selection", "selection.bed", None),
        (
            "selection.top",
            "selection.top.bed",
            Some("--keep-bed-strand"),
        ),
        (
            "selection.bottom",
            "selection.bottom.bed",
            Some("--keepStrand"),
        ),
    ] {
        let prefix = directory.path().join(label);
        let mut extract = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"));
        extract.args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--bed",
            fixture.join(bed).to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ]);
        if let Some(argument) = keep_strand {
            extract.arg(argument);
        }
        let result = extract.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let observed = std::fs::read_to_string(format!("{}_CpG.bedGraph", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_eq!(
            observed,
            std::fs::read_to_string(fixture.join(format!("expected.{label}.bedGraph"))).unwrap()
        );

        let mut mbias = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"));
        mbias.args([
            "mbias",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "-l",
            fixture.join(bed).to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ]);
        if let Some(argument) = keep_strand {
            mbias.arg(argument);
        }
        let result = mbias.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            std::fs::read(format!("{}_mbias.tsv", prefix.display())).unwrap(),
            std::fs::read(fixture.join(format!("expected.{label}.mbias.tsv"))).unwrap()
        );
    }
}

#[test]
fn exhaustive_report_obeys_bed_selection() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("selected");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--bed",
            fixture.join("selection.bed").to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--cytosine-report",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        std::fs::read(format!("{}.cytosine_report.txt", prefix.display())).unwrap(),
        std::fs::read(fixture.join("expected.selection.cytosine-report.tsv")).unwrap()
    );
}

#[test]
fn per_read_bed_selection_uses_alignment_spans_and_bisulfite_strands() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    for (label, bed, keep_strand) in [
        ("top", "selection.top.bed", true),
        ("bottom", "selection.bottom.bed", true),
        ("empty", "selection.empty.bed", false),
    ] {
        let output = directory.path().join(format!("{label}.tsv"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"));
        command.args([
            "per-read",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--bed",
            fixture.join(bed).to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ]);
        if keep_strand {
            command.arg("--keep-bed-strand");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        if label == "empty" {
            assert!(std::fs::read(output).unwrap().is_empty());
        } else {
            assert_eq!(
                std::fs::read(output).unwrap(),
                std::fs::read(fixture.join(format!("expected.selection.{label}.per-read.tsv")))
                    .unwrap()
            );
        }
    }
}

#[test]
fn invalid_bed_preserves_per_read_output() {
    let fixture = extract_fixture();
    let directory = tempfile::tempdir().unwrap();
    let bed = directory.path().join("invalid.bed");
    let output = directory.path().join("reads.tsv");
    std::fs::write(&bed, b"chrSynthetic -1 4\n").unwrap();
    std::fs::write(&output, b"keep\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "per-read",
            fixture.join("synthetic.fa").to_str().unwrap(),
            fixture.join("synthetic.bam").to_str().unwrap(),
            "--bed",
            bed.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"keep\n");
}

#[test]
fn conversion_filter_matches_live_extract_and_mbias_goldens() {
    let fixture =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/conversion");
    let directory = tempfile::tempdir().unwrap();
    for threshold in ["0.5", "0.75"] {
        let prefix = directory.path().join(threshold);
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args([
                "extract",
                fixture.join("reference.fa").to_str().unwrap(),
                fixture.join("input.bam").to_str().unwrap(),
                "--minConversionEfficiency",
                threshold,
                "--output-prefix",
                prefix.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let observed = std::fs::read_to_string(format!("{}_CpG.bedGraph", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_eq!(
            observed,
            std::fs::read_to_string(fixture.join(format!("expected.{threshold}.bedGraph")))
                .unwrap()
        );

        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args([
                "mbias",
                fixture.join("reference.fa").to_str().unwrap(),
                fixture.join("input.bam").to_str().unwrap(),
                "--minimum-conversion-efficiency",
                threshold,
                "--output-prefix",
                prefix.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            std::fs::read(format!("{}_mbias.tsv", prefix.display())).unwrap(),
            std::fs::read(fixture.join(format!("expected.{threshold}.mbias.tsv"))).unwrap()
        );
    }
}

#[test]
fn opposite_strand_variant_filter_matches_live_goldens() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/variant");
    let directory = tempfile::tempdir().unwrap();
    for (fraction, expected, alias) in [
        ("0.35", "expected.filtered.bedGraph", true),
        ("0.6", "expected.retained.bedGraph", false),
    ] {
        let prefix = directory.path().join(fraction);
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"));
        command.args([
            "extract",
            fixture.join("reference.fa").to_str().unwrap(),
            fixture.join("input.bam").to_str().unwrap(),
            "--chg",
            "--minimum-opposite-depth",
            "4",
            "--output-prefix",
            prefix.to_str().unwrap(),
        ]);
        if alias {
            command.args(["--maxVariantFrac", fraction]);
        } else {
            command.args(["--maximum-variant-fraction", fraction]);
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let observed = std::fs::read_to_string(format!("{}_CpG.bedGraph", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_eq!(
            observed,
            std::fs::read_to_string(fixture.join(expected)).unwrap()
        );
        let observed = std::fs::read_to_string(format!("{}_CHG.bedGraph", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_eq!(
            observed,
            std::fs::read_to_string(fixture.join(expected.replace(".bedGraph", ".CHG.bedGraph")))
                .unwrap()
        );
    }
}

#[test]
fn variant_filter_matches_all_context_goldens() {
    let fixture =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/conversion");
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("contexts");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("reference.fa").to_str().unwrap(),
            fixture.join("input.bam").to_str().unwrap(),
            "--chg",
            "--chh",
            "--minimum-opposite-depth",
            "4",
            "--maximum-variant-fraction",
            "0.35",
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    for context in ["CpG", "CHG", "CHH"] {
        let observed = std::fs::read_to_string(format!("{}_{context}.bedGraph", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_eq!(
            observed,
            std::fs::read_to_string(fixture.join(format!("expected.variant.{context}.bedGraph")))
                .unwrap()
        );
    }
}

#[test]
fn variant_filter_applies_documented_boundaries_and_merged_exclusion() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/variant");
    let directory = tempfile::tempdir().unwrap();
    for (label, depth, fraction, expected) in [
        ("equal", "4", "0.25", "expected.filtered.bedGraph"),
        ("unknown", "4", "0.45", "expected.filtered.bedGraph"),
        ("depth", "5", "0.35", "expected.retained.bedGraph"),
    ] {
        let prefix = directory.path().join(label);
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
            .args([
                "extract",
                fixture.join("reference.fa").to_str().unwrap(),
                fixture.join("input.bam").to_str().unwrap(),
                "--minOppositeDepth",
                depth,
                "--maximum-variant-fraction",
                fraction,
                "--output-prefix",
                prefix.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(result.status.success());
        let observed = std::fs::read_to_string(format!("{}_CpG.bedGraph", prefix.display()))
            .unwrap()
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_eq!(
            observed,
            std::fs::read_to_string(fixture.join(expected)).unwrap()
        );
    }

    let merged = directory.path().join("merged");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("reference.fa").to_str().unwrap(),
            fixture.join("input.bam").to_str().unwrap(),
            "--chg",
            "--minimum-opposite-depth",
            "4",
            "--maximum-variant-fraction",
            "0.35",
            "--merge-context",
            "--output-prefix",
            merged.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(result.status.success());
    for context in ["CpG", "CHG"] {
        assert_eq!(
            std::fs::read_to_string(format!("{}_{context}.bedGraph", merged.display()))
                .unwrap()
                .lines()
                .count(),
            1
        );
    }

    let report = directory.path().join("report");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-methyl"))
        .args([
            "extract",
            fixture.join("reference.fa").to_str().unwrap(),
            fixture.join("input.bam").to_str().unwrap(),
            "--minimum-opposite-depth",
            "4",
            "--maximum-variant-fraction",
            "0.35",
            "--cytosine-report",
            "--output-prefix",
            report.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(result.status.success());
    assert_eq!(
        std::fs::read(format!("{}.cytosine_report.txt", report.display())).unwrap(),
        std::fs::read(fixture.join("expected.cytosine-report.tsv")).unwrap()
    );
}
