use std::process::Command;

fn extract_fixture() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/extract")
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
