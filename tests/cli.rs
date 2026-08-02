use std::process::Command;

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
