use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

const READ_LENGTH: usize = 100;

fn main() -> io::Result<()> {
    let mut args = env::args_os().skip(1);
    let reference = args
        .next()
        .expect("usage: generate_fixture <reference.fa> <records> [single|paired]");
    let records: usize = args
        .next()
        .expect("usage: generate_fixture <reference.fa> <records> [single|paired]")
        .to_string_lossy()
        .parse()
        .expect("record count must be an integer");
    let mode = args.next().unwrap_or_else(|| "single".into());
    assert!(args.next().is_none(), "too many arguments");

    let reference_length = match mode.to_str() {
        Some("single") => records.saturating_mul(12).saturating_add(READ_LENGTH),
        Some("paired") => records.saturating_mul(25).saturating_add(READ_LENGTH + 10),
        _ => panic!("mode must be single or paired"),
    }
    .max(100_000);
    write_reference(Path::new(&reference), reference_length)?;
    match mode.to_str() {
        Some("single") => write_single(records, reference_length),
        Some("paired") => write_paired(records, reference_length),
        _ => unreachable!(),
    }
}

fn write_reference(path: &Path, length: usize) -> io::Result<()> {
    let mut output = BufWriter::with_capacity(1024 * 1024, File::create(path)?);
    writeln!(output, ">chrWgbs")?;
    for start in (0..length).step_by(80) {
        let end = length.min(start + 80);
        for position in start..end {
            output.write_all(&[reference_base(position)])?;
        }
        output.write_all(b"\n")?;
    }
    output.flush()
}

fn write_single(records: usize, reference_length: usize) -> io::Result<()> {
    let mut output = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    write_header(&mut output, reference_length)?;
    let quality = "I".repeat(READ_LENGTH);
    for record in 0..records {
        let start = record * 12;
        let top = record % 2 == 0;
        let mut flag = if top { 0 } else { 16 };
        if record % 97 == 0 {
            flag |= 1024;
        }
        let mapq = if record % 997 == 0 { 0 } else { 60 };
        writeln!(
            output,
            "r{record:08}\t{flag}\tchrWgbs\t{}\t{mapq}\t100M\t*\t0\t0\t{}\t{}\tXG:Z:{}",
            start + 1,
            sequence(start, record, top),
            quality,
            if top { "CT" } else { "GA" }
        )?;
    }
    output.flush()
}

fn write_paired(fragments: usize, reference_length: usize) -> io::Result<()> {
    let mut output = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    write_header(&mut output, reference_length)?;
    let quality = "I".repeat(READ_LENGTH);
    for fragment in 0..fragments {
        let first = fragment * 25;
        let second = first + 10;
        let top = fragment % 2 == 0;
        let duplicate = if fragment % 97 == 0 { 1024 } else { 0 };
        let (first_flag, second_flag, tag) = if top {
            (99 | duplicate, 147 | duplicate, "CT")
        } else {
            (83 | duplicate, 163 | duplicate, "GA")
        };
        let mapq = if fragment % 997 == 0 { 0 } else { 60 };
        writeln!(
            output,
            "p{fragment:08}\t{first_flag}\tchrWgbs\t{}\t{mapq}\t100M\t=\t{}\t110\t{}\t{quality}\tXG:Z:{tag}",
            first + 1,
            second + 1,
            sequence(first, fragment, top)
        )?;
        writeln!(
            output,
            "p{fragment:08}\t{second_flag}\tchrWgbs\t{}\t{mapq}\t100M\t=\t{}\t-110\t{}\t{quality}\tXG:Z:{tag}",
            second + 1,
            first + 1,
            sequence(second, fragment, top)
        )?;
    }
    output.flush()
}

fn write_header(output: &mut impl Write, reference_length: usize) -> io::Result<()> {
    writeln!(output, "@HD\tVN:1.6\tSO:coordinate")?;
    writeln!(output, "@SQ\tSN:chrWgbs\tLN:{reference_length}")
}

fn sequence(start: usize, molecule: usize, top: bool) -> String {
    let mut sequence = String::with_capacity(READ_LENGTH);
    for offset in 0..READ_LENGTH {
        let position = start + offset;
        let reference = reference_base(position);
        let methylated = (molecule + position / 20) % 3 != 0;
        let base = match (top, reference, methylated) {
            (true, b'C', false) => b'T',
            (false, b'G', false) => b'A',
            _ => reference,
        };
        sequence.push(char::from(base));
    }
    sequence
}

fn reference_base(position: usize) -> u8 {
    match position % 20 {
        0 => b'C',
        1 => b'G',
        _ if splitmix64(position as u64) & 1 == 0 => b'A',
        _ => b'T',
    }
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}
