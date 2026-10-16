use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use flate2::write::GzEncoder;
use flate2::Compression;
use mutcode::{decode, encode, DecodeOptions, EncodeOptions};

static NEXT: AtomicU64 = AtomicU64::new(0);

struct Fixture(PathBuf);

impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "mutcode-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, text).unwrap();
        path
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn vcf(records: &str) -> String {
    format!(
        "##fileformat=VCFv4.2\n##contig=<ID=g,length=12>\n##contig=<ID=unchanged,length=6>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n{records}"
    )
}

fn options(reference: &Path, vcf: &Path, output: &Path) -> EncodeOptions {
    EncodeOptions {
        reference: reference.into(),
        vcfs: vec![vcf.into()],
        snippy_dir: None,
        output: output.into(),
        sample: None,
        haplotype: None,
        threads: 1,
        check_consensus: None,
    }
}

#[test]
fn all_snippy_literal_types_roundtrip_and_independent_consensus() {
    let fixture = Fixture::new();
    let reference = fixture.file("ref.fa", ">g\nAACCGGTTACGT\n>unchanged\natgNRY\n");
    let input = fixture.file(
        "sample.vcf",
        &vcf("g\t1\t.\tA\tT\t60\t.\tTYPE=snp\tGT\t1/1\n\
         g\t3\t.\tCC\tAA\t60\tPASS\tTYPE=mnp\tGT\t1\n\
         g\t5\t.\tG\tGAC\t60\t.\tTYPE=ins\tGT\t1/1\n\
         g\t7\t.\tTT\tT\t60\t.\tTYPE=del\tGT\t1/1\n\
         g\t10\t.\tCGT\tTANN\t60\t.\tTYPE=complex\tGT\t1/1\n"),
    );
    let expected = fixture.file("expected.fa", ">g\nTAAAGACGTATANN\n>unchanged\natgNRY\n");
    let output = fixture.0.join("codes.tsv");
    let mut config = options(&reference, &input, &output);
    config.check_consensus = Some(expected);
    let summary = encode(&config).unwrap();
    assert_eq!(
        (
            summary.samples,
            summary.genes,
            summary.rows,
            summary.applied_records
        ),
        (1, 2, 2, 5)
    );
    let decoded = fixture.0.join("decoded");
    decode(&DecodeOptions {
        reference,
        codes: output,
        output: Some(decoded.clone()),
        sample: None,
        threads: 2,
    })
    .unwrap();
    assert_eq!(
        fs::read_to_string(decoded.join("s000001.fa")).unwrap(),
        ">g\nTAAAGACGTATANN\n>unchanged\natgNRY\n"
    );
    assert_eq!(
        fs::read_to_string(decoded.join("samples.tsv")).unwrap(),
        "sample\tfasta\nsample\ts000001.fa\n"
    );
}

#[test]
fn gzip_discovery_and_deterministic_threads() {
    use std::io::Write;
    let fixture = Fixture::new();
    let reference = fixture.file("ref.fa", ">g\nAACCGGTTACGT\n>unchanged\nATGNRY\n");
    let directory = fixture.0.join("snippy");
    fs::create_dir(&directory).unwrap();
    let final_vcf = vcf("g\t2\t.\tA\tC,G\t60\t.\tTYPE=snp\tDP:GT\t20:2/2\n");
    fs::write(directory.join("sample.vcf"), &final_vcf).unwrap();
    fs::write(directory.join("sample.raw.vcf"), "must not read").unwrap();
    fs::write(directory.join("sample.filt.vcf"), "must not read").unwrap();
    fs::write(directory.join("sample.subs.vcf"), "must not read").unwrap();
    let mut compressed = GzEncoder::new(
        fs::File::create(directory.join("sample.vcf.gz")).unwrap(),
        Compression::fast(),
    );
    compressed.write_all(final_vcf.as_bytes()).unwrap();
    compressed.finish().unwrap();
    let mut config = options(
        &reference,
        &directory.join("sample.vcf"),
        &fixture.0.join("one.tsv"),
    );
    config.vcfs.clear();
    config.snippy_dir = Some(directory.clone());
    encode(&config).unwrap();
    config.output = fixture.0.join("two.tsv");
    config.threads = 2;
    encode(&config).unwrap();
    assert_eq!(
        fs::read(fixture.0.join("one.tsv")).unwrap(),
        fs::read(&config.output).unwrap()
    );
    fs::remove_file(directory.join("sample.vcf")).unwrap();
    config.output = fixture.0.join("compressed.tsv.gz");
    encode(&config).unwrap();
    let summary = decode(&DecodeOptions {
        reference,
        codes: config.output,
        output: None,
        sample: None,
        threads: 2,
    })
    .unwrap();
    assert_eq!(summary.rows, 2);
}

#[test]
fn fail_explicitly_on_nonreconstructible_inputs_and_leave_no_codes() {
    let fixture = Fixture::new();
    let reference = fixture.file("ref.fa", ">g\nAACCGGTTACGT\n>unchanged\nATGNRY\n");
    for (index, records) in [
        "g\t1\t.\tC\tT\t60\t.\t.\tGT\t1/1\n",
        "g\t1\t.\tA\tT\t60\tLowQual\t.\tGT\t1/1\n",
        "g\t1\t.\tA\tT\t60\t.\t.\tGT\t0/1\n",
        "g\t1\t.\tA\tT\t60\t.\t.\tGT\t./.\n",
        "g\t1\t.\tA\t<DEL>\t60\t.\tEND=4\tGT\t1/1\n",
        "g\t1\t.\tA\t*\t60\t.\t.\tGT\t1/1\n",
        "g\t1\t.\tAA\tT\t60\t.\t.\tGT\t1/1\n\
         g\t2\t.\tA\tG\t60\t.\t.\tGT\t1/1\n",
        "g\t1\t.\tA\tT\t60\t.\t.\tGT\t1/1\n\
         g\t1\t.\tA\tT\t60\t.\t.\tGT\t1/1\n",
        "g\t13\t.\tA\tT\t60\t.\t.\tGT\t1/1\n",
        "other\t1\t.\tA\tT\t60\t.\t.\tGT\t1/1\n",
    ]
    .iter()
    .enumerate()
    {
        let input = fixture.file(&format!("{index}.vcf"), &vcf(records));
        let output = fixture.0.join(format!("{index}.tsv"));
        assert!(
            encode(&options(&reference, &input, &output)).is_err(),
            "{records}"
        );
        assert!(!output.exists());
    }
}

#[test]
fn empty_vcf_encodes_every_gene_and_refuses_wrong_reference_or_truncated_codes() {
    let fixture = Fixture::new();
    let reference = fixture.file("ref.fa", ">g\nAACCGGTTACGT\n>unchanged\nATGNRY\n");
    let input = fixture.file("sample.vcf", &vcf(""));
    let codes = fixture.0.join("codes.tsv");
    encode(&options(&reference, &input, &codes)).unwrap();
    assert!(fs::read_to_string(&codes)
        .unwrap()
        .lines()
        .skip(2)
        .all(|row| row.ends_with("\t[]")));
    let other = fixture.file("other.fa", ">g\nAACCGGTTACGA\n>unchanged\nATGNRY\n");
    assert!(decode(&DecodeOptions {
        reference: other,
        codes: codes.clone(),
        output: None,
        sample: None,
        threads: 2,
    })
    .is_err());
    let text = fs::read_to_string(&codes).unwrap();
    let truncated = fixture.file(
        "truncated.tsv",
        &text.lines().take(3).collect::<Vec<_>>().join("\n"),
    );
    assert!(decode(&DecodeOptions {
        reference: reference.clone(),
        codes: truncated,
        output: None,
        sample: None,
        threads: 2,
    })
    .is_err());
    let damaged = fixture.file("damaged.tsv", &text.replace("\t[]", "\t[[0,1,\"T\"]]"));
    assert!(decode(&DecodeOptions {
        reference,
        codes: damaged,
        output: None,
        sample: None,
        threads: 2,
    })
    .is_err());
}

#[test]
fn multisample_phased_selection_and_output_overwrite_protection() {
    let fixture = Fixture::new();
    let reference = fixture.file("ref.fa", ">g\nAACCGGTTACGT\n>unchanged\nATGNRY\n");
    let text = vcf("g\t1\t.\tA\tC,T\t60\t.\t.\tGT\t1|2\n")
        .replace("\tsample\n", "\tsample\tother\n")
        .replace("\t1|2\n", "\t1|2\t0/0\n");
    let input = fixture.file("sample.vcf", &text);
    let codes = fixture.0.join("codes.tsv");
    let mut config = options(&reference, &input, &codes);
    config.haplotype = Some(2);
    assert_eq!(encode(&config).unwrap().samples, 2);
    assert!(encode(&config).is_err());
    let output = fixture.0.join("decoded");
    let config = DecodeOptions {
        reference,
        codes,
        output: Some(output.clone()),
        sample: Some("sample".into()),
        threads: 2,
    };
    decode(&config).unwrap();
    assert!(fs::read_to_string(output.join("s000002.fa"))
        .unwrap()
        .contains("TACCGGTTACGT"));
    assert!(decode(&config).is_err());
}

#[test]
fn parallel_gene_batches_and_samples_preserve_exact_gzip_and_fasta_output() {
    use std::io::Read;
    let fixture = Fixture::new();
    let dna = "AACCGGTTACGT".repeat(80);
    let reference_text: String = (0..1025).map(|gene| format!(">g{gene}\n{dna}\n")).collect();
    let reference = fixture.file("reference.fa", &reference_text);
    let inputs: Vec<_> = (0..2)
        .map(|sample| {
            let alt = if sample == 0 { "T" } else { "AAC" };
            let records: String = (0..1025)
                .map(|gene| format!("g{gene}\t1\t.\tA\t{alt}\t60\t.\t.\tGT\t1/1\n"))
                .collect();
            fixture.file(
                &format!("sample{sample}.vcf"),
                &format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample{sample}\n{records}"),
            )
        })
        .collect();
    let mut encode_config = options(&reference, &inputs[0], &fixture.0.join("one.tsv.gz"));
    encode_config.vcfs = inputs;
    assert_eq!(encode(&encode_config).unwrap().rows, 2050);
    encode_config.output = fixture.0.join("four.tsv.gz");
    encode_config.threads = 4;
    encode(&encode_config).unwrap();
    assert_eq!(
        fs::read(fixture.0.join("one.tsv.gz")).unwrap(),
        fs::read(&encode_config.output).unwrap()
    );

    let mut decode_config = DecodeOptions {
        reference,
        codes: encode_config.output.clone(),
        output: Some(fixture.0.join("serial")),
        sample: None,
        threads: 1,
    };
    decode(&decode_config).unwrap();
    decode_config.output = Some(fixture.0.join("parallel"));
    decode_config.threads = 4;
    assert_eq!(decode(&decode_config).unwrap().rows, 2050);
    for file in ["s000001.fa", "s000002.fa", "samples.tsv"] {
        assert_eq!(
            fs::read(fixture.0.join("serial").join(file)).unwrap(),
            fs::read(fixture.0.join("parallel").join(file)).unwrap()
        );
    }
    for (index, prefix) in ["T", "AAC"].iter().enumerate() {
        let target = format!("{prefix}{}", &dna[1..]);
        let wrapped: String = target
            .as_bytes()
            .chunks(80)
            .map(|chunk| format!("{}\n", std::str::from_utf8(chunk).unwrap()))
            .collect();
        let expected: String = (0..1025)
            .map(|gene| format!(">g{gene}\n{wrapped}"))
            .collect();
        assert!(
            fs::read_to_string(
                fixture
                    .0
                    .join("parallel")
                    .join(format!("s{:06}.fa", index + 1))
            )
            .unwrap()
                == expected,
            "reconstructed FASTA differs from literal expected alleles for sample {index}"
        );
    }
    decode_config.output = None;
    assert_eq!(decode(&decode_config).unwrap().rows, 2050);
    encode_config.output = fixture.0.join("zero.tsv");
    encode_config.threads = 0;
    assert!(encode(&encode_config).is_err());
    decode_config.threads = 0;
    assert!(decode(&decode_config).is_err());

    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(&decode_config.codes).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let mut rows: Vec<_> = text.lines().map(str::to_owned).collect();
    rows[514] = rows[514].replace("\tg512\t", "\tg511\t");
    decode_config.codes = fixture.file("corrupt.tsv", &rows.join("\n"));
    decode_config.output = Some(fixture.0.join("must-not-publish"));
    decode_config.threads = 4;
    assert!(decode(&decode_config).is_err());
    assert!(!decode_config.output.unwrap().exists());
}
