"""Fast, mostly data-free unit tests for the variant_calling_nf helpers.

These cover the pure logic (variant typing, Ts/Tv, GT-dosage parsing, GTF
interval classification, VCF parsing, the planted-variant self-check, and
demo-data generation) without invoking the heavy aligner / GATK toolchain —
that path is exercised by `run_local.sh --demo`.
Run from the repo root:  pytest variant_calling_nf/tests/
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"


def load(stem: str):
    """Import a bin/ script as a module by file path."""
    path = BIN / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- annotate_variants: typing, Ts/Tv, GT dosage ---------------------------

def test_variant_kind():
    av = load("annotate_variants")
    assert av.variant_kind("A", "G") == "snv"
    assert av.variant_kind("A", "ATG") == "ins"
    assert av.variant_kind("ATG", "A") == "del"
    assert av.variant_kind("AT", "GC") == "mnv"


def test_ts_or_tv():
    av = load("annotate_variants")
    assert av.ts_or_tv("A", "G") == "ts"   # purine <-> purine
    assert av.ts_or_tv("C", "T") == "ts"   # pyrimidine <-> pyrimidine
    assert av.ts_or_tv("A", "C") == "tv"
    assert av.ts_or_tv("G", "T") == "tv"
    assert av.ts_or_tv("A", "ATG") == ""   # not a SNV


def test_gt_dosage():
    av = load("annotate_variants")
    assert av.gt_dosage("0/0") == 0
    assert av.gt_dosage("0/1") == 1
    assert av.gt_dosage("1/1") == 2
    assert av.gt_dosage("1|1") == 2        # phased
    assert av.gt_dosage("0/1:35,12:47") == 1   # extra FORMAT subfields ignored
    assert av.gt_dosage("./.") is None


def test_classify_region(tmp_path):
    av = load("annotate_variants")
    gtf = tmp_path / "g.gtf"
    # one gene with two exons (1-50 and 151-200); 51-150 is intronic
    gtf.write_text(
        'chr1\td\texon\t1\t50\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";\n'
        'chr1\td\texon\t151\t200\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";\n')
    exons, genes = av.read_intervals(gtf)
    assert av.classify_region("chr1", 25, exons, genes) == "exonic"
    assert av.classify_region("chr1", 100, exons, genes) == "intronic"
    assert av.classify_region("chr1", 500, exons, genes) == "intergenic"
    assert av.classify_region("chr2", 25, exons, genes) == "intergenic"


def test_parse_vcf(tmp_path):
    av = load("annotate_variants")
    vcf = tmp_path / "c.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t1/1:28\n")
    recs = list(av.parse_vcf(vcf))
    assert len(recs) == 1
    chrom, pos, ref, alt, filt, samples, dosages = recs[0]
    assert (chrom, pos, ref, alt, filt) == ("chr1", 10, "A", "G", "PASS")
    assert samples == ["S1", "S2"]
    assert dosages == [1, 2]


# ---- check_truth: recall / precision / concordance --------------------------

def test_check_truth_dosage_and_has_alt():
    ct = load("check_truth")
    assert ct.dosage("0/1") == 1
    assert ct.dosage("./.") is None
    assert ct.has_alt({"a": 0, "b": 1}) is True
    assert ct.has_alt({"a": 0, "b": None}) is False


def test_check_truth_perfect_recovery(tmp_path):
    ct = load("check_truth")
    body = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n")
    truth = tmp_path / "t.vcf"; truth.write_text(body)
    called = tmp_path / "c.vcf"; called.write_text(body)
    _, t = ct.load_vcf(truth)
    _, c = ct.load_vcf(called, pass_only=True)
    tsites = {k: v for k, v in t.items() if ct.has_alt(v)}
    csites = {k: v for k, v in c.items() if ct.has_alt(v)}
    assert set(tsites) == set(csites)            # perfect recall + precision


def test_check_truth_respects_pass_only(tmp_path):
    ct = load("check_truth")
    vcf = tmp_path / "c.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tG\t50\tQD2\t.\tGT\t0/1\n")     # filtered out
    _, sites = ct.load_vcf(vcf, pass_only=True)
    assert sites == {}


# ---- make_demo_data: haplotype construction + generation --------------------

def test_build_haplotype_applies_variants():
    dd = load("make_demo_data")
    ref = "AAAACCCCGGGGTTTT"
    snv = [dict(pos0=4, ref="C", alt="T")]
    assert dd.build_haplotype(ref, snv) == "AAAATCCCGGGGTTTT"
    ins = [dict(pos0=4, ref="C", alt="CGG")]
    assert dd.build_haplotype(ref, ins) == "AAAACGGCCCGGGGTTTT"
    # anchor C kept, the next two C's deleted -> one C of the original run remains
    dele = [dict(pos0=4, ref="CCC", alt="C")]
    assert dd.build_haplotype(ref, dele) == "AAAACCGGGGTTTT"


def test_demo_data_generation(tmp_path):
    pytest.importorskip("numpy")
    dd = load("make_demo_data")
    import sys
    argv = ["make_demo_data", "--outdir", str(tmp_path), "--chroms", "1",
            "--chrom-len", "6000", "--reps", "2", "--depth", "10", "--seed", "1"]
    old = sys.argv
    sys.argv = argv
    try:
        dd.main()
    finally:
        sys.argv = old
    assert (tmp_path / "genome.fa").exists()
    assert (tmp_path / "genes.gtf").exists()
    assert (tmp_path / "truth.vcf").exists()
    sheet = (tmp_path / "samplesheet.csv").read_text().splitlines()
    assert sheet[0] == "sample_id,group,fastq_1,fastq_2"
    assert len(sheet) == 1 + 4               # 2 groupA + 2 groupB
    assert (tmp_path / "reads" / "A_s1_R1.fastq.gz").exists()
    # truth VCF carries the sample genotype columns
    truth = [l for l in (tmp_path / "truth.vcf").read_text().splitlines()
             if l.startswith("#CHROM")][0].split("\t")
    assert truth[9:] == ["A_s1", "A_s2", "B_s1", "B_s2"]


# ---- make_overview / make_geno_heatmap: pure helpers ------------------------

def test_tstv_ratio():
    pytest.importorskip("pandas")
    import pandas as pd
    mo = load("make_overview")
    df = pd.DataFrame({"kind": ["snv", "snv", "snv", "ins"],
                       "tstv": ["ts", "ts", "tv", ""]})
    assert mo.tstv_ratio(df) == pytest.approx(2.0)   # 2 ts / 1 tv
