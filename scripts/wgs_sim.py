#!/usr/bin/env python3
"""
wgs_sim.py — a small, dependency-light simulator for haploid bacterial WGS.

Sub-commands
------------
synthetic-ref   Write a random reference genome (offline fallback only).
simulate        Plant known SNPs/indels into a reference, write the truth VCF,
                the mutated genome, and paired-end Illumina-like FASTQ files.

Why simulate?  Real reads have no answer key.  Planting mutations at known
positions gives an exact truth set, so precision / recall of the variant
caller can be *measured* rather than assumed (the same idea as GIAB
benchmarking, scaled down to a bacterial genome).

Read model
----------
* fragment length ~ Normal(insert_mean, insert_sd), both mates read inwards
* per-base substitution error grows linearly along the read (Illumina-like),
  scaled by a per-read log-normal factor (some reads are worse than others)
* base qualities are the Phred value of the error probability actually used,
  so qualities are calibrated — the caller can trust them, as with real BQSR
Optional confounders (used by the variant-ml-filter repo to create realistic false positives):
* --plasmid-paralogs  diverged copies of chromosomal segments on an extra
                      contig that is NOT in the reference -> collapsed-repeat
                      mismapping, the classic bacterial WGS false positive
* --contam-frac       a fraction of reads from a second, unrelated strain
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i
_COMP = np.frombuffer(b"TGCAN", dtype=np.uint8)  # complement by index A,C,G,T,N
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


# ----------------------------------------------------------------- FASTA I/O
def read_fasta(path: str | Path) -> tuple[str, str]:
    """Return (name, sequence) of the first record; sequence upper-cased."""
    name, chunks = None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                chunks.append(line)
    if name is None:
        sys.exit(f"no FASTA record in {path}")
    return name, "".join(chunks).upper()


def write_fasta(path: str | Path, records: list[tuple[str, str]], width: int = 80) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


# ------------------------------------------------------- synthetic reference
def synthetic_reference(length: int, gc: float, rng: np.random.Generator) -> str:
    """Random genome with E. coli-like GC and three copies of a 1.3 kb
    IS-like element, so that repeats (and their mapping problems) exist."""
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    seq = BASES[rng.choice(4, size=length, p=p)].tobytes().decode()
    is_elem = BASES[rng.choice(4, size=1300, p=p)].tobytes().decode()
    for frac in (0.23, 0.51, 0.78):
        pos = int(length * frac)
        seq = seq[:pos] + is_elem + seq[pos + len(is_elem):]
    return seq[:length]


# ------------------------------------------------------------ mutation model
def plant_variants(seq: str, n_snps: int, n_indels: int, rng: np.random.Generator,
                   min_gap: int = 25, edge: int = 600, ts_frac: float = 0.67):
    """Return (mutant_seq, events) where events = [(pos1, ref, alt, type)].
    Positions are 1-based, VCF-style (indels anchored on the preceding base)."""
    L = len(seq)
    n_total = n_snps + n_indels
    candidates = np.arange(edge, L - edge)
    chosen: list[int] = []
    # rejection sampling keeps variants >= min_gap apart (no overlapping events)
    for pos in rng.permutation(candidates):
        if len(chosen) == n_total:
            break
        if seq[pos] == "N":
            continue
        if all(abs(pos - c) >= min_gap for c in chosen):
            chosen.append(int(pos))
    chosen.sort()
    kinds = np.array(["SNP"] * n_snps + ["INDEL"] * n_indels)
    rng.shuffle(kinds)

    events, out, last = [], [], 0
    for pos0, kind in zip(chosen, kinds):
        ref_base = seq[pos0]
        if kind == "SNP":
            if rng.random() < ts_frac:
                alt = TRANSITION[ref_base]
            else:
                alt = rng.choice([b for b in "ACGT" if b not in (ref_base, TRANSITION[ref_base])])
            out.append(seq[last:pos0]); out.append(alt); last = pos0 + 1
            events.append((pos0 + 1, ref_base, alt, "SNP"))
        else:
            k = int(rng.integers(1, 7))
            if rng.random() < 0.5:  # insertion after pos0
                ins = BASES[rng.integers(0, 4, size=k)].tobytes().decode()
                out.append(seq[last:pos0 + 1]); out.append(ins); last = pos0 + 1
                events.append((pos0 + 1, ref_base, ref_base + ins, "INS"))
            else:                    # deletion of k bases after pos0
                ref = seq[pos0:pos0 + 1 + k]
                out.append(seq[last:pos0 + 1]); last = pos0 + 1 + k
                events.append((pos0 + 1, ref, ref_base, "DEL"))
    out.append(seq[last:])
    return "".join(out), events


def diverge(seq: str, rate: float, rng: np.random.Generator) -> str:
    """Substitute a fraction `rate` of positions (used for paralogs/contaminant)."""
    arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    idx = np.flatnonzero(rng.random(len(arr)) < rate)
    cur = _LUT[arr[idx]]
    ok = cur < 4
    new = (cur[ok] + rng.integers(1, 4, size=ok.sum())) % 4
    arr[idx[ok]] = BASES[new]
    return arr.tobytes().decode()


def write_truth_vcf(path: Path, chrom: str, length: int, sample: str, events) -> None:
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=wgs_sim.py (planted truth set)\n")
        fh.write(f"##contig=<ID={chrom},length={length}>\n")
        fh.write('##INFO=<ID=VT,Number=1,Type=String,Description="Variant type: SNP, INS or DEL">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (haploid)">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for pos, ref, alt, vt in events:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tVT={vt}\tGT\t1\n")


# --------------------------------------------------------------- read model
def _encode(seq: str) -> np.ndarray:
    return _LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]


def simulate_pairs(contigs: list[tuple[str, str, float]], n_pairs: int, read_len: int,
                   ins_mean: float, ins_sd: float, err: float, rng: np.random.Generator):
    """Yield (r1_codes, r1_qual, r2_codes, r2_qual, origin) batches.
    contigs = [(name, seq, relative_weight)]; weight multiplies length."""
    weights = np.array([len(s) * w for _, s, w in contigs], dtype=float)
    weights /= weights.sum()
    per_contig = rng.multinomial(n_pairs, weights)
    # error rises towards the 3' end, as on Illumina instruments (mean ~= err)
    pos_err = err * (0.3 + 2.1 * (np.arange(read_len) / read_len) ** 2)
    for (name, seq, _), n in zip(contigs, per_contig):
        if n == 0:
            continue
        g = _encode(seq)
        L = len(g)
        flen = np.clip(rng.normal(ins_mean, ins_sd, n).astype(int), read_len, L - 1)
        start = (rng.random(n) * (L - flen)).astype(int)
        cols = np.arange(read_len)
        r1 = g[start[:, None] + cols]
        r2 = _COMP[g[(start + flen - read_len)[:, None] + cols]][:, ::-1]
        r2 = _LUT[r2]  # back to codes (complement table returned ASCII)
        flip = rng.random(n) < 0.5          # fragment from the reverse strand
        r1[flip], r2[flip] = r2[flip].copy(), r1[flip].copy()
        reads = []
        for r in (r1, r2):
            read_factor = rng.lognormal(0.0, 0.5, size=(n, 1))
            p = np.clip(pos_err[None, :] * read_factor, 1e-4, 0.3)
            hit = (rng.random(r.shape) < p) & (r < 4)
            r = r.copy()
            r[hit] = (r[hit] + rng.integers(1, 4, size=hit.sum())) % 4
            q = np.clip(np.round(-10 * np.log10(p)), 2, 41).astype(np.uint8)
            reads.append((r, q))
        yield reads[0][0], reads[0][1], reads[1][0], reads[1][1], name


def write_fastq_pair(path1: Path, path2: Path, batches, prefix: str) -> int:
    alphabet = np.frombuffer(b"ACGTN", dtype=np.uint8)
    n = 0
    with gzip.open(path1, "wb", compresslevel=3) as f1, gzip.open(path2, "wb", compresslevel=3) as f2:
        for r1, q1, r2, q2, origin in batches:
            s1, s2 = alphabet[r1], alphabet[r2]
            qq1, qq2 = (q1 + 33), (q2 + 33)
            buf1, buf2 = [], []
            for i in range(r1.shape[0]):
                rid = f"@{prefix}:{origin}:{n + i}"
                buf1.append(b"%s/1\n%s\n+\n%s\n" % (rid.encode(), s1[i].tobytes(), qq1[i].tobytes()))
                buf2.append(b"%s/2\n%s\n+\n%s\n" % (rid.encode(), s2[i].tobytes(), qq2[i].tobytes()))
            f1.write(b"".join(buf1)); f2.write(b"".join(buf2))
            n += r1.shape[0]
    return n


# --------------------------------------------------------------------- CLI
def cmd_synthetic_ref(a) -> None:
    rng = np.random.default_rng(a.seed)
    write_fasta(a.out, [(a.name, synthetic_reference(a.length, a.gc, rng))])
    print(f"synthetic reference: {a.length:,} bp -> {a.out}", file=sys.stderr)


def cmd_simulate(a) -> None:
    rng = np.random.default_rng(a.seed)
    chrom, ref = read_fasta(a.reference)
    mutant, events = plant_variants(ref, a.n_snps, a.n_indels, rng)
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    write_truth_vcf(outdir / f"{a.prefix}.truth.vcf", chrom, len(ref), a.sample, events)
    write_fasta(outdir / f"{a.prefix}.genome.fa", [(f"{a.sample}_chr", mutant)])

    contigs: list[tuple[str, str, float]] = [("chr", mutant, 1.0)]
    # extra plasmid made of diverged copies of chromosomal segments
    if a.plasmid_paralogs > 0:
        parts = []
        for _ in range(a.plasmid_paralogs):
            s = int(rng.integers(1000, len(ref) - a.paralog_len - 1000))
            parts.append(diverge(ref[s:s + a.paralog_len], a.paralog_div, rng))
        contigs.append(("plasmid", "".join(parts), a.plasmid_copies))
    # reads from a second, unrelated strain (lab contamination / mixed culture)
    if a.contam_frac > 0:
        contam = diverge(ref, a.contam_div, rng)
        w = a.contam_frac / (1 - a.contam_frac)
        contigs.append(("contam", contam, w))

    total_len = sum(len(s) * w for _, s, w in contigs)
    n_pairs = int(round(a.coverage * total_len / (2 * a.read_len)))
    batches = simulate_pairs(contigs, n_pairs, a.read_len, a.insert_mean, a.insert_sd, a.error_rate, rng)
    n = write_fastq_pair(outdir / f"{a.prefix}_R1.fastq.gz", outdir / f"{a.prefix}_R2.fastq.gz",
                         batches, a.sample)
    n_snp = sum(e[3] == "SNP" for e in events)
    print(f"{a.sample}: {n_snp} SNPs + {len(events) - n_snp} indels planted; "
          f"{n:,} read pairs ({a.coverage}x, err={a.error_rate}) -> {outdir}", file=sys.stderr)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("synthetic-ref", help="random reference (offline fallback)")
    s.add_argument("--out", required=True)
    s.add_argument("--name", default="synthetic_chr")
    s.add_argument("--length", type=int, default=200_000)
    s.add_argument("--gc", type=float, default=0.508)  # E. coli K-12 GC content
    s.add_argument("--seed", type=int, default=1)
    s.set_defaults(func=cmd_synthetic_ref)

    m = sub.add_parser("simulate", help="plant variants + simulate PE reads")
    m.add_argument("--reference", required=True)
    m.add_argument("--outdir", required=True)
    m.add_argument("--prefix", required=True, help="output file prefix")
    m.add_argument("--sample", required=True)
    m.add_argument("--n-snps", type=int, default=150)
    m.add_argument("--n-indels", type=int, default=25)
    m.add_argument("--coverage", type=float, default=40)
    m.add_argument("--read-len", type=int, default=150)
    m.add_argument("--insert-mean", type=float, default=350)
    m.add_argument("--insert-sd", type=float, default=50)
    m.add_argument("--error-rate", type=float, default=0.004)
    m.add_argument("--plasmid-paralogs", type=int, default=0)
    m.add_argument("--paralog-len", type=int, default=1500)
    m.add_argument("--paralog-div", type=float, default=0.03)
    m.add_argument("--plasmid-copies", type=float, default=1.0)
    m.add_argument("--contam-frac", type=float, default=0.0)
    m.add_argument("--contam-div", type=float, default=0.01)
    m.add_argument("--seed", type=int, default=42)
    m.set_defaults(func=cmd_simulate)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
