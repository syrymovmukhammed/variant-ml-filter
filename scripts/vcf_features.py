#!/usr/bin/env python3
"""
vcf_features.py — turn a bcftools VCF into a feature table for ML.

One row per ALT allele (input must be split with `bcftools norm -m -any`).
Features are the same metrics a human analyst looks at when hard-filtering
(QUAL, depth, mapping quality, allele fraction, strand bias ...) plus three
derived ones that hard filters usually ignore:

  DP_RATIO   site depth / median depth of the sample's calls. Collapsed repeats
             (paralog reads piling onto one locus) show up as excess depth.
  QD         QUAL / DP — GATK's QualByDepth, robust to coverage differences.
  NEIGHBORS  number of other calls within +-200 bp. True bacterial SNPs are
             sparse; artefacts from mismapped paralogs come in dense clusters.

With --truth, a `label` column is added: 1 = call matches the truth set
(true variant), 0 = artefact (false positive).
"""
from __future__ import annotations

import argparse
import gzip
import sys
from bisect import bisect_left, bisect_right

import numpy as np
import pandas as pd

FEATURES = [
    "QUAL", "DP", "AF", "MQ", "MQ0F", "SP", "VDB", "SGB", "RPBZ", "MQBZ",
    "BQBZ", "MQSBZ", "SCBZ", "IS_INDEL", "INDEL_LEN", "IMF",
    "DP_RATIO", "QD", "NEIGHBORS",
]
NEIGHBOR_WINDOW = 200


def _info_dict(info: str) -> dict[str, str]:
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = "1"  # flags such as INDEL
    return d


def _num(d: dict, key: str, default: float = 0.0) -> float:
    v = d.get(key)
    if v is None or v == ".":
        return default
    try:
        return float(v.split(",")[0])
    except ValueError:
        return default


def parse_vcf(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = f[0], int(f[1]), f[3].upper(), f[4].upper()
            if "," in alt:
                sys.exit(f"{path}: multi-allelic record at {chrom}:{pos}; run bcftools norm -m -any first")
            qual = float(f[5]) if f[5] != "." else 0.0
            info = _info_dict(f[7])
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            ad = [int(x) for x in fmt.get("AD", "0,0").split(",") if x not in (".", "")] or [0, 0]
            ad_ref, ad_alt = (ad + [0, 0])[:2]
            is_indel = int(len(ref) != len(alt))
            rows.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "QUAL": qual,
                "DP": _num(info, "DP"),
                "AF": ad_alt / (ad_ref + ad_alt) if ad_ref + ad_alt else 0.0,
                "MQ": _num(info, "MQ"),
                "MQ0F": _num(info, "MQ0F"),
                "SP": _num(fmt, "SP"),
                "VDB": _num(info, "VDB", 0.5),
                "SGB": _num(info, "SGB"),
                "RPBZ": _num(info, "RPBZ"),
                "MQBZ": _num(info, "MQBZ"),
                "BQBZ": _num(info, "BQBZ"),
                "MQSBZ": _num(info, "MQSBZ"),
                "SCBZ": _num(info, "SCBZ"),
                "IS_INDEL": is_indel,
                "INDEL_LEN": abs(len(ref) - len(alt)),
                "IMF": _num(info, "IMF", 1.0 if not is_indel else 0.0),
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt"] + FEATURES)
    med = float(np.median(df["DP"])) or 1.0
    df["DP_RATIO"] = df["DP"] / med
    df["QD"] = df["QUAL"] / df["DP"].clip(lower=1)
    # calls within +-NEIGHBOR_WINDOW bp on the same contig (excluding itself)
    df["NEIGHBORS"] = 0
    for chrom, idx in df.groupby("chrom").groups.items():
        p = df.loc[idx, "pos"].to_numpy()
        order = np.sort(p)
        df.loc[idx, "NEIGHBORS"] = [
            bisect_right(order, x + NEIGHBOR_WINDOW) - bisect_left(order, x - NEIGHBOR_WINDOW) - 1
            for x in p
        ]
    return df


def truth_keys(path: str) -> set[tuple]:
    opener = gzip.open if path.endswith(".gz") else open
    keys = set()
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                f = line.split("\t")
                keys.add((f[0], int(f[1]), f[3].upper(), f[4].upper()))
    return keys


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--split", default="unknown")
    ap.add_argument("--truth", help="normalised truth VCF -> adds `label`")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = parse_vcf(a.vcf)
    df.insert(0, "split", a.split)
    df.insert(0, "sample", a.sample)
    if a.truth:
        t = truth_keys(a.truth)
        df["label"] = [int((c, p, r, x) in t) for c, p, r, x in
                       zip(df["chrom"], df["pos"], df["ref"], df["alt"])]
        df.attrs["n_truth"] = len(t)
    df.to_csv(a.out, sep="\t", index=False)
    lab = f", {int(df['label'].sum())} true / {int((1 - df['label']).sum())} artefacts" if a.truth else ""
    print(f"{a.sample}: {len(df)} calls{lab}", file=sys.stderr)


if __name__ == "__main__":
    main()
