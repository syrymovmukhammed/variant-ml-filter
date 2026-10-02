# variant-ml-filter — ML filter for false-positive variant calls

[![CI](https://github.com/syrymovmukhammed/variant-ml-filter/actions/workflows/main.yml/badge.svg)](https://github.com/syrymovmukhammed/variant-ml-filter/actions/workflows/main.yml)
![python](https://img.shields.io/badge/python-3.12-blue.svg)
![scikit-learn](https://img.shields.io/badge/scikit--learn-RF%20%7C%20LogReg-orange.svg)

A Random Forest and a Logistic Regression learn to tell **real variants from
artefacts** in bacterial WGS calls, using the metrics bcftools writes into the
VCF (QUAL, depth, mapping quality, allele fraction, strand and position bias…).
The models are compared head-to-head with the GATK-style hard filter of
[wgs-variant-pipeline](https://github.com/syrymovmukhammed/wgs-variant-pipeline) on isolates they never saw, and the winner is saved as a reusable
VCF filter.

## Where the labels come from

Machine learning on variants needs labelled calls, and labels need a truth set.
Twelve *E. coli* isolates are simulated (`config/samples.tsv`), each with 150 SNPs +
25 indels planted at known positions — and with the confounders that make
bacterial WGS hard:

| confounder | how it is simulated | the artefact it causes |
|---|---|---|
| low coverage (8–40×), higher error rate | per-isolate `coverage`, `error_rate` | sequencing-error calls |
| contamination / mixed culture | 5–15 % of reads from a second strain | low allele-fraction calls |
| plasmid paralogs | diverged (3 %) copies of chromosomal genes on a contig absent from the reference | collapsed-repeat mismapping — dense clusters of mid-AF calls |

Calls are made **permissively** (`bcftools call -P 0.01`, no filters), then each call is
labelled 1 if it matches the truth set and 0 otherwise.

## Features (`scripts/vcf_features.py`)

`QUAL, DP, AF, MQ, MQ0F, SP, VDB, SGB, RPBZ, MQBZ, BQBZ, MQSBZ, SCBZ, IS_INDEL, INDEL_LEN, IMF`
straight from the VCF, plus three derived ones:

- `DP_RATIO` — depth relative to the isolate's median: collapsed repeats have excess depth;
- `QD` — QUAL / DP (GATK QualByDepth): comparable across coverages;
- `NEIGHBORS` — other calls within ±200 bp: real bacterial SNPs are sparse, paralog artefacts cluster.

## Validation

- **GroupKFold by isolate** on the 10 training isolates — every call of an isolate stays in one fold, so isolate-level properties (coverage, contamination) cannot leak from train to validation.
- **Held-out test isolates** (`T01`, `T02`) — touched once, for the final numbers.
- **Recall is counted against the full truth set**, so a filter is also charged for variants the caller never produced — numbers are directly comparable with the wgs-variant-pipeline benchmark.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 00 | `steps/00_get_reference.sh` | MG1655 reference (NCBI, first 200 kb) |
| 01 | `steps/01_simulate_cohort.sh` | 12 isolates: reads + truth VCF |
| 02 | `steps/02_align_and_call.sh` | bwa mem → markdup → permissive bcftools calls |
| 03 | `steps/03_extract_features.sh` | labelled `features.tsv` |
| 04 | `steps/04_train_models.sh` | CV + test metrics, ROC/PR curves, feature importance, `variant_filter_model.joblib` |
| 05 | `steps/05_apply_filter.sh` | test VCFs annotated with `INFO/ML_PROB`, `FILTER=ML_FAIL`; `results/summary.md` |

## Running

```bash
sudo apt-get install -y bwa samtools bcftools tabix
pip install -r requirements.txt
bash run_all.sh               # ~1 min
bash tests/check_outputs.sh   # CI acceptance test: model must beat the hard filter on held-out isolates
```

## Results (reference local run, offline fallback reference)

Held-out isolates T01 + T02 — 424 raw calls, 350 true variants:

| filter | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|
| none (raw permissive calls) | 342 | 82 | 8 | 0.807 | 0.977 | 0.884 |
| GATK-style hard filter | 307 | 1 | 43 | 0.997 | 0.877 | 0.933 |
| Logistic Regression | 340 | 1 | 10 | 0.997 | 0.971 | **0.984** |
| Random Forest | 339 | 1 | 11 | 0.997 | 0.969 | 0.983 |

Both models remove the same artefacts as the hard filter, but keep ~33 more true variants.
The hard filter's fixed `DP ≥ 10` throws away real variants in low-coverage isolates, while
the models judge depth relative to the isolate (`DP_RATIO`) and combine it with QUAL, AF and
SNP clustering (`NEIGHBORS`). Logistic Regression is selected (best grouped-CV PR-AUC); the
task is close to linearly separable here, so the forest gains nothing — an honest outcome,
not a failure. Real data with messier artefacts is where the forest is expected to pull ahead.

## Using the model on any VCF

```bash
python3 scripts/apply_filter.py \
    --vcf ../wgs-variant-pipeline/results/05_variants/ECOLI_SIM01.raw.vcf.gz \
    --model results/04_models/variant_filter_model.joblib \
  | bcftools view -Oz -o ECOLI_SIM01.ml.vcf.gz
```

(Input must be a bcftools VCF split with `bcftools norm -m -any`. On the raw calls of
[wgs-variant-pipeline](https://github.com/syrymovmukhammed/wgs-variant-pipeline) the model reaches F1 0.988 vs 0.983 for that
pipeline's hard filter.)

## Full-size data

- **Whole chromosome:** `REGION_END=4641652 THREADS=8 bash run_all.sh`
  (reference: [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)).
- **More / different isolates:** edit `config/samples.tsv` (one row per isolate).
- **Real reads with a real truth set:** sequence a strain that also has a finished
  (long-read) assembly, derive its truth VCF by aligning the assembly to the reference
  (the `minimap2 --cs` → `paf_to_vcf.py` route of [popgen-vcf-analysis](https://github.com/syrymovmukhammed/popgen-vcf-analysis)), and label the Illumina calls against it.
  Example pair: ENA run [ERR022075](https://www.ebi.ac.uk/ena/browser/view/ERR022075)
  (Illumina, MG1655) vs the MG1655 reference itself — every PASS call there is a false positive.

## Layout

```
variant-ml-filter/
├── .github/workflows/main.yml   # CI: train + test on every push to main
├── config.sh, lib/common.sh, run_all.sh
├── config/samples.tsv        # cohort design: coverage, error, contamination, paralogs, split
├── steps/00…05_*.sh
├── scripts/
│   ├── wgs_sim.py            # read simulator with confounders (same as in wgs-variant-pipeline)
│   ├── vcf_features.py       # VCF -> feature table (+ labels)
│   ├── train_models.py       # grouped CV, test evaluation, plots, model export
│   ├── apply_filter.py       # model as a VCF filter (ML_PROB / ML_FAIL)
│   └── make_report.py
├── tests/check_outputs.sh
└── requirements.txt / environment.yml
```

## Related repositories

Part of a three-repository set on bacterial WGS:

- [wgs-variant-pipeline](https://github.com/syrymovmukhammed/wgs-variant-pipeline) — reads → BAM → VCF, bcftools + GATK, benchmarked against a truth set
- [variant-ml-filter](https://github.com/syrymovmukhammed/variant-ml-filter) — Random Forest / Logistic Regression filter for false-positive calls
- [popgen-vcf-analysis](https://github.com/syrymovmukhammed/popgen-vcf-analysis) — multi-sample VCF → SNP distances, PCA, clustering, NJ tree
