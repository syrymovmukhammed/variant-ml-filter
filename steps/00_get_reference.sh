#!/usr/bin/env bash
# Step 00 — reference genome: E. coli K-12 MG1655 (RefSeq NC_000913.3).
# Test mode fetches only the first ${REGION_END} bp through NCBI E-utilities
# (a ~200 kb download). If NCBI is unreachable, a synthetic genome with the
# same GC content is generated instead and the run is clearly marked as such.
source "$(dirname "$0")/../lib/common.sh"
need_cmd samtools bwa curl python3
mkdir -p "${DATA_DIR}"

if [[ -n "${REF_FASTA}" ]]; then
  log "using local reference: ${REF_FASTA}"
  cp "${REF_FASTA}" "${REF}"
  echo "LOCAL:${REF_FASTA}" > "${DATA_DIR}/reference.source"
else
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}&rettype=fasta&retmode=text&seq_start=1&seq_stop=${REGION_END}"
  tmp="${DATA_DIR:?}/.ref.download.fa"
  log "downloading ${REF_ACC}:1-${REGION_END} from NCBI"
  if curl -fsSL --retry 4 --retry-delay 5 --max-time 300 "${url}" -o "${tmp}" \
     && head -c1 "${tmp}" | grep -q '>' \
     && [[ "$(grep -v '>' "${tmp}" | tr -d '\nACGTNacgtn' | wc -c)" -eq 0 ]]; then
    # rename the record to the bare accession so VCF contig names are clean
    { echo ">${REF_ACC}"; grep -v '>' "${tmp}"; } > "${REF}"
    echo "NCBI:${REF_ACC}:1-${REGION_END}" > "${DATA_DIR}/reference.source"
  else
    warn "NCBI download failed — generating a SYNTHETIC reference (offline fallback)"
    python3 "${PROJECT_DIR}/scripts/wgs_sim.py" synthetic-ref \
      --out "${REF}" --name synthetic_chr --length "${REGION_END}" --seed "${SEED}"
    echo "SYNTHETIC:offline-fallback" > "${DATA_DIR}/reference.source"
  fi
  rm -f "${tmp:?}"
fi

log "indexing reference"
samtools faidx "${REF}"
bwa index "${REF}" 2> "${DATA_DIR}/bwa_index.log"
log "reference ready: $(cut -f1,2 "${REF}.fai" | tr '\t' ' ') bp  [source: $(cat "${DATA_DIR}/reference.source")]"
