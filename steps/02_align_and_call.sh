#!/usr/bin/env bash
# Step 02 — per isolate: bwa mem -> fixmate/sort/markdup -> bcftools call.
# Calling is deliberately PERMISSIVE (raised prior -P, no filters) so that the
# raw VCF contains both true variants and artefacts; deciding between them is
# the job of the classifier. All annotations a model can learn from are kept.
source "$(dirname "$0")/../lib/common.sh"
need_cmd bwa samtools bcftools
need_file "${REF}.bwt"
mkdir -p "${CALL_DIR}"

while read -r s _; do
  bam="${CALL_DIR}/${s}.bam"
  rg="@RG\tID:${s}\tSM:${s}\tPL:ILLUMINA\tLB:${s}_lib1"
  bwa mem -t "${THREADS}" -R "${rg}" "${REF}" \
      "${SIM_DIR}/${s}_R1.fastq.gz" "${SIM_DIR}/${s}_R2.fastq.gz" 2> /dev/null \
    | samtools fixmate -m - - \
    | samtools sort -@ "${THREADS}" -o - - 2> /dev/null \
    | samtools markdup - "${bam}" 2> /dev/null
  samtools index "${bam}"

  bcftools mpileup -f "${REF}" -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
      --min-MQ 0 --min-BQ 0 --max-depth 2000 -Ou "${bam}" 2> /dev/null \
    | bcftools call -m -v --ploidy 1 -P "${CALL_PRIOR}" -Ou \
    | bcftools norm -f "${REF}" -m -any -Oz -o "${CALL_DIR}/${s}.raw.vcf.gz" 2> /dev/null
  bcftools index -t "${CALL_DIR}/${s}.raw.vcf.gz"
  bcftools norm -f "${REF}" -m -any "${SIM_DIR}/${s}.truth.vcf" -Oz \
      -o "${CALL_DIR}/${s}.truth.norm.vcf.gz" 2> /dev/null
  log "${s}: $(bcftools view -H "${CALL_DIR}/${s}.raw.vcf.gz" | wc -l) raw calls, $(grep -vc '^#' "${SIM_DIR}/${s}.truth.vcf") true variants"
done < <(cohort)
