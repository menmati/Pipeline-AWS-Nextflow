#!/usr/bin/env bash
# ------------------------------------------------------------------
# Download each SRA run, subsample 1,000,000 read pairs and compress.
# Raw FASTQ files are deleted before the next run is downloaded, so
# the disk footprint stays below ~1 GB instead of >70 GB.
#
# Usage (from the repository root):
#   bash scripts/download_and_subsample.sh [accession_list] [n_pairs] [out_dir]
# Requires: sra-tools (fasterq-dump), seqtk, gzip
# ------------------------------------------------------------------
set -euo pipefail

ACC_LIST="${1:-data/SRR_Acc_List.txt}"
N_PAIRS="${2:-1000000}"
OUT_DIR="${3:-data/fastq}"
TMP_DIR="data/tmp"
SEED=42
THREADS=4

mkdir -p "$OUT_DIR" "$TMP_DIR"

while read -r SRR; do
    [[ -z "$SRR" ]] && continue
    if [[ -s "${OUT_DIR}/${SRR}_1.fastq.gz" && -s "${OUT_DIR}/${SRR}_2.fastq.gz" ]]; then
        echo "[skip] ${SRR} already subsampled"
        continue
    fi

    echo "[1/3] Downloading ${SRR}..."
    fasterq-dump --split-files --threads "$THREADS" --outdir "$TMP_DIR" "$SRR"

    echo "[2/3] Subsampling ${N_PAIRS} read pairs (seed ${SEED})..."
    # Same seed and same number on both files keeps the read pairs in sync
    seqtk sample -s"$SEED" "${TMP_DIR}/${SRR}_1.fastq" "$N_PAIRS" | gzip > "${OUT_DIR}/${SRR}_1.fastq.gz"
    seqtk sample -s"$SEED" "${TMP_DIR}/${SRR}_2.fastq" "$N_PAIRS" | gzip > "${OUT_DIR}/${SRR}_2.fastq.gz"

    echo "[3/3] Removing raw files for ${SRR}..."
    rm -f "${TMP_DIR}/${SRR}"_*.fastq
done < "$ACC_LIST"

echo "Done. Subsampled files are in ${OUT_DIR}/"
