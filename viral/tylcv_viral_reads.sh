#!/usr/bin/env bash
# tylcv_viral_reads.sh - TYLCV UMIs per Visium spot, two routes:
#   A (recommended): merged reference + re-run -> use spaceranger_tylcv.sh (this script only points there).
#                    The matrix then contains two viral features: TYLCV_V and TYLCV_C.
#   B (quick check): take the reads that did not map to tomato from the existing BAMs, align them to TYLCV
#                    and count UMIs per barcode -> one <sample>.csv per section for the notebook (cell 14).
#
# Assumptions to verify before trusting strand-specific numbers:
#   * TYLCV_FA is the infectious-clone sequence in virion-sense orientation (standard GenBank numbering);
#   * Visium v1 3' libraries are sense-stranded: reads on '+' come from virion-sense transcripts (V1/V2),
#     reads on '-' from complementary-sense transcripts (C1-C4). Total viral UMIs do not depend on this.
#
# Usage:  bash tylcv_viral_reads.sh A     or     bash tylcv_viral_reads.sh B
set -euo pipefail

ROUTE="${1:-B}"
TYLCV_FA="TYLCV_IL.fa"            # from make_tylcv_reference.py (or your clone), single record, virion sense
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# ---- route A inputs: the SAME FASTA/GTF used to build your original Space Ranger reference. Inside a Space Ranger
#      reference folder they are <ref>/fasta/genome.fa and <ref>/genes/genes.gtf(.gz); gzipped files are fine.
GENOME_FA="/path/to/your_reference/fasta/genome.fa"
GENES_GTF="/path/to/your_reference/genes/genes.gtf.gz"
REF_NAME="SL_plus_TYLCV"

# ---- route B inputs: tab-separated "<section_name> <path/to/possorted_genome_bam.bam>", one line per section.
#      Section names must match the notebook (e.g. infected_7dpi_1). Pair each BAM with the section whose loom
#      (made from that BAM) matched it in the notebook's 02_loom_mapping.csv - the A1/A2 7-dpi files are crossed.
#      Cell 14 double-checks this: it warns if a CSV's viral UMIs sit on another section's tissue.
BAM_LIST="bams.tsv"
OUTDIR="viral_counts"

route_a() {
  echo "Route A (merged tomato + TYLCV reference) now lives in spaceranger_tylcv.sh:"
  echo "  bash spaceranger_tylcv.sh check | reference | commands"
}

find_tool() {   # if a tool is installed but not on PATH (common with Homebrew on Apple Silicon), add its folder
  command -v "$1" >/dev/null 2>&1 && return 0
  local d
  for d in /opt/homebrew/bin /usr/local/bin "$HOME/homebrew/bin" "${CONDA_PREFIX:-/nonexistent}/bin" \
           "$HOME/miniforge3/bin" "$HOME/mambaforge/bin" "$HOME/miniconda3/bin" "$HOME/anaconda3/bin"; do
    if [[ -x "$d/$1" ]]; then
      export PATH="$d:$PATH"
      echo "   note: using $d/$1 (that folder is not on your PATH; see the fix printed by: bash $0 --path-help)"
      return 0
    fi
  done
  return 1
}

path_help() {
  cat <<'TXT'
samtools/minimap2 installed with Homebrew but "not found"? Homebrew on Apple Silicon lives in /opt/homebrew/bin,
which must be added to your shell once:
  echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
  eval "$(/opt/homebrew/bin/brew shellenv)"
Then check:  samtools --version | head -1 ;  minimap2 --version
(On Intel Macs Homebrew uses /usr/local/bin, which is normally on PATH already.)
Simplest alternative if you have conda: one environment with all three tools (Apple Silicon builds exist):
  conda create -y -n tylcv_tools -c conda-forge -c bioconda samtools minimap2 pysam
  conda activate tylcv_tools
TXT
}

align_unmapped() {   # $1 = index | scan  (uses SAMPLE, BAM, THREADS, OUTDIR, TYLCV_FA)
  local SRC
  if [[ "$1" == index ]]; then SRC=(samtools view -u -@ "$THREADS" "$BAM" '*')
  else SRC=(samtools view -u -f 4 -@ "$THREADS" "$BAM"); fi
  "${SRC[@]}" \
    | samtools fastq -T CB,UB - 2> "$OUTDIR/$SAMPLE.log" \
    | minimap2 -ax sr -y -t "$THREADS" "$TYLCV_FA" - 2>> "$OUTDIR/$SAMPLE.log" \
    | samtools view -b -F 0x904 -q 1 - > "$OUTDIR/$SAMPLE.tylcv.bam"
}

reads_done() { grep -o 'processed [0-9]* reads' "$OUTDIR/$SAMPLE.log" 2>/dev/null | tail -1 | awk '{print $2}' || true; }

route_b() {
  find_tool samtools || { echo "samtools not found."; path_help; exit 1; }
  find_tool minimap2 || { echo "minimap2 not found (brew install minimap2)."; path_help; exit 1; }
  python3 -c "import pysam" 2>/dev/null || {
    echo "python3 ($(command -v python3)) cannot import pysam: activate the environment where you installed it"
    echo "(e.g. source ~/tylcv/bin/activate) or run: python3 -m pip install pysam"; exit 1; }
  echo "tools: $(command -v samtools) | $(command -v minimap2) | $(command -v python3)"
  mkdir -p "$OUTDIR"
  while IFS=$'\t' read -r SAMPLE BAM <&3; do
    [[ -z "${SAMPLE:-}" || "$SAMPLE" == \#* ]] && continue
    [[ -s "${BAM:-}" ]] || { echo "BAM not found for '$SAMPLE': '${BAM:-}' (bams.tsv columns must be TAB-separated)"; exit 1; }
    echo ">> $SAMPLE  ($BAM)"
    # Unmapped reads, streamed. With an index, '*' jumps straight to the unplaced reads at the end of a sorted BAM.
    # The index is trusted only if the number of reads it promises equals the number actually read; otherwise
    # (or without an index) the whole BAM is scanned with -f 4.
    mode=scan
    if [[ -s "$BAM.bai" || -s "${BAM%.bam}.bai" ]]; then
      expected=$(samtools idxstats "$BAM" 2>/dev/null | awk '$1 == "*" {print $4}' || true)
      if align_unmapped index; then
        got=$(reads_done)
        if [[ -n "$expected" && "$got" == "$expected" ]]; then
          mode=index; echo "   index verified: $got unmapped reads read directly"
        else
          echo "   index does not match this BAM (index: ${expected:-?} unmapped reads, read: ${got:-?}) -> full scan"
        fi
      else
        echo "   reading through the index failed -> full scan"
      fi
    fi
    if [[ "$mode" == scan ]]; then
      align_unmapped scan
      echo "   full scan: $(reads_done) unmapped reads"
    fi
    python3 count_viral_umis.py "$OUTDIR/$SAMPLE.tylcv.bam" --out "$OUTDIR/$SAMPLE.csv"
  done 3< "$BAM_LIST"
  echo "Copy $OUTDIR/*.csv to MyDrive/spatial_transcriptomics/viral_counts/ and run the notebook (cell 14)."
  echo "Mock sections should give ~0 viral UMIs: they define the background threshold."
}

[[ "$ROUTE" == "--path-help" ]] && { path_help; exit 0; }
case "$ROUTE" in
  A) route_a ;;
  B) route_b ;;
  *) echo "usage: bash $0 A|B"; exit 1 ;;
esac
