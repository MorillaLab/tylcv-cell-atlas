#!/usr/bin/env bash
# make_bams_tsv.sh - prepare the BAMs for `bash tylcv_viral_reads.sh B`.
#
# Usage: bash make_bams_tsv.sh "/path/to/folder holding the BAMs (searched recursively)" [bams.tsv]
#
# 1. Extracts the .bai index files from any .zip under that folder (into <folder>/_bai_extracted/).
# 2. Finds each section's position-sorted Space Ranger BAM by facility sample code (N0522...), else by name.
#    velocyto's "cellsorted_*.bam" copies are ignored (used only if a section has no original BAM).
# 3. Pairs each BAM with its index: exact file name first, else the same sample code; no index -> full scan.
# 4. Links BAM and index into ./bam_links/<section>.bam(.bai) (nothing is copied or moved) and writes bams.tsv.
#
# Sections follow your AnnData. The 7-dpi infected pair is CROSSED relative to the file names:
# N0523 -> infected_7dpi_1 and N0522 -> infected_7dpi_2 (spot counts 212/290 in your velocity step; confirm in
# 02_loom_mapping.csv). tylcv_viral_reads.sh re-checks every index, and notebook cell 14 every pairing.
set -euo pipefail
ROOT="${1:?usage: bash make_bams_tsv.sh /path/to/folder_with_the_BAMs [bams.tsv]}"
OUT="${2:-bams.tsv}"
[[ -d "$ROOT" ]] || { echo "Not a folder: $ROOT"; exit 1; }
abspath() { (cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")"); }

code_of() {             # facility sample code of each section
  case "$1" in
    control_7dpi_1) echo N0524 ;;   control_7dpi_2) echo N0525 ;;
    control_24dpi_1) echo N0036 ;;  control_24dpi_2) echo N0037 ;;
    infected_7dpi_1) echo N0523 ;;  infected_7dpi_2) echo N0522 ;;    # crossed pair
    infected_24dpi_1) echo N0034 ;; infected_24dpi_2) echo N0035 ;;
  esac
}
fallback_of() {         # "<file-name regex>|<path regex for the time point>" if the code is not in the name
  case "$1" in
    control_7dpi_1) echo 'A3[-_ ]*mock|(^|[^0-9])7[-_ ]?dpi' ;;   control_7dpi_2) echo 'A4[-_ ]*mock|(^|[^0-9])7[-_ ]?dpi' ;;
    control_24dpi_1) echo 'A3[-_ ]*mock|24[-_ ]?dpi' ;;           control_24dpi_2) echo 'A4[-_ ]*mock|24[-_ ]?dpi' ;;
    infected_7dpi_1) echo 'A2[-_ ]*tylcv|(^|[^0-9])7[-_ ]?dpi' ;;  infected_7dpi_2) echo 'A1[-_ ]*tylcv|(^|[^0-9])7[-_ ]?dpi' ;;
    infected_24dpi_1) echo 'A1[-_ ]*tylcv|24[-_ ]?dpi' ;;          infected_24dpi_2) echo 'A2[-_ ]*tylcv|24[-_ ]?dpi' ;;
  esac
}
rel() { printf '%s\n' "${1#"$ROOT"/}"; }      # path below ROOT
by_code() {             # $1 = newline-separated paths, $2 = sample code (searched in folder and file names below ROOT)
  printf '%s\n' "$1" | while IFS= read -r b; do
    [[ -z "$b" ]] && continue
    rel "$b" | grep -i -E "(^|[^0-9])$2([^0-9]|$)" >/dev/null && echo "$b"
  done || true
}
codes_in() { printf '%s\n' "$1" | grep -o -E 'N0[0-9]{3}' | sort -u | tr '\n' ' ' || true; }
by_name() {             # $1 = paths, $2 = "<name regex>|<path regex>"
  printf '%s\n' "$1" | while IFS= read -r b; do
    [[ -z "$b" ]] && continue
    basename "$b" | grep -i -E "${2%%|*}" >/dev/null && printf '%s\n' "$b" | grep -i -E "${2#*|}" >/dev/null && echo "$b"
  done || true
}
count() { printf '%s' "$1" | grep -c . || true; }

echo "== 1. Index files inside zips"
EXTRACT="$ROOT/_bai_extracted"; skipped=0
rm -rf "$EXTRACT"                       # this script's own folder: always rebuilt from the zips
while IFS= read -r z; do
  [[ -z "$z" ]] && continue
  if unzip -l "$z" 2>/dev/null | grep -i '\.bai$' >/dev/null; then
    dest="$EXTRACT/$(basename "$z" .zip)"; mkdir -p "$dest"
    unzip -o -q "$z" '*.bai' -d "$dest"          # keep folders: several indexes can share one file name
    echo "  ok  $(basename "$z"): $(find "$dest" -type f -name '*.bai' | grep -c . || true) index file(s)"
  else
    skipped=$((skipped + 1))
  fi
done < <(find -L "$ROOT" -type f -iname '*.zip' -not -path '*/_bai_extracted/*' 2>/dev/null)
echo "  ($skipped other zip files hold no index files and were left alone)"

ALL=$(find -L "$ROOT" -type f -name '*.bam' 2>/dev/null | grep -v -E '\.tylcv\.bam$' | grep -v '/bam_links/' || true)
CELLSORTED=$(printf '%s\n' "$ALL" | grep -E '/cellsorted_[^/]*$' || true)
BAMS=$(printf '%s\n' "$ALL" | grep -v -E '/cellsorted_[^/]*$' || true)
BAIS=$(find -L "$ROOT" -type f -name '*.bai' 2>/dev/null | grep -v '/bam_links/' || true)
[[ -n "$ALL" ]] || { echo "No .bam files found under $ROOT"; exit 1; }
echo "== 2. Found $(count "$BAMS") original BAMs, $(count "$CELLSORTED") velocyto cellsorted copies (ignored) and $(count "$BAIS") index files"
printf '%s\n' "$BAMS" | while IFS= read -r b; do [[ -n "$b" ]] && echo "      $(basename "$b")"; done

echo "== 3. Pairing (section <- BAM, index)"
mkdir -p bam_links
{ echo "# section<TAB>BAM - written by make_bams_tsv.sh; links in bam_links/ point to your original files"
  echo "# 7-dpi infected pair crossed on purpose: N0523 -> infected_7dpi_1, N0522 -> infected_7dpi_2"; } > "$OUT"
ok=0
for s in control_7dpi_1 control_7dpi_2 control_24dpi_1 control_24dpi_2 \
         infected_7dpi_1 infected_7dpi_2 infected_24dpi_1 infected_24dpi_2; do
  code=$(code_of "$s"); fb=$(fallback_of "$s")
  hits=$(by_code "$BAMS" "$code"); how="sample code $code"
  if [[ -z "$hits" ]]; then hits=$(by_name "$BAMS" "$fb"); how="name pattern"; fi
  if [[ -z "$hits" && -n "$CELLSORTED" ]]; then hits=$(by_code "$CELLSORTED" "$code"); how="cellsorted copy, full scan"; fi
  n=$(count "$hits")
  if [[ "$n" -ne 1 ]]; then
    tag=AMBIGUOUS; [[ "$n" -eq 0 ]] && tag=MISSING
    printf '# %s %s: %s candidate BAMs - add the right line by hand\n' "$tag" "$s" "$n" >> "$OUT"
    printf '  !!  %-17s %s matching BAMs\n' "$s" "$n"
    printf '%s\n' "$hits" | while IFS= read -r h; do [[ -n "$h" ]] && echo "        $h"; done
    continue
  fi
  bam="$hits"; bb=$(basename "$bam"); bai=""; ihow=""
  if [[ "$how" != cellsorted* ]]; then
    bai=$(printf '%s\n' "$BAIS" | while IFS= read -r i; do
            [[ -z "$i" ]] && continue
            ib=$(basename "$i"); c=$(codes_in "$(rel "$i")")
            if [[ "$ib" == "$bb.bai" || "$ib" == "${bb%.bam}.bai" ]] && [[ -z "$c" || "$c" == "$code " ]]; then echo "$i"; fi
          done | head -1 || true)
    ihow="exact name"
    if [[ -z "$bai" ]]; then
      cand=$(by_code "$BAIS" "$code")
      if [[ $(count "$cand") -eq 1 ]]; then bai="$cand"; ihow="sample code in path"; fi
    fi
  fi
  ln -sf "$(abspath "$bam")" "bam_links/$s.bam"
  rm -f "bam_links/$s.bam.bai"
  if [[ -n "$bai" ]]; then ln -sf "$(abspath "$bai")" "bam_links/$s.bam.bai"; idx="index: $(basename "$bai") ($ihow)"
  else idx="no index -> full scan"; fi
  hdr="header not read (samtools not on PATH)"
  if command -v samtools >/dev/null 2>&1; then
    hc=$(codes_in "$(samtools view -H "$bam" 2>/dev/null | grep '^@RG' || true)")   # read groups only: @PG lines echo file names
    if [[ -z "$hc" ]]; then hdr="read groups name no sample code"
    elif [[ "$hc" == "$code " ]]; then hdr="read group confirms $code"
    else hdr="READ GROUP SAYS ${hc% }, file name says $code: check this file"; fi
  fi
  printf '%s\t%s\n' "$s" "$(abspath "bam_links/$s.bam")" >> "$OUT"; ok=$((ok + 1))
  printf '  ok  %-17s <- %s  [%s; %s; %s]\n' "$s" "$(rel "$bam")" "$how" "$idx" "$hdr"
done
echo "Wrote $OUT ($ok of 8 sections). Next: bash tylcv_viral_reads.sh B"
