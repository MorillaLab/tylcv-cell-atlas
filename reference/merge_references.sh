#!/usr/bin/env bash
# merge_references.sh - merge the tomato and TYLCV references (FASTA + GTF) for `spaceranger mkref`.
#
# Usage:
#   bash merge_references.sh \
#     --tomato-fasta /path/to/ref/fasta/genome.fa \
#     --tomato-gtf   /path/to/ref/genes/genes.gtf.gz \
#     --tylcv-fasta  tylcv_ref/TYLCV_IL.fa \
#     --tylcv-gtf    tylcv_ref/tylcv_spaceranger.gtf \
#     --out          merged_ref_inputs \
#     [--mkref SL_TYLCV] [--threads 16] [--memgb 64]
#
# Writes <out>/genome.fa, <out>/genes.gtf and <out>/MANIFEST.txt (inputs, checksums, counts).
# Plain or gzipped inputs both work; the inputs are never modified. With --mkref NAME it also runs
# `spaceranger mkref` and checks that the TYLCV features made it into the new reference.
set -euo pipefail

usage() { sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
die()  { echo "ERROR: $*" >&2; exit 1; }
ok()   { echo "  ok  $*"; }
cat_any() { gzip -dcf "$1"; }                       # prints plain or gzipped files
md5of() { md5sum "$1" 2>/dev/null | cut -d' ' -f1 || md5 -q "$1"; }

TF="" TG="" VF="" VG="" OUT="" MKREF="" THREADS=16 MEMGB=64
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tomato-fasta) TF=$2; shift 2 ;;
    --tomato-gtf)   TG=$2; shift 2 ;;
    --tylcv-fasta)  VF=$2; shift 2 ;;
    --tylcv-gtf)    VG=$2; shift 2 ;;
    --out)          OUT=$2; shift 2 ;;
    --mkref)        MKREF=$2; shift 2 ;;
    --threads)      THREADS=$2; shift 2 ;;
    --memgb)        MEMGB=$2; shift 2 ;;
    -h|--help)      usage ;;
    *) die "unknown option '$1' (see --help)" ;;
  esac
done
for v in TF TG VF VG OUT; do [[ -n "${!v}" ]] || usage; done
for f in "$TF" "$TG" "$VF" "$VG"; do [[ -s "$f" ]] || die "missing or empty input: $f"; done
[[ -e "$OUT/genome.fa" || -e "$OUT/genes.gtf" ]] && die "$OUT already holds merged files; pick a new --out"
mkdir -p "$OUT"

echo "== 1/4 TYLCV inputs"
VSEQ=$(cat_any "$VF" | tr -d '\r' | awk 'NF')
n=$(printf '%s\n' "$VSEQ" | awk '/^>/' | wc -l)
[[ "$n" -eq 1 ]] || die "$VF must hold exactly one FASTA record (found $n)"
VNAME=$(printf '%s\n' "$VSEQ" | awk '/^>/ {sub(/^>/, ""); print $1}')
VBODY=$(printf '%s\n' "$VSEQ" | awk '!/^>/' | tr -d ' \t\n')
VLEN=${#VBODY}
BAD=$(printf '%s' "$VBODY" | tr -d 'ACGTNRYKMSWBDHVacgtnrykmswbdhv' | wc -c)
[[ "$BAD" -eq 0 ]] || die "$VF contains $BAD characters that are not nucleotides"
ok "FASTA record '$VNAME', $VLEN nt"
VGTF=$(cat_any "$VG" | tr -d '\r' | awk 'NF')
NL=$(printf '%s\n' "$VGTF" | awk '!/^#/' | wc -l)
[[ "$NL" -ge 1 ]] || die "$VG has no feature lines"
printf '%s\n' "$VGTF" | awk -F'\t' -v name="$VNAME" -v len="$VLEN" '
  /^#/ {next}
  NF != 9 {printf "    line %d: %d columns (GTF needs 9, tab-separated)\n", NR, NF; bad=1; next}
  $1 != name {printf "    line %d: seqname \"%s\" differs from the FASTA record \"%s\"\n", NR, $1, name; bad=1}
  ($4+0) < 1 || ($5+0) > (len+0) || ($4+0) > ($5+0) {printf "    line %d: %s-%s outside 1-%s\n", NR, $4, $5, len; bad=1}
  $9 !~ /gene_id "/ || $9 !~ /transcript_id "/ {printf "    line %d: needs gene_id and transcript_id\n", NR; bad=1}
  END {exit bad}' || die "TYLCV GTF problems listed above"
VIDS=$(printf '%s\n' "$VGTF" | awk '!/^#/' | grep -o 'gene_id "[^"]*"' | sed 's/gene_id "//; s/"$//' | sort -u | tr '\n' ' ' || true)
ok "GTF: $NL feature lines on '$VNAME'; gene_id: $VIDS"

echo "== 2/4 Tomato inputs (reads the genome once; may take a minute)"
cat_any "$TF" | awk '/^>/ {sub(/^>/, ""); print $1}' > "$OUT/.tomato_contigs"
NC=$(wc -l < "$OUT/.tomato_contigs")
[[ "$NC" -ge 1 ]] || die "no FASTA headers found in $TF"
DUPS=$(sort "$OUT/.tomato_contigs" | uniq -d | head -3 | tr '\n' ' ' || true)
[[ -z "$DUPS" ]] || die "duplicate contig names in the tomato FASTA: $DUPS"
if grep -qx -- "$VNAME" "$OUT/.tomato_contigs"; then
  die "the tomato FASTA already has a contig named '$VNAME' (an already-merged reference?)"
fi
ok "$NC tomato contigs; none named '$VNAME'"
cat_any "$TG" | awk -F'\t' '!/^#/ && NF >= 9 {print $1}' | sort -u > "$OUT/.gtf_seqnames"
MISSING=$(comm -23 "$OUT/.gtf_seqnames" <(sort -u "$OUT/.tomato_contigs") | head -5 | tr '\n' ' ' || true)
[[ -z "$MISSING" ]] || die "the tomato GTF uses contigs absent from the FASTA ($MISSING) - FASTA and GTF from different assemblies?"
for id in $VIDS; do
  c=$(cat_any "$TG" | grep -c "gene_id \"$id\"" || true)
  [[ "${c:-0}" -eq 0 ]] || die "gene_id \"$id\" already exists in the tomato GTF"
done
ok "tomato GTF matches the FASTA ($(wc -l < "$OUT/.gtf_seqnames") contigs used); no gene_id collision"

echo "== 3/4 Writing merged files"
{ cat_any "$TF"; echo; printf '%s\n' "$VSEQ"; } | awk 'NF' > "$OUT/genome.fa"
{ cat_any "$TG"; echo; printf '%s\n' "$VGTF"; } | awk 'NF' > "$OUT/genes.gtf"
MC=$(awk '/^>/' "$OUT/genome.fa" | wc -l)
[[ "$MC" -eq $((NC + 1)) ]] || die "merged FASTA has $MC contigs, expected $((NC + 1))"
GOT=$(awk -v n=">$VNAME" '/^>/ {p = ($1 == n); next} p' "$OUT/genome.fa" | tr -d '\n')
[[ "$GOT" == "$VBODY" ]] || die "the $VNAME sequence in the merged FASTA differs from the input"
VL=$(awk -F'\t' -v n="$VNAME" '$1 == n' "$OUT/genes.gtf" | wc -l)
[[ "$VL" -eq "$NL" ]] || die "merged GTF holds $VL '$VNAME' lines, expected $NL"
ok "genome.fa: $MC contigs ($NC tomato + $VNAME, sequence verified)"
ok "genes.gtf: $(awk '!/^#/' "$OUT/genes.gtf" | wc -l) feature lines ($NL from $VNAME)"
rm -f "$OUT/.tomato_contigs" "$OUT/.gtf_seqnames"
{
  echo "merge_references.sh run $(date -u +%Y-%m-%dT%H:%MZ)"
  echo "tomato FASTA : $TF   md5 $(md5of "$TF")"
  echo "tomato GTF   : $TG   md5 $(md5of "$TG")"
  echo "TYLCV FASTA  : $VF   md5 $(md5of "$VF")   record $VNAME, $VLEN nt"
  echo "TYLCV GTF    : $VG   md5 $(md5of "$VG")   gene_id $VIDS"
  echo "merged genome.fa md5 $(md5of "$OUT/genome.fa")   contigs $MC"
  echo "merged genes.gtf md5 $(md5of "$OUT/genes.gtf")"
} > "$OUT/MANIFEST.txt"
ok "checksums written to $OUT/MANIFEST.txt"

echo "== 4/4 spaceranger mkref"
if [[ -z "$MKREF" ]]; then
  echo "  Not run (no --mkref). To build the reference:"
  echo "  spaceranger mkref --genome=SL_TYLCV --fasta=$OUT/genome.fa --genes=$OUT/genes.gtf --nthreads=$THREADS --memgb=$MEMGB"
  exit 0
fi
command -v spaceranger >/dev/null || die "spaceranger is not on PATH"
[[ -e "$MKREF" ]] && die "$MKREF already exists; remove it or choose another --mkref name"
spaceranger mkref --genome="$MKREF" --fasta="$OUT/genome.fa" --genes="$OUT/genes.gtf" \
                  --nthreads="$THREADS" --memgb="$MEMGB"
G=$(ls "$MKREF"/genes/genes.gtf* 2>/dev/null | head -1 || true)
[[ -n "$G" ]] || die "mkref finished but $MKREF/genes/genes.gtf is missing"
R=$(cat_any "$G" | awk -F'\t' -v n="$VNAME" '$1 == n' | wc -l)
[[ "$R" -eq "$NL" ]] || die "$MKREF holds $R '$VNAME' GTF lines, expected $NL"
echo "$MKREF built with spaceranger mkref ($(spaceranger --version 2>&1 | head -1))" >> "$OUT/MANIFEST.txt"
ok "reference ready: $(cd "$MKREF" && pwd)  (use it as --transcriptome)"
