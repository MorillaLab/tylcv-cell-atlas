#!/usr/bin/env python3
"""Count unique TYLCV UMIs per Visium barcode and strand from a BAM/SAM of reads aligned to the virus.

Reads must carry the CB (corrected barcode) and UB (UMI) tags that Space Ranger writes; route B of
tylcv_viral_reads.sh keeps them through samtools fastq -T CB,UB and minimap2 -y.
Forward alignments -> TYLCV_V (virion-sense transcripts), reverse -> TYLCV_C (complementary-sense),
assuming sense-stranded Visium 3' reads and a virion-sense reference. Output: barcode,TYLCV_V,TYLCV_C,total
"""
import argparse
import collections
import csv

import pysam


def count(bam_path, min_mapq=1):
    umis = collections.defaultdict(set)
    used = skipped = 0
    with pysam.AlignmentFile(bam_path, "r", check_sq=False) as fh:
        for read in fh:
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
                skipped += 1
                continue
            if not (read.has_tag("CB") and read.has_tag("UB")):
                skipped += 1
                continue
            feature = "TYLCV_C" if read.is_reverse else "TYLCV_V"
            umis[(read.get_tag("CB"), feature)].add(read.get_tag("UB"))
            used += 1
    table = collections.defaultdict(lambda: {"TYLCV_V": 0, "TYLCV_C": 0})
    for (barcode, feature), s in umis.items():
        table[barcode][feature] = len(s)
    return table, used, skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=1)
    args = ap.parse_args()
    table, used, skipped = count(args.bam, args.min_mapq)
    with open(args.out, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["barcode", "TYLCV_V", "TYLCV_C", "total"])
        for bc in sorted(table):
            c = table[bc]
            w.writerow([bc, c["TYLCV_V"], c["TYLCV_C"], c["TYLCV_V"] + c["TYLCV_C"]])
    print(f"{args.bam}: {used} reads counted, {skipped} skipped, {len(table)} barcodes with viral UMIs")


if __name__ == "__main__":
    main()
