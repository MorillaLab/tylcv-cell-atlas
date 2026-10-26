#!/usr/bin/env python3
"""make_tylcv_reference.py - TYLCV reference files for Space Ranger and the audit notebook.

Default source: GenBank AJ489258.1, TYLCV-IL [ES:Alm:Pep:99], the Almeria isolate described by
Morilla et al. 2005 (IHSM). If your infectious clone was sequenced, use it instead (--genbank / --fasta):
reads that differ from the reference by more than a few percent map poorly.

Writes into --outdir:
  TYLCV_IL.fa                 one record named 'TYLCV', virion-sense orientation
  tylcv_spaceranger.gtf       two strand-specific features spanning the genome: TYLCV_V (+), TYLCV_C (-).
                              Use THIS one for `spaceranger mkref` (overlapping ORFs would make reads ambiguous).
  tylcv_orfs.gtf              V2, V1 (CP), C1 (Rep), C2 (TrAP), C3 (REn), C4 from the GenBank CDS features.
                              Annotation/visualisation only - do not use it for counting.
  tylcv_reference_report.txt  provenance (accession, definition, length, MD5) and checks

Standard library only. Examples:
  python3 make_tylcv_reference.py                         # downloads AJ489258.1 from NCBI
  python3 make_tylcv_reference.py --accession NC_004005.1 # the NCBI RefSeq for TYLCV
  python3 make_tylcv_reference.py --genbank my_clone.gb   # your own annotated clone
  python3 make_tylcv_reference.py --fasta my_clone.fa     # FASTA only (no ORF GTF)
"""
import argparse
import hashlib
import re
import sys
import urllib.request
from pathlib import Path

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text"
COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")
ORF_RULES = [  # (canonical name, regex on /gene, /product, /note, /locus_tag)
    ("V2", r"\b(A?V2|pre-?coat|precoat|AV2)\b"),
    ("V1", r"\b(A?V1|coat protein|CP)\b"),
    ("C1", r"\b(A?C1|Rep|replication[- ]associated)\b"),
    ("C2", r"\b(A?C2|TrAP|transcriptional activator)\b"),
    ("C3", r"\b(A?C3|REn|replication enhancer)\b"),
    ("C4", r"\b(A?C4)\b"),
]


def revcomp(s):
    return s.translate(COMP)[::-1]


def download(acc):
    with urllib.request.urlopen(EFETCH.format(acc=acc), timeout=60) as r:
        text = r.read().decode()
    if "LOCUS" not in text[:200]:
        sys.exit(f"NCBI did not return a GenBank record for {acc}:\n{text[:300]}")
    return text


def parse_genbank(text):
    rec = {"locus_len": None, "definition": "", "accession": "", "features": [], "seq": ""}
    m = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", text, re.M)
    rec["locus_len"] = int(m.group(1)) if m else None
    m = re.search(r"^DEFINITION\s+(.+?)\n(?=\S)", text, re.M | re.S)
    rec["definition"] = " ".join(m.group(1).split()) if m else ""
    m = re.search(r"^VERSION\s+(\S+)", text, re.M) or re.search(r"^ACCESSION\s+(\S+)", text, re.M)
    rec["accession"] = m.group(1) if m else ""
    origin = text.split("\nORIGIN", 1)
    if len(origin) == 2:
        rec["seq"] = re.sub(r"[^A-Za-z]", "", origin[1].split("//")[0]).upper()
    feats = text.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0] if "\nFEATURES" in text else ""
    cur = None
    for line in feats.splitlines()[1:]:
        key, body = line[5:21].strip(), line[21:].rstrip()
        if key:                                  # new feature
            cur = {"type": key, "location": body, "quals": {}}
            rec["features"].append(cur)
            last_q = None
        elif cur is not None:
            if body.startswith("/"):
                q, _, v = body[1:].partition("=")
                cur["quals"][q] = v.strip('"')
                last_q = q
            elif last_q is None:                 # location continues on the next line
                cur["location"] += body.strip()
            else:                                # qualifier value continues
                cur["quals"][last_q] = (cur["quals"][last_q] + " " + body.strip().strip('"')).strip()
    return rec


def parse_location(loc):
    strand = "-" if "complement" in loc else "+"
    parts = [(int(a), int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", loc)]
    return strand, parts


def name_cds(quals):
    text = " ".join(quals.get(k, "") for k in ("gene", "product", "note", "locus_tag"))
    for name, rx in ORF_RULES:
        if re.search(rx, text, re.I):
            return name
    return quals.get("gene") or quals.get("product", "orf").replace(" ", "_")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--accession", default="AJ489258.1")
    src.add_argument("--genbank")
    src.add_argument("--fasta")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--name", default="TYLCV", help="contig name used in FASTA and GTF")
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    report = []

    if a.fasta:
        lines = Path(a.fasta).read_text().splitlines()
        headers = [l for l in lines if l.startswith(">")]
        if len(headers) != 1:
            sys.exit(f"{a.fasta} must contain exactly one record (found {len(headers)})")
        rec = {"accession": headers[0][1:].split()[0], "definition": headers[0][1:], "features": [],
               "seq": "".join(l.strip() for l in lines if not l.startswith(">")).upper(), "locus_len": None}
        source = f"local FASTA {a.fasta}"
    else:
        text = Path(a.genbank).read_text() if a.genbank else download(a.accession)
        rec = parse_genbank(text)
        source = f"local GenBank {a.genbank}" if a.genbank else f"NCBI efetch {a.accession}"
        (out / f"{(rec['accession'] or 'tylcv').replace('/', '_')}.gb").write_text(text)

    seq = re.sub(r"[^ACGTRYKMSWBDHVN]", "", rec["seq"])
    if not seq:
        sys.exit("No sequence found.")
    if rec["locus_len"] and rec["locus_len"] != len(seq):
        sys.exit(f"Length mismatch: LOCUS says {rec['locus_len']} bp, sequence has {len(seq)}")
    L = len(seq)
    report += [f"source: {source}", f"accession: {rec['accession']}", f"definition: {rec['definition']}",
               f"length: {L} nt", f"GC: {100 * (seq.count('G') + seq.count('C')) / L:.1f}%"]
    if not 2600 <= L <= 2900:
        report.append(f"WARNING: {L} nt is unusual for a monopartite begomovirus (~2.78 kb).")

    cds = []
    for f in rec["features"]:
        if f["type"] != "CDS":
            continue
        strand, parts = parse_location(f["location"])
        if parts:
            cds.append({"name": name_cds(f["quals"]), "strand": strand, "parts": parts,
                        "product": f["quals"].get("product", "")})

    flipped = False
    v_strands = {c["strand"] for c in cds if c["name"] in ("V1", "V2")}
    if v_strands == {"-"}:                       # record is complementary-sense: flip to virion sense
        seq = revcomp(seq)
        for c in cds:
            c["strand"] = "+" if c["strand"] == "-" else "-"
            c["parts"] = sorted((L - e + 1, L - s + 1) for s, e in c["parts"])
        flipped = True
    report.append("orientation: " + ("reverse-complemented to virion sense" if flipped else
                                     "virion sense (V1/V2 on +)" if v_strands == {"+"} else
                                     "UNVERIFIED (no V1/V2 annotation) - check that the record is virion sense"))

    circ = seq + seq[:8]
    nona = [m.start() + 1 for m in re.finditer("TAATATTAC", circ) if m.start() < L]
    report.append(f"conserved nonanucleotide TAATATTAC (origin of replication) at: {nona or 'not found'}")

    for c in cds:
        nt = "".join(seq[s - 1:e] for s, e in c["parts"])
        nt = nt if c["strand"] == "+" else revcomp(nt)
        c["aa"] = (len(nt) // 3) - 1
        c["start_ok"] = nt[:3] == "ATG"
        c["stop_ok"] = nt[-3:] in ("TAA", "TAG", "TGA")
    expected = {"V1": "+", "V2": "+", "C1": "-", "C2": "-", "C3": "-", "C4": "-"}
    for c in cds:
        flag = "" if expected.get(c["name"], c["strand"]) == c["strand"] else "  <-- unexpected strand"
        report.append(f"CDS {c['name']:4s} {c['strand']} {','.join(f'{s}..{e}' for s, e in c['parts'])}  "
                      f"{c['aa']} aa  ATG={c['start_ok']} stop={c['stop_ok']}  {c['product']}{flag}")
    if not cds:
        report.append("no CDS features: tylcv_orfs.gtf not written (FASTA input or unannotated record)")

    md5 = hashlib.md5(seq.encode()).hexdigest()
    report.append(f"sequence MD5 (as written, uppercase): {md5}")
    fasta = f">{a.name} {rec['accession']} {rec['definition']}\n" + "\n".join(seq[i:i + 70] for i in range(0, L, 70)) + "\n"
    (out / "TYLCV_IL.fa").write_text(fasta)

    attrs = lambda g: f'gene_id "{g}"; transcript_id "{g}.1"; gene_name "{g}"; gene_biotype "protein_coding";'
    sr = (f"{a.name}\tcustom\texon\t1\t{L}\t.\t+\t.\t{attrs(a.name + '_V')}\n"
          f"{a.name}\tcustom\texon\t1\t{L}\t.\t-\t.\t{attrs(a.name + '_C')}\n")
    (out / "tylcv_spaceranger.gtf").write_text(sr)

    if cds:
        lines = []
        for c in sorted(cds, key=lambda c: (c["strand"], c["parts"][0][0])):
            g = f"{a.name}_{c['name']}"
            extra = f' product "{c["product"]}";' if c["product"] else ""
            lo, hi = min(s for s, _ in c["parts"]), max(e for _, e in c["parts"])
            lines.append(f"{a.name}\tGenBank\ttranscript\t{lo}\t{hi}\t.\t{c['strand']}\t.\t{attrs(g)}{extra}")
            for s, e in c["parts"]:
                lines.append(f"{a.name}\tGenBank\texon\t{s}\t{e}\t.\t{c['strand']}\t.\t{attrs(g)}{extra}")
                lines.append(f"{a.name}\tGenBank\tCDS\t{s}\t{e}\t.\t{c['strand']}\t0\t{attrs(g)}{extra}")
        (out / "tylcv_orfs.gtf").write_text("\n".join(lines) + "\n")

    (out / "tylcv_reference_report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    print(f"\nwritten to {out.resolve()}: TYLCV_IL.fa, tylcv_spaceranger.gtf" + (", tylcv_orfs.gtf" if cds else ""))


if __name__ == "__main__":
    main()
