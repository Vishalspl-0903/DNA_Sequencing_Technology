"""
Phase 0 / step 1 -- turn the raw 6.9 GB NCBI plasmid query into a defensible
plasmid reference set for the simulation arm.

The raw file is NOT PLSDB: it is an unfiltered nucleotide query. It contains
1,581 records >=500 kb (up to 11.85 Mb) which the proposal's own appendix rule
(contigs >=500 kb are chromosome) would classify as chromosomes. Feeding those
to Badread would silently simulate chromosomes as plasmids.

Pass 1 (this script) streams the raw file once and writes a manifest of every
record with the fields needed to select on. No sequence is written here.

Filters, all fixed before any selection:
  F1  header declares 'plasmid'
  F2  header declares 'complete sequence'   (drop partial/draft records)
  F3  1,000 bp <= length < 500,000 bp       (matches the >=500 kb chromosome rule)
  F4  ambiguous-base fraction <= 0.1%       (Badread needs clean templates)
  F5  first occurrence of an accession wins (deduplicate)
  F6  first occurrence of an identical sequence wins (deduplicate replicons)
"""
import os, csv, hashlib, collections

RAW = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\sequences.fasta"
OUT = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work\refs\plasmid_manifest.csv"

MIN_LEN, MAX_LEN = 1_000, 500_000
MAX_N_FRAC = 0.001

ACGT = set(b"ACGTacgt")


def records(path):
    """Stream (header, seq_bytes) without holding the file in memory."""
    hdr, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if hdr is not None:
                    yield hdr, b"".join(chunks)
                hdr = line[1:].decode("utf-8", "replace").rstrip()
                chunks = []
            else:
                chunks.append(line.strip())
    if hdr is not None:
        yield hdr, b"".join(chunks)


def parse_header(h):
    parts = h.split()
    acc = parts[0] if parts else ""
    organism = " ".join(parts[1:3]) if len(parts) >= 3 else ""
    genus = parts[1] if len(parts) >= 2 else ""
    return acc, genus, organism


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    seen_acc, seen_seq = set(), set()
    reasons = collections.Counter()
    kept = 0
    total = 0

    with open(OUT, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo)
        w.writerow(["accession", "genus", "organism", "length", "gc",
                    "n_frac", "verdict", "reason"])

        for hdr, seq in records(RAW):
            total += 1
            acc, genus, organism = parse_header(hdr)
            low = hdr.lower()
            L = len(seq)

            gc = (seq.count(b"G") + seq.count(b"C") +
                  seq.count(b"g") + seq.count(b"c"))
            n_amb = L - sum(seq.count(bytes([c])) for c in b"ACGTacgt")
            n_frac = (n_amb / L) if L else 1.0
            gc_frac = (gc / L) if L else 0.0

            reason = ""
            if "plasmid" not in low:
                reason = "F1_not_plasmid"
            elif "complete sequence" not in low:
                reason = "F2_not_complete"
            elif L < MIN_LEN:
                reason = "F3_too_short"
            elif L >= MAX_LEN:
                reason = "F3_chromosome_sized"
            elif n_frac > MAX_N_FRAC:
                reason = "F4_ambiguous_bases"
            elif acc in seen_acc:
                reason = "F5_dup_accession"
            else:
                h = hashlib.md5(seq.upper()).hexdigest()
                if h in seen_seq:
                    reason = "F6_dup_sequence"
                else:
                    seen_seq.add(h)

            if reason:
                reasons[reason] += 1
                verdict = "reject"
            else:
                seen_acc.add(acc)
                kept += 1
                verdict = "keep"

            w.writerow([acc, genus, organism, L, "%.4f" % gc_frac,
                        "%.6f" % n_frac, verdict, reason])

            if total % 10000 == 0:
                print("  ...%d scanned, %d kept" % (total, kept), flush=True)

    print()
    print("=" * 62)
    print("RAW records scanned : %d" % total)
    print("KEPT                : %d" % kept)
    print("REJECTED            : %d" % (total - kept))
    for r, c in reasons.most_common():
        print("    %-22s %d" % (r, c))
    print("manifest -> %s" % OUT)


if __name__ == "__main__":
    main()
