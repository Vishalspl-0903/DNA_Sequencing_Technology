"""
demux_wick_by_barcode.py -- split a pooled tech_rep_*_reads.fastq.gz into one
FASTQ per isolate, using the barcode_arrangement column in the matching
sequencing_summary.txt.gz.

Verified against Wick's own method.md for this dataset (rrwick/Small-plasmid-
Nanopore): the summary file's column 21 is `barcode_arrangement`; reads
tagged `unclassified` or a barcode not in BARCODE_TO_ISOLATE are dropped, the
same way the original pipeline did (`rm barcode06.fastq.gz barcode09.fastq.gz`
in method.md means only 7 of 12 barcodes carry real isolates here).

Usage:
    python demux_wick_by_barcode.py \
        --reads Dataset/tech_rep_1_rapid_reads.fastq.gz \
        --summary Dataset/tech_rep_1_rapid_sequencing_summary.txt.gz \
        --out-dir Dataset/wick/tech_rep_1_rapid
"""
import argparse
import csv
import gzip
import os


# From method.md's final "assemblies/" renaming step -- confirmed mapping,
# not guessed. Barcodes 06 and 09 were dropped by the original authors too.
BARCODE_TO_ISOLATE = {
    "barcode01": "Acinetobacter_baumannii_J9",
    "barcode02": "Citrobacter_koseri_MINF_9D",
    "barcode03": "Enterobacter_kobei_MSB1_1B",
    "barcode04": "Haemophilus_unknown_M1C132_1",
    "barcode05": "Klebsiella_oxytoca_MSB1_2C",
    "barcode07": "Klebsiella_variicola_INF345",
    "barcode08": "Serratia_marcescens_17-147-1671",
}


def load_read_to_isolate(summary_path):
    """Map read_id -> isolate name (or None) from the sequencing summary."""
    opener = gzip.open if summary_path.endswith(".gz") else open
    read_to_isolate = {}
    with opener(summary_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "read_id" not in reader.fieldnames or "barcode_arrangement" not in reader.fieldnames:
            raise SystemExit(
                f"Expected columns 'read_id' and 'barcode_arrangement' in "
                f"{summary_path}, found: {reader.fieldnames}. Guppy/MinKNOW "
                f"versions occasionally rename these -- check the header "
                f"and adjust COLUMN names below if needed."
            )
        for row in reader:
            bc = row["barcode_arrangement"]
            read_to_isolate[row["read_id"]] = BARCODE_TO_ISOLATE.get(bc)
    return read_to_isolate


def demux(reads_path, summary_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    read_to_isolate = load_read_to_isolate(summary_path)
    print(f"loaded {len(read_to_isolate)} read->barcode assignments")

    handles = {}
    for isolate in set(BARCODE_TO_ISOLATE.values()):
        handles[isolate] = open(os.path.join(out_dir, f"{isolate}.fastq"), "w")

    counts = {isolate: 0 for isolate in handles}
    dropped = 0

    opener = gzip.open if reads_path.endswith(".gz") else open
    with opener(reads_path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            # ONT read IDs in the header look like "@<uuid> runid=... ..."
            read_id = header[1:].split()[0].strip()
            isolate = read_to_isolate.get(read_id)
            if isolate is None:
                dropped += 1
                continue
            handles[isolate].writelines([header, seq, plus, qual])
            counts[isolate] += 1

    for h in handles.values():
        h.close()

    print()
    for isolate, n in sorted(counts.items()):
        print(f"  {isolate:35s} {n:8d} reads")
    print(f"  {'(dropped: unclassified/other barcode)':35s} {dropped:8d} reads")
    print(f"\nwrote FASTQ files to {out_dir}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    demux(args.reads, args.summary, args.out_dir)
