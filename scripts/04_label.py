"""
Phase 0 / step 4 -- label construction from an assembly graph.

Implements Sec. III of the proposal exactly:
  * reference plasmid is the QUERY, GFA segments are the TARGET
  * minimap2 -x asm5, secondary alignments retained
  * identity threshold applied PER ALIGNMENT BLOCK, before merging
  * alignment intervals projected onto the plasmid and MERGED before any
    coverage is computed (never summed -- that would double-count repeats and
    push coverage past 100%)

Label thresholds (fixed before modelling):
  recovered  : >=95% of the plasmid to a SINGLE contig at >=95% identity
  fragmented : >=95% cumulative across >=2 contigs
  absorbed   : containment within a contig of >= 3x the plasmid length
  absent     : <50% anywhere
  50-95%     : held out, reported separately

IMPLEMENTATION DECISION (flagged, not silently absorbed): the proposal does not
state precedence between 'recovered' and 'absorbed', and a plasmid can satisfy
both (>=95% into one contig that is also >=3x its length). We treat absorbed as
a refinement of recovered -- if the single covering contig is >=3x the plasmid
length the event is absorption, not recovery -- because physically the plasmid
has been merged into a larger replicon. This needs confirming before Phase 0
results are interpreted.

Support set S(p): segment enters if its merged coverage of p is >=10% of the
plasmid length AND >=500 bp, at >=95% identity. Deliberately looser than the
label thresholds so fragmented plasmids keep all their pieces.
"""
import os, sys, csv, json, subprocess, collections

MIN_IDENT      = 0.95
REC_COV        = 0.95
ABSENT_COV     = 0.50
ABSORB_RATIO   = 3.0
SUP_FRAC       = 0.10
SUP_BP         = 500


# ---------------------------------------------------------------- fasta / gfa
def read_fasta(path):
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


def parse_gfa(path):
    """Return segments {name: {'len','seq','depth'}} and raw link list."""
    segs, links = {}, []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f:
                continue
            if f[0] == "S":
                name, seq = f[1], f[2]
                depth = None
                for tag in f[3:]:
                    if tag.startswith("dp:i:"):
                        depth = float(tag[5:])
                    elif tag.startswith("DP:f:"):
                        depth = float(tag[5:])
                    elif tag.startswith("dp:f:"):
                        depth = float(tag[5:])
                L = len(seq) if seq != "*" else 0
                segs[name] = {"len": L, "seq": seq, "depth": depth}
            elif f[0] == "L":
                links.append((f[1], f[3]))
    return segs, links


def merge(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def total(iv):
    return sum(e - s for s, e in iv)


# ---------------------------------------------------------------- main
def label_one(asm_dir, ref_fasta, out_dir, mm2):
    gfa = os.path.join(asm_dir, "assembly_graph.gfa")
    if not os.path.exists(gfa):
        print("  no GFA in %s -- skipping" % asm_dir)
        return None
    os.makedirs(out_dir, exist_ok=True)

    segs, links = parse_gfa(gfa)

    # segments -> fasta (alignment target)
    seg_fa = os.path.join(out_dir, "segments.fasta")
    with open(seg_fa, "w") as fo:
        for n, d in segs.items():
            if d["len"]:
                fo.write(">%s\n%s\n" % (n, d["seq"]))

    # reference plasmids (query) = every record of the sim reference bar the chromosome
    pl_fa = os.path.join(out_dir, "plasmids.fasta")
    plen = {}
    with open(pl_fa, "w") as fo:
        for hdr, seq in read_fasta(ref_fasta):
            name = hdr.split()[0]
            if name == "chromosome":
                continue
            plen[name] = len(seq)
            fo.write(">%s\n%s\n" % (name, seq.decode()))

    if not plen:
        print("  no plasmids in %s" % ref_fasta)
        return None

    # minimap2: target=segments, query=plasmids; keep secondaries
    paf = os.path.join(out_dir, "aln.paf")
    cmd = [mm2, "-x", "asm5", "-c", "--secondary=yes", "-N", "50",
           "-t", "8", seg_fa, pl_fa]
    with open(paf, "w") as fo, open(os.path.join(out_dir, "mm2.log"), "w") as fe:
        subprocess.run(cmd, stdout=fo, stderr=fe, check=True)

    # PAF: qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq
    per_seg = collections.defaultdict(lambda: collections.defaultdict(list))
    nblocks = nkept = 0
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12:
                continue
            nblocks += 1
            q, qs, qe = f[0], int(f[2]), int(f[3])
            t = f[5]
            nmatch, alen = int(f[9]), int(f[10])
            if alen == 0 or (nmatch / alen) < MIN_IDENT:   # identity per block
                continue
            nkept += 1
            per_seg[q][t].append((qs, qe))

    rows, node_pos = [], collections.defaultdict(set)
    for p, L in sorted(plen.items()):
        seg_cov = {}
        for t, iv in per_seg.get(p, {}).items():
            seg_cov[t] = total(merge(iv))
        all_iv = [x for iv in per_seg.get(p, {}).values() for x in iv]
        cum = total(merge(all_iv))
        cum_frac = cum / L

        support = [t for t, c in seg_cov.items()
                   if c >= SUP_FRAC * L and c >= SUP_BP]

        best_t, best_c = None, 0
        for t, c in seg_cov.items():
            if c > best_c:
                best_t, best_c = t, c
        best_frac = best_c / L if L else 0.0

        if best_frac >= REC_COV:
            tlen = segs[best_t]["len"]
            if tlen >= ABSORB_RATIO * L:
                lab = "absorbed"
            else:
                lab = "recovered"
        elif cum_frac >= REC_COV and len(support) >= 2:
            lab = "fragmented"
        elif cum_frac < ABSENT_COV:
            lab = "absent"
        else:
            lab = "gray_50_95"

        if lab in ("fragmented", "absorbed"):
            for t in support:
                node_pos[t].add(lab)

        rows.append({
            "plasmid": p, "plasmid_len": L, "label": lab,
            "best_segment": best_t or "", "best_cov_frac": round(best_frac, 4),
            "best_segment_len": segs[best_t]["len"] if best_t else 0,
            "cumulative_cov_frac": round(cum_frac, 4),
            "n_support": len(support), "support": ";".join(sorted(support)),
        })

    with open(os.path.join(out_dir, "plasmid_labels.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    with open(os.path.join(out_dir, "node_labels.csv"), "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["segment", "seg_len", "node_positive", "mechanisms"])
        for n, d in segs.items():
            m = sorted(node_pos.get(n, []))
            w.writerow([n, d["len"], 1 if m else 0, ";".join(m)])

    counts = collections.Counter(r["label"] for r in rows)
    summary = {
        "asm_dir": asm_dir, "n_segments": len(segs), "n_links": len(links),
        "n_plasmids": len(rows), "labels": dict(counts),
        "n_missing_for_graph_head": counts["absent"],
        "paf_blocks": nblocks, "paf_blocks_passing_identity": nkept,
    }
    json.dump(summary, open(os.path.join(out_dir, "label_summary.json"), "w"),
              indent=2)
    return summary


def default_root():
    """minimap2 lives in the WSL env, so this script normally runs under WSL
    against /mnt/d; keep the Windows path working for local inspection."""
    return os.environ.get(
        "PIPE_ROOT",
        r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work" if os.name == "nt"
        else "/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work")


if __name__ == "__main__":
    root = default_root()
    mm2 = sys.argv[1] if len(sys.argv) > 1 else "minimap2"
    asm_root = os.path.join(root, "asm")
    out_root = os.path.join(root, "graphs")
    todo = sorted(d for d in os.listdir(asm_root)
                  if os.path.isdir(os.path.join(asm_root, d)))
    allsum = []
    for tag in todo:
        iso = tag.split("_")[0]          # sim01_lig_d30 -> sim01
        ref = os.path.join(root, "sim", iso + "_ref.fasta")
        if not os.path.exists(ref):
            continue
        print("labelling", tag)
        s = label_one(os.path.join(asm_root, tag), ref,
                      os.path.join(out_root, tag), mm2)
        if s:
            s["tag"] = tag
            allsum.append(s)
            print("   ", s["labels"], "segments=%d" % s["n_segments"])
    json.dump(allsum, open(os.path.join(out_root, "all_label_summaries.json"), "w"),
              indent=2)
    print("\nlabelled %d assemblies" % len(allsum))
