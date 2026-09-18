"""
Phase 0 / step 5b -- GFA + labels -> PyTorch Geometric Data objects.

This is the object the GNN consumes, so it is built exactly as Sec. IV
specifies, including the parts that are easy to get wrong:

  * a VIRTUAL GLOBAL NODE connected bidirectionally to every segment. Without
    it, message passing over an isolated node aggregates across an empty
    neighbourhood and the network degenerates to an MLP on x_v -- for precisely
    the positive class. The global node is appended as the LAST node index and
    carries a flag feature so it is distinguishable from real segments.
  * edges are the UNDIRECTED SIMPLE graph: orientation discarded, multi-edges
    collapsed, SELF-LOOPS DROPPED (a self-loop is how a circular contig is
    represented -- counting it would reintroduce the circularity leak).
  * S(p) is stored as a pooling selector for the plasmid head, NOT as a node
    feature. No node carries any mark of which plasmid it supports.
  * y_node marks fragmented/absorbed support segments. Absence has no node, so
    it can never appear here; it lives only in y_graph.

Saved as a plain torch file so it can be loaded with or without a GPU.
"""
import os, csv, json, collections
import torch
from torch_geometric.data import Data

# PATCHED: the original hardcoded a Windows dev path
# (D:\Plasmid-GNN\DNA_Sequencing_Technology\work) that does not exist on any
# other machine and does not match this repo's actual layout (asm/, graphs/,
# results/ live at the repo root, not under a "work" subfolder). Derive ROOT
# from the script's own location instead, so this runs unmodified on macOS,
# Linux, or Windows, straight out of `git clone`.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATS = ["log10_length", "relative_depth", "delta_gc",
         "degree", "component_size", "n_components", "depth_x_length"]
POSITIVE_LABELS = ("fragmented", "absorbed", "absent")


def parse_links(gfa):
    edges = set()
    with open(gfa, "r", errors="replace") as fh:
        for line in fh:
            if line[:2] == "L\t":
                f = line.split("\t")
                a, b = f[1], f[3]
                if a != b:                       # self-loops dropped
                    edges.add(frozenset((a, b)))  # multi-edges collapsed
    return edges


def build(tag, node_feats):
    asm = os.path.join(ROOT, "asm", tag)
    gdir = os.path.join(ROOT, "graphs", tag)
    gfa = os.path.join(asm, "assembly_graph.gfa")
    lp = os.path.join(gdir, "plasmid_labels.csv")
    npth = os.path.join(gdir, "node_labels.csv")
    if not (os.path.exists(gfa) and os.path.exists(lp)):
        return None

    segs = sorted(node_feats.keys())
    idx = {s: i for i, s in enumerate(segs)}
    n = len(segs)
    if n == 0:
        return None

    # node features + one extra column flagging the virtual global node
    x = torch.zeros((n + 1, len(FEATS) + 1), dtype=torch.float)
    for s in segs:
        x[idx[s], :len(FEATS)] = torch.tensor(
            [node_feats[s][f] for f in FEATS], dtype=torch.float)
    x[n, len(FEATS)] = 1.0                       # global node flag

    src, dst = [], []
    for e in parse_links(gfa):
        a, b = tuple(e)
        if a in idx and b in idx:
            src += [idx[a], idx[b]]
            dst += [idx[b], idx[a]]
    # virtual global node <-> every segment, bidirectional
    for i in range(n):
        src += [n, i]
        dst += [i, n]
    edge_index = torch.tensor([src, dst], dtype=torch.long) if src else \
        torch.zeros((2, 0), dtype=torch.long)

    y_node = torch.zeros(n + 1, dtype=torch.float)
    if os.path.exists(npth):
        for r in csv.DictReader(open(npth)):
            if r["segment"] in idx and r["node_positive"] == "1":
                y_node[idx[r["segment"]]] = 1.0

    plasmids, y_plasmid, support_sets, gray = [], [], [], 0
    plasmid_labels4, plasmid_lens = [], []   # PATCHED: see note below
    n_absent = 0
    for r in csv.DictReader(open(lp)):
        if r["label"] == "gray_50_95":
            gray += 1
            continue
        if r["label"] == "absent":
            n_absent += 1
        plasmids.append(r["plasmid"])
        y_plasmid.append(1.0 if r["label"] in POSITIVE_LABELS else 0.0)
        sup = [idx[s] for s in r["support"].split(";") if s and s in idx]
        support_sets.append(sup)
        # PATCHED: 13_train_gnn.py's plasmid-pooled and node heads need the
        # full 4-way label (recovered/fragmented/absorbed/absent) and the
        # true plasmid length -- both already sit in plasmid_labels.csv but
        # were being dropped here, collapsed into the binary y_plasmid used
        # by the 06_baselines.py heads. Carry them through as parallel lists
        # instead of changing y_plasmid, so 06_baselines.py is unaffected.
        plasmid_labels4.append(r["label"])
        plasmid_lens.append(float(r["plasmid_len"]))

    d = Data(x=x, edge_index=edge_index)
    d.y_node = y_node
    d.node_mask = torch.tensor([1] * n + [0], dtype=torch.bool)   # exclude global
    d.y_plasmid = torch.tensor(y_plasmid, dtype=torch.float)
    d.y_graph = torch.tensor([float(n_absent)], dtype=torch.float)  # Poisson count
    d.support = support_sets                     # pooling selector, not a feature
    d.plasmid_names = plasmids
    d.plasmid_label4 = plasmid_labels4            # PATCHED: added
    d.plasmid_len = plasmid_lens                  # PATCHED: added
    d.global_index = n
    d.tag = tag
    d.isolate = tag.split("_")[0]
    d.n_segments = n
    d.gray_held_out = gray
    return d


def main():
    nf = collections.defaultdict(dict)
    fp = os.path.join(ROOT, "results", "node_features.csv")
    for r in csv.DictReader(open(fp)):
        nf[r["tag"]][r["segment"]] = {k: float(r[k]) for k in FEATS}

    graphs = []
    for tag in sorted(nf):
        d = build(tag, nf[tag])
        if d is None:
            continue
        graphs.append(d)
        empty = sum(1 for s in d.support if not s)
        print("  %-18s nodes=%-4d edges=%-5d plasmids=%-3d pos=%-2d "
              "absent=%-2d empty_S(p)=%d"
              % (tag, d.num_nodes, d.edge_index.shape[1], len(d.plasmid_names),
                 int(d.y_plasmid.sum()), int(d.y_graph.item()), empty))

    if not graphs:
        print("no graphs to export"); return
    out = os.path.join(ROOT, "results", "pyg_dataset.pt")
    torch.save(graphs, out)
    print()
    print("  graphs        :", len(graphs))
    print("  total nodes   :", sum(g.num_nodes for g in graphs))
    print("  total edges   :", sum(g.edge_index.shape[1] for g in graphs))
    print("  plasmid units :", sum(len(g.plasmid_names) for g in graphs))
    print("  node positives:", int(sum(g.y_node.sum() for g in graphs)))
    print("  isolates      :", len({g.isolate for g in graphs}))
    print("-> %s" % out)


if __name__ == "__main__":
    main()
