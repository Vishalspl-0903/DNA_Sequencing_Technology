"""
Phase 0 / step 5 -- node features from the GFA, per Sec. IV + Appendix.

Six node features, all computable at inference from long reads alone:
  log10_length, relative_depth, delta_gc, degree, component stats, depth x length

Appendix rules implemented literally:
  * chromosome set = all segments >=500 kb; reference depth is their
    LENGTH-WEIGHTED MEDIAN depth, reference GC their LENGTH-WEIGHTED MEAN.
  * if no segment reaches 500 kb the isolate falls back to the length-weighted
    median over all segments and is FLAGGED in the manifest (flagged isolates
    are reported separately, never pooled).
  * degree is computed on the undirected SIMPLE graph: orientation signs
    discarded, multi-edges collapsed, SELF-LOOPS DROPPED. A circular contig is
    represented with a self-loop, so counting them would reintroduce the
    circularity leak through the graph instead of the contig report.
  * self-loop count is retained as a manifest diagnostic and NEVER as a feature.

Flye's circularity flag is likewise excluded from X and kept only as a
diagnostic -- the label definition itself invokes circularity.
"""
import os, csv, json, collections, statistics

ROOT = r"D:\DNA-Sequencing\work"
CHROM_MIN = 500_000


def parse_gfa(path):
    segs, raw_links, selfloops = {}, [], 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                name, seq = f[1], f[2]
                depth = None
                for tag in f[3:]:
                    if tag[:5] in ("dp:i:", "dp:f:", "DP:f:"):
                        depth = float(tag[5:])
                gc = 0
                if seq != "*":
                    u = seq.upper()
                    gc = (u.count("G") + u.count("C")) / len(u) if u else 0.0
                segs[name] = {"len": len(seq) if seq != "*" else 0,
                              "depth": depth if depth is not None else 0.0,
                              "gc": gc}
            elif f[0] == "L":
                a, b = f[1], f[3]
                if a == b:
                    selfloops += 1
                else:
                    raw_links.append(frozenset((a, b)))
    simple_edges = set(raw_links)          # multi-edges collapsed
    return segs, simple_edges, selfloops


def weighted_median(pairs):
    """pairs = [(value, weight)] -> length-weighted median."""
    pairs = sorted(pairs)
    tot = sum(w for _, w in pairs)
    if tot == 0:
        return 0.0
    acc = 0
    for v, w in pairs:
        acc += w
        if acc >= tot / 2:
            return v
    return pairs[-1][0]


def components(nodes, edges):
    adj = collections.defaultdict(set)
    for e in edges:
        a, b = tuple(e)
        adj[a].add(b); adj[b].add(a)
    seen, comp_of, sizes = set(), {}, []
    for n in nodes:
        if n in seen:
            continue
        stack, members = [n], []
        seen.add(n)
        while stack:
            x = stack.pop(); members.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        cid = len(sizes)
        for m in members:
            comp_of[m] = cid
        sizes.append(len(members))
    return comp_of, sizes


def features_for(tag):
    gfa = os.path.join(ROOT, "asm", tag, "assembly_graph.gfa")
    if not os.path.exists(gfa):
        return None
    segs, edges, selfloops = parse_gfa(gfa)
    if not segs:
        return None

    chrom = [n for n, d in segs.items() if d["len"] >= CHROM_MIN]
    fallback = False
    if chrom:
        ref_depth = weighted_median([(segs[n]["depth"], segs[n]["len"]) for n in chrom])
        tw = sum(segs[n]["len"] for n in chrom)
        ref_gc = sum(segs[n]["gc"] * segs[n]["len"] for n in chrom) / tw
    else:
        fallback = True
        ref_depth = weighted_median([(d["depth"], d["len"]) for d in segs.values()])
        tw = sum(d["len"] for d in segs.values()) or 1
        ref_gc = sum(d["gc"] * d["len"] for d in segs.values()) / tw
    if ref_depth <= 0:
        ref_depth = 1.0

    deg = collections.Counter()
    for e in edges:
        a, b = tuple(e)
        deg[a] += 1; deg[b] += 1
    comp_of, sizes = components(list(segs), edges)

    import math
    rows = []
    for n, d in segs.items():
        rel_depth = d["depth"] / ref_depth
        log_len = math.log10(d["len"]) if d["len"] > 0 else 0.0
        rows.append({
            "tag": tag, "segment": n,
            "log10_length": round(log_len, 5),
            "relative_depth": round(rel_depth, 5),
            "delta_gc": round(d["gc"] - ref_gc, 5),
            "degree": deg.get(n, 0),
            "component_size": sizes[comp_of[n]],
            "n_components": len(sizes),
            "depth_x_length": round(rel_depth * log_len, 5),
            # diagnostics only, never features:
            "_raw_len": d["len"], "_raw_depth": d["depth"], "_gc": round(d["gc"], 5),
        })

    diag = {"tag": tag, "n_segments": len(segs), "n_simple_edges": len(edges),
            "self_loops": selfloops, "n_components": len(sizes),
            "chrom_set_size": len(chrom), "chrom_fallback": fallback,
            "ref_depth": round(ref_depth, 3), "ref_gc": round(ref_gc, 5),
            "isolated_nodes": sum(1 for n in segs if deg.get(n, 0) == 0)}
    return rows, diag


def main():
    asm_root = os.path.join(ROOT, "asm")
    tags = sorted(d for d in os.listdir(asm_root)
                  if os.path.isdir(os.path.join(asm_root, d)))
    allrows, diags = [], []
    for t in tags:
        r = features_for(t)
        if not r:
            print("  %s: no graph" % t); continue
        rows, diag = r
        allrows.extend(rows); diags.append(diag)
        print("  %-14s segments=%3d edges=%3d selfloops=%3d comps=%3d "
              "isolated=%3d chromset=%d%s"
              % (t, diag["n_segments"], diag["n_simple_edges"], diag["self_loops"],
                 diag["n_components"], diag["isolated_nodes"], diag["chrom_set_size"],
                 "  [FALLBACK]" if diag["chrom_fallback"] else ""))

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    fp = os.path.join(ROOT, "results", "node_features.csv")
    with open(fp, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(allrows[0].keys()))
        w.writeheader(); w.writerows(allrows)
    json.dump(diags, open(os.path.join(ROOT, "results", "graph_diagnostics.json"), "w"),
              indent=2)
    print("\n-> %s  (%d nodes over %d graphs)" % (fp, len(allrows), len(diags)))
    print("-> graph_diagnostics.json")


if __name__ == "__main__":
    main()
