from pathlib import Path
import argparse
import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from scipy.linalg import eigvals, svdvals

ROOT = Path(__file__).resolve().parents[1]
CURRENCY = {
    "ATP", "ADP", "AMP", "NAD", "NADH", "NADP", "NADPH",
    "FAD", "FADH2", "COA", "CO-A", "WATER", "H2O",
    "PHOSPHATE", "PI", "H+", "PROTON", "CO2"
}

def clean_name(value):
    return re.sub(r"\s+", " ", str(value).strip())

def is_currency(name):
    # BiGG SBML species identifiers commonly look like M_atp_c or atp_c.
    raw = str(name).strip().lower()
    raw = raw[2:] if raw.startswith("m_") else raw
    base = raw.rsplit("_", 1)[0] if raw.endswith(("_c", "_e", "_p")) else raw
    key = re.sub(r"[^A-Za-z0-9+]", "", base).upper()
    aliases = {
        "COENZYMEA": "COA", "COA": "COA", "H2O": "H2O",
        "WATER": "WATER", "PHOSPHATE": "PHOSPHATE",
        "INORGANICPHOSPHATE": "PHOSPHATE", "PROTON": "PROTON",
        "H": "H+", "CARBONDIOXIDE": "CO2", "ATPC": "ATP"
    }
    key = aliases.get(key, key)
    return key in CURRENCY or key in {
        "ATP", "ADP", "AMP", "NAD", "NADH", "NADP", "NADPH",
        "FAD", "FADH2", "COA", "H2O", "PI", "H", "CO2"
    }

def read_csv_reactions(path):
    df = pd.read_csv(path).fillna("")
    required = {"reaction_id", "substrates", "products", "reversible"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    reactions = []
    for _, row in df.iterrows():
        subs = [clean_name(x) for x in str(row["substrates"]).split(";") if clean_name(x)]
        prods = [clean_name(x) for x in str(row["products"]).split(";") if clean_name(x)]
        reactions.append((subs, prods, str(row["reversible"]).lower() in {"true", "1", "yes"}))
    return reactions

def read_sbml_reactions(path):
    root = ET.parse(path).getroot()
    reactions = []
    for rxn in (el for el in root.iter() if el.tag.split("}")[-1] == "reaction"):
        subs, prods = [], []
        for child in rxn:
            tag = child.tag.split("}")[-1]
            if tag == "listOfReactants":
                subs = [x.attrib.get("species", "") for x in child]
            elif tag == "listOfProducts":
                prods = [x.attrib.get("species", "") for x in child]
        reversible = rxn.attrib.get("reversible", "false").lower() == "true"
        if subs and prods:
            reactions.append((subs, prods, reversible))
    if not reactions:
        raise ValueError("No SBML reactions found. Check the downloaded model file.")
    return reactions

def build_edges(reactions, remove_currency=False):
    edges = set()
    for subs, prods, reversible in reactions:
        if remove_currency:
            subs = [s for s in subs if not is_currency(s)]
            prods = [p for p in prods if not is_currency(p)]
        for s in subs:
            for p in prods:
                if s != p:
                    edges.add((s, p))
                    if reversible:
                        edges.add((p, s))
    return sorted(edges)

def exact_controllability(adj, eig_tol=1e-7, rank_tol=1e-8):
    n = adj.shape[0]
    if n == 0:
        return 0, 0.0, None
    values = eigvals(adj.astype(float))
    # Cluster numerically close eigenvalues; conjugate eigenvalues remain separate.
    clusters = []
    for val in values:
        for cluster in clusters:
            if abs(val - cluster[0]) <= eig_tol * max(1.0, abs(val), abs(cluster[0])):
                cluster.append(val)
                break
        else:
            clusters.append([val])
    max_geom = 0
    max_eval = None
    for cluster in clusters:
        lam = sum(cluster) / len(cluster)
        singular = svdvals(lam * np.eye(n, dtype=complex) - adj.astype(complex))
        threshold = rank_tol * max(1.0, singular[0] if len(singular) else 1.0)
        rank = int(np.sum(singular > threshold))
        geom = n - rank
        if geom > max_geom:
            max_geom, max_eval = geom, lam
    return max_geom, max_geom / n, max_eval

def make_matrix(edges):
    nodes = sorted({x for edge in edges for x in edge})
    index = {node: i for i, node in enumerate(nodes)}
    A = np.zeros((len(nodes), len(nodes)), dtype=float)
    for source, target in edges:
        A[index[source], index[target]] = 1.0
    return nodes, A

def run(reactions, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    records = []
    edge_sets = {}
    for label, reduced in [("Full participant network", False), ("Without currency metabolites", True)]:
        edges = build_edges(reactions, remove_currency=reduced)
        nodes, A = make_matrix(edges)
        nd, fraction, dominant = exact_controllability(A)
        records.append({
            "network": label, "nodes_N": len(nodes), "edges_E": len(edges),
            "driver_nodes_ND": nd, "driver_fraction_nD": fraction,
            "dominant_eigenvalue_real": float(np.real(dominant)) if dominant is not None else np.nan,
            "dominant_eigenvalue_imag": float(np.imag(dominant)) if dominant is not None else np.nan,
            "eigenvalue_tolerance": 1e-7, "rank_tolerance": 1e-8
        })
        edge_sets[label] = (nodes, edges)
        pd.DataFrame(edges, columns=["source_metabolite", "target_metabolite"]).to_csv(
            outdir / ("metabolite_edges_full.csv" if not reduced else "metabolite_edges_without_currency.csv"),
            index=False
        )
    result = pd.DataFrame(records)
    result.to_csv(outdir / "controllability_comparison.csv", index=False)

    # Sensitivity to SVD rank tolerance on both network variants.
    sens = []
    for label, (nodes, edges) in edge_sets.items():
        _, A = make_matrix(edges)
        for tol in [1e-6, 1e-8, 1e-10, 1e-12]:
            nd, frac, _ = exact_controllability(A, rank_tol=tol)
            sens.append({"network": label, "rank_tolerance": tol, "ND": nd, "nD": frac})
    pd.DataFrame(sens).to_csv(outdir / "tolerance_sensitivity.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (label, (nodes, edges)) in zip(axes, edge_sets.items()):
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)
        pos = nx.spring_layout(G, seed=7)
        nx.draw_networkx(G, pos=pos, ax=ax, node_size=450, font_size=6,
                         arrows=True, with_labels=True, width=0.8)
        ax.set_title(label)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "network_comparison.png", dpi=220)
    plt.close(fig)
    print(result.to_string(index=False))
    print(f"\nFiles saved in: {outdir}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(ROOT / "data/ecoli_central_metabolism_subset.csv"))
    parser.add_argument("--sbml", default=None, help="Optional SBML model file; overrides --input")
    parser.add_argument("--output", default=str(ROOT / "results"))
    args = parser.parse_args()
    reactions = read_sbml_reactions(args.sbml) if args.sbml else read_csv_reactions(args.input)
    run(reactions, Path(args.output))

if __name__ == "__main__":
    main()
