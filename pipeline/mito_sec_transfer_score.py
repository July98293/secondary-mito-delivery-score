

import os, sys, io, warnings
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

import h5py
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT  = os.path.join(ROOT, "otput")
CACHE_PATH = os.path.join(ROOT, "_gene_expr_cache.npz")
H5AD_PATH  = os.path.join(DATA, "immuno-data.h5ad")

os.makedirs(OUT, exist_ok=True)

MODULES = {
    "Exosome":    {
        "weight": 0.20, "role": "functional",
        "genes": ["TSG101","PDCD6IP","VPS4A","VPS4B","CD63","CD9","CD81","SDCBP"],
    },
    "Mito Health":{
        "weight": 0.18, "role": "functional",
        "genes": ["TFAM","OPA1","MFN1","MFN2","DNM1L","PINK1","PRKN","PPARGC1A","ATP5F1A"],
    },
    "Rab Traffic":{
        "weight": 0.12, "role": "functional",
        "genes": ["RAB27A","RAB27B","RAB11A","RAB11B","RAB35","RAB5A","RAB7A"],
    },
    "MDV":        {
        "weight": 0.10, "role": "functional",
        "genes": ["VPS35","SNX9","ATG9A","MUL1","RAB9A","OPTN"],
    },
    "TNT":        {
        "weight": 0.08, "role": "functional",
        "genes": ["RHOA","CDC42","ACTA2","INF2","FMNL1","FMNL2","ARPC2","ARPC3"],
    },
    "Adhesion":   {
        "weight": 0.08, "role": "functional",
        "genes": ["CXCR4","CCR7","ITGA4","ITGB1","ICAM1"],
    },
    "Metabolism": {
        "weight": 0.08, "role": "functional",
        "genes": ["PFKFB3","SLC2A1","LDHA","CPT1A","HIF1A"],
    },
    "MV Form.":   {
        "weight": 0.08, "role": "functional",
        "genes": ["ARF6","ROCK1","ROCK2","MYH9","ARRDC1","FLOT1","FLOT2","PLD2"],
    },
    "Shedding":   {
        "weight": None, "role": "bonus+",
        "genes": ["CD63","LAMP1","TSG101"],
    },
    "Activation": {
        "weight": None, "role": "penalty-",
        "genes": ["IL2RA","CD69","TNF","IL1B","IFNG","HLA-DRA"],
    },
    "Viability":  {
        "weight": None, "role": "viability",
        "anti": ["BCL2","MCL1","BIRC2","XIAP"],
        "pro":  ["BAX","BAK1"],
        "genes": ["BCL2","MCL1","BIRC2","XIAP","BAX","BAK1"],
    },
}

FUNCTIONAL_KEYS = ["Exosome","Mito Health","Rab Traffic","MDV","TNT",
                   "Adhesion","Metabolism","MV Form."]
FEATURE_KEYS    = FUNCTIONAL_KEYS + ["Shedding","Activation","Viability"]

ALL_GENES = sorted({g for m in MODULES.values() for g in m["genes"]})
ALIASES   = {"HLA-DRA": ["HLA-DRA","HLADRA","HLA_DRA"]}

def read_categorical(grp):
    cats  = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["categories"][:]]
    codes = grp["codes"][:]
    return [cats[c] if 0 <= c < len(cats) else "NA" for c in codes]


def load_expression():
    if os.path.exists(CACHE_PATH):
        print(f"  NPZ cache: {CACHE_PATH}")
        c = np.load(CACHE_PATH, allow_pickle=True)
        expr         = c["expr"]
        cell_types   = list(c["cell_types"])
        tissues      = list(c["tissues"])
        cached_genes = list(c["gene_cols"])
        cached_upper = {g.upper(): i for i, g in enumerate(cached_genes)}

        gene_idx, missing = {}, []
        for g in ALL_GENES:
            for variant in [g.upper(), g.upper().replace("-","_"), g.upper().replace("-","")]:
                if variant in cached_upper:
                    gene_idx[g] = cached_upper[variant]; break
            else:
                found = False
                for alias in ALIASES.get(g, []):
                    if alias.upper() in cached_upper:
                        gene_idx[g] = cached_upper[alias.upper()]; found = True; break
                if not found:
                    missing.append(g)

        if missing:
            print(f"  Missing from cache ({len(missing)}): {missing} — reading h5ad …")
            extra_expr, extra_names, cell_types, tissues = _read_h5ad(missing, cell_types, tissues)
            start = expr.shape[1]
            for i, g in enumerate(extra_names):
                gene_idx[g] = start + i
            expr = np.hstack([expr, extra_expr])
        else:
            print(f"  All {len(ALL_GENES)} genes in cache")

        return expr, np.array(cell_types), np.array(tissues), gene_idx

    print("  No cache — reading h5ad …")
    expr, found_names, cell_types, tissues = _read_h5ad(ALL_GENES, None, None)
    return expr, np.array(cell_types), np.array(tissues), {g: i for i, g in enumerate(found_names)}


def _read_h5ad(target_genes, cell_types_in, tissues_in):
    CHUNK = 20_000
    with h5py.File(H5AD_PATH, "r") as f:
        if cell_types_in is None:
            cell_types_in = read_categorical(f["obs"]["cell_type"])
            tissues_in    = read_categorical(f["obs"]["tissue"])

        var_grp = f["var"]
        raw = var_grp["_index"][:] if "_index" in var_grp else var_grp[list(var_grp.keys())[0]][:]
        gene_names = [x.decode() if isinstance(x, bytes) else str(x) for x in raw]
        gene_upper = {g.upper(): i for i, g in enumerate(gene_names)}

        target_cols, found_names = [], []
        for g in target_genes:
            for c in [g.upper()] + [a.upper() for a in ALIASES.get(g, [])]:
                if c in gene_upper:
                    target_cols.append(gene_upper[c]); found_names.append(g); break

        target_col_arr = np.array(target_cols, dtype=np.int32)
        n_cells  = len(cell_types_in)
        expr_out = np.zeros((n_cells, len(target_cols)), dtype=np.float32)

        try:
            total_counts = f["obs"]["total_counts"][:].astype(np.float32)
        except Exception:
            total_counts = None

        lyr = f["layers"].get("decontXcounts", f["layers"].get("counts"))
        if lyr is None:
            raise RuntimeError("No counts layer in h5ad")
        indptr = lyr["indptr"][:]

        print(f"  Extracting {len(target_cols)} genes from {n_cells:,} cells …")
        for start in range(0, n_cells, CHUNK):
            end = min(start + CHUNK, n_cells)
            chunk_idx  = lyr["indices"][indptr[start]:indptr[end]].astype(np.int32)
            chunk_data = lyr["data"][indptr[start]:indptr[end]].astype(np.float32)
            row_lens   = np.diff(indptr[start:end+1])
            row_ids    = np.repeat(np.arange(end - start), row_lens)
            mask = np.isin(chunk_idx, target_col_arr)
            if mask.any():
                np.add.at(expr_out[start:end],
                          (row_ids[mask], np.searchsorted(target_col_arr, chunk_idx[mask])),
                          chunk_data[mask])
            if (start // CHUNK) % 5 == 0:
                print(f"    {end:,}/{n_cells:,}")

        if total_counts is not None:
            safe = np.where(total_counts > 0, total_counts, 1.0)
            expr_out = np.log1p(expr_out / safe[:, None] * 10_000)
        else:
            rs = expr_out.sum(axis=1, keepdims=True)
            expr_out = np.log1p(expr_out / np.where(rs > 0, rs, 1.0) * 10_000)

        return expr_out, found_names, cell_types_in, tissues_in
#score pe rcell
def compute_scores(expr, cell_types, tissues, gene_idx):
    groups = {}
    for i, (ct, ti) in enumerate(zip(cell_types, tissues)):
        groups.setdefault((ct, ti), []).append(i)

    rows = []
    for (ct, ti), idxs in groups.items():
        idxs = np.array(idxs)
        row  = {"cell_type": ct, "tissue": ti, "n_cells": len(idxs)}

        for mod in FUNCTIONAL_KEYS + ["Shedding", "Activation"]:
            cols = [gene_idx[g] for g in MODULES[mod]["genes"] if g in gene_idx]
            row[mod] = float(expr[np.ix_(idxs, cols)].mean()) if cols else 0.0

        anti = [gene_idx[g] for g in MODULES["Viability"]["anti"] if g in gene_idx]
        pro  = [gene_idx[g] for g in MODULES["Viability"]["pro"]  if g in gene_idx]
        row["Viability_raw"] = (
            float(expr[np.ix_(idxs, anti)].mean()) if anti else 0.0
        ) - (
            float(expr[np.ix_(idxs, pro)].mean())  if pro  else 0.0
        )
        rows.append(row)

    return pd.DataFrame(rows)


def _minmax(s):
    mn, mx = s.min(), s.max()
    return pd.Series(0.0, index=s.index) if mx == mn else (s - mn) / (mx - mn)


def compute_uber(df):
    df = df.copy()
    for col in FUNCTIONAL_KEYS + ["Shedding", "Activation", "Viability_raw"]:
        df[col + "_norm"] = _minmax(df[col])
    df["Viability"] = df["Viability_raw_norm"]
    df["Functional"] = sum(MODULES[k]["weight"] * df[k + "_norm"] for k in FUNCTIONAL_KEYS)
    df["Uber Score"] = (
        0.6 * (df["Functional"] + 0.10 * df["Shedding_norm"] - 0.25 * df["Activation_norm"])
        + 0.4 * df["Viability"]
    )
    df = df.sort_values("Uber Score", ascending=False).reset_index(drop=True)
    df.insert(0, "Rank", range(1, len(df) + 1))
    return df

def main():
    print("=" * 72)
    print("Mito Secondary Transfer Score Pipeline")
    print("=" * 72)

    print("\n[1] Loading expression data …")
    expr, cell_types, tissues, gene_idx = load_expression()
    print(f"  {expr.shape[0]:,} cells × {expr.shape[1]} genes | "
          f"{len(set(cell_types))} cell types | {len(set(tissues))} tissues")

    missing_g = [g for g in ALL_GENES if g not in gene_idx]
    print(f"  Genes matched: {len(gene_idx)}/{len(ALL_GENES)}"
          + (f"  MISSING: {missing_g}" if missing_g else ""))

    print("\n[2] Computing per-(cell_type × tissue) feature scores …")
    df_raw = compute_scores(expr, cell_types, tissues, gene_idx)
    print(f"  Groups: {len(df_raw)}")

    print("\n[3] Normalising (0→1) and computing Uber Score …")
    df = compute_uber(df_raw)

    norm_cols = [k + "_norm" for k in FUNCTIONAL_KEYS] + \
                ["Shedding_norm", "Activation_norm", "Viability"]
    save_cols = (["Rank","cell_type","tissue","n_cells","Uber Score","Functional","Viability"]
                 + FUNCTIONAL_KEYS + ["Shedding","Activation"]
                 + norm_cols)
    csv_path = os.path.join(OUT, "mito_sec_transfer_scores.csv")
    df[save_cols].to_csv(csv_path, index=False, encoding="utf-8")
    print(f"\n[4] Saved: {csv_path}  ({len(df)} rows)")

    print("\nGLOBAL TOP 20")
    print(f"  {'#':<4} {'Cell Type':<28} {'Tissue':<22} {'Uber':>6}  {'Funct':>6}  {'Viab':>6}")
    print("  " + "-" * 75)
    for _, r in df.head(20).iterrows():
        print(f"  {int(r['Rank']):<4} {r['cell_type']:<28} {r['tissue']:<22} "
              f"{r['Uber Score']:>6.3f}  {r['Functional']:>6.3f}  {r['Viability']:>6.3f}")

    print("\nPER-TISSUE TOP 3:")
    for tissue in sorted(df["tissue"].unique()):
        sub = df[df["tissue"] == tissue].head(3)
        best = sub.iloc[0]
        others = ", ".join(sub.iloc[1:]["cell_type"].tolist())
        print(f"  {tissue:<35} #1 {best['cell_type']:<25} {best['Uber Score']:.3f}")

    print("\nDone. Run visualize/visualize_sec_transfer.py for figures.")


if __name__ == "__main__":
    main()
