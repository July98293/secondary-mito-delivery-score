

import os, sys, io, warnings
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm, ListedColormap

ROOT    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT     = os.path.join(ROOT, "otput")
CSV     = os.path.join(OUT, "mito_sec_transfer_scores.csv")

FUNCTIONAL_KEYS = ["Exosome","Mito Health","Rab Traffic","MDV","TNT",
                   "Adhesion","Metabolism","MV Form."]
FEATURE_KEYS    = FUNCTIONAL_KEYS + ["Shedding","Activation","Viability"]
NORM_COLS       = [k + "_norm" for k in FUNCTIONAL_KEYS] + \
                  ["Shedding_norm","Activation_norm","Viability"]
FEATURE_LABELS  = FUNCTIONAL_KEYS + ["Shedding","Activation ↓","Viability"]

KEY_TISSUES = [
    "blood", "bone marrow", "lymph node", "lung",
    "spleen", "thymus", "liver", "adipose tissue",
    "large intestine", "kidney", "skin of body", "heart",
]


def _minmax(s):
    mn, mx = s.min(), s.max()
    return pd.Series(0.0, index=s.index) if mx == mn else (s - mn) / (mx - mn)


def load_data():
    df = pd.read_csv(CSV, encoding="utf-8")
    # Re-normalise raw feature scores globally (same as pipeline)
    for col in FUNCTIONAL_KEYS + ["Shedding", "Activation"]:
        if col + "_norm" not in df.columns:
            df[col + "_norm"] = _minmax(df[col])
    # Viability is already normalised in CSV
    if "Viability" not in df.columns:
        df["Viability"] = 0.0
    return df


def abbrev_ct(name, maxlen=18):
    replacements = {
        "classical monocyte": "cl. monocyte",
        "intermediate monocyte": "int. monocyte",
        "natural killer cell": "NK cell",
        "hematopoietic precursor cell": "HSC",
        "innate lymphoid cell": "ILC",
        "regulatory T cell": "Treg",
        "tissue-resident macrophage": "tissue macrophage",
        "CD4-positive, alpha-beta T cell": "CD4+ T",
        "CD8-positive, alpha-beta T cell": "CD8+ T",
        "mature NK T cell": "NKT cell",
        "mononuclear phagocyte": "mono. phago.",
        "hematopoietic cell": "hema. cell",
    }
    name = replacements.get(name, name)
    return name[:maxlen] + "…" if len(name) > maxlen else name



def viz_best_per_tissue(df, out_path):
    best = (
        df.loc[df.groupby("tissue")["Uber Score"].idxmax()]
          .copy()
          .sort_values("Uber Score", ascending=True)
    )

    # assign a colour to each cell type
    cell_types = best["cell_type"].unique()
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(cell_types), 1)))
    ct_color = {ct: palette[i % len(palette)] for i, ct in enumerate(cell_types)}

    fig_h = max(10, len(best) * 0.28)
    fig, ax = plt.subplots(figsize=(12, fig_h))

    bars = ax.barh(
        best["tissue"],
        best["Uber Score"],
        color=[ct_color[ct] for ct in best["cell_type"]],
        edgecolor="white", linewidth=0.4, height=0.7,
    )

    for bar, (_, row) in zip(bars, best.iterrows()):
        w = bar.get_width()
        ax.text(w + 0.004, bar.get_y() + bar.get_height() / 2,
                f"{abbrev_ct(row['cell_type'])}  {w:.3f}",
                va="center", ha="left", fontsize=7.5)

    handles = [mpatches.Patch(color=ct_color[ct], label=abbrev_ct(ct))
               for ct in sorted(ct_color)]
    ax.legend(handles=handles, title="Cell Type", loc="lower right",
              fontsize=7, title_fontsize=8, ncol=2,
              framealpha=0.85, edgecolor="gray")

    ax.set_xlim(0, best["Uber Score"].max() * 1.35)
    ax.set_xlabel("Über Score", fontsize=11)
    ax.set_title("Best Cell Type per Tissue  —  Mito Secondary Transfer Score",
                 fontsize=13, fontweight="bold", pad=10)
    ax.tick_params(axis="y", labelsize=7.5)
    ax.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def viz_celltype_x_tissue(df, out_path):
    pivot = df.pivot_table(
        index="cell_type", columns="tissue",
        values="Uber Score", aggfunc="max"
    )

    ct_cover  = pivot.notna().sum(axis=1)
    tis_cover = pivot.notna().sum(axis=0)
    pivot = pivot.loc[ct_cover >= 8, tis_cover >= 5]

    row_means = pivot.mean(axis=1, skipna=True)
    col_means = pivot.mean(axis=0, skipna=True)
    pivot = pivot.loc[row_means.sort_values(ascending=False).index,
                      col_means.sort_values(ascending=False).head(30).index]

    mat = pivot.values.astype(float)
    row_labels = [abbrev_ct(ct, 22) for ct in pivot.index]
    col_labels  = list(pivot.columns)
    n_rows, n_cols = mat.shape

    fig_w = max(14, n_cols * 0.45)
    fig_h = max(6,  n_rows * 0.40)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd",
                   vmin=np.nanmin(mat), vmax=np.nanmax(mat),
                   interpolation="nearest")

    for i in range(n_rows):
        for j in range(n_cols):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=5.5, color="black" if v < 0.55 else "white",
                        fontweight="bold")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title("Über Score by Cell Type × Tissue  —  Mito Secondary Transfer",
                 fontsize=12, fontweight="bold", pad=10)
    ax.set_xlabel("Tissue", fontsize=10)
    ax.set_ylabel("Cell Type", fontsize=10)

    plt.colorbar(im, ax=ax, fraction=0.015, pad=0.01, label="Über Score")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def viz_feature_x_celltype_by_organ(df, out_path):
    # Keep only key tissues that exist in data
    available = set(df["tissue"].str.lower().unique())
    tissues_to_plot = [t for t in KEY_TISSUES if t in available]

    n_tissues = len(tissues_to_plot)
    n_cols    = 3
    n_rows    = int(np.ceil(n_tissues / n_cols))
    fig_w     = n_cols * 7.5
    fig_h     = n_rows * 5.5

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h),
                             squeeze=False)

    for idx, tissue in enumerate(tissues_to_plot):
        ax  = axes[idx // n_cols][idx % n_cols]
        sub = (df[df["tissue"].str.lower() == tissue]
               .sort_values("Uber Score", ascending=False)
               .head(8))

        if sub.empty:
            ax.set_visible(False); continue

        ct_labels = [abbrev_ct(ct, 16) for ct in sub["cell_type"]]

        
        mat = np.zeros((len(FEATURE_KEYS), len(sub)))
        for fi, feat in enumerate(FEATURE_KEYS):
            if feat == "Viability":
                col = "Viability"
            elif feat == "Activation":
                col = "Activation_norm"
            elif feat == "Shedding":
                col = "Shedding_norm"
            else:
                col = feat + "_norm"

            if col in sub.columns:
                vals = sub[col].values.astype(float)
            else:
                vals = np.zeros(len(sub))
            mn, mx = vals.min(), vals.max()
            mat[fi] = (vals - mn) / (mx - mn) if mx > mn else vals

        act_fi = FEATURE_KEYS.index("Activation")

        im = ax.imshow(mat, aspect="auto", cmap="RdYlGn",
                       vmin=0, vmax=1, interpolation="nearest")

        act_row = (1 - mat[act_fi]).reshape(1, -1)
        ax.imshow(act_row, aspect="auto", cmap="RdYlGn",
                  vmin=0, vmax=1, interpolation="nearest",
                  extent=(-0.5, len(sub) - 0.5, act_fi + 0.5, act_fi - 0.5))

        for fi in range(len(FEATURE_KEYS)):
            for ci in range(len(sub)):
                v = mat[fi, ci]
                ax.text(ci, fi, f"{v:.2f}", ha="center", va="center",
                        fontsize=6, color="black" if 0.25 < v < 0.75 else "white",
                        fontweight="bold")

        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(ct_labels, rotation=40, ha="right", fontsize=7)
        ax.set_yticks(range(len(FEATURE_LABELS)))
        ax.set_yticklabels(FEATURE_LABELS, fontsize=8)
        ax.set_title(tissue.title(), fontsize=10, fontweight="bold")

        uber_vals = sub["Uber Score"].values
        ax2 = ax.inset_axes([0, -0.22, 1, 0.14])
        ax2.bar(range(len(sub)), uber_vals, color="steelblue", width=0.7)
        ax2.set_xlim(-0.5, len(sub) - 0.5)
        ax2.set_ylim(0, max(uber_vals) * 1.2)
        ax2.set_xticks([]); ax2.set_yticks([0, round(max(uber_vals), 2)])
        ax2.tick_params(labelsize=6)
        ax2.set_ylabel("Über", fontsize=6)
        ax2.spines[["top","right"]].set_visible(False)

    for idx in range(n_tissues, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    fig.suptitle(
        "Feature Scores by Cell Type per Tissue  —  Mito Secondary Transfer\n"
        "(green = high, red = low;  Activation ↓ is inverted)",
        fontsize=13, fontweight="bold", y=1.01
    )
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# main
def main():
    print("=" * 60)
    print("Mito Secondary Transfer — Visualizations")
    print("=" * 60)

    if not os.path.exists(CSV):
        print(f"ERROR: {CSV} not found.")
        print("Run pipeline/mito_sec_transfer_score.py first.")
        sys.exit(1)

    print(f"\nLoading {CSV} …")
    df = load_data()
    print(f"  {len(df)} rows | {df['cell_type'].nunique()} cell types | "
          f"{df['tissue'].nunique()} tissues")

    print("\n[1] Best cell type per tissue …")
    viz_best_per_tissue(df, os.path.join(OUT, "sec_best_cell_per_tissue.png"))

    print("\n[2] Cell type × tissue heatmap (Uber Score) …")
    viz_celltype_x_tissue(df, os.path.join(OUT, "sec_heatmap_celltype_x_tissue.png"))

    print("\n[3] Feature × cell type by organ …")
    viz_feature_x_celltype_by_organ(
        df, os.path.join(OUT, "sec_heatmap_feature_x_celltype_by_organ.png")
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
