
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import pandas as pd
import numpy as np

ROOT    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT     = os.path.join(ROOT, "output")
CSV_IN  = os.path.join(OUT, "mito_sec_transfer_scores.csv")

VESICLE_COLS = ["Exosome_norm", "Rab Traffic_norm", "MDV_norm", "MV Form._norm", "Shedding_norm"]
TNT_COL      = "TNT_norm"

def minmax(s):
    mn, mx = s.min(), s.max()
    return pd.Series(0.0, index=s.index) if mx == mn else (s - mn) / (mx - mn)


def load(path):
    df = pd.read_csv(path, encoding="utf-8")
    functional = ["Exosome","Mito Health","Rab Traffic","MDV","TNT",
                  "Adhesion","Metabolism","MV Form."]
    for col in functional + ["Shedding", "Activation"]:
        if col + "_norm" not in df.columns:
            df[col + "_norm"] = minmax(df[col])
    return df


def score_modes(df):
    df = df.copy()
    available_v = [c for c in VESICLE_COLS if c in df.columns]
    df["vesicle_score"] = df[available_v].mean(axis=1)
    df["tnt_score"]     = df[TNT_COL] if TNT_COL in df.columns else 0.0
    df["margin"]        = df["vesicle_score"] - df["tnt_score"]
    df["winner"]        = df["margin"].apply(lambda m: "Vesicle" if m > 0 else "TNT")
    df["dominance"]     = df["margin"].abs().apply(
        lambda m: "strong" if m > 0.15 else ("moderate" if m > 0.06 else "marginal")
    )
    return df


def build_report(top_per_tissue):
    lines = []
    lines.append("=" * 78)
    lines.append("VESICLE vs TNT TRANSFER — Top Cell Type per Tissue")
    lines.append("=" * 78)
    lines.append("")
    lines.append(
        "Vesicle score = mean(Exosome, Rab Traffic, MDV, MV Form., Shedding)  [0-1 norm]"
    )
    lines.append("TNT score     = TNT_norm  [0-1 norm]")
    lines.append("")
    lines.append(f"  {'Tissue':<32} {'Cell Type':<28} {'Vesl':>5}  {'TNT':>5}  {'Winner':<8}  {'Margin':>7}  Strength")
    lines.append("  " + "-" * 98)

    tally = {"Vesicle": 0, "TNT": 0}
    for _, r in top_per_tissue.sort_values("tissue").iterrows():
        tally[r["winner"]] += 1
        ct = r["cell_type"]
        if len(ct) > 26: ct = ct[:25] + "…"
        ti = r["tissue"]
        if len(ti) > 30: ti = ti[:29] + "…"
        lines.append(
            f"  {ti:<32} {ct:<28} {r['vesicle_score']:>5.3f}  {r['tnt_score']:>5.3f}"
            f"  {r['winner']:<8}  {r['margin']:>+7.3f}  {r['dominance']}"
        )

    lines.append("")
    lines.append("─" * 78)
    lines.append(f"SUMMARY:  Vesicle preferred in {tally['Vesicle']} tissues | "
                 f"TNT preferred in {tally['TNT']} tissues")

    strong   = top_per_tissue[top_per_tissue["dominance"] == "strong"]
    moderate = top_per_tissue[top_per_tissue["dominance"] == "moderate"]
    marginal = top_per_tissue[top_per_tissue["dominance"] == "marginal"]
    lines.append(f"           Strong ({len(strong)}) | Moderate ({len(moderate)}) | Marginal ({len(marginal)})")
    lines.append("")

    lines.append("TOP 5 TISSUES FOR VESICLE TRANSFER:")
    top_v = top_per_tissue[top_per_tissue["winner"] == "Vesicle"].nlargest(5, "vesicle_score")
    for _, r in top_v.iterrows():
        lines.append(f"  {r['tissue']:<35}  {r['cell_type']:<28}  vesicle={r['vesicle_score']:.3f}")

    lines.append("")
    lines.append("TOP 5 TISSUES FOR TNT TRANSFER:")
    top_t = top_per_tissue[top_per_tissue["winner"] == "TNT"].nlargest(5, "tnt_score")
    for _, r in top_t.iterrows():
        lines.append(f"  {r['tissue']:<35}  {r['cell_type']:<28}  TNT={r['tnt_score']:.3f}")

    return "\n".join(lines)


def main():
    if not os.path.exists(CSV_IN):
        print(f"ERROR: {CSV_IN} not found. Run mito_sec_transfer_score.py first.")
        sys.exit(1)

    df = score_modes(load(CSV_IN))

    top_idx = df.groupby("tissue")["Uber Score"].idxmax()
    top = df.loc[top_idx].copy().reset_index(drop=True)

    csv_out = os.path.join(OUT, "vescicle_vs_tnt.csv")
    top[["tissue","cell_type","Uber Score","vesicle_score","tnt_score",
         "winner","margin","dominance"]].to_csv(csv_out, index=False, encoding="utf-8")

    txt_out = os.path.join(OUT, "vescicle_vs_tnt.txt")
    report  = build_report(top)
    with open(txt_out, "w", encoding="utf-8") as fh:
        fh.write(report)

    print(report)
    print(f"\nSaved:\n  {csv_out}\n  {txt_out}")


if __name__ == "__main__":
    main()
