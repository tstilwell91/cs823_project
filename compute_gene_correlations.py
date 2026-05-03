#!/usr/bin/env python3
"""
Compute AlphaGenome vs Evo2 correlations for one gene.

Example:
  python compute_gene_correlations.py \
    --gene-symbol BACH2 \
    --alphagenome-csv BACH2_alphagenome_variant_summary_bio_ranked.csv \
    --evo2-csv BACH2_evo2_outputs_top100/BACH2_evo2_scores.csv \
    --output-dir BACH2_final_correlation

If your Evo2 output file is still named FOXP3_evo2_scores.csv inside a gene-specific
folder, that is okay. Just point --evo2-csv to that file.
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def require_columns(df, cols, file_label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{file_label} is missing required columns: {missing}")


def choose_alphagenome_score_column(df):
    if "abs_score" in df.columns:
        return "abs_score"
    if "alphagenome_abs_score" in df.columns:
        return "alphagenome_abs_score"
    raise ValueError(
        "Could not find AlphaGenome score column. Expected 'abs_score' "
        "or 'alphagenome_abs_score'."
    )


def safe_correlation(x, y, label):
    if len(x) < 3:
        print(f"WARNING: Only {len(x)} variants for {label}; correlation may not be meaningful.")
    if x.nunique() < 2 or y.nunique() < 2:
        print(f"WARNING: Constant input for {label}; returning NaN correlations.")
        return {
            "pearson_r": float("nan"),
            "pearson_p": float("nan"),
            "spearman_rho": float("nan"),
            "spearman_p": float("nan"),
        }

    pearson = pearsonr(x, y)
    spearman = spearmanr(x, y)

    return {
        "pearson_r": pearson.statistic,
        "pearson_p": pearson.pvalue,
        "spearman_rho": spearman.statistic,
        "spearman_p": spearman.pvalue,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compute correlations between AlphaGenome scores and Evo2 delta scores."
    )
    parser.add_argument("--gene-symbol", required=True, help="Gene symbol, e.g. FOXP3, BACH2, SATB1.")
    parser.add_argument("--alphagenome-csv", required=True, help="AlphaGenome ranked variant CSV.")
    parser.add_argument("--evo2-csv", required=True, help="Evo2 scores CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for final correlation outputs.")
    parser.add_argument(
        "--top-n-label",
        default="top100",
        help="Label to include in output filenames, e.g. top20, top100. Default: top100.",
    )
    parser.add_argument(
        "--top-table-n",
        type=int,
        default=10,
        help="Number of rows for top AlphaGenome and top Evo2 tables. Default: 10.",
    )
    parser.add_argument(
        "--print-columns",
        action="store_true",
        help="Print input CSV columns for debugging.",
    )

    args = parser.parse_args()

    gene = args.gene_symbol.upper()
    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(args.alphagenome_csv):
        raise FileNotFoundError(f"AlphaGenome CSV not found: {args.alphagenome_csv}")
    if not os.path.exists(args.evo2_csv):
        raise FileNotFoundError(f"Evo2 CSV not found: {args.evo2_csv}")

    ag = pd.read_csv(args.alphagenome_csv)
    evo = pd.read_csv(args.evo2_csv)

    if args.print_columns:
        print("AlphaGenome columns:")
        print(ag.columns.tolist())
        print("\nEvo2 columns:")
        print(evo.columns.tolist())

    require_columns(ag, ["variant_id_custom"], "AlphaGenome CSV")
    require_columns(evo, ["variant_id_custom", "evo2_ref_score", "evo2_alt_score", "evo2_delta_score"], "Evo2 CSV")

    ag_score_col = choose_alphagenome_score_column(ag)

    # Keep one AlphaGenome row per variant.
    # If AlphaGenome contains multiple tracks per variant, retain the highest absolute regulatory score.
    ag_one = (
        ag.sort_values(ag_score_col, ascending=False)
          .drop_duplicates("variant_id_custom", keep="first")
          .copy()
          .rename(columns={ag_score_col: "alphagenome_abs_score"})
    )

    # Keep one Evo2 row per variant.
    evo_one = (
        evo.sort_values("variant_id_custom")
           .drop_duplicates("variant_id_custom", keep="first")
           .copy()
    )

    if "evo2_abs_delta_score" not in evo_one.columns:
        evo_one["evo2_abs_delta_score"] = evo_one["evo2_delta_score"].abs()

    evo_cols = [
        "variant_id_custom",
        "evo2_ref_score",
        "evo2_alt_score",
        "evo2_delta_score",
        "evo2_abs_delta_score",
    ]

    # Preserve any useful Evo2/AlphaGenome columns if they exist, but avoid duplicates.
    merged = ag_one.merge(
        evo_one[evo_cols],
        on="variant_id_custom",
        how="inner",
    )

    if len(merged) == 0:
        raise ValueError(
            "No variants overlapped between AlphaGenome and Evo2 files. "
            "Check variant_id_custom formatting."
        )

    x = merged["alphagenome_abs_score"]
    y_signed = merged["evo2_delta_score"]
    y_abs = merged["evo2_abs_delta_score"]

    signed_corr = safe_correlation(x, y_signed, "signed Evo2 delta")
    abs_corr = safe_correlation(x, y_abs, "absolute Evo2 delta")

    summary = pd.DataFrame(
        [
            {
                "gene": gene,
                "comparison": "AlphaGenome abs score vs signed Evo2 delta",
                "n": len(merged),
                **signed_corr,
            },
            {
                "gene": gene,
                "comparison": "AlphaGenome abs score vs absolute Evo2 delta",
                "n": len(merged),
                **abs_corr,
            },
        ]
    )

    desc = pd.DataFrame(
        [
            {
                "gene": gene,
                "n_variants": len(merged),
                "alphagenome_min": x.min(),
                "alphagenome_median": x.median(),
                "alphagenome_max": x.max(),
                "evo2_delta_min": y_signed.min(),
                "evo2_delta_median": y_signed.median(),
                "evo2_delta_max": y_signed.max(),
                "evo2_abs_delta_median": y_abs.median(),
                "evo2_abs_delta_max": y_abs.max(),
                "n_alt_lower_likelihood": int((y_signed < 0).sum()),
                "n_alt_higher_likelihood": int((y_signed > 0).sum()),
                "n_alt_equal_likelihood": int((y_signed == 0).sum()),
            }
        ]
    )

    label = args.top_n_label

    merged_path = os.path.join(args.output_dir, f"{gene}_alphagenome_evo2_merged_{label}.csv")
    summary_path = os.path.join(args.output_dir, f"{gene}_correlation_summary_{label}.csv")
    desc_path = os.path.join(args.output_dir, f"{gene}_descriptive_summary_{label}.csv")
    top_ag_path = os.path.join(args.output_dir, f"{gene}_top{args.top_table_n}_alphagenome_with_evo2.csv")
    top_evo2_path = os.path.join(args.output_dir, f"{gene}_top{args.top_table_n}_evo2_disruptive.csv")

    merged.to_csv(merged_path, index=False)
    summary.to_csv(summary_path, index=False)
    desc.to_csv(desc_path, index=False)

    top_cols = ["variant_id_custom", "alphagenome_abs_score", "evo2_delta_score", "evo2_abs_delta_score"]
    extra_cols = [
        c for c in ["output_type", "track_name", "biosample_name", "region_id", "chrom", "pos", "ref", "alt"]
        if c in merged.columns and c not in top_cols
    ]
    top_cols = top_cols + extra_cols

    top_ag = merged.sort_values("alphagenome_abs_score", ascending=False).head(args.top_table_n)[top_cols]
    top_evo2 = merged.sort_values("evo2_abs_delta_score", ascending=False).head(args.top_table_n)[top_cols]

    top_ag.to_csv(top_ag_path, index=False)
    top_evo2.to_csv(top_evo2_path, index=False)

    # Scatter plot: signed Evo2 delta.
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y_signed)
    plt.axhline(0, linestyle="--")
    plt.xlabel("AlphaGenome absolute regulatory score")
    plt.ylabel("Evo2 delta score (ALT - REF)")
    plt.title(f"{gene}: AlphaGenome-prioritized variants scored by Evo2")

    text = (
        f"n = {len(merged)}\n"
        f"Pearson r = {signed_corr['pearson_r']:.3f}, p = {signed_corr['pearson_p']:.3g}\n"
        f"Spearman ρ = {signed_corr['spearman_rho']:.3f}, p = {signed_corr['spearman_p']:.3g}"
    )
    plt.text(
        0.02,
        0.98,
        text,
        transform=plt.gca().transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", alpha=0.2),
    )
    plt.tight_layout()
    signed_plot_path = os.path.join(args.output_dir, f"{gene}_alphagenome_vs_evo2_signed_delta_{label}.png")
    plt.savefig(signed_plot_path, dpi=200)
    plt.close()

    # Scatter plot: absolute Evo2 delta.
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y_abs)
    plt.xlabel("AlphaGenome absolute regulatory score")
    plt.ylabel("|Evo2 delta score|")
    plt.title(f"{gene}: AlphaGenome score vs absolute Evo2 disruption")

    text = (
        f"n = {len(merged)}\n"
        f"Pearson r = {abs_corr['pearson_r']:.3f}, p = {abs_corr['pearson_p']:.3g}\n"
        f"Spearman ρ = {abs_corr['spearman_rho']:.3f}, p = {abs_corr['spearman_p']:.3g}"
    )
    plt.text(
        0.02,
        0.98,
        text,
        transform=plt.gca().transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", alpha=0.2),
    )
    plt.tight_layout()
    abs_plot_path = os.path.join(args.output_dir, f"{gene}_alphagenome_vs_evo2_abs_delta_{label}.png")
    plt.savefig(abs_plot_path, dpi=200)
    plt.close()

    print("\n=== Correlation summary ===")
    print(summary.to_string(index=False))

    print("\n=== Descriptive summary ===")
    print(desc.to_string(index=False))

    print("\nSaved outputs:")
    for path in [
        merged_path,
        summary_path,
        desc_path,
        top_ag_path,
        top_evo2_path,
        signed_plot_path,
        abs_plot_path,
    ]:
        print(f"  {path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
