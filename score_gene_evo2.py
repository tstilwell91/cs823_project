#!/usr/bin/env python3
"""
Score FOXP3 / regulatory SNVs with Evo2 from the command line.

Default example:
    chrX:49269631:T>C

This script adapts the BRCA1 zero-shot VEP notebook pattern:
    reference sequence -> alternate sequence -> Evo2 likelihood score -> delta score

It can run either one variant from command-line arguments or the top N rows from an
AlphaGenome ranked CSV, such as FOXP3_alphagenome_variant_summary_bio_ranked.csv.

Coordinate convention:
    - Input positions are 1-based GRCh38/hg38 coordinates.
    - UCSC sequence API uses 0-based half-open intervals internally.

Outputs:
    - <outdir>/input_variants_for_evo2.csv
    - <outdir>/sequence_metadata.csv
    - <outdir>/<gene>_evo2_scores.csv
    - <outdir>/<gene>_alphagenome_vs_evo2.png, if matplotlib is available
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd

VALID_BASES = set("ACGT")
UCSC_SEQUENCE_URL = "https://api.genome.ucsc.edu/getData/sequence"


def parse_variant_string(variant: str) -> Tuple[str, int, str, str]:
    """Parse chrX:49269631:T>C or chrX:49269631:T:C."""
    try:
        chrom, pos, change = variant.strip().split(":", 2)
        if ">" in change:
            ref, alt = change.split(">", 1)
        elif ":" in change:
            ref, alt = change.split(":", 1)
        else:
            raise ValueError
        return chrom, int(pos), ref.upper(), alt.upper()
    except Exception as exc:
        raise argparse.ArgumentTypeError(
            "Variant must look like chrX:49269631:T>C or chrX:49269631:T:C"
        ) from exc


def fetch_reference_sequence_ucsc(chrom: str, start0: int, end0: int, genome: str = "hg38") -> str:
    """Fetch reference DNA sequence from UCSC using a 0-based, half-open interval."""
    params = urlencode({
        "genome": genome,
        "chrom": chrom,
        "start": int(start0),
        "end": int(end0),
    })
    url = f"{UCSC_SEQUENCE_URL}?{params}"
    with urlopen(url, timeout=90) as response:
        payload = response.read().decode("utf-8")
    data = json.loads(payload)
    if "dna" not in data:
        raise RuntimeError(f"UCSC response did not contain DNA. URL={url}. Response={data}")
    seq = data["dna"].upper()
    expected_len = int(end0) - int(start0)
    if len(seq) != expected_len:
        raise RuntimeError(f"Expected {expected_len} bp from UCSC, got {len(seq)} bp. URL={url}")
    return seq


def parse_ref_alt_sequences(
    chrom: str,
    pos_1based: int,
    ref: str,
    alt: str,
    window_size: int,
    genome: str,
) -> Dict[str, object]:
    """Return REF and ALT sequence windows centered on a 1-based SNV coordinate."""
    ref = str(ref).upper()
    alt = str(alt).upper()

    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("Only SNVs are supported: one reference base and one alternate base.")
    if ref not in VALID_BASES or alt not in VALID_BASES:
        raise ValueError(f"REF and ALT must be A/C/G/T. Got ref={ref}, alt={alt}")
    if ref == alt:
        raise ValueError(f"REF and ALT are identical for {chrom}:{pos_1based}:{ref}>{alt}")

    p0 = int(pos_1based) - 1
    start0 = max(0, p0 - window_size // 2)
    end0 = start0 + window_size
    snv_index = p0 - start0

    ref_seq = fetch_reference_sequence_ucsc(chrom, start0, end0, genome=genome)
    observed_ref = ref_seq[snv_index]
    if observed_ref != ref:
        raise AssertionError(
            f"Reference base mismatch at {chrom}:{pos_1based}. "
            f"Expected {ref}, observed {observed_ref}. "
            f"Check genome build, coordinate system, and REF/ALT orientation."
        )

    alt_seq = ref_seq[:snv_index] + alt + ref_seq[snv_index + 1:]
    if len(ref_seq) != window_size or len(alt_seq) != window_size:
        raise AssertionError("Internal sequence length error.")

    return {
        "chromosome": chrom,
        "position": int(pos_1based),
        "start0": int(start0),
        "end0": int(end0),
        "snv_index": int(snv_index),
        "observed_ref": observed_ref,
        "ref_seq": ref_seq,
        "alt_seq": alt_seq,
    }


def make_single_variant_df(args: argparse.Namespace) -> pd.DataFrame:
    chrom, pos, ref, alt = args.variant
    return pd.DataFrame([
        {
            "gene_symbol": args.gene_symbol,
            "chromosome": chrom,
            "position": pos,
            "reference_base": ref,
            "alternate_base": alt,
            "variant_id_custom": f"{chrom}:{pos}:{ref}>{alt}",
            "alphagenome_abs_score": args.alphagenome_abs_score,
            "alphagenome_output_type": args.alphagenome_output_type,
            "alphagenome_track_name": args.alphagenome_track_name,
            "alphagenome_biosample_name": args.alphagenome_biosample_name,
        }
    ])


def load_variants_from_csv(path: str, top_n: int) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = ["gene_symbol", "chromosome", "position", "reference_base", "alternate_base", "variant_id_custom"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    rename_map = {
        "abs_score": "alphagenome_abs_score",
        "output_type": "alphagenome_output_type",
        "track_name": "alphagenome_track_name",
        "biosample_name": "alphagenome_biosample_name",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    df = df.drop_duplicates(subset=["variant_id_custom"]).head(top_n).copy()
    return df.reset_index(drop=True)


def import_evo2_class():
    try:
        from evo2.models import Evo2  # official BRCA1 notebook path
        return Evo2
    except Exception:
        from evo2 import Evo2  # fallback used in some installs
        return Evo2


def score_in_batches(model, sequences: List[str], batch_size: int) -> List[float]:
    scores: List[float] = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        print(f"Scoring batch {i // batch_size + 1}: {i + 1}-{i + len(batch)} of {len(sequences)}", flush=True)
        batch_scores = model.score_sequences(batch)
        scores.extend([float(x) for x in batch_scores])
    return scores


def save_plot(results_df: pd.DataFrame, out_png: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Skipping plot because matplotlib could not be imported: {exc}")
        return

    plt.figure(figsize=(6, 4))
    if "alphagenome_abs_score" in results_df.columns and results_df["alphagenome_abs_score"].notna().any():
        x = results_df["alphagenome_abs_score"].astype(float)
        y = results_df["evo2_delta_score"].astype(float)
        plt.scatter(x, y, s=60)
        for _, row in results_df.iterrows():
            plt.text(float(row["alphagenome_abs_score"]), float(row["evo2_delta_score"]), str(row["variant_id_custom"]), fontsize=7)
        plt.xlabel("AlphaGenome absolute regulatory score")
    else:
        x = np.arange(len(results_df))
        y = results_df["evo2_delta_score"].astype(float)
        plt.scatter(x, y, s=60)
        plt.xticks(x, results_df["variant_id_custom"], rotation=45, ha="right", fontsize=8)
        plt.xlabel("Variant")

    plt.axhline(0, linestyle="--", linewidth=1)
    plt.ylabel("Evo2 delta score (ALT - REF)")
    plt.title("AlphaGenome-prioritized variants scored by Evo2")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved plot: {out_png}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Score SNVs with Evo2 using REF vs ALT sequence likelihood deltas."
    )
    parser.add_argument(
        "--variant",
        type=parse_variant_string,
        default=parse_variant_string("chrX:49269631:T>C"),
        help="Single variant, e.g. chrX:49269631:T>C. Ignored if --alphagenome-csv is used.",
    )
    parser.add_argument("--gene-symbol", default="FOXP3")
    parser.add_argument("--model-name", default="evo2_7b", help="Model name, e.g. evo2_7b or evo2_1b_base.")
    parser.add_argument("--genome", default="hg38", help="UCSC genome name. Use hg38 for GRCh38 coordinates.")
    parser.add_argument("--window-size", type=int, default=8192)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--output-dir", default="FOXP3_evo2_outputs")
    parser.add_argument("--alphagenome-csv", default=None, help="Optional AlphaGenome ranked CSV to score multiple variants.")
    parser.add_argument("--top-n", type=int, default=20, help="Top unique variants to score from --alphagenome-csv.")
    parser.add_argument("--skip-model", action="store_true", help="Build REF/ALT sequences only, do not load Evo2 or score.")
    parser.add_argument("--no-plot", action="store_true")

    # Metadata defaults for the top FOXP3 AlphaGenome hit.
    parser.add_argument("--alphagenome-abs-score", type=float, default=155104.0)
    parser.add_argument("--alphagenome-output-type", default="CHIP_HISTONE")
    parser.add_argument("--alphagenome-track-name", default="CL:0000792 Histone ChIP-seq H3K27ac")
    parser.add_argument("--alphagenome-biosample-name", default="CD4-positive, CD25-positive, alpha-beta regulatory T cell")

    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=== Evo2 SNV scoring ===")
    print(f"Output directory: {outdir.resolve()}")
    print(f"Model: {args.model_name}")
    print(f"Genome: {args.genome}")
    print(f"Window size: {args.window_size}")

    if args.alphagenome_csv:
        variants_df = load_variants_from_csv(args.alphagenome_csv, args.top_n)
        print(f"Loaded {len(variants_df)} unique variants from {args.alphagenome_csv}")
    else:
        variants_df = make_single_variant_df(args)
        print(f"Using single variant: {variants_df.loc[0, 'variant_id_custom']}")

    variants_csv = outdir / "input_variants_for_evo2.csv"
    variants_df.to_csv(variants_csv, index=False)
    print(f"Saved: {variants_csv}")

    records = []
    ref_sequences: List[str] = []
    alt_sequences: List[str] = []

    print("Fetching reference sequences and creating ALT sequences...")
    for idx, row in variants_df.iterrows():
        parsed = parse_ref_alt_sequences(
            chrom=str(row["chromosome"]),
            pos_1based=int(row["position"]),
            ref=str(row["reference_base"]),
            alt=str(row["alternate_base"]),
            window_size=args.window_size,
            genome=args.genome,
        )
        records.append(parsed)
        ref_sequences.append(str(parsed["ref_seq"]))
        alt_sequences.append(str(parsed["alt_seq"]))
        print(
            f"  {row['variant_id_custom']}: window {parsed['chromosome']}:{parsed['start0']}-{parsed['end0']}, "
            f"index={parsed['snv_index']}, REF={parsed['observed_ref']} ALT={row['alternate_base']}",
            flush=True,
        )

    sequence_metadata = pd.DataFrame([{k: v for k, v in rec.items() if k not in ("ref_seq", "alt_seq")} for rec in records])
    meta_csv = outdir / "sequence_metadata.csv"
    sequence_metadata.to_csv(meta_csv, index=False)
    print(f"Saved: {meta_csv}")

    if args.skip_model:
        print("--skip-model was provided. Stopping before Evo2 load/scoring.")
        return 0

    # Deduplicate identical REF windows. Variants near each other may share the same REF sequence.
    ref_seq_to_index: Dict[str, int] = {}
    unique_ref_sequences: List[str] = []
    ref_seq_indexes: List[int] = []

    for seq in ref_sequences:
        if seq not in ref_seq_to_index:
            ref_seq_to_index[seq] = len(unique_ref_sequences)
            unique_ref_sequences.append(seq)
        ref_seq_indexes.append(ref_seq_to_index[seq])

    ref_seq_indexes_arr = np.array(ref_seq_indexes)

    print(f"Unique REF sequences: {len(unique_ref_sequences)}")
    print(f"ALT sequences: {len(alt_sequences)}")

    Evo2 = import_evo2_class()
    print(f"Loading Evo2 model: {args.model_name}", flush=True)
    start_time = time.time()
    model = Evo2(args.model_name)
    print(f"Evo2 model loaded in {(time.time() - start_time) / 60:.2f} minutes.", flush=True)

    print("Scoring reference sequences...")
    ref_scores_unique = score_in_batches(model, unique_ref_sequences, batch_size=args.batch_size)
    print("Scoring alternate sequences...")
    alt_scores = score_in_batches(model, alt_sequences, batch_size=args.batch_size)

    ref_scores = np.array(ref_scores_unique, dtype=float)[ref_seq_indexes_arr]
    alt_scores_arr = np.array(alt_scores, dtype=float)
    delta_scores = alt_scores_arr - ref_scores

    results_df = variants_df.copy()
    results_df["window_size"] = args.window_size
    results_df["evo2_model_name"] = args.model_name
    results_df["evo2_ref_score"] = ref_scores
    results_df["evo2_alt_score"] = alt_scores_arr
    results_df["evo2_delta_score"] = delta_scores
    results_df["evo2_abs_delta_score"] = np.abs(delta_scores)
    results_df["evo2_interpretation"] = np.where(
        results_df["evo2_delta_score"] < 0,
        "ALT lower likelihood than REF",
        "ALT higher likelihood than REF",
    )

    gene_for_name = str(results_df["gene_symbol"].iloc[0]) if "gene_symbol" in results_df.columns else args.gene_symbol
    out_csv = outdir / f"{gene_for_name}_evo2_scores.csv"
    results_df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")
    print(results_df[["variant_id_custom", "evo2_ref_score", "evo2_alt_score", "evo2_delta_score", "evo2_interpretation"]].to_string(index=False))

    if not args.no_plot:
        save_plot(results_df, outdir / f"{gene_for_name}_alphagenome_vs_evo2.png")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
