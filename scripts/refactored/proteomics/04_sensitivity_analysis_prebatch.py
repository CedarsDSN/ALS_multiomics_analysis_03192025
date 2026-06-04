# ============================================================
# Script: sensitivity_analysis_prebatch.py
# Project: ALS multi-omics analysis
# Category: Proteomics / batch correction sensitivity
#
# Purpose:
#   Runs a simple protein-level C9orf72-vs-control sensitivity analysis on pre-batch-corrected and post-batch-corrected proteomics matrices, using independent t-tests per protein, FDR correction, and overlap summaries.
#
# Inputs:
#   - Pre-batch proteomics matrix: `AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt`.
#   - Post-batch proteomics matrix: `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
#   - `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
#
# Outputs:
#   - `DE_results_C9orf72_vs_Control_PREBATCH.tsv`.
#   - `DE_results_C9orf72_vs_Control_POSTBATCH.tsv`.
#   - Console summary of sample sizes, significant protein counts, top proteins, and proteins gained/lost after correction.
#
# Upstream dependencies:
#   - No upstream script required beyond proteomics matrices and mutation metadata.
#
# How to run:
#   - Current: edit `BASE_DIR`, then run `python sensitivity_analysis_prebatch.py`.
#   - Refactored target included in this package: `python scripts_refactored/sensitivity_analysis_prebatch_argparse.py --prebatch <file> --postbatch <file> --subtypes <csv> --outdir analysis/batch_05_proteomics/results/prebatch_sensitivity`.
# ============================================================

#!/usr/bin/env python3
"""
Sensitivity analysis: C9orf72 vs Control proteomics before and after batch correction.

This is the command-line refactor of the uploaded script. It keeps the same scientific
logic: identify C9orf72 and healthy-control samples from mutation metadata, run an
independent t-test per protein, apply FDR correction, write pre- and post-batch results,
and summarize significant-protein overlap.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
try:
    from scipy.stats import false_discovery_control
except Exception:  # pragma: no cover
    false_discovery_control = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="C9orf72 vs Control proteomics pre/post batch sensitivity analysis")
    parser.add_argument("--prebatch", required=True, type=Path, help="Pre-batch-corrected protein matrix TSV")
    parser.add_argument("--postbatch", required=True, type=Path, help="Post-batch-corrected protein matrix TSV")
    parser.add_argument("--subtypes", required=True, type=Path, help="Mutation subtype metadata CSV")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--min-n", type=int, default=3, help="Minimum nonmissing, positive observations per group")
    return parser.parse_args()


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    if false_discovery_control is not None:
        return false_discovery_control(pvalues)
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def participant_from_sample(sample_name: str) -> str | None:
    match = re.match(r"(CASE|CTRL)[_-]([A-Z0-9]+)", sample_name)
    if not match:
        return None
    prefix, code = match.groups()
    return f"{prefix}-{code}"


def run_analysis(protein_file: Path, subtype_file: Path, label: str, min_n: int) -> tuple[pd.DataFrame, int, int]:
    print("=" * 80)
    print(f"Analysis: {label}")
    print("=" * 80)
    protein_df = pd.read_csv(protein_file, sep="\t")
    mutation_df = pd.read_csv(subtype_file)

    required_cols = {"Participant_ID", "Mutation_Subtype", "SUBJECT_GROUP"}
    missing = required_cols.difference(mutation_df.columns)
    if missing:
        raise ValueError(f"Subtype metadata missing columns: {sorted(missing)}")
    if "Protein" not in protein_df.columns:
        raise ValueError("Protein matrix must contain a 'Protein' column")

    sample_cols = [c for c in protein_df.columns if c not in {"Protein", "nFragment", "nPeptide"}]
    sample_info = pd.DataFrame({
        "Sample_ID": sample_cols,
        "Participant_ID": [participant_from_sample(c) for c in sample_cols],
    })
    sample_info = sample_info.merge(
        mutation_df[["Participant_ID", "Mutation_Subtype", "SUBJECT_GROUP"]],
        on="Participant_ID",
        how="left",
    )
    c9_samples = sample_info.loc[sample_info["Mutation_Subtype"] == "C9orf72", "Sample_ID"].tolist()
    ctrl_samples = sample_info.loc[sample_info["SUBJECT_GROUP"] == "Healthy Control", "Sample_ID"].tolist()
    c9_samples = [c for c in c9_samples if c in protein_df.columns]
    ctrl_samples = [c for c in ctrl_samples if c in protein_df.columns]
    print(f"C9orf72 samples: {len(c9_samples)}")
    print(f"Control samples: {len(ctrl_samples)}")

    results: list[dict[str, float | str | int]] = []
    for _, row in protein_df.iterrows():
        protein = row["Protein"]
        c9_vals = pd.to_numeric(row[c9_samples], errors="coerce").to_numpy(dtype=float)
        ctrl_vals = pd.to_numeric(row[ctrl_samples], errors="coerce").to_numpy(dtype=float)
        c9_clean = c9_vals[np.isfinite(c9_vals) & (c9_vals > 0)]
        ctrl_clean = ctrl_vals[np.isfinite(ctrl_vals) & (ctrl_vals > 0)]
        if len(c9_clean) < min_n or len(ctrl_clean) < min_n:
            continue
        c9_mean = float(np.mean(c9_clean))
        ctrl_mean = float(np.mean(ctrl_clean))
        log2fc = float(np.log2(c9_mean + 1.0) - np.log2(ctrl_mean + 1.0))
        t_stat, p_val = stats.ttest_ind(c9_clean, ctrl_clean, equal_var=False, nan_policy="omit")
        if not np.isfinite(p_val):
            continue
        results.append({
            "Protein": protein,
            "C9_mean": c9_mean,
            "Control_mean": ctrl_mean,
            "log2FC": log2fc,
            "pvalue": float(p_val),
            "t_statistic": float(t_stat),
            "n_C9": int(len(c9_clean)),
            "n_Control": int(len(ctrl_clean)),
        })
    results_df = pd.DataFrame(results)
    if results_df.empty:
        raise RuntimeError(f"No proteins passed filtering for {label}")
    results_df["FDR"] = bh_fdr(results_df["pvalue"].to_numpy(dtype=float))
    results_df = results_df.sort_values("FDR")
    return results_df, len(c9_samples), len(ctrl_samples)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for f in [args.prebatch, args.postbatch, args.subtypes]:
        if not f.exists():
            raise FileNotFoundError(f)

    results_pre, n_c9_pre, n_ctrl_pre = run_analysis(args.prebatch, args.subtypes, "PRE-batch corrected", args.min_n)
    results_post, n_c9_post, n_ctrl_post = run_analysis(args.postbatch, args.subtypes, "POST-batch corrected", args.min_n)

    out_pre = args.outdir / "DE_results_C9orf72_vs_Control_PREBATCH.tsv"
    out_post = args.outdir / "DE_results_C9orf72_vs_Control_POSTBATCH.tsv"
    results_pre.to_csv(out_pre, sep="\t", index=False)
    results_post.to_csv(out_post, sep="\t", index=False)

    sig_pre_05 = set(results_pre.loc[results_pre["FDR"] < 0.05, "Protein"])
    sig_post_05 = set(results_post.loc[results_post["FDR"] < 0.05, "Protein"])
    summary = pd.DataFrame([
        {"metric": "n_C9_pre", "value": n_c9_pre},
        {"metric": "n_Control_pre", "value": n_ctrl_pre},
        {"metric": "n_C9_post", "value": n_c9_post},
        {"metric": "n_Control_post", "value": n_ctrl_post},
        {"metric": "proteins_tested_pre", "value": len(results_pre)},
        {"metric": "proteins_tested_post", "value": len(results_post)},
        {"metric": "sig_FDR05_pre", "value": len(sig_pre_05)},
        {"metric": "sig_FDR05_post", "value": len(sig_post_05)},
        {"metric": "sig_FDR05_both", "value": len(sig_pre_05 & sig_post_05)},
        {"metric": "sig_FDR05_pre_only", "value": len(sig_pre_05 - sig_post_05)},
        {"metric": "sig_FDR05_post_only", "value": len(sig_post_05 - sig_pre_05)},
    ])
    summary_file = args.outdir / "pre_vs_post_batch_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)
    print(f"Wrote: {out_pre}")
    print(f"Wrote: {out_post}")
    print(f"Wrote: {summary_file}")


if __name__ == "__main__":
    main()
