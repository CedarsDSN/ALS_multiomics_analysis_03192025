# ============================================================
# Script: sensitivity_analysis_prebatch.py
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from docs/SCRIPT_SPECS_FULL.md
# 
# ## `sensitivity_analysis_prebatch.py`
# 
# **Category:** Proteomics / batch correction sensitivity
# 
# **Purpose**
# 
# Runs a simple protein-level C9orf72-vs-control sensitivity analysis on pre-batch-corrected and post-batch-corrected proteomics matrices, using independent t-tests per protein, FDR correction, and overlap summaries.
# 
# **Inputs**
# 
# - Pre-batch proteomics matrix: `AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt`.
# - Post-batch proteomics matrix: `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
# - `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
# 
# **Outputs**
# 
# - `DE_results_C9orf72_vs_Control_PREBATCH.tsv`.
# - `DE_results_C9orf72_vs_Control_POSTBATCH.tsv`.
# - Console summary of sample sizes, significant protein counts, top proteins, and proteins gained/lost after correction.
# 
# **Upstream dependencies**
# 
# - No upstream script required beyond proteomics matrices and mutation metadata.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# python scripts/documented/proteomics/sensitivity_analysis_prebatch.py
# ```
# 
# Refactored command-line version:
# 
# ```bash
# python scripts/refactored/proteomics/04_sensitivity_analysis_prebatch.py \
#   --prebatch /path/to/prebatch.txt \
#   --postbatch /path/to/postbatch.txt \
#   --subtypes /path/to/ALS_Mutation_Subtypes_with_ALS_Control.csv \
#   --outdir results/proteomics/prebatch_sensitivity
# ```
# 
# **Refactor notes**
# 
# - Use `argparse` instead of hard-coded paths.
# - Write a comparison summary TSV in addition to console output.
# - Optionally replace t-tests with limma-style modeling if covariate adjustment is required.
# ============================================================

#!/usr/bin/env python
"""
Sensitivity Analysis: C9orf72 vs Control on PRE-batch corrected data
Comparing to post-batch corrected results
"""
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import sys

print("="*80)
print("SENSITIVITY ANALYSIS: Pre-Batch Corrected Data")
print("C9orf72 vs Control Differential Expression")
print("="*80)

# Paths
BASE_DIR = Path("/path/to/private/ALS")
PROTEOMICS_DIR = BASE_DIR / "proteomics" / "4_matrix"
OUTPUT_DIR = BASE_DIR / "proteomics" / "4_matrix" / "outputs" / "DE_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# PRE-batch corrected data
PROTEIN_MATRIX_PREBATCH = PROTEOMICS_DIR / "AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt"
PROTEIN_MATRIX_POSTBATCH = PROTEOMICS_DIR / "AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt"
MUTATION_SUBTYPES = BASE_DIR / "ALS_Mutation_Subtypes_with_ALS_Control.csv"

def run_analysis(protein_file, label):
    """Run DE analysis on given protein matrix"""
    print(f"\n{'='*80}")
    print(f"Analysis: {label}")
    print(f"{'='*80}")
    
    print(f"\n1. Loading data...")
    sys.stdout.flush()
    
    # Load protein matrix
    print(f"   Loading: {protein_file.name}")
    protein_df = pd.read_csv(protein_file, sep='\t')
    print(f"   ✓ {protein_df.shape[0]} proteins × {protein_df.shape[1]} columns")
    
    # Load mutation subtypes
    mutation_df = pd.read_csv(MUTATION_SUBTYPES)
    
    print(f"\n2. Identifying sample groups...")
    sys.stdout.flush()
    
    # Get sample columns
    all_cols = protein_df.columns.tolist()
    sample_cols = [c for c in all_cols if c not in ['Protein', 'nFragment', 'nPeptide']]
    
    # Extract participant IDs
    import re
    participant_map = {}
    for col in sample_cols:
        match = re.match(r'(CASE|CTRL)[_-]([A-Z0-9]+)', col)
        if match:
            prefix, code = match.groups()
            participant_id = f"{prefix}-{code}"
            participant_map[col] = participant_id
    
    # Create sample info
    sample_info = pd.DataFrame([
        {'Sample_ID': col, 'Participant_ID': participant_map.get(col)}
        for col in sample_cols
    ])
    
    # Merge with mutation data
    sample_info = sample_info.merge(
        mutation_df[['Participant_ID', 'Mutation_Subtype', 'SUBJECT_GROUP']], 
        on='Participant_ID', 
        how='left'
    )
    
    # Identify groups
    c9_samples = sample_info[sample_info['Mutation_Subtype'] == 'C9orf72']['Sample_ID'].tolist()
    ctrl_samples = sample_info[sample_info['SUBJECT_GROUP'] == 'Healthy Control']['Sample_ID'].tolist()
    
    print(f"   ✓ C9orf72: {len(c9_samples)} samples")
    print(f"   ✓ Control: {len(ctrl_samples)} samples")
    sys.stdout.flush()
    
    print(f"\n3. Running t-tests...")
    sys.stdout.flush()
    
    # Filter to selected samples
    selected_samples = c9_samples + ctrl_samples
    selected_cols = ['Protein'] + selected_samples
    analysis_df = protein_df[selected_cols].copy()
    
    # Get data matrices
    c9_cols = [c for c in c9_samples if c in analysis_df.columns]
    ctrl_cols = [c for c in ctrl_samples if c in analysis_df.columns]
    
    c9_data = analysis_df[c9_cols].values
    ctrl_data = analysis_df[ctrl_cols].values
    protein_names = analysis_df['Protein'].values
    
    print(f"   Testing {len(protein_names)} proteins...")
    results = []
    
    for i in range(len(protein_names)):
        if i % 500 == 0:
            print(f"   Progress: {i}/{len(protein_names)}", end='\r')
            sys.stdout.flush()
        
        protein = protein_names[i]
        c9_vals = c9_data[i, :]
        ctrl_vals = ctrl_data[i, :]
        
        # Remove NaN and zeros
        c9_clean = c9_vals[~np.isnan(c9_vals) & (c9_vals > 0)]
        ctrl_clean = ctrl_vals[~np.isnan(ctrl_vals) & (ctrl_vals > 0)]
        
        # Need at least 3 observations per group
        if len(c9_clean) < 3 or len(ctrl_clean) < 3:
            continue
        
        # Calculate means
        c9_mean = np.mean(c9_clean)
        ctrl_mean = np.mean(ctrl_clean)
        
        # Log2 fold change
        log2fc = np.log2(c9_mean + 1) - np.log2(ctrl_mean + 1)
        
        # T-test
        try:
            t_stat, p_val = stats.ttest_ind(c9_clean, ctrl_clean)
            
            results.append({
                'Protein': protein,
                'C9_mean': c9_mean,
                'Control_mean': ctrl_mean,
                'log2FC': log2fc,
                'pvalue': p_val,
                't_statistic': t_stat,
                'n_C9': len(c9_clean),
                'n_Control': len(ctrl_clean)
            })
        except:
            pass
    
    print(f"\n   ✓ Tested {len(results)} proteins with sufficient data")
    sys.stdout.flush()
    
    print(f"\n4. Calculating FDR...")
    sys.stdout.flush()
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    from scipy.stats import false_discovery_control
    results_df['FDR'] = false_discovery_control(results_df['pvalue'].values)
    
    # Sort by FDR
    results_df = results_df.sort_values('FDR')
    
    return results_df, len(c9_samples), len(ctrl_samples)


# Run both analyses
print("\n" + "="*80)
print("RUNNING BOTH ANALYSES")
print("="*80)

# Pre-batch corrected
results_pre, n_c9_pre, n_ctrl_pre = run_analysis(PROTEIN_MATRIX_PREBATCH, "PRE-batch corrected")
output_pre = OUTPUT_DIR / "DE_results_C9orf72_vs_Control_PREBATCH.tsv"
results_pre.to_csv(output_pre, sep='\t', index=False)
print(f"\n   ✓ Saved to: {output_pre}")

# Post-batch corrected
results_post, n_c9_post, n_ctrl_post = run_analysis(PROTEIN_MATRIX_POSTBATCH, "POST-batch corrected")
output_post = OUTPUT_DIR / "DE_results_C9orf72_vs_Control_POSTBATCH.tsv"
results_post.to_csv(output_post, sep='\t', index=False)
print(f"\n   ✓ Saved to: {output_post}")

# Compare results
print(f"\n{'='*80}")
print("COMPARISON: Pre-Batch vs Post-Batch Correction")
print(f"{'='*80}")

print(f"\nSample sizes:")
print(f"  Pre-batch:  {n_c9_pre} C9orf72, {n_ctrl_pre} Controls")
print(f"  Post-batch: {n_c9_post} C9orf72, {n_ctrl_post} Controls")

print(f"\nProteins tested:")
print(f"  Pre-batch:  {len(results_pre)}")
print(f"  Post-batch: {len(results_post)}")

sig_pre_05 = results_pre[results_pre['FDR'] < 0.05]
sig_pre_01 = results_pre[results_pre['FDR'] < 0.01]
sig_post_05 = results_post[results_post['FDR'] < 0.05]
sig_post_01 = results_post[results_post['FDR'] < 0.01]

print(f"\nSignificant proteins (FDR < 0.05):")
print(f"  Pre-batch:  {len(sig_pre_05)}")
print(f"  Post-batch: {len(sig_post_05)}")

print(f"\nSignificant proteins (FDR < 0.01):")
print(f"  Pre-batch:  {len(sig_pre_01)}")
print(f"  Post-batch: {len(sig_post_01)}")

print(f"\n{'-'*80}")
print("TOP 20 PROTEINS - PRE-BATCH CORRECTED")
print(f"{'-'*80}")
print(results_pre[['Protein', 'log2FC', 'pvalue', 'FDR']].head(20).to_string(index=False))

print(f"\n{'-'*80}")
print("TOP 20 PROTEINS - POST-BATCH CORRECTED")
print(f"{'-'*80}")
print(results_post[['Protein', 'log2FC', 'pvalue', 'FDR']].head(20).to_string(index=False))

# Find proteins significant in pre but not post
if len(sig_pre_05) > 0:
    pre_sig_proteins = set(sig_pre_05['Protein'].tolist())
    post_sig_proteins = set(sig_post_05['Protein'].tolist())
    
    lost_in_correction = pre_sig_proteins - post_sig_proteins
    gained_in_correction = post_sig_proteins - pre_sig_proteins
    
    print(f"\n{'-'*80}")
    print("BATCH CORRECTION IMPACT")
    print(f"{'-'*80}")
    print(f"Proteins significant in pre-batch only: {len(lost_in_correction)}")
    print(f"Proteins significant in post-batch only: {len(gained_in_correction)}")
    print(f"Proteins significant in both: {len(pre_sig_proteins & post_sig_proteins)}")
    
    if lost_in_correction:
        print(f"\nProteins lost after batch correction:")
        for p in list(lost_in_correction)[:10]:
            pre_row = results_pre[results_pre['Protein'] == p].iloc[0]
            print(f"  {p}: log2FC={pre_row['log2FC']:.3f}, pre-FDR={pre_row['FDR']:.4f}")

print(f"\n{'='*80}")
print("✓ Sensitivity analysis complete!")
print(f"{'='*80}")
