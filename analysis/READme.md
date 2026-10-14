# analysis/

Results for the sporadic ALS subtyping project. Every analysis batch has
its own subfolder holding its outputs, logs and a short README. This file
indexes all of them.

Code lives in `scripts/`. What the code produced lives here. Each batch
below names the scripts that produce it, so the two READMEs can be read
side by side.

## Contents

- [Folder layout](#folder-layout)
- [Results index](#results-index)
- [01 Backbone: axes and clustering](#01-backbone-axes-and-clustering)
- [02 Patient level and CCT2 classifier](#02-patient-level-and-cct2-classifier)
- [03 Culture composition sensitivity](#03-culture-composition-sensitivity)
- [04 LASSO panels](#04-lasso-panels)
- [05 Genotype DE and convergence](#05-genotype-de-and-convergence)
- [06 Mechanism: IPA, GO and axis components](#06-mechanism-ipa-go-and-axis-components)
- [07 Germline](#07-germline)
- [08 Clinical](#08-clinical)
- [09 Pathway robustness](#09-pathway-robustness)
- [10 Spatial validation](#10-spatial-validation)
- [11 Figures](#11-figures)
- [12 Tables](#12-tables)
- [Superseded results](#superseded-results)
- [Adding a new batch](#adding-a-new-batch)

## Folder layout

```
analysis/
├── README.md                        this file
├── 01_backbone/
├── 02_patient_level/
├── 03_composition_sensitivity/
├── 04_lasso/
├── 05_genotype_de_convergence/
├── 06_mechanism/
├── 07_germline/
├── 08_clinical/
├── 09_pathway_robustness/
├── 10_spatial/
├── 11_figures/
├── 12_tables/
└── archive/                         superseded results, kept for the record
```

Each batch folder holds:

- `README.md`: question, scripts, date run, landmarks, status
- `results/`: the `.rds`, `.csv` and `.tsv` outputs
- `plots/`: diagnostic plots not in the paper
- `logs/`: console output from the audited run (`sink()` to a `.txt`)

Large pinned `.rds` files (anything over 50 MB, for example
`backbone_FINAL.rds`) are not committed. Their README lists the path on
the cluster and the MD5 checksum instead.

## Results index

One line per batch. Status is `final` (feeds the manuscript),
`supporting` (in supplementary material or methods only) or
`superseded` (moved to `archive/`).

| Batch | Question | Main output | Headline result | Status |
|---|---|---|---|---|
| 01 | Do sporadic patients separate on three pathway axes? | `backbone_FINAL.rds` | Three states at k = 3 (PAC 0.299) | final |
| 02 | Does one gene mark S-ALS-3? | `patient_level_FINAL.rds` | CCT2 AUC 0.892 | final |
| 03 | Is the partition a culture artefact? | `diff_quality_sensitivity_FINAL.rds` | Partly for S-ALS-1; CCT2 effect holds | final |
| 04 | How small a gene panel recovers each state? | `lasso_FINAL.rds` | 12 to 16 genes; CCT2 for S-ALS-3 only | final |
| 05 | How do C9orf72, SOD1 and sporadic relate? | `genotype_de_corrected_FINAL.rds` | C9orf72 closer to sporadic at pathway level | final |
| 06 | What biology sits behind S-ALS-3? | IPA exports, GO tables | MYC and mTORC1 translation program | final |
| 07 | Is the state encoded in germline DNA? | `genomic_burden_FINAL.rds` | No test survives correction | final |
| 08 | Does CCT2 track clinical decline? | `downstream_analyses_FINAL.rds` | ALSFRS-R link weakens after composition adjustment | final |
| 09 | Are the axes Hallmark-specific? | `database_independence_FINAL.rds` | GO and Reactome reproduce them | supporting |
| 10 | Does CCT2 track the cell cycle axis in tissue? | `spatial_results_FINAL.rds` | rho 0.62 in frontal cortex | final |
| 11 | Figures 1 to 9 | PDF + 300 dpi TIFF | | final |
| 12 | Tables 1, 2, S1 to S4 | CSVs | | final |

Cohort throughout: 644 RNA-seq samples from 619 participants (185 / 46 /
13 / 400 samples; 173 control, 43 C9orf72, 12 SOD1, 391 sporadic
participants), after removing 98 replicated reference-line samples.

---

## 01 Backbone: axes and clustering

**Scripts:** `01_rebuild_backbone_AUDITED.R`
**Output:** `backbone_FINAL.rds`

Three Hallmark-derived axes scored by ssGSEA on the sex- and
SV-residualised matrix, then ConsensusClusterPlus on the 400 sporadic
samples.

| Landmark | Value |
|---|---|
| Axis sizes (metabolic / cell cycle / inflammation) | 512 / 610 / 538 genes |
| PAC, k = 2 / 3 / 4 / 5 | 0.310 / **0.299** / 0.368 / 0.356 |
| Subtype sizes, samples | 89 / 204 / 107 |
| Cross-axis Spearman (met-infl / met-cc / infl-cc) | 0.24 / 0.61 / 0.60 |
| PCA variance, PC1 / PC2 | 66.3% / 25.4% |
| SVs used | SV1 to SV6 from the 742-sample fit |

Round 2 clustering (marker-residualised) was dropped on 13 August 2026.
See [Superseded results](#superseded-results).

## 02 Patient level and CCT2 classifier

**Scripts:** `02_patient_level_AUDITED.R`, `24_figure5_cct2.R`
**Outputs:** `patient_level_FINAL.rds`, `figure5_stats.rds`

Earliest timepoint per participant, 391 patients.

| Landmark | Value |
|---|---|
| Patients per subtype | 85 / 204 / 102 |
| CCT2 AUC, S-ALS-3 vs S-ALS-2 | 0.892 (95% CI 0.851 to 0.927) |
| t-test / Wilcoxon | P = 2.7e-32 / P = 5.6e-29 |
| S-ALS-3 above the S-ALS-2 maximum | 28.4% |
| Dip test within S-ALS-3 | P = 0.679 |
| AUC without the top 15 patients | 0.873 |
| CCT2 vs control, VST (S1 / S2 / S3) | -0.060 (P 0.022) / -0.106 (P 7.0e-10) / +0.138 (P 3.6e-12) |
| All sporadic vs control | -0.033 (P 0.044) |

Protein (`08_cct2_protein_test.R`, `cct2_protein_by_subtype.csv`):

| Landmark | Value |
|---|---|
| Patients with proteomics | 265 (59 / 144 / 62) |
| Kruskal-Wallis / S3 vs S2 Wilcoxon | P = 0.005 / P = 0.033 |
| RNA-protein Spearman | rho 0.21, P = 8e-4 |
| TRiC 8-subunit composite | P = 0.13 (not specific to the complex) |

## 03 Culture composition sensitivity

**Scripts:** `03_differentiation_sensitivity.R`, `03b_marker_rescue_rerun.R`
**Output:** `diff_quality_sensitivity_FINAL.rds`

| Landmark | Value |
|---|---|
| Differentiation quality across subtypes | eta squared 0.382, driven by S-ALS-1 |
| ARI after S100B residualisation | 0.565 |
| ARI after full glial + DQ panel | 0.273 |
| CCT2 S3 effect: unadjusted / +3 composites / +full panel | +0.244 / +0.217 / +0.231 (all P < 1e-23) |

Conclusion: S-ALS-1 is composition-confounded and treated as
exploratory. S-ALS-3 vs S-ALS-2 is the primary contrast.

To confirm on the corrected cohort: Tukey S-ALS-2 vs S-ALS-3 difference
in differentiation quality (0.37, P = 0.0009 in an earlier draft).

## 04 LASSO panels

**Scripts:** `04_lasso_AUDITED.R`, `onevsrestLASSO.R`
**Outputs:** `lasso_FINAL.rds`, Supplementary Table S1

| Landmark | Value |
|---|---|
| Candidate genes | 1,470 Hallmark axis genes |
| Panel size per subtype | 12 to 16 genes |
| Nested AUC | 0.936 to 0.970 |
| CCT2 bootstrap selection, S-ALS-3 | 74.7%; never selected for S-ALS-1 or S-ALS-2 |

## 05 Genotype DE and convergence

**Scripts:** `18_full_recovery_and_de_AUDITED.R`, `14_tals_rebuild.R`,
`16a_pin_proteomics_inputs.R`, `16a_preamble.R`,
`16_convergence_corrected.R`
**Outputs:** `genotype_de_corrected_FINAL.rds`,
`convergence_results_corrected/`, `tals_de_results.rds`,
`proteomics_inputs_FINAL.rds`

Design `~SEX + 6 clean SVs + Mutation_Subtype` on the 644-sample refit.
SV3, SV4 and SV5 were excluded (genotype ANOVA P = 0.004, 0.022, 0.023).
SV10 passed (P = 0.25) but was dropped to keep six SVs, matching the
backbone. Kept: SV1, SV2, SV6, SV7, SV8, SV9.

| Comparison | C9-Spor | C9-SOD1 | SOD1-Spor |
|---|---|---|---|
| Gene logFC Spearman, 24,253 genes | 0.326 | 0.026 | 0.263 |
| RNA Hallmark NES Spearman | 0.402 | -0.070 | 0.045 |
| Protein Hallmark NES Spearman | 0.319 | -0.282 | 0.108 |
| Target ALS bulk motor cortex | 0.751 | 0.487 | 0.751 |

Top-500 C9orf72 anchor: 94.6% same direction in sporadic (rho 0.803),
56.6% in SOD1 (rho 0.295).

GSE56500 laser-capture motor neurons: C9 vs sporadic gene rho 0.59,
Hallmark NES rho 0.96. Directional only (3 per group).

Data-driven GO on the three contrasts (padj < 0.05, |log2FC| > 0.3):
38 / 25 / 231 genes, no significant GO term in any contrast
(Supplementary Note 2).

## 06 Mechanism: IPA, GO and axis components

**Scripts:** `07_ipa_prep_AUDITED.R`, final IPA DE script, GO and GSVA
code (see batch README)
**Outputs:** IPA noSexChr Upstream Regulators and Canonical Pathways
exports, GO tables, `figure6_stats.rds`

DE: S-ALS-3 vs S-ALS-2, 311 samples (107 / 204), design
`~SEX + SV1..SV6 + subtype`. IPA version 01-23-01, Core Analysis,
Ingenuity Knowledge Base (genes only), direct and indirect,
experimentally observed. Universe 24,021 genes without sex chromosomes.
Cutoff FDR < 0.001 and |log2FC| >= 0.3 gives 2,865 genes (771 up,
2,094 down).

| Regulator or pathway | z | P |
|---|---|---|
| MYC | +2.14 | 6.0e-28 |
| MYCN | +5.58 | 2.2e-19 |
| RPTOR | +2.88 | 1.4e-5 |
| torin1 / sirolimus | -3.06 / -3.93 | 2.5e-28 / 1.3e-21 |
| LARP1 | -4.82 | 1.9e-23 |
| RICTOR | -4.70 | 3.1e-9 |
| Translation elongation (pathway) | +6.56 | |
| GCN2 / EIF2AK4 amino acid deficiency (pathway) | +6.40 | 7.3e-15 |

Use the pathway z for GCN2 (+6.40), not the upstream regulator z for
EIF2AK4 (+1.30).

GO-BP (clusterProfiler 4.18.4, same list): up in S-ALS-3, cytoplasmic
translation (54 genes, adj P 2.8e-21), sterol biosynthesis, aerobic
respiration, glycolysis. Down: miRNA silencing (adj P 3.8e-19),
mitotic division and spindle checkpoint, cilium assembly.

Cell cycle axis components (GSVA 2.4.9, unresidualised VST,
S-ALS-3 vs S-ALS-2):

| Hallmark set | Direction | Wilcoxon P |
|---|---|---|
| MYC targets V1 | up | 2.1e-11 |
| DNA repair | up | 0.002 |
| E2F targets | none | 0.90 |
| G2M checkpoint | down | 3.6e-4 |

The cell cycle signal in S-ALS-3 is a MYC target signal, not mitosis.

## 07 Germline

**Scripts:** `genomic_preamble.R`, `genomic_analysis.R`,
`tric_burden.R`, `TRICBurdenMYCBurden.R`, `ciseQTLTesting.R`
**Outputs:** `genomic_burden_FINAL.rds`, `tric_burden_results/`,
`myc_burden_results/`, `cct2_cis_eqtl_results.csv`, Supplementary Table S2

360 sporadic patients with WGS (280 in S-ALS-3 vs S-ALS-2).

| Test | Result |
|---|---|
| Rare damaging burden, ANOVA / Kruskal-Wallis | P = 0.864 / P = 0.05 |
| Rare LoF burden | P = 0.30 |
| 34-gene ALS panel carriers | 242 / 360 (67.2%), chi-square P = 0.95 |
| CCT2 rare damaging carriers | 0 of 360 |
| CCT2 cis-eQTL, 3,782 variants | minimum FDR q = 0.712 |
| Top 50 cis variants in ALS GWAS | 42 matched, minimum P = 0.116 |
| TRiC carriers / burden | P = 0.19 / P = 0.35 |
| MYC target burden / carriers | P = 0.85 / P = 0.73 |
| MYC per-gene scan, top hit XPOT | 5/88 vs 0/192, P = 0.003, q = 0.362 |
| Genomic k-means / blind consensus vs subtype | ARI -0.005 / 0.003 |

## 08 Clinical

**Scripts:** `05_downstream_ladders_AUDITED.R`,
`06_clinical_screen_AUDITED.R`, `27_figure8_clinical.R`
**Outputs:** `downstream_analyses_FINAL.rds`,
`clinical_screen_combined_FINAL.rds`, `figure8_stats.rds`,
Supplementary Table S3

All models adjust for age at onset, sex and SV1 to SV6. CCT2 effects
are per SD (SD of CCT2 VST = 0.193).

ALSFRS-R slope, points per month, n = 244:

| Model | beta | P |
|---|---|---|
| Unadjusted | -0.135 | 0.002 |
| Adjusted | -0.137 | 0.005 |
| Excluding slope < -2 (n = 233) | -0.092 | 0.018 |
| Trimmed 1st/99th pct (n = 238) | -0.115 | 0.008 |
| + subtype | -0.122 | 0.039 |
| + subtype + composition | -0.084 | 0.21 |
| + subtype + composition + glia | -0.076 | 0.27 |

log10 NfL, n = 297: CCT2 +0.039 (P 0.045), falls to +0.008 (P 0.75)
with subtype. Subtype omnibus P = 0.015; with both terms, subtype
P = 0.11 and CCT2 P = 0.75. Not separable.

With NfL in the ALSFRS-R model (n = 188): adjusted R squared 0.058 to
0.262; CCT2 -0.142 (P 0.010) to -0.099 (P 0.047).

Disease duration, 182 deceased patients: ANOVA P = 0.20,
Kruskal-Wallis P = 0.08.

Endpoint screen: 25 screened, 23 modeled, none survives FDR (minimum
q = 0.12). CCT2 nominal hits are all decline rates: grip strength
(P 0.003), portal ALSFRS-R (P 0.006), vital capacity (P 0.039).
Pathway by endpoint matrix, 30 tests, minimum FDR 0.213.

## 09 Pathway robustness

**Scripts:** `30_pathway_preamble.R`, `31_pathway_endpoint_screen.R`,
`32_database_independence.R`
**Outputs:** `database_independence_FINAL.rds`, `phase*_*.csv`

| Landmark | Value |
|---|---|
| GO-BP axis correlation with Hallmark | 0.93 / 0.96 / 0.96 |
| Reactome axis correlation with Hallmark | 0.88 / 0.97 / 0.79 |
| Partition ARI vs Hallmark, GO-BP / Reactome | 0.494 / 0.594 |
| GO-BP vs Reactome partition ARI | 0.614 |

Supplementary Note 2.

## 10 Spatial validation

**Scripts:** `spatialdata.R`, `28_figure9_spatial.R`
**Outputs:** `spatial_results_FINAL.rds`, `figure9_stats.rds`

Target ALS Visium, frontal cortex only, top quartile of layer 5b
projection neuron spots pseudobulked per subject. 28 ALS, 11 control.

| Landmark | Value |
|---|---|
| CCT2 ALS vs control | lower in ALS, Wilcoxon P = 0.024 |
| CCT2 vs cell cycle axis, within ALS | rho 0.62, P = 5.2e-4; Pearson r 0.53, P = 0.004 |
| Dropping top 1 / top 3 CCT2 patients | rho 0.70 / 0.64 |
| CCT2 removed from the axis | rho 0.62 (ranks unchanged) |
| CCT2 vs metabolic / inflammation | rho 0.15 (P 0.43) / -0.38 (P 0.045, ns after correction) |

Motor cortex (4 samples) is descriptive only.

## 11 Figures

**Scripts:** `00_figure_common.R`, `20_figure1_cohort.R` to
`28_figure9_spatial.R`

Each figure is saved as a font-embedded vector PDF and a 300 dpi LZW
TIFF, 180 mm wide. Subtype palette `#7E57C2 / #90A4AE / #F28E2B`.
Anything rendered before 5 October 2026 uses the old palette and should
be regenerated.

| Figure | Content | File |
|---|---|---|
| 1 | Study design and cohort | `Figure1.pdf/.tiff` |
| 2 | Genotype similarity | `Figure2.pdf/.tiff` |
| 3 | Axes | `Figure3.pdf/.tiff` |
| 4 | Consensus clustering | `Figure4.pdf/.tiff` |
| 5 | CCT2 classifier and protein | `Figure5.pdf/.tiff` |
| 6 | Upstream regulators and pathways | `Figure6.pdf/.tiff` |
| 7 | Germline | `Figure7.pdf/.tiff` |
| 8 | Clinical screen and robustness | `Figure8.pdf/.tiff` |
| 9 | Spatial | `Figure9.pdf/.tiff` |

`Figure6_panelD_mechanism` (schematic) is not in the paper.

## 12 Tables

**Script:** `make_tables.R`

| Table | Content | File |
|---|---|---|
| 1 | Cohort characteristics | `Table1_cohort.csv` |
| 2 | Subtype characteristics | `Table2_subtypes.csv` |
| S1 | LASSO panels | `TableS1_lasso_summary.csv`, `TableS1_lasso_panels.csv` |
| S2 | Germline tests | `TableS2_genomic_tests.csv` |
| S3 | Clinical endpoint screen, 23 endpoints | `TableS3_clinical_screen.csv` |
| S4 | Axis composition | `TableS4_axis_composition.csv`, `TableS4b_axis_genes.csv` |

Supplementary tables are numbered in order of first citation. Any
console line starting `TYPED:` marks a value with no saved object
behind it; check those by hand.

---

## Superseded results

Kept in `archive/` so the record is complete. Do not cite any of these.

| Result | Why it was replaced | Replacement |
|---|---|---|
| Convergence on the pre-correction cohort (RNA NES 0.761 / -0.077 / 0.104) | Reference line counted as 98 controls | Batch 05 |
| RNA NES 0.804 / +0.170 / 0.348 | Not reproducible from any object on disk; retired Aug 2026 | Batch 05 |
| Round 2 clustering (269 / 130 / 1) | Degenerate partition, dropped 13 Aug 2026 | Batch 01 round 1 |
| PRSet pathway scores | GWAS on GRCh37, genotypes on GRCh38; about 1% of positions matched | Not reported |
| Spatial on all 43 samples (P 0.099, rho 0.55) | Mixed regions, motor cortex n = 4 | Batch 10 frontal only |
| GO drill-down from June (226 vs 76) | Built on `IPAPrepare.R`, pre-correction | Batch 06 rerun |
| Data-driven GO from June (5 developmental terms) | Pre-correction contrasts | Batch 05 rerun, 0 terms |
| Fig. 8c without SVs (+ subtype P = 0.056) | Missing SV covariates | Batch 08 |
| GCN2 z = +1.30 | Upstream regulator value, not the pathway | Batch 06, +6.40 |
| Earlier "77% above S-ALS-2 maximum" | Draft error | 28.4% |
| MOFA, NMF, GMM, DDRTree, WGCNA explorations | Early approaches, replaced by the axis framework | Not reported |

## Adding a new batch

1. Create `analysis/NN_short_name/` with the next free number.
2. Copy this template into its `README.md`:

```markdown
# NN short name

**Question:** one sentence
**Scripts:** scripts/...
**Run:** YYYY-MM-DD, fresh session, by initials
**Inputs:** pinned files read
**Outputs:** files in results/
**Status:** exploratory | supporting | final | superseded

## Landmarks
| Landmark | Value |
|---|---|

## Notes
Decisions, caveats, anything a reviewer would ask.
```

3. Add one row to the [Results index](#results-index) and a section
   below it.
4. If the batch replaces an earlier one, move the old folder to
   `archive/` and add a row to [Superseded results](#superseded-results).
