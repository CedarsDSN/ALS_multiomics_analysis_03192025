# scripts/

All analysis code for the sporadic ALS subtyping project. This README
indexes every script, says what it reads and writes, and gives the run
order.

Code lives here. What the code produced lives in `analysis/`.

## Contents

- [Run order](#run-order)
- [Master index](#master-index)
- [Core pipeline](#core-pipeline)
- [Genotype DE and convergence](#genotype-de-and-convergence)
- [Clinical](#clinical)
- [Germline](#germline)
- [Pathway](#pathway)
- [Spatial](#spatial)
- [Figures and tables](#figures-and-tables)
- [Conventions](#conventions)
- [Pinned inputs](#pinned-inputs)
- [Environment](#environment)
- [Archive](#archive)

## Run order

```
01, 02                       required first, in order
  |
  +-- 03, 03b, 04            core, independent of each other
  +-- 05, 06                 clinical
  +-- 18 -> 16a -> 16        genotype DE and convergence (14 feeds 16)
  +-- 07                     IPA prep
  +-- 30 -> 31, 32           pathway
  +-- germline arm           own working directory, see below
  +-- spatialdata.R          spatial
        |
        +-- 00_figure_common -> 20 to 28     figures, last
        +-- make_tables.R                    tables, last
```

Scripts 01 and 02 build the objects everything else reads. Nothing
runs before them. Figures and tables run last, because they only
consume pinned objects.

## Master index

| Script | Group | Reads | Writes |
|---|---|---|---|
| `01_rebuild_backbone_AUDITED.R` | core | `cohort_inputs_FINAL.rds` | `backbone_FINAL.rds` |
| `02_patient_level_AUDITED.R` | core | backbone, clinical, NfL | `patient_level_FINAL.rds` |
| `03_differentiation_sensitivity.R` | core | backbone, patient level | `diff_quality_sensitivity_FINAL.rds` |
| `03b_marker_rescue_rerun.R` | core | as 03, plus counts | marker-rescue adjudication |
| `04_lasso_AUDITED.R` | core | backbone, patient level | `lasso_FINAL.rds` |
| `05_downstream_ladders_AUDITED.R` | clinical | patient level, clinical, PRS | `downstream_analyses_FINAL.rds` |
| `06_clinical_screen_AUDITED.R` | clinical | patient level, portal CSVs | `clinical_screen_combined_FINAL.rds` |
| `AllClinicalVariableTesting_newlabels.R` | clinical | sourced by 06 | the five screen waves |
| `07_ipa_prep_AUDITED.R` | pathway | backbone, counts | `ipa_inputs_FINAL/`, `*_noXY.tsv` |
| `08_cct2_protein_test.R` | convergence | protein matrix, patient level | `cct2_protein_by_subtype.csv` |
| `14_tals_rebuild.R` | convergence | Target ALS bulk, metadata | `tals_de_results.rds` |
| `16_convergence_corrected.R` | convergence | 18, 16a, 14 | `convergence_results_corrected/` |
| `16a_pin_proteomics_inputs.R` | convergence | proteomics TSVs | `proteomics_inputs_FINAL.rds` |
| `16a_preamble.R` | convergence | 18, convergence pin | six objects for 16 |
| `18_full_recovery_and_de_AUDITED.R` | convergence | `cohort_inputs_FINAL.rds` | `genotype_de_corrected_FINAL.rds` |
| `genomic_preamble.R` | germline | backbone, patient level | supplies objects |
| `genomic_analysis.R` | germline | ANNOVAR reports | `genomic_burden_FINAL.rds` |
| `tric_burden.R` | germline | preamble, burden | `tric_burden_results/tric_models.rds` |
| `TRICBurdenMYCBurden.R` | germline | preamble, burden | `myc_burden_results/myc_models.rds` |
| `ciseQTLTesting.R` | germline | genotypes, GWAS base | `cct2_cis_eqtl_results.csv` |
| `30_pathway_preamble.R` | pathway | pinned files only | supplies objects |
| `31_pathway_endpoint_screen.R` | pathway | 30 | `phase2_*`, `phase6_*`, `phase7_*` CSVs |
| `32_database_independence.R` | pathway | 30 | `database_independence_FINAL.rds` |
| `spatialdata.R` | spatial | Target ALS spatial, backbone | `spatial_results_FINAL.rds` |
| `00_figure_common.R` | figures | nothing | loader, palette, path guards |
| `20_figure1_cohort.R` | figures | backbone, patient level | Figure 1 |
| `21_figure2_genotype_de.R` | figures | 18, 16, 16a | Figure 2 |
| `22_figure3_axes.R` | figures | backbone, patient level | Figure 3 |
| `23_figure4_consensus.R` | figures | backbone | Figure 4 |
| `24_figure5_cct2.R` | figures | patient level | Figure 5 |
| `25_figure6_mechanism.R` | figures | IPA exports | Figure 6 |
| `26_figure7_genomic.R` | figures | germline outputs | Figure 7 |
| `27_figure8_clinical.R` | figures | 05, 06, patient level | Figure 8 |
| `28_figure9_spatial.R` | figures | `spatial_results_FINAL.rds` | Figure 9 |
| `make_tables.R` | tables | all pins | Tables 1, 2, S1 to S4 |

35 scripts. Seven more are in `archive/` and must not be run.

## Core pipeline

```r
source("01_rebuild_backbone_AUDITED.R")
source("02_patient_level_AUDITED.R")
source("03_differentiation_sensitivity.R")   # optional
source("03b_marker_rescue_rerun.R")          # optional, follows 03
source("04_lasso_AUDITED.R")
```

**01** scores three Hallmark-derived axes (metabolic, cell cycle,
inflammation) by ssGSEA on residualised expression, then partitions the
sporadic samples with ConsensusClusterPlus at k = 3. Settings are
maxK 3, reps 1000, pItem 0.8, km on euclidean, seed 42.

Two things to know before editing it. `axis_sets`, the three gene sets
the whole paper rests on, is defined here and saved into the backbone.
It used to exist only as a session object, so the core construct was
never under version control. Change it here and you change the paper.

A second residualisation round was dropped on 13 August 2026 because it
produced a degenerate 269 / 130 / 1 split. The script keeps `d2 <- d1`
so the shape of the old code stays legible, but round 2 does not run.
Do not reinstate it.

**02** collapses 400 sporadic samples to 391 patients by earliest
timepoint, then joins clinical variables, NfL and CCT2. Earliest
timepoint is the canonical deduplication rule project-wide. The NfL
file is required; an earlier version skipped it silently and produced
a patient table that looked complete.

**03** tests whether the subtyping is a differentiation-quality
artefact. The answer it prints is partly yes for S-ALS-1 and no for the
S-ALS-3 against S-ALS-2 contrast where CCT2 lives. **03b** follows up
on the marker panel; the glial markers fall below the expression
filter.

**04** runs one-vs-rest LASSO with bootstrap stability selection and
selects CCT2 for S-ALS-3. `axis_sets` comes from the backbone here,
not from the session.

## Genotype DE and convergence

```r
source("18_full_recovery_and_de_AUDITED.R")
source("14_tals_rebuild.R")             # independent of 18
source("16a_pin_proteomics_inputs.R")   # independent of both
source("16a_preamble.R")
source("16_convergence_corrected.R")
source("08_cct2_protein_test.R")
```

**The surrogate variable question.** This is the decision most likely
to look like an inconsistency, so it is written down.

Script 01 uses SV1 to SV6 from an SVA fit on the full 742 samples,
including the batch-replicated reference control line. Script 18 uses
a separate refit on the corrected 644 with the genotype-correlated
surrogate variables (3, 4, 5) excluded, leaving 1, 2, 6, 7, 8 and 9.

The two are different on purpose. Refitting on the corrected cohort
alone gave surrogate variables correlated 0.50 to 0.65 with the
biological axes, which dissolved the subtype structure, so the backbone
keeps the 742-sample fit. The DE layer has the opposite requirement,
surrogate variables carrying no genotype signal, which that fit cannot
guarantee. Do not try to reconcile them.

Design for the DE layer is `~SEX + 6 clean SVs + Mutation_Subtype`.
Two sensitivity fits bracket the canonical values, one without SVA and
one on the anchor SVs. Reference only, never primary.

**The metric.** Gene-level cross-genotype similarity is the Spearman
correlation of `log2FoldChange` across all shared genes. A
Wald-statistic version exists in the codebase and is secondary. It must
not appear in the paper or on slides. Ordering is always C9 to
Sporadic, C9 to SOD1, SOD1 to Sporadic.

**Object case, which costs an hour if missed.** `res_C9`, `res_SOD1`,
`res_Sporadic` uppercase are RNA from DESeq2 with `log2FoldChange` and
`padj`. `res_c9`, `res_sod1`, `res_spor` lowercase are protein from
limma with `gene`, `logFC`, `t`, `P.Value`, `adj.P.Val`. They differ
only by letter case. Check the columns, not the name.

**CCT2 in the protein matrix.** CCT2 is absent from the
symbol-collapsed matrix. Its row is `"1/sp|P78371|TCPB_HUMAN"`, and
`TCPB_HUMAN` defeats symbol mapping. `16a_pin_proteomics_inputs.R`
rescues it by UniProt accession, so take CCT2 from the pin rather than
re-deriving it.

**Output directories.** `convergence_results_corrected/` is current.
`~/Documents/convergence_results/` holds a superseded pre-correction
run and is blocked by absolute path in `00_figure_common.R`, with one
named exception for `genotype_de_corrected_FINAL.rds`, which
legitimately lives there. A script resolving a convergence file into
the unsuffixed directory is a bug, not a fallback.

`16a_preamble.R` is the one live script never audited end to end. Read
its path-resolution block before running it from a clean session.

## Clinical

```r
source("05_downstream_ladders_AUDITED.R")
source("06_clinical_screen_AUDITED.R")   # sources the wave monolith itself
```

**05** fits the NfL and ALSFRS-R outcome ladders. Scope, recorded in
the header because the analysis is defensible either way: ladders, BIC
and head-to-head run on all 391 sporadic patients with three-level
subtype, while the PRS sections keep the S2 and S3 clean subset per the
PRSet specification. The script prints a three-specification
robustness table so a reader can see the surrogate variables stabilise
the ALSFRS-R result rather than create it.

**06** is a thin runner. It pins the working directory, sources
`AllClinicalVariableTesting_newlabels.R`, combines the five waves and
gates the result. Do not source the monolith directly.

Two details in the 06 gate are not obvious. Endpoint labels are not
unique, since "Creatinine baseline" appears in both wave 1 and wave 4,
so the gate orders and compares rather than matching on label. And
`setequal` resolves into the wrong namespace with `conflicted`
attached, so the script calls `base::setequal` explicitly. Leave both
alone.

## Germline

These read the ANNOVAR exon reports with `list.files(".")`, so **set
the working directory to the folder holding the
`ALS.chrN.anno_report_exon.txt` files** before starting.

```r
setwd("/path/to/annovar/reports")
source("genomic_preamble.R")      # must run first
source("genomic_analysis.R")
source("tric_burden.R")
source("TRICBurdenMYCBurden.R")
source("ciseQTLTesting.R")
```

`genomic_preamble.R` bridges the current backbone labels to the object
names the older genomic code expects, so nothing downstream runs on a
stale label set.

Rare means MAF below 0.01 on the ANNOVAR "Max Freq of 3" column.
Damaging classes are `nonsynonymous_SNV`, `stopgain`, `stoploss`,
`frameshift_deletion`, `frameshift_insertion`. Both are constants at
the top of the burden scripts.

**Zero-carrier tests are not estimable.** With no carriers in a
subgroup, `chisq.test` silently degrades to a goodness-of-fit test on a
single-column table and returns a meaningless p-value. CCT2 has zero
carriers in some strata. Both burden scripts guard against this. Do not
remove the guards.

**Genome build.** The van Rheenen summary statistics (GCST90027164)
ship as GRCh37 while AnswerALS genotypes are GRCh38, a +393,780 bp
offset on chr12 that leaves about 1 percent of positions matching by
coincidence. Use the harmonised GRCh38 release, the `.h.tsv.gz` file,
and take its `hm_*` columns as a set. The lookup also checks allele
orientation. This mismatch is why the PRSet analysis was dropped.

**Deliberate overwrite.** The final `saveRDS` in `genomic_analysis.R`
moves `burden_subtype` from 387 rows to 360. The comment on that line
says so. It is a pseudoreplication fix, not a lost subset.

## Pathway

```r
source("07_ipa_prep_AUDITED.R")     # independent

source("30_pathway_preamble.R")     # required before 31 and 32
source("31_pathway_endpoint_screen.R")
source("32_database_independence.R")
```

**07** writes the DE tables IPA consumes. The sex-chromosome exclusion
happens inside IPA through the analysis settings, not by filtering the
saved data. The `_noXY.tsv` files exist to make that setup
reproducible. The gate runs at tolerance rather than bit-identity,
because the annotation round trip is not exactly reversible.

**32** is the reviewer-facing robustness check: rebuild the three axes
from GO-BP and from Reactome, recluster, compare to the Hallmark
partition. `RUN_CLEANED` stays `FALSE`; it guards the abandoned
round-2 residualisation. The script prints cluster sizes before any
ARI, because a near-zero ARI against a partition with one singleton is
degenerate rather than a negative result. `ConsensusClusterPlus` is
called with `plot = "png"` after `graphics.off()`, since `plot = NULL`
inherits whatever device is current and throws once anything has
written a PDF. `rrvgo` and `GOSemSim` are optional.

**31** runs the pathway by clinical endpoint matrix and three further
phases. Patient-level end to end; the predecessor joined a 391-row
patient frame to sample-level ssGSEA scores and pseudoreplicated
silently. Scores are computed on `expr_resid`, matching how the axes
were built. Phase 6's centroid distances are descriptive and should
not be reported as a test.

## Spatial

```r
source("spatialdata.R")
```

Frontal cortex only. Motor cortex contributes four samples, two per
group, too few for a comparison, and pooling regions mixes different
baseline expression. An earlier all-43-sample version is underpowered
and its numbers should not be quoted.

The circularity check rebuilds the cell cycle axis without CCT2 and
repeats the correlation. The result is unchanged, which means CCT2
contributes negligibly to the composite. It is not two independent
estimates agreeing, since they share nearly all their input genes.

`spatialdata.R` and the archived `11_spatial.R` both write
`spatial_results_FINAL.rds`. The published numbers come from
`spatialdata.R`. If both end up on the path, whichever runs second
silently wins.

## Figures and tables

```r
source("00_figure_common.R")   # first, every session
source("21_figure2_genotype_de.R")   # or any other
source("make_tables.R")
```

`00_figure_common.R` provides `p_()` and `load_pin()` for paths,
`get_backbone()`, `get_patients()`, `get_patients_full()` and
`get_group()` for pinned objects, `fmt_p()`, `save_fig()`,
`check_fonts()`, `figure_readiness()`, `report_structure()` and
`ALL_INPUTS`.

**Palette.**

```r
SCOL <- c("S-ALS-1" = "#7E57C2", "S-ALS-2" = "#90A4AE", "S-ALS-3" = "#F28E2B")
```

Plus `NCOL` for neutral categories. This replaced
`#E64B35` / `#4DBBD5` / `#00A087`, which survives only as a fallback in
`26_figure7_genomic.R`. Any figure rendered before 5 October 2026 uses
the old palette and needs regenerating.

**Font embedding.** `ggsave` to `.pdf` uses R's default device, which
does not embed fonts. `pdffonts` shows Helvetica with `emb: no`,
non-Mac viewers substitute and blur the text, and Nature requires
embedded fonts. `device = cairo_pdf` does not fix it here: cairo fails
to load because R's `cairo.so` needs XQuartz libraries that are not
installed, and R falls back to the default device silently. `save_fig()`
therefore drives the native macOS device,
`grDevices::quartz(file = ..., type = "pdf", ...)` then `print(plot)`
then `dev.off()`, which embeds font subsets with nothing to install,
and also writes a 300 dpi TIFF. Verify with `pdffonts` showing
`emb: yes` and `sub: yes`. Output is 510 by 316 pt, which is 180 by
111 mm at double column.

**Per-figure notes.** 21 resolves its own inputs rather than using
`p_()`, with a check that stops if it lands in the superseded
convergence directory; that is deliberate and should not be
refactored. 22 reads gene-set sizes from `backbone$axis_sets` and
states the cross-axis correlations as computed, because the hardcoded
sizes and the "all below 0.5" caption were both wrong. 24 consumes the
earliest-timepoint patient table rather than averaging across a
patient's samples. 25 needs the two IPA export CSVs from the noSexChr
analysis and stops rather than returning NA, because the old version
produced a blank schematic that looked like a rendering bug. 26 panel b
must read the post-fix eQTL CSV. 27 takes survival rows from
`downstream_analyses_FINAL.rds` rather than the session.

**make_tables.R** builds Table 1, Table 2 and Supplementary S1 to S4,
numbered by order of first citation as Nature Communications requires.
Where a value has no saved object behind it, the console prints
`TYPED:`. Those are the only places a number can drift out of sync with
the analysis, so check them by hand before submission.

## Conventions

Three rules hold in every script.

1. **Inputs are pinned.** Each script starts with `readRDS()` of a
   named file. No script depends on an object that happens to be in
   the session.
2. **Outputs are saved.** Each script ends with `saveRDS()` to a named
   file.
3. **Reruns are gated.** Before overwriting a canonical `_FINAL` file,
   a script recomputes its results and compares them to what is on
   disk. On a mismatch it writes `*_AUDITCHECK.rds` and leaves the
   canonical file untouched.

All `_FINAL` files are treated as immutable. If a script says it wrote
an `_AUDITCHECK` file, stop and work out why before running anything
downstream.

To run anything: fresh R session with no workspace restore, `source()`
top to bottom with no manual steps and no console detours between the
seed and any stochastic step, then compare the printed landmarks to
the values in the script header.

Scripts with `_AUDITED` in the name have been verified end to end from
a clean session. Before sourcing one you just pulled, confirm you have
the right copy:

```r
grep("AUDITED", readLines("01_rebuild_backbone_AUDITED.R"), value = TRUE)[1]
```

Near-identical filenames caused two incidents during the audit, so the
check is worth the two seconds.

## Pinned inputs

None are in the repository. They are expected alongside the project
data directory.

| File | Written by | Consumed by |
|---|---|---|
| `cohort_inputs_FINAL.rds` | not yet scripted | 01, 18 |
| `backbone_FINAL.rds` | 01 | everything |
| `patient_level_FINAL.rds` | 02 | everything |
| `genotype_de_corrected_FINAL.rds` | 18 | 16, 21 |
| `proteomics_inputs_FINAL.rds` | 16a | 16, 21 |
| `tals_de_results.rds` | 14 | 16 |
| `genomic_burden_FINAL.rds` | `genomic_analysis.R` | burden scripts, 26 |

`cohort_inputs_FINAL.rds` has no build script. Writing
`00_build_inputs.R` to produce it from the raw release files, with an
`identical()` cross-check against the current file, is the one
outstanding gap in the provenance chain.

Useful crosswalks: recover `Participant_ID` with
`sub("-[0-9]+-T$", "", sample_id)`, and apply
`distinct(Participant_ID, .keep_all = TRUE)` to `cct2_df` before any
merge, because of longitudinal duplicates.

## Environment

R 4.3 or later with DESeq2, sva, limma, GSVA, ConsensusClusterPlus,
fgsea, msigdbr, glmnet, survival, org.Hs.eg.db, AnnotationDbi, Seurat,
data.table, dplyr, ggplot2, patchwork. The germline arm also needs
ANNOVAR, bcftools and PLINK. Figure 6 needs IPA exports.

Four things bite regardless of which script you are running.

- `ConsensusClusterPlus` calls `layout()` internally even with
  `plot = NULL` and fails headless. Wrap calls in `pdf(NULL)` and
  `dev.off()`. A failed `pdf()` leaves the device stack broken.
- `conflicted` is attached. Declare preferences before sourcing,
  `base::unname`, `base::setequal` and `BiocParallel::register` in
  particular.
- `fgsea` is Monte Carlo and is not reproducible from a seed alone
  when a parallel backend is registered. Set the seed immediately
  before each call and pass `nproc = 1` where exact reproduction
  matters.
- macOS filesystems are case-insensitive, so `backbone_final.rds` and
  `backbone_FINAL.rds` are the same file.

Target ALS data comes through the dnastack CLI.
`InvalidExistingConfigurationError` is fixed by deleting
`~/.dnastack/config.yaml` and re-running
`dnastack use dataportal.targetals.org`. Queries need
`gs_uri LIKE '%filename_suffix'` patterns, and `subject_id` is
unreliable for joins, so extract the Library ID from the filename
instead.

## Archive

Seven superseded scripts live in `archive/`. Do not run them. Nothing
here sources anything there.

| Script | Why |
|---|---|
| `11_spatial.R` | filename collision with `spatialdata.R` |
| `12_rna_protein_convergence.R` | pre-correction convergence layer |
| `13_convergence_addenda.R` | addenda to 12, same layer |
| `15_genotype_de_rebuild.R` | superseded by the fit-3 specification in 18 |
| `17_recover_inputs.R` | one-off workspace recovery after a crash |
| `18_full_recovery_and_de.R` | saves under the canonical filename **before** its checksum runs |
| `ipa_rerun_FINAL.R` | session dependency and an undefined `old_f`, cannot run clean |

The unsuffixed `18_full_recovery_and_de.R` is the one that caused real
damage. It writes the anchor-SV fit to
`genotype_de_corrected_FINAL.rds` before its checksum can object. To
check whether this has happened, look for the string "742-sample
anchor" in `x$design` and a `NULL` in `x$clean_svs_used`. The fix is to
rerun `18_full_recovery_and_de_AUDITED.R`.
