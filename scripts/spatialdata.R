library(Seurat)
library(dplyr)
library(stringr)
library(readr)
library(Matrix)
library(hdf5r)
library(DESeq2)

mainDir <- "~/Documents/DSN/ALS/targetALS/"
setwd(mainDir)
download_root  <- "./spatial_data"     # wherever batch_download.py wrote to
processed_root <- "./processed"
meta <- read_csv("collections.postmortem_tissue_core.columbia_rnaseq_spatial_metadata.csv")

all_files <- list.files(download_root, recursive = TRUE, full.names = TRUE)
cat("Total files downloaded:", length(all_files), "\n")
print(head(all_files, 10))

meta <- meta %>%
  mutate(filename_key = `Library ID`)

dir.create(processed_root, showWarnings = FALSE)

for (i in seq_len(nrow(meta))) {
  key <- meta$filename_key[i]
  lib <- meta$`Library ID`[i]
  area_dir <- file.path(processed_root, lib)
  spatial_dir <- file.path(area_dir, "spatial")
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Match on the library ID bounded by dots (matches the confirmed filename
  # pattern "{Subject}.{LibraryID}.{TissueTag}.{filetype}"), avoiding any
  # partial-substring collisions between similarly-prefixed library IDs.
  matches <- all_files[str_detect(all_files, fixed(paste0(".", key, ".")))]
  
  h5  <- matches[str_detect(matches, "filtered_feature_bc_matrix\\.h5$")]
  pos <- matches[str_detect(matches, "tissue_positions\\.csv$")]
  scf <- matches[str_detect(matches, "scalefactors_json\\.json$")]
  img <- matches[str_detect(matches, "tissue_hires_image\\.png$")]
  
  if (length(h5)  == 1) file.copy(h5,  file.path(area_dir, "filtered_feature_bc_matrix.h5"))
  if (length(pos) == 1) file.copy(pos, file.path(spatial_dir, "tissue_positions.csv"))
  if (length(scf) == 1) file.copy(scf, file.path(spatial_dir, "scalefactors_json.json"))
  if (length(img) == 1) file.copy(img, file.path(spatial_dir, "tissue_hires_image.png"))
  
  if (length(h5) != 1 || length(pos) != 1 || length(scf) != 1) {
    warning(paste("Incomplete file set for", lib, "-- check manually"))
  }
}

layer5b_markers <- c("FEZF2", "CRYM")
betz_markers    <- c("BCL11B", "POU3F1", "NEFH")

pseudobulk_one_area <- function(library_id, q = 0.75) {
  area_dir <- file.path(processed_root, library_id)
  
  # Load10X_Spatial defaults to tissue_lowres_image.png, which we didn't
  # download (only hires, to keep the pull small). Build the image object
  # explicitly pointing at the hires file instead.
  img_obj <- Read10X_Image(image.dir = file.path(area_dir, "spatial"),
                           image.name = "tissue_hires_image.png",
                           filter.matrix = TRUE)
  obj <- Load10X_Spatial(data.dir = area_dir, image = img_obj)
  
  # AddModuleScore needs log-normalized data, not raw counts. Without this,
  # the "data" layer doesn't exist yet and AddModuleScore fails.
  obj <- NormalizeData(obj)
  
  obj <- AddModuleScore(obj, features = list(c(layer5b_markers, betz_markers)),
                        name = "Layer5b_Betz_Score")
  thresh <- quantile(obj$Layer5b_Betz_Score1, q)
  keep_spots <- colnames(obj)[obj$Layer5b_Betz_Score1 >= thresh]
  
  # Pseudobulking still uses raw counts (correct -- normalize once, at the
  # pseudobulk level, via DESeq2 in Step 4), not the log-normalized layer.
  counts <- GetAssayData(obj, layer = "counts")[, keep_spots, drop = FALSE]
  Matrix::rowSums(counts)
}

pb_list <- list()
for (lib in meta$`Library ID`) {
  message("Processing ", lib)
  pb_list[[lib]] <- tryCatch(
    pseudobulk_one_area(lib),
    error = function(e) { message("  failed: ", e$message); NULL }
  )
}
pb_list <- pb_list[!sapply(pb_list, is.null)]


pb_matrix <- do.call(cbind, pb_list)
colnames(pb_matrix) <- names(pb_list)

lib_to_subject <- setNames(meta$`Subject ID`, meta$`Library ID`)
lib_to_tissue  <- setNames(meta$`Tissue Source`, meta$`Library ID`)
subject_tissue <- paste(lib_to_subject[colnames(pb_matrix)],
                        lib_to_tissue[colnames(pb_matrix)], sep = "__")

pb_by_subject <- sapply(unique(subject_tissue), function(st) {
  cols <- colnames(pb_matrix)[subject_tissue == st]
  if (length(cols) == 1) pb_matrix[, cols]
  else Matrix::rowSums(pb_matrix[, cols, drop = FALSE])
})
colnames(pb_by_subject) <- unique(subject_tissue)

subject_ids <- sapply(strsplit(colnames(pb_by_subject), "__"), `[`, 1)
tissue_ids  <- sapply(strsplit(colnames(pb_by_subject), "__"), `[`, 2)
clinical_group <- meta$`Clinical Group`[match(subject_ids, meta$`Subject ID`)]

col_data <- data.frame(row.names = colnames(pb_by_subject),
                       subject = subject_ids, tissue = tissue_ids,
                       clinical_group = clinical_group)

dds <- DESeqDataSetFromMatrix(countData = round(pb_by_subject),
                              colData = col_data, design = ~ 1)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

library(msigdbr)
library(GSVA)

hallmark_tbl <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()
hallmark_list <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)

axis_sets <- list(
  AXIS_METABOLIC = unique(c(
    hallmark_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
    hallmark_list[["HALLMARK_FATTY_ACID_METABOLISM"]],
    hallmark_list[["HALLMARK_GLYCOLYSIS"]]
  )),
  AXIS_CELL_CYCLE = unique(c(
    hallmark_list[["HALLMARK_E2F_TARGETS"]],
    hallmark_list[["HALLMARK_G2M_CHECKPOINT"]],
    hallmark_list[["HALLMARK_MYC_TARGETS_V1"]],
    hallmark_list[["HALLMARK_DNA_REPAIR"]]
  )),
  AXIS_INFLAMMATION = unique(c(
    hallmark_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
    hallmark_list[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
    hallmark_list[["HALLMARK_COMPLEMENT"]]
  ))
)

# Check rownames(vst_mat) first -- if it's already gene symbols (likely,
# given Space Ranger's standard 10x reference), skip straight to scoring.
# Only run the ENSEMBL-to-SYMBOL mapping block from the bulk script if
# rownames(vst_mat) come back as ENSG IDs instead.

vst_mat_num <- as.matrix(vst_mat)
storage.mode(vst_mat_num) <- "double"

param <- ssgseaParam(vst_mat_num, axis_sets, normalize = TRUE)
axis_scores <- gsva(param)

axis_df <- as.data.frame(t(axis_scores))   # samples x axes
axis_df <- as.data.frame(scale(axis_df))   # z-score each axis, matching
# the bulk pipeline convention
axis_df$sample_id <- rownames(axis_df)

cct2 <- vst_mat["CCT2", ]
# If this errors, check whether rownames are gene symbols or ENSEMBL IDs;
# if ENSEMBL, strip the version suffix first with your usual
# sub("\\..*$", "", ...) convention before matching.

results_df <- col_data
results_df$sample_id <- rownames(col_data)
results_df$CCT2 <- cct2
results_df <- merge(results_df, axis_df, by = "sample_id")

# Frontal cortex: the real test (n=28 ALS vs 11 control)
frontal <- results_df %>% filter(tissue == "Frontal Cortex")
wilcox.test(CCT2 ~ clinical_group, data = frontal)
wilcox.test(AXIS_CELL_CYCLE ~ clinical_group, data = frontal)
wilcox.test(AXIS_METABOLIC ~ clinical_group, data = frontal)

# Motor cortex: descriptive only (n=2 vs 2). Report values and n, never a
# p-value framed as a finding.
motor <- results_df %>% filter(str_detect(tissue, "Motor"))
print(motor[, c("subject", "tissue", "clinical_group", "CCT2",
                "AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")])

tapply(frontal$CCT2, frontal$clinical_group, median)


frontal_als <- frontal %>% filter(clinical_group == "ALS")
cat("within-ALS frontal n =", nrow(frontal_als), "\n")

## The validation-relevant readout: does the CCT2-metabolic coupling
## that defines S-ALS-3 exist in vivo, across ALS subjects?
print(cor.test(frontal_als$CCT2, frontal_als$AXIS_METABOLIC,    method = "spearman"))
print(cor.test(frontal_als$CCT2, frontal_als$AXIS_CELL_CYCLE,   method = "spearman"))
print(cor.test(frontal_als$CCT2, frontal_als$AXIS_INFLAMMATION, method = "spearman"))

## Dispersion: a high-CCT2 subgroup predicts wider spread in ALS
print(fligner.test(CCT2 ~ clinical_group, data = frontal))
print(tapply(frontal$CCT2, frontal$clinical_group,
             function(x) c(median = median(x), IQR = IQR(x), n = length(x))))

## For completeness alongside the other two axes
print(wilcox.test(AXIS_INFLAMMATION ~ clinical_group, data = frontal))
print(tapply(frontal$AXIS_METABOLIC, frontal$clinical_group, median))

cc_set_noCCT2 <- setdiff(axis_sets$AXIS_CELL_CYCLE, "CCT2")
param_cc <- ssgseaParam(vst_mat_num, list(AXIS_CC_noCCT2 = cc_set_noCCT2),
                        normalize = TRUE)
cc_no <- as.numeric(gsva(param_cc))
names(cc_no) <- colnames(vst_mat_num)

frontal_als$AXIS_CC_noCCT2 <- scale(cc_no)[frontal_als$sample_id, 1]
print(cor.test(frontal_als$CCT2, frontal_als$AXIS_CC_noCCT2, method = "spearman"))

saveRDS(list(results_df = results_df,
             frontal_tests = list(
               cct2_cc      = wilcox.test(CCT2 ~ clinical_group, data = frontal),
               inflam_cc    = wilcox.test(AXIS_INFLAMMATION ~ clinical_group, data = frontal),
               cct2_ccaxis  = cor.test(frontal_als$CCT2, frontal_als$AXIS_CELL_CYCLE, method = "spearman"),
               cct2_ccaxis_noCCT2 = cor.test(frontal_als$CCT2, frontal_als$AXIS_CC_noCCT2, method = "spearman"),
               cct2_metab   = cor.test(frontal_als$CCT2, frontal_als$AXIS_METABOLIC, method = "spearman"),
               cct2_inflam  = cor.test(frontal_als$CCT2, frontal_als$AXIS_INFLAMMATION, method = "spearman"),
               fligner      = fligner.test(CCT2 ~ clinical_group, data = frontal)),
             motor_descriptive = motor,
             axis_cc_noCCT2 = cc_no),
        "spatial_results_FINAL.rds")
