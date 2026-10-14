# 08_cct2_protein_test.R   (AUDITED 2026-08-22)
# CCT2 protein by subtype -- codifies the verified console analysis.
# This test previously existed ONLY as a chat/console block; this file
# is its durable form.
#
# Requires the collapsed proteomics matrix from your proteomics script
# (prot_expr_collapsed: proteins x 436 participants, columns named by
# Participant_ID). Point PROT_RDS at wherever that script saves it, or
# run this in the proteomics session.
#
# Verified landmarks (2026-08-21): n = 265 (59/144/62);
# S-ALS-3 vs S-ALS-2 Wilcoxon W = 3629, p = 0.03347;
# Kruskal-Wallis p = 0.004838; CCT2 RNA-protein rho = 0.2051535,
# p = 7.985e-4; TRiC composite S3 vs S2 p = 0.1324 (ns, by design the
# specificity check). CCT2 is quantified as "1/sp|P78371|TCPB_HUMAN"
# and must be rescued manually: TCPB_HUMAN defeats symbol mapping.

suppressPackageStartupMessages({ library(dplyr) })
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               base::setdiff, base::intersect, base::union,
                               base::unname, .quiet = TRUE)
}
set.seed(42)
setwd("~/Documents/DSN/ALS")

## ---- 0. Inputs ----
if (!exists("prot_expr_collapsed")) {
  PROT_RDS <- "proteomics/proteomics_matrix_collapsed.rds"   # EDIT if needed
  if (file.exists(PROT_RDS)) {
    prot_expr_collapsed <- readRDS(PROT_RDS)
  } else stop("prot_expr_collapsed not in session and ", PROT_RDS,
              " not found. Run the proteomics script first or fix the path.")
}
if (!exists("prot_gene"))
  stop("prot_gene (symbol-mapped matrix) not in session; run the ",
       "proteomics script sections 1-5 first.")

pl  <- readRDS("patient_level_FINAL.rds")
pat <- if (is.data.frame(pl)) pl else pl$pat
sub_col <- intersect(c("subtype", "cluster_final", "cluster"), names(pat))[1]
rna_col <- intersect(c("CCT2", "CCT2_vst"), names(pat))[1]
stopifnot(!is.na(sub_col))

## ---- 1. Rescue CCT2 (UniProt P78371; symbol mapping fails on TCPB) ----
cct2_row <- grep("P78371", rownames(prot_expr_collapsed), value = TRUE)[1]
stopifnot(!is.na(cct2_row))
cct2_vec <- as.numeric(prot_expr_collapsed[cct2_row, ])
names(cct2_vec) <- colnames(prot_expr_collapsed)
cat("CCT2 rescued from row:", cct2_row, "|",
    sum(is.finite(cct2_vec)), "participants\n")

## ---- 2. Test frame ----
df <- tibble::tibble(Participant_ID = names(cct2_vec),
                     CCT2_prot = cct2_vec) %>%
  filter(is.finite(CCT2_prot)) %>%
  inner_join(pat %>% distinct(Participant_ID, .keep_all = TRUE) %>%
               select(Participant_ID, subtype = all_of(sub_col),
                      CCT2_rna = any_of(rna_col)),
             by = "Participant_ID")
cat("n =", nrow(df), "(expect 265)\n")
print(table(df$subtype))          # expect 59 / 144 / 62

## ---- 3. Tests ----
d23 <- df %>% filter(grepl("2$|3$", as.character(subtype)))
w  <- wilcox.test(CCT2_prot ~ droplevels(factor(subtype)), data = d23)
kw <- kruskal.test(CCT2_prot ~ subtype, data = df)
print(w); print(kw)
print(df %>% group_by(subtype) %>%
        summarise(n = n(), median = median(CCT2_prot), IQR = IQR(CCT2_prot)))
ct <- if (!is.na(rna_col))
  cor.test(df$CCT2_rna, df$CCT2_prot, method = "spearman") else NULL
if (!is.null(ct)) print(ct)

## ---- 4. TRiC composite (specificity check; expected ns) ----
tric_syms <- c("TCP1", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8")
stopifnot(all(tric_syms %in% rownames(prot_gene)))
common_ids <- intersect(colnames(prot_gene), names(cct2_vec))
tric_mat <- rbind(prot_gene[tric_syms, common_ids, drop = FALSE],
                  CCT2 = cct2_vec[common_ids])
tric_score <- colMeans(t(scale(t(tric_mat))), na.rm = TRUE)
df$TRiC_score  <- tric_score[df$Participant_ID]
d23$TRiC_score <- tric_score[d23$Participant_ID]
wt <- wilcox.test(TRiC_score ~ droplevels(factor(subtype)),
                  data = d23 %>% filter(is.finite(TRiC_score)))
print(wt)

## ---- 5. Verification gate vs the saved analysis ----
out_csv <- "proteomics/cct2_protein_by_subtype.csv"
gate <- c(n_ok = nrow(df) == 265,
          w_ok = abs(w$p.value - 0.03347) < 1e-3,
          kw_ok = abs(kw$p.value - 0.004838) < 1e-4,
          rho_ok = is.null(ct) || abs(unname(ct$estimate) - 0.2051535) < 1e-4,
          tric_ok = abs(wt$p.value - 0.1324) < 1e-3)
print(gate)
if (all(gate)) {
  readr::write_csv(df, out_csv)
  cat("EXACT MATCH: reproduces the verified analysis. Script VERIFIED;",
      "saved", out_csv, "\n")
} else {
  readr::write_csv(df, sub("\\.csv$", "_AUDITCHECK.csv", out_csv))
  cat("MISMATCH: saved _AUDITCHECK.csv; existing csv untouched.\n")
}
sessionInfo()
