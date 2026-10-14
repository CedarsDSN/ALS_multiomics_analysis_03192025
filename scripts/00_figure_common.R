# =====================================================================
# 00_figure_common.R   (rewritten clean 2026-08-26; updated 2026-10-05)
#
# Shared loader, palette, theme and helpers for the figure scripts.
#
# THE ONE THING THIS FILE MUST GET RIGHT
# ~/Documents/convergence_results/ holds the CONTAMINATED script-12
# outputs. Its convergence_FINAL.rds gives RNA pathway NES
# 0.761 / -0.077 / 0.104, computed on the pre-correction DE layer whose
# control group contains the replicated reference line. The corrected
# outputs live in convergence_results_corrected/.
#
# That directory is BLOCKED by absolute path. The block is applied after
# EVERY lookup, direct or recursive, so no search route can reach it.
# The single exception is genotype_de_corrected_FINAL.rds, which is
# current and happens to live there; it is returned by explicit name.
#
# 2026-10-05 changes (search "NEW" / "CHANGED"):
#   - save_fig() now also writes a 300 dpi TIFF next to every PDF and
#     checks that the PDF fonts are embedded.
#   - new check_fonts() and check_all_figures() for a pre-submission pass.
#   - new fmt_p() so every panel formats P values the same way.
#   - Fig 8 now depends on downstream_analyses_FINAL.rds (panel c uses the
#     canonical SV-adjusted models); readiness check updated to match.
#   - fixed the doubled "save_fig <- save_fig <-" assignment.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})
set.seed(42)

## ---- blocked paths ---------------------------------------------------
BLOCKED_DIRS <- path.expand("~/Documents/convergence_results")
ALLOW_ANYWAY <- "genotype_de_corrected_FINAL.rds"

is_blocked <- function(paths) {
  if (!length(paths)) return(logical(0))
  norm    <- normalizePath(paths, mustWork = FALSE)
  blocked <- normalizePath(BLOCKED_DIRS, mustWork = FALSE)
  vapply(norm, function(p)
    any(startsWith(p, paste0(blocked, .Platform$file.sep))),
    logical(1), USE.NAMES = FALSE)
}

## ---- search path -----------------------------------------------------
PIN_DIRS <- getOption("als.pin_dirs", c(
  ".",
  "~/Documents/DSN/ALS",
  "~/Documents/convergence_results_corrected",
  "~/Documents/DSN/ALS/convergence_results_corrected",
  "~/Documents/DSN/ALS/clinical_endpoints",
  "~/Documents/DSN/ALS/proteomics",
  "~/Documents/DSN/ALS/prset",
  "~/Documents/DSN/ALS/cct2_eqtl",
  "~/Documents/DSN/ALS/tric_burden_results",
  "~/Documents/DSN/ALS/myc_burden_results",
  "~/Documents/DSN/ALS/IPA_inputs_FINAL",
  "~/Documents/DSN/ALS/targetALS",
  "~/Documents/targetALS",
  "~/Documents/IPA_inputs_FINAL",
  "~/Documents"
))
PIN_DIRS <- path.expand(PIN_DIRS)
PIN <- PIN_DIRS[1]

## ---- resolve a file --------------------------------------------------
p_ <- function(f) {
  base <- basename(f)
  if (base %in% ALLOW_ANYWAY) {
    ex <- file.path(BLOCKED_DIRS, base)
    ex <- ex[file.exists(ex)]
    if (length(ex)) return(ex[1])
  }
  direct <- file.path(PIN_DIRS, f)
  direct <- direct[file.exists(direct)]
  if (length(direct)) direct <- direct[!is_blocked(direct)]
  if (length(direct)) return(direct[1])
  for (d in PIN_DIRS) {
    if (!dir.exists(d)) next
    g <- list.files(d, pattern = paste0("^", base, "$"),
                    recursive = TRUE, full.names = TRUE)
    if (length(g)) g <- g[!is_blocked(g)]
    if (length(g)) return(g[1])
  }
  file.path(PIN_DIRS[1], f)
}

need_file <- function(f) {
  path <- p_(f)
  if (!file.exists(path))
    stop("required pinned file not found: ", basename(f),
         "\n  searched: ", paste(PIN_DIRS, collapse = ", "))
  path
}

load_pin <- function(f) {
  path <- need_file(f)
  cat("  loaded", basename(path), "from", dirname(path), "\n")
  readRDS(path)
}

locate_inputs <- function() {
  want <- ALL_INPUTS
  cat("blocked:", BLOCKED_DIRS, "\n\n")
  for (f in want) {
    path <- p_(f)
    cat(sprintf("  %-40s %s\n", f,
                if (file.exists(path)) dirname(path) else "NOT FOUND"))
  }
  cat("\nIPA exports (Figure 6 needs the noSexChr pair):\n")
  ipa <- unlist(lapply(PIN_DIRS[dir.exists(PIN_DIRS)], function(d)
    list.files(d, pattern = "Upstream|Canonical", ignore.case = TRUE,
               full.names = TRUE, recursive = TRUE)))
  ipa <- unique(ipa[!is_blocked(ipa)])
  if (length(ipa)) for (i in ipa) cat("  ", i, "\n") else
    cat("   NONE FOUND -- Figure 6 will stop\n")
  invisible(NULL)
}

## Every input any figure script reads, and which figure needs it.
ALL_INPUTS <- c(
  "backbone_FINAL.rds",                 # fig 1,3,4,5
  "patient_level_FINAL.rds",            # fig 1,3,5,8
  "genotype_de_corrected_FINAL.rds",    # fig 2
  "convergence_FINAL.rds",              # fig 2 (RNA NES row)
  "proteomics_inputs_FINAL.rds",        # fig 2 (protein NES row)
  "proteomics_cct2_FINAL.rds",          # fig 5 e/f
  "cct2_protein_by_subtype.csv",        # fig 5 e/f (csv route)
  "lasso_FINAL.rds",                    # fig 5 (stability, if used)
  "genomic_burden_FINAL.rds",           # fig 7
  "cct2_cis_eqtl_results.csv",          # fig 7
  "tric_burden_full.csv",               # fig 7
  "clinical_screen_combined_FINAL.rds", # fig 8 a/b
  "downstream_analyses_FINAL.rds",      # fig 8 c (CHANGED: now required; SV + composition models)
  "spatial_results_FINAL.rds"           # fig 9
)

check_convergence_source <- function() {
  path <- p_("convergence_FINAL.rds")
  cat("convergence_FINAL.rds resolves to:\n  ", path, "\n")
  if (isTRUE(is_blocked(path)))
    stop("resolved to the CONTAMINATED convergence file. This must not happen.")
  cv <- readRDS(path)
  if (!is.null(cv$pattern_table)) {
    pt <- as.data.frame(cv$pattern_table)
    cat("RNA rows in this file:\n")
    print(pt[pt$omic == "RNA", c("pair", "n", "spearman")])
    cat("(contaminated signature 0.761 / -0.077 / 0.104;",
        "corrected ~0.40 / -0.06 / 0.05)\n")
  }
  invisible(path)
}

## ---- palettes and theme ----------------------------------------------
## Subtype and genotype palettes are deliberately disjoint, so a reader
## never links S-ALS-3 with SOD1 or S-ALS-1 with C9orf72.
SCOL <- c("S-ALS-1" = "#7E57C2", "S-ALS-2" = "#90A4AE", "S-ALS-3" = "#F28E2B")
GCOL <- c(Control = "#7F7F7F", Sporadic = "#4DBBD5",
          C9orf72 = "#E64B35", SOD1 = "#00A087")
## NEW: neutral colours for non-subtype, non-genotype categories
## (case vs control in Fig 9, highlighted points in Fig 8, etc.)
NCOL <- c(neutral = "grey20", muted = "grey60", case = "#5B6C8F", control = "#7F7F7F")
AX   <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")
ALAB <- c(AXIS_METABOLIC = "Metabolic", AXIS_INFLAMMATION = "Inflammation",
          AXIS_CELL_CYCLE = "Cell cycle")
W    <- 180 / 25.4     # Nature Communications double-column width, inches

thm <- theme_classic(base_size = 8) +
  theme(axis.text     = element_text(size = 7, colour = "black"),
        axis.title    = element_text(size = 8),
        legend.text   = element_text(size = 7),
        legend.title  = element_text(size = 8),
        legend.key.size = unit(3, "mm"),
        plot.title    = element_text(size = 9, hjust = 0),
        plot.subtitle = element_text(size = 7),
        strip.background = element_blank(),
        strip.text    = element_text(size = 7),
        plot.margin   = margin(3, 3, 3, 3))

void_thm <- theme_void(base_size = 8) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 9, hjust = 0),
        plot.subtitle = element_text(size = 7),
        plot.margin   = margin(3, 3, 3, 3))

## ---- standard accessors ----------------------------------------------
get_backbone <- function() {
  b <- load_pin("backbone_FINAL.rds")
  stopifnot(!is.null(b$ax_round1), !is.null(b$meta))
  d <- as.data.frame(b$ax_round1)
  if (is.null(d$Participant_ID))
    d$Participant_ID <- b$meta$Participant_ID[match(d$sample_id, b$meta$sample_id)]
  if (is.null(d$subtype)) d$subtype <- d$cluster_r1
  d$subtype <- factor(as.character(d$subtype), levels = names(SCOL))
  b$ax_round1 <- d
  b
}

get_patients <- function() {
  pl  <- load_pin("patient_level_FINAL.rds")
  pat <- as.data.frame(pl$pat)
  pat$subtype <- factor(as.character(pat$subtype), levels = names(SCOL))
  stopifnot(!any(duplicated(pat$Participant_ID)))
  attr(pat, "pl") <- pl
  pat
}

## NEW: the patient table with SVs and composition scores (Fig 8 c, and any
## future panel that needs the canonical covariates). Same 391 patients.
get_patients_full <- function() {
  ds <- load_pin("downstream_analyses_FINAL.rds")
  d  <- as.data.frame(ds$df)
  d$subtype <- factor(as.character(d$subtype), levels = names(SCOL))
  stopifnot(nrow(d) == 391, !any(duplicated(d$Participant_ID)))
  attr(d, "covariates") <- ds$covariates   # "AGE_AT_SYMPTOM_ONSET + SEX + SV1 + ... + SV6"
  d
}

## ---- genotype group -----------------------------------------------------
## meta carries SUBJECT_GROUP.x / SUBJECT_SUBGROUP rather than a tidy
## "group" column, and may or may not carry Mutation_Subtype. Resolve it
## once here so the figures do not each guess.
get_group <- function(meta) {
  meta <- as.data.frame(meta)
  direct <- intersect(c("Mutation_Subtype", "group", "Group",
                        "genotype", "Genotype"), names(meta))
  if (length(direct)) {
    g <- as.character(meta[[direct[1]]])
    cat("  group column:", direct[1], "\n")
  } else {
    sg  <- grep("^SUBJECT_GROUP", names(meta), value = TRUE)[1]
    ssg <- grep("^SUBJECT_SUBGROUP", names(meta), value = TRUE)[1]
    if (is.na(sg))
      stop("no genotype column found in meta. Columns: ",
           paste(names(meta), collapse = ", "))
    cat("  deriving group from", sg,
        if (!is.na(ssg)) paste("+", ssg) else "", "\n")
    grp <- as.character(meta[[sg]])
    sub <- if (!is.na(ssg)) as.character(meta[[ssg]]) else rep("", nrow(meta))
    g <- ifelse(grepl("control", grp, ignore.case = TRUE), "Control",
                ifelse(grepl("c9", sub, ignore.case = TRUE), "C9orf72",
                       ifelse(grepl("sod1", sub, ignore.case = TRUE), "SOD1",
                              ifelse(grepl("als|mnd", grp, ignore.case = TRUE), "Sporadic",
                                     NA_character_))))
  }
  g <- factor(g, levels = c("Control", "Sporadic", "C9orf72", "SOD1"))
  cat("  group counts:", paste(names(table(g)), table(g), sep = "=",
                               collapse = "  "),
      "| unmapped:", sum(is.na(g)), "\n")
  cat("  (expect Control 185 / Sporadic 400 / C9orf72 46 / SOD1 13 at",
      "sample level)\n")
  g
}

## ---- helpers ----------------------------------------------------------
auc_fn <- function(score, pos) {
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

lfc_rho <- function(a, b) {
  sh <- intersect(rownames(a), rownames(b))
  suppressWarnings(cor(a[sh, "log2FoldChange"], b[sh, "log2FoldChange"],
                       method = "spearman", use = "complete.obs"))
}

fresh_device <- function() {
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
}

## NEW: one P-value format for every panel, matching the manuscript text:
##   P >= 0.001 -> "P = 0.033" (2 significant digits, 3 decimals max)
##   P <  0.001 -> "P = 5.2e-4" style, or "P < 0.001" if short = TRUE
## Capital P throughout, as in the text.
fmt_p <- function(p, short = FALSE, prefix = "P = ") {
  vapply(p, function(x) {
    if (is.na(x)) return("P = NA")
    if (x >= 0.001) return(paste0(prefix, formatC(signif(x, 2), format = "fg", digits = 2)))
    if (short) return("P < 0.001")
    paste0(prefix, formatC(x, format = "e", digits = 1))
  }, character(1))
}

## NEW: confirm a PDF has every font embedded (Nature requirement).
## Uses the poppler `pdffonts` tool if it is installed (brew install poppler);
## if not, it says so and moves on rather than failing.
check_fonts <- function(file, quiet = FALSE) {
  if (!nzchar(Sys.which("pdffonts"))) {
    if (!quiet) cat("  (pdffonts not found; install poppler to check font embedding)\n")
    return(invisible(NA))
  }
  out <- system2("pdffonts", file, stdout = TRUE, stderr = TRUE)
  rows <- out[-(1:2)]
  if (!length(rows)) { if (!quiet) cat("  ", basename(file), ": no fonts (pure vector)\n"); return(invisible(TRUE)) }
  ## 'emb' is the 5th field from the right in pdffonts output
  emb <- vapply(strsplit(trimws(rows), "\\s+"),
                function(f) f[length(f) - 4], character(1))
  ok <- all(emb == "yes")
  if (!quiet) cat(sprintf("  %-32s fonts embedded: %s\n", basename(file),
                          if (ok) "yes" else paste("NO  <--", paste(rows[emb != "yes"], collapse = "; "))))
  invisible(ok)
}

## CHANGED: save_fig() writes the vector PDF (for the manuscript) and, by
## default, a 300 dpi TIFF (for the journal / collaborators) from the same
## plot object, then checks the PDF fonts. Both use the native macOS quartz
## device, which embeds fonts without XQuartz.
##   tiff = FALSE   skip the TIFF
##   dpi  = 600     higher resolution if a journal asks for it
save_fig <- function(plot, file, height, width = W, tiff = TRUE, dpi = 300) {
  grDevices::quartz(file = file, type = "pdf", width = width, height = height)
  print(plot)
  grDevices::dev.off()
  cat("wrote", file, sprintf("(%.0f x %.0f mm)\n", width * 25.4, height * 25.4))
  check_fonts(file)
  if (tiff) {
    tf <- sub("\\.pdf$", ".tiff", file, ignore.case = TRUE)
    if (identical(tf, file)) tf <- paste0(file, ".tiff")
    grDevices::quartz(file = tf, type = "tiff", width = width, height = height,
                      dpi = dpi, bg = "white")
    print(plot)
    grDevices::dev.off()
    cat("wrote", tf, sprintf("(%d dpi)\n", dpi))
  }
  invisible(file)
}

## NEW: pre-submission pass over every figure PDF in the working directory.
check_all_figures <- function(pattern = "^Figure.*\\.pdf$") {
  fs <- list.files(".", pattern = pattern)
  if (!length(fs)) { cat("no figure PDFs found matching", pattern, "\n"); return(invisible(NULL)) }
  cat("Font embedding:\n")
  res <- vapply(fs, check_fonts, logical(1))
  tifs <- sub("\\.pdf$", ".tiff", fs)
  cat("\n300 dpi TIFF present for:",
      paste(fs[file.exists(tifs)], collapse = ", "), "\n")
  miss <- fs[!file.exists(tifs)]
  if (length(miss)) cat("TIFF missing for:", paste(miss, collapse = ", "), "\n")
  invisible(res)
}

## ---- structure report --------------------------------------------------
report_structure <- function() {
  files <- grep("\\.rds$", ALL_INPUTS, value = TRUE)
  for (f in files) {
    path <- p_(f)
    cat("\n==================================================\n", f, "\n")
    if (!file.exists(path)) { cat("  ABSENT\n"); next }
    cat("  dir:", dirname(path), "\n")
    x <- readRDS(path)
    cat("  class:", paste(class(x), collapse = "/"), "\n")
    if (is.list(x) && !is.data.frame(x)) {
      for (nm in names(x)) {
        el <- x[[nm]]
        cat(sprintf("  $%-22s %-16s", nm, paste(class(el), collapse = "/")))
        if (is.data.frame(el) || is.matrix(el))
          cat(" dim ", paste(dim(el), collapse = " x "))
        if (is.data.frame(el))
          cat("\n      cols: ", paste(utils::head(names(el), 25), collapse = ", "))
        cat("\n")
      }
    } else if (is.data.frame(x)) {
      cat("  dim ", paste(dim(x), collapse = " x "), "\n  cols: ",
          paste(names(x), collapse = ", "), "\n")
    }
  }
  invisible(NULL)
}

## ---- readiness: can each figure run right now? -----------------------
figure_readiness <- function() {
  need <- list(
    "Fig 1 cohort"      = c("backbone_FINAL.rds", "patient_level_FINAL.rds"),
    "Fig 2 convergence" = c("genotype_de_corrected_FINAL.rds",
                            "convergence_FINAL.rds", "proteomics_inputs_FINAL.rds"),
    "Fig 3 axes"        = c("backbone_FINAL.rds", "patient_level_FINAL.rds"),
    "Fig 4 consensus"   = "backbone_FINAL.rds",
    "Fig 5 CCT2"        = c("patient_level_FINAL.rds", "cct2_protein_by_subtype.csv"),
    "Fig 6 mechanism"   = character(0),          # IPA csv, checked separately
    "Fig 7 genomic"     = c("genomic_burden_FINAL.rds", "cct2_cis_eqtl_results.csv",
                            "tric_burden_full.csv"),
    "Fig 8 clinical"    = c("clinical_screen_combined_FINAL.rds",
                            "downstream_analyses_FINAL.rds"),   # CHANGED
    "Fig 9 spatial"     = "spatial_results_FINAL.rds")
  for (fg in names(need)) {
    miss <- need[[fg]][!file.exists(vapply(need[[fg]], p_, ""))]
    if (fg == "Fig 6 mechanism") {
      ipa <- unlist(lapply(PIN_DIRS[dir.exists(PIN_DIRS)], function(d)
        list.files(d, pattern = "Upstream|Canonical", ignore.case = TRUE,
                   full.names = TRUE, recursive = TRUE)))
      miss <- if (length(ipa[!is_blocked(ipa)]) >= 2) character(0) else "IPA exports"
    }
    cat(sprintf("  %-20s %s\n", fg,
                if (!length(miss)) "READY" else
                  paste("missing:", paste(miss, collapse = ", "))))
  }
  invisible(NULL)
}

cat("figure common loaded.\n")
cat("  locate_inputs()            where each pinned file resolves to\n")
cat("  check_convergence_source() confirm the corrected convergence file\n")
cat("  report_structure()         inspect the pinned objects\n")
cat("  figure_readiness()         which of the nine figures can run now\n")
cat("  check_all_figures()        fonts embedded + 300 dpi TIFF present, all figures\n")