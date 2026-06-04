# ============================================================
# ALS project IO/bootstrap helpers
# ============================================================
# These helpers make legacy-style scripts safer to run from GitHub:
# - parse --config, --base-dir, --outdir, and --objects
# - load YAML config when available
# - load .rds object bundles into the script environment
# - redirect relative output filenames into --outdir
# - save selected objects as a handoff bundle

als_parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    config = "config/config.yml",
    base_dir = NULL,
    outdir = NULL,
    objects = character(0),
    redirect_outputs = TRUE,
    help = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      out$help <- TRUE; i <- i + 1
    } else if (key == "--config") {
      out$config <- args[[i + 1]]; i <- i + 2
    } else if (key == "--base-dir") {
      out$base_dir <- args[[i + 1]]; i <- i + 2
    } else if (key == "--outdir") {
      out$outdir <- args[[i + 1]]; i <- i + 2
    } else if (key == "--objects") {
      out$objects <- c(out$objects, unlist(strsplit(args[[i + 1]], ",", fixed = TRUE))); i <- i + 2
    } else if (key == "--no-output-redirect") {
      out$redirect_outputs <- FALSE; i <- i + 1
    } else {
      warning("Ignoring unknown argument: ", key)
      i <- i + 1
    }
  }
  out$objects <- unique(trimws(out$objects[out$objects != ""]))
  out
}

als_print_help <- function(script_id = "script") {
  cat("Usage:\n")
  cat("  Rscript ", script_id, " --config config/config.yml --base-dir /path/to/private/data --outdir results/<analysis> [--objects bundle1.rds,bundle2.rds]\n", sep = "")
  cat("\nOptions:\n")
  cat("  --config              YAML config file. Default: config/config.yml\n")
  cat("  --base-dir            Private ALS project/data root. Overrides project.base_dir in config.\n")
  cat("  --outdir              Output directory for figures/tables/object handoff.\n")
  cat("  --objects             Comma-separated .rds object bundles to load before analysis.\n")
  cat("  --no-output-redirect  Preserve legacy relative output locations.\n")
}

als_load_config <- function(path = "config/config.yml") {
  if (!file.exists(path)) return(list())
  if (!requireNamespace("yaml", quietly = TRUE)) {
    warning("Config file exists but R package 'yaml' is not installed. Returning empty config.")
    return(list())
  }
  yaml::read_yaml(path)
}

als_cfg_get <- function(cfg, key, default = NULL) {
  if (is.null(cfg) || length(cfg) == 0 || is.null(key) || key == "") return(default)
  parts <- strsplit(key, "\\.")[[1]]
  val <- cfg
  for (p in parts) {
    if (is.null(val[[p]])) return(default)
    val <- val[[p]]
  }
  val
}

als_choose <- function(arg_value, cfg, cfg_key, default) {
  if (!is.null(arg_value) && !is.na(arg_value) && nzchar(arg_value)) return(path.expand(arg_value))
  cfg_val <- als_cfg_get(cfg, cfg_key, NULL)
  if (!is.null(cfg_val) && !is.na(cfg_val) && nzchar(cfg_val)) return(path.expand(cfg_val))
  path.expand(default)
}

als_is_abs_or_home <- function(path) {
  grepl("^/|^~|^[A-Za-z]:", path)
}

als_out_path <- function(path) {
  if (!is.character(path) || length(path) == 0 || is.na(path[1])) return(path)
  if (!isTRUE(getOption("als_project.redirect_outputs", TRUE))) return(path)
  p <- path[1]
  if (!nzchar(p) || als_is_abs_or_home(p)) return(path)
  root <- getOption("als_project.outdir", "results")
  fp <- file.path(root, p)
  dir.create(dirname(fp), recursive = TRUE, showWarnings = FALSE)
  if (length(path) == 1) fp else c(fp, path[-1])
}

als_project_path <- function(path) {
  if (!is.character(path) || length(path) == 0 || is.na(path[1])) return(path)
  p <- path[1]
  if (als_is_abs_or_home(p)) return(path.expand(p))
  file.path(getOption("als_project.base_dir", getwd()), p)
}

als_init <- function(script_id, base_dir, outdir, redirect_outputs = TRUE) {
  if (!dir.exists(base_dir)) warning("Base directory does not exist yet: ", base_dir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  options(
    als_project.script_id = script_id,
    als_project.base_dir = normalizePath(base_dir, mustWork = FALSE),
    als_project.outdir = normalizePath(outdir, mustWork = FALSE),
    als_project.redirect_outputs = redirect_outputs
  )
  message("ALS project script: ", script_id)
  message("Base directory: ", getOption("als_project.base_dir"))
  message("Output directory: ", getOption("als_project.outdir"))
}

als_load_objects <- function(paths, envir = .GlobalEnv) {
  if (length(paths) == 0) return(invisible(character(0)))
  loaded <- character(0)
  for (p in paths) {
    p <- path.expand(p)
    if (!file.exists(p)) stop("Object bundle not found: ", p, call. = FALSE)
    obj <- readRDS(p)
    nm <- basename(p)
    if (is.list(obj) && !is.data.frame(obj) && !is.null(names(obj))) {
      for (name in names(obj)) assign(name, obj[[name]], envir = envir)
      loaded <- c(loaded, names(obj))
      message("Loaded object bundle: ", p, " [", length(names(obj)), " objects]")
    } else {
      assign(tools::file_path_sans_ext(nm), obj, envir = envir)
      loaded <- c(loaded, tools::file_path_sans_ext(nm))
      message("Loaded object: ", p)
    }
  }
  invisible(unique(loaded))
}

als_save_handoff <- function(script_id, object_names, envir = .GlobalEnv, outdir = getOption("als_project.outdir", "results")) {
  existing <- object_names[vapply(object_names, exists, logical(1), envir = envir, inherits = TRUE)]
  if (length(existing) == 0) {
    message("No configured handoff objects found for ", script_id, "; skipping .rds handoff.")
    return(invisible(NULL))
  }
  bundle <- stats::setNames(lapply(existing, get, envir = envir, inherits = TRUE), existing)
  f <- file.path(outdir, paste0(script_id, "_objects.rds"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, f)
  message("Saved handoff bundle: ", f, " [", length(existing), " objects]")
  invisible(f)
}

als_write_session_info <- function(script_id, outdir = getOption("als_project.outdir", "results")) {
  f <- file.path(outdir, paste0(script_id, "_sessionInfo.txt"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = f)
  invisible(f)
}

als_register_output_wrappers <- function(envir = .GlobalEnv) {
  # Preserve originals.
  .base_pdf <- grDevices::pdf
  .base_png <- grDevices::png
  .base_jpeg <- grDevices::jpeg
  .base_tiff <- grDevices::tiff
  .base_write_csv <- utils::write.csv
  .base_write_table <- utils::write.table
  .base_saveRDS <- base::saveRDS
  .base_writeLines <- base::writeLines
  assign("pdf", function(file, ...) .base_pdf(als_out_path(file), ...), envir = envir)
  assign("png", function(filename, ...) .base_png(als_out_path(filename), ...), envir = envir)
  assign("jpeg", function(filename, ...) .base_jpeg(als_out_path(filename), ...), envir = envir)
  assign("tiff", function(filename, ...) .base_tiff(als_out_path(filename), ...), envir = envir)
  assign("write.csv", function(x, file, ...) .base_write_csv(x, als_out_path(file), ...), envir = envir)
  assign("write.table", function(x, file, ...) .base_write_table(x, als_out_path(file), ...), envir = envir)
  assign("saveRDS", function(object, file = "", ...) .base_saveRDS(object, file = als_out_path(file), ...), envir = envir)
  assign("writeLines", function(text, con = stdout(), ...) {
    if (is.character(con) && length(con) == 1) con <- als_out_path(con)
    .base_writeLines(text, con = con, ...)
  }, envir = envir)
  assign("ggsave", function(filename, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for ggsave")
    ggplot2::ggsave(filename = als_out_path(filename), ...)
  }, envir = envir)
  assign("write_csv", function(x, file, ...) {
    if (!requireNamespace("readr", quietly = TRUE)) stop("readr is required for write_csv")
    readr::write_csv(x, als_out_path(file), ...)
  }, envir = envir)
  assign("write_tsv", function(x, file, ...) {
    if (!requireNamespace("readr", quietly = TRUE)) stop("readr is required for write_tsv")
    readr::write_tsv(x, als_out_path(file), ...)
  }, envir = envir)
  assign("fwrite", function(x, file, ...) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required for fwrite")
    data.table::fwrite(x, file = als_out_path(file), ...)
  }, envir = envir)
  invisible(TRUE)
}
