# -----------------------------------------------------------------------------
# Script: mmrm_subgroup_analysis.R
# Goal:
#   Re-run the MMRM pipeline from scripts/mmrm_main_analysis.R separately in
#   baseline-defined subgroups:
#     1) visit == bl and edeqpurging > 0
#     2) visit == bl and edeq15obe > 0
#     3) visit == bl and edeqpurging = 0 and edeq15obe > 0 (binge only)
#
# Approach:
#   - Read the same long-format analysis dataset used by mmrm_main_analysis.R.
#   - Create subject-level subgroup membership using baseline values.
#   - For each subgroup, write a temporary filtered dataset (all visits kept for
#     eligible participants) and execute the main MMRM pipeline with subgroup-
#     specific input/output paths.
# -----------------------------------------------------------------------------

source_script <- file.path("scripts", "mmrm_main_analysis.R")
input_path <- file.path("data", "analysis_dataset.csv")
output_root <- file.path("outputs", "mmrm_subgroup_analysis")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(source_script))
stopifnot(file.exists(input_path))

raw <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!("visit" %in% names(raw))) stop("Column 'visit' not found in input data.")
if (!("id" %in% names(raw))) stop("Column 'id' not found in input data.")

raw[["visit"]] <- tolower(as.character(raw[["visit"]]))
for (v in c("edeqpurging", "edeq15obe")) {
  if (!(v %in% names(raw))) stop(sprintf("Column '%s' not found in input data.", v))
  raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
}

baseline <- raw[raw[["visit"]] == "bl", , drop = FALSE]
if (nrow(baseline) == 0) stop("No baseline rows found (visit == 'bl').")

subgroups <- list(
  list(
    label = "bl_anypurging",
    desc = "Baseline edeqpurging > 0",
    id_filter = function(bl) {
      !is.na(bl[["edeqpurging"]]) & bl[["edeqpurging"]] > 0
    }
  ),
  list(
    label = "bl_anybinge",
    desc = "Baseline edeq15obe > 0",
    id_filter = function(bl) {
      !is.na(bl[["edeq15obe"]]) & bl[["edeq15obe"]] > 0
    }
  ),
  list(
    label = "bl_binge_only",
    desc = "Baseline edeqpurging = 0 and edeq15obe > 0 (binge only)",
    id_filter = function(bl) {
      !is.na(bl[["edeqpurging"]]) & bl[["edeqpurging"]] == 0 &
        !is.na(bl[["edeq15obe"]]) & bl[["edeq15obe"]] > 0
    }
  )
)

patch_once <- function(text, pattern, replacement, label) {
  out <- sub(pattern, replacement, text, perl = TRUE)
  if (identical(out, text)) {
    stop(sprintf("Could not patch '%s' in %s", label, source_script))
  }
  out
}

script_text <- paste(readLines(source_script, warn = FALSE), collapse = "\n")
status <- list()

for (sg in subgroups) {
  message("Preparing subgroup: ", sg$label, " (", sg$desc, ")")

  keep_idx <- sg$id_filter(baseline)
  if (!is.logical(keep_idx) || length(keep_idx) != nrow(baseline)) {
    stop("Subgroup id_filter must return a logical vector with length nrow(baseline).")
  }
  ids <- unique(baseline[["id"]][keep_idx])
  ids <- ids[!is.na(ids)]

  if (length(ids) == 0) {
    warning("No participants found for subgroup ", sg$label, "; skipping.")
    status[[length(status) + 1]] <- data.frame(
      subgroup = sg$label,
      description = sg$desc,
      n_ids = 0,
      n_rows = 0,
      run_status = "skipped_no_ids",
      stringsAsFactors = FALSE
    )
    next
  }

  dat_sub <- raw[raw[["id"]] %in% ids, , drop = FALSE]
  out_dir <- file.path(output_root, sg$label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tmp_csv <- tempfile(pattern = paste0("mmrm_", sg$label, "_"),
                      tmpdir = tempdir(), fileext = ".csv")
  write.csv(dat_sub, tmp_csv, row.names = FALSE)

  run_text <- script_text
  run_text <- patch_once(
    run_text,
    'input_path <- file\\.path\\("data", "analysis_dataset\\.csv"\\)',
    sprintf('input_path <- "%s"', gsub("\\\\", "/", tmp_csv)),
    "input_path"
  )
  run_text <- patch_once(
    run_text,
    'results_dir <- file\\.path\\("outputs", "mmrm_main_analysis"\\)',
    sprintf('results_dir <- file.path("outputs", "mmrm_subgroup_analysis", "%s")', sg$label),
    "results_dir"
  )

  tmp_script <- tempfile(pattern = paste0("mmrm_main_analysis_", sg$label, "_"),
                         tmpdir = tempdir(), fileext = ".R")
  writeLines(strsplit(run_text, "\n", fixed = TRUE)[[1]], tmp_script)

  ok <- TRUE
  err_msg <- NA_character_
  tryCatch({
    run_env <- new.env(parent = globalenv())
    sys.source(tmp_script, envir = run_env)
  }, error = function(e) {
    ok <<- FALSE
    err_msg <<- conditionMessage(e)
  })

  status[[length(status) + 1]] <- data.frame(
    subgroup = sg$label,
    description = sg$desc,
    n_ids = length(unique(dat_sub[["id"]])),
    n_rows = nrow(dat_sub),
    run_status = if (ok) "completed" else paste0("failed: ", err_msg),
    stringsAsFactors = FALSE
  )
}

status_df <- do.call(rbind, status)
write.csv(status_df,
          file.path(output_root, "mmrm_subgroup_run_status.csv"),
          row.names = FALSE)

print(status_df)
cat("Subgroup MMRM finished. Outputs under", output_root, "\n")
