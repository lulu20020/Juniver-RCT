# -----------------------------------------------------------------------------
# Script: mmrm_main_analysis.R
# Goal
#   1) Main models (randomised comparison period: wk4/wk8/wk12)
#        Base: outcome ~ GROUP * visit + baseline_outcome
#        cov_base (baseline covariate set) is defined in code as:
#          age, gender, bmi, ethnicity_whitevsother, country
#        Enhanced variants:
#          - main_enhanced: + baseline age, gender, bmi, ethnicity_whitevsother, country
#          - main_enh_ready: main_enhanced + baseline readythoughts + readybehs
#          - main_enh_tx:    main_enhanced + baseline edtx_binary + mhtx_binary
#          - main_enh_all:   main_enhanced + readiness + prior treatment (readythoughts,
#                            readybehs, edtx_binary, mhtx_binary)
#        edeqglobal additionally runs covariance sensitivity (UN→AR(1)→CS);
#        diagnostics saved for edeqglobal base.
#   2) Exploratory models (extended follow-up: wk4/wk8/wk12/wk16/wk20/wk24)
#        - explore_piece (piecewise delayed-start):
#            outcome ~ group + visit_week + post12_week +
#                      group:post12_week + baseline_outcome + cov_base
#            visit_week = numeric week; post12_week = pmax(0, visit_week - 12)
#            group:post12_week tests delayed-start divergence post wk12.
#        - explore_treatweeks (treatment exposure weeks):
#            outcome ~ group + visit_week + treatment_weeks +
#                      group:treatment_weeks + baseline_outcome + cov_base
#            treatment_weeks: Juniver 4/8/12/16/20/24; Control 0/0/0/4/8/12.
#        Exploratory models use baseline covariates cov_base =
#        age, gender, bmi, ethnicity_whitevsother, country; covariance fallback
#        AR(1) → CS → UN. Predictions/contrasts for exploratory models are
#        computed via the model matrix using reference covariates (mean for
#        baseline/bmi; modal level for categorical covariates).
#   3) Save visit-specific contrasts, predicted means, fixed effects,
#      sample sizes, and selected covariance structures. Exploratory outputs use labels
#      explore_piece and explore_treatweeks.
#
# Input:
#   - Long-format analysis dataset. The confidential participant-level dataset
#     is not included in the public repository.
#   - Required columns: id, group, visit, age, gender, country,
#     ethnicity_whitevsother, bmi, and the outcomes listed below.
#
# Outcomes analysed:
#   edeqglobal; EDE-Q subscores: edeq_restraint, edeq_eatconc,
#   edeq_shape, edeq_weight; phq9total, gad7total, cia, pss,
#   edeq15obe, edeq16v, edeq17lax, edeqpurging, bmi
#
# Repeated structure:
#   - Main visits: wk4, wk8, wk12
#   - Subject: id (factor)
#   - Covariance tried in order: UN (~us), AR(1) (~ar1), CS (~cs); first to converge kept
#
# Conventions:
#   - age, gender treated as categorical; covariates use baseline values.
#   - baseline is a covariate, not a repeated level.
#   - If outcome == bmi, bmi is automatically removed from covariates in both
#     fit_one() (main models) and fit_explore() (exploratory models), to avoid
#     adjusting BMI for itself.
#   - No imputation; rows with missing outcome/baseline are dropped per model.
# Dependencies: mmrm, emmeans, dplyr, tidyr, broom
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(mmrm)
library(emmeans)
library(broom)
emm_options(rg.limit = 3000000)

input_path <- file.path("data", "analysis_dataset.csv")
results_dir <- file.path("outputs", "mmrm_main_analysis")
contr_dir <- file.path(results_dir, "contrasts")
pred_dir <- file.path(results_dir, "predicted_means")
ss_dir <- file.path(results_dir, "sample_size")
fe_dir <- file.path(results_dir, "fixed_effects")
cov_dir <- file.path(results_dir, "covariance_selected")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ss_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fe_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cov_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_path))

# Outcomes of interest
outcomes <- c("edeqglobal",
              "edeq_restraint", "edeq_eatconc", "edeq_shape", "edeq_weight",
              "phq9total", "gad7total", "cia", "pss",
              "edeq15obe", "edeq16v", "edeq17lax", "edeqpurging", "bmi")

# Load data
raw <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
raw[outcomes] <- lapply(raw[outcomes], function(x) suppressWarnings(as.numeric(x)))
raw[["visit"]] <- tolower(raw[["visit"]])

# Standardise group labels (some exports use 1/2 instead of names)
if ("group" %in% names(raw)) {
  raw[["group"]] <- trimws(as.character(raw[["group"]]))
  raw[["group"]] <- dplyr::recode(raw[["group"]],
                                  "1" = "Juniver",
                                  "2" = "Control",
                                  .default = raw[["group"]])
}
keep_visits_main <- c("wk4", "wk8", "wk12")
keep_visits_exploratory <- c("wk4", "wk8", "wk12", "wk16", "wk20", "wk24")
cov_base_default <- c("age", "gender", "bmi",
                      "ethnicity_whitevsother", "country")

# Helper: add baseline covariate for outcome
add_baseline <- function(df, outcome) {
  bl_vals <- df %>%
    filter(visit == "bl") %>%
    group_by(id) %>%
    summarise(
      baseline = {
        x <- .data[[outcome]]
        x <- x[!is.na(x)]
        if (length(x) == 0) NA_real_ else x[[1]]
      },
      .groups = "drop"
    )
  df %>% left_join(bl_vals, by = "id")
}

add_baseline_covars <- function(df, covars){
  if (length(covars) == 0) return(df)
  bl_cov <- df %>%
    filter(visit == "bl") %>%
    group_by(id) %>%
    summarise(
      across(all_of(covars), ~ {
        x <- .x
        x <- x[!is.na(x)]
        if (length(x) == 0) NA else x[[1]]
      }),
      .groups = "drop"
    )
  names(bl_cov) <- c("id", paste0(covars, "_bl"))
  out <- df %>% left_join(bl_cov, by = "id")
  for (v in covars){
    blv <- paste0(v, "_bl")
    out[[v]] <- out[[blv]]
    out[[blv]] <- NULL
  }
  out
}

fit_explore <- function(df, outcome,
                        model = c("piece", "treatweeks"),
                        cov_order = c("ar1", "cs", "us"),
                        visits_use = keep_visits_exploratory,
                        covar_set = cov_base_default,
                        return_fit = FALSE) {
  model <- match.arg(model)

  df_o <- add_baseline(df, outcome)
  df_o <- df_o %>%
    filter(visit %in% c(visits_use, "bl")) %>%
    filter(!is.na(.data[[outcome]]), !is.na(baseline))
  if (nrow(df_o) == 0) return(NULL)

  covars_use <- if (is.null(covar_set)) cov_base_default else covar_set
  if (outcome == "bmi") {
    covars_use <- setdiff(covars_use, "bmi")
  }
  df_o <- add_baseline_covars(df_o, covars_use)

  df_o <- df_o %>%
    mutate(
      id = factor(.data[["id"]]),
      visit = factor(visit, levels = visits_use, ordered = FALSE),
    group = factor(.data[["group"]], levels = c("Juniver", "Control")),
      age = factor(.data[["age"]]),
      gender = factor(.data[["gender"]]),
      country = factor(.data[["country"]]),
      ethnicity_whitevsother = factor(.data[["ethnicity_whitevsother"]])
    ) %>% droplevels()
  if (n_distinct(df_o[["group"]]) < 2) return(NULL)

  visit_week_map <- setNames(c(4, 8, 12, 16, 20, 24), visits_use)
  df_o <- df_o %>%
    mutate(
      visit_week = as.numeric(visit_week_map[as.character(visit)]),
      post12_week = pmax(0, visit_week - 12),
      treatment_weeks = ifelse(group == "Juniver",
                               visit_week,
                               pmax(0, visit_week - 12))
    )

  cat_cov <- intersect(c("age", "gender", "country",
                         "ethnicity_whitevsother"), covars_use)
  num_cov <- intersect(c("bmi"), covars_use)
  keep_cat <- cat_cov[sapply(cat_cov, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]
  keep_num <- if (length(num_cov) == 0) character(0) else
    num_cov[sapply(num_cov, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]

  rhs <- switch(model,
                piece = c("group", "visit_week", "post12_week",
                          "group:post12_week", "baseline"),
                treatweeks = c("group", "visit_week", "treatment_weeks",
                               "group:treatment_weeks", "baseline"))
  rhs <- c(rhs, keep_cat, keep_num)

  rhs_covars <- c(keep_cat, keep_num)
  if (length(rhs_covars) > 0) {
    df_o <- df_o %>%
      filter(if_all(all_of(rhs_covars), ~ !is.na(.x))) %>%
      droplevels()
    if (n_distinct(df_o[["group"]]) < 2) return(NULL)
    keep_cat <- keep_cat[sapply(keep_cat, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]
    keep_num <- if (length(keep_num) == 0) character(0) else
      keep_num[sapply(keep_num, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]
    rhs <- switch(model,
                  piece = c("group", "visit_week", "post12_week",
                            "group:post12_week", "baseline", keep_cat, keep_num),
                  treatweeks = c("group", "visit_week", "treatment_weeks",
                                 "group:treatment_weeks", "baseline", keep_cat, keep_num))
  }

  form <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  df_fit <- df_o %>% filter(visit %in% visits_use,
                            !is.na(visit_week))
  if (nrow(df_fit) == 0) return(NULL)

  cov_forms <- list(
    ar1 = ~ ar1(visit | id),
    cs = ~ cs(visit | id),
    us = ~ us(visit | id)
  )

  fit <- NULL
  cov_used <- NA_character_
  for (cov in cov_order) {
    fit <- try(mmrm::mmrm(formula = form,
                          data = df_fit,
                          covariance = cov_forms[[cov]],
                          reml = TRUE),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cov_used <- cov
      break
    }
  }
  if (inherits(fit, "try-error")) return(NULL)

  vis_weeks <- as.numeric(visit_week_map[visits_use])
  visit_label_map <- setNames(visits_use, as.character(vis_weeks))
  grp_levels <- levels(df_fit$group)

  # Reference covariate values: mean for numeric, modal category for factors
  mode_level <- function(x) {
    xt <- table(x)
    names(xt)[which.max(xt)]
  }
  ref_base <- list(
    baseline = mean(df_fit$baseline, na.rm = TRUE),
    bmi = mean(df_fit$bmi, na.rm = TRUE),
    age = factor(mode_level(df_fit$age), levels = levels(df_fit$age)),
    gender = factor(mode_level(df_fit$gender), levels = levels(df_fit$gender)),
    country = factor(mode_level(df_fit$country), levels = levels(df_fit$country)),
    ethnicity_whitevsother = factor(mode_level(df_fit$ethnicity_whitevsother),
                                    levels = levels(df_fit$ethnicity_whitevsother))
  )

  newdata <- expand.grid(
    group = grp_levels,
    visit_week = vis_weeks,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (model == "piece") {
    newdata$post12_week <- pmax(0, newdata$visit_week - 12)
  } else {
    newdata$treatment_weeks <- ifelse(newdata$group == "Juniver",
                                      newdata$visit_week,
                                      pmax(0, newdata$visit_week - 12))
  }
  newdata$baseline <- ref_base$baseline
  newdata$age <- ref_base$age
  newdata$gender <- ref_base$gender
  newdata$country <- ref_base$country
  newdata$ethnicity_whitevsother <- ref_base$ethnicity_whitevsother
  newdata$bmi <- ref_base$bmi
  newdata$visit <- factor(visit_label_map[as.character(newdata$visit_week)],
                          levels = visits_use)
  newdata$group <- factor(newdata$group, levels = grp_levels)

  mm_terms <- delete.response(terms(form))
  Xref <- model.matrix(mm_terms, data = newdata,
                       contrasts.arg = attr(model.matrix(mm_terms, data = df_fit),
                                            "contrasts"))
  beta <- fit$beta_est
  V <- fit$beta_vcov
  # Sparse fits can be rank-deficient, so align the reference
  # design matrix to the estimable coefficient set actually returned by mmrm.
  if (!is.null(names(beta)) && !is.null(colnames(Xref))) {
    miss_cols <- setdiff(names(beta), colnames(Xref))
    if (length(miss_cols) > 0) {
      Xmiss <- matrix(0, nrow = nrow(Xref), ncol = length(miss_cols))
      colnames(Xmiss) <- miss_cols
      Xref <- cbind(Xref, Xmiss)
    }
    Xref <- Xref[, names(beta), drop = FALSE]
  }
  pred_vals <- as.vector(Xref %*% beta)
  pred_se <- sqrt(diag(Xref %*% V %*% t(Xref)))
  preds <- newdata %>%
    mutate(pred = pred_vals, se = pred_se) %>%
    transmute(visit = as.character(visit), group, pred, se)

  # Group contrasts per visit
  cont_list <- lapply(vis_weeks, function(vv) {
    rows <- which(newdata$visit_week == vv)
    if (length(rows) < 2) return(NULL)
    r_ctrl <- rows[newdata$group[rows] == "Control"]
    r_jun <- rows[newdata$group[rows] == "Juniver"]
    if (length(r_ctrl) != 1 || length(r_jun) != 1) return(NULL)
    cvec <- Xref[r_jun, , drop = FALSE] - Xref[r_ctrl, , drop = FALSE]
    est <- as.numeric(cvec %*% beta)
    se <- sqrt(as.numeric(cvec %*% V %*% t(cvec)))
    df_val <- nrow(df_fit) - ncol(Xref)
    tval <- est / se
    pval <- 2 * pt(-abs(tval), df = df_val)
    ci_half <- qt(0.975, df = df_val) * se
    data.frame(visit = visit_label_map[as.character(vv)],
               contrast = est,
               se = se,
               statistic = tval,
               df = df_val,
               conf.low = est - ci_half,
               conf.high = est + ci_half,
               p.value = pval,
               stringsAsFactors = FALSE)
  })
  cont <- bind_rows(cont_list)

  n_total <- n_distinct(df_fit$id)
  n_by_group <- df_fit %>% distinct(id, group) %>% count(group, name = "n")

  fe <- tryCatch({
    tidy(fit, effects = "fixed", conf.int = TRUE)
  }, error = function(e) NULL)

  list(contrasts = cont,
       preds = preds,
       fit = if (return_fit) fit else NULL,
       n_total = n_total,
       n_by_group = n_by_group,
       fe = fe,
       cov_structure = cov_used)
}

fit_one <- function(df, outcome, enhanced = FALSE,
                    cov_order = c("us", "ar1", "cs"),
                    return_fit = FALSE,
                    visits_use = keep_visits_main,
                    covar_set = NULL){
  df_o <- add_baseline(df, outcome)
  df_o <- df_o %>%
    filter(visit %in% c(visits_use, "bl")) %>%
    filter(!is.na(.data[[outcome]]), !is.na(baseline))
  if (nrow(df_o) == 0) return(NULL)

  base_covars_default <- c("age", "gender", "bmi",
                           "ethnicity_whitevsother", "country")
  covars_use <- if (is.null(covar_set)) base_covars_default else covar_set
  if (outcome == "bmi") {
    covars_use <- setdiff(covars_use, "bmi")
  }
  df_o <- add_baseline_covars(df_o, covars_use)

  df_o <- df_o %>% mutate(
    id = factor(.data[["id"]]),
    visit = factor(visit, levels = visits_use, ordered = FALSE),
    group = factor(.data[["group"]], levels = c("Juniver", "Control")),
    age = factor(.data[["age"]]),
    gender = factor(.data[["gender"]]),
    country = factor(.data[["country"]]),
    ethnicity_whitevsother = factor(.data[["ethnicity_whitevsother"]])
  ) %>% droplevels()
  if (n_distinct(df_o[["group"]]) < 2) return(NULL)

  cat_cov <- intersect(c("age", "gender", "country", "ethnicity_whitevsother",
                         "readythoughts", "readybehs",
                         "edtx_binary", "mhtx_binary"), covars_use)
  num_cov <- intersect(c("bmi"), covars_use)
  keep_cat <- cat_cov[sapply(cat_cov, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]
  keep_num <- if (length(num_cov) == 0) character(0) else
    num_cov[sapply(num_cov, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]

  rhs_base <- c("group*visit", "baseline")
  rhs <- rhs_base
  if (enhanced) {
    rhs <- c(rhs, keep_cat, keep_num)
    if (length(rhs) > 0) {
      df_o <- df_o %>%
        filter(if_all(all_of(setdiff(rhs, c(rhs_base))), ~ !is.na(.x))) %>%
        droplevels()
      if (n_distinct(df_o[["group"]]) < 2) return(NULL)
      keep_cat <- keep_cat[sapply(keep_cat, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]
      keep_num <- if (length(keep_num) == 0) character(0) else
        keep_num[sapply(keep_num, function(v) n_distinct(df_o[[v]], na.rm = TRUE) > 1)]
      rhs <- c(rhs_base, keep_cat, keep_num)
    }
  }

  form <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))

  # Fit with covariance fallbacks
  fit <- NULL
  cov_used <- NA_character_
  cov_forms <- list(
    us = ~ us(visit | id),
    ar1 = ~ ar1(visit | id),
    cs = ~ cs(visit | id)
  )
  for (cov in cov_order) {
    fit <- try(mmrm::mmrm(formula = form,
                          data = df_o %>% filter(visit %in% visits_use),
                          covariance = cov_forms[[cov]],
                          reml = TRUE),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cov_used <- cov
      break
    }
  }
  if (inherits(fit, "try-error")) return(NULL)

  preds <- data.frame()
  cont <- data.frame()

  # Estimated marginal means and contrasts
  emm_try <- try(emmeans(fit, ~ group | visit, nesting = NULL), silent = TRUE)
  if (inherits(emm_try, "try-error") || length(levels(df_o$group)) < 2) return(NULL)
  cont <- contrast(emm_try, method = list(diff = c(1, -1))) %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    transmute(visit = as.character(visit),
              contrast = estimate,
              se = SE,
              statistic = t.ratio,
              df,
              conf.low = lower.CL,
              conf.high = upper.CL,
              p.value)

  preds <- emm_try %>% summary() %>% as.data.frame() %>%
    transmute(visit = as.character(visit), group, pred = emmean, se = SE)

  if (is.null(preds) || is.null(cont)) return(NULL)
  if ("visit" %in% names(cont)) cont$visit <- as.character(cont$visit)
  if ("visit" %in% names(preds)) preds$visit <- as.character(preds$visit)

  df_used <- df_o %>% filter(visit %in% visits_use)
  n_total <- n_distinct(df_used$id)
  n_by_group <- df_used %>% distinct(id, group) %>% count(group, name = "n")

  fe <- tryCatch({
    tidy(fit, effects = "fixed", conf.int = TRUE)
  }, error = function(e) NULL)

  list(contrasts = cont,
       preds = preds,
       fit = if (return_fit) fit else NULL,
       n_total = n_total,
       n_by_group = n_by_group,
       fe = fe,
       cov_structure = cov_used)
}

all_contrasts <- list()
all_preds <- list()
all_sizes <- list()
all_fixed <- list()
all_cov <- list()
primary_base <- NULL

for (outcome in outcomes) {
  message("Main base model: ", outcome)
  res_base <- fit_one(raw, outcome, enhanced = FALSE,
                      return_fit = (outcome == "edeqglobal"),
                      visits_use = keep_visits_main)
  if (!is.null(res_base)) {
    res_base$contrasts$outcome <- outcome
    res_base$contrasts$model <- "main_base"
    all_contrasts[[length(all_contrasts) + 1]] <- res_base$contrasts
    if (!is.null(res_base$preds)) {
      res_base$preds$outcome <- outcome
      res_base$preds$model <- "main_base"
      all_preds[[length(all_preds) + 1]] <- res_base$preds
    }
    if (!is.null(res_base$fe)) {
      res_base$fe$outcome <- outcome
      res_base$fe$model <- "main_base"
      all_fixed[[length(all_fixed) + 1]] <- res_base$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "main_base",
                 n_total = res_base$n_total) %>%
      bind_rows(res_base$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "main_base",
                            group = as.character(group),
                            n = n))
    if (outcome == "edeqglobal") primary_base <- res_base
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "main_base",
      cov_structure = res_base$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  # Enhanced variants
  cov_base <- c("age", "gender", "bmi",
                "ethnicity_whitevsother", "country")
  cov_ready <- c(cov_base, "readythoughts", "readybehs")
  cov_tx <- c(cov_base, "edtx_binary", "mhtx_binary")
  cov_all <- c(cov_base, "readythoughts", "readybehs",
               "edtx_binary", "mhtx_binary")

  message("Main enhanced model: ", outcome)
  res_enh <- fit_one(raw, outcome, enhanced = TRUE,
                     visits_use = keep_visits_main,
                     covar_set = cov_base)
  if (!is.null(res_enh)) {
    res_enh$contrasts$outcome <- outcome
    res_enh$contrasts$model <- "main_enhanced"
    all_contrasts[[length(all_contrasts) + 1]] <- res_enh$contrasts
    if (!is.null(res_enh$fe)) {
      res_enh$fe$outcome <- outcome
      res_enh$fe$model <- "main_enhanced"
      all_fixed[[length(all_fixed) + 1]] <- res_enh$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "main_enhanced",
                 n_total = res_enh$n_total) %>%
      bind_rows(res_enh$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "main_enhanced",
                            group = as.character(group),
                            n = n))
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "main_enhanced",
      cov_structure = res_enh$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  message("Enhanced + readiness: ", outcome)
  res_ready <- fit_one(raw, outcome, enhanced = TRUE,
                       visits_use = keep_visits_main,
                       covar_set = cov_ready)
  if (!is.null(res_ready)) {
    res_ready$contrasts$outcome <- outcome
    res_ready$contrasts$model <- "main_enh_ready"
    all_contrasts[[length(all_contrasts) + 1]] <- res_ready$contrasts
    if (!is.null(res_ready$fe)) {
      res_ready$fe$outcome <- outcome
      res_ready$fe$model <- "main_enh_ready"
      all_fixed[[length(all_fixed) + 1]] <- res_ready$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "main_enh_ready",
                 n_total = res_ready$n_total) %>%
      bind_rows(res_ready$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "main_enh_ready",
                            group = as.character(group),
                            n = n))
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "main_enh_ready",
      cov_structure = res_ready$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  message("Enhanced + prior tx: ", outcome)
  res_tx <- fit_one(raw, outcome, enhanced = TRUE,
                    visits_use = keep_visits_main,
                    covar_set = cov_tx)
  if (!is.null(res_tx)) {
    res_tx$contrasts$outcome <- outcome
    res_tx$contrasts$model <- "main_enh_tx"
    all_contrasts[[length(all_contrasts) + 1]] <- res_tx$contrasts
    if (!is.null(res_tx$fe)) {
      res_tx$fe$outcome <- outcome
      res_tx$fe$model <- "main_enh_tx"
      all_fixed[[length(all_fixed) + 1]] <- res_tx$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "main_enh_tx",
                 n_total = res_tx$n_total) %>%
      bind_rows(res_tx$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "main_enh_tx",
                            group = as.character(group),
                            n = n))
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "main_enh_tx",
      cov_structure = res_tx$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  message("Enhanced + readiness + prior tx: ", outcome)
  res_all <- fit_one(raw, outcome, enhanced = TRUE,
                     visits_use = keep_visits_main,
                     covar_set = cov_all)
  if (!is.null(res_all)) {
    res_all$contrasts$outcome <- outcome
    res_all$contrasts$model <- "main_enh_all"
    all_contrasts[[length(all_contrasts) + 1]] <- res_all$contrasts
    if (!is.null(res_all$fe)) {
      res_all$fe$outcome <- outcome
      res_all$fe$model <- "main_enh_all"
      all_fixed[[length(all_fixed) + 1]] <- res_all$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "main_enh_all",
                 n_total = res_all$n_total) %>%
      bind_rows(res_all$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "main_enh_all",
                            group = as.character(group),
                            n = n))
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "main_enh_all",
      cov_structure = res_all$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  # Exploratory: piecewise delayed-start (wk4–wk24)
  message("Exploratory piecewise (wk4-24): ", outcome)
  res_piece <- fit_explore(raw, outcome,
                           model = "piece",
                           cov_order = c("ar1", "cs", "us"),
                           visits_use = keep_visits_exploratory,
                           covar_set = cov_base_default)
  if (!is.null(res_piece)) {
    res_piece$contrasts$outcome <- outcome
    res_piece$contrasts$model <- "explore_piece"
    all_contrasts[[length(all_contrasts) + 1]] <- res_piece$contrasts
    if (!is.null(res_piece$preds)) {
      res_piece$preds$outcome <- outcome
      res_piece$preds$model <- "explore_piece"
      all_preds[[length(all_preds) + 1]] <- res_piece$preds
    }
    if (!is.null(res_piece$fe)) {
      res_piece$fe$outcome <- outcome
      res_piece$fe$model <- "explore_piece"
      all_fixed[[length(all_fixed) + 1]] <- res_piece$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "explore_piece",
                 n_total = res_piece$n_total) %>%
      bind_rows(res_piece$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "explore_piece",
                            group = as.character(group),
                            n = n))
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "explore_piece",
      cov_structure = res_piece$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  # Exploratory: treatment exposure weeks (wk4–wk24)
  message("Exploratory treatment-weeks (wk4-24): ", outcome)
  res_treat <- fit_explore(raw, outcome,
                           model = "treatweeks",
                           cov_order = c("ar1", "cs", "us"),
                           visits_use = keep_visits_exploratory,
                           covar_set = cov_base_default)
  if (!is.null(res_treat)) {
    res_treat$contrasts$outcome <- outcome
    res_treat$contrasts$model <- "explore_treatweeks"
    all_contrasts[[length(all_contrasts) + 1]] <- res_treat$contrasts
    if (!is.null(res_treat$preds)) {
      res_treat$preds$outcome <- outcome
      res_treat$preds$model <- "explore_treatweeks"
      all_preds[[length(all_preds) + 1]] <- res_treat$preds
    }
    if (!is.null(res_treat$fe)) {
      res_treat$fe$outcome <- outcome
      res_treat$fe$model <- "explore_treatweeks"
      all_fixed[[length(all_fixed) + 1]] <- res_treat$fe
    }
    all_sizes[[length(all_sizes) + 1]] <-
      data.frame(outcome = outcome, model = "explore_treatweeks",
                 n_total = res_treat$n_total) %>%
      bind_rows(res_treat$n_by_group %>%
                  transmute(outcome = outcome,
                            model = "explore_treatweeks",
                            group = as.character(group),
                            n = n))
    all_cov[[length(all_cov) + 1]] <- data.frame(
      outcome = outcome,
      model = "explore_treatweeks",
      cov_structure = res_treat$cov_structure,
      stringsAsFactors = FALSE
    )
  }

  # Sensitivity covariance for primary outcome (main)
  if (outcome == "edeqglobal") {
    for (covx in c("ar1", "cs")) {
      message("Sensitivity (main ", covx, "): ", outcome)
      res_sens <- fit_one(raw, outcome, enhanced = FALSE,
                          cov_order = covx,
                          visits_use = keep_visits_main)
      if (!is.null(res_sens)) {
      res_sens$contrasts$outcome <- outcome
      res_sens$contrasts$model <- paste0("main_base_", covx)
      all_contrasts[[length(all_contrasts) + 1]] <- res_sens$contrasts
      if (!is.null(res_sens$fe)) {
        res_sens$fe$outcome <- outcome
        res_sens$fe$model <- paste0("main_base_", covx)
        all_fixed[[length(all_fixed) + 1]] <- res_sens$fe
      }
      all_sizes[[length(all_sizes) + 1]] <-
        data.frame(outcome = outcome, model = paste0("main_base_", covx),
                   n_total = res_sens$n_total) %>%
        bind_rows(res_sens$n_by_group %>%
                    transmute(outcome = outcome,
                                model = paste0("main_base_", covx),
                                group = as.character(group),
                                n = n))
      all_cov[[length(all_cov) + 1]] <- data.frame(
        outcome = outcome,
        model = paste0("main_base_", covx),
        cov_structure = res_sens$cov_structure,
        stringsAsFactors = FALSE
      )
      }
    }
  }
}

cov_out <- bind_rows(all_cov)
if (nrow(cov_out) > 0) {
  cov_out <- cov_out %>% relocate(outcome)
  write.csv(cov_out,
            file.path(cov_dir, "mmrm_covariance_selected.csv"),
            row.names = FALSE)
}

contr_out <- bind_rows(all_contrasts)
if (nrow(contr_out) > 0) {
  contr_out <- contr_out %>% relocate(outcome)
  write.csv(contr_out,
            file.path(contr_dir, "mmrm_group_contrasts.csv"),
            row.names = FALSE)
  # Split by model for clarity.
  for (m in unique(contr_out$model)) {
    dir_m <- file.path(contr_dir, m)
    dir.create(dir_m, showWarnings = FALSE, recursive = TRUE)
    write.csv(contr_out %>% filter(model == m),
              file.path(dir_m, paste0("mmrm_group_contrasts_", m, ".csv")),
              row.names = FALSE)
    contr_out %>% filter(model == m) %>%
      group_by(outcome) %>%
      group_walk(~ write.csv(.x,
                             file.path(dir_m,
                                       paste0("contrasts_", .y$outcome, ".csv")),
                                       row.names = FALSE))
  }
}

size_out <- bind_rows(all_sizes)
if (nrow(size_out) > 0) {
  size_out <- size_out %>% relocate(outcome)
  write.csv(size_out, file.path(ss_dir, "mmrm_sample_sizes.csv"),
            row.names = FALSE)
  for (m in unique(size_out$model)) {
    write.csv(size_out %>% filter(model == m),
              file.path(ss_dir, paste0("mmrm_sample_sizes_", m, ".csv")),
              row.names = FALSE)
    dir_m <- file.path(ss_dir, m)
    dir.create(dir_m, showWarnings = FALSE, recursive = TRUE)
    size_out %>% filter(model == m) %>%
      group_by(outcome) %>%
      group_walk(~ write.csv(.x,
                             file.path(dir_m,
                                       paste0("sample_size_", .y$outcome, ".csv")),
                             row.names = FALSE))
  }
}

fe_out <- bind_rows(all_fixed)
if (nrow(fe_out) > 0) {
  fe_out <- fe_out %>% relocate(outcome)
  write.csv(fe_out, file.path(fe_dir, "mmrm_fixed_effects.csv"),
            row.names = FALSE)
  for (m in unique(fe_out$model)) {
    write.csv(fe_out %>% filter(model == m),
              file.path(fe_dir, paste0("mmrm_fixed_effects_", m, ".csv")),
              row.names = FALSE)
    # Per-outcome files for this model.
    dir_m <- file.path(fe_dir, m)
    dir.create(dir_m, showWarnings = FALSE, recursive = TRUE)
    fe_out %>% filter(model == m) %>%
      group_by(outcome) %>%
      group_walk(~ write.csv(.x, file.path(dir_m,
                                          paste0("fixedeffect_", .y$outcome, ".csv")),
                             row.names = FALSE))
  }
}

pred_out <- bind_rows(all_preds)
if (nrow(pred_out) > 0) {
  pred_out <- pred_out %>% relocate(outcome)
  write.csv(pred_out, file.path(pred_dir, "mmrm_predicted_means.csv"),
            row.names = FALSE)
  for (m in unique(pred_out$model)) {
    dir_m <- file.path(pred_dir, m)
    dir.create(dir_m, showWarnings = FALSE, recursive = TRUE)
    write.csv(pred_out %>% filter(model == m),
              file.path(dir_m, paste0("mmrm_predicted_means_", m, ".csv")),
              row.names = FALSE)
    pred_out %>% filter(model == m) %>%
      group_by(outcome) %>%
      group_walk(~ write.csv(.x,
                             file.path(dir_m,
                                       paste0("predicted_means_", .y$outcome, ".csv")),
                                       row.names = FALSE))
  }
}

cat("MMRM finished. Outputs in", results_dir, "\n")
