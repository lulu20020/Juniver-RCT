# -----------------------------------------------------------------------------
# Script: mnar_sensitivity_analysis.R
# Purpose:
#   MNAR sensitivity analyses for the primary outcome (EDEQglobal) during the
#   randomised period (wk4/wk8/wk12), run under two model specifications:
#     1) No-covariate model:
#          outcome ~ group * visit + baseline
#     2) Covariate-adjusted model:
#          outcome ~ group * visit + baseline +
#                    age + gender + bmi + ethnicity_whitevsother + country
#
#   For each specification, this script runs:
#     A) Pattern-mixture delta tipping-point analysis (group-specific deltas)
#        on unobserved wk12 outcomes only.
#     B) Reference-based MI (rbmi; J2R / CR / CIR).
#
# Notes:
#   - Input data are not included in the public repository.
#   - Visits used: bl, wk4, wk8, wk12.
#   - Observed outcomes are not changed.
#   - Only missing wk12 outcomes are replaced by MAR predictions + delta.
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(mmrm)
library(emmeans)
library(rbmi)

input_path <- file.path("data", "analysis_dataset.csv")
stopifnot(file.exists(input_path))

res_dir <- file.path("outputs", "mnar_sensitivity_analysis")
tab_dir <- file.path(res_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)

raw$group <- trimws(as.character(raw$group))
raw$group <- dplyr::recode(raw$group,
                           "1" = "Juniver", "2" = "Control",
                           .default = raw$group)
raw$visit <- tolower(raw$visit)

keep_visits <- c("wk4", "wk8", "wk12")
outcome <- "edeqglobal"
covars_all <- c("age", "gender", "bmi", "ethnicity_whitevsother", "country")

add_baseline <- function(df, vars) {
  bl <- df %>%
    filter(visit == "bl") %>%
    group_by(id) %>%
    summarise(
      baseline = {
        x <- .data[[outcome]]
        x <- x[!is.na(x)]
        if (length(x) == 0) NA_real_ else as.numeric(x[[1]])
      },
      across(all_of(vars), ~ {
        x <- .x
        x <- x[!is.na(x)]
        if (length(x) == 0) NA else x[[1]]
      }),
      .groups = "drop"
    )

  out <- df %>% left_join(bl, by = "id")
  for (v in vars) {
    out[[v]] <- dplyr::coalesce(out[[v]], bl[[v]][match(out$id, bl$id)])
  }
  out
}

dat <- raw %>%
  filter(visit %in% c("bl", keep_visits)) %>%
  add_baseline(vars = covars_all) %>%
  mutate(
    id = factor(id),
    visit = factor(visit, levels = keep_visits, ordered = FALSE),
    group = factor(group, levels = c("Juniver", "Control")),
    baseline = as.numeric(baseline),
    !!outcome := as.numeric(.data[[outcome]]),
    across(all_of(covars_all), ~ ifelse(is.character(.x), trimws(.x), .x)),
    age = factor(age),
    gender = factor(gender),
    country = factor(country),
    ethnicity_whitevsother = factor(ethnicity_whitevsother),
    bmi = as.numeric(bmi)
  )

delta_grid <- expand.grid(
  delta_j = c(0, 0.25, 0.5, 0.75, 1.0),
  delta_c = c(0, 0.25, 0.5, 0.75, 1.0)
)

run_delta_grid <- function(dat, spec_name, use_covariates) {
  covars_use <- if (isTRUE(use_covariates)) covars_all else character(0)

  df_fit_base <- dat %>%
    filter(visit %in% keep_visits,
           !is.na(.data[[outcome]]),
           !is.na(baseline))

  if (length(covars_use) > 0) {
    cov_keep <- covars_use[sapply(covars_use, function(v) {
      dplyr::n_distinct(df_fit_base[[v]], na.rm = TRUE) > 1
    })]
    cov_keep <- as.character(cov_keep)
  } else {
    cov_keep <- character(0)
  }

  if (length(cov_keep) > 0) {
    df_fit_base <- df_fit_base %>%
      tidyr::drop_na(all_of(cov_keep))
  }

  if (nrow(df_fit_base) == 0) {
    stop("No analysable rows for spec: ", spec_name)
  }

  primary_form <- as.formula(
    paste0(outcome, " ~ group*visit + baseline",
           if (length(cov_keep) > 0) {
             paste0(" + ", paste(cov_keep, collapse = " + "))
           } else {
             ""
           })
  )

  fit_mar <- mmrm::mmrm(
    formula = primary_form,
    data = df_fit_base,
    covariance = ~ us(visit | id),
    reml = TRUE
  )

  wk12_missing <- dat %>%
    filter(visit == "wk12", is.na(.data[[outcome]]))

  if (nrow(wk12_missing) > 0) {
    pred_tmp <- try(predict(fit_mar, newdata = wk12_missing), silent = TRUE)
    fallback <- mean(df_fit_base[[outcome]], na.rm = TRUE)
    if (inherits(pred_tmp, "try-error")) {
      wk12_missing$pred_mar <- fallback
    } else {
      pred_num <- as.numeric(pred_tmp)
      pred_num[!is.finite(pred_num)] <- fallback
      wk12_missing$pred_mar <- pred_num
    }
  }

  run_one_delta <- function(dj, dc) {
    df_imp <- dat

    if (nrow(wk12_missing) > 0) {
      df_imp$outcome_adj <- df_imp[[outcome]]
      wk12_missing2 <- wk12_missing %>%
        mutate(
          delta = ifelse(group == "Juniver", dj, dc),
          outcome_mnar = pred_mar + delta
        )

      df_imp <- df_imp %>%
        left_join(wk12_missing2 %>% select(id, outcome_mnar), by = "id") %>%
        mutate(
          outcome_adj = ifelse(
            visit == "wk12" & is.na(.data[[outcome]]),
            outcome_mnar,
            outcome_adj
          )
        ) %>%
        select(-outcome_mnar)
    } else {
      df_imp$outcome_adj <- df_imp[[outcome]]
    }

    df_imp2 <- df_imp %>%
      filter(visit %in% keep_visits, !is.na(outcome_adj), !is.na(baseline))

    cov_keep_imp <- as.character(unlist(intersect(cov_keep, names(df_imp2))))
    if (length(cov_keep_imp) > 0) {
      keep_idx <- sapply(cov_keep_imp, function(v) {
        dplyr::n_distinct(df_imp2[[v]], na.rm = TRUE) > 1
      })
      cov_keep_imp <- cov_keep_imp[as.logical(keep_idx)]
      if (length(cov_keep_imp) > 0) {
        df_imp2 <- df_imp2 %>% tidyr::drop_na(all_of(cov_keep_imp))
      }
    }

    if (nrow(df_imp2) == 0) {
      return(tibble(
        spec = spec_name,
        delta_j = dj,
        delta_c = dc,
        estimate = NA_real_,
        se = NA_real_,
        df = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_
      ))
    }

    form_imp <- as.formula(
      paste0("outcome_adj ~ group*visit + baseline",
             if (length(cov_keep_imp) > 0) {
               paste0(" + ", paste(cov_keep_imp, collapse = " + "))
             } else {
               ""
             })
    )

    fit_imp <- try(mmrm::mmrm(
      formula = form_imp,
      data = df_imp2,
      covariance = ~ us(visit | id),
      reml = TRUE
    ), silent = TRUE)

    if (inherits(fit_imp, "try-error")) {
      return(tibble(
        spec = spec_name,
        delta_j = dj,
        delta_c = dc,
        estimate = NA_real_,
        se = NA_real_,
        df = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_
      ))
    }

    emm <- emmeans(fit_imp, ~ group | visit)
    cont <- contrast(emm, method = list(diff = c(1, -1))) %>%
      summary(infer = TRUE)
    row_wk12 <- cont[as.character(cont$visit) == "wk12", , drop = FALSE]

    if (nrow(row_wk12) != 1) {
      return(tibble(
        spec = spec_name,
        delta_j = dj,
        delta_c = dc,
        estimate = NA_real_,
        se = NA_real_,
        df = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_
      ))
    }

    tibble(
      spec = spec_name,
      delta_j = dj,
      delta_c = dc,
      estimate = row_wk12$estimate,
      se = row_wk12$SE,
      df = row_wk12$df,
      conf.low = row_wk12$lower.CL,
      conf.high = row_wk12$upper.CL,
      p.value = row_wk12$p.value
    )
  }

  tp_out <- bind_rows(lapply(seq_len(nrow(delta_grid)), function(i) {
    run_one_delta(delta_grid$delta_j[i], delta_grid$delta_c[i])
  }))

  out_csv <- file.path(tab_dir, paste0("delta_tipping_", spec_name, ".csv"))
  write.csv(tp_out, out_csv, row.names = FALSE)

  tp_out
}

rbmi_run_one <- function(dat, strategy_active, spec_name, use_covariates,
                         out_csv, n_imp = 200) {
  strategy_map <- c("J2R" = "JR", "CR" = "CR", "CIR" = "CIR")
  strategy_rbmi <- strategy_map[[strategy_active]]
  if (is.null(strategy_rbmi)) {
    stop("Unknown strategy: ", strategy_active)
  }

  covars_use <- if (isTRUE(use_covariates)) covars_all else character(0)

  rbmi_dat <- dat %>%
    dplyr::filter(visit %in% keep_visits) %>%
    dplyr::transmute(
      id = factor(id),
      visit_num = factor(match(as.character(visit), keep_visits), levels = 1:3),
      group = factor(as.character(group), levels = c("Control", "Juniver")),
      outcome_bl = as.numeric(baseline),
      outcome = as.numeric(.data[[outcome]]),
      change = outcome - outcome_bl,
      age = age,
      gender = gender,
      bmi = as.numeric(bmi),
      ethnicity_whitevsother = ethnicity_whitevsother,
      country = country
    ) %>%
    dplyr::arrange(id, visit_num, dplyr::desc(!is.na(outcome))) %>%
    dplyr::group_by(id, visit_num) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tidyr::complete(id, visit_num) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      group = {
        x <- group[!is.na(group)]
        if (length(x) == 0) NA else x[[1]]
      },
      outcome_bl = {
        x <- outcome_bl[is.finite(outcome_bl)]
        if (length(x) == 0) NA_real_ else as.numeric(x[[1]])
      },
      age = {
        x <- age[!is.na(age)]
        if (length(x) == 0) NA else x[[1]]
      },
      gender = {
        x <- gender[!is.na(gender)]
        if (length(x) == 0) NA else x[[1]]
      },
      bmi = {
        x <- bmi[is.finite(bmi)]
        if (length(x) == 0) NA_real_ else as.numeric(x[[1]])
      },
      ethnicity_whitevsother = {
        x <- ethnicity_whitevsother[!is.na(ethnicity_whitevsother)]
        if (length(x) == 0) NA else x[[1]]
      },
      country = {
        x <- country[!is.na(country)]
        if (length(x) == 0) NA else x[[1]]
      }
    ) %>%
    dplyr::ungroup()

  keep_flag <- rbmi_dat %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      has_bl = any(is.finite(outcome_bl)),
      n_obs_visits = dplyr::n_distinct(visit_num[is.finite(change)]),
      .groups = "drop"
    )

  removed_bl <- keep_flag %>% dplyr::filter(!has_bl)
  write.csv(
    removed_bl,
    file.path(tab_dir, paste0("rbmi_removed_missing_baseline_", spec_name, ".csv")),
    row.names = FALSE
  )

  removed_single <- keep_flag %>% dplyr::filter(has_bl, n_obs_visits < 2)
  write.csv(
    removed_single,
    file.path(tab_dir, paste0("rbmi_removed_single_visit_", spec_name, ".csv")),
    row.names = FALSE
  )

  keep_ids <- keep_flag %>%
    dplyr::filter(has_bl, n_obs_visits >= 2) %>%
    dplyr::pull(id)

  rbmi_dat <- rbmi_dat %>%
    dplyr::filter(id %in% keep_ids) %>%
    droplevels()

  if (length(covars_use) > 0) {
    cov_keep_rbmi <- covars_use[sapply(covars_use, function(v) {
      dplyr::n_distinct(rbmi_dat[[v]], na.rm = TRUE) > 1
    })]
    cov_keep_rbmi <- as.character(cov_keep_rbmi)
  } else {
    cov_keep_rbmi <- character(0)
  }

  if (length(cov_keep_rbmi) > 0) {
    rbmi_dat <- rbmi_dat %>% tidyr::drop_na(dplyr::all_of(cov_keep_rbmi))
  }

  if (length(covars_use) > 0) {
    cov_keep_rbmi <- covars_use[sapply(covars_use, function(v) {
      dplyr::n_distinct(rbmi_dat[[v]], na.rm = TRUE) > 1
    })]
    cov_keep_rbmi <- as.character(cov_keep_rbmi)
  } else {
    cov_keep_rbmi <- character(0)
  }

  before_base <- dat %>%
    dplyr::filter(visit == "bl") %>%
    dplyr::transmute(
      id = factor(id),
      group = factor(as.character(group), levels = c("Control", "Juniver")),
      outcome_bl = as.numeric(.data[[outcome]]),
      age = age,
      gender = gender,
      bmi = as.numeric(bmi),
      ethnicity_whitevsother = ethnicity_whitevsother,
      country = country
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(set = "before")

  after_base <- rbmi_dat %>%
    dplyr::select(id, group, outcome_bl, age, gender, bmi,
                  ethnicity_whitevsother, country) %>%
    dplyr::distinct() %>%
    dplyr::mutate(set = "after")

  base_summary <- dplyr::bind_rows(before_base, after_base)
  write.csv(
    base_summary,
    file.path(tab_dir, paste0("rbmi_baseline_before_after_", spec_name, ".csv")),
    row.names = FALSE
  )

  data_ice <- rbmi_dat %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      group = {
        x <- group[!is.na(group)]
        if (length(x) == 0) NA else x[[1]]
      },
      ice_visit = {
        v <- as.integer(as.character(visit_num))
        y <- change
        vv <- v[is.na(y)]
        if (length(vv) == 0) NA_integer_ else min(vv)
      },
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(ice_visit)) %>%
    dplyr::transmute(
      id = id,
      visit_num = factor(ice_visit, levels = 1:3),
      strategy = ifelse(group == "Juniver", strategy_rbmi, "MAR")
    )

  covariates_imp <- c("visit_num*outcome_bl", "visit_num*group", cov_keep_rbmi)
  covariates_an <- c("outcome_bl", cov_keep_rbmi)

  vars_imp <- rbmi::set_vars(
    subjid = "id",
    visit = "visit_num",
    outcome = "change",
    group = "group",
    covariates = covariates_imp,
    strategy = "strategy"
  )
  vars_an <- vars_imp
  vars_an$covariates <- covariates_an

  method <- rbmi::method_approxbayes(n_samples = n_imp)
  draw_obj <- rbmi::draws(
    data = rbmi_dat,
    data_ice = data_ice,
    vars = vars_imp,
    method = method,
    quiet = TRUE,
    ncores = 1
  )
  imp_obj <- rbmi::impute(
    draw_obj,
    references = c("Control" = "Control", "Juniver" = "Control")
  )
  ana_obj <- rbmi::analyse(imp_obj, vars = vars_an)
  pool_obj <- rbmi::pool(ana_obj)

  extract_trt <- function(pool_obj) {
    if (!is.null(pool_obj$results)) {
      df <- as.data.frame(pool_obj$results)
      row <- df[df$parameter == "trt_3", , drop = FALSE]
      if (nrow(row) != 1) return(NULL)
      pval <- if (!is.null(row$pvalue)) row$pvalue else row$pval
      return(list(
        est = row$est,
        se = row$se,
        ci = c(row$lci, row$uci),
        pvalue = pval
      ))
    }
    if (!is.null(pool_obj$pars)) {
      row <- pool_obj$pars[["trt_3"]]
      if (is.null(row)) return(NULL)
      pval <- if (!is.null(row$pvalue)) row$pvalue else row$pval
      return(list(
        est = row$est,
        se = row$se,
        ci = row$ci,
        pvalue = pval
      ))
    }
    NULL
  }

  trt_row <- extract_trt(pool_obj)
  if (is.null(trt_row)) {
    stop("rbmi pooled results did not contain trt_3 estimate (", spec_name, ")")
  }

  out <- tibble::tibble(
    spec = spec_name,
    strategy = strategy_active,
    visit = "wk12",
    contrast = "Juniver - Control",
    estimate = trt_row$est,
    se = trt_row$se,
    conf.low = trt_row$ci[[1]],
    conf.high = trt_row$ci[[2]],
    p.value = trt_row$pvalue,
    n_imputations = n_imp
  )

  write.csv(out, out_csv, row.names = FALSE)
  out
}

specs <- data.frame(
  spec = c("no_covariates", "with_covariates"),
  use_covariates = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

all_delta <- list()
all_rbmi <- list()

for (i in seq_len(nrow(specs))) {
  spec_name <- specs$spec[i]
  use_covs <- isTRUE(specs$use_covariates[i])
  message("Running delta tipping-point: ", spec_name)
  all_delta[[length(all_delta) + 1]] <- run_delta_grid(
    dat = dat,
    spec_name = spec_name,
    use_covariates = use_covs
  )
}

delta_all <- bind_rows(all_delta)
write.csv(delta_all,
          file.path(tab_dir, "delta_tipping_all_specs.csv"),
          row.names = FALSE)

for (i in seq_len(nrow(specs))) {
  spec_name <- specs$spec[i]
  use_covs <- isTRUE(specs$use_covariates[i])
  for (strategy_i in c("J2R", "CR", "CIR")) {
    message("Running rbmi ", strategy_i, ": ", spec_name)
    out_csv <- file.path(
      tab_dir,
      paste0("rbi_", tolower(strategy_i), "_", spec_name, ".csv")
    )
    rbmi_out <- rbmi_run_one(
      dat = dat,
      strategy_active = strategy_i,
      spec_name = spec_name,
      use_covariates = use_covs,
      out_csv = out_csv,
      n_imp = 200
    )
    all_rbmi[[length(all_rbmi) + 1]] <- rbmi_out
  }
}

rbmi_all <- bind_rows(all_rbmi)
write.csv(rbmi_all,
          file.path(tab_dir, "rbi_all_specs.csv"),
          row.names = FALSE)

# Main-result comparison table for manuscript reporting:
# no-covariates vs with-covariates at wk12.
main_cmp <- delta_all %>%
  filter(delta_j == 0, delta_c == 0) %>%
  transmute(
    spec = spec,
    method = "delta_mar_reference",
    estimate = estimate,
    conf.low = conf.low,
    conf.high = conf.high,
    p.value = p.value
  )

main_cmp <- bind_rows(
  main_cmp,
  rbmi_all %>%
    transmute(
      spec = spec,
      method = paste0("rbmi_", tolower(strategy)),
      estimate = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value
    )
)

main_cmp <- main_cmp %>%
  mutate(
    spec = factor(spec, levels = c("no_covariates", "with_covariates")),
    method = factor(method, levels = c("delta_mar_reference",
                                       "rbmi_j2r", "rbmi_cr", "rbmi_cir"))
  ) %>%
  arrange(spec, method) %>%
  mutate(spec = as.character(spec), method = as.character(method))

write.csv(main_cmp,
          file.path(tab_dir,
                    "main_result_no_covariates_vs_with_covariates.csv"),
          row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(res_dir, "session_info.txt"))
