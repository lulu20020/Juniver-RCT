# -----------------------------------------------------------------------------
# Script: missingness_diagnostics.R
# Purpose:
#   - Run missingness diagnostics on the long-format analysis dataset using a
#     baseline-complete subsample.
#   - Models:
#       * Model A1: I(wk4 missing)  ~ GROUP + baseline EDEQglobal + baseline BMI +
#                   baseline age + gender + baseline readythoughts + baseline readybehs
#       * Model A2: I(wk8 missing)  ~ same as A1
#       * Model A3: I(wk12 missing) ~ same as A1
#       * Model B1: I(wk12 missing) ~ GROUP + baseline EDEQglobal + I(wk4 missing) +
#                   I(wk8 missing) + baseline BMI + baseline age + gender +
#                   baseline readythoughts + baseline readybehs +
#                   country + ethnicity_whitevsother + baseline PHQ9/GAD7/CIA/PSS
#   - Exploratory (post wk12):
#       * Model E (wk24): I(wk24 missing) ~ treatment_weeks (Juniver=24, Control=12)
#         + early missing (wk4/wk8/wk12) + same baseline covariates as B1.
#         Interpretation: exploratory – reflects association between cumulative
#         weeks on treatment and wk24 missingness (early vs delayed treatment).
#   - Sample: keep IDs with no missing at baseline on key vars
#             (EDEQglobal, age, gender, BMI, country, ethnicity_whitevsother,
#              readythoughts, readybehs, PHQ9, GAD7, CIA, PSS),
#             then use all visits for those IDs.
#   - Key outputs:
#       * OR tables (csv): model_wk4_missing / model_wk8_missing / model_wk12_missing
#       * Visit-level retention table:
#         randomised N, observed N, missing N, missing %, retention %
# Notes:
#   - Input data are not included in the public repository.
#   - No imputation; baseline-complete IDs only.
#   - To change required baseline vars or covariates, edit required_bl and formulas.
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(forcats)
library(broom)
library(readxl)

input_path <- file.path("data", "analysis_dataset.csv")
results_dir <- file.path("outputs", "missingness_diagnostics")
analysis_varlist_path <- file.path("data", "analysis_variable_list.xlsx")
analysis_missing_dir <- file.path("outputs", "missingness_summary")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(analysis_missing_dir)) dir.create(analysis_missing_dir, recursive = TRUE)

stopifnot(file.exists(input_path))

# --- Load and clean ------------------------------------------------------
dat <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
dat$visit <- tolower(dat$visit)

# Outcome variable name
edeq_var <- "edeqglobal"

# --- Additional missingness output using Analysis-sheet variables --------
read_analysis_vars <- function(path, sheet = "Analysis") {
  if (!file.exists(path)) return(character(0))
  x <- readxl::read_xlsx(path, sheet = sheet)
  vals <- unlist(x, use.names = FALSE)
  vals <- trimws(as.character(vals))
  vals <- vals[!is.na(vals) & vals != ""]
  unique(vals)
}

analysis_vars <- read_analysis_vars(analysis_varlist_path, "Analysis")
analysis_vars <- intersect(analysis_vars, names(dat))

if (length(analysis_vars) > 0) {
  dat_miss <- dat
  dat_miss$group <- as.character(dat_miss$group)
  dat_miss$group[dat_miss$group == "1"] <- "Juniver"
  dat_miss$group[dat_miss$group == "2"] <- "Control"
  dat_miss$group <- factor(dat_miss$group, levels = c("Juniver", "Control"))

  groups_avail <- levels(dat_miss$group)
  groups_avail <- groups_avail[!is.na(groups_avail)]

  # Visit x variable missingness (overall + by group)
  visit_levels_all <- c("bl", "wk4", "wk8", "wk12", "wk16", "wk20", "wk24")
  dat_miss$visit <- factor(dat_miss$visit, levels = visit_levels_all, ordered = TRUE)

  # Long-format cell table (id x visit x variable), robust to duplicate records.
  cell_dat <- dat_miss %>%
    select(id, visit, group, all_of(analysis_vars)) %>%
    pivot_longer(cols = all_of(analysis_vars), names_to = "variable", values_to = "value") %>%
    group_by(id, visit, group, variable) %>%
    summarise(
      value = {
        nv <- value[!is.na(value)]
        if (length(nv) == 0) NA else nv[1]
      },
      .groups = "drop"
    )

  # Variable-level missingness across all id-visit cells (long-format denominator).
  miss_sum <- cell_dat %>%
    group_by(variable) %>%
    summarise(
      total_cells_total = n(),
      unique_ids_total = dplyr::n_distinct(id),
      visits_total = dplyr::n_distinct(visit),
      missing_count_total = sum(is.na(value)),
      observed_count_total = sum(!is.na(value)),
      missing_rate_total = round(mean(is.na(value)), 4),
      .groups = "drop"
    )

  for (g in groups_avail) {
    tmp_g <- cell_dat %>%
      filter(group == g) %>%
      group_by(variable) %>%
      summarise(
        total_cells = n(),
        unique_ids = dplyr::n_distinct(id),
        visits_n = dplyr::n_distinct(visit),
        missing_count = sum(is.na(value)),
        observed_count = sum(!is.na(value)),
        missing_rate = round(mean(is.na(value)), 4),
        .groups = "drop"
      )
    miss_sum <- miss_sum %>%
      left_join(
        tmp_g %>% rename(
          !!sprintf("total_cells_%s", g) := total_cells,
          !!sprintf("unique_ids_%s", g) := unique_ids,
          !!sprintf("visits_%s", g) := visits_n,
          !!sprintf("missing_count_%s", g) := missing_count,
          !!sprintf("observed_count_%s", g) := observed_count,
          !!sprintf("missing_rate_%s", g) := missing_rate
        ),
        by = "variable"
      )
  }

  write.csv(
    miss_sum,
    file.path(analysis_missing_dir, "analysis_missing_summary.csv"),
    row.names = FALSE
  )

  miss_by_visit <- cell_dat %>%
    group_by(variable, visit) %>%
    summarise(
      total_cells_total = n(),
      unique_ids_total = dplyr::n_distinct(id),
      missing_count_total = sum(is.na(value)),
      observed_count_total = sum(!is.na(value)),
      missing_rate_total = round(mean(is.na(value)), 4),
      .groups = "drop"
    )

  for (g in groups_avail) {
    tmp <- cell_dat %>%
      filter(group == g) %>%
      group_by(variable, visit) %>%
      summarise(
        n_cells = n(),
        n_ids = dplyr::n_distinct(id),
        mc = sum(is.na(value)),
        obs = sum(!is.na(value)),
        mr = round(mean(is.na(value)), 4),
        .groups = "drop"
      )

    miss_by_visit <- miss_by_visit %>%
      left_join(
        tmp %>%
          rename(
            !!sprintf("total_cells_%s", g) := n_cells,
            !!sprintf("unique_ids_%s", g) := n_ids,
            !!sprintf("missing_count_%s", g) := mc,
            !!sprintf("observed_count_%s", g) := obs,
            !!sprintf("missing_rate_%s", g) := mr
          ),
        by = c("variable", "visit")
      )
  }

  write.csv(
    miss_by_visit,
    file.path(analysis_missing_dir, "analysis_missing_by_visit.csv"),
    row.names = FALSE
  )

  # Visit-level summary over all Analysis variables
  visit_summary <- do.call(
    rbind,
    lapply(c("total", groups_avail), function(glabel) {
      df <- if (glabel == "total") {
        cell_dat
      } else {
        cell_dat[cell_dat$group == glabel, , drop = FALSE]
      }
      do.call(
        rbind,
        lapply(levels(cell_dat$visit), function(v) {
          if (is.na(v)) return(NULL)
          sub <- df[df$visit == v, , drop = FALSE]
          if (nrow(sub) == 0) return(NULL)
          total_cells <- nrow(sub)
          missing_cells <- sum(is.na(sub$value))
          id_miss <- sub %>%
            group_by(id) %>%
            summarise(miss_prop = mean(is.na(value)), .groups = "drop")
          data.frame(
            group = glabel,
            visit = v,
            total_rows = dplyr::n_distinct(sub$id),
            unique_ids = dplyr::n_distinct(sub$id),
            variables_n = dplyr::n_distinct(sub$variable),
            total_cells = total_cells,
            missing_cells = missing_cells,
            observed_cells = total_cells - missing_cells,
            missing_rate_cells = round(missing_cells / total_cells, 4),
            mean_missing_prop_per_id = round(mean(id_miss$miss_prop), 4),
            median_missing_prop_per_id = round(stats::median(id_miss$miss_prop), 4),
            complete_id_rate = round(mean(id_miss$miss_prop == 0), 4),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  write.csv(
    visit_summary,
    file.path(analysis_missing_dir, "analysis_missing_by_visit_summary.csv"),
    row.names = FALSE
  )
}

# Baseline required variables
required_bl <- c(
  edeq_var, "age", "gender", "bmi", "country", "ethnicity_whitevsother",
  "readythoughts", "readybehs",
  "phq9total", "gad7total", "cia", "pss"
)

# Keep IDs with complete baseline on required vars
dat_bl <- dat %>% filter(visit == "bl")
complete_ids <- dat_bl %>%
  filter(if_all(all_of(required_bl), ~ !is.na(.x))) %>%
  pull(id)

dat_cc <- dat %>% filter(id %in% complete_ids)

# Missing indicators per visit for target
miss_flags <- dat_cc %>%
  select(id, visit, all_of(edeq_var)) %>%
  mutate(miss = is.na(.data[[edeq_var]])) %>%
  select(-all_of(edeq_var)) %>%
  pivot_wider(names_from = visit, values_from = miss, names_prefix = "miss_")

# Assemble analysis baseline
baseline <- dat_cc %>% filter(visit == "bl") %>%
  select(id, group, age, gender, bmi_bl = bmi, country, ethnicity_whitevsother,
         edeq_bl = all_of(edeq_var),
         readythoughts_bl = readythoughts,
         readybehs_bl = readybehs,
         phq9_bl = phq9total,
         gad7_bl = gad7total,
         cia_bl = cia,
         pss_bl = pss)

# Merge
analysis <- baseline %>% left_join(miss_flags, by = "id")

# Factorize categorical vars
analysis <- analysis %>%
  mutate(
    group = trimws(as.character(group)),
    group = dplyr::case_when(
      tolower(group) %in% c("1", "juniver") ~ "Juniver",
      tolower(group) %in% c("2", "control") ~ "Control",
      TRUE ~ group
    ),
    group = factor(group, levels = c("Juniver", "Control")),
    age = factor(age),
    gender = factor(gender),
    country = factor(country),
    ethnicity_whitevsother = factor(ethnicity_whitevsother)
  )

# Relevel gender so female ("1") is the reference, enabling
# female-vs-male and female-vs-non-binary reporting.
if ("1" %in% levels(analysis$gender)) {
  analysis$gender <- forcats::fct_relevel(analysis$gender, "1")
}

# Ensure missing flags exist
for (v in c("miss_wk4", "miss_wk8", "miss_wk12", "miss_wk16", "miss_wk20", "miss_wk24")) {
  if (!v %in% names(analysis)) analysis[[v]] <- NA
}

# Cumulative treatment weeks at wk24 (exploratory): Juniver ~24, Control ~12
analysis <- analysis %>%
  mutate(treatment_weeks_wk24 = ifelse(group == "Juniver", 24, 12))

# --- Models --------------------------------------------------------------
fit_wk4 <- glm(miss_wk4 ~ group + edeq_bl + age + gender + bmi_bl +
                 readythoughts_bl + readybehs_bl,
               data = analysis, family = binomial)
fit_wk8 <- glm(miss_wk8 ~ group + edeq_bl + age + gender + bmi_bl +
                 readythoughts_bl + readybehs_bl,
               data = analysis, family = binomial)
fit_wk4  <- glm(miss_wk4  ~ group + edeq_bl + age + gender + bmi_bl +
                  readythoughts_bl + readybehs_bl,
                data = analysis, family = binomial)
fit_wk8  <- glm(miss_wk8  ~ group + edeq_bl + age + gender + bmi_bl +
                  readythoughts_bl + readybehs_bl,
                data = analysis, family = binomial)
fit_wk12A <- glm(miss_wk12 ~ group + edeq_bl + age + gender + bmi_bl +
                   readythoughts_bl + readybehs_bl,
                 data = analysis, family = binomial)

fit_wk12B <- glm(miss_wk12 ~ group + edeq_bl + miss_wk4 + miss_wk8 +
                   age + gender + bmi_bl + country + ethnicity_whitevsother +
                   readythoughts_bl + readybehs_bl +
                   phq9_bl + gad7_bl + cia_bl + pss_bl,
                 data = analysis, family = binomial)

# Exploratory post-wk12 model (wk24)
fit_wk24E <- glm(miss_wk24 ~ treatment_weeks_wk24 +
                   miss_wk4 + miss_wk8 + miss_wk12 +
                   age + gender + bmi_bl + country + ethnicity_whitevsother +
                   readythoughts_bl + readybehs_bl +
                   phq9_bl + gad7_bl + cia_bl + pss_bl,
                 data = analysis, family = binomial)

# Helper to tidy odds ratios
tidy_or <- function(fit) {
  tab <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    transmute(
      term,
      OR = estimate,
      conf.low,
      conf.high,
      z.value = statistic,
      p.value
    )

  # Helper for reciprocal contrast (A vs B from B vs A).
  invert_row <- function(df_row, new_term) {
    df_row %>%
      mutate(
        OR_old = OR,
        low_old = conf.low,
        high_old = conf.high,
        z_old = z.value
      ) %>%
      transmute(
        term = new_term,
        OR = 1 / OR_old,
        conf.low = 1 / high_old,
        conf.high = 1 / low_old,
        z.value = -z_old,
        p.value = p.value
      )
  }

  # Gender: report female vs male and female vs non-binary.
  # With female ("1") as reference, glm returns:
  #   gender2 = male vs female
  #   gender0 = non-binary vs female
  if ("gender2" %in% tab$term) {
    tab <- bind_rows(
      tab %>% filter(term != "gender2"),
      invert_row(tab %>% filter(term == "gender2"), "gender_female_vs_male")
    )
  }
  if ("gender0" %in% tab$term) {
    tab <- bind_rows(
      tab %>% filter(term != "gender0"),
      invert_row(
        tab %>% filter(term == "gender0"),
        "gender_female_vs_non_binary"
      )
    )
  }

  # Ethnicity: report white vs other.
  # With level "1" as white and "2" as other, glm returns:
  #   ethnicity_whitevsother2 = other vs white
  if ("ethnicity_whitevsother2" %in% tab$term) {
    tab <- bind_rows(
      tab %>% filter(term != "ethnicity_whitevsother2"),
      invert_row(
        tab %>% filter(term == "ethnicity_whitevsother2"),
        "ethnicity_white_vs_other"
      )
    )
  }

  # Generic explicit labels for remaining categorical contrasts:
  # <variable>_<level>_vs_<reference>.
  sanitize <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[^a-z0-9]+", "_", x)
    gsub("^_|_$", "", x)
  }

  for (v in names(fit$xlevels)) {
    lv <- fit$xlevels[[v]]
    if (length(lv) < 2) next
    ref <- lv[[1]]
    for (k in lv[-1]) {
      raw_term <- paste0(v, k)
      new_term <- paste0(v, "_", sanitize(k), "_vs_", sanitize(ref))
      tab$term[tab$term == raw_term] <- new_term
    }
  }

  tab
}

write.csv(tidy_or(fit_wk4),   file.path(results_dir, "model_wk4_missing.csv"),    row.names = FALSE)
write.csv(tidy_or(fit_wk8),   file.path(results_dir, "model_wk8_missing.csv"),    row.names = FALSE)
write.csv(tidy_or(fit_wk12A), file.path(results_dir, "model_wk12_missing_A.csv"), row.names = FALSE)
write.csv(tidy_or(fit_wk12B), file.path(results_dir, "model_wk12_missing_B.csv"), row.names = FALSE)
write.csv(tidy_or(fit_wk24E), file.path(results_dir, "model_wk24_missing_exploratory.csv"), row.names = FALSE)

# --- Visit-level randomised / missing / retention table ------------------
visits_all <- c("bl", "wk4", "wk8", "wk12", "wk16", "wk20", "wk24")
visit_labels <- c(
  bl = "Baseline",
  wk4 = "Week 4",
  wk8 = "Week 8",
  wk12 = "Week 12",
  wk16 = "Week 16",
  wk20 = "Week 20",
  wk24 = "Week 24"
)

dat_ret <- dat %>%
  mutate(
    group_std = case_when(
      tolower(trimws(as.character(group))) %in% c("1", "juniver") ~ "Juniver",
      tolower(trimws(as.character(group))) %in% c("2", "control") ~ "Control",
      TRUE ~ as.character(group)
    ),
    visit = tolower(as.character(visit)),
    outcome_value = .data[[edeq_var]]
  ) %>%
  filter(!is.na(id), visit %in% visits_all)

# One row per id x visit; if duplicates exist, keep one observed value if present.
cell_ret <- dat_ret %>%
  group_by(id, visit) %>%
  summarise(
    group_std = {
      g <- group_std[!is.na(group_std)]
      if (length(g) == 0) NA_character_ else g[[1]]
    },
    outcome_value = {
      x <- outcome_value[!is.na(outcome_value)]
      if (length(x) == 0) NA_real_ else as.numeric(x[[1]])
    },
    .groups = "drop"
  )

id_group <- dat_ret %>%
  group_by(id) %>%
  summarise(
    group_std = {
      g <- group_std[!is.na(group_std)]
      if (length(g) == 0) NA_character_ else g[[1]]
    },
    .groups = "drop"
  ) %>%
  filter(group_std %in% c("Juniver", "Control"))

group_levels <- c("Total", "Juniver", "Control")

build_row <- function(v, gname) {
  ids_use <- if (gname == "Total") {
    id_group$id
  } else {
    id_group$id[id_group$group_std == gname]
  }
  n_rand <- length(ids_use)
  sub <- cell_ret %>% filter(id %in% ids_use, visit == v)
  n_obs <- dplyr::n_distinct(sub$id[!is.na(sub$outcome_value)])
  n_miss <- n_rand - n_obs
  miss_pct <- if (n_rand > 0) round(100 * n_miss / n_rand, 1) else NA_real_
  ret_pct <- if (n_rand > 0) round(100 * n_obs / n_rand, 1) else NA_real_

  data.frame(
    Visit = visit_labels[[v]],
    Group = gname,
    check.names = FALSE,
    `Randomised, N` = n_rand,
    `Observed, No` = n_obs,
    `Missing, No.` = n_miss,
    `Missing, %` = miss_pct,
    `Retention, %` = ret_pct
  )
}

retention_tbl <- do.call(
  rbind,
  lapply(visits_all, function(v) {
    do.call(rbind, lapply(group_levels, function(g) build_row(v, g)))
  })
)

write.csv(
  retention_tbl,
  file.path(results_dir, "visit_randomised_missing_retention.csv"),
  row.names = FALSE
)

cat("Done. Results in:", results_dir, "\n")
