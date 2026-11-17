# ============================================================
# 04_summary_tables_simulation_updated.R
# 功能：
# 1. 彙整 GEE / LME / DL 之 simulation 結果
# 2. 依 Overall / Control / Statin 分別計算 RMSE、MAE、MAPE、Bias 平均與 SD
# 3. 輸出 CSV 與 Excel 工作簿
# ============================================================

read_master_safely <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

standardize_sim_result_cols <- function(dat) {

  if (nrow(dat) == 0) return(dat)

  dat <- dat %>%
    dplyr::mutate(
      outcome = dplyr::recode(
        as.character(outcome),
        HbA1C = "HbA1c",
        Glucose_AC = "Glu_AC",
        Triglycerides = "TG",
        .default = as.character(outcome)
      ),
      model = as.character(model)
    )

  if (!("outcome_label" %in% names(dat))) {
    dat <- dat %>% dplyr::mutate(outcome_label = unname(SIM_OUTCOME_LABELS[outcome]))
  } else {
    dat <- dat %>% dplyr::mutate(outcome_label = as.character(outcome_label))
  }

  if (!("Group_eval" %in% names(dat))) {
    dat <- dat %>% dplyr::mutate(Group_eval = "Overall")
  } else {
    dat <- dat %>% dplyr::mutate(Group_eval = as.character(Group_eval))
  }

  if (!("Group_label" %in% names(dat))) {
    dat <- dat %>% dplyr::mutate(Group_label = unname(SIM_GROUP_LABELS[Group_eval]))
  } else {
    dat <- dat %>% dplyr::mutate(Group_label = as.character(Group_label))
  }

  dat %>%
    dplyr::mutate(
      outcome = factor(outcome, levels = SIM_OUTCOMES),
      outcome_label = factor(outcome_label, levels = SIM_OUTCOME_LABELS[SIM_OUTCOMES]),
      Group_eval = factor(Group_eval, levels = c("Overall", "Control", "Statin")),
      Group_label = factor(Group_label, levels = unname(SIM_GROUP_LABELS[c("Overall", "Control", "Statin")]))
    )
}

summarise_one_metric_set <- function(dat) {

  dat %>%
    dplyr::group_by(scenario, outcome, outcome_label, model, Group_eval, Group_label) %>%
    dplyr::summarise(
      RMSE_mean = mean(RMSE, na.rm = TRUE),
      RMSE_sd = stats::sd(RMSE, na.rm = TRUE),
      MAE_mean = mean(MAE, na.rm = TRUE),
      MAE_sd = stats::sd(MAE, na.rm = TRUE),
      MAPE_mean = mean(MAPE, na.rm = TRUE),
      MAPE_sd = stats::sd(MAPE, na.rm = TRUE),
      Bias_mean = mean(Bias, na.rm = TRUE),
      Bias_sd = stats::sd(Bias, na.rm = TRUE),
      n_test_mean = mean(n_test, na.rm = TRUE),
      n_test_sd = stats::sd(n_test, na.rm = TRUE),
      n_run_success = dplyr::n_distinct(
        run[
          is.finite(RMSE) &
            is.finite(MAE) &
            is.finite(MAPE) &
            is.finite(Bias)
        ]
      ),
      n_run_total = dplyr::n_distinct(run),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_run_failed = n_run_total - n_run_success,
      success_rate = n_run_success / n_run_total,
      RMSE_MCSE = RMSE_sd / sqrt(n_run_success),
      MAE_MCSE = MAE_sd / sqrt(n_run_success),
      MAPE_MCSE = MAPE_sd / sqrt(n_run_success),
      Bias_MCSE = Bias_sd / sqrt(n_run_success)
    ) %>%
    dplyr::arrange(
      scenario,
      factor(outcome, levels = SIM_OUTCOMES),
      factor(model, levels = c("GEE", "LME", "LSTM-NARX", "TCN", "ODE-LSTM")),
      Group_eval
    )
}


make_paper_metric_model_wide <- function(summary_df, group_filter = "Overall") {

  if (nrow(summary_df) == 0) return(tibble::tibble())

  dat <- summary_df %>%
    dplyr::filter(as.character(Group_eval) == group_filter) %>%
    dplyr::mutate(
      RMSE_fmt = sprintf("%.3f ± %.3f", RMSE_mean, RMSE_sd),
      MAE_fmt = sprintf("%.3f ± %.3f", MAE_mean, MAE_sd),
      MAPE_fmt = sprintf("%.3f ± %.3f", MAPE_mean, MAPE_sd),
      Bias_fmt = sprintf("%.3f ± %.3f", Bias_mean, Bias_sd)
    ) %>%
    dplyr::select(
      scenario, Group_eval, Group_label, model, outcome, outcome_label,
      RMSE_fmt, MAE_fmt, MAPE_fmt, Bias_fmt
    ) %>%
    tidyr::pivot_longer(
      cols = c(RMSE_fmt, MAE_fmt, MAPE_fmt, Bias_fmt),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        RMSE_fmt = "RMSE",
        MAE_fmt = "MAE",
        MAPE_fmt = "MAPE",
        Bias_fmt = "Bias"
      ),
      metric = factor(metric, levels = c("RMSE", "MAE", "MAPE", "Bias")),
      model = factor(model, levels = c("GEE", "LME", "LSTM-NARX", "TCN", "ODE-LSTM")),
      outcome_label = factor(outcome_label, levels = SIM_OUTCOME_LABELS[SIM_OUTCOMES])
    ) %>%
    dplyr::arrange(metric, model)

  dat %>%
    dplyr::select(metric, model, outcome_label, value) %>%
    tidyr::pivot_wider(
      names_from = outcome_label,
      values_from = value
    ) %>%
    dplyr::arrange(metric, model)
}

make_best_model_summary <- function(summary_df) {

  if (nrow(summary_df) == 0) return(list(by_metric = tibble::tibble(), count = tibble::tibble()))

  metric_long <- summary_df %>%
    dplyr::select(
      scenario, outcome, outcome_label, model, Group_eval, Group_label,
      RMSE_mean, MAE_mean, MAPE_mean, Bias_mean
    ) %>%
    tidyr::pivot_longer(
      cols = c(RMSE_mean, MAE_mean, MAPE_mean, Bias_mean),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        RMSE_mean = "RMSE",
        MAE_mean = "MAE",
        MAPE_mean = "MAPE",
        Bias_mean = "Bias"
      ),
      criterion_value = dplyr::if_else(metric == "Bias", abs(value), value)
    )

  by_metric <- metric_long %>%
    dplyr::filter(is.finite(criterion_value)) %>%
    dplyr::group_by(scenario, Group_eval, Group_label, outcome, outcome_label, metric) %>%
    dplyr::slice_min(order_by = criterion_value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      scenario,
      Group_eval,
      Group_label,
      outcome,
      outcome_label,
      metric,
      best_model = as.character(model),
      best_value = value
    ) %>%
    dplyr::arrange(scenario, Group_eval, outcome, factor(metric, levels = c("RMSE", "MAE", "MAPE", "Bias")))

  count <- by_metric %>%
    dplyr::filter(metric %in% c("RMSE", "MAE", "MAPE")) %>%
    dplyr::count(scenario, Group_eval, Group_label, outcome, outcome_label, best_model, name = "n_best_error_metrics") %>%
    dplyr::group_by(scenario, Group_eval, Group_label, outcome, outcome_label) %>%
    dplyr::arrange(dplyr::desc(n_best_error_metrics), best_model, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(overall_best_model_by_error_metrics = best_model)

  list(by_metric = by_metric, count = count)
}

summarise_simulation_outputs <- function(
    out_dir,
    scenario_name = NULL
) {

  lme_path <- file.path(out_dir, "lme_results", "lme_results_master.csv")
  gee_path <- file.path(out_dir, "gee_results", "gee_results_master.csv")
  dl_path <- file.path(out_dir, "dl_results", "dl_results_master.csv")

  lme_res <- read_master_safely(lme_path)
  gee_res <- read_master_safely(gee_path)
  dl_res <- read_master_safely(dl_path)

  lme_res <- standardize_sim_result_cols(lme_res)
  gee_res <- standardize_sim_result_cols(gee_res)
  dl_res <- standardize_sim_result_cols(dl_res)

  stat_pred <- dplyr::bind_rows(gee_res, lme_res)
  dl_pred <- dl_res

  summary_stat_pred <- if (nrow(stat_pred) > 0) {
    summarise_one_metric_set(stat_pred)
  } else {
    tibble::tibble()
  }

  summary_dl <- if (nrow(dl_pred) > 0) {
    summarise_one_metric_set(dl_pred)
  } else {
    tibble::tibble()
  }

  all_model_summary <- dplyr::bind_rows(
    summary_stat_pred,
    summary_dl
  ) %>%
    dplyr::arrange(
      scenario,
      factor(outcome, levels = SIM_OUTCOMES),
      factor(model, levels = c("GEE", "LME", "LSTM-NARX", "TCN", "ODE-LSTM")),
      Group_eval
    )

  mcse_summary <- all_model_summary %>%
    dplyr::summarise(
      max_RMSE_MCSE = max(RMSE_MCSE, na.rm = TRUE),
      max_MAE_MCSE = max(MAE_MCSE, na.rm = TRUE),
      max_MAPE_MCSE = max(MAPE_MCSE, na.rm = TRUE),
      max_Bias_MCSE = max(Bias_MCSE, na.rm = TRUE),
      minimum_success_rate = min(success_rate, na.rm = TRUE)
    )

  summary_dir <- file.path(out_dir, "summary_tables")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- ifelse(is.null(scenario_name), "simulation", scenario_name)

  readr::write_csv(
    summary_stat_pred,
    file.path(summary_dir, paste0(prefix, "_GEE_LME_summary_group.csv"))
  )

  if (nrow(summary_dl) > 0) {
    readr::write_csv(
      summary_dl,
      file.path(summary_dir, paste0(prefix, "_DL_summary_group.csv"))
    )
  }

  readr::write_csv(
    all_model_summary,
    file.path(summary_dir, paste0(prefix, "_five_model_summary_group.csv"))
  )

  readr::write_csv(
    mcse_summary,
    file.path(summary_dir, paste0(prefix, "_MCSE_and_success_summary.csv"))
  )

  # 論文主表常用寬格式：每個模型與組別一列，各指標以 mean (SD) 呈現。
  table_for_paper <- all_model_summary %>%
    dplyr::mutate(
      RMSE_fmt = sprintf("%.3f (%.3f)", RMSE_mean, RMSE_sd),
      MAE_fmt = sprintf("%.3f (%.3f)", MAE_mean, MAE_sd),
      MAPE_fmt = sprintf("%.3f (%.3f)", MAPE_mean, MAPE_sd),
      Bias_fmt = sprintf("%.3f (%.3f)", Bias_mean, Bias_sd)
    ) %>%
    dplyr::select(
      scenario, outcome, outcome_label, model, Group_eval, Group_label,
      RMSE_fmt, MAE_fmt, MAPE_fmt, Bias_fmt,
      n_test_mean,
      n_run_success, n_run_failed, n_run_total, success_rate,
      RMSE_MCSE, MAE_MCSE, MAPE_MCSE, Bias_MCSE
    )

  readr::write_csv(
    table_for_paper,
    file.path(summary_dir, paste0(prefix, "_five_model_summary_group_for_paper.csv"))
  )


  paper_table_overall <- make_paper_metric_model_wide(all_model_summary, "Overall")
  paper_table_control <- make_paper_metric_model_wide(all_model_summary, "Control")
  paper_table_statin <- make_paper_metric_model_wide(all_model_summary, "Statin")

  readr::write_csv(
    paper_table_overall,
    file.path(summary_dir, paste0(prefix, "_paper_table_Overall_wide.csv"))
  )
  readr::write_csv(
    paper_table_control,
    file.path(summary_dir, paste0(prefix, "_paper_table_Control_wide.csv"))
  )
  readr::write_csv(
    paper_table_statin,
    file.path(summary_dir, paste0(prefix, "_paper_table_Statin_wide.csv"))
  )

  best_model_summary <- make_best_model_summary(all_model_summary)

  readr::write_csv(
    best_model_summary$by_metric,
    file.path(summary_dir, paste0(prefix, "_best_model_by_metric.csv"))
  )
  readr::write_csv(
    best_model_summary$count,
    file.path(summary_dir, paste0(prefix, "_best_model_summary_for_paper.csv"))
  )

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(
        GEE_LME_summary_group = summary_stat_pred,
        DL_summary_group = summary_dl,
        five_model_summary_group = all_model_summary,
        MCSE_success_summary = mcse_summary,
        table_for_paper = table_for_paper,
        paper_table_Overall = paper_table_overall,
        paper_table_Control = paper_table_control,
        paper_table_Statin = paper_table_statin,
        best_model_by_metric = best_model_summary$by_metric,
        best_model_summary = best_model_summary$count
      ),
      file.path(summary_dir, paste0(prefix, "_simulation_summary_tables.xlsx"))
    )
  }

  list(
    summary_stat_pred = summary_stat_pred,
    summary_dl = summary_dl,
    all_model_summary = all_model_summary,
    mcse_summary = mcse_summary,
    table_for_paper = table_for_paper,
    paper_table_overall = paper_table_overall,
    paper_table_control = paper_table_control,
    paper_table_statin = paper_table_statin,
    best_model_by_metric = best_model_summary$by_metric,
    best_model_summary = best_model_summary$count
  )
}

consolidate_all_scenarios <- function(sim_root, scenarios) {

  summary_files <- file.path(
    sim_root,
    scenarios,
    "summary_tables",
    paste0(scenarios, "_five_model_summary_group.csv")
  )

  existing <- summary_files[file.exists(summary_files)]

  if (length(existing) == 0) {
    warning("尚未找到任何情境的 five_model_summary_group.csv，略過全情境彙整。")
    return(invisible(NULL))
  }

  all_summary <- purrr::map_dfr(existing, readr::read_csv, show_col_types = FALSE) %>%
    standardize_sim_result_cols()

  out_dir <- file.path(sim_root, "summary_all_scenarios")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(
    all_summary,
    file.path(out_dir, "all_scenarios_five_model_summary_group.csv")
  )

  paper_long <- all_summary %>%
    dplyr::mutate(
      RMSE_fmt = sprintf("%.3f ± %.3f", RMSE_mean, RMSE_sd),
      MAE_fmt = sprintf("%.3f ± %.3f", MAE_mean, MAE_sd),
      MAPE_fmt = sprintf("%.3f ± %.3f", MAPE_mean, MAPE_sd),
      Bias_fmt = sprintf("%.3f ± %.3f", Bias_mean, Bias_sd)
    ) %>%
    dplyr::select(
      scenario, outcome, outcome_label, model, Group_eval, Group_label,
      RMSE_fmt, MAE_fmt, MAPE_fmt, Bias_fmt,
      n_test_mean,
      n_run_success, n_run_failed, n_run_total, success_rate,
      RMSE_MCSE, MAE_MCSE, MAPE_MCSE, Bias_MCSE
    )

  readr::write_csv(
    paper_long,
    file.path(out_dir, "all_scenarios_table_for_paper_long.csv")
  )

  # 依論文表 5.3--5.7：預設輸出 Overall 寬表；Control / Statin 另存供補充或新版老師要求使用。
  for (sc in unique(as.character(all_summary$scenario))) {
    sc_df <- all_summary %>% dplyr::filter(as.character(scenario) == sc)
    for (grp in c("Overall", "Control", "Statin")) {
      wide <- make_paper_metric_model_wide(sc_df, grp)
      readr::write_csv(
        wide,
        file.path(out_dir, paste0(sc, "_paper_table_", grp, "_wide.csv"))
      )
    }
  }

  best <- make_best_model_summary(all_summary)

  readr::write_csv(
    best$by_metric,
    file.path(out_dir, "all_scenarios_best_model_by_metric.csv")
  )
  readr::write_csv(
    best$count,
    file.path(out_dir, "all_scenarios_best_model_summary_for_table_5_8.csv")
  )

  if (requireNamespace("writexl", quietly = TRUE)) {
    sheets <- list(
      all_summary = all_summary,
      table_for_paper_long = paper_long,
      best_model_by_metric = best$by_metric,
      best_model_summary = best$count
    )

    for (sc in unique(as.character(all_summary$scenario))) {
      sc_df <- all_summary %>% dplyr::filter(as.character(scenario) == sc)
      sheets[[paste0(sc, "_Overall")]] <- make_paper_metric_model_wide(sc_df, "Overall")
      sheets[[paste0(sc, "_Control")]] <- make_paper_metric_model_wide(sc_df, "Control")
      sheets[[paste0(sc, "_Statin")]] <- make_paper_metric_model_wide(sc_df, "Statin")
    }

    writexl::write_xlsx(
      sheets,
      file.path(out_dir, "all_scenarios_simulation_summary_tables.xlsx")
    )
  }

  invisible(all_summary)
}
