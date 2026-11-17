# ============================================================
# 03_fit_lme_gee_simulation_updated.R
# 功能：
# 1. 以 train set 配適 GEE 與 LME-RI+RS
# 2. 以 test set 每位病患最後一次可評估觀測值作為預測目標
# 3. 於原始尺度計算 RMSE、MAE、MAPE、Bias
# 4. 指標分 Overall / Control / Statin 呈現
# ============================================================


calculate_sim_smearing_factor <- function(observed_y, fitted_log, label) {
  keep <- (
    is.finite(observed_y) &
      observed_y > 0 &
      is.finite(fitted_log)
  )

  if (!any(keep)) {
    stop("無法估計 empirical smearing factor：", label, call. = FALSE)
  }

  sf <- mean(
    exp(log(observed_y[keep]) - fitted_log[keep]),
    na.rm = TRUE
  )

  if (!is.finite(sf) || sf <= 0) {
    stop("empirical smearing factor 無效：", label, call. = FALSE)
  }

  sf
}

calc_pred_metrics <- function(y_true, y_pred) {

  keep <- is.finite(y_true) & is.finite(y_pred) & y_true > 0

  y_true <- y_true[keep]
  y_pred <- y_pred[keep]

  if (length(y_true) == 0) {
    return(tibble::tibble(
      RMSE = NA_real_,
      MAE = NA_real_,
      MAPE = NA_real_,
      Bias = NA_real_,
      n_test = 0L
    ))
  }

  error <- y_pred - y_true

  tibble::tibble(
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    MAPE = mean(abs(error / y_true)) * 100,
    Bias = mean(error),
    n_test = length(y_true)
  )
}

calc_metrics_by_group <- function(pred_df) {

  pred_overall <- pred_df %>%
    dplyr::mutate(Group_eval = "Overall")

  pred_group <- pred_df %>%
    dplyr::mutate(Group_eval = as.character(Group))

  dplyr::bind_rows(pred_overall, pred_group) %>%
    dplyr::group_by(scenario, run, outcome, outcome_label, model, model_detail, Group_eval) %>%
    dplyr::summarise(
      calc_pred_metrics(y_true = y_true, y_pred = y_pred),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Group_eval = factor(Group_eval, levels = c("Overall", "Control", "Statin")),
      Group_label = SIM_GROUP_LABELS[as.character(Group_eval)]
    ) %>%
    dplyr::arrange(outcome, model, Group_eval)
}

build_sim_model_formula <- function(outcome_name) {
  stats::as.formula(
    paste0(
      "log(", outcome_name, ") ~ ",
      "GroupStatin * Time_in_Years + ",
      "Age_baseline + Sex1 + BMI"
    )
  )
}

prepare_sim_model_data <- function(sim_data, outcome_name) {

  sim_data %>%
    dplyr::mutate(
      ID = as.factor(ID),
      Group = factor(Group, levels = c("Control", "Statin")),
      split = factor(split, levels = c("train", "validation", "test")),
      GroupStatin = as.integer(Group == "Statin"),
      Sex1 = as.integer(Sex1)
    ) %>%
    dplyr::filter(
      is.finite(.data[[outcome_name]]),
      .data[[outcome_name]] > 0,
      is.finite(Time_in_Years),
      is.finite(Age_baseline),
      is.finite(BMI),
      !is.na(Sex1),
      !is.na(split)
    ) %>%
    dplyr::arrange(ID, Time_in_Years)
}

get_test_last_observation <- function(dat) {

  dat %>%
    dplyr::filter(split == "test") %>%
    dplyr::group_by(ID) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::slice_max(Time_in_Years, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

make_failed_metric_rows <- function(scenario, run, outcome_name, model, model_detail) {

  tibble::tibble(
    scenario = scenario,
    run = run,
    outcome = outcome_name,
    outcome_label = SIM_OUTCOME_LABELS[[outcome_name]],
    model = model,
    model_detail = model_detail,
    Group_eval = factor(c("Overall", "Control", "Statin"), levels = c("Overall", "Control", "Statin")),
    Group_label = SIM_GROUP_LABELS[as.character(Group_eval)],
    RMSE = NA_real_,
    MAE = NA_real_,
    MAPE = NA_real_,
    Bias = NA_real_,
    n_test = 0L
  )
}

fit_all_gee_pred <- function(
    sim_data,
    base_params,
    scenario,
    run,
    corstr = "exchangeable"
) {

  all_results <- list()

  for (outcome_name in names(base_params)) {

    dat <- prepare_sim_model_data(sim_data, outcome_name)
    train_data <- dat %>% dplyr::filter(split == "train")
    test_last <- get_test_last_observation(dat)

    form <- build_sim_model_formula(outcome_name)

    fit <- tryCatch({
      geepack::geeglm(
        formula = form,
        id = ID,
        data = train_data,
        family = gaussian(),
        corstr = corstr,
        std.err = "san.se",
        na.action = na.exclude
      )
    }, error = function(e) {
      message("GEE failed | outcome = ", outcome_name,
              " | scenario = ", scenario,
              " | run = ", run,
              " | ", e$message)
      NULL
    })

    if (is.null(fit) || nrow(test_last) == 0) {
      tmp <- make_failed_metric_rows(
        scenario = scenario,
        run = run,
        outcome_name = outcome_name,
        model = "GEE",
        model_detail = "GEE_exchangeable"
      )
    } else {
      train_fitted_log <- as.numeric(
        predict(fit, newdata = train_data, type = "response")
      )

      smearing_factor <- calculate_sim_smearing_factor(
        observed_y = train_data[[outcome_name]],
        fitted_log = train_fitted_log,
        label = paste("GEE", outcome_name, scenario, run)
      )

      pred_log <- as.numeric(
        predict(fit, newdata = test_last, type = "response")
      )

      pred_df <- test_last %>%
        dplyr::transmute(
          scenario = scenario,
          run = run,
          ID = as.character(ID),
          Group = as.character(Group),
          Time_in_Years = Time_in_Years,
          outcome = outcome_name,
          outcome_label = SIM_OUTCOME_LABELS[[outcome_name]],
          model = "GEE",
          model_detail = "GEE_exchangeable_empirical_smearing",
          target_definition = "last_simulated_test_target_with_history",
          smearing_factor = smearing_factor,
          y_true = .data[[outcome_name]],
          y_pred = exp(pred_log) * smearing_factor
        )

      tmp <- calc_metrics_by_group(pred_df)
    }

    all_results[[outcome_name]] <- tmp
  }

  dplyr::bind_rows(all_results)
}

fit_all_lme_pred <- function(
    sim_data,
    base_params,
    scenario,
    run,
    random_slope = TRUE
) {

  all_results <- list()

  for (outcome_name in names(base_params)) {

    dat <- prepare_sim_model_data(sim_data, outcome_name)
    train_data <- dat %>% dplyr::filter(split == "train")
    test_last <- get_test_last_observation(dat)

    form <- build_sim_model_formula(outcome_name)

    model_detail <- ifelse(random_slope, "LME_RI_RS", "LME_RI")

    fit <- tryCatch({
      nlme::lme(
        fixed = form,
        random = if (random_slope) ~ Time_in_Years | ID else ~ 1 | ID,
        data = train_data,
        method = "REML",
        na.action = na.exclude,
        control = nlme::lmeControl(
          opt = "optim",
          maxIter = 150,
          msMaxIter = 150,
          niterEM = 30,
          returnObject = TRUE
        )
      )
    }, error = function(e) {
      message(model_detail, " failed | outcome = ", outcome_name,
              " | scenario = ", scenario,
              " | run = ", run,
              " | ", e$message)
      NULL
    })

    if (is.null(fit) || nrow(test_last) == 0) {
      tmp <- make_failed_metric_rows(
        scenario = scenario,
        run = run,
        outcome_name = outcome_name,
        model = "LME",
        model_detail = model_detail
      )
    } else {
      # 測試集病患視為新病患，因此 level = 0，只使用固定效果預測。
      # 這與目前資料分析.R 及 DL 的最後一次觀測目標預測邏輯一致。
      train_fitted_log <- as.numeric(
        predict(fit, newdata = train_data, level = 0)
      )

      smearing_factor <- calculate_sim_smearing_factor(
        observed_y = train_data[[outcome_name]],
        fitted_log = train_fitted_log,
        label = paste(model_detail, outcome_name, scenario, run)
      )

      pred_log <- as.numeric(
        predict(fit, newdata = test_last, level = 0)
      )

      pred_df <- test_last %>%
        dplyr::transmute(
          scenario = scenario,
          run = run,
          ID = as.character(ID),
          Group = as.character(Group),
          Time_in_Years = Time_in_Years,
          outcome = outcome_name,
          outcome_label = SIM_OUTCOME_LABELS[[outcome_name]],
          model = "LME",
          model_detail = paste0(model_detail, "_level0_empirical_smearing"),
          target_definition = "last_simulated_test_target_with_history",
          smearing_factor = smearing_factor,
          y_true = .data[[outcome_name]],
          y_pred = exp(pred_log) * smearing_factor
        )

      tmp <- calc_metrics_by_group(pred_df)
    }

    all_results[[outcome_name]] <- tmp
  }

  dplyr::bind_rows(all_results)
}
