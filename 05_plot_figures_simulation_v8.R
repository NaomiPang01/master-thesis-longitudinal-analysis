# ============================================================
# 05_plot_figures_simulation_v8.R
#
# 功能：
# 1. 繪製 S0 固定效果參考軌跡圖
# 2. 由五情境整合摘要檔繪製 S0–S4 五模型預測誤差比較圖
# 3. 繪製附錄 E1、E2 所需之 S0 代表性模擬資料 LME 殘差診斷圖
#
# 使用前提：
# - 工作目錄設定於 Simulation 資料夾
# - simulation_outputs/S0 至 S4 已存在
# - summary_all_scenarios 整合輸出已完成
# ============================================================


# ============================================================
# 0. 載入套件
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "ggh4x",
  "patchwork",
  "purrr",
  "tibble",
  "nlme"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "尚未安裝以下套件：",
    paste(missing_packages, collapse = ", "),
    "\n請先執行：install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(purrr)
library(tibble)
library(nlme)

# ============================================================
# 1. 全域設定
# ============================================================

SIM_OUTCOMES <- c(
  "HbA1c",
  "Glu_AC",
  "LDL",
  "TG"
)

SIM_OUTCOME_LABELS <- c(
  HbA1c = "HbA1c",
  Glu_AC = "Glu-AC",
  LDL = "LDL",
  TG = "TG"
)

SIM_GROUP_LABELS <- c(
  Overall = "整體",
  Control = "控制組",
  Statin = "史他汀組"
)

SIM_MODEL_LEVELS <- c(
  "GEE",
  "LME",
  "LSTM-NARX",
  "TCN",
  "ODE-LSTM"
)

SIM_METRIC_LEVELS <- c(
  "RMSE",
  "MAE",
  "MAPE",
  "Bias"
)

SIM_SCENARIOS <- paste0("S", 0:4)

# 模擬輸出根目錄
simulation_root <- "simulation_outputs"

# 五情境整合摘要檔
all_scenarios_summary_file <- file.path(
  simulation_root,
  "summary_all_scenarios",
  "all_scenarios_five_model_summary_group.csv"
)

# 五情境比較圖輸出目錄
all_scenarios_figure_dir <- file.path(
  simulation_root,
  "summary_all_scenarios",
  "figures"
)

dir.create(
  all_scenarios_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. 固定效果參考軌跡圖
# ============================================================

plot_reference_trajectory <- function(
    out_dir,
    scenario_name = "S0"
) {
  
  ref_path <- file.path(
    out_dir,
    "summary_tables",
    "reference_trajectory.csv"
  )
  
  if (!file.exists(ref_path)) {
    warning(
      "找不到 reference_trajectory.csv：",
      ref_path
    )
    return(NULL)
  }
  
  ref <- readr::read_csv(
    ref_path,
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      outcome = dplyr::recode(
        as.character(outcome),
        HbA1C = "HbA1c",
        Glucose_AC = "Glu_AC",
        Triglycerides = "TG",
        .default = as.character(outcome)
      ),
      outcome = factor(
        outcome,
        levels = SIM_OUTCOMES
      ),
      outcome_label = factor(
        unname(SIM_OUTCOME_LABELS[as.character(outcome)]),
        levels = unname(SIM_OUTCOME_LABELS[SIM_OUTCOMES])
      ),
      Group = factor(
        Group,
        levels = c("Control", "Statin")
      ),
      Group_label = factor(
        unname(
          SIM_GROUP_LABELS[
            as.character(Group)
          ]
        ),
        levels = unname(
          SIM_GROUP_LABELS[
            c("Control", "Statin")
          ]
        )
      )
    )
  
  if (any(is.na(ref$outcome_label))) {
    stop(
      "reference_trajectory.csv 中含有無法辨識的 outcome。"
    )
  }
  
  if (any(is.na(ref$Group_label))) {
    stop(
      "reference_trajectory.csv 中含有無法辨識的 Group。"
    )
  }
  
  p <- ggplot2::ggplot(
    ref,
    ggplot2::aes(
      x = Time_in_Years,
      y = fit,
      color = Group_label,
      group = Group_label
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(
      ~ outcome_label,
      scales = "free_y",
      nrow = 2
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::labs(
      title = paste0(
        scenario_name,
        " 情境下四項主要代謝指標之固定效果參考軌跡"
      ),
      x = "追蹤時間（年）",
      y = "模擬平均值",
      color = "組別"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )
  
  fig_dir <- file.path(
    out_dir,
    "figures"
  )
  
  dir.create(
    fig_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  out_file <- file.path(
    fig_dir,
    paste0(
      scenario_name,
      "_reference_trajectory.png"
    )
  )
  
  ggplot2::ggsave(
    filename = out_file,
    plot = p,
    width = 10,
    height = 7,
    dpi = 600,
    bg = "white"
  )
  
  message(
    "參考軌跡圖已儲存至：",
    normalizePath(out_file)
  )
  
  return(p)
}


# ============================================================
# 3. 五模型分組預測誤差比較圖
# ============================================================

plot_simulation_metrics_by_group <- function(
    summary_file,
    scenario_name = "S0",
    figure_dir
) {
  
  if (!scenario_name %in% SIM_SCENARIOS) {
    stop(
      "scenario_name 必須為：",
      paste(SIM_SCENARIOS, collapse = ", ")
    )
  }
  
  if (!file.exists(summary_file)) {
    stop(
      "找不到整合輸出檔：",
      summary_file
    )
  }
  
  dat <- readr::read_csv(
    summary_file,
    show_col_types = FALSE
  ) %>%
    dplyr::filter(
      scenario == scenario_name
    )
  
  if (nrow(dat) == 0) {
    stop(
      "整合輸出檔中找不到情境：",
      scenario_name
    )
  }
  
  required_columns <- c(
    "scenario",
    "outcome",
    "outcome_label",
    "model",
    "Group_eval",
    "RMSE_mean",
    "RMSE_sd",
    "MAE_mean",
    "MAE_sd",
    "MAPE_mean",
    "MAPE_sd",
    "Bias_mean",
    "Bias_sd"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "整合輸出檔缺少欄位：",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  dat <- dat %>%
    dplyr::mutate(
      outcome = factor(
        outcome,
        levels = SIM_OUTCOMES
      ),
      outcome_label = factor(
        outcome_label,
        levels = unname(
          SIM_OUTCOME_LABELS[
            SIM_OUTCOMES
          ]
        )
      ),
      model = factor(
        model,
        levels = SIM_MODEL_LEVELS
      ),
      Group_eval = factor(
        Group_eval,
        levels = c(
          "Overall",
          "Control",
          "Statin"
        )
      ),
      Group_label_plot = factor(
        Group_eval,
        levels = c(
          "Overall",
          "Control",
          "Statin"
        ),
        labels = c(
          "整體",
          "控制組",
          "史他汀組"
        )
      )
    )
  
  if (any(is.na(dat$outcome))) {
    stop("資料中含有無法辨識的 outcome。")
  }
  
  if (any(is.na(dat$model))) {
    stop("資料中含有無法辨識的 model。")
  }
  
  if (any(is.na(dat$Group_label_plot))) {
    stop("資料中含有無法辨識的 Group_eval。")
  }
  
  metric_long <- dat %>%
    dplyr::select(
      scenario,
      outcome,
      outcome_label,
      model,
      Group_eval,
      Group_label_plot,
      RMSE_mean,
      RMSE_sd,
      MAE_mean,
      MAE_sd,
      MAPE_mean,
      MAPE_sd,
      Bias_mean,
      Bias_sd
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        RMSE_mean,
        RMSE_sd,
        MAE_mean,
        MAE_sd,
        MAPE_mean,
        MAPE_sd,
        Bias_mean,
        Bias_sd
      ),
      names_to = c(
        "metric",
        ".value"
      ),
      names_pattern = "^(RMSE|MAE|MAPE|Bias)_(mean|sd)$"
    ) %>%
    dplyr::mutate(
      metric = factor(
        metric,
        levels = SIM_METRIC_LEVELS
      ),
      errorbar_lower = dplyr::if_else(
        as.character(metric) == "Bias",
        mean - sd,
        pmax(mean - sd, 0)
      ),
      errorbar_upper = mean + sd
    )
  
  dodge_width <- 0.55
  
  bias_reference <- tidyr::expand_grid(
    metric = factor(
      "Bias",
      levels = SIM_METRIC_LEVELS
    ),
    outcome_label = factor(
      unname(
        SIM_OUTCOME_LABELS[
          SIM_OUTCOMES
        ]
      ),
      levels = unname(
        SIM_OUTCOME_LABELS[
          SIM_OUTCOMES
        ]
      )
    ),
    yintercept = 0
  )
  
  p <- ggplot2::ggplot(
    metric_long
  ) +
    
    ggplot2::geom_hline(
      data = bias_reference,
      mapping = ggplot2::aes(
        yintercept = yintercept
      ),
      linetype = "dashed",
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    
    ggplot2::geom_errorbar(
      mapping = ggplot2::aes(
        x = model,
        ymin = errorbar_lower,
        ymax = errorbar_upper,
        group = Group_label_plot
      ),
      position = ggplot2::position_dodge(
        width = dodge_width
      ),
      width = 0.16,
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = model,
        y = mean,
        shape = Group_label_plot,
        group = Group_label_plot
      ),
      position = ggplot2::position_dodge(
        width = dodge_width
      ),
      size = 2.4,
      na.rm = TRUE
    ) +
    
    ggh4x::facet_grid2(
      metric ~ outcome_label,
      scales = "free_y",
      independent = "y"
    ) +
    
    ggplot2::scale_shape_manual(
      name = "組別",
      breaks = c(
        "整體",
        "控制組",
        "史他汀組"
      ),
      labels = c(
        "整體",
        "控制組",
        "史他汀組"
      ),
      values = c(
        "整體" = 16,
        "控制組" = 17,
        "史他汀組" = 15
      ),
      drop = FALSE
    ) +
    
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 3
        )
      )
    ) +
    
    ggplot2::theme_bw(
      base_size = 11
    ) +
    
    ggplot2::labs(
      title = paste0(
        scenario_name,
        " 情境下五模型預測誤差比較"
      ),
      x = "模型",
      y = "評估指標平均值",
      caption = paste0(
        "點表示各次成功模擬之平均值；",
        "誤差棒表示平均值 ± 1 個標準差。"
      )
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 10
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 9
      ),
      panel.spacing = grid::unit(
        0.6,
        "lines"
      )
    )
  
  dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  out_file <- file.path(
    figure_dir,
    paste0(
      scenario_name,
      "_five_model_metrics_by_group.png"
    )
  )
  
  ggplot2::ggsave(
    filename = out_file,
    plot = p,
    width = 14,
    height = 9,
    dpi = 600,
    bg = "white"
  )
  
  message(
    "五模型比較圖已儲存至：",
    normalizePath(out_file)
  )
  
  return(p)
}


# ============================================================
# 4. 附錄 E1、E2：S0 代表性模擬資料 LME 殘差診斷
# ============================================================

plot_S0_residual_diagnostics <- function(
    representative_data_path,
    out_dir,
    scenario_name = "S0"
) {
  
  if (!file.exists(representative_data_path)) {
    warning(
      "找不到代表性 S0 模擬資料：",
      representative_data_path
    )
    return(NULL)
  }
  
  sim_data <- readr::read_csv(
    representative_data_path,
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      ID = as.factor(ID),
      Group = factor(
        Group,
        levels = c(
          "Control",
          "Statin"
        )
      ),
      GroupStatin = as.integer(
        Group == "Statin"
      ),
      Sex1 = as.integer(Sex1)
    )
  
  fit_lme_and_extract_diag <- function(
    data,
    outcome_name
  ) {
    
    form <- stats::as.formula(
      paste0(
        "log(",
        outcome_name,
        ") ~ ",
        "GroupStatin * Time_in_Years + ",
        "Age_baseline + Sex1 + BMI"
      )
    )
    
    fit <- tryCatch(
      {
        nlme::lme(
          fixed = form,
          random = ~ Time_in_Years | ID,
          data = data,
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
      },
      error = function(e) {
        warning(
          "S0 residual diagnostic LME failed | outcome = ",
          outcome_name,
          " | ",
          e$message
        )
        NULL
      }
    )
    
    if (is.null(fit)) {
      return(
        tibble::tibble()
      )
    }
    
    dat_fit <- nlme::getData(fit)
    
    tibble::tibble(ID = as.character(dat_fit$ID),
      outcome = outcome_name,
      outcome_label = SIM_OUTCOME_LABELS[[outcome_name]],
      Time_in_Years = dat_fit$Time_in_Years,
      fitted = as.numeric(fitted(fit,level = 1)),
      residual = as.numeric(residuals(fit,type = "normalized")),
      raw_residual = as.numeric(residuals(fit,type = "response"))
    )
  }
  
  diag_df <- purrr::map_dfr(
    SIM_OUTCOMES,
    ~ fit_lme_and_extract_diag(
      sim_data,
      .x
    )
  ) %>%
    dplyr::mutate(
      outcome = factor(
        outcome,
        levels = SIM_OUTCOMES
      ),
      outcome_label = factor(
        outcome_label,
        levels = unname(
          SIM_OUTCOME_LABELS[
            SIM_OUTCOMES
          ]
        )
      )
    )
  
  if (nrow(diag_df) == 0) {
    warning(
      "S0 residual diagnostics 沒有可用資料。"
    )
    return(NULL)
  }
  
  fig_dir <- file.path(
    out_dir,
    "figures"
  )
  
  summary_dir <- file.path(
    out_dir,
    "summary_tables"
  )
  
  dir.create(
    fig_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    summary_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  readr::write_csv(
    diag_df,
    file.path(
      summary_dir,
      paste0(
        scenario_name,
        "_lme_residual_diagnostics_data.csv"
      )
    )
  )
  
  qq_df <- diag_df %>%
    dplyr::group_by(
      outcome,
      outcome_label
    ) %>%
    dplyr::arrange(
      residual,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      n = dplyr::n(),
      theoretical = stats::qnorm(
        stats::ppoints(n)
      ),
      sample = sort(residual)
    ) %>%
    dplyr::ungroup()
  
  qq_line_df <- qq_df %>%
    dplyr::group_by(
      outcome,
      outcome_label
    ) %>%
    dplyr::summarise(
      x_q1 = stats::quantile(
        theoretical,
        0.25,
        na.rm = TRUE
      ),
      x_q3 = stats::quantile(
        theoretical,
        0.75,
        na.rm = TRUE
      ),
      y_q1 = stats::quantile(
        sample,
        0.25,
        na.rm = TRUE
      ),
      y_q3 = stats::quantile(
        sample,
        0.75,
        na.rm = TRUE
      ),
      slope = (
        y_q3 - y_q1
      ) / (
        x_q3 - x_q1
      ),
      intercept = y_q1 -
        slope * x_q1,
      x_min = min(
        theoretical,
        na.rm = TRUE
      ),
      x_max = max(
        theoretical,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_min = intercept +
        slope * x_min,
      y_max = intercept +
        slope * x_max
    )
  
  base_theme <- ggplot2::theme_bw(
    base_size = 11
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  p_res_fit <- ggplot2::ggplot(
    diag_df,
    ggplot2::aes(
      x = fitted,
      y = residual
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.25,
      size = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~ outcome_label,
      scales = "free_x",
      ncol = 1
    ) +
    ggplot2::labs(
      x = "配適值",
      y = "標準化殘差",
      title = "殘差對配適值圖"
    ) +
    base_theme
  
  p_qq <- ggplot2::ggplot(
    qq_df,
    ggplot2::aes(
      x = theoretical,
      y = sample
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.35,
      size = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      data = qq_line_df,
      ggplot2::aes(
        x = x_min,
        xend = x_max,
        y = y_min,
        yend = y_max
      ),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    ggplot2::facet_wrap(
      ~ outcome_label,
      scales = "free",
      ncol = 1
    ) +
    ggplot2::labs(
      x = "理論常態分位數",
      y = "樣本分位數",
      title = "常態 Q-Q 圖"
    ) +
    base_theme
  
  fig_E1 <- p_res_fit +
    p_qq +
    patchwork::plot_layout(
      ncol = 2
    )
  
  E1_file <- file.path(
    fig_dir,
    paste0(
      scenario_name,
      "_Appendix_Figure_E1_residual_fitted_QQ.png"
    )
  )
  
  ggplot2::ggsave(
    filename = E1_file,
    plot = fig_E1,
    width = 12,
    height = 9,
    dpi = 600,
    bg = "white"
  )
  
  p_hist <- ggplot2::ggplot(
    diag_df,
    ggplot2::aes(
      x = residual
    )
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(
        y = ggplot2::after_stat(
          density
        )
      ),
      bins = 35,
      color = "white",
      na.rm = TRUE
    ) +
    ggplot2::geom_density(
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~ outcome_label,
      scales = "free",
      ncol = 1
    ) +
    ggplot2::labs(
      x = "標準化殘差",
      y = "密度",
      title = "殘差直方圖"
    ) +
    base_theme
  
  p_time <- ggplot2::ggplot(
    diag_df,
    ggplot2::aes(
      x = Time_in_Years,
      y = residual
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.25,
      size = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~ outcome_label,
      scales = "free_x",
      ncol = 1
    ) +
    ggplot2::labs(
      x = "追蹤時間（年）",
      y = "標準化殘差",
      title = "殘差對追蹤時間圖"
    ) +
    base_theme
  
  fig_E2 <- p_hist +
    p_time +
    patchwork::plot_layout(
      ncol = 2
    )
  
  E2_file <- file.path(
    fig_dir,
    paste0(
      scenario_name,
      "_Appendix_Figure_E2_residual_hist_time.png"
    )
  )
  
  ggplot2::ggsave(
    filename = E2_file,
    plot = fig_E2,
    width = 12,
    height = 9,
    dpi = 600,
    bg = "white"
  )
  
  message(
    "附錄圖 E1 已儲存至：",
    normalizePath(E1_file)
  )
  
  message(
    "附錄圖 E2 已儲存至：",
    normalizePath(E2_file)
  )
  
  invisible(
    list(
      E1 = fig_E1,
      E2 = fig_E2,
      diagnostics = diag_df
    )
  )
}

# ============================================================
# 5. 主執行區
# ============================================================

# 檢查整合摘要檔是否存在
if (!file.exists(all_scenarios_summary_file)) {
  stop(
    "找不到五情境整合摘要檔：",
    all_scenarios_summary_file
  )
}

# 建立圖檔輸出資料夾
dir.create(
  all_scenarios_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 一次產生 S0–S4 五模型預測誤差比較圖
plot_list_metrics_s0_s4 <- lapply(
  SIM_SCENARIOS,
  function(s) {
    plot_simulation_metrics_by_group(
      summary_file = all_scenarios_summary_file,
      scenario_name = s,
      figure_dir = all_scenarios_figure_dir
    )
  }
)

# 為圖形清單加上情境名稱
names(plot_list_metrics_s0_s4) <- SIM_SCENARIOS

message(
  "S0–S4 五模型預測誤差比較圖已全部完成。"
)

# ============================================================
# 6. 完成訊息
# ============================================================

message(
  "05_plot_figures_simulation_v8.R 執行完成。"
)
