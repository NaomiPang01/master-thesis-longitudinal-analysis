############################################################
# 完整 cohort 多重插補推論 + prediction-safe 單一資料預測
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(geepack)
  library(nlme)
  library(mitools)
  library(tibble)
  library(patchwork)
  library(writexl)
})

set.seed(5014)

SCRIPT_VERSION <- "scheme-A-final-v18"

# 第一次執行請維持 FALSE。
# 檢查 AIC 輸出後，若四項 outcome 均支持 LME-RI+RS，
# 再改為 TRUE，並重新執行完整 RScript。
AIC_SELECTION_REVIEWED <- TRUE

OBSERVED_FLAG_MAP <- c(
  HbA1c = "HbA1c_observed",
  Glu_AC = "Glu_AC_observed",
  LDL = "LDL_observed",
  TG = "TG_observed"
)

EXPECTED_DL_MODELS <- c("LSTM-NARX", "TCN", "ODE-LSTM")
OUTCOMES <- c("HbA1c", "Glu_AC", "LDL", "TG")

OUTCOME_LABELS <- c(
  HbA1c = "HbA1c",
  Glu_AC = "Glu-AC",
  LDL = "LDL",
  TG = "TG"
)

MODEL_LEVELS_ALL <- c("GEE", "LME_RI", "LME_RI_RS")

MODEL_LABELS <- c(
  GEE = "GEE",
  LME_RI = "LME-RI",
  LME_RI_RS = "LME-RI+RS"
)

MODEL_LEVELS_FIVE <- c(
  "GEE", "LME", "LSTM-NARX", "TCN", "ODE-LSTM"
)

GROUP_LEVELS <- c("Overall", "Control", "Statin")

GROUP_LABELS_ZH <- c(
  Overall = "整體",
  Control = "控制組",
  Statin = "史他汀組"
)

GROUP_LABELS_EN <- c(
  Overall = "Overall",
  Control = "Control group",
  Statin = "Statin group"
)

MODEL_COVARS <- c("Age_baseline", "Sex", "BMI")
GEE_CORSTR <- "exchangeable"
MIN_PLOT_TARGET_N <- 5L
BOOTSTRAP_B <- 1000L
BOOTSTRAP_SEED <- 5014L

# ============================================================
# 0. 路徑
# ============================================================

pick_existing_path <- function(candidates, error_message) {
  hit <- candidates[file.exists(candidates)][1]

  if (length(hit) == 0L || is.na(hit)) {
    stop(error_message, call. = FALSE)
  }

  hit
}

DATA_PATH <- pick_existing_path(
  c(
    file.path(
      "analysis_outputs",
      "dl_data",
      "imputed_long_with_split_type2.csv"
    ),
    file.path(
      "dl_data",
      "imputed_long_with_split_type2.csv"
    ),
    "imputed_long_with_split_type2.csv"
  ),
  "找不到 imputed_long_with_split_type2.csv。請先執行最新版表格.R。"
)

PREDICTION_DATA_PATH <- pick_existing_path(
  c(
    file.path(
      "analysis_outputs",
      "dl_data",
      "prediction_long_with_split_type2.csv"
    ),
    file.path(
      "dl_data",
      "prediction_long_with_split_type2.csv"
    ),
    "prediction_long_with_split_type2.csv"
  ),
  "找不到 prediction_long_with_split_type2.csv。請先執行最新版表格.R。"
)

OUT_DIR <- file.path("analysis_outputs", "model_analysis")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")
RDS_DIR <- file.path(OUT_DIR, "rds")

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# 清除 v16 已停用或重新命名之舊輸出，避免新舊結果混用。
legacy_output_files <- c(
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_backtransform_parameters.csv"
  ),
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_smearing_factors.csv"
  ),
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_metrics_by_imp_group.csv"
  ),
  file.path(
    TABLE_DIR,
    "inference_prediction_row_alignment_summary_qc.csv"
  )
)

invisible(
  file.remove(
    legacy_output_files[
      file.exists(legacy_output_files)
    ]
  )
)

# ============================================================
# 1. 共用函數
# ============================================================

parse_logical_safe <- function(x) {
  if (is.logical(x)) return(x)

  x_chr <- tolower(trimws(as.character(x)))

  dplyr::case_when(
    x_chr %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x_chr %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

all_n_imp_equal <- function(x, expected) {
  length(x) > 0L &&
    all(!is.na(x)) &&
    all(as.integer(x) == as.integer(expected))
}

qc_row <- function(check, value, expected, fatal = TRUE) {
  value_num <- as.numeric(value)
  expected_num <- as.numeric(expected)

  tibble(
    check = check,
    value = value_num,
    expected = expected_num,
    fatal = fatal,
    passed = isTRUE(
      length(value_num) == 1L &&
        is.finite(value_num) &&
        value_num == expected_num
    )
  )
}

standardize_input <- function(path) {
  readr::read_csv(
    path,
    show_col_types = FALSE
  ) %>%
    mutate(
      row_id = as.integer(row_id),
      ID = as.character(ID),
      Group = factor(
        Group,
        levels = c("Control", "Statin")
      ),
      Sex = factor(
        Sex,
        levels = c("Female", "Male")
      ),
      split = factor(
        split,
        levels = c("train", "validation", "test")
      ),
      .imp = as.integer(.imp),
      Time_Category = factor(
        Time_Category,
        levels = paste0("T", 1:19),
        ordered = TRUE
      ),
      across(
        all_of(unname(OBSERVED_FLAG_MAP)),
        parse_logical_safe
      )
    ) %>%
    arrange(
      .imp,
      ID,
      Time_in_Years,
      row_id
    )
}

check_id_contiguous <- function(df) {
  if (nrow(df) == 0L) return(FALSE)

  id_runs <- rle(as.character(df$ID))$values

  length(id_runs) == dplyr::n_distinct(df$ID)
}

build_model_formula <- function(outcome_name) {
  rhs <- paste(
    c("Group * Time_in_Years", MODEL_COVARS),
    collapse = " + "
  )

  as.formula(
    paste0(
      "log(",
      outcome_name,
      ") ~ ",
      rhs
    )
  )
}

prepare_model_data <- function(
    df,
    outcome_name,
    use_split = c("all", "train", "test"),
    require_observed_outcome = FALSE
) {
  use_split <- match.arg(use_split)

  observed_col <- unname(
    OBSERVED_FLAG_MAP[[outcome_name]]
  )

  out <- df %>%
    mutate(
      ID = factor(ID),
      Group = factor(
        Group,
        levels = c("Control", "Statin")
      ),
      Sex = factor(
        Sex,
        levels = c("Female", "Male")
      )
    )

  if (use_split != "all") {
    out <- out %>%
      filter(as.character(split) == use_split)
  }

  if (isTRUE(require_observed_outcome)) {
    out <- out %>%
      filter(.data[[observed_col]] %in% TRUE)
  }

  out %>%
    filter(
      is.finite(.data[[outcome_name]]),
      .data[[outcome_name]] > 0,
      is.finite(Time_in_Years),
      is.finite(Age_baseline),
      is.finite(BMI),
      !is.na(Sex),
      !is.na(Group)
    ) %>%
    arrange(
      ID,
      Time_in_Years,
      row_id
    )
}

get_last_observed_target_data <- function(
    df,
    outcome_name,
    split_name = c("train", "validation", "test"),
    min_history = 1L
) {
  split_name <- match.arg(split_name)

  observed_col <- unname(
    OBSERVED_FLAG_MAP[[outcome_name]]
  )

  df_split <- df %>%
    mutate(
      ID = as.character(ID),
      Group = factor(
        Group,
        levels = c("Control", "Statin")
      ),
      Sex = factor(
        Sex,
        levels = c("Female", "Male")
      )
    ) %>%
    filter(
      as.character(split) == split_name
    ) %>%
    arrange(
      ID,
      Time_in_Years,
      row_id
    )

  df_split %>%
    group_by(ID) %>%
    group_modify(
      function(.x, .y) {
        g <- .x %>%
          arrange(
            Time_in_Years,
            row_id
          ) %>%
          mutate(
            target_row_position = row_number(),
            n_split_rows = n(),
            n_originally_observed_outcome = sum(
              .data[[observed_col]] %in% TRUE,
              na.rm = TRUE
            )
          )

        valid_target <- (
          g[[observed_col]] %in% TRUE &
            is.finite(g[[outcome_name]]) &
            g[[outcome_name]] > 0 &
            is.finite(g$Time_in_Years) &
            is.finite(g$Age_baseline) &
            is.finite(g$BMI) &
            !is.na(g$Sex) &
            !is.na(g$Group) &
            g$target_row_position >= (min_history + 1L)
        )

        candidate_pos <- which(valid_target)

        if (length(candidate_pos) == 0L) {
          return(g[0, , drop = FALSE])
        }

        target_pos <- max(candidate_pos)
        target_row <- g[target_pos, , drop = FALSE]

        history_rows <- g[
          seq_len(target_pos - 1L),
          ,
          drop = FALSE
        ]

        valid_history <- (
          history_rows[[observed_col]] %in% TRUE &
            is.finite(history_rows[[outcome_name]]) &
            history_rows[[outcome_name]] > 0 &
            is.finite(history_rows$Time_in_Years)
        )

        history_n_available <- sum(
          valid_history,
          na.rm = TRUE
        )

        if (history_n_available < min_history) {
          return(g[0, , drop = FALSE])
        }

        target_row %>%
          mutate(
            history_length_available = target_pos - 1L,
            history_length_used_for_eligibility =
              history_n_available,
            target_originally_observed = TRUE,
            target_definition = paste0(
              "last_originally_observed_",
              split_name,
              "_target_with_observed_history"
            )
          )
      }
    ) %>%
    ungroup() %>%
    mutate(ID = factor(ID)) %>%
    arrange(
      ID,
      Time_in_Years,
      row_id
    )
}

get_test_last_target_data <- function(df, outcome_name) {
  get_last_observed_target_data(
    df = df,
    outcome_name = outcome_name,
    split_name = "test"
  ) %>%
    rename(
      n_test_rows = n_split_rows
    )
}

calc_pred_metrics <- function(y_true, y_pred) {
  keep <- (
    is.finite(y_true) &
      is.finite(y_pred) &
      y_true > 0
  )

  y_true <- y_true[keep]
  y_pred <- y_pred[keep]

  if (length(y_true) == 0L) {
    return(
      tibble(
        RMSE = NA_real_,
        MAE = NA_real_,
        MAPE = NA_real_,
        Bias = NA_real_,
        n_test = 0L
      )
    )
  }

  error <- y_pred - y_true

  tibble(
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    MAPE = mean(abs(error / y_true)) * 100,
    Bias = mean(error),
    n_test = length(y_true)
  )
}

calc_metrics_by_group <- function(pred_df) {
  bind_rows(
    pred_df %>% mutate(Group_eval = "Overall"),
    pred_df %>% mutate(Group_eval = as.character(Group))
  ) %>%
    group_by(
      .imp,
      outcome,
      model,
      Group_eval
    ) %>%
    summarise(
      calc_pred_metrics(y_true, y_pred),
      .groups = "drop"
    ) %>%
    mutate(
      Group_eval = factor(
        Group_eval,
        levels = GROUP_LEVELS
      ),
      Group_label = GROUP_LABELS_ZH[
        as.character(Group_eval)
      ]
    )
}

safe_fit_list <- function(x) {
  x[
    !vapply(
      x,
      is.null,
      logical(1)
    )
  ]
}

calculate_smearing_factor <- function(
    observed_y,
    fitted_log,
    label
) {
  valid <- (
    is.finite(observed_y) &
      observed_y > 0 &
      is.finite(fitted_log)
  )

  if (!any(valid)) {
    stop(
      "無法計算 smearing factor：沒有有效訓練資料。 ",
      label,
      call. = FALSE
    )
  }

  sf <- mean(
    exp(
      log(observed_y[valid]) -
        fitted_log[valid]
    ),
    na.rm = TRUE
  )

  if (!is.finite(sf) || sf <= 0) {
    stop(
      "無法計算 smearing factor：",
      label,
      call. = FALSE
    )
  }

  sf
}

# ============================================================
# 2. 讀取資料與輸入 QC
# ============================================================

dt_long <- standardize_input(DATA_PATH)
dt_prediction <- standardize_input(PREDICTION_DATA_PATH)

required_cols <- c(
  "row_id",
  "ID",
  "Group",
  "split",
  ".imp",
  "Time_in_Years",
  "Time_Category",
  OUTCOMES,
  MODEL_COVARS,
  unname(OBSERVED_FLAG_MAP)
)

missing_cols_full <- setdiff(
  required_cols,
  names(dt_long)
)

missing_cols_pred <- setdiff(
  required_cols,
  names(dt_prediction)
)

if (length(missing_cols_full) > 0L) {
  stop(
    "完整補值資料缺少欄位：",
    paste(missing_cols_full, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_cols_pred) > 0L) {
  stop(
    "prediction-safe 資料缺少欄位：",
    paste(missing_cols_pred, collapse = ", "),
    call. = FALSE
  )
}

input_qc <- tibble(
  check = c(
    "number_of_imputations",
    "duplicate_row_id_within_imputation",
    "patient_multiple_split_within_imputation",
    "patient_split_inconsistent_across_imputations",
    "patient_group_inconsistent_across_imputations",
    "prediction_dataset_count",
    "prediction_BMI_complete"
  ),
  value = c(
    n_distinct(dt_long$.imp),
    nrow(
      dt_long %>%
        count(.imp, row_id) %>%
        filter(n > 1L)
    ),
    nrow(
      dt_long %>%
        distinct(.imp, ID, split) %>%
        count(.imp, ID) %>%
        filter(n != 1L)
    ),
    nrow(
      dt_long %>%
        distinct(ID, split) %>%
        count(ID) %>%
        filter(n != 1L)
    ),
    nrow(
      dt_long %>%
        distinct(ID, Group) %>%
        count(ID) %>%
        filter(n != 1L)
    ),
    n_distinct(dt_prediction$.imp),
    sum(is.na(dt_prediction$BMI))
  ),
  expected = c(
    10L,
    0L,
    0L,
    0L,
    0L,
    1L,
    0L
  )
) %>%
  mutate(
    passed = value == expected
  )

observed_flag_input_qc <- map_dfr(
  names(OBSERVED_FLAG_MAP),
  function(outcome_name) {
    flag_name <- unname(
      OBSERVED_FLAG_MAP[[outcome_name]]
    )

    inconsistent_n <- dt_long %>%
      group_by(row_id) %>%
      summarise(
        n_flag_value = n_distinct(
          .data[[flag_name]],
          na.rm = FALSE
        ),
        .groups = "drop"
      ) %>%
      filter(
        n_flag_value != 1L
      ) %>%
      nrow()

    tibble(
      check = c(
        paste0(flag_name, "_parse_missing"),
        paste0(
          flag_name,
          "_inconsistent_across_imputations"
        )
      ),
      value = c(
        sum(is.na(dt_long[[flag_name]])),
        inconsistent_n
      ),
      expected = c(0L, 0L),
      passed = value == expected
    )
  }
)

row_count_by_imp <- dt_long %>%
  count(
    .imp,
    name = "row_n"
  )

row_id_imp_qc <- dt_long %>%
  group_by(row_id) %>%
  summarise(
    n_imp = n_distinct(.imp),
    .groups = "drop"
  ) %>%
  filter(
    n_imp != 10L
  )

row_structure_input_qc <- tibble(
  check = c(
    "row_count_identical_across_imputations",
    "row_id_present_in_all_imputations"
  ),
  value = c(
    n_distinct(row_count_by_imp$row_n),
    nrow(row_id_imp_qc)
  ),
  expected = c(1L, 0L),
  passed = value == expected
)

prediction_flag_missing_n <- sum(
  is.na(
    as.matrix(
      dt_prediction %>%
        select(
          all_of(unname(OBSERVED_FLAG_MAP))
        )
    )
  )
)

inference_same_time_duplicate_n <- dt_long %>%
  count(
    .imp,
    ID,
    Time_in_Years,
    name = "n"
  ) %>%
  filter(n > 1L) %>%
  nrow()

prediction_same_time_duplicate_n <- dt_prediction %>%
  count(
    ID,
    Time_in_Years,
    name = "n"
  ) %>%
  filter(n > 1L) %>%
  nrow()

prediction_structure_qc <- tibble(
  check = c(
    "prediction_duplicate_row_id",
    "prediction_patient_multiple_split",
    "prediction_patient_multiple_group",
    "prediction_group_parse_missing",
    "prediction_sex_parse_missing",
    "prediction_split_parse_missing",
    "prediction_observed_flag_parse_missing",
    "inference_same_patient_same_time_duplicate",
    "prediction_same_patient_same_time_duplicate"
  ),
  value = c(
    nrow(
      dt_prediction %>%
        count(row_id) %>%
        filter(n > 1L)
    ),
    nrow(
      dt_prediction %>%
        distinct(ID, split) %>%
        count(ID) %>%
        filter(n != 1L)
    ),
    nrow(
      dt_prediction %>%
        distinct(ID, Group) %>%
        count(ID) %>%
        filter(n != 1L)
    ),
    sum(is.na(dt_prediction$Group)),
    sum(is.na(dt_prediction$Sex)),
    sum(is.na(dt_prediction$split)),
    prediction_flag_missing_n,
    inference_same_time_duplicate_n,
    prediction_same_time_duplicate_n
  ),
  expected = rep(0L, 9L)
) %>%
  mutate(
    passed = value == expected
  )

inference_row_manifest <- dt_long %>%
  filter(.imp == 1L) %>%
  transmute(
    row_id,
    ID = as.character(ID),
    Time_in_Years
  ) %>%
  distinct()

prediction_row_manifest <- dt_prediction %>%
  transmute(
    row_id,
    ID = as.character(ID),
    Time_in_Years
  ) %>%
  distinct()

row_manifest_alignment_qc <- full_join(
  inference_row_manifest %>%
    rename(
      ID_inference = ID,
      Time_in_Years_inference = Time_in_Years
    ),
  prediction_row_manifest %>%
    rename(
      ID_prediction = ID,
      Time_in_Years_prediction = Time_in_Years
    ),
  by = "row_id"
) %>%
  mutate(
    passed = (
      !is.na(ID_inference) &
        !is.na(ID_prediction) &
        ID_inference == ID_prediction &
        is.finite(Time_in_Years_inference) &
        is.finite(Time_in_Years_prediction) &
        abs(
          Time_in_Years_inference -
            Time_in_Years_prediction
        ) <= 1e-12
    )
  )

row_manifest_alignment_summary_qc <- tibble(
  check = "inference_prediction_row_alignment",
  value = sum(!row_manifest_alignment_qc$passed),
  expected = 0L,
  passed = value == expected
)

write_csv(
  row_manifest_alignment_qc,
  file.path(
    TABLE_DIR,
    "inference_prediction_row_alignment_qc.csv"
  )
)

input_qc <- bind_rows(
  input_qc,
  observed_flag_input_qc,
  row_structure_input_qc,
  prediction_structure_qc,
  row_manifest_alignment_summary_qc
)

write_csv(
  input_qc,
  file.path(
    TABLE_DIR,
    "model_analysis_input_qc.csv"
  )
)

if (any(!input_qc$passed)) {
  print(
    input_qc %>%
      filter(!passed)
  )

  stop(
    "資料分析輸入 QC 未通過。",
    call. = FALSE
  )
}

imp_list <- split(
  dt_long,
  dt_long$.imp
)

imp_list_prediction <- split(
  dt_prediction,
  dt_prediction$.imp
)

split_check_patient <- dt_long %>%
  distinct(
    .imp,
    ID,
    Group,
    split
  ) %>%
  count(
    .imp,
    Group,
    split,
    name = "patient_n"
  ) %>%
  arrange(
    .imp,
    Group,
    split
  )

split_check_observation <- dt_long %>%
  count(
    .imp,
    Group,
    split,
    name = "observation_n"
  ) %>%
  arrange(
    .imp,
    Group,
    split
  )

write_csv(
  split_check_patient,
  file.path(
    TABLE_DIR,
    "split_check_patient_by_imp.csv"
  )
)

write_csv(
  split_check_observation,
  file.path(
    TABLE_DIR,
    "split_check_observation_by_imp.csv"
  )
)

prediction_evaluation_scope <- tibble(
  analysis_part = "GEE/LME prediction evaluation",
  training_data =
    "train set; all originally observed outcome rows",
  evaluation_data = "test set",
  evaluation_unit =
    "last originally observed target per evaluable test patient",
  note = paste0(
    "推論分析使用完整多重插補資料；預測比較以 train set 中原始觀測 outcome 列配適 GEE/LME，",
    "並於 test set 每位病患最後一次原始觀測 target 評估。prediction-safe 資料之 outcome 不插補，",
    "BMI 僅使用 train set 資訊補值。GEE 與 LME 均使用訓練資料殘差估計之 empirical Duan ",
    "smearing factor 進行原尺度反轉換；LME 對訓練資料及測試資料均採固定效果層級預測 ",
    "（level = 0），因此測試集預測代表未納入測試病患個體隨機效果之 population-level prediction；",
    "validation set 僅供 DL 調參。"
  )
)

write_csv(
  prediction_evaluation_scope,
  file.path(
    TABLE_DIR,
    "prediction_evaluation_scope.csv"
  )
)

# ============================================================
# 3. Target 樣本數
# ============================================================

target_sample_size_by_split <- map_dfr(
  imp_list_prediction,
  function(df_i) {
    map_dfr(
      OUTCOMES,
      function(outcome_name) {
        map_dfr(
          c("train", "validation", "test"),
          function(split_name) {
            target_df <- get_last_observed_target_data(
              df_i,
              outcome_name,
              split_name
            )

            tibble(
              .imp = unique(df_i$.imp)[1],
              outcome = outcome_name,
              outcome_label =
                OUTCOME_LABELS[[outcome_name]],
              split = split_name,
              target_patient_n =
                n_distinct(target_df$ID),
              target_row_n = nrow(target_df)
            )
          }
        )
      }
    )
  }
)

target_sample_size_summary <- target_sample_size_by_split %>%
  group_by(
    outcome,
    outcome_label,
    split
  ) %>%
  summarise(
    target_n_min = min(target_patient_n),
    target_n_max = max(target_patient_n),
    target_n = first(target_patient_n),
    dataset_n = n_distinct(.imp),
    single_dataset_consistent = (
      target_n_min == target_n_max &
        dataset_n == length(imp_list_prediction)
    ),
    .groups = "drop"
  )

gee_lme_training_observation_count <- map_dfr(
  imp_list_prediction,
  function(df_i) {
    map_dfr(
      OUTCOMES,
      function(outcome_name) {
        train_df <- prepare_model_data(
          df_i,
          outcome_name,
          use_split = "train",
          require_observed_outcome = TRUE
        )

        tibble(
          .imp = unique(df_i$.imp)[1],
          outcome = outcome_name,
          outcome_label =
            OUTCOME_LABELS[[outcome_name]],
          train_observation_n = nrow(train_df),
          train_patient_n =
            n_distinct(train_df$ID)
        )
      }
    )
  }
)

write_csv(
  target_sample_size_by_split,
  file.path(
    TABLE_DIR,
    "prediction_target_sample_size_by_split_imp.csv"
  )
)

write_csv(
  target_sample_size_summary,
  file.path(
    TABLE_DIR,
    "prediction_target_sample_size_by_split.csv"
  )
)

write_csv(
  gee_lme_training_observation_count,
  file.path(
    TABLE_DIR,
    "GEE_LME_training_observation_count.csv"
  )
)

if (
  any(
    !target_sample_size_summary$
      single_dataset_consistent
  )
) {
  stop(
    "prediction-safe target 樣本數 QC 未通過。",
    call. = FALSE
  )
}

# ============================================================
# 4. GEE / LME 推論
# ============================================================

fit_gee_inference_one_imp <- function(
    df,
    outcome_name
) {
  dat <- prepare_model_data(
    df,
    outcome_name,
    "all"
  )

  if (!check_id_contiguous(dat)) {
    stop(
      "GEE inference 資料中相同 ID 未連續排列：",
      outcome_name,
      call. = FALSE
    )
  }

  tryCatch(
    geepack::geeglm(
      formula = build_model_formula(outcome_name),
      id = ID,
      data = dat,
      family = gaussian(),
      corstr = GEE_CORSTR,
      std.err = "san.se",
      na.action = na.exclude
    ),
    error = function(e) {
      message(
        "GEE inference failed: ",
        outcome_name,
        " | ",
        e$message
      )
      NULL
    }
  )
}

fit_lme_inference_one_imp <- function(
    df,
    outcome_name,
    random_slope
) {
  dat <- prepare_model_data(
    df,
    outcome_name,
    "all"
  )

  tryCatch(
    nlme::lme(
      fixed = build_model_formula(outcome_name),
      random = if (random_slope) {
        ~ Time_in_Years | ID
      } else {
        ~ 1 | ID
      },
      data = dat,
      method = "REML",
      na.action = na.exclude,
      control = nlme::lmeControl(
        opt = "optim",
        maxIter = 150,
        msMaxIter = 150,
        niterEM = 30,
        returnObject = FALSE
      )
    ),
    error = function(e) {
      message(
        "LME inference failed: ",
        outcome_name,
        " | ",
        e$message
      )
      NULL
    }
  )
}

run_inference <- function(
    model_type,
    outcome_name
) {
  map(
    imp_list,
    function(df_i) {
      if (model_type == "GEE") {
        fit_gee_inference_one_imp(
          df_i,
          outcome_name
        )
      } else if (model_type == "LME_RI") {
        fit_lme_inference_one_imp(
          df_i,
          outcome_name,
          FALSE
        )
      } else {
        fit_lme_inference_one_imp(
          df_i,
          outcome_name,
          TRUE
        )
      }
    }
  )
}

gee_inference_fits <- setNames(
  map(
    OUTCOMES,
    ~ run_inference("GEE", .x)
  ),
  OUTCOMES
)

lme_ri_inference_fits <- setNames(
  map(
    OUTCOMES,
    ~ run_inference("LME_RI", .x)
  ),
  OUTCOMES
)

lme_rirs_inference_fits <- setNames(
  map(
    OUTCOMES,
    ~ run_inference("LME_RI_RS", .x)
  ),
  OUTCOMES
)

extract_inference_fit_qc <- function(
    fit_nested_list,
    model_name
) {
  map_dfr(
    names(fit_nested_list),
    function(outcome_name) {
      fits <- fit_nested_list[[outcome_name]]

      map_dfr(
        seq_along(fits),
        function(j) {
          fit <- fits[[j]]
          is_lme <- inherits(fit, "lme")

          tibble(
            .imp = j,
            outcome = outcome_name,
            outcome_label =
              OUTCOME_LABELS[[outcome_name]],
            model = model_name,
            model_label =
              MODEL_LABELS[[model_name]],
            fit_ok = !is.null(fit),
            num_iter = if (
              is_lme &&
                !is.null(fit)
            ) {
              suppressWarnings(
                as.integer(fit$numIter)[1]
              )
            } else {
              NA_integer_
            },
            apvar_ok = if (
              is_lme &&
                !is.null(fit)
            ) {
              is.matrix(fit$apVar) &&
                all(is.finite(fit$apVar))
            } else {
              NA
            },
            apvar_note = if (
              is_lme &&
                !is.null(fit) &&
                !is.matrix(fit$apVar)
            ) {
              as.character(fit$apVar)[1]
            } else {
              NA_character_
            }
          )
        }
      )
    }
  )
}

inference_fit_qc <- bind_rows(
  extract_inference_fit_qc(
    gee_inference_fits,
    "GEE"
  ),
  extract_inference_fit_qc(
    lme_ri_inference_fits,
    "LME_RI"
  ),
  extract_inference_fit_qc(
    lme_rirs_inference_fits,
    "LME_RI_RS"
  )
)

write_csv(
  inference_fit_qc,
  file.path(
    TABLE_DIR,
    "inference_model_fit_qc.csv"
  )
)

if (any(!inference_fit_qc$fit_ok)) {
  print(
    inference_fit_qc %>%
      filter(!fit_ok)
  )

  stop(
    "至少一組推論模型配適失敗。",
    call. = FALSE
  )
}

saveRDS(
  list(
    gee = gee_inference_fits,
    lme_ri = lme_ri_inference_fits,
    lme_rirs = lme_rirs_inference_fits
  ),
  file.path(
    RDS_DIR,
    "inference_model_fits_GEE_LME_type2.rds"
  )
)

pool_coef_table <- function(
    fit_list,
    model_name,
    outcome_name
) {
  fit_list <- safe_fit_list(fit_list)

  if (length(fit_list) == 0L) {
    return(
      tibble(
        model = model_name,
        model_label =
          MODEL_LABELS[[model_name]],
        outcome = outcome_name,
        outcome_label =
          OUTCOME_LABELS[[outcome_name]],
        term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p_value = NA_real_,
        p_value_fmt = NA_character_,
        lower = NA_real_,
        upper = NA_real_,
        miss_info = NA_real_,
        n_imp = 0L
      )
    )
  }

  coef_list <- lapply(
    fit_list,
    function(fit) {
      if (inherits(fit, "lme")) {
        nlme::fixed.effects(fit)
      } else {
        coef(fit)
      }
    }
  )

  vcov_list <- lapply(
    fit_list,
    vcov
  )

  term_list <- lapply(
    coef_list,
    names
  )

  reference_terms <- term_list[[1]]

  term_sets_identical <- all(
    vapply(
      term_list,
      function(x) {
        identical(
          x,
          reference_terms
        )
      },
      logical(1)
    )
  )

  if (!term_sets_identical) {
    print(
      tibble(
        imputation_index =
          seq_along(term_list),
        terms = vapply(
          term_list,
          paste,
          collapse = " | ",
          character(1)
        )
      )
    )

    stop(
      "Rubin 彙整失敗：不同插補資料的係數名稱或順序不一致。 ",
      model_name,
      " / ",
      outcome_name,
      call. = FALSE
    )
  }

  common_terms <- reference_terms

  coef_mat <- do.call(
    rbind,
    lapply(
      coef_list,
      function(x) {
        as.numeric(
          x[common_terms]
        )
      }
    )
  )

  colnames(coef_mat) <- common_terms

  within_var_mat <- do.call(
    rbind,
    lapply(
      vcov_list,
      function(vmat) {
        vmat <- as.matrix(vmat)

        as.numeric(
          diag(
            vmat[
              common_terms,
              common_terms,
              drop = FALSE
            ]
          )
        )
      }
    )
  )

  colnames(within_var_mat) <- common_terms

  m <- nrow(coef_mat)
  q_bar <- colMeans(coef_mat)
  u_bar <- colMeans(within_var_mat)

  b_var <- if (m > 1L) {
    apply(
      coef_mat,
      2,
      var
    )
  } else {
    rep(
      0,
      length(common_terms)
    )
  }

  total_var <- (
    u_bar +
      (1 + 1 / m) * b_var
  )

  std_error <- sqrt(total_var)

  df_old <- ifelse(
    is.finite(b_var) &
      b_var > 0 &
      is.finite(u_bar) &
      u_bar > 0,
    (m - 1) *
      (
        1 +
          u_bar /
          ((1 + 1 / m) * b_var)
      )^2,
    Inf
  )

  statistic <- q_bar / std_error

  p_value <- ifelse(
    is.finite(statistic),
    2 * pt(
      abs(statistic),
      df = df_old,
      lower.tail = FALSE
    ),
    NA_real_
  )

  crit <- ifelse(
    is.finite(df_old),
    qt(
      0.975,
      df = df_old
    ),
    qnorm(0.975)
  )

  miss_lookup <- tryCatch(
    {
      mi_res <- mitools::MIcombine(
        results = coef_list,
        variances = vcov_list
      )

      s_df <- as.data.frame(
        summary(mi_res)
      )

      s_df$term <- rownames(s_df)
      rownames(s_df) <- NULL

      names(s_df) <- gsub(
        "[()]",
        "",
        names(s_df)
      )

      names(s_df) <- gsub(
        " ",
        "_",
        names(s_df)
      )

      if ("missInfo" %in% names(s_df)) {
        s_df %>%
          select(
            term,
            miss_info = missInfo
          )
      } else {
        tibble(
          term = common_terms,
          miss_info = NA_real_
        )
      }
    },
    error = function(e) {
      tibble(
        term = common_terms,
        miss_info = NA_real_
      )
    }
  )

  tibble(
    model = model_name,
    model_label =
      MODEL_LABELS[[model_name]],
    outcome = outcome_name,
    outcome_label =
      OUTCOME_LABELS[[outcome_name]],
    term = common_terms,
    estimate = as.numeric(q_bar),
    std_error = as.numeric(std_error),
    statistic = as.numeric(statistic),
    df = as.numeric(df_old),
    p_value = as.numeric(p_value),
    p_value_fmt = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    lower = as.numeric(
      q_bar - crit * std_error
    ),
    upper = as.numeric(
      q_bar + crit * std_error
    ),
    n_imp = m
  ) %>%
    left_join(
      miss_lookup,
      by = "term"
    ) %>%
    relocate(
      miss_info,
      .after = upper
    )
}

pooled_parameter_estimates <- bind_rows(
  map_dfr(
    OUTCOMES,
    ~ pool_coef_table(
      gee_inference_fits[[.x]],
      "GEE",
      .x
    )
  ),
  map_dfr(
    OUTCOMES,
    ~ pool_coef_table(
      lme_ri_inference_fits[[.x]],
      "LME_RI",
      .x
    )
  ),
  map_dfr(
    OUTCOMES,
    ~ pool_coef_table(
      lme_rirs_inference_fits[[.x]],
      "LME_RI_RS",
      .x
    )
  )
)

write_csv(
  pooled_parameter_estimates,
  file.path(
    TABLE_DIR,
    "pooled_parameter_estimates_Rubin_GEE_LME.csv"
  )
)

# ============================================================
# 5. LME AIC 比較與人工確認
# ============================================================

extract_lme_aic_by_imp <- function(
    fit_nested_list,
    model_name
) {
  map_dfr(
    names(fit_nested_list),
    function(outcome_name) {
      fits <- fit_nested_list[[outcome_name]]

      map_dfr(
        seq_along(fits),
        function(j) {
          fit <- fits[[j]]

          tibble(
            .imp = as.integer(j),
            outcome = outcome_name,
            outcome_label =
              OUTCOME_LABELS[[outcome_name]],
            model = model_name,
            model_label =
              MODEL_LABELS[[model_name]],
            AIC = if (!is.null(fit)) {
              suppressWarnings(
                as.numeric(stats::AIC(fit))
              )
            } else {
              NA_real_
            }
          )
        }
      )
    }
  )
}

lme_aic_by_imp <- bind_rows(
  extract_lme_aic_by_imp(
    lme_ri_inference_fits,
    "LME_RI"
  ),
  extract_lme_aic_by_imp(
    lme_rirs_inference_fits,
    "LME_RI_RS"
  )
) %>%
  arrange(
    factor(outcome, levels = OUTCOMES),
    .imp,
    factor(
      model,
      levels = c(
        "LME_RI",
        "LME_RI_RS"
      )
    )
  )

lme_aic_pairwise_by_imp <- lme_aic_by_imp %>%
  select(
    .imp,
    outcome,
    outcome_label,
    model,
    AIC
  ) %>%
  pivot_wider(
    names_from = model,
    values_from = AIC
  ) %>%
  mutate(
    pair_complete =
      is.finite(LME_RI) &
      is.finite(LME_RI_RS),
    delta_AIC_RI_minus_RIRS =
      if_else(
        pair_complete,
        LME_RI - LME_RI_RS,
        NA_real_
      ),
    preferred_model = case_when(
      !pair_complete ~ NA_character_,
      LME_RI_RS < LME_RI ~
        "LME_RI_RS",
      LME_RI < LME_RI_RS ~
        "LME_RI",
      TRUE ~ "Tie"
    )
  ) %>%
  arrange(
    factor(outcome, levels = OUTCOMES),
    .imp
  )

lme_aic_comparison <- lme_aic_by_imp %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label
  ) %>%
  summarise(
    AIC_mean = safe_mean(AIC),
    AIC_sd = safe_sd(AIC),
    AIC_min = safe_min(AIC),
    AIC_max = safe_max(AIC),
    n_imp = sum(is.finite(AIC)),
    .groups = "drop"
  ) %>%
  arrange(
    factor(outcome, levels = OUTCOMES),
    factor(
      model,
      levels = c(
        "LME_RI",
        "LME_RI_RS"
      )
    )
  )

lme_aic_selection_summary <- lme_aic_pairwise_by_imp %>%
  group_by(
    outcome,
    outcome_label
  ) %>%
  summarise(
    n_imp_total = n(),
    n_imp_complete =
      sum(pair_complete),
    n_RI_preferred = sum(
      preferred_model == "LME_RI",
      na.rm = TRUE
    ),
    n_RIRS_preferred = sum(
      preferred_model == "LME_RI_RS",
      na.rm = TRUE
    ),
    n_tie = sum(
      preferred_model == "Tie",
      na.rm = TRUE
    ),
    delta_AIC_mean = safe_mean(
      delta_AIC_RI_minus_RIRS
    ),
    delta_AIC_sd = safe_sd(
      delta_AIC_RI_minus_RIRS
    ),
    delta_AIC_min = safe_min(
      delta_AIC_RI_minus_RIRS
    ),
    delta_AIC_max = safe_max(
      delta_AIC_RI_minus_RIRS
    ),
    mean_AIC_preferred_model =
      case_when(
        !is.finite(delta_AIC_mean) ~
          NA_character_,
        delta_AIC_mean > 0 ~
          "LME_RI_RS",
        delta_AIC_mean < 0 ~
          "LME_RI",
        TRUE ~ "Tie"
      ),
    all_imputations_prefer_RIRS = (
      n_imp_complete ==
        length(imp_list) &
        n_RIRS_preferred ==
        length(imp_list)
    ),
    .groups = "drop"
  ) %>%
  arrange(
    factor(outcome, levels = OUTCOMES)
  )

write_csv(
  lme_aic_by_imp,
  file.path(
    TABLE_DIR,
    "lme_aic_by_imputation.csv"
  )
)

write_csv(
  lme_aic_pairwise_by_imp,
  file.path(
    TABLE_DIR,
    "lme_aic_pairwise_by_imputation.csv"
  )
)

write_csv(
  lme_aic_comparison,
  file.path(
    TABLE_DIR,
    "lme_aic_comparison_for_paper.csv"
  )
)

write_csv(
  lme_aic_selection_summary,
  file.path(
    TABLE_DIR,
    "lme_aic_selection_summary.csv"
  )
)

cat(
  "\n============================================================\n",
  "LME 隨機結構 AIC 比較結果\n",
  "delta_AIC_RI_minus_RIRS > 0：LME-RI+RS AIC 較小\n",
  "delta_AIC_RI_minus_RIRS < 0：LME-RI AIC 較小\n",
  "============================================================\n",
  sep = ""
)

print(
  lme_aic_selection_summary,
  n = Inf,
  width = Inf
)

aic_all_complete <- all(
  lme_aic_selection_summary$
    n_imp_complete ==
    length(imp_list)
)

if (!aic_all_complete) {
  print(
    lme_aic_selection_summary %>%
      filter(
        n_imp_complete !=
          length(imp_list)
      )
  )

  stop(
    paste0(
      "部分 outcome 的 LME-RI 與 LME-RI+RS ",
      "未在十組插補資料中取得完整配對 AIC。"
    ),
    call. = FALSE
  )
}

if (!isTRUE(AIC_SELECTION_REVIEWED)) {
  stop(
    paste0(
      "程式已完成 LME AIC 輸出，並暫停於模型選擇階段。\n",
      "請檢查：\n",
      "1. lme_aic_by_imputation.csv\n",
      "2. lme_aic_pairwise_by_imputation.csv\n",
      "3. lme_aic_selection_summary.csv\n\n",
      "確認四項 outcome 是否皆支持 LME-RI+RS 後，",
      "將程式最前方 AIC_SELECTION_REVIEWED 改為 TRUE，",
      "再重新執行完整 RScript。"
    ),
    call. = FALSE
  )
}

all_outcomes_mean_AIC_prefer_RIRS <- all(
  !is.na(
    lme_aic_selection_summary$
      mean_AIC_preferred_model
  ) &
    lme_aic_selection_summary$
      mean_AIC_preferred_model ==
      "LME_RI_RS"
)

if (!all_outcomes_mean_AIC_prefer_RIRS) {
  print(
    lme_aic_selection_summary %>%
      filter(
        is.na(mean_AIC_preferred_model) |
          mean_AIC_preferred_model !=
          "LME_RI_RS"
      )
  )

  stop(
    paste0(
      "並非四項代謝指標均由平均 REML AIC 支持 LME-RI+RS。\n",
      "目前程式不得繼續固定使用 LME-RI+RS；",
      "需依 outcome 選擇主要 LME 模型後再調整後續程式。"
    ),
    call. = FALSE
  )
}

MAIN_LME_MODEL <- setNames(
  rep(
    "LME_RI_RS",
    length(OUTCOMES)
  ),
  OUTCOMES
)

MODEL_LEVELS_MAIN <- c(
  "GEE",
  "LME_RI_RS"
)

cat(
  "\n四項代謝指標之平均 REML AIC 均支持 LME-RI+RS。\n",
  "後續 LME 診斷、預測圖與五模型比較將使用 LME-RI+RS。\n\n",
  sep = ""
)

# ============================================================
# 6. LME variance components 與 GEE correlation
# ============================================================

extract_lme_varcomp_one <- function(
    fit,
    outcome_name,
    model_name,
    imp_id
) {
  if (is.null(fit)) {
    return(
      tibble(
        .imp = imp_id,
        outcome = outcome_name,
        outcome_label =
          OUTCOME_LABELS[[outcome_name]],
        model = model_name,
        model_label =
          MODEL_LABELS[[model_name]],
        sigma_b0_sq = NA_real_,
        sigma_b1_sq = NA_real_,
        rho_b0b1 = NA_real_,
        sigma_sq = NA_real_
      )
    )
  }

  D_mat <- tryCatch(
    as.matrix(
      getVarCov(
        fit,
        type = "random.effects"
      )
    ),
    error = function(e) {
      NULL
    }
  )

  sigma_b0_sq <- if (
    !is.null(D_mat) &&
      nrow(D_mat) >= 1L
  ) {
    as.numeric(D_mat[1, 1])
  } else {
    NA_real_
  }

  sigma_b1_sq <- if (
    !is.null(D_mat) &&
      nrow(D_mat) >= 2L
  ) {
    as.numeric(D_mat[2, 2])
  } else {
    NA_real_
  }

  rho_b0b1 <- if (
    !is.null(D_mat) &&
      nrow(D_mat) >= 2L &&
      is.finite(sigma_b0_sq) &&
      is.finite(sigma_b1_sq) &&
      sigma_b0_sq > 0 &&
      sigma_b1_sq > 0
  ) {
    as.numeric(
      D_mat[1, 2] /
        sqrt(
          sigma_b0_sq *
            sigma_b1_sq
        )
    )
  } else {
    NA_real_
  }

  tibble(
    .imp = imp_id,
    outcome = outcome_name,
    outcome_label =
      OUTCOME_LABELS[[outcome_name]],
    model = model_name,
    model_label =
      MODEL_LABELS[[model_name]],
    sigma_b0_sq = sigma_b0_sq,
    sigma_b1_sq = sigma_b1_sq,
    rho_b0b1 = rho_b0b1,
    sigma_sq = fit$sigma^2
  )
}

extract_lme_varcomp <- function(
    fit_nested_list,
    model_name
) {
  map_dfr(
    names(fit_nested_list),
    function(outcome_name) {
      map_dfr(
        seq_along(
          fit_nested_list[[outcome_name]]
        ),
        function(j) {
          extract_lme_varcomp_one(
            fit_nested_list[[outcome_name]][[j]],
            outcome_name,
            model_name,
            j
          )
        }
      )
    }
  )
}

lme_variance_components_by_imp <- bind_rows(
  extract_lme_varcomp(
    lme_ri_inference_fits,
    "LME_RI"
  ),
  extract_lme_varcomp(
    lme_rirs_inference_fits,
    "LME_RI_RS"
  )
)

lme_variance_components_summary_for_paper <-
  lme_variance_components_by_imp %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label
  ) %>%
  summarise(
    sigma_b0_sq_mean =
      safe_mean(sigma_b0_sq),
    sigma_b0_sq_sd =
      safe_sd(sigma_b0_sq),
    sigma_b1_sq_mean =
      safe_mean(sigma_b1_sq),
    sigma_b1_sq_sd =
      safe_sd(sigma_b1_sq),
    rho_b0b1_mean =
      safe_mean(rho_b0b1),
    rho_b0b1_sd =
      safe_sd(rho_b0b1),
    sigma_sq_mean =
      safe_mean(sigma_sq),
    sigma_sq_sd =
      safe_sd(sigma_sq),
    n_imp = sum(is.finite(sigma_sq)),
    .groups = "drop"
  )

lme_random_structure_qc <-
  lme_variance_components_by_imp %>%
  filter(
    model == "LME_RI_RS"
  ) %>%
  mutate(
    random_slope_near_zero =
      is.finite(sigma_b1_sq) &
      sigma_b1_sq < 1e-6,
    random_correlation_near_boundary =
      is.finite(rho_b0b1) &
      abs(rho_b0b1) > 0.95
  ) %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label
  ) %>%
  summarise(
    n_imp = n(),
    n_valid_slope_variance =
      sum(is.finite(sigma_b1_sq)),
    n_valid_intercept_slope_correlation =
      sum(is.finite(rho_b0b1)),
    n_slope_near_zero =
      sum(
        random_slope_near_zero,
        na.rm = TRUE
      ),
    n_corr_near_boundary =
      sum(
        random_correlation_near_boundary,
        na.rm = TRUE
      ),
    minimum_random_slope_variance =
      safe_min(sigma_b1_sq),
    maximum_random_slope_variance =
      safe_max(sigma_b1_sq),
    maximum_absolute_intercept_slope_correlation =
      safe_max(abs(rho_b0b1)),
    .groups = "drop"
  ) %>%
  mutate(
    slope_variance_qc_passed =
      n_valid_slope_variance ==
      n_imp,
    correlation_qc_passed =
      n_valid_intercept_slope_correlation ==
      n_imp,
    any_slope_near_zero =
      n_slope_near_zero > 0L,
    any_correlation_near_boundary =
      n_corr_near_boundary > 0L
  )

extract_gee_alpha_safe <- function(fit) {
  if (
    is.null(fit) ||
      !inherits(fit, "geeglm")
  ) {
    return(NA_real_)
  }

  alpha_direct <- tryCatch(
    as.numeric(fit$geese$alpha)[1],
    error = function(e) {
      NA_real_
    }
  )

  if (
    length(alpha_direct) == 1L &&
      is.finite(alpha_direct)
  ) {
    return(alpha_direct)
  }

  tryCatch(
    {
      corr_table <- summary(fit)$corr

      if (
        !is.null(corr_table) &&
          nrow(corr_table) >= 1L &&
          "Estimate" %in%
          names(corr_table)
      ) {
        as.numeric(
          corr_table$Estimate[1]
        )
      } else {
        NA_real_
      }
    },
    error = function(e) {
      NA_real_
    }
  )
}

extract_gee_scale_safe <- function(fit) {
  if (
    is.null(fit) ||
      !inherits(fit, "geeglm")
  ) {
    return(NA_real_)
  }

  scale_summary <- tryCatch(
    {
      dispersion_table <-
        summary(fit)$dispersion

      if (
        !is.null(dispersion_table) &&
          nrow(dispersion_table) >= 1L &&
          "Estimate" %in%
          names(dispersion_table)
      ) {
        as.numeric(
          dispersion_table$Estimate[1]
        )
      } else {
        NA_real_
      }
    },
    error = function(e) {
      NA_real_
    }
  )

  if (is.finite(scale_summary)) {
    return(scale_summary)
  }

  tryCatch(
    as.numeric(fit$geese$gamma)[1],
    error = function(e) {
      NA_real_
    }
  )
}

gee_working_correlation_by_imp <- map_dfr(
  names(gee_inference_fits),
  function(outcome_name) {
    fits <- gee_inference_fits[[outcome_name]]

    map_dfr(
      seq_along(fits),
      function(j) {
        fit <- fits[[j]]

        tibble(
          .imp = as.integer(j),
          outcome = outcome_name,
          outcome_label =
            OUTCOME_LABELS[[outcome_name]],
          model = "GEE",
          model_label = "GEE",
          corstr = tryCatch(
            as.character(
              fit$geese$model$corstr
            ),
            error = function(e) {
              NA_character_
            }
          ),
          alpha =
            extract_gee_alpha_safe(fit),
          scale =
            extract_gee_scale_safe(fit)
        )
      }
    )
  }
) %>%
  distinct(
    outcome,
    .imp,
    .keep_all = TRUE
  ) %>%
  arrange(
    factor(outcome, levels = OUTCOMES),
    .imp
  )

gee_working_correlation_summary_for_paper <-
  gee_working_correlation_by_imp %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label,
    corstr
  ) %>%
  summarise(
    alpha_mean = safe_mean(alpha),
    alpha_sd = safe_sd(alpha),
    scale_mean = safe_mean(scale),
    scale_sd = safe_sd(scale),
    n_imp_alpha =
      sum(is.finite(alpha)),
    n_imp_scale =
      sum(is.finite(scale)),
    .groups = "drop"
  ) %>%
  arrange(
    factor(outcome, levels = OUTCOMES)
  )

gee_alpha_failed <-
  gee_working_correlation_by_imp %>%
  filter(
    is.na(corstr) |
      corstr != GEE_CORSTR |
      !is.finite(alpha)
  )

if (nrow(gee_alpha_failed) > 0L) {
  print(gee_alpha_failed)

  stop(
    "部分 GEE working correlation parameter 無法擷取。",
    call. = FALSE
  )
}

write_csv(
  lme_variance_components_by_imp,
  file.path(
    TABLE_DIR,
    "lme_variance_components_by_imp.csv"
  )
)

write_csv(
  lme_variance_components_summary_for_paper,
  file.path(
    TABLE_DIR,
    "lme_variance_components_summary_for_paper.csv"
  )
)

write_csv(
  lme_random_structure_qc,
  file.path(
    TABLE_DIR,
    "lme_random_structure_qc.csv"
  )
)

write_csv(
  gee_working_correlation_by_imp,
  file.path(
    TABLE_DIR,
    "gee_working_correlation_by_imp.csv"
  )
)

write_csv(
  gee_working_correlation_summary_for_paper,
  file.path(
    TABLE_DIR,
    "gee_working_correlation_summary_for_paper.csv"
  )
)

# ============================================================
# 7. LME 診斷
# ============================================================

get_main_lme_fit_imp <- function(
    outcome_name,
    imp_index
) {
  selected_model <- unname(
    MAIN_LME_MODEL[[outcome_name]]
  )

  if (
    is.null(selected_model) ||
      is.na(selected_model)
  ) {
    stop(
      "尚未指定 LME 主模型：",
      outcome_name,
      call. = FALSE
    )
  }

  fit <- if (
    selected_model ==
      "LME_RI_RS"
  ) {
    lme_rirs_inference_fits[[outcome_name]][[imp_index]]
  } else {
    lme_ri_inference_fits[[outcome_name]][[imp_index]]
  }

  if (is.null(fit)) {
    stop(
      "指定的 LME 主模型無法取得：",
      outcome_name,
      " / ",
      selected_model,
      " / imputation ",
      imp_index,
      call. = FALSE
    )
  }

  list(
    fit = fit,
    model = selected_model,
    .imp = imp_index
  )
}

get_main_lme_fit <- function(
    outcome_name
) {
  get_main_lme_fit_imp(
    outcome_name,
    1L
  )
}

diagnostic_data <- function(fit) {
  if (is.null(fit)) {
    stop(
      "LME 診斷收到空模型物件。",
      call. = FALSE
    )
  }

  dat_fit <- tryCatch(
    getData(fit),
    error = function(e) {
      NULL
    }
  )

  tibble(
    fitted = as.numeric(
      fitted(
        fit,
        level = 1
      )
    ),
    resid_norm = as.numeric(
      residuals(
        fit,
        type = "normalized"
      )
    ),
    Time_in_Years = if (
      !is.null(dat_fit) &&
        "Time_in_Years" %in%
        names(dat_fit)
    ) {
      dat_fit$Time_in_Years
    } else {
      NA_real_
    }
  )
}

diagnostic_summary_one_fit <- function(
    fit,
    outcome_name,
    model_name,
    imp_index
) {
  d <- diagnostic_data(fit)

  x <- d$resid_norm[
    is.finite(d$resid_norm)
  ]

  x_sd <- if (length(x) >= 2L) {
    sd(x)
  } else {
    NA_real_
  }

  tibble(
    .imp = imp_index,
    outcome = outcome_name,
    outcome_label =
      OUTCOME_LABELS[[outcome_name]],
    model = model_name,
    model_label =
      MODEL_LABELS[[model_name]],
    n = length(x),
    resid_mean = safe_mean(x),
    resid_sd = safe_sd(x),
    resid_min = safe_min(x),
    resid_q1 = if (length(x) > 0L) {
      as.numeric(
        quantile(
          x,
          0.25,
          names = FALSE
        )
      )
    } else {
      NA_real_
    },
    resid_median = if (
      length(x) > 0L
    ) {
      median(x)
    } else {
      NA_real_
    },
    resid_q3 = if (length(x) > 0L) {
      as.numeric(
        quantile(
          x,
          0.75,
          names = FALSE
        )
      )
    } else {
      NA_real_
    },
    resid_max = safe_max(x),
    resid_skewness = if (
      length(x) >= 3L &&
        is.finite(x_sd) &&
        x_sd > 0
    ) {
      mean(
        (
          (x - mean(x)) /
            x_sd
        )^3
      )
    } else {
      NA_real_
    },
    resid_excess_kurtosis = if (
      length(x) >= 4L &&
        is.finite(x_sd) &&
        x_sd > 0
    ) {
      mean(
        (
          (x - mean(x)) /
            x_sd
        )^4
      ) - 3
    } else {
      NA_real_
    },
    abs_resid_gt3_n =
      sum(abs(x) > 3),
    abs_resid_gt3_pct = if (
      length(x) > 0L
    ) {
      mean(abs(x) > 3) * 100
    } else {
      NA_real_
    },
    spearman_resid_fitted =
      suppressWarnings(
        cor(
          d$resid_norm,
          d$fitted,
          method = "spearman",
          use = "complete.obs"
        )
      ),
    spearman_abs_resid_fitted =
      suppressWarnings(
        cor(
          abs(d$resid_norm),
          d$fitted,
          method = "spearman",
          use = "complete.obs"
        )
      ),
    spearman_resid_time =
      suppressWarnings(
        cor(
          d$resid_norm,
          d$Time_in_Years,
          method = "spearman",
          use = "complete.obs"
        )
      )
  )
}

lme_diagnostics_by_imp <- map_dfr(
  OUTCOMES,
  function(outcome_name) {
    map_dfr(
      seq_along(imp_list),
      function(j) {
        main <- get_main_lme_fit_imp(
          outcome_name,
          j
        )

        diagnostic_summary_one_fit(
          fit = main$fit,
          outcome_name = outcome_name,
          model_name = main$model,
          imp_index = j
        )
      }
    )
  }
)

lme_diagnostics_summary <-
  lme_diagnostics_by_imp %>%
  filter(
    .imp == 1L
  ) %>%
  select(
    -.imp
  )

lme_diagnostics_all_imputations_summary <-
  lme_diagnostics_by_imp %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label
  ) %>%
  summarise(
    n_imp = n_distinct(.imp),
    resid_mean_mean =
      safe_mean(resid_mean),
    resid_mean_min =
      safe_min(resid_mean),
    resid_mean_max =
      safe_max(resid_mean),
    resid_sd_mean =
      safe_mean(resid_sd),
    resid_sd_min =
      safe_min(resid_sd),
    resid_sd_max =
      safe_max(resid_sd),
    skewness_mean =
      safe_mean(resid_skewness),
    skewness_min =
      safe_min(resid_skewness),
    skewness_max =
      safe_max(resid_skewness),
    excess_kurtosis_mean =
      safe_mean(
        resid_excess_kurtosis
      ),
    excess_kurtosis_min =
      safe_min(
        resid_excess_kurtosis
      ),
    excess_kurtosis_max =
      safe_max(
        resid_excess_kurtosis
      ),
    abs_resid_gt3_pct_mean =
      safe_mean(abs_resid_gt3_pct),
    abs_resid_gt3_pct_min =
      safe_min(abs_resid_gt3_pct),
    abs_resid_gt3_pct_max =
      safe_max(abs_resid_gt3_pct),
    abs_spearman_resid_fitted_max =
      safe_max(
        abs(spearman_resid_fitted)
      ),
    abs_spearman_abs_resid_fitted_max =
      safe_max(
        abs(
          spearman_abs_resid_fitted
        )
      ),
    abs_spearman_resid_time_max =
      safe_max(
        abs(spearman_resid_time)
      ),
    .groups = "drop"
  )

write_csv(
  lme_diagnostics_summary,
  file.path(
    TABLE_DIR,
    "lme_diagnostics_summary_first_imputation.csv"
  )
)

write_csv(
  lme_diagnostics_by_imp,
  file.path(
    TABLE_DIR,
    "lme_diagnostics_by_imputation.csv"
  )
)

write_csv(
  lme_diagnostics_all_imputations_summary,
  file.path(
    TABLE_DIR,
    "lme_diagnostics_all_imputations_summary.csv"
  )
)

make_lme_diagnostic_components <- function(
    fit,
    outcome_name,
    lang = c("zh", "en")
) {
  lang <- match.arg(lang)

  d <- diagnostic_data(fit)

  outcome_lab <-
    OUTCOME_LABELS[[outcome_name]]

  if (lang == "zh") {
    titles <- c(
      "殘差對配適值",
      "常態 Q-Q 圖",
      "殘差分布",
      "殘差對追蹤時間"
    )

    xlp <- c(
      "對數尺度配適值",
      "理論分位數",
      "標準化殘差",
      "追蹤時間（年）"
    )

    ylp <- c(
      "標準化殘差",
      "樣本分位數",
      "次數",
      "標準化殘差"
    )
  } else {
    titles <- c(
      "Residuals versus fitted values",
      "Normal Q-Q plot",
      "Residual distribution",
      "Residuals versus follow-up time"
    )

    xlp <- c(
      "Fitted values on the log scale",
      "Theoretical quantiles",
      "Normalized residuals",
      "Follow-up time (years)"
    )

    ylp <- c(
      "Normalized residuals",
      "Sample quantiles",
      "Count",
      "Normalized residuals"
    )
  }

  base_theme <- theme_bw(
    base_size = 8
  ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 9
      ),
      axis.title = element_text(
        size = 8
      ),
      axis.text = element_text(
        size = 7
      )
    )

  p1 <- ggplot(
    d,
    aes(
      fitted,
      resid_norm
    )
  ) +
    geom_point(
      alpha = 0.20,
      size = 0.45,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    base_theme +
    labs(
      title = paste0(
        outcome_lab,
        if (lang == "zh") "：" else ": ",
        titles[1]
      ),
      x = xlp[1],
      y = ylp[1]
    )

  p2 <- ggplot(
    d,
    aes(
      sample = resid_norm
    )
  ) +
    stat_qq(
      alpha = 0.30,
      size = 0.45,
      na.rm = TRUE
    ) +
    stat_qq_line(
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    base_theme +
    labs(
      title = paste0(
        outcome_lab,
        if (lang == "zh") "：" else ": ",
        titles[2]
      ),
      x = xlp[2],
      y = ylp[2]
    )

  p3 <- ggplot(
    d,
    aes(resid_norm)
  ) +
    geom_histogram(
      bins = 40,
      color = "white",
      na.rm = TRUE
    ) +
    base_theme +
    labs(
      title = paste0(
        outcome_lab,
        if (lang == "zh") "：" else ": ",
        titles[3]
      ),
      x = xlp[3],
      y = ylp[3]
    )

  p4 <- ggplot(
    d,
    aes(
      Time_in_Years,
      resid_norm
    )
  ) +
    geom_point(
      alpha = 0.20,
      size = 0.45,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    base_theme +
    labs(
      title = paste0(
        outcome_lab,
        if (lang == "zh") "：" else ": ",
        titles[4]
      ),
      x = xlp[4],
      y = ylp[4]
    )

  list(p1, p2, p3, p4)
}

for (lang in c("zh", "en")) {
  component_list <- list()

  for (outcome_name in OUTCOMES) {
    main <- get_main_lme_fit(
      outcome_name
    )

    component_list <- c(
      component_list,
      make_lme_diagnostic_components(
        main$fit,
        outcome_name,
        lang
      )
    )
  }

  panel_grid <-
    patchwork::wrap_plots(
      component_list,
      ncol = 4
    )

  title_text <- if (lang == "zh") {
    paste0(
      "四項主要代謝指標之 LME 主模型診斷圖\n",
      "（第一組多重插補資料）"
    )
  } else {
    paste0(
      "LME main-model diagnostics for four metabolic indicators\n",
      "(First multiply imputed dataset)"
    )
  }

  title_panel <-
    patchwork::wrap_elements(
      full = grid::textGrob(
        label = title_text,
        x = 0.5,
        y = 0.5,
        hjust = 0.5,
        vjust = 0.5,
        gp = grid::gpar(
          fontsize = if (
            lang == "zh"
          ) {
            13
          } else {
            12.5
          },
          fontface = "bold",
          lineheight = 1.05
        )
      )
    )

  p_all <- (
    title_panel /
      panel_grid
  ) +
    patchwork::plot_layout(
      heights = c(0.075, 1)
    )

  ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(
        "LME_diagnostic_all_outcomes_4x4_",
        lang,
        ".png"
      )
    ),
    plot = p_all,
    width = 14,
    height = 10.8,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    bg = "white"
  )
}

extract_random_effects_first_imp <- function(
    outcome_name
) {
  main <- get_main_lme_fit(
    outcome_name
  )

  re_df <- as.data.frame(
    ranef(main$fit)
  )

  re_df$ID <- rownames(re_df)

  re_df %>%
    pivot_longer(
      cols = -ID,
      names_to = "effect",
      values_to = "random_effect"
    ) %>%
    mutate(
      outcome = outcome_name,
      outcome_label =
        OUTCOME_LABELS[[outcome_name]],
      model = main$model,
      model_label =
        MODEL_LABELS[[main$model]]
    )
}

lme_random_effects_first_imp <- map_dfr(
  OUTCOMES,
  extract_random_effects_first_imp
) %>%
  group_by(
    outcome,
    effect
  ) %>%
  mutate(
    effect_z = {
      s <- sd(
        random_effect,
        na.rm = TRUE
      )

      if (
        is.finite(s) &&
          s > 0
      ) {
        (
          random_effect -
            mean(
              random_effect,
              na.rm = TRUE
            )
        ) / s
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()

lme_random_effects_summary <-
  lme_random_effects_first_imp %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label,
    effect
  ) %>%
  summarise(
    n_patient = sum(
      is.finite(random_effect)
    ),
    mean = safe_mean(random_effect),
    sd = safe_sd(random_effect),
    min = safe_min(random_effect),
    q1 = if (
      any(is.finite(random_effect))
    ) {
      as.numeric(
        quantile(
          random_effect,
          0.25,
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    },
    median = if (
      any(is.finite(random_effect))
    ) {
      median(
        random_effect,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    q3 = if (
      any(is.finite(random_effect))
    ) {
      as.numeric(
        quantile(
          random_effect,
          0.75,
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    },
    max = safe_max(random_effect),
    n_abs_z_gt_3 = sum(
      abs(effect_z) > 3,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  lme_random_effects_first_imp,
  file.path(
    TABLE_DIR,
    "lme_random_effects_first_imputation.csv"
  )
)

write_csv(
  lme_random_effects_summary,
  file.path(
    TABLE_DIR,
    "lme_random_effects_summary_first_imputation.csv"
  )
)

plot_random_effect_qq <- function(
    re_df,
    lang = c("zh", "en")
) {
  lang <- match.arg(lang)

  plot_df <- re_df %>%
    mutate(
      outcome = factor(
        as.character(outcome),
        levels = OUTCOMES
      ),
      outcome_label_plot = factor(
        OUTCOME_LABELS[
          as.character(outcome)
        ],
        levels = unname(
          OUTCOME_LABELS[OUTCOMES]
        )
      ),
      effect_label = if (
        lang == "zh"
      ) {
        recode(
          as.character(effect),
          `(Intercept)` = "隨機截距",
          Time_in_Years =
            "隨機時間斜率",
          .default =
            as.character(effect)
        )
      } else {
        recode(
          as.character(effect),
          `(Intercept)` =
            "Random intercept",
          Time_in_Years =
            "Random time slope",
          .default =
            as.character(effect)
        )
      },
      effect_label = factor(
        effect_label,
        levels = if (
          lang == "zh"
        ) {
          c(
            "隨機截距",
            "隨機時間斜率"
          )
        } else {
          c(
            "Random intercept",
            "Random time slope"
          )
        }
      )
    ) %>%
    filter(
      !is.na(outcome_label_plot),
      !is.na(effect_label),
      is.finite(effect_z)
    )

  ggplot(
    plot_df,
    aes(
      sample = effect_z
    )
  ) +
    stat_qq(
      alpha = 0.30,
      size = 0.55,
      na.rm = TRUE
    ) +
    stat_qq_line(
      linewidth = 0.50,
      na.rm = TRUE
    ) +
    facet_grid(
      rows = vars(outcome_label_plot),
      cols = vars(effect_label),
      scales = "fixed",
      drop = FALSE
    ) +
    theme_bw(base_size = 10) +
    labs(
      title = if (lang == "zh") {
        "四項主要代謝指標之 LME 標準化隨機效果常態 Q-Q 圖"
      } else {
        "Normal Q-Q plots of standardized LME random effects"
      },
      subtitle = if (lang == "zh") {
        "第一組多重插補資料"
      } else {
        "First multiply imputed dataset"
      },
      x = if (lang == "zh") {
        "理論分位數"
      } else {
        "Theoretical quantiles"
      },
      y = if (lang == "zh") {
        "標準化隨機效果樣本分位數"
      } else {
        "Sample quantiles of standardized random effects"
      }
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9.5,
        margin = margin(b = 7)
      ),
      strip.text = element_text(
        face = "bold",
        size = 9.5
      ),
      axis.title = element_text(
        size = 9.5
      ),
      axis.text = element_text(
        size = 8.5
      ),
      panel.spacing = grid::unit(
        0.12,
        "cm"
      ),
      plot.margin = margin(
        t = 6,
        r = 6,
        b = 6,
        l = 6
      )
    )
}

for (lang in c("zh", "en")) {
  p_re <- plot_random_effect_qq(
    lme_random_effects_first_imp,
    lang
  )

  ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(
        "LME_random_effect_QQ_all_outcomes_",
        lang,
        ".png"
      )
    ),
    plot = p_re,
    width = 11,
    height = 9,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    bg = "white"
  )
}

# ============================================================
# 8. GEE / LME 預測
# ============================================================

fit_gee_prediction_one_imp <- function(
    df,
    outcome_name
) {
  train_data <- prepare_model_data(
    df,
    outcome_name,
    "train",
    TRUE
  )

  test_target <- get_test_last_target_data(
    df,
    outcome_name
  )

  if (
    nrow(train_data) == 0L ||
      nrow(test_target) == 0L
  ) {
    return(
      list(
        fit = NULL,
        pred = NULL,
        metrics = NULL,
        smearing = NULL
      )
    )
  }

  if (!check_id_contiguous(train_data)) {
    stop(
      "GEE prediction 資料中相同 ID 未連續排列：",
      outcome_name,
      call. = FALSE
    )
  }

  imp_id <- unique(df$.imp)[1]

  fit <- tryCatch(
    geepack::geeglm(
      formula =
        build_model_formula(outcome_name),
      id = ID,
      data = train_data,
      family = gaussian(),
      corstr = GEE_CORSTR,
      std.err = "san.se",
      na.action = na.exclude
    ),
    error = function(e) {
      message(
        "GEE prediction failed: ",
        outcome_name,
        " | ",
        e$message
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(
      list(
        fit = NULL,
        pred = NULL,
        metrics = NULL,
        smearing = NULL
      )
    )
  }

  train_fitted_log <- as.numeric(
    predict(
      fit,
      newdata = train_data,
      type = "response"
    )
  )

  smearing_factor <-
    calculate_smearing_factor(
      observed_y =
        train_data[[outcome_name]],
      fitted_log =
        train_fitted_log,
      label = paste(
        "GEE",
        outcome_name
      )
    )

  pred_log <- as.numeric(
    predict(
      fit,
      newdata = test_target,
      type = "response"
    )
  )

  pred_df <- test_target %>%
    transmute(
      .imp = imp_id,
      row_id = as.integer(row_id),
      ID = as.character(ID),
      Group = as.character(Group),
      Time_Category =
        as.character(Time_Category),
      Time_in_Years,
      n_test_rows,
      n_originally_observed_outcome,
      target_row_position,
      history_length_available,
      history_length_used_for_eligibility,
      target_originally_observed,
      target_definition,
      outcome = outcome_name,
      outcome_label =
        OUTCOME_LABELS[[outcome_name]],
      model = "GEE",
      model_label = "GEE",
      prediction_level =
        "population_fixed_effects",
      backtransform_method =
        "empirical_Duan_smearing_train_residuals",
      smearing_factor =
        smearing_factor,
      y_true =
        .data[[outcome_name]],
      y_pred =
        exp(pred_log) *
        smearing_factor
    )

  list(
    fit = fit,
    pred = pred_df,
    metrics =
      calc_metrics_by_group(pred_df),
    smearing = tibble(
      .imp = imp_id,
      outcome = outcome_name,
      outcome_label =
        OUTCOME_LABELS[[outcome_name]],
      model = "GEE",
      model_label = "GEE",
      prediction_level =
        "population_fixed_effects",
      backtransform_method =
        "empirical_Duan_smearing_train_residuals",
      train_observation_n =
        nrow(train_data),
      train_patient_n =
        n_distinct(train_data$ID),
      smearing_factor =
        smearing_factor
    )
  )
}

fit_lme_prediction_one_imp <- function(
    df,
    outcome_name,
    random_slope = TRUE
) {
  train_data <- prepare_model_data(
    df,
    outcome_name,
    "train",
    TRUE
  )

  test_target <- get_test_last_target_data(
    df,
    outcome_name
  )

  if (
    nrow(train_data) == 0L ||
      nrow(test_target) == 0L
  ) {
    return(
      list(
        fit = NULL,
        pred = NULL,
        metrics = NULL,
        smearing = NULL
      )
    )
  }

  imp_id <- unique(df$.imp)[1]

  model_name <- if (random_slope) {
    "LME_RI_RS"
  } else {
    "LME_RI"
  }

  fit <- tryCatch(
    nlme::lme(
      fixed =
        build_model_formula(outcome_name),
      random = if (random_slope) {
        ~ Time_in_Years | ID
      } else {
        ~ 1 | ID
      },
      data = train_data,
      method = "REML",
      na.action = na.exclude,
      control = nlme::lmeControl(
        opt = "optim",
        maxIter = 150,
        msMaxIter = 150,
        niterEM = 30,
        returnObject = FALSE
      )
    ),
    error = function(e) {
      message(
        model_name,
        " prediction failed: ",
        outcome_name,
        " | ",
        e$message
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(
      list(
        fit = NULL,
        pred = NULL,
        metrics = NULL,
        smearing = NULL
      )
    )
  }

  # 使用訓練資料固定效果層級預測估計 empirical smearing factor。
  train_fitted_log <- as.numeric(
    predict(
      fit,
      newdata = train_data,
      level = 0
    )
  )

  smearing_factor <-
    calculate_smearing_factor(
      observed_y =
        train_data[[outcome_name]],
      fitted_log =
        train_fitted_log,
      label = paste(
        model_name,
        outcome_name
      )
    )

  # 測試病患採 level = 0，不帶入個體隨機效果。
  pred_log <- as.numeric(
    predict(
      fit,
      newdata = test_target,
      level = 0
    )
  )

  if (
    length(pred_log) != nrow(test_target) ||
      any(!is.finite(pred_log))
  ) {
    stop(
      "LME 測試資料對數尺度預測無效：",
      model_name,
      " / ",
      outcome_name,
      call. = FALSE
    )
  }

  y_pred_original <- exp(pred_log) *
    smearing_factor

  if (
    any(!is.finite(y_pred_original)) ||
      any(y_pred_original <= 0)
  ) {
    stop(
      "LME empirical smearing 原尺度預測產生非有限值或非正值：",
      model_name,
      " / ",
      outcome_name,
      call. = FALSE
    )
  }

  pred_df <- test_target %>%
    transmute(
      .imp = imp_id,
      row_id = as.integer(row_id),
      ID = as.character(ID),
      Group = as.character(Group),
      Time_Category =
        as.character(Time_Category),
      Time_in_Years,
      n_test_rows,
      n_originally_observed_outcome,
      target_row_position,
      history_length_available,
      history_length_used_for_eligibility,
      target_originally_observed,
      target_definition,
      outcome = outcome_name,
      outcome_label =
        OUTCOME_LABELS[[outcome_name]],
      model = model_name,
      model_label =
        MODEL_LABELS[[model_name]],
      prediction_level =
        "population_fixed_effects_level_0",
      backtransform_method =
        "empirical_Duan_smearing_train_residuals",
      smearing_factor =
        smearing_factor,
      y_true =
        .data[[outcome_name]],
      y_pred =
        y_pred_original
    )

  list(
    fit = fit,
    pred = pred_df,
    metrics =
      calc_metrics_by_group(pred_df),
    smearing = tibble(
      .imp = imp_id,
      outcome = outcome_name,
      outcome_label =
        OUTCOME_LABELS[[outcome_name]],
      model = model_name,
      model_label =
        MODEL_LABELS[[model_name]],
      prediction_level =
        "population_fixed_effects_level_0",
      backtransform_method =
        "empirical_Duan_smearing_train_residuals",
      train_observation_n =
        nrow(train_data),
      train_patient_n =
        n_distinct(train_data$ID),
      smearing_factor =
        smearing_factor
    )
  )
}

run_prediction <- function(
    model_type,
    outcome_name
) {
  map(
    imp_list_prediction,
    function(df_i) {
      if (model_type == "GEE") {
        fit_gee_prediction_one_imp(
          df_i,
          outcome_name
        )
      } else if (
        model_type == "LME_RI"
      ) {
        fit_lme_prediction_one_imp(
          df_i,
          outcome_name,
          FALSE
        )
      } else {
        fit_lme_prediction_one_imp(
          df_i,
          outcome_name,
          TRUE
        )
      }
    }
  )
}

gee_prediction_results <- setNames(
  map(
    OUTCOMES,
    ~ run_prediction("GEE", .x)
  ),
  OUTCOMES
)

lme_ri_prediction_results <- setNames(
  map(
    OUTCOMES,
    ~ run_prediction("LME_RI", .x)
  ),
  OUTCOMES
)

lme_rirs_prediction_results <- setNames(
  map(
    OUTCOMES,
    ~ run_prediction("LME_RI_RS", .x)
  ),
  OUTCOMES
)

extract_nested <- function(
    res_list,
    element_name
) {
  map_dfr(
    res_list,
    function(outcome_res) {
      map_dfr(
        outcome_res,
        function(x) {
          out <- x[[element_name]]

          if (is.null(out)) {
            tibble()
          } else {
            out
          }
        }
      )
    }
  )
}

pred_GEE_LME_individual <- bind_rows(
  extract_nested(
    gee_prediction_results,
    "pred"
  ),
  extract_nested(
    lme_ri_prediction_results,
    "pred"
  ),
  extract_nested(
    lme_rirs_prediction_results,
    "pred"
  )
)

prediction_metrics_GEE_LME_by_imp_group <-
  bind_rows(
    extract_nested(
      gee_prediction_results,
      "metrics"
    ),
    extract_nested(
      lme_ri_prediction_results,
      "metrics"
    ),
    extract_nested(
      lme_rirs_prediction_results,
      "metrics"
    )
  )

prediction_smearing_factors <- bind_rows(
  extract_nested(
    gee_prediction_results,
    "smearing"
  ),
  extract_nested(
    lme_ri_prediction_results,
    "smearing"
  ),
  extract_nested(
    lme_rirs_prediction_results,
    "smearing"
  )
)

prediction_metrics_GEE_LME_summary_group <-
  prediction_metrics_GEE_LME_by_imp_group %>%
  mutate(
    outcome_label = OUTCOME_LABELS[
      as.character(outcome)
    ],
    model_label = MODEL_LABELS[
      as.character(model)
    ],
    dataset_n = 1L
  ) %>%
  arrange(
    factor(outcome, levels = OUTCOMES),
    factor(
      model,
      levels = MODEL_LEVELS_ALL
    ),
    factor(
      Group_eval,
      levels = GROUP_LEVELS
    )
  )

write_csv(
  pred_GEE_LME_individual,
  file.path(
    TABLE_DIR,
    "GEE_LME_test_predictions_individual.csv"
  )
)

write_csv(
  prediction_metrics_GEE_LME_summary_group,
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_metrics_summary_group.csv"
  )
)

write_csv(
  prediction_smearing_factors,
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_smearing_factors.csv"
  )
)

saveRDS(
  list(
    gee = gee_prediction_results,
    lme_ri = lme_ri_prediction_results,
    lme_rirs =
      lme_rirs_prediction_results
  ),
  file.path(
    RDS_DIR,
    "prediction_model_results_GEE_LME_type2.rds"
  )
)

prediction_fit_qc <- map_dfr(
  OUTCOMES,
  function(outcome_name) {
    tibble(
      outcome = outcome_name,
      model = c(
        "GEE",
        "LME_RI",
        "LME_RI_RS"
      ),
      fit_ok = c(
        !is.null(
          gee_prediction_results[[outcome_name]][[1]]$fit
        ),
        !is.null(
          lme_ri_prediction_results[[outcome_name]][[1]]$fit
        ),
        !is.null(
          lme_rirs_prediction_results[[outcome_name]][[1]]$fit
        )
      ),
      prediction_ok = c(
        !is.null(
          gee_prediction_results[[outcome_name]][[1]]$pred
        ) &&
          nrow(
            gee_prediction_results[[outcome_name]][[1]]$pred
          ) > 0L,
        !is.null(
          lme_ri_prediction_results[[outcome_name]][[1]]$pred
        ) &&
          nrow(
            lme_ri_prediction_results[[outcome_name]][[1]]$pred
          ) > 0L,
        !is.null(
          lme_rirs_prediction_results[[outcome_name]][[1]]$pred
        ) &&
          nrow(
            lme_rirs_prediction_results[[outcome_name]][[1]]$pred
          ) > 0L
      )
    )
  }
)

stat_target_model_qc <-
  pred_GEE_LME_individual %>%
  mutate(
    outcome = as.character(outcome),
    model = as.character(model),
    Group = as.character(Group)
  ) %>%
  group_by(
    .imp,
    outcome,
    row_id,
    ID,
    Group
  ) %>%
  summarise(
    n_models = n_distinct(model),
    y_true_min = min(y_true),
    y_true_max = max(y_true),
    all_target_originally_observed =
      all(
        target_originally_observed %in%
          TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    y_true_range =
      y_true_max - y_true_min,
    passed = (
      n_models ==
        length(MODEL_LEVELS_ALL) &
        is.finite(y_true_range) &
        abs(y_true_range) <= 1e-10 &
        all_target_originally_observed
    )
  )

stat_prediction_duplicate_qc <-
  pred_GEE_LME_individual %>%
  count(
    .imp,
    outcome,
    model,
    row_id,
    name = "n"
  ) %>%
  filter(n != 1L)

write_csv(
  prediction_fit_qc,
  file.path(
    TABLE_DIR,
    "prediction_model_fit_qc.csv"
  )
)

write_csv(
  stat_target_model_qc,
  file.path(
    TABLE_DIR,
    "GEE_LME_target_model_consistency_qc.csv"
  )
)

write_csv(
  stat_prediction_duplicate_qc,
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_duplicate_qc.csv"
  )
)

smearing_factor_qc <-
  prediction_smearing_factors %>%
  mutate(
    passed = (
      is.finite(smearing_factor) &
        smearing_factor > 0 &
        is.finite(train_observation_n) &
        train_observation_n > 0 &
        is.finite(train_patient_n) &
        train_patient_n > 0 &
        backtransform_method ==
          "empirical_Duan_smearing_train_residuals"
    )
  )

write_csv(
  smearing_factor_qc,
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_smearing_factor_qc.csv"
  )
)

if (
  any(!prediction_fit_qc$fit_ok) ||
    any(!prediction_fit_qc$prediction_ok) ||
    any(!stat_target_model_qc$passed) ||
    nrow(stat_prediction_duplicate_qc) > 0L ||
    any(!smearing_factor_qc$passed)
) {
  stop(
    "GEE/LME 預測、target、重複列或 empirical smearing factor QC 未通過。",
    call. = FALSE
  )
}

# ============================================================
# 9. 預測圖
# ============================================================

pred_mean_plot_df <-
  pred_GEE_LME_individual %>%
  mutate(
    outcome = factor(
      outcome,
      levels = OUTCOMES
    ),
    outcome_label = OUTCOME_LABELS[
      as.character(outcome)
    ],
    model = factor(
      model,
      levels = MODEL_LEVELS_ALL
    ),
    model_label = MODEL_LABELS[
      as.character(model)
    ],
    Group = factor(
      Group,
      levels = c(
        "Control",
        "Statin"
      )
    ),
    Group_label = factor(
      GROUP_LABELS_ZH[
        as.character(Group)
      ],
      levels = c(
        "控制組",
        "史他汀組"
      )
    ),
    Time_Category = factor(
      Time_Category,
      levels = paste0("T", 1:19),
      ordered = TRUE
    ),
    Time_Num = as.integer(
      str_remove(
        as.character(Time_Category),
        "T"
      )
    )
  ) %>%
  group_by(
    outcome,
    outcome_label,
    model,
    model_label,
    Group,
    Group_label,
    Time_Category,
    Time_Num
  ) %>%
  summarise(
    true_mean =
      mean(y_true, na.rm = TRUE),
    pred_mean =
      mean(y_pred, na.rm = TRUE),
    true_sd = safe_sd(y_true),
    pred_sd = safe_sd(y_pred),
    Time_mid =
      mean(Time_in_Years, na.rm = TRUE),
    n_target = n_distinct(row_id),
    dataset_n = n_distinct(.imp),
    .groups = "drop"
  ) %>%
  mutate(
    true_se =
      true_sd / sqrt(n_target),
    pred_se =
      pred_sd / sqrt(n_target),
    t_crit = qt(
      0.975,
      df = pmax(
        n_target - 1L,
        1L
      )
    ),
    true_ci_lower =
      true_mean - t_crit * true_se,
    true_ci_upper =
      true_mean + t_crit * true_se,
    pred_ci_lower =
      pred_mean - t_crit * pred_se,
    pred_ci_upper =
      pred_mean + t_crit * pred_se,
    plot_eligible =
      n_target >= MIN_PLOT_TARGET_N,
    true_mean_plot = if_else(
      plot_eligible,
      true_mean,
      NA_real_
    ),
    pred_mean_plot = if_else(
      plot_eligible,
      pred_mean,
      NA_real_
    )
  ) %>%
  select(-t_crit)

write_csv(
  pred_mean_plot_df,
  file.path(
    TABLE_DIR,
    "GEE_LME_prediction_plot_cell_counts.csv"
  )
)

pred_mean_plot_long <-
  pred_mean_plot_df %>%
  select(
    outcome,
    outcome_label,
    model,
    model_label,
    Group,
    Group_label,
    Time_Category,
    Time_Num,
    Time_mid,
    n_target,
    dataset_n,
    plot_eligible,
    true_mean_plot,
    pred_mean_plot
  ) %>%
  pivot_longer(
    cols = c(
      true_mean_plot,
      pred_mean_plot
    ),
    names_to = "series",
    values_to = "mean_value"
  ) %>%
  mutate(
    series_label = recode(
      series,
      true_mean_plot =
        "測試集實際平均值",
      pred_mean_plot =
        "模型預測平均值"
    )
  )

plot_traj_panel <- function(
    model_to_plot,
    outcome_to_plot,
    continuous_time = FALSE,
    lang = c("zh", "en")
) {
  lang <- match.arg(lang)

  plot_df <- pred_mean_plot_long %>%
    filter(
      as.character(model) ==
        model_to_plot,
      as.character(outcome) ==
        outcome_to_plot
    ) %>%
    mutate(
      group_plot = if (
        lang == "zh"
      ) {
        factor(
          Group_label,
          levels = c(
            "控制組",
            "史他汀組"
          )
        )
      } else {
        factor(
          GROUP_LABELS_EN[
            as.character(Group)
          ],
          levels = c(
            "Control group",
            "Statin group"
          )
        )
      },
      series_plot = if (
        lang == "zh"
      ) {
        factor(
          series_label,
          levels = c(
            "測試集實際平均值",
            "模型預測平均值"
          )
        )
      } else {
        factor(
          recode(
            as.character(series_label),
            `測試集實際平均值` =
              "Observed test mean",
            `模型預測平均值` =
              "Model-predicted mean"
          ),
          levels = c(
            "Observed test mean",
            "Model-predicted mean"
          )
        )
      }
    )

  x_var <- if (
    continuous_time
  ) {
    "Time_mid"
  } else {
    "Time_Num"
  }

  p <- ggplot(
    plot_df,
    aes(
      x = .data[[x_var]],
      y = mean_value,
      color = group_plot,
      shape = series_plot
    )
  ) +
    geom_point(
      size = 1.85,
      alpha = 0.90,
      stroke = 0.70,
      na.rm = TRUE
    ) +
    scale_shape_manual(
      values = if (lang == "zh") {
        c(
          "測試集實際平均值" = 1,
          "模型預測平均值" = 15
        )
      } else {
        c(
          "Observed test mean" = 1,
          "Model-predicted mean" = 15
        )
      }
    ) +
    theme_bw(base_size = 9.5) +
    labs(
      title = paste0(
        MODEL_LABELS[[model_to_plot]],
        if (lang == "zh") "：" else ": ",
        OUTCOME_LABELS[[outcome_to_plot]]
      ),
      x = if (continuous_time) {
        if (lang == "zh") {
          "追蹤區間平均時間（年）"
        } else {
          "Mean follow-up time (years)"
        }
      } else {
        if (lang == "zh") {
          "追蹤時間區間"
        } else {
          "Follow-up interval"
        }
      },
      y = if (lang == "zh") {
        paste0(
          OUTCOME_LABELS[[outcome_to_plot]],
          " 平均值"
        )
      } else {
        paste0(
          "Mean ",
          OUTCOME_LABELS[[outcome_to_plot]]
        )
      },
      color = if (lang == "zh") {
        "組別"
      } else {
        "Group"
      },
      shape = if (lang == "zh") {
        "資料類型"
      } else {
        "Series"
      }
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10.5
      ),
      legend.position = "bottom",
      axis.text.x =
        element_text(size = 7.3),
      axis.text.y =
        element_text(size = 8)
    )

  if (!continuous_time) {
    p <- p +
      scale_x_continuous(
        breaks = seq(
          1,
          19,
          by = 2
        ),
        labels = paste0(
          "T",
          seq(1, 19, by = 2)
        ),
        limits = c(0.5, 19.5)
      )
  }

  p
}

save_panel_plot <- function(
    models,
    file_name,
    continuous_time = FALSE,
    lang = c("zh", "en")
) {
  lang <- match.arg(lang)
  
  panel_list <- list()
  
  for (m in models) {
    for (o in OUTCOMES) {
      panel_list[[paste(m, o, sep = "_")]] <-
        plot_traj_panel(
          model_to_plot = m,
          outcome_to_plot = o,
          continuous_time = continuous_time,
          lang = lang
        )
    }
  }
  
  # 八格主圖
  panel_grid <- patchwork::wrap_plots(
    panel_list,
    ncol = 4,
    guides = "collect"
  )
  
  # 將共同圖例放在下方
  panel_grid <- panel_grid &
    ggplot2::theme(
      legend.position = "bottom"
    )
  
  # 標題文字
  title_text <- if (lang == "zh") {
    paste0(
      "GEE 與 LME：測試集最後一次原始觀測目標之",
      "實際平均值與預測平均值\n",
      if (continuous_time) {
        "（依區間平均追蹤時間呈現；僅顯示 n ≥ 5 之區間）"
      } else {
        "（依追蹤時間區間呈現；僅顯示 n ≥ 5 之區間）"
      }
    )
  } else {
    paste0(
      "Observed and predicted means at the last originally observed test target\n",
      if (continuous_time) {
        "(By mean follow-up time; intervals with n ≥ 5 only)"
      } else {
        "(By follow-up interval; intervals with n ≥ 5 only)"
      }
    )
  }
  
  # 使用 grid grob 建立獨立標題列，
  # 完全避開 plot_annotation() 的版本相容警告
  title_panel <- patchwork::wrap_elements(
    full = grid::textGrob(
      label = title_text,
      x = 0.015,
      y = 0.50,
      hjust = 0,
      vjust = 0.5,
      gp = grid::gpar(
        fontsize = if (lang == "zh") 13 else 12.5,
        fontface = "plain",
        lineheight = 1.05
      )
    )
  )
  
  # 標題列與八格圖上下組合
  p_final <- (
    title_panel /
      panel_grid
  ) +
    patchwork::plot_layout(
      heights = c(0.09, 1)
    )
  
  ggsave(
    filename = file.path(
      FIGURE_DIR,
      file_name
    ),
    plot = p_final,
    width = 11.4,
    height = 7.6,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    bg = "white"
  )
  
  invisible(p_final)
}

for (lang in c("zh", "en")) {
  save_panel_plot(
    models = MODEL_LEVELS_MAIN,
    file_name = paste0(
      "GEE_LME_observed_vs_predicted_8panel_interval_",
      lang,
      ".png"
    ),
    continuous_time = FALSE,
    lang = lang
  )
  
  save_panel_plot(
    models = MODEL_LEVELS_MAIN,
    file_name = paste0(
      "GEE_LME_observed_vs_predicted_8panel_continuous_time_",
      lang,
      ".png"
    ),
    continuous_time = TRUE,
    lang = lang
  )
}

# ============================================================
# 10. 讀取 DL 結果與共同 target QC
# ============================================================

DL_INDIVIDUAL_PATH <-
  pick_existing_path(
    c(
      file.path(
        "analysis_outputs",
        "dl_output",
        "tables",
        "dl_test_predictions_individual.csv"
      ),
      file.path(
        "dl_output",
        "tables",
        "dl_test_predictions_individual.csv"
      ),
      "dl_test_predictions_individual.csv"
    ),
    "找不到最新版 DL individual prediction 檔。"
  )

DL_SAMPLE_SIZE_PATH <-
  pick_existing_path(
    c(
      file.path(
        "analysis_outputs",
        "dl_output",
        "tables",
        "dl_sequence_sample_size.csv"
      ),
      file.path(
        "dl_output",
        "tables",
        "dl_sequence_sample_size.csv"
      ),
      "dl_sequence_sample_size.csv"
    ),
    "找不到最新版 DL sequence sample size 檔。"
  )

dl_prediction_individual <- read_csv(
  DL_INDIVIDUAL_PATH,
  show_col_types = FALSE
)

dl_sequence_sample_size <- read_csv(
  DL_SAMPLE_SIZE_PATH,
  show_col_types = FALSE
)

required_dl_individual_cols <- c(
  ".imp",
  "row_id",
  "ID",
  "Group",
  "outcome",
  "model",
  "y_true",
  "y_pred",
  "target_originally_observed"
)

required_dl_sample_cols <- c(
  ".imp",
  "model",
  "outcome",
  "train_target_n",
  "validation_target_n",
  "test_target_n"
)

missing_dl_individual_cols <- setdiff(
  required_dl_individual_cols,
  names(dl_prediction_individual)
)

missing_dl_sample_cols <- setdiff(
  required_dl_sample_cols,
  names(dl_sequence_sample_size)
)

if (
  length(missing_dl_individual_cols) >
  0L
) {
  stop(
    "DL individual prediction 缺少必要欄位：",
    paste(
      missing_dl_individual_cols,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  length(missing_dl_sample_cols) >
  0L
) {
  stop(
    "DL sequence sample size 缺少必要欄位：",
    paste(
      missing_dl_sample_cols,
      collapse = ", "
    ),
    call. = FALSE
  )
}

dl_sequence_sample_size <-
  dl_sequence_sample_size %>%
  mutate(
    .imp = as.integer(.imp),
    model = as.character(model),
    outcome = as.character(outcome)
  )

unknown_dl_models_in_sample_size <-
  setdiff(
    unique(
      dl_sequence_sample_size$model
    ),
    EXPECTED_DL_MODELS
  )

if (
  length(
    unknown_dl_models_in_sample_size
  ) > 0L
) {
  stop(
    "DL sequence sample size 含未預期模型：",
    paste(
      unknown_dl_models_in_sample_size,
      collapse = ", "
    ),
    call. = FALSE
  )
}

dl_target_sample_model_qc <-
  dl_sequence_sample_size %>%
  select(
    all_of(required_dl_sample_cols)
  ) %>%
  pivot_longer(
    cols = ends_with("_target_n"),
    names_to = "split",
    values_to = "target_n_DL"
  ) %>%
  mutate(
    split = str_remove(
      split,
      "_target_n"
    ),
    target_n_DL =
      as.integer(target_n_DL)
  ) %>%
  group_by(
    .imp,
    outcome,
    split
  ) %>%
  summarise(
    n_models = n_distinct(model),
    target_n_min =
      min(target_n_DL),
    target_n_max =
      max(target_n_DL),
    passed = (
      n_models ==
        length(EXPECTED_DL_MODELS) &
        target_n_min ==
        target_n_max
    ),
    .groups = "drop"
  )

target_sample_alignment_all_splits <-
  target_sample_size_by_split %>%
  transmute(
    .imp = as.integer(.imp),
    outcome = as.character(outcome),
    split = as.character(split),
    target_n_R =
      as.integer(target_patient_n)
  ) %>%
  full_join(
    dl_target_sample_model_qc %>%
      transmute(
        .imp,
        outcome,
        split,
        target_n_DL =
          as.integer(target_n_min),
        dl_model_count = n_models,
        dl_models_consistent = passed
      ),
    by = c(
      ".imp",
      "outcome",
      "split"
    )
  ) %>%
  mutate(
    passed = (
      !is.na(target_n_R) &
        !is.na(target_n_DL) &
        dl_models_consistent %in% TRUE &
        target_n_R == target_n_DL
    )
  )

write_csv(
  dl_target_sample_model_qc,
  file.path(
    TABLE_DIR,
    "DL_target_sample_size_model_qc.csv"
  )
)

write_csv(
  target_sample_alignment_all_splits,
  file.path(
    TABLE_DIR,
    "GEE_LME_DL_target_sample_size_alignment_all_splits.csv"
  )
)

if (
  any(!dl_target_sample_model_qc$passed) ||
    any(
      !target_sample_alignment_all_splits$
        passed
    )
) {
  stop(
    "R 與 DL 的 Train/Validation/Test target 樣本數不一致。",
    call. = FALSE
  )
}

dl_prediction_individual <-
  dl_prediction_individual %>%
  mutate(
    .imp = as.integer(.imp),
    row_id = as.integer(row_id),
    ID = as.character(ID),
    Group = as.character(Group),
    outcome = as.character(outcome),
    model = as.character(model),
    y_true = as.numeric(y_true),
    y_pred = as.numeric(y_pred),
    target_originally_observed =
      parse_logical_safe(
        target_originally_observed
      )
  )

unknown_dl_models <- setdiff(
  unique(dl_prediction_individual$model),
  EXPECTED_DL_MODELS
)

if (length(unknown_dl_models) > 0L) {
  stop(
    "DL individual prediction 含未預期模型：",
    paste(
      unknown_dl_models,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  any(
    !(
      dl_prediction_individual$
        target_originally_observed %in%
        TRUE
    )
  )
) {
  stop(
    "DL individual prediction 含非原始觀測 target。",
    call. = FALSE
  )
}

if (
  any(
    !is.finite(
      dl_prediction_individual$y_true
    )
  ) ||
    any(
      !is.finite(
        dl_prediction_individual$y_pred
      )
    ) ||
    any(
      dl_prediction_individual$y_true <=
        0
    )
) {
  stop(
    "DL individual prediction 含非有限真值、非有限預測值或非正真值。",
    call. = FALSE
  )
}

dl_target_model_consistency_qc <-
  dl_prediction_individual %>%
  group_by(
    .imp,
    outcome,
    row_id,
    ID,
    Group
  ) %>%
  summarise(
    n_models = n_distinct(model),
    y_true_min = min(y_true),
    y_true_max = max(y_true),
    all_targets_observed =
      all(
        target_originally_observed %in%
          TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    y_true_range =
      y_true_max - y_true_min,
    passed = (
      n_models ==
        length(EXPECTED_DL_MODELS) &
        is.finite(y_true_range) &
        abs(y_true_range) <= 1e-10 &
        all_targets_observed
    )
  )

dl_target_model_count_qc <-
  dl_prediction_individual %>%
  count(
    .imp,
    outcome,
    model,
    name = "target_n"
  ) %>%
  group_by(
    .imp,
    outcome
  ) %>%
  summarise(
    n_models = n_distinct(model),
    target_n_min = min(target_n),
    target_n_max = max(target_n),
    passed = (
      n_models ==
        length(EXPECTED_DL_MODELS) &
        target_n_min ==
        target_n_max
    ),
    .groups = "drop"
  )

dl_prediction_duplicate_qc <-
  dl_prediction_individual %>%
  count(
    .imp,
    outcome,
    model,
    row_id,
    name = "n"
  ) %>%
  filter(n != 1L)

write_csv(
  dl_target_model_consistency_qc,
  file.path(
    TABLE_DIR,
    "DL_target_row_truth_consistency_qc.csv"
  )
)

write_csv(
  dl_target_model_count_qc,
  file.path(
    TABLE_DIR,
    "DL_target_model_count_consistency_qc.csv"
  )
)

write_csv(
  dl_prediction_duplicate_qc,
  file.path(
    TABLE_DIR,
    "DL_prediction_duplicate_qc.csv"
  )
)

if (
  any(
    !dl_target_model_consistency_qc$
      passed
  ) ||
    any(
      !dl_target_model_count_qc$passed
    ) ||
    nrow(dl_prediction_duplicate_qc) >
    0L
) {
  stop(
    "LSTM-NARX、TCN 與 ODE-LSTM 使用的 test target、病患、真值或預測列不一致。",
    call. = FALSE
  )
}

r_targets <- pred_GEE_LME_individual %>%
  filter(
    as.character(model) == "GEE"
  ) %>%
  transmute(
    .imp = as.integer(.imp),
    outcome = as.character(outcome),
    row_id = as.integer(row_id),
    ID = as.character(ID),
    Group = as.character(Group),
    y_true_R = as.numeric(y_true)
  ) %>%
  distinct()

dl_targets_all_models <-
  dl_prediction_individual %>%
  transmute(
    .imp = as.integer(.imp),
    outcome = as.character(outcome),
    model = as.character(model),
    row_id = as.integer(row_id),
    ID = as.character(ID),
    Group = as.character(Group),
    y_true_DL = as.numeric(y_true)
  ) %>%
  distinct()

target_alignment_qc <-
  dl_targets_all_models %>%
  full_join(
    r_targets,
    by = c(
      ".imp",
      "outcome",
      "row_id",
      "ID",
      "Group"
    )
  ) %>%
  mutate(
    alignment_status = case_when(
      is.na(y_true_R) ~ "DL_only",
      is.na(y_true_DL) ~ "R_only",
      TRUE ~ "matched"
    ),
    absolute_truth_difference =
      abs(y_true_R - y_true_DL),
    tolerance = 1e-6 * pmax(
      1,
      abs(y_true_R),
      abs(y_true_DL),
      na.rm = TRUE
    ),
    passed = (
      alignment_status == "matched" &
        is.finite(
          absolute_truth_difference
        ) &
        absolute_truth_difference <=
        tolerance
    )
  )

target_alignment_summary <-
  target_alignment_qc %>%
  count(
    outcome,
    model,
    alignment_status,
    passed,
    name = "n"
  )

write_csv(
  target_alignment_qc,
  file.path(
    TABLE_DIR,
    "GEE_LME_DL_target_alignment_qc.csv"
  )
)

write_csv(
  target_alignment_summary,
  file.path(
    TABLE_DIR,
    "GEE_LME_DL_target_alignment_summary.csv"
  )
)

if (any(!target_alignment_qc$passed)) {
  stop(
    "GEE/LME 與至少一個 DL 模型的 test target、病患或真值不一致。",
    call. = FALSE
  )
}

# ============================================================
# 11. 五模型 bootstrap
# ============================================================

bootstrap_metric_ci_one <- function(
    pred_df,
    group_eval,
    B = BOOTSTRAP_B,
    seed = BOOTSTRAP_SEED
) {
  d <- if (
    group_eval == "Overall"
  ) {
    pred_df
  } else {
    pred_df %>%
      filter(
        Group == group_eval
      )
  }

  d <- d %>%
    filter(
      is.finite(y_true),
      is.finite(y_pred),
      y_true > 0
    ) %>%
    distinct(
      ID,
      row_id,
      .keep_all = TRUE
    )

  point <- calc_pred_metrics(
    d$y_true,
    d$y_pred
  )

  if (
    nrow(d) < 2L ||
      B < 2L
  ) {
    return(
      tibble(
        RMSE = point$RMSE,
        RMSE_lower = NA_real_,
        RMSE_upper = NA_real_,
        MAE = point$MAE,
        MAE_lower = NA_real_,
        MAE_upper = NA_real_,
        MAPE = point$MAPE,
        MAPE_lower = NA_real_,
        MAPE_upper = NA_real_,
        Bias = point$Bias,
        Bias_lower = NA_real_,
        Bias_upper = NA_real_,
        n_test = point$n_test,
        bootstrap_B = B
      )
    )
  }

  set.seed(seed)

  group_indices <- if (
    group_eval == "Overall"
  ) {
    split(
      seq_len(nrow(d)),
      d$Group
    )
  } else {
    list(
      group = seq_len(nrow(d))
    )
  }

  boot_mat <- matrix(
    NA_real_,
    nrow = B,
    ncol = 4L,
    dimnames = list(
      NULL,
      c(
        "RMSE",
        "MAE",
        "MAPE",
        "Bias"
      )
    )
  )

  for (b in seq_len(B)) {
    idx <- unlist(
      lapply(
        group_indices,
        function(ii) {
          sample(
            ii,
            length(ii),
            replace = TRUE
          )
        }
      ),
      use.names = FALSE
    )

    boot_mat[b, ] <- unlist(
      calc_pred_metrics(
        d$y_true[idx],
        d$y_pred[idx]
      )[
        c(
          "RMSE",
          "MAE",
          "MAPE",
          "Bias"
        )
      ]
    )
  }

  q <- apply(
    boot_mat,
    2,
    quantile,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 6
  )

  tibble(
    RMSE = point$RMSE,
    RMSE_lower = q[1, "RMSE"],
    RMSE_upper = q[2, "RMSE"],
    MAE = point$MAE,
    MAE_lower = q[1, "MAE"],
    MAE_upper = q[2, "MAE"],
    MAPE = point$MAPE,
    MAPE_lower = q[1, "MAPE"],
    MAPE_upper = q[2, "MAPE"],
    Bias = point$Bias,
    Bias_lower = q[1, "Bias"],
    Bias_upper = q[2, "Bias"],
    n_test = point$n_test,
    bootstrap_B = B
  )
}

summarise_individual_predictions_with_bootstrap <-
  function(pred_df) {
    pred_df %>%
      mutate(
        outcome = as.character(outcome),
        model = as.character(model),
        Group = as.character(Group)
      ) %>%
      group_split(
        outcome,
        model,
        .keep = TRUE
      ) %>%
      map_dfr(
        function(dm) {
          outcome_name <-
            unique(dm$outcome)

          model_name <-
            unique(dm$model)

          map_dfr(
            GROUP_LEVELS,
            function(g) {
              seed_offset <- sum(
                utf8ToInt(
                  paste(
                    outcome_name,
                    model_name,
                    g
                  )
                )
              )

              bootstrap_metric_ci_one(
                dm,
                g,
                seed =
                  BOOTSTRAP_SEED +
                  seed_offset
              ) %>%
                mutate(
                  outcome =
                    outcome_name,
                  outcome_label =
                    OUTCOME_LABELS[[outcome_name]],
                  model =
                    model_name,
                  Group_eval = g,
                  Group_label =
                    GROUP_LABELS_ZH[[g]],
                  dataset_n = 1L,
                  .before = 1
                )
            }
          )
        }
      ) %>%
      mutate(
        outcome = factor(
          outcome,
          levels = OUTCOMES
        ),
        model = factor(
          model,
          levels = MODEL_LEVELS_FIVE
        ),
        Group_eval = factor(
          Group_eval,
          levels = GROUP_LEVELS
        )
      ) %>%
      arrange(
        outcome,
        model,
        Group_eval
      )
  }

stat_predictions_for_five_models <-
  pred_GEE_LME_individual %>%
  filter(
    as.character(model) %in%
      c(
        "GEE",
        "LME_RI_RS"
      )
  ) %>%
  mutate(
    model = recode(
      as.character(model),
      GEE = "GEE",
      LME_RI_RS = "LME"
    )
  ) %>%
  select(
    .imp,
    row_id,
    ID,
    Group,
    outcome,
    model,
    y_true,
    y_pred,
    target_originally_observed
  )

dl_predictions_for_five_models <-
  dl_prediction_individual %>%
  select(
    .imp,
    row_id,
    ID,
    Group,
    outcome,
    model,
    y_true,
    y_pred,
    target_originally_observed
  )

five_model_predictions_individual <-
  bind_rows(
    stat_predictions_for_five_models,
    dl_predictions_for_five_models
  ) %>%
  mutate(
    .imp = as.integer(.imp),
    row_id = as.integer(row_id),
    ID = as.character(ID),
    Group = as.character(Group),
    outcome = as.character(outcome),
    model = as.character(model),
    y_true = as.numeric(y_true),
    y_pred = as.numeric(y_pred),
    target_originally_observed =
      parse_logical_safe(
        target_originally_observed
      )
  )

five_model_target_qc <-
  five_model_predictions_individual %>%
  group_by(
    .imp,
    outcome,
    row_id,
    ID,
    Group
  ) %>%
  summarise(
    n_models = n_distinct(model),
    y_true_min = min(y_true),
    y_true_max = max(y_true),
    all_target_originally_observed =
      all(
        target_originally_observed %in%
          TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    y_true_range =
      y_true_max - y_true_min,
    passed = (
      n_models ==
        length(MODEL_LEVELS_FIVE) &
        is.finite(y_true_range) &
        abs(y_true_range) <= 1e-10 &
        all_target_originally_observed
    )
  )

five_model_duplicate_qc <-
  five_model_predictions_individual %>%
  count(
    .imp,
    outcome,
    model,
    row_id,
    name = "n"
  ) %>%
  filter(n != 1L)

five_model_model_count_qc <-
  five_model_predictions_individual %>%
  count(
    .imp,
    outcome,
    model,
    name = "target_n"
  ) %>%
  group_by(
    .imp,
    outcome
  ) %>%
  summarise(
    n_models = n_distinct(model),
    target_n_min = min(target_n),
    target_n_max = max(target_n),
    passed = (
      n_models ==
        length(MODEL_LEVELS_FIVE) &
        target_n_min ==
        target_n_max
    ),
    .groups = "drop"
  )

write_csv(
  five_model_target_qc,
  file.path(
    TABLE_DIR,
    "five_model_target_consistency_qc.csv"
  )
)

write_csv(
  five_model_duplicate_qc,
  file.path(
    TABLE_DIR,
    "five_model_prediction_duplicate_qc.csv"
  )
)

write_csv(
  five_model_model_count_qc,
  file.path(
    TABLE_DIR,
    "five_model_target_count_consistency_qc.csv"
  )
)

if (
  any(!five_model_target_qc$passed) ||
    nrow(five_model_duplicate_qc) >
    0L ||
    any(
      !five_model_model_count_qc$passed
    )
) {
  stop(
    "五模型 target、真值、樣本數或個體預測列不一致。",
    call. = FALSE
  )
}

five_model_prediction_summary_group <-
  summarise_individual_predictions_with_bootstrap(
    five_model_predictions_individual
  )

five_model_prediction_table_for_paper <-
  five_model_prediction_summary_group %>%
  transmute(
    outcome = as.character(outcome),
    outcome_label,
    model = as.character(model),
    Group_eval =
      as.character(Group_eval),
    Group_label,
    n_test,
    RMSE,
    RMSE_95CI = sprintf(
      "%.3f (%.3f, %.3f)",
      RMSE,
      RMSE_lower,
      RMSE_upper
    ),
    MAE,
    MAE_95CI = sprintf(
      "%.3f (%.3f, %.3f)",
      MAE,
      MAE_lower,
      MAE_upper
    ),
    MAPE,
    MAPE_95CI = sprintf(
      "%.3f (%.3f, %.3f)",
      MAPE,
      MAPE_lower,
      MAPE_upper
    ),
    Bias,
    Bias_95CI = sprintf(
      "%.3f (%.3f, %.3f)",
      Bias,
      Bias_lower,
      Bias_upper
    ),
    bootstrap_B
  )

write_csv(
  five_model_predictions_individual,
  file.path(
    TABLE_DIR,
    "five_model_test_predictions_individual.csv"
  )
)

write_csv(
  five_model_prediction_summary_group,
  file.path(
    TABLE_DIR,
    "five_model_prediction_metrics_with_bootstrap_CI.csv"
  )
)

write_csv(
  five_model_prediction_table_for_paper,
  file.path(
    TABLE_DIR,
    "five_model_prediction_table_for_paper.csv"
  )
)

# ============================================================
# 12. 最終 QC 與 Excel
# ============================================================

final_analysis_qc <- bind_rows(
  qc_row(
    "input_qc_all_passed",
    as.integer(all(input_qc$passed)),
    1L
  ),
  qc_row(
    "inference_prediction_row_alignment",
    sum(!row_manifest_alignment_qc$passed),
    0L
  ),
  qc_row(
    "inference_same_patient_same_time_no_duplicates",
    inference_same_time_duplicate_n,
    0L
  ),
  qc_row(
    "prediction_same_patient_same_time_no_duplicates",
    prediction_same_time_duplicate_n,
    0L
  ),
  qc_row(
    "inference_fit_count",
    nrow(inference_fit_qc),
    length(OUTCOMES) *
      length(MODEL_LEVELS_ALL) *
      length(imp_list)
  ),
  qc_row(
    "inference_fit_all_successful",
    as.integer(
      all(inference_fit_qc$fit_ok)
    ),
    1L
  ),
  qc_row(
    "lme_apvar_all_valid",
    as.integer(
      all(
        inference_fit_qc$apvar_ok[
          inference_fit_qc$model !=
            "GEE"
        ] %in% TRUE
      )
    ),
    1L
  ),
  qc_row(
    "pooled_parameters_all_10_imputations",
    as.integer(
      all_n_imp_equal(
        pooled_parameter_estimates$n_imp,
        length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "lme_aic_all_10_imputations",
    as.integer(
      all_n_imp_equal(
        lme_aic_comparison$n_imp,
        length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "lme_aic_pair_complete_all_10_imputations",
    as.integer(
      all(
        lme_aic_selection_summary$
          n_imp_complete ==
          length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "lme_aic_selection_reviewed",
    as.integer(
      isTRUE(AIC_SELECTION_REVIEWED)
    ),
    1L
  ),
  qc_row(
    "all_outcomes_mean_AIC_prefer_RIRS",
    as.integer(
      all_outcomes_mean_AIC_prefer_RIRS
    ),
    1L
  ),
  qc_row(
    "lme_variance_components_all_10_imputations",
    as.integer(
      all_n_imp_equal(
        lme_variance_components_summary_for_paper$
          n_imp,
        length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "lme_random_slope_variance_complete",
    as.integer(
      all(
        lme_random_structure_qc$
          slope_variance_qc_passed
      )
    ),
    1L
  ),
  qc_row(
    "lme_random_correlation_complete",
    as.integer(
      all(
        lme_random_structure_qc$
          correlation_qc_passed
      )
    ),
    1L
  ),
  qc_row(
    "lme_random_slope_not_near_zero",
    as.integer(
      all(
        !lme_random_structure_qc$
          any_slope_near_zero
      )
    ),
    1L,
    fatal = FALSE
  ),
  qc_row(
    "lme_random_correlation_not_near_boundary",
    as.integer(
      all(
        !lme_random_structure_qc$
          any_correlation_near_boundary
      )
    ),
    1L,
    fatal = FALSE
  ),
  qc_row(
    "gee_alpha_all_10_imputations",
    as.integer(
      all(
        gee_working_correlation_summary_for_paper$
          n_imp_alpha ==
          length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "gee_scale_all_10_imputations",
    as.integer(
      all(
        gee_working_correlation_summary_for_paper$
          n_imp_scale ==
          length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "lme_diagnostics_all_10_imputations",
    as.integer(
      all(
        lme_diagnostics_all_imputations_summary$
          n_imp ==
          length(imp_list)
      )
    ),
    1L
  ),
  qc_row(
    "prediction_fit_all_successful",
    as.integer(
      all(
        prediction_fit_qc$fit_ok &
          prediction_fit_qc$prediction_ok
      )
    ),
    1L
  ),
  qc_row(
    "empirical_smearing_all_models_valid",
    as.integer(
      all(smearing_factor_qc$passed)
    ),
    1L
  ),
  qc_row(
    "empirical_smearing_expected_row_count",
    nrow(prediction_smearing_factors),
    length(OUTCOMES) *
      length(MODEL_LEVELS_ALL) *
      length(imp_list_prediction)
  ),
  qc_row(
    "empirical_smearing_method_consistent",
    as.integer(
      all(
        prediction_smearing_factors$
          backtransform_method ==
          "empirical_Duan_smearing_train_residuals"
      )
    ),
    1L
  ),
  qc_row(
    "statistical_target_model_consistency",
    as.integer(
      all(stat_target_model_qc$passed)
    ),
    1L
  ),
  qc_row(
    "statistical_prediction_no_duplicates",
    as.integer(
      nrow(
        stat_prediction_duplicate_qc
      ) == 0L
    ),
    1L
  ),
  qc_row(
    "dl_target_sample_size_consistent_across_models",
    as.integer(
      all(dl_target_sample_model_qc$passed)
    ),
    1L
  ),
  qc_row(
    "gee_lme_dl_target_sample_size_alignment_all_splits",
    as.integer(
      all(
        target_sample_alignment_all_splits$
          passed
      )
    ),
    1L
  ),
  qc_row(
    "dl_target_row_truth_consistency",
    as.integer(
      all(
        dl_target_model_consistency_qc$
          passed
      )
    ),
    1L
  ),
  qc_row(
    "dl_target_count_consistency",
    as.integer(
      all(dl_target_model_count_qc$passed)
    ),
    1L
  ),
  qc_row(
    "dl_prediction_no_duplicates",
    as.integer(
      nrow(
        dl_prediction_duplicate_qc
      ) == 0L
    ),
    1L
  ),
  qc_row(
    "gee_lme_dl_target_alignment_all_models",
    as.integer(
      all(target_alignment_qc$passed)
    ),
    1L
  ),
  qc_row(
    "five_model_target_consistency",
    as.integer(
      all(five_model_target_qc$passed)
    ),
    1L
  ),
  qc_row(
    "five_model_target_count_consistency",
    as.integer(
      all(
        five_model_model_count_qc$passed
      )
    ),
    1L
  ),
  qc_row(
    "five_model_prediction_no_duplicates",
    as.integer(
      nrow(
        five_model_duplicate_qc
      ) == 0L
    ),
    1L
  ),
  qc_row(
    "five_model_summary_row_count",
    nrow(
      five_model_prediction_summary_group
    ),
    length(OUTCOMES) *
      length(MODEL_LEVELS_FIVE) *
      length(GROUP_LEVELS)
  ),
  qc_row(
    "five_model_bootstrap_ci_complete",
    as.integer(
      all(
        is.finite(
          as.matrix(
            five_model_prediction_summary_group %>%
              select(
                RMSE_lower,
                RMSE_upper,
                MAE_lower,
                MAE_upper,
                MAPE_lower,
                MAPE_upper,
                Bias_lower,
                Bias_upper
              )
          )
        )
      )
    ),
    1L
  )
)

write_csv(
  final_analysis_qc,
  file.path(
    TABLE_DIR,
    "final_model_analysis_qc.csv"
  )
)

fatal_qc_failed <-
  final_analysis_qc %>%
  filter(
    fatal,
    !passed
  )

if (nrow(fatal_qc_failed) > 0L) {
  print(fatal_qc_failed)

  stop(
    "資料分析最終必要 QC 未通過。請先處理上述項目，不可直接定稿。",
    call. = FALSE
  )
}

nonfatal_qc_flagged <-
  final_analysis_qc %>%
  filter(
    !fatal,
    !passed
  )

if (nrow(nonfatal_qc_flagged) > 0L) {
  warning(
    "部分 LME 隨機結構診斷出現警示；請檢查 lme_random_structure_qc.csv。"
  )

  print(nonfatal_qc_flagged)
}

output_sheets <- list(
  input_qc = input_qc,
  row_alignment_summary =
    row_manifest_alignment_summary_qc,
  row_alignment_detail =
    row_manifest_alignment_qc,
  split_patient = split_check_patient,
  split_observation =
    split_check_observation,
  prediction_scope =
    prediction_evaluation_scope,
  target_n_by_split =
    target_sample_size_summary,
  target_n_by_split_imp =
    target_sample_size_by_split,
  gee_lme_train_n =
    gee_lme_training_observation_count,
  inference_fit_qc =
    inference_fit_qc,
  pooled_parameters =
    pooled_parameter_estimates,
  gee_working_corr =
    gee_working_correlation_summary_for_paper,
  gee_working_corr_imp =
    gee_working_correlation_by_imp,
  lme_varcomp =
    lme_variance_components_summary_for_paper,
  lme_varcomp_imp =
    lme_variance_components_by_imp,
  lme_random_structure =
    lme_random_structure_qc,
  lme_aic_by_imp =
    lme_aic_by_imp,
  lme_aic_pairwise =
    lme_aic_pairwise_by_imp,
  lme_aic_summary =
    lme_aic_comparison,
  lme_aic_selection =
    lme_aic_selection_summary,
  lme_diagnostics_imp1 =
    lme_diagnostics_summary,
  lme_diagnostics_by_imp =
    lme_diagnostics_by_imp,
  lme_diagnostics_all =
    lme_diagnostics_all_imputations_summary,
  lme_random_effects =
    lme_random_effects_summary,
  smearing_factors =
    prediction_smearing_factors,
  smearing_qc =
    smearing_factor_qc,
  plot_cell_counts =
    pred_mean_plot_df,
  prediction_fit_qc =
    prediction_fit_qc,
  stat_target_qc =
    stat_target_model_qc,
  stat_duplicate_qc =
    stat_prediction_duplicate_qc,
  pred_metrics_single =
    prediction_metrics_GEE_LME_summary_group,
  dl_target_n_qc =
    dl_target_sample_model_qc,
  target_n_alignment =
    target_sample_alignment_all_splits,
  dl_target_row_qc =
    dl_target_model_consistency_qc,
  dl_target_count_qc =
    dl_target_model_count_qc,
  dl_duplicate_qc =
    dl_prediction_duplicate_qc,
  target_alignment =
    target_alignment_qc,
  target_align_summary =
    target_alignment_summary,
  five_model_target_qc =
    five_model_target_qc,
  five_model_count_qc =
    five_model_model_count_qc,
  five_model_dup_qc =
    five_model_duplicate_qc,
  five_model_bootstrap =
    five_model_prediction_summary_group,
  five_model_paper =
    five_model_prediction_table_for_paper,
  final_qc =
    final_analysis_qc
)

write_xlsx(
  output_sheets,
  file.path(
    TABLE_DIR,
    "GEE_LME_model_outputs_type2.xlsx"
  )
)

message("資料分析.R 已執行完成。")
message("程式版本：", SCRIPT_VERSION)
message("輸出資料夾：", OUT_DIR)
message("請優先檢查：")
message("1. model_analysis_input_qc.csv")
message("2. final_model_analysis_qc.csv")
message("3. lme_aic_selection_summary.csv")
message("4. lme_random_structure_qc.csv")
message("5. GEE_LME_prediction_smearing_factors.csv")
message("6. GEE_LME_prediction_metrics_summary_group.csv")
message("7. five_model_prediction_metrics_with_bootstrap_CI.csv")
message("8. five_model_prediction_table_for_paper.csv")
message("9. GEE_LME_DL_target_alignment_qc.csv")
message("10. GEE_LME_model_outputs_type2.xlsx")
