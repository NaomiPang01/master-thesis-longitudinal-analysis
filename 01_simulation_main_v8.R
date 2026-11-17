# ============================================================
# 模擬主程式（正式版）
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(MASS)
  library(readr)
  library(purrr)
  library(nlme)
  library(geepack)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(tibble)
})

source("02_functions_simulation_v8.R")
source("03_fit_lme_gee_simulation_v8.R")
source("04_summary_tables_simulation_v8.R")
source("05_plot_figures_simulation_v8.R")

# ============================================================
# 0. 使用者設定
# ============================================================

SEED_BASE <- 5014L
RUN_MODE <- "formal_all"

# pilot_all 與 formal_all 使用完全相同的模型及訓練設定，
# 並共用 simulation_outputs；正式執行時可由 S0 run 4 接續。
if (RUN_MODE == "pilot_all") {
  SCENARIOS_TO_RUN <- "S0"
  N_RUNS <- 3L
  RUN_GEE <- TRUE
  RUN_LME <- TRUE
  RUN_DL <- TRUE
  DL_EPOCHS <- 50L
  DL_PATIENCE <- 5L
} else if (RUN_MODE == "formal_all") {
  SCENARIOS_TO_RUN <- c("S0", "S1", "S2", "S3", "S4")
  N_RUNS <- 1000L
  RUN_GEE <- TRUE
  RUN_LME <- TRUE
  RUN_DL <- TRUE
  DL_EPOCHS <- 50L
  DL_PATIENCE <- 5L
} else {
  stop("未知的 RUN_MODE：", RUN_MODE, call. = FALSE)
}

RUN_START <- 1L
DL_MAX_LEN <- 7L
DL_BATCH_SIZE <- 32L
TARGET_N <- 1500L
SAVE_SIM_DATA <- FALSE
KEEP_S0_RUN1_FOR_DIAGNOSTICS <- TRUE
SKIP_EXISTING <- TRUE
USE_CURRENT_FIXED_EFFECTS <- FALSE
USE_CURRENT_VARCOMP <- FALSE

# TRUE：每完成一個情境後詢問是否繼續；FALSE：自動跑完全部情境。
ASK_BEFORE_NEXT_SCENARIO <- TRUE

SIM_ROOT <- "simulation_outputs"
FIXED_SKELETON_PATH <- file.path(
  SIM_ROOT,
  "design",
  "fixed_skeleton_1500_patients.csv"
)
SAMPLING_ALLOCATION_PATH <- file.path(
  SIM_ROOT,
  "design",
  "Group_split_sampling_allocation.csv"
)

# ------------------------------------------------------------
# Python 設定
# ------------------------------------------------------------

Sys.setenv(SIM_PYTHON = "C:/Users/NaomiP/anaconda3/python.exe")

find_python_executable <- function() {
  env_py <- Sys.getenv("SIM_PYTHON", unset = "")
  if (nzchar(env_py) && file.exists(env_py)) {
    return(list(
      cmd = normalizePath(env_py, winslash = "/", mustWork = TRUE),
      prefix_args = character(0),
      source = "SIM_PYTHON"
    ))
  }
  
  is_bad_python <- function(x) {
    !nzchar(x) || grepl("WindowsApps", x, ignore.case = TRUE)
  }
  
  py <- Sys.which("python")
  if (!is_bad_python(py)) {
    return(list(cmd = unname(py), prefix_args = character(0), source = "python"))
  }
  
  py3 <- Sys.which("python3")
  if (!is_bad_python(py3)) {
    return(list(cmd = unname(py3), prefix_args = character(0), source = "python3"))
  }
  
  py_launcher <- Sys.which("py")
  if (nzchar(py_launcher)) {
    return(list(cmd = unname(py_launcher), prefix_args = "-3", source = "py -3"))
  }
  
  list(cmd = NA_character_, prefix_args = character(0), source = "not_found")
}

PYTHON_INFO <- find_python_executable()
PYTHON_EXE <- PYTHON_INFO$cmd
PYTHON_PREFIX_ARGS <- PYTHON_INFO$prefix_args
PYTHON_SCRIPT <- "run_dl_models_simulation_v8.py"

message("RUN_MODE：", RUN_MODE)
message("SIM_ROOT：", SIM_ROOT)
message("Python source：", PYTHON_INFO$source)
message("Python executable：", ifelse(is.na(PYTHON_EXE), "not found", PYTHON_EXE))

if (RUN_DL) {
  if (is.na(PYTHON_EXE) || !nzchar(PYTHON_EXE)) {
    stop("RUN_DL = TRUE，但找不到可用的 Python。")
  }
  if (!file.exists(PYTHON_SCRIPT)) {
    stop("RUN_DL = TRUE，但找不到 Python 腳本：", PYTHON_SCRIPT)
  }
  
  version_out <- system2(
    PYTHON_EXE,
    args = c(PYTHON_PREFIX_ARGS, "--version"),
    stdout = TRUE,
    stderr = TRUE
  )
  message(paste(version_out, collapse = "\n"))
  
  package_check <- system2(
    PYTHON_EXE,
    args = c(
      PYTHON_PREFIX_ARGS,
      "-c",
      shQuote("import tensorflow, torch, pandas, numpy, sklearn; print('Python packages OK')")
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  package_status <- attr(package_check, "status")
  message(paste(package_check, collapse = "\n"))
  
  if (!is.null(package_status) && package_status != 0L) {
    stop("Python 套件檢查失敗，請確認目前 Anaconda environment。")
  }
}

# ============================================================
# 1. 尋找資料
# ============================================================

DATA_CANDIDATES <- c(
  file.path("analysis_outputs", "dl_data", "prediction_long_with_split_type2.csv"),
  file.path("dl_data", "prediction_long_with_split_type2.csv"),
  "prediction_long_with_split_type2.csv"
)

POOLED_PARAM_PATH <- NA_character_
VARCOMP_PATH <- NA_character_

DATA_PATH <- find_first_existing(DATA_CANDIDATES)

if (is.na(DATA_PATH)) {
  stop("找不到 prediction_long_with_split_type2.csv。")
}

# ============================================================
# 2. empirical skeleton 與參數
# ============================================================

skeleton_all <- load_empirical_skeleton(
  data_path = DATA_PATH,
  imputation_id = 1,
  seed = SEED_BASE,
  create_split_if_missing = FALSE
)

dir.create(
  file.path(SIM_ROOT, "design"),
  recursive = TRUE,
  showWarnings = FALSE
)

# pilot 與 formal 必須使用完全相同的 1500 位病患。
# 已存在固定 skeleton 時直接讀取，不再重新抽樣。
if (file.exists(FIXED_SKELETON_PATH)) {
  message("讀取既有固定 skeleton：", FIXED_SKELETON_PATH)
  
  skeleton <- read_csv(
    FIXED_SKELETON_PATH,
    show_col_types = FALSE
  ) %>%
    mutate(
      row_id = as.integer(row_id),
      ID = as.character(ID),
      Group = as.character(Group),
      split = as.character(split),
      Time_in_Years = as.numeric(Time_in_Years),
      GroupStatin = ifelse(Group == "Statin", 1L, 0L),
      Sex1 = as.integer(Sex1)
    ) %>%
    arrange(ID, Time_in_Years, row_id)
  
  sampling_allocation <- if (
    file.exists(SAMPLING_ALLOCATION_PATH)
  ) {
    read_csv(
      SAMPLING_ALLOCATION_PATH,
      show_col_types = FALSE
    )
  } else {
    skeleton %>%
      distinct(ID, Group, split) %>%
      count(Group, split, name = "n_sample")
  }
} else {
  skeleton <- select_patients_by_group_split(
    skeleton = skeleton_all,
    target_n = TARGET_N,
    seed = SEED_BASE
  )
  
  sampling_allocation <- attr(
    skeleton,
    "sampling_allocation"
  )
  
  skeleton <- skeleton %>%
    arrange(ID, Time_in_Years, row_id) %>%
    mutate(
      ID = as.character(ID),
      Group = as.character(Group),
      split = as.character(split),
      GroupStatin = ifelse(Group == "Statin", 1L, 0L),
      Sex1 = as.integer(Sex1)
    )
  
  if (is.null(sampling_allocation)) {
    sampling_allocation <- skeleton %>%
      distinct(ID, Group, split) %>%
      count(Group, split, name = "n_sample")
  }
  
  write_csv(
    skeleton,
    FIXED_SKELETON_PATH
  )
  
  write_csv(
    sampling_allocation,
    SAMPLING_ALLOCATION_PATH
  )
}

SKELETON_PATIENT_N <- n_distinct(skeleton$ID)
SKELETON_ROW_N <- nrow(skeleton)

if (SKELETON_PATIENT_N != TARGET_N) {
  stop(
    "固定 skeleton 病患數不等於 TARGET_N：",
    SKELETON_PATIENT_N,
    " versus ",
    TARGET_N,
    call. = FALSE
  )
}

if (
  nrow(
    skeleton %>%
    count(row_id) %>%
    filter(n != 1L)
  ) > 0L
) {
  stop("固定 skeleton 的 row_id 不唯一。", call. = FALSE)
}

if (
  nrow(
    skeleton %>%
    distinct(ID, Group, split) %>%
    count(ID) %>%
    filter(n != 1L)
  ) > 0L
) {
  stop(
    "至少一位病患出現在多個 Group × split 階層。",
    call. = FALSE
  )
}

skeleton_qc <- summarise_skeleton(skeleton)

# T1-T19 是追蹤時間區間，不代表每位病患最多只能有 19 筆。
# 同一病患可在同一 Time_Category 內有多筆實際觀測。
time_category_qc <- skeleton %>%
  mutate(
    Time_Category_chr = as.character(Time_Category),
    Time_Category_num = suppressWarnings(
      as.integer(
        str_remove(
          Time_Category_chr,
          "^T"
        )
      )
    )
  )

if (
  any(
    is.na(time_category_qc$Time_Category_num) |
    time_category_qc$Time_Category_num < 1L |
    time_category_qc$Time_Category_num > 19L
  )
) {
  stop(
    "固定 skeleton 含 T1-T19 以外的 Time_Category。",
    call. = FALSE
  )
}

patient_time_category_qc <- time_category_qc %>%
  count(
    ID,
    Time_Category_chr,
    name = "observation_n_in_interval"
  )

time_category_summary <- tibble(
  patient_n = SKELETON_PATIENT_N,
  observation_n_total = SKELETON_ROW_N,
  maximum_observations_per_patient = max(
    count(skeleton, ID)$n
  ),
  maximum_unique_intervals_per_patient = max(
    time_category_qc %>%
      distinct(ID, Time_Category_chr) %>%
      count(ID) %>%
      pull(n)
  ),
  patient_n_with_more_than_19_observations = sum(
    count(skeleton, ID)$n > 19L
  ),
  patient_interval_cells_with_multiple_observations = sum(
    patient_time_category_qc$
      observation_n_in_interval > 1L
  ),
  maximum_observations_within_one_interval = max(
    patient_time_category_qc$
      observation_n_in_interval
  )
)

write_csv(
  sampling_allocation,
  SAMPLING_ALLOCATION_PATH
)
write_csv(
  skeleton_qc$patient_by_group,
  file.path(
    SIM_ROOT,
    "design",
    "skeleton_patient_by_group_split.csv"
  )
)
write_csv(
  skeleton_qc$observation_by_group,
  file.path(
    SIM_ROOT,
    "design",
    "skeleton_observation_by_group_split.csv"
  )
)
write_csv(
  skeleton_qc$repeat_summary %>%
    mutate(
      observation_n_total = SKELETON_ROW_N,
      target_patient_n = SKELETON_PATIENT_N
    ),
  file.path(
    SIM_ROOT,
    "design",
    "skeleton_repeated_measure_summary.csv"
  )
)
write_csv(
  time_category_summary,
  file.path(
    SIM_ROOT,
    "design",
    "skeleton_time_category_qc.csv"
  )
)

message("模擬 skeleton 病患數：", SKELETON_PATIENT_N)
message("模擬 skeleton 觀測筆數：", SKELETON_ROW_N)
message(
  "每位病患最大觀測筆數：",
  time_category_summary$
    maximum_observations_per_patient
)
message(
  "每位病患最大不同追蹤區間數：",
  time_category_summary$
    maximum_unique_intervals_per_patient
)
print(
  skeleton %>%
    distinct(ID, Group) %>%
    count(Group)
)
print(
  skeleton %>%
    distinct(ID, Group, split) %>%
    count(split, Group)
)
print(time_category_summary)

BASE_PARAMS <- get_base_parameters(
  pooled_path = POOLED_PARAM_PATH,
  varcomp_path = VARCOMP_PATH,
  use_current_fixed = USE_CURRENT_FIXED_EFFECTS,
  use_current_varcomp = USE_CURRENT_VARCOMP
)

base_param_table <- base_parameters_to_table(BASE_PARAMS)
print(base_param_table)

# ============================================================
# 3. 情境設定
# ============================================================

scenario_table <- tribble(
  ~scenario, ~scenario_label, ~sb_mult, ~sigma_mult, ~rho_mode,
  "S0", "基準情境",           1.0, 1.0, "original",
  "S1", "個體間變異較大",     1.5, 1.0, "original",
  "S2", "個體間變異較小",     0.5, 1.0, "original",
  "S3", "殘差變異較大",       1.0, 1.5, "original",
  "S4", "隨機效果正相關情境", 1.0, 1.0, "positive"
)

# ============================================================
# 4. 資料夾、彙整與 QA
# ============================================================

make_scenario_folders <- function(out_dir) {
  folders <- c(
    file.path(out_dir, "simulated_data"),
    file.path(out_dir, "gee_results"),
    file.path(out_dir, "lme_results"),
    file.path(out_dir, "dl_results"),
    file.path(out_dir, "summary_tables"),
    file.path(out_dir, "figures"),
    file.path(out_dir, "logs")
  )
  walk(folders, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))
}

append_failed_run <- function(out_dir, scenario, run, component, reason) {
  failure_path <- file.path(out_dir, "logs", "failed_runs.csv")
  new_row <- tibble(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    scenario = scenario,
    run = as.integer(run),
    component = component,
    reason = as.character(reason)
  )
  
  if (file.exists(failure_path)) {
    old <- read_csv(failure_path, show_col_types = FALSE)
    bind_rows(old, new_row) %>%
      distinct(scenario, run, component, .keep_all = TRUE) %>%
      write_csv(failure_path)
  } else {
    write_csv(new_row, failure_path)
  }
}

validate_stat_result_file <- function(
    path,
    expected_model
) {
  if (!file.exists(path)) {
    return(
      list(
        ok = FALSE,
        reason = "result file does not exist"
      )
    )
  }
  
  x <- tryCatch(
    read_csv(
      path,
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(x)) {
    return(
      list(
        ok = FALSE,
        reason = "result file cannot be read"
      )
    )
  }
  
  required_cols <- c(
    "scenario",
    "run",
    "outcome",
    "model",
    "Group_eval",
    "RMSE",
    "MAE",
    "MAPE",
    "Bias",
    "n_test"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(x)
  )
  
  if (length(missing_cols) > 0L) {
    return(
      list(
        ok = FALSE,
        reason = paste(
          "missing columns:",
          paste(
            missing_cols,
            collapse = ", "
          )
        )
      )
    )
  }
  
  expected_grid <- expand_grid(
    outcome = c(
      "HbA1c",
      "Glu_AC",
      "LDL",
      "TG"
    ),
    model = expected_model,
    Group_eval = c(
      "Overall",
      "Control",
      "Statin"
    )
  )
  
  observed_grid <- x %>%
    distinct(
      outcome,
      model,
      Group_eval
    )
  
  missing_grid <- anti_join(
    expected_grid,
    observed_grid,
    by = c(
      "outcome",
      "model",
      "Group_eval"
    )
  )
  
  metric_cols <- c(
    "RMSE",
    "MAE",
    "MAPE",
    "Bias"
  )
  
  invalid_rows <- x %>%
    filter(
      if_any(
        all_of(metric_cols),
        ~ !is.finite(.x)
      ) |
        !is.finite(n_test) |
        n_test <= 0
    )
  
  expected_row_n <- nrow(expected_grid)
  
  ok <- (
    nrow(x) == expected_row_n &&
      nrow(missing_grid) == 0L &&
      nrow(invalid_rows) == 0L
  )
  
  list(
    ok = ok,
    reason = if (ok) {
      "OK"
    } else {
      paste0(
        "rows=",
        nrow(x),
        "; missing combinations=",
        nrow(missing_grid),
        "; invalid rows=",
        nrow(invalid_rows)
      )
    }
  )
}

validate_dl_result_file <- function(path) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, reason = "DL output file does not exist"))
  }
  
  x <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(x)) {
    return(list(ok = FALSE, reason = "DL output cannot be read"))
  }
  
  required_cols <- c(
    "scenario", "run", "outcome", "model", "Group_eval",
    "RMSE", "MAE", "MAPE", "Bias", "n_test"
  )
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0L) {
    return(list(
      ok = FALSE,
      reason = paste("missing columns:", paste(missing_cols, collapse = ", "))
    ))
  }
  
  complete_grid <- expand_grid(
    outcome = c("HbA1c", "Glu_AC", "LDL", "TG"),
    model = c("LSTM-NARX", "TCN", "ODE-LSTM"),
    Group_eval = c("Overall", "Control", "Statin")
  )
  
  observed_grid <- x %>% distinct(outcome, model, Group_eval)
  missing_grid <- anti_join(
    complete_grid, observed_grid,
    by = c("outcome", "model", "Group_eval")
  )
  
  metric_cols <- c("RMSE", "MAE", "MAPE", "Bias")
  bad_metric <- x %>%
    filter(if_any(all_of(metric_cols), ~ !is.finite(.x)) | n_test <= 0)
  
  if (nrow(x) != 36L || nrow(missing_grid) > 0L || nrow(bad_metric) > 0L) {
    return(list(
      ok = FALSE,
      reason = paste0(
        "incomplete DL result: rows=", nrow(x),
        ", missing combinations=", nrow(missing_grid),
        ", invalid metric rows=", nrow(bad_metric)
      )
    ))
  }
  
  list(ok = TRUE, reason = "OK")
}

rebuild_master_results <- function(out_dir) {
  read_and_write_master <- function(subdir, output_name) {
    files <- list.files(
      file.path(out_dir, subdir),
      pattern = "^S[0-9]_run_.*\\.csv$",
      full.names = TRUE
    )
    if (length(files) == 0L) return(invisible(NULL))
    
    map_dfr(files, read_csv, show_col_types = FALSE) %>%
      arrange(scenario, run, outcome, model, Group_eval) %>%
      write_csv(file.path(out_dir, subdir, output_name))
  }
  
  read_and_write_master("gee_results", "gee_results_master.csv")
  read_and_write_master("lme_results", "lme_results_master.csv")
  read_and_write_master("dl_results", "dl_results_master.csv")
  invisible(TRUE)
}

# ============================================================
# 5. 執行單一情境
# ============================================================

run_one_scenario <- function(
    scenario_name,
    n_runs = N_RUNS,
    run_start = RUN_START,
    run_gee = RUN_GEE,
    run_lme = RUN_LME,
    run_dl = RUN_DL,
    save_sim_data = SAVE_SIM_DATA,
    skip_existing = SKIP_EXISTING
) {
  if (run_start > n_runs) {
    stop("RUN_START 不可大於 N_RUNS。")
  }
  
  sc_setting <- scenario_table %>% filter(scenario == scenario_name)
  if (nrow(sc_setting) != 1L) {
    stop("scenario_table 中找不到情境：", scenario_name)
  }
  
  out_dir <- file.path(SIM_ROOT, scenario_name)
  make_scenario_folders(out_dir)
  
  write_csv(base_param_table,
            file.path(out_dir, "summary_tables", "base_parameters_used.csv"))
  write_csv(base_param_table,
            file.path(out_dir, "summary_tables",
                      "Appendix_Table_E2_S0_data_generation_parameters.csv"))
  
  dl_training_settings <- tibble(
    item = c(
      "模型", "預測目標", "序列長度上限", "最大訓練回合",
      "batch size", "early stopping patience", "評估指標", "評估尺度"
    ),
    setting = c(
      "LSTM-NARX、TCN、ODE-LSTM",
      "測試集每位病患最後一次可評估觀測值",
      as.character(DL_MAX_LEN),
      as.character(DL_EPOCHS),
      as.character(DL_BATCH_SIZE),
      as.character(DL_PATIENCE),
      "RMSE、MAE、MAPE、Bias；分 Overall、Control、Statin",
      "模型於 log/standardized 尺度訓練，轉回原始尺度後計算誤差"
    )
  )
  
  write_csv(dl_training_settings,
            file.path(out_dir, "summary_tables",
                      "Appendix_Table_E1_DL_training_settings.csv"))
  write_csv(scenario_table,
            file.path(out_dir, "summary_tables", "scenario_table.csv"))
  write_csv(skeleton_qc$patient_by_group,
            file.path(out_dir, "summary_tables",
                      "skeleton_patient_by_group_split.csv"))
  write_csv(skeleton_qc$observation_by_group,
            file.path(out_dir, "summary_tables",
                      "skeleton_observation_by_group_split.csv"))
  write_csv(skeleton_qc$repeat_summary,
            file.path(out_dir, "summary_tables",
                      "skeleton_repeated_measure_summary.csv"))
  
  ref_all <- create_reference_trajectories(
    skeleton = skeleton,
    base_params = BASE_PARAMS
  )
  write_csv(ref_all,
            file.path(out_dir, "summary_tables", "reference_trajectory.csv"))
  
  message("\n==============================")
  message("Running scenario: ", scenario_name, " - ", sc_setting$scenario_label)
  message("==============================")
  
  for (run in seq.int(run_start, n_runs)) {
    message("Scenario: ", scenario_name, " | run: ", run)
    run_tag <- sprintf("%04d", run)
    
    sim_path <- file.path(
      out_dir, "simulated_data",
      paste0(scenario_name, "_run_", run_tag, ".csv")
    )
    gee_path <- file.path(
      out_dir, "gee_results",
      paste0(scenario_name, "_run_", run_tag, ".csv")
    )
    lme_path <- file.path(
      out_dir, "lme_results",
      paste0(scenario_name, "_run_", run_tag, ".csv")
    )
    dl_path <- file.path(
      out_dir, "dl_results",
      paste0(scenario_name, "_run_", run_tag, ".csv")
    )
    
    gee_check <- if (
      run_gee &&
      file.exists(gee_path)
    ) {
      validate_stat_result_file(
        gee_path,
        expected_model = "GEE"
      )
    } else {
      list(
        ok = FALSE,
        reason = "not available"
      )
    }
    
    lme_check <- if (
      run_lme &&
      file.exists(lme_path)
    ) {
      validate_stat_result_file(
        lme_path,
        expected_model = "LME"
      )
    } else {
      list(
        ok = FALSE,
        reason = "not available"
      )
    }
    
    dl_check_existing <- if (
      run_dl &&
      file.exists(dl_path)
    ) {
      validate_dl_result_file(dl_path)
    } else {
      list(
        ok = FALSE,
        reason = "not available"
      )
    }
    
    gee_done <- (
      !run_gee ||
        (
          skip_existing &&
            isTRUE(gee_check$ok)
        )
    )
    
    lme_done <- (
      !run_lme ||
        (
          skip_existing &&
            isTRUE(lme_check$ok)
        )
    )
    
    dl_done <- (
      !run_dl ||
        (
          skip_existing &&
            isTRUE(
              dl_check_existing$ok
            )
        )
    )
    
    if (gee_done && lme_done && dl_done) {
      message("  skip existing completed run: ", run)
      next
    }
    
    sim_data <- simulate_full_dataset(
      skeleton = skeleton,
      base_params = BASE_PARAMS,
      scenario_setting = sc_setting,
      seed = SEED_BASE + run
    ) %>%
      mutate(
        scenario = scenario_name,
        run = run,
        ID_num = as.integer(factor(ID))
      ) %>%
      arrange(ID, Time_in_Years, row_id)
    
    if (
      nrow(sim_data) != SKELETON_ROW_N ||
      n_distinct(sim_data$ID) != SKELETON_PATIENT_N ||
      !identical(
        as.integer(sim_data$row_id),
        as.integer(skeleton$row_id)
      )
    ) {
      stop(
        "模擬後資料列數、病患數或 row_id 與固定 skeleton 不一致。",
        call. = FALSE
      )
    }
    
    write_csv(sim_data, sim_path)
    
    if (run_gee && !gee_done) {
      tryCatch({
        gee_res <- fit_all_gee_pred(
          sim_data = sim_data,
          base_params = BASE_PARAMS,
          scenario = scenario_name,
          run = run
        )
        write_csv(gee_res, gee_path)
      }, error = function(e) {
        append_failed_run(out_dir, scenario_name, run, "GEE", conditionMessage(e))
        warning("GEE failed: ", scenario_name, " run ", run)
      })
    }
    
    if (run_lme && !lme_done) {
      tryCatch({
        lme_res <- fit_all_lme_pred(
          sim_data = sim_data,
          base_params = BASE_PARAMS,
          scenario = scenario_name,
          run = run,
          random_slope = TRUE
        )
        write_csv(lme_res, lme_path)
      }, error = function(e) {
        append_failed_run(out_dir, scenario_name, run, "LME", conditionMessage(e))
        warning("LME failed: ", scenario_name, " run ", run)
      })
    }
    
    if (run_dl && !dl_done) {
      args <- c(
        PYTHON_PREFIX_ARGS,
        PYTHON_SCRIPT,
        "--input", sim_path,
        "--output", dl_path,
        "--seed", as.character(SEED_BASE + run),
        "--scenario", scenario_name,
        "--run", as.character(run),
        "--max_len", as.character(DL_MAX_LEN),
        "--epochs", as.character(DL_EPOCHS),
        "--batch_size", as.character(DL_BATCH_SIZE),
        "--patience", as.character(DL_PATIENCE)
      )
      
      py_log <- file.path(
        out_dir, "logs",
        paste0("python_", scenario_name, "_run_", run_tag, ".log")
      )
      
      status <- system2(
        PYTHON_EXE,
        args = args,
        stdout = py_log,
        stderr = py_log
      )
      
      dl_check <- validate_dl_result_file(dl_path)
      
      if (status != 0L || !isTRUE(dl_check$ok)) {
        append_failed_run(
          out_dir, scenario_name, run, "DL",
          paste0("status=", status, "; ", dl_check$reason)
        )
        if (file.exists(dl_path)) {
          invalid_path <- sub("\\.csv$", "_invalid.csv", dl_path)
          file.rename(dl_path, invalid_path)
        }
        warning("Python DL failed or incomplete: ", scenario_name, " run ", run)
      }
    }
    
    keep_representative <- (
      KEEP_S0_RUN1_FOR_DIAGNOSTICS &&
        scenario_name == "S0" &&
        run == 1L
    )
    
    if (!save_sim_data && !keep_representative && file.exists(sim_path)) {
      file.remove(sim_path)
    }
    
    if (run %% 10L == 0L) {
      rebuild_master_results(out_dir)
      invisible(gc())
    }
  }
  
  rebuild_master_results(out_dir)
  
  summary_res <- summarise_simulation_outputs(
    out_dir = out_dir,
    scenario_name = scenario_name
  )
  
  plot_reference_trajectory(
    out_dir = out_dir,
    scenario_name = scenario_name
  )
  
  if (scenario_name == "S0") {
    representative_path <- file.path(
      out_dir, "simulated_data", "S0_run_0001.csv"
    )
    
    if (file.exists(representative_path)) {
      plot_S0_residual_diagnostics(
        representative_data_path = representative_path,
        out_dir = out_dir,
        scenario_name = "S0"
      )
    } else {
      warning("找不到 S0_run_0001.csv，略過附錄圖 E1/E2。")
    }
  }
  
  invisible(summary_res)
}

# ============================================================
# 6. 開始執行
# ============================================================
ask_continue_scenario <- function(current_scenario, next_scenario) {
  question_text <- paste0("情境 ", current_scenario, 
                          " 已完成。\n\n是否繼續執行情境 ",
                          next_scenario, "？")
  
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
    isTRUE(rstudioapi::isAvailable())) {
    answer <- rstudioapi::showQuestion(
      title = "模擬情境執行確認",
      message = question_text,
      ok = "繼續",
      cancel = "停止"
    )
    
    return(isTRUE(answer))
  }
  
  if (!interactive()) {
    message(
      "目前不是互動式 R 工作階段，無法顯示詢問介面；",
      "為避免未經確認自動執行下一情境，程式將暫停。"
    )
    return(FALSE)
  }
  
  repeat {
    answer <- readline(prompt = paste0("\n", question_text, "\n請輸入 y 繼續或 n 停止："))
    
    answer <- tolower(trimws(answer))
    
    if (answer %in% c("y", "yes")) {
      return(TRUE)
    }
    
    if (answer %in% c("n", "no")) {
      return(FALSE)
    }
    
    message("未偵測到有效輸入。空白輸入不會停止，", "請明確輸入 y 或 n。")
    
    Sys.sleep(0.5)
  }
}

start_time <- Sys.time()
user_stopped <- FALSE
completed_scenarios <- character(0)

for (i in seq_along(SCENARIOS_TO_RUN)) {
  
  scenario_name <- SCENARIOS_TO_RUN[[i]]
  
  run_one_scenario(
    scenario_name = scenario_name,
    n_runs = N_RUNS,
    run_start = RUN_START,
    run_gee = RUN_GEE,
    run_lme = RUN_LME,
    run_dl = RUN_DL,
    save_sim_data = SAVE_SIM_DATA,
    skip_existing = SKIP_EXISTING
  )
  
  completed_scenarios <- c(completed_scenarios, scenario_name)
  
  # 每完成一個情境後，立即彙整目前已完成的情境。
  consolidate_all_scenarios(
    sim_root = SIM_ROOT,
    scenarios = completed_scenarios
  )
  
  # 最後一個情境完成後不再詢問。
  if (
    ASK_BEFORE_NEXT_SCENARIO &&
    i < length(SCENARIOS_TO_RUN)
  ) {
    next_scenario <- SCENARIOS_TO_RUN[[i + 1L]]
    
    continue_next <-
      ask_continue_scenario(
        current_scenario = scenario_name,
        next_scenario = next_scenario
      )
    
    if (!continue_next) {
      user_stopped <- TRUE
      
      message(
        "使用者選擇暫停。已完成情境：",
        paste(completed_scenarios, collapse = ", ")
      )
      
      break
    }
  }
}

end_time <- Sys.time()

if (user_stopped) {
  message(
    "模擬已暫停，尚未完成全部指定情境。"
  )
} else {
  message(
    "全部指定情境執行完成。"
  )
}

message("本次已處理情境：", paste(completed_scenarios, collapse = ", "))

message("總執行時間：")
print(end_time - start_time)

message("主要輸出資料夾：", SIM_ROOT)
