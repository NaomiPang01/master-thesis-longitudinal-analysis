# ============================================================
# 研究流程：
# 1. 讀取最新版「糖尿病患回診資料_with_type.csv／.xlsx」；若兩者皆通過最新版 schema，
#    必須再通過分析關鍵欄位一致性檢查，避免誤讀舊版或不同步檔案。
# 2. 主分析限定第二型糖尿病（Type 2）且基準年齡 >= 18 歲。
# 3. 依符合主分析追蹤事件之就診層級用藥紀錄建立固定控制組／史他汀組。
# 4. 建立 T1–T19 三個月時間區間，並保留至少兩個不同追蹤日期者。
# 5. 表 4.1 使用每位病患第一筆合格觀測之原始 baseline 資料。
# 6. 同一病患於同一報告日期之多筆紀錄，在進入 MICE、GEE/LME 與 DL 前先彙整為單一時間點；
#    連續型檢驗值與 BMI 取非缺失平均，固定變項需保持一致，並保留來源就診追溯欄位。
# 7. Outcome 與 BMI 僅建立 QA flags，不以臨床範圍、IQR 或相鄰突變自動刪值。
# 8. 推論資料流程：完整 cohort 以 MICE 產生 10 組完成資料，僅供 GEE/LME
#    參數推論、模型配適與診斷；Age、Sex、血壓、HDL、TC 不插補。
# 9. 預測資料流程：先固定病患層級 train/validation/test，再建立單一
#    prediction-safe dataset。四項 outcome 不插補，歷史序列與測試目標僅使用
#    原始觀測值；BMI 缺失僅使用 train set 估計之分組／整體中位數補值。
# 10. GEE/LME 與 DL 實證預測共用同一 prediction-safe dataset；
#    validation 僅供 DL 早期停止與模型選擇，test 每位病患取最後一次原始觀測 target。
# 11. 所有主要圖同時輸出全中文與全英文版本。
# 12. 先完整統計原始時間類別（目前資料可見 T1–T21），
#     再依預先設定排除 T20、T21，主分析保留 T1–T19。
# 13. 分別輸出 inference 用 10 組 MICE 資料與 prediction-safe 單一資料，
#     不得將完整 cohort MICE 資料用於實證預測比較。
# 14. 全體病患軌跡圖預設不輸出，避免嚴重 overplotting；
#     正文主圖使用分組平均值、95% CI，並僅於指定時間點標示 n。
#
# 圖形命名：
# - 中文版：*_zh.png；組別為「控制組」「史他汀組」
# - 英文版：*_en.png；組別為「Control group」「Statin group」
#
# 主要 outcome：HbA1c、Glu-AC、LDL、TG
# HDL、TC：僅作補充描述
# Systolic_BP、Diastolic_BP：若上游檔案有提供，僅作附錄資料可用性描述
# ============================================================

# ============================================================
# 0. 套件、路徑與固定設定
# ============================================================

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "lubridate",
  "purrr", "ggplot2", "patchwork", "mice",
  "writexl", "readr", "tibble"
)

invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "尚未安裝套件：", pkg,
      "。請先執行 install.packages('", pkg, "')"
    )
  }
  library(pkg, character.only = TRUE)
}))

set.seed(5014)

LOCAL_TZ <- "Asia/Taipei"
RANDOM_SEED <- 5014L
PLOT_DPI <- 600L
MICE_M <- 10L
MICE_MAXIT <- 10L
SCRIPT_VERSION <- "scheme-A-final-v10"

raw_root <- "raw_data"
preprocess_folder <- "糖尿病患抽血資料整理"
output_root <- "analysis_outputs"
table_folder <- file.path(output_root, "tables")
figure_folder <- file.path(output_root, "figures")
rds_folder <- file.path(output_root, "rds")
dl_folder <- file.path(output_root, "dl_data")
qc_folder <- file.path(output_root, "qc")

dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
dir.create(table_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(dl_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_folder, showWarnings = FALSE, recursive = TRUE)

# 清除本版已停止輸出的舊檔，避免與本次結果混淆。
obsolete_output_paths <- c(
  file.path(rds_folder, "imp_type2_inference.rds"),
  file.path(rds_folder, "imp_long_with_split_type2_inference.rds"),
  file.path(rds_folder, "prediction_long_with_split_type2.rds"),
  file.path(dl_folder, "dl_raw_with_split_type2.csv"),
  file.path(dl_folder, "dl_raw_with_split_type2_full_qc.csv"),
  file.path(dl_folder, "dl_raw_train_type2.csv"),
  file.path(dl_folder, "dl_raw_validation_type2.csv"),
  file.path(dl_folder, "dl_raw_test_type2.csv"),
  file.path(dl_folder, "imputed_long_train_type2.csv"),
  file.path(dl_folder, "imputed_long_validation_type2.csv"),
  file.path(dl_folder, "imputed_long_test_type2.csv")
)

obsolete_output_paths <- c(
  obsolete_output_paths,
  list.files(
    dl_folder,
    pattern = "^imputed_dataset_[0-9]+_type2[.]csv$",
    full.names = TRUE
  ),
  list.files(
    figure_folder,
    pattern = "^(mean_trend_ci_|trajectory_sampled_|trajectory_all_).*[.]png$",
    full.names = TRUE
  )
)

obsolete_output_paths <- unique(
  obsolete_output_paths[file.exists(obsolete_output_paths)]
)

if (length(obsolete_output_paths) > 0L) {
  unlink(obsolete_output_paths)
}

# CSV 與 XLSX 應由同一次最新版資料串檔產生。
# 若兩者皆通過最新版 schema，本程式會比較分析關鍵欄位；不一致時直接停止，
# 避免因檔名相同但版本不同而靜默誤用舊資料。XLSX 作為第一順位來源。
main_xlsx_path <- file.path(
  preprocess_folder,
  "糖尿病患回診資料_with_type.xlsx"
)
main_csv_path <- file.path(
  preprocess_folder,
  "糖尿病患回診資料_with_type.csv"
)

main_data_candidates <- unique(
  c(
    main_xlsx_path,
    main_csv_path,
    "糖尿病患回診資料_with_type.xlsx",
    "糖尿病患回診資料_with_type.csv"
  )
)

main_data_candidates <- main_data_candidates[
  file.exists(main_data_candidates)
]

if (length(main_data_candidates) == 0L) {
  stop(
    "找不到糖尿病患回診資料_with_type.csv 或 .xlsx。",
    "請確認檔案路徑。"
  )
}

med_path <- file.path(raw_root, "病人藥品明細.xlsx")
diag_path <- file.path(raw_root, "病人疾病診斷明細.xlsx")

if (!file.exists(med_path)) {
  stop("找不到藥品明細檔：", med_path)
}
if (!file.exists(diag_path)) {
  stop("找不到疾病診斷明細檔：", diag_path)
}

# 全體病患軌跡圖在大樣本下會嚴重重疊，僅作資料品質檢查，不建議放正文。
# 正文以分組平均趨勢與 95% CI 為主；附錄另輸出各組隨機抽樣軌跡圖。
make_all_trajectory <- FALSE
make_sampled_trajectory <- TRUE
sample_per_group_trajectory <- 100L
N_LABEL_TIMEPOINTS <- c(1L, 4L, 7L, 10L, 13L, 16L, 19L)

# T1–T19：每 3 個月一區間，共 19 個區間。
MAX_TIME_CATEGORY_INDEX <- 18L

# 僅作 QA 標記，不自動刪值。
clinical_review_ranges <- list(
  HbA1c = c(3, 20),
  Glu_AC = c(40, 600),
  HDL = c(10, 120),
  LDL = c(20, 400),
  TC = c(70, 500),
  TG = c(20, 1000),
  BMI = c(10, 80)
)

primary_outcomes <- c("HbA1c", "Glu_AC", "LDL", "TG")
secondary_outcomes <- c("HDL", "TC")
response_vars <- c(primary_outcomes, secondary_outcomes)

outcome_labels <- c(
  HbA1c = "HbA1c",
  Glu_AC = "Glu-AC",
  LDL = "LDL",
  TG = "TG",
  HDL = "HDL",
  TC = "TC"
)

# ============================================================
# 1. 工具函數
# ============================================================

assert_has_cols <- function(dat, required_cols, data_name = "資料") {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      data_name, "缺少必要欄位：",
      paste(missing_cols, collapse = ", ")
    )
  }
  invisible(TRUE)
}

read_xlsx_hospital <- function(path, drop_first_col = TRUE) {
  if (!file.exists(path)) stop("找不到檔案：", path)

  dat <- readxl::read_excel(path)

  if (drop_first_col && ncol(dat) >= 1L) {
    first_col_name <- names(dat)[1L]
    if (
      is.na(first_col_name) ||
      first_col_name == "" ||
      grepl("^\\.\\.\\.\\d+$", first_col_name)
    ) {
      dat <- dat[, -1, drop = FALSE]
    }
  }

  dat %>%
    dplyr::select(-dplyr::matches("^\\.\\.\\.\\d+$"))
}

read_main_data <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    return(readxl::read_excel(path))
  }

  if (ext != "csv") {
    stop("不支援的主資料格式：", ext)
  }

  # 先讀取欄名，再只對實際存在的日期時間欄指定 character parser。
  # 如此可避免上游未提供 BP測量日期時間 時出現 named parser 警告。
  header_names <- names(
    readr::read_csv(
      path,
      n_max = 0,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
  )

  spec_args <- list(
    .default = readr::col_guess()
  )

  datetime_candidates <- intersect(
    c("測量日期時間", "BP測量日期時間"),
    header_names
  )

  for (v in datetime_candidates) {
    spec_args[[v]] <- readr::col_character()
  }

  csv_col_spec <- do.call(
    readr::cols,
    spec_args
  )

  readr::read_csv(
    path,
    show_col_types = FALSE,
    na = c("", "NA", "NaN", "NULL"),
    col_types = csv_col_spec
  )
}

format_datetime_for_export <- function(x) {
  if (inherits(x, c("POSIXct", "POSIXt"))) {
    out <- rep(NA_character_, length(x))
    ok <- !is.na(x)
    out[ok] <- format(
      x[ok],
      format = "%Y-%m-%d %H:%M:%S",
      tz = LOCAL_TZ
    )
    return(out)
  }

  x_chr <- as.character(x)
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x_chr
}

prepare_export_data <- function(dat) {
  out <- dat

  datetime_cols <- intersect(
    c("測量日期時間", "BP測量日期時間"),
    names(out)
  )

  for (v in datetime_cols) {
    out[[v]] <- format_datetime_for_export(out[[v]])
  }

  out
}

parse_date_from_excel <- function(x) {
  if (inherits(x, "Date")) return(x)

  if (inherits(x, c("POSIXct", "POSIXt"))) {
    return(as.Date(x, tz = LOCAL_TZ))
  }

  if (is.numeric(x)) {
    return(as.Date(x, origin = "1899-12-30"))
  }

  x_chr <- stringr::str_squish(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_

  parsed <- suppressWarnings(
    lubridate::parse_date_time(
      x_chr,
      orders = c(
        "Ymd HMS", "Ymd HM", "Ymd",
        "ymd HMS", "ymd HM", "ymd",
        "mdy HMS", "mdy HM", "mdy"
      ),
      tz = LOCAL_TZ
    )
  )

  as.Date(parsed, tz = LOCAL_TZ)
}

standardize_id <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

standardize_sex_binary <- function(x) {
  x_chr <- stringr::str_squish(as.character(x))

  dplyr::case_when(
    x_chr %in% c(
      "1", "男", "M", "m",
      "Male", "male", "MALE"
    ) ~ 1,
    x_chr %in% c(
      "0", "女", "F", "f",
      "Female", "female", "FEMALE"
    ) ~ 0,
    TRUE ~ NA_real_
  )
}

as_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

mean_or_na <- function(x) {
  x_num <- as_numeric_safely(x)
  x_valid <- x_num[!is.na(x_num) & is.finite(x_num)]

  if (length(x_valid) == 0L) {
    NA_real_
  } else {
    mean(x_valid)
  }
}

first_nonmissing <- function(x) {
  idx <- which(!is.na(x))

  if (length(idx) == 0L) {
    return(x[NA_integer_])
  }

  x[idx[1L]]
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) NA_real_ else stats::sd(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else stats::median(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    NA_real_
  } else {
    as.numeric(stats::quantile(x, prob, names = FALSE))
  }
}

critical_t_975 <- function(n) {
  if (is.na(n) || n <= 1L) NA_real_ else stats::qt(0.975, df = n - 1L)
}

format_p_value <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.0001 ~ "<0.0001",
    TRUE ~ formatC(p, format = "f", digits = 4)
  )
}

make_missing_summary <- function(dat, vars, group_vars = NULL) {
  vars <- intersect(vars, names(dat))

  if (length(vars) == 0L) {
    return(tibble::tibble())
  }

  if (is.null(group_vars) || length(group_vars) == 0L) {
    return(
      purrr::map_dfr(vars, function(v) {
        tibble::tibble(
          variable = v,
          total_n = nrow(dat),
          missing_n = sum(is.na(dat[[v]])),
          observed_n = sum(!is.na(dat[[v]])),
          missing_rate_percent = round(
            mean(is.na(dat[[v]])) * 100,
            2
          )
        )
      })
    )
  }

  dat %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_vars))
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      summary = purrr::map(
        data,
        ~ make_missing_summary(.x, vars)
      )
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(summary) %>%
    dplyr::ungroup()
}

detect_adjacent_jump <- function(x, threshold = 3) {
  out <- rep(FALSE, length(x))
  valid_idx <- which(!is.na(x) & is.finite(x))

  if (length(valid_idx) < 3L) return(out)

  sd_x <- stats::sd(x[valid_idx])
  if (is.na(sd_x) || sd_x == 0) return(out)

  for (k in seq.int(2L, length(valid_idx))) {
    current_idx <- valid_idx[k]
    previous_idx <- valid_idx[k - 1L]

    if (
      abs(x[current_idx] - x[previous_idx]) >
        threshold * sd_x
    ) {
      out[current_idx] <- TRUE
    }
  }

  out
}

add_quality_flags <- function(
  dat,
  id_var = "ID",
  date_var = "報告日期",
  visit_var = "就診號"
) {
  out <- dat
  vars_to_check <- intersect(
    names(clinical_review_ranges),
    names(out)
  )

  for (v in vars_to_check) {
    range_v <- clinical_review_ranges[[v]]

    clinical_flag <- paste0(v, "_clinical_review_flag")
    adjacent_flag <- paste0(v, "_adjacent_jump_flag")
    iqr_flag <- paste0(v, "_iqr_outlier_flag")

    out[[clinical_flag]] <- (
      !is.na(out[[v]]) &
        (
          out[[v]] < range_v[1] |
            out[[v]] > range_v[2]
        )
    )

    out <- out %>%
      dplyr::group_by(.data[[id_var]]) %>%
      dplyr::arrange(
        .data[[date_var]],
        .data[[visit_var]],
        .by_group = TRUE
      ) %>%
      dplyr::mutate(
        "{adjacent_flag}" := detect_adjacent_jump(.data[[v]])
      ) %>%
      dplyr::ungroup()

    x_valid <- out[[v]][
      !is.na(out[[v]]) &
        is.finite(out[[v]])
    ]

    if (length(x_valid) > 0L) {
      q1 <- safe_quantile(x_valid, 0.25)
      q3 <- safe_quantile(x_valid, 0.75)
      iqr_value <- q3 - q1
      lower <- q1 - 1.5 * iqr_value
      upper <- q3 + 1.5 * iqr_value

      out[[iqr_flag]] <- (
        !is.na(out[[v]]) &
          (
            out[[v]] < lower |
              out[[v]] > upper
          )
      )
    } else {
      out[[iqr_flag]] <- FALSE
    }
  }

  out
}

summarise_quality_flags <- function(dat) {
  flag_vars <- names(dat)[
    stringr::str_detect(
      names(dat),
      "_clinical_review_flag$|_adjacent_jump_flag$|_iqr_outlier_flag$"
    )
  ]

  if (length(flag_vars) == 0L) {
    return(
      tibble::tibble(
        flag_variable = character(),
        flagged_n = integer()
      )
    )
  }

  dat %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(flag_vars),
        ~ sum(.x %in% TRUE, na.rm = TRUE)
      )
    ) %>%
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "flag_variable",
      values_to = "flagged_n"
    ) %>%
    dplyr::arrange(flag_variable)
}

categorical_test_info <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)

  if (
    length(unique(x[ok])) < 2L ||
      length(unique(g[ok])) < 2L
  ) {
    return(
      tibble::tibble(
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        test = NA_character_
      )
    )
  }

  tab <- table(x[ok], g[ok])

  chi_obj <- tryCatch(
    suppressWarnings(stats::chisq.test(tab, correct = FALSE)),
    error = function(e) NULL
  )

  use_fisher <- (
    is.null(chi_obj) ||
      any(chi_obj$expected < 5)
  )

  if (use_fisher) {
    p <- tryCatch(
      stats::fisher.test(tab)$p.value,
      error = function(e) NA_real_
    )
    test_name <- "Fisher's exact test"
  } else {
    p <- chi_obj$p.value
    test_name <- "Pearson chi-square test"
  }

  tibble::tibble(
    p_value = p,
    p_value_formatted = format_p_value(p),
    test = test_name
  )
}

continuous_test_info <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)

  if (length(unique(g[ok])) < 2L) {
    return(
      tibble::tibble(
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        test = "Wilcoxon rank-sum test"
      )
    )
  }

  p <- tryCatch(
    stats::wilcox.test(
      x[ok] ~ g[ok],
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )

  tibble::tibble(
    p_value = p,
    p_value_formatted = format_p_value(p),
    test = "Wilcoxon rank-sum test"
  )
}

split_ids_by_group <- function(
  df_group,
  train_prop = 0.70,
  validation_prop = 0.15
) {
  n <- nrow(df_group)
  if (n == 0L) return(df_group)

  idx <- sample(seq_len(n))
  n_train <- floor(train_prop * n)
  n_validation <- floor(validation_prop * n)

  split_value <- rep("test", n)

  if (n_train > 0L) {
    split_value[idx[seq_len(n_train)]] <- "train"
  }

  if (n_validation > 0L) {
    validation_index <- seq.int(
      n_train + 1L,
      n_train + n_validation
    )
    split_value[idx[validation_index]] <- "validation"
  }

  df_group$split <- split_value
  df_group
}

# ============================================================
# 2. 讀取主資料、糖尿病類型與成人資格
# ============================================================

main_required_cols <- c(
  "病歷號",
  "就診號",
  "報告日期",
  "Repeat",
  "diabetes_type",
  "性別",
  "Age",
  "身高",
  "體重",
  "BMI",
  "basic_time_gap_days",
  "HbA1C",
  "Glucose_AC",
  "LDL-C",
  "Triglycerides",
  "HDL-C",
  "Total_cholesterol",
  "HbA1C_left_censored",
  "HbA1C_right_censored"
)

# 最新版資料串檔主檔不再強制輸出下列欄位。
# Age_source／sex_source 僅供追溯；BP 若存在才做附錄缺失率，
# 不得因舊版表格.R 需求而迫使上游主檔重新加入。
main_optional_cols <- c(
  "Age_source",
  "sex_source",
  "BP_time_gap_days",
  "Systolic_BP",
  "Diastolic_BP",
  "BP測量日期",
  "BP測量日期時間"
)

normalize_main_for_compare <- function(dat) {
  numeric_cols <- intersect(
    c(
      "Repeat", "Age", "身高", "體重", "BMI",
      "basic_time_gap_days", "HbA1C", "Glucose_AC",
      "LDL-C", "Triglycerides", "HDL-C",
      "Total_cholesterol", "HbA1C_left_censored",
      "HbA1C_right_censored"
    ),
    names(dat)
  )

  out <- dat %>%
    dplyr::transmute(
      病歷號 = standardize_id(病歷號),
      就診號 = standardize_id(就診號),
      報告日期 = as.character(
        parse_date_from_excel(報告日期)
      ),
      diabetes_type = as.character(diabetes_type),
      性別 = as.character(性別),
      dplyr::across(
        dplyr::all_of(numeric_cols),
        ~ round(as_numeric_safely(.x), 8)
      )
    ) %>%
    dplyr::arrange(
      病歷號,
      報告日期,
      就診號,
      Repeat
    )

  out <- as.data.frame(
    out,
    stringsAsFactors = FALSE
  )
  attr(out, "spec") <- NULL
  attr(out, "problems") <- NULL
  rownames(out) <- NULL
  out
}

main_read_attempts <- list()
main_valid_data <- list()

for (candidate_path in main_data_candidates) {
  candidate_data <- tryCatch(
    read_main_data(candidate_path),
    error = function(e) {
      main_read_attempts[[candidate_path]] <<- tibble::tibble(
        path = candidate_path,
        file_type = tolower(tools::file_ext(candidate_path)),
        status = "read_error",
        row_n = NA_integer_,
        column_n = NA_integer_,
        missing_required_columns = conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(candidate_data)) {
    next
  }

  missing_candidate_cols <- setdiff(
    main_required_cols,
    names(candidate_data)
  )

  main_read_attempts[[candidate_path]] <- tibble::tibble(
    path = candidate_path,
    file_type = tolower(tools::file_ext(candidate_path)),
    status = if (
      length(missing_candidate_cols) == 0L
    ) {
      "accepted"
    } else {
      "rejected_schema_mismatch"
    },
    row_n = nrow(candidate_data),
    column_n = ncol(candidate_data),
    missing_required_columns = if (
      length(missing_candidate_cols) == 0L
    ) {
      ""
    } else {
      paste(missing_candidate_cols, collapse = " | ")
    }
  )

  if (length(missing_candidate_cols) == 0L) {
    main_valid_data[[candidate_path]] <- candidate_data
  } else {
    warning(
      "主資料檔欄位不符合最新版 schema，略過：",
      candidate_path,
      "；缺少欄位：",
      paste(missing_candidate_cols, collapse = ", ")
    )
  }
}

main_source_summary <- dplyr::bind_rows(
  main_read_attempts
)

accepted_main_paths <- names(main_valid_data)

if (length(accepted_main_paths) == 0L) {
  stop(
    "找到的 CSV／XLSX 均無法通過最新版欄位檢查。",
    "請重新執行最新版資料串檔.R。"
  )
}

reference_main_path <- accepted_main_paths[[1L]]
reference_main_normalized <- normalize_main_for_compare(
  main_valid_data[[reference_main_path]]
)

main_source_consistency_qc <- purrr::map_dfr(
  accepted_main_paths,
  function(path_i) {
    normalized_i <- normalize_main_for_compare(
      main_valid_data[[path_i]]
    )

    tibble::tibble(
      reference_path = reference_main_path,
      compared_path = path_i,
      row_n = nrow(normalized_i),
      column_n = ncol(main_valid_data[[path_i]]),
      analysis_critical_columns_identical = identical(
        reference_main_normalized,
        normalized_i
      )
    )
  }
)

if (
  any(
    !main_source_consistency_qc$
      analysis_critical_columns_identical
  )
) {
  inconsistent_paths <- main_source_consistency_qc %>%
    dplyr::filter(
      !analysis_critical_columns_identical
    ) %>%
    dplyr::pull(compared_path)

  stop(
    "同名 CSV 與 XLSX 的分析關鍵欄位不一致，為避免誤用舊版資料已停止執行。\n",
    "請先重新執行最新版資料串檔.R，使 CSV 與 XLSX 同步，",
    "或移除明確不再使用的舊版檔案。\n",
    "不一致檔案：",
    paste(inconsistent_paths, collapse = " | ")
  )
}

main_data_path <- reference_main_path
dt_raw <- main_valid_data[[main_data_path]]

assert_has_cols(
  dt_raw,
  main_required_cols,
  "最新版串檔主資料"
)

message("主資料讀取來源：", main_data_path)

main_optional_column_qc <- tibble::tibble(
  variable = main_optional_cols,
  present_in_input = main_optional_cols %in%
    names(dt_raw),
  role = c(
    "optional provenance",
    "optional provenance",
    "optional BP QA",
    "optional BP appendix",
    "optional BP appendix",
    "optional BP provenance",
    "optional BP provenance"
  )
)

if (!("Age_source" %in% names(dt_raw))) {
  dt_raw$Age_source <-
    "Not provided by preprocessing file"
}

if (!("sex_source" %in% names(dt_raw))) {
  dt_raw$sex_source <-
    "Not provided by preprocessing file"
}

dt_raw <- dt_raw %>%
  dplyr::mutate(
    病歷號 = standardize_id(病歷號),
    就診號 = standardize_id(就診號),
    報告日期 = parse_date_from_excel(報告日期),
    diabetes_type = as.character(diabetes_type),
    Age_source = as.character(Age_source),
    sex_source = as.character(sex_source)
  )

type_count_raw_patient <- dt_raw %>%
  dplyr::distinct(病歷號, diabetes_type) %>%
  dplyr::count(diabetes_type, name = "patient_n") %>%
  dplyr::arrange(diabetes_type)

type_count_raw_observation <- dt_raw %>%
  dplyr::count(diabetes_type, name = "observation_n") %>%
  dplyr::arrange(diabetes_type)

dt_type2 <- dt_raw %>%
  dplyr::filter(diabetes_type == "Type 2") %>%
  dplyr::arrange(病歷號, 報告日期, 就診號)

type2_baseline_age <- dt_type2 %>%
  dplyr::group_by(病歷號) %>%
  dplyr::summarise(
    Age_baseline_eligibility = dplyr::first(Age),
    Age_source_baseline = dplyr::first(Age_source),
    observation_n = dplyr::n(),
    distinct_report_date_n = dplyr::n_distinct(報告日期),
    first_report_date = min(報告日期),
    last_report_date = max(報告日期),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    adult_status = dplyr::case_when(
      is.na(Age_baseline_eligibility) ~ "Baseline age missing",
      Age_baseline_eligibility < 18 ~ "Age < 18",
      TRUE ~ "Age >= 18"
    )
  )

adult_eligibility_summary <- type2_baseline_age %>%
  dplyr::group_by(adult_status) %>%
  dplyr::summarise(
    patient_n = dplyr::n(),
    observation_n = sum(observation_n),
    .groups = "drop"
  ) %>%
  dplyr::arrange(adult_status)

under18_patient_detail <- type2_baseline_age %>%
  dplyr::filter(adult_status == "Age < 18") %>%
  dplyr::arrange(Age_baseline_eligibility, 病歷號)

age_missing_patient_detail <- type2_baseline_age %>%
  dplyr::filter(adult_status == "Baseline age missing") %>%
  dplyr::arrange(病歷號)

adult_type2_ids <- type2_baseline_age %>%
  dplyr::filter(
    !is.na(Age_baseline_eligibility),
    Age_baseline_eligibility >= 18
  ) %>%
  dplyr::select(病歷號)

dt_adult_type2 <- dt_type2 %>%
  dplyr::semi_join(
    adult_type2_ids,
    by = "病歷號"
  )

message("主分析限定成人 Type 2 diabetes。")
message(
  "Type 2 病患數：",
  dplyr::n_distinct(dt_type2$病歷號)
)
# Type 2 病患數：3483
message(
  "成人 Type 2 病患數：",
  dplyr::n_distinct(dt_adult_type2$病歷號)
)
# 成人 Type 2 病患數：3474
message(
  "未滿 18 歲病患數：",
  nrow(under18_patient_detail)
)
# 未滿 18 歲病患數：2
message(
  "基準年齡缺失病患數：",
  nrow(age_missing_patient_detail)
)
# 基準年齡缺失病患數：7

# ============================================================
# 3. 定義固定史他汀組與控制組
# ============================================================

dt_med <- read_xlsx_hospital(med_path)

assert_has_cols(
  dt_med,
  c("病歷號", "就診號", "藥品名稱"),
  "病人藥品明細"
)

dt_med <- dt_med %>%
  dplyr::mutate(
    病歷號 = standardize_id(病歷號),
    就診號 = standardize_id(就診號),
    藥品名稱 = stringr::str_squish(as.character(藥品名稱)),
    drug_name_lower = stringr::str_to_lower(藥品名稱)
  )

dt_adult_key <- dt_adult_type2 %>%
  dplyr::select(病歷號, 就診號) %>%
  dplyr::distinct()

dt_med <- dt_med %>%
  dplyr::semi_join(
    dt_adult_key,
    by = c("病歷號", "就診號")
  )

# 史他汀成分與院內實際出現之商品名。
# Rosis 為 furosemide，不列入史他汀。
# Caduet、Pravafen、Atorcal、Agitin、Dehypotin 均含史他汀成分。
statin_terms <- c(
  "atorvastatin",
  "rosuvastatin",
  "simvastatin",
  "pravastatin",
  "lovastatin",
  "fluvastatin",
  "pitavastatin",
  "crestor",
  "lipitor",
  "zulitor",
  "rosulator",
  "delipic",
  "mevalotin",
  "livalo",
  "tulip",
  "vytorin",
  "atozet",
  "rosuzet",
  "cretrol",
  "caduet",
  "pravafen",
  "atorcal",
  "agitin",
  "dehypotin",
  "lescol",
  "zocor",
  "pravachol"
)

statin_pattern <- paste(statin_terms, collapse = "|")
ezetimibe_single_pattern <- "ezetimibe|ezetor|ezetrol"
rosis_pattern <- "(^|[^a-z])rosis([^a-z]|$)"

dt_med_classified <- dt_med %>%
  dplyr::mutate(
    has_statin = dplyr::coalesce(
      stringr::str_detect(
        drug_name_lower,
        statin_pattern
      ),
      FALSE
    ),
    has_ezetimibe_keyword = dplyr::coalesce(
      stringr::str_detect(
        drug_name_lower,
        ezetimibe_single_pattern
      ),
      FALSE
    ),
    is_rosis = dplyr::coalesce(
      stringr::str_detect(
        drug_name_lower,
        rosis_pattern
      ),
      FALSE
    )
  )

statin_drug_dictionary <- tibble::tibble(
  statin_term = statin_terms
)

statin_drug_check <- dt_med_classified %>%
  dplyr::distinct(
    藥品名稱,
    has_statin,
    has_ezetimibe_keyword,
    is_rosis
  ) %>%
  dplyr::arrange(
    dplyr::desc(has_statin),
    藥品名稱
  )

critical_statin_brand_terms <- c(
  "caduet",
  "pravafen",
  "atorcal",
  "agitin",
  "dehypotin"
)

critical_statin_brand_pattern <- paste(
  critical_statin_brand_terms,
  collapse = "|"
)

known_statin_brand_qc <- dt_med_classified %>%
  dplyr::filter(
    stringr::str_detect(
      drug_name_lower,
      critical_statin_brand_pattern
    )
  ) %>%
  dplyr::distinct(
    藥品名稱,
    has_statin
  ) %>%
  dplyr::arrange(藥品名稱)

# 利用院內命名標記與常見成分字根，列出仍未被分類的疑似降血脂藥品。
# 此表只作人工 QA，不直接變更分組。
lipid_name_review_pattern <- paste(
  c(
    "血脂",
    "ator",
    "rosu",
    "simva",
    "prava",
    "pitava",
    "lova",
    "fluva",
    "ezeti",
    "crestor",
    "lipitor",
    "livalo",
    "tulip",
    "zulitor",
    "rosulator",
    "delipic",
    "mevalotin",
    "vytorin",
    "atozet",
    "cretrol",
    "caduet",
    "pravafen",
    "atorcal",
    "agitin",
    "dehypotin",
    "rosis"
  ),
  collapse = "|"
)

unclassified_lipid_name_qc <- dt_med_classified %>%
  dplyr::filter(
    stringr::str_detect(
      drug_name_lower,
      lipid_name_review_pattern
    ),
    !has_statin,
    !has_ezetimibe_keyword,
    !is_rosis
  ) %>%
  dplyr::distinct(藥品名稱) %>%
  dplyr::arrange(藥品名稱)

rosis_exclusion <- dt_med_classified %>%
  dplyr::filter(is_rosis) %>%
  dplyr::summarise(
    rosis_record_n = dplyr::n(),
    rosis_patient_n = dplyr::n_distinct(病歷號),
    rosis_visit_n = dplyr::n_distinct(
      paste(病歷號, 就診號)
    ),
    rosis_misclassified_as_statin_n = sum(
      has_statin,
      na.rm = TRUE
    )
  )

if (nrow(rosis_exclusion) == 0L) {
  rosis_exclusion <- tibble::tibble(
    rosis_record_n = 0L,
    rosis_patient_n = 0L,
    rosis_visit_n = 0L,
    rosis_misclassified_as_statin_n = 0L
  )
}

if (
  nrow(known_statin_brand_qc) > 0L &&
  any(!known_statin_brand_qc$has_statin)
) {
  stop(
    "已知含史他汀商品未被正確辨識，",
    "請檢查 known_statin_brand_qc 與 statin_pattern。"
  )
}

if (nrow(unclassified_lipid_name_qc) > 0L) {
  warning(
    "仍有疑似降血脂藥品名稱未被分類；",
    "請檢查 unclassified_lipid_qc 工作表。"
  )
}

if (
  rosis_exclusion$rosis_misclassified_as_statin_n[1L] > 0L
) {
  stop("Rosis 仍被誤判為史他汀，請停止分析。")
}

ezetimibe_drug_check <- dt_med_classified %>%
  dplyr::filter(has_ezetimibe_keyword) %>%
  dplyr::distinct(
    藥品名稱,
    has_statin
  ) %>%
  dplyr::arrange(藥品名稱)

visit_med_class <- dt_med_classified %>%
  dplyr::group_by(病歷號, 就診號) %>%
  dplyr::summarise(
    any_statin = any(has_statin, na.rm = TRUE),
    any_ezetimibe_single = any(
      has_ezetimibe_keyword,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Group_Visit = dplyr::case_when(
      any_statin ~ "Statin",
      !any_statin & any_ezetimibe_single ~ "Ezetimibe",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Group_Visit)) %>%
  dplyr::select(
    病歷號,
    就診號,
    Group_Visit
  )

dt_exposure <- dt_adult_type2 %>%
  dplyr::left_join(
    visit_med_class,
    by = c("病歷號", "就診號"),
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    Group_Visit = dplyr::if_else(
      is.na(Group_Visit),
      "Control",
      Group_Visit
    )
  )

dt_exposure_check <- dt_exposure %>%
  dplyr::group_by(病歷號) %>%
  dplyr::summarise(
    Unique_Group_Count = dplyr::n_distinct(
      Group_Visit,
      na.rm = TRUE
    ),
    Group_Set = paste(
      sort(unique(Group_Visit)),
      collapse = "+"
    ),
    observation_n = dplyr::n(),
    .groups = "drop"
  )

fixed_exposure_ids <- dt_exposure_check %>%
  dplyr::filter(
    Unique_Group_Count == 1L,
    Group_Set %in% c("Control", "Statin")
  ) %>%
  dplyr::select(病歷號)

excluded_exposure_summary <- dt_exposure_check %>%
  dplyr::mutate(
    exposure_status = dplyr::case_when(
      病歷號 %in% fixed_exposure_ids$病歷號 ~
        "Included: fixed Control/Statin",
      Group_Set == "Ezetimibe" ~
        "Excluded: fixed Ezetimibe only",
      stringr::str_detect(Group_Set, "\\+") ~
        "Excluded: mixed or time-varying exposure",
      TRUE ~
        "Excluded: other"
    )
  ) %>%
  dplyr::group_by(exposure_status) %>%
  dplyr::summarise(
    patient_n = dplyr::n(),
    observation_n = sum(observation_n),
    .groups = "drop"
  ) %>%
  dplyr::arrange(exposure_status)

ezetimibe_exposure_check <- dt_exposure_check %>%
  dplyr::filter(
    stringr::str_detect(
      Group_Set,
      "Ezetimibe"
    )
  )

dt_fixed_exposure <- dt_exposure %>%
  dplyr::semi_join(
    fixed_exposure_ids,
    by = "病歷號"
  ) %>%
  dplyr::group_by(病歷號) %>%
  dplyr::mutate(
    Group = dplyr::first(Group_Visit)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-Group_Visit) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c("Control", "Statin")
    )
  )

fixed_exposure_group_patient <- dt_fixed_exposure %>%
  dplyr::distinct(病歷號, Group) %>%
  dplyr::count(Group, name = "patient_n")

fixed_exposure_group_observation <- dt_fixed_exposure %>%
  dplyr::count(Group, name = "observation_n")

message("成人 Type 2 固定暴露樣本（病患數）：")
print(fixed_exposure_group_patient)
# A tibble: 2 × 2
# Group   patient_n
# <fct>       <int>
# 1 Control      1423
# 2 Statin       1191
message("成人 Type 2 固定暴露樣本（觀測筆數）：")
print(fixed_exposure_group_observation)
# A tibble: 2 × 2
# Group   observation_n
# <fct>           <int>
# 1 Control          6877
# 2 Statin          11256

# ============================================================
# 4. 欄位名稱／型別一致化與共病
# ============================================================

dt_harmonized <- dt_fixed_exposure %>%
  dplyr::rename(
    ID = 病歷號,
    Sex_raw = 性別,
    HbA1c = HbA1C,
    Glu_AC = Glucose_AC,
    HDL = `HDL-C`,
    LDL = `LDL-C`,
    TG = Triglycerides,
    TC = Total_cholesterol,
    Weight = 體重,
    Height = 身高,
    HbA1c_left_censored = HbA1C_left_censored,
    HbA1c_right_censored = HbA1C_right_censored
  ) %>%
  dplyr::mutate(
    ID = standardize_id(ID),
    就診號 = standardize_id(就診號),
    報告日期 = parse_date_from_excel(報告日期),
    Sex_numeric = standardize_sex_binary(Sex_raw),
    Sex = factor(
      Sex_numeric,
      levels = c(0, 1),
      labels = c("Female", "Male")
    ),
    Group = factor(
      Group,
      levels = c("Control", "Statin")
    ),
    dplyr::across(
      dplyr::any_of(
        c(
          "HbA1c", "Glu_AC", "HDL", "LDL",
          "TG", "TC", "BMI",
          "Systolic_BP", "Diastolic_BP",
          "Age", "Weight", "Height",
          "BP_time_gap_days",
          "basic_time_gap_days"
        )
      ),
      as_numeric_safely
    )
  )

if (any(is.na(dt_harmonized$Sex))) {
  stop(
    "成人 Type 2 固定暴露資料中仍有 Sex 缺失；",
    "請回查資料串檔 QC。"
  )
}

if (any(is.na(dt_harmonized$Age))) {
  stop(
    "成人資格篩選後仍有 Age 缺失；",
    "請檢查成人篩選邏輯。"
  )
}

dt_diag <- read_xlsx_hospital(diag_path) %>%
  dplyr::mutate(
    ID = standardize_id(病歷號),
    就診號 = standardize_id(就診號)
  )

diag_code_candidates <- c(
  "ICD10診斷碼",
  "ICD-10診斷碼",
  "ICD10",
  "診斷碼",
  "疾病碼"
)

diag_code_col <- intersect(
  diag_code_candidates,
  names(dt_diag)
)[1L]

if (is.na(diag_code_col)) {
  stop("疾病診斷明細中找不到 ICD 診斷碼欄位。")
}

# 依既定論文定義保留，不擅自改變。
icd_patterns <- list(
  CVD_History = "^I(2[0-5]|6[0-9])",
  HTN_History = "^I1[0-5]",
  CKD_History = "^(N1[7-9]|I12)",
  Liver_History = "^K7[0-7]"
)

comorbidity_definition <- tibble::tribble(
  ~variable, ~icd_pattern, ~description,
  "CVD_History", icd_patterns$CVD_History,
  "I20–I25 或 I60–I69",
  "HTN_History", icd_patterns$HTN_History,
  "I10–I15",
  "CKD_History", icd_patterns$CKD_History,
  "依目前研究既定定義：N17–N19 或 I12",
  "Liver_History", icd_patterns$Liver_History,
  "K70–K77"
)

dt_key_final_patient <- dt_harmonized %>%
  dplyr::distinct(ID)

dt_comorbidities <- dt_diag %>%
  dplyr::semi_join(
    dt_key_final_patient,
    by = "ID"
  ) %>%
  dplyr::mutate(
    icd_code = stringr::str_to_upper(
      stringr::str_replace_all(
        stringr::str_squish(
          as.character(.data[[diag_code_col]])
        ),
        "\\.",
        ""
      )
    ),
    CVD_flag = stringr::str_detect(
      icd_code,
      icd_patterns$CVD_History
    ),
    HTN_flag = stringr::str_detect(
      icd_code,
      icd_patterns$HTN_History
    ),
    CKD_flag = stringr::str_detect(
      icd_code,
      icd_patterns$CKD_History
    ),
    Liver_flag = stringr::str_detect(
      icd_code,
      icd_patterns$Liver_History
    )
  ) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    CVD_History = as.integer(
      any(CVD_flag, na.rm = TRUE)
    ),
    HTN_History = as.integer(
      any(HTN_flag, na.rm = TRUE)
    ),
    CKD_History = as.integer(
      any(CKD_flag, na.rm = TRUE)
    ),
    Liver_History = as.integer(
      any(Liver_flag, na.rm = TRUE)
    ),
    .groups = "drop"
  )

dt_harmonized <- dt_harmonized %>%
  dplyr::left_join(
    dt_comorbidities,
    by = "ID",
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        CVD_History,
        HTN_History,
        CKD_History,
        Liver_History
      ),
      ~ tidyr::replace_na(.x, 0L)
    ),
    dplyr::across(
      c(
        CVD_History,
        HTN_History,
        CKD_History,
        Liver_History
      ),
      ~ factor(.x, levels = c(0, 1))
    )
  )

# ============================================================
# 5. 同一病患、同一報告日期之多筆紀錄彙整
# ============================================================
# 同一病患在同一報告日期可能因不同就診號或資料來源而出現多列。
# 為避免相同實際時間被 GEE/LME 或序列模型視為多個時間步，
# 在 MICE 與所有後續模型之前先彙整為單一病患－日期紀錄。
#
# 彙整原則：
# - HbA1c、Glu-AC、HDL、LDL、TG、TC、BMI、身高、體重與血壓：
#   取該日非缺失值平均。
# - HbA1c 左／右設限旗標：任一來源列為 1 即記為 1。
# - Group、Sex、Age 與共病等固定變項必須在同一病患同日內一致；
#   若不一致則停止執行。
# - 就診號保留排序後第一個值，另以 source_visit_ids 保留完整來源。
# - 此步驟發生在 MICE 前，因此推論與預測資料使用相同的時間點定義。

same_date_duplicate_before_qc <- dt_harmonized %>%
  dplyr::group_by(ID, 報告日期) %>%
  dplyr::summarise(
    source_row_n = dplyr::n(),
    source_visit_n = dplyr::n_distinct(就診號, na.rm = TRUE),
    source_visit_ids = paste(
      sort(unique(stats::na.omit(as.character(就診號)))),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(source_row_n > 1L) %>%
  dplyr::arrange(ID, 報告日期)

same_date_fixed_vars <- intersect(
  c(
    "Group", "Sex", "Sex_numeric", "Sex_raw", "Age",
    "diabetes_type", "CVD_History", "HTN_History",
    "CKD_History", "Liver_History"
  ),
  names(dt_harmonized)
)

same_date_fixed_constancy_qc <- purrr::map_dfr(
  same_date_fixed_vars,
  function(v) {
    dt_harmonized %>%
      dplyr::group_by(ID, 報告日期) %>%
      dplyr::summarise(
        source_row_n = dplyr::n(),
        distinct_value_n = dplyr::n_distinct(.data[[v]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        source_row_n > 1L,
        distinct_value_n > 1L
      ) %>%
      dplyr::mutate(variable = v) %>%
      dplyr::select(
        variable,
        ID,
        報告日期,
        source_row_n,
        distinct_value_n
      )
  }
)

if (nrow(same_date_fixed_constancy_qc) > 0L) {
  print(same_date_fixed_constancy_qc)
  stop(
    "同一病患同一報告日期的固定變項不一致，",
    "不可直接彙整；請檢查 same_date_fixed_constancy_qc。"
  )
}

same_date_mean_vars <- intersect(
  c(
    "HbA1c", "Glu_AC", "HDL", "LDL", "TG", "TC",
    "BMI", "Weight", "Height", "Systolic_BP", "Diastolic_BP"
  ),
  names(dt_harmonized)
)

same_date_any_vars <- intersect(
  c(
    "HbA1c_left_censored",
    "HbA1c_right_censored"
  ),
  names(dt_harmonized)
)

same_date_key_vars <- c("ID", "報告日期")

same_date_first_vars <- setdiff(
  names(dt_harmonized),
  c(
    same_date_key_vars,
    same_date_mean_vars,
    same_date_any_vars,
    "就診號"
  )
)

dt_same_date_aggregated <- dt_harmonized %>%
  dplyr::arrange(ID, 報告日期, 就診號) %>%
  dplyr::group_by(ID, 報告日期) %>%
  dplyr::summarise(
    source_row_n = dplyr::n(),
    source_visit_n = dplyr::n_distinct(就診號, na.rm = TRUE),
    source_visit_ids = paste(
      sort(unique(stats::na.omit(as.character(就診號)))),
      collapse = " | "
    ),
    就診號 = first_nonmissing(就診號),
    dplyr::across(
      dplyr::all_of(same_date_mean_vars),
      mean_or_na
    ),
    dplyr::across(
      dplyr::all_of(same_date_any_vars),
      ~ as.numeric(any(as_numeric_safely(.x) == 1, na.rm = TRUE))
    ),
    dplyr::across(
      dplyr::all_of(same_date_first_vars),
      first_nonmissing
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Control", "Statin")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    dplyr::across(
      c(CVD_History, HTN_History, CKD_History, Liver_History),
      ~ factor(as.character(.x), levels = c("0", "1"))
    )
  ) %>%
  dplyr::arrange(ID, 報告日期, 就診號)

same_date_duplicate_after_qc <- dt_same_date_aggregated %>%
  dplyr::count(ID, 報告日期, name = "duplicate_n") %>%
  dplyr::filter(duplicate_n > 1L)

same_date_aggregation_summary <- tibble::tibble(
  item = c(
    "rows_before_aggregation",
    "rows_after_aggregation",
    "rows_removed_by_aggregation",
    "patient_date_cells_aggregated",
    "patients_with_aggregation",
    "remaining_same_date_duplicates"
  ),
  value = c(
    nrow(dt_harmonized),
    nrow(dt_same_date_aggregated),
    nrow(dt_harmonized) - nrow(dt_same_date_aggregated),
    nrow(same_date_duplicate_before_qc),
    dplyr::n_distinct(same_date_duplicate_before_qc$ID),
    nrow(same_date_duplicate_after_qc)
  )
)

if (nrow(same_date_duplicate_after_qc) > 0L) {
  stop("同日彙整後仍存在同一病患同一報告日期的重複列。")
}

# ============================================================
# 6. 異常值 QA：只標記，不改值
# ============================================================

dt_flagged <- add_quality_flags(dt_same_date_aggregated)
quality_flag_summary <- summarise_quality_flags(dt_flagged)

# ============================================================
# 7. 建立 T1–T19 與最終縱向分析 cohort
# ============================================================

dt_model_pre_filter <- dt_flagged %>%
  dplyr::arrange(ID, 報告日期, 就診號) %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(
    Base_Report_Date = dplyr::first(報告日期),
    Time_in_Days = as.numeric(
      difftime(
        報告日期,
        Base_Report_Date,
        units = "days"
      )
    ),
    Time_in_Years = Time_in_Days / 365.25,
    Time_in_Months = Time_in_Years * 12,
    Time_Cat_Index_3M = floor(Time_in_Months / 3),
    Time_Category_raw = paste0(
      "T",
      Time_Cat_Index_3M + 1L
    )
  ) %>%
  dplyr::ungroup()

time_category_count_before <- dt_model_pre_filter %>%
  dplyr::group_by(Time_Category_raw) %>%
  dplyr::summarise(
    observation_n = dplyr::n(),
    patient_n = dplyr::n_distinct(ID),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    time_num = as.integer(
      stringr::str_remove(
        Time_Category_raw,
        "T"
      )
    )
  ) %>%
  dplyr::arrange(time_num)

time_category_count_by_group_before <- dt_model_pre_filter %>%
  dplyr::group_by(
    Time_Category_raw,
    Group
  ) %>%
  dplyr::summarise(
    observation_n = dplyr::n(),
    patient_n = dplyr::n_distinct(ID),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    time_num = as.integer(
      stringr::str_remove(
        Time_Category_raw,
        "T"
      )
    )
  ) %>%
  dplyr::arrange(time_num, Group)

# 原始追蹤資料的時間類別範圍與 T20/T21 排除紀錄。
time_category_range_qc <- tibble::tibble(
  min_time_num = min(
    time_category_count_before$time_num,
    na.rm = TRUE
  ),
  max_time_num = max(
    time_category_count_before$time_num,
    na.rm = TRUE
  ),
  observed_time_categories = paste(
    time_category_count_before$Time_Category_raw,
    collapse = " | "
  )
)

time_category_decision <- time_category_count_before %>%
  dplyr::mutate(
    analysis_status = dplyr::if_else(
      time_num <= MAX_TIME_CATEGORY_INDEX + 1L,
      "Included: T1-T19",
      "Excluded: after T19"
    ),
    exclusion_reason = dplyr::case_when(
      time_num == 20L ~
        "T20 observation count is sparse",
      time_num == 21L ~
        "T21 observation count is sparse",
      time_num > 21L ~
        "Beyond the prespecified T1-T19 analysis window",
      TRUE ~
        ""
    )
  )

time_category_excluded_after_t19 <- time_category_decision %>%
  dplyr::filter(
    time_num > MAX_TIME_CATEGORY_INDEX + 1L
  )

time_category_exclusion_totals <- time_category_decision %>%
  dplyr::group_by(analysis_status) %>%
  dplyr::summarise(
    observation_n = sum(observation_n),
    patient_n_sum_across_periods = sum(patient_n),
    .groups = "drop"
  )

if (
  max(
    time_category_count_before$time_num,
    na.rm = TRUE
  ) > 21L
) {
  warning(
    "原始資料出現 T21 之後的時間類別；",
    "目前仍依預先設定只保留 T1-T19，",
    "請檢查 time_category_decision 工作表。"
  )
}

dt_model_t19 <- dt_model_pre_filter %>%
  dplyr::filter(
    Time_Cat_Index_3M >= 0,
    Time_Cat_Index_3M <= MAX_TIME_CATEGORY_INDEX
  ) %>%
  dplyr::mutate(
    Time_Category = factor(
      paste0("T", Time_Cat_Index_3M + 1L),
      levels = paste0("T", 1:19),
      ordered = TRUE
    )
  )

longitudinal_eligibility <- dt_model_t19 %>%
  dplyr::group_by(ID, Group) %>%
  dplyr::summarise(
    observation_n = dplyr::n(),
    distinct_report_date_n = dplyr::n_distinct(報告日期),
    first_report_date = min(報告日期),
    last_report_date = max(報告日期),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    longitudinal_status = dplyr::if_else(
      distinct_report_date_n >= 2L,
      "Included: at least 2 distinct dates",
      "Excluded: fewer than 2 distinct dates"
    )
  )

longitudinal_eligibility_summary <- longitudinal_eligibility %>%
  dplyr::group_by(longitudinal_status) %>%
  dplyr::summarise(
    patient_n = dplyr::n(),
    observation_n = sum(observation_n),
    .groups = "drop"
  )

eligible_longitudinal_ids <- longitudinal_eligibility %>%
  dplyr::filter(distinct_report_date_n >= 2L) %>%
  dplyr::select(ID)

dt_model <- dt_model_t19 %>%
  dplyr::semi_join(
    eligible_longitudinal_ids,
    by = "ID"
  ) %>%
  dplyr::arrange(ID, 報告日期, 就診號) %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(
    Age_baseline = dplyr::first(Age),
    BMI_baseline_raw = dplyr::first(BMI)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    row_id = dplyr::row_number(),
    HbA1c_observed = !is.na(HbA1c),
    Glu_AC_observed = !is.na(Glu_AC),
    LDL_observed = !is.na(LDL),
    TG_observed = !is.na(TG),
    BMI_observed = !is.na(BMI),
    Group = factor(
      Group,
      levels = c("Control", "Statin")
    ),
    Sex = factor(
      Sex,
      levels = c("Female", "Male")
    ),
    Time_Category = factor(
      Time_Category,
      levels = paste0("T", 1:19),
      ordered = TRUE
    )
  )

if (any(is.na(dt_model$Age_baseline))) {
  stop("最終分析 cohort 仍有 Age_baseline 缺失。")
}

if (any(dt_model$Age_baseline < 18, na.rm = TRUE)) {
  stop("最終分析 cohort 仍含基準年齡未滿 18 歲者。")
}

repeat_sequence_qc <- dt_model %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    observation_n = dplyr::n(),
    distinct_report_date_n = dplyr::n_distinct(報告日期),
    .groups = "drop"
  ) %>%
  dplyr::filter(distinct_report_date_n < 2L)

final_event_duplicate_qc <- dt_model %>%
  dplyr::count(
    ID,
    報告日期,
    就診號,
    name = "duplicate_n"
  ) %>%
  dplyr::filter(duplicate_n > 1L)

same_patient_same_time_qc <- dt_model %>%
  dplyr::count(
    ID,
    Time_in_Years,
    name = "duplicate_n"
  ) %>%
  dplyr::filter(duplicate_n > 1L)

same_patient_date_multi_visit_qc <- dt_model %>%
  dplyr::group_by(ID, 報告日期) %>%
  dplyr::summarise(
    observation_n = dplyr::n(),
    distinct_visit_n = dplyr::n_distinct(就診號),
    visit_ids = paste(
      sort(unique(as.character(就診號))),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(observation_n > 1L)

age_baseline_constancy_qc <- dt_model %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    age_baseline_n_distinct = dplyr::n_distinct(
      Age_baseline
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(age_baseline_n_distinct != 1L)

patient_flow <- tibble::tibble(
  step = c(
    "資料串檔後全部糖尿病相關檢驗紀錄",
    "限定 Type 2 diabetes",
    "排除基準年齡 <18 歲或年齡缺失",
    "限定固定 Control/Statin 暴露",
    "合併同一病患同一報告日期之多筆紀錄",
    "保留 T1–T19",
    "保留至少兩個不同追蹤日期"
  ),
  patient_n = c(
    dplyr::n_distinct(dt_raw$病歷號),
    dplyr::n_distinct(dt_type2$病歷號),
    dplyr::n_distinct(dt_adult_type2$病歷號),
    dplyr::n_distinct(dt_fixed_exposure$病歷號),
    dplyr::n_distinct(dt_same_date_aggregated$ID),
    dplyr::n_distinct(dt_model_t19$ID),
    dplyr::n_distinct(dt_model$ID)
  ),
  observation_n = c(
    nrow(dt_raw),
    nrow(dt_type2),
    nrow(dt_adult_type2),
    nrow(dt_fixed_exposure),
    nrow(dt_same_date_aggregated),
    nrow(dt_model_t19),
    nrow(dt_model)
  )
)

group_count_final_patient <- dt_model %>%
  dplyr::distinct(ID, Group) %>%
  dplyr::count(Group, name = "patient_n")

group_count_final_observation <- dt_model %>%
  dplyr::count(Group, name = "observation_n")

repeat_patient_summary_final <- dt_model %>%
  dplyr::count(ID, name = "n_obs") %>%
  dplyr::summarise(
    patient_n = dplyr::n(),
    mean_n_obs = mean(n_obs),
    sd_n_obs = stats::sd(n_obs),
    median_n_obs = stats::median(n_obs),
    q1_n_obs = safe_quantile(n_obs, 0.25),
    q3_n_obs = safe_quantile(n_obs, 0.75),
    min_n_obs = min(n_obs),
    max_n_obs = max(n_obs)
  )

time_category_count_final <- dt_model %>%
  dplyr::group_by(Time_Category, Group) %>%
  dplyr::summarise(
    observation_n = dplyr::n(),
    patient_n = dplyr::n_distinct(ID),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Time_Category = as.character(Time_Category)
  )

message("最終成人 Type 2 主分析樣本：")
message("病患數：", dplyr::n_distinct(dt_model$ID))
# 病患數：1790
message("觀測筆數：", nrow(dt_model))
# 觀測筆數：17252
print(table(dt_model$Group))
# Control  Statin 
# 6238   11014 

# ------------------------------------------------------------
# 7.1 後續模型轉換規則與 log 可用性 QC
# ------------------------------------------------------------

# 在這保留原始尺度與原始單位，不在此處建立永久 log 或 z-score 欄位。
# 本程式中的 standardize／standardized 僅指欄位格式與編碼一致化，不代表 z-score。
# GEE / LME：於後續模型程式中對 outcome 取自然對數；共變數維持原始尺度。
# DL：先依既定 split 區分 train/validation/test，再僅用 train set 估計 log 後
#     之平均值與標準差，並將相同參數套用至 validation / test；評估前轉回原始尺度。
transformation_rules <- tibble::tribble(
  ~analysis_stage, ~outcome_scale, ~continuous_covariate_scale, ~parameter_source, ~evaluation_scale,
  "Descriptive tables and figures",
  "Original scale",
  "Original scale",
  "Not applicable",
  "Original scale",
  "GEE / LME",
  "Natural-log outcome",
  "Original covariate scale unless the model script explicitly specifies otherwise",
  "Not applicable",
  "Predictions transformed back to original scale",
  "Deep learning",
  "Natural-log outcome followed by z-score",
  "Time, baseline age and BMI standardized by training-set parameters",
  "Training set only",
  "Predictions inverse-standardized and exponentiated to original scale"
)

model_data_usage <- tibble::tribble(
  ~analysis_component, ~model_fitting_data, ~validation_role, ~evaluation_data, ~inference_or_metric,
  "GEE / LME statistical inference",
  "full eligible cohort across all split labels",
  "not used",
  "not applicable",
  "Coefficients, standard errors, 95% CI and p-values from full-cohort models; pooled across 10 imputations by Rubin's rules",
  "GEE / LME prediction",
  "train set",
  "not used",
  "test set",
  "RMSE, MAE, MAPE and Bias on the original scale; reported overall and by Control / Statin group",
  "Deep learning prediction",
  "train set",
  "validation set for tuning and early stopping",
  "test set",
  "RMSE, MAE, MAPE and Bias on the original scale; reported overall and by Control / Statin group"
)

imputation_split_scope <- tibble::tribble(
  ~item, ~setting, ~interpretation,
  
  "Inference imputation scope",
  "The final eligible cohort is multiply imputed into 10 completed datasets",
  "These 10 imputed datasets are used only for full-cohort GEE and LME statistical inference",
  
  "Prediction data construction",
  "A single prediction-safe dataset is created after fixing the patient-level train, validation and test split",
  "HbA1c, Glu-AC, LDL and TG are not imputed across data splits; only observed historical outcome values are used",
  
  "BMI preprocessing for prediction",
  "Missing time-varying BMI values are filled using group-specific medians estimated from the training set only",
  "Validation and test data do not contribute to the BMI imputation parameters",
  
  "Patient separation",
  "Each patient belongs to exactly one of train, validation or test",
  "No patient appears in more than one split",
  
  "Thesis wording",
  "Full-cohort multiple imputation is used for GEE and LME inference, whereas prediction models use a separate prediction-safe dataset",
  "Do not describe the prediction analysis as using the 10 full-cohort imputed datasets"
)

log_eligibility_qc <- purrr::map_dfr(
  primary_outcomes,
  function(v) {
    x <- dt_model[[v]]

    tibble::tibble(
      variable = v,
      observed_n = sum(!is.na(x)),
      zero_n = sum(x == 0, na.rm = TRUE),
      negative_n = sum(x < 0, na.rm = TRUE),
      nonpositive_n = sum(x <= 0, na.rm = TRUE),
      min_observed = ifelse(
        any(!is.na(x)),
        min(x, na.rm = TRUE),
        NA_real_
      )
    )
  }
)

if (any(log_eligibility_qc$nonpositive_n > 0L)) {
  stop(
    "主要代謝指標中存在 <= 0 的觀測值，無法直接取自然對數；",
    "請檢查 log_eligibility_qc。"
  )
}

# ============================================================
# 8. 敘述統計、缺失率與論文表格
# ============================================================

safe_wilcox_p <- function(x, g) {
  ok <- !is.na(x) & is.finite(x) & !is.na(g)

  if (
    sum(ok) == 0L ||
      length(unique(g[ok])) < 2L
  ) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      x[ok] ~ g[ok],
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

format_mean_sd <- function(mean_value, sd_value, digits = 1L) {
  if (
    is.na(mean_value) ||
      is.na(sd_value)
  ) {
    return(NA_character_)
  }

  paste0(
    formatC(
      mean_value,
      format = "f",
      digits = digits
    ),
    " ± ",
    formatC(
      sd_value,
      format = "f",
      digits = digits
    )
  )
}

format_n_percent <- function(n, denominator, digits = 1L) {
  if (
    is.na(n) ||
      is.na(denominator) ||
      denominator <= 0
  ) {
    return(NA_character_)
  }

  paste0(
    n,
    " (",
    formatC(
      n / denominator * 100,
      format = "f",
      digits = digits
    ),
    "%)"
  )
}

format_missing_n_percent <- function(
    missing_n,
    total_n,
    digits = 1L
) {
  if (
    is.na(missing_n) ||
      is.na(total_n) ||
      total_n <= 0
  ) {
    return(NA_character_)
  }

  paste0(
    missing_n,
    " (",
    formatC(
      missing_n / total_n * 100,
      format = "f",
      digits = digits
    ),
    "%)"
  )
}

format_p_for_thesis <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p > 0.99 ~ ">0.99",
    TRUE ~ formatC(
      p,
      format = "f",
      digits = 3
    )
  )
}

# ------------------------------------------------------------
# 7.1 事件層級與病患－時間區間層級缺失率
# ------------------------------------------------------------

missing_response_overall <- make_missing_summary(
  dt_model,
  response_vars
)

missing_response_by_group <- make_missing_summary(
  dt_model,
  response_vars,
  group_vars = "Group"
)

missing_covariate_overall <- make_missing_summary(
  dt_model,
  c(
    "Age_baseline",
    "Sex",
    "BMI_baseline_raw",
    "BMI",
    "Systolic_BP",
    "Diastolic_BP"
  )
)

missing_covariate_by_group <- make_missing_summary(
  dt_model,
  c(
    "Age_baseline",
    "Sex",
    "BMI_baseline_raw",
    "BMI",
    "Systolic_BP",
    "Diastolic_BP"
  ),
  group_vars = "Group"
)

# 同一病患於同一三個月區間若有多筆檢驗，先取病患內平均，
# 避免同一病患在同一區間重複貢獻多次。
patient_period_outcomes <- dt_model %>%
  dplyr::select(
    ID,
    Group,
    Time_Category,
    dplyr::all_of(response_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(response_vars),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  dplyr::group_by(
    ID,
    Group,
    Time_Category,
    outcome
  ) %>%
  dplyr::summarise(
    value = if (
      all(is.na(value))
    ) {
      NA_real_
    } else {
      mean(value, na.rm = TRUE)
    },
    .groups = "drop"
  )

missing_response_by_group_time <- patient_period_outcomes %>%
  dplyr::filter(
    outcome %in% primary_outcomes
  ) %>%
  dplyr::group_by(
    Group,
    Time_Category,
    outcome
  ) %>%
  dplyr::summarise(
    total_patient_n = dplyr::n(),
    missing_patient_n = sum(is.na(value)),
    observed_patient_n = sum(!is.na(value)),
    missing_rate_percent = round(
      missing_patient_n / total_patient_n * 100,
      2
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Time_Category = as.character(Time_Category),
    Time_Num = as.integer(
      stringr::str_remove(
        Time_Category,
        "T"
      )
    )
  ) %>%
  dplyr::arrange(outcome, Time_Num, Group)

# ------------------------------------------------------------
# 7.2 圖 4.1 所需平均值、SD、SE 與 95% CI
# ------------------------------------------------------------

make_mean_sd_ci_patient_period <- function(
    dat_long,
    vars
) {
  dat_long %>%
    dplyr::filter(
      outcome %in% vars
    ) %>%
    dplyr::group_by(
      Time_Category,
      Group,
      outcome
    ) %>%
    dplyr::summarise(
      patient_n = sum(!is.na(value)),
      mean = safe_mean(value),
      sd = safe_sd(value),
      se = dplyr::if_else(
        patient_n > 1L,
        sd / sqrt(patient_n),
        NA_real_
      ),
      t_critical = critical_t_975(patient_n),
      ci_lower = dplyr::if_else(
        patient_n > 1L,
        mean - t_critical * se,
        NA_real_
      ),
      ci_upper = dplyr::if_else(
        patient_n > 1L,
        mean + t_critical * se,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Time_Category = as.character(Time_Category),
      Time_Num = as.integer(
        stringr::str_remove(
          Time_Category,
          "T"
        )
      )
    ) %>%
    dplyr::arrange(
      outcome,
      Time_Num,
      Group
    )
}

mean_sd_ci_figure4_1 <- make_mean_sd_ci_patient_period(
  patient_period_outcomes,
  primary_outcomes
)

# ------------------------------------------------------------
# 7.3 表 4.2－4.5 與附錄 A1－A2
# ------------------------------------------------------------

make_outcome_time_results <- function(
    dat_long,
    outcome_name
) {
  outcome_label <- unname(
    outcome_labels[outcome_name]
  )

  if (is.na(outcome_label)) {
    outcome_label <- outcome_name
  }

  dat_one <- dat_long %>%
    dplyr::filter(
      outcome == outcome_name
    )

  p_tbl <- dat_one %>%
    dplyr::group_by(Time_Category) %>%
    dplyr::summarise(
      p_value = safe_wilcox_p(
        value,
        Group
      ),
      p_value_formatted = format_p_for_thesis(
        p_value
      ),
      test = "Wilcoxon rank-sum test",
      .groups = "drop"
    )

  overall_tbl <- dat_one %>%
    dplyr::group_by(Time_Category) %>%
    dplyr::summarise(
      Group = "Overall",
      total_patient_n = dplyr::n(),
      observed_patient_n = sum(!is.na(value)),
      missing_patient_n = sum(is.na(value)),
      missing_rate_percent = round(
        missing_patient_n / total_patient_n * 100,
        2
      ),
      mean = safe_mean(value),
      sd = safe_sd(value),
      se = dplyr::if_else(
        observed_patient_n > 1L,
        sd / sqrt(observed_patient_n),
        NA_real_
      ),
      t_critical = critical_t_975(
        observed_patient_n
      ),
      ci_lower = dplyr::if_else(
        observed_patient_n > 1L,
        mean - t_critical * se,
        NA_real_
      ),
      ci_upper = dplyr::if_else(
        observed_patient_n > 1L,
        mean + t_critical * se,
        NA_real_
      ),
      median = safe_median(value),
      q1 = safe_quantile(value, 0.25),
      q3 = safe_quantile(value, 0.75),
      .groups = "drop"
    )

  group_tbl <- dat_one %>%
    dplyr::group_by(
      Time_Category,
      Group
    ) %>%
    dplyr::summarise(
      total_patient_n = dplyr::n(),
      observed_patient_n = sum(!is.na(value)),
      missing_patient_n = sum(is.na(value)),
      missing_rate_percent = round(
        missing_patient_n / total_patient_n * 100,
        2
      ),
      mean = safe_mean(value),
      sd = safe_sd(value),
      se = dplyr::if_else(
        observed_patient_n > 1L,
        sd / sqrt(observed_patient_n),
        NA_real_
      ),
      t_critical = critical_t_975(
        observed_patient_n
      ),
      ci_lower = dplyr::if_else(
        observed_patient_n > 1L,
        mean - t_critical * se,
        NA_real_
      ),
      ci_upper = dplyr::if_else(
        observed_patient_n > 1L,
        mean + t_critical * se,
        NA_real_
      ),
      median = safe_median(value),
      q1 = safe_quantile(value, 0.25),
      q3 = safe_quantile(value, 0.75),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Group = as.character(Group)
    )

  long_tbl <- dplyr::bind_rows(
    overall_tbl,
    group_tbl
  ) %>%
    dplyr::left_join(
      p_tbl,
      by = "Time_Category",
      relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      outcome = outcome_name,
      outcome_label = outcome_label,
      Time_Category = as.character(Time_Category),
      Time_Num = as.integer(
        stringr::str_remove(
          Time_Category,
          "T"
        )
      ),
      Group = factor(
        Group,
        levels = c(
          "Overall",
          "Statin",
          "Control"
        )
      )
    ) %>%
    dplyr::arrange(
      Time_Num,
      Group
    )

  summary_wide <- long_tbl %>%
    dplyr::mutate(
      mean_sd = purrr::map2_chr(
        mean,
        sd,
        ~ format_mean_sd(
          .x,
          .y,
          digits = 1L
        )
      )
    ) %>%
    dplyr::select(
      Time_Category,
      Group,
      mean_sd
    ) %>%
    tidyr::pivot_wider(
      names_from = Group,
      values_from = mean_sd
    )

  overall_missing <- long_tbl %>%
    dplyr::filter(
      Group == "Overall"
    ) %>%
    dplyr::transmute(
      Time_Category,
      Time_Num,
      available_n = observed_patient_n,
      missing_n = missing_patient_n,
      missing_percent = missing_rate_percent,
      missing_n_percent = paste0(
        missing_patient_n,
        " (",
        formatC(
          missing_rate_percent,
          format = "f",
          digits = 1
        ),
        "%)"
      ),
      p_value,
      p_value_formatted,
      test
    )

  ready_tbl <- overall_missing %>%
    dplyr::left_join(
      summary_wide,
      by = "Time_Category",
      relationship = "one-to-one"
    ) %>%
    dplyr::transmute(
      outcome = outcome_name,
      outcome_label = outcome_label,
      Time_Category,
      Time_Num,
      available_n,
      missing_n,
      missing_percent,
      missing_n_percent,
      Overall,
      Statin,
      Control,
      p_value,
      p_value_formatted,
      test
    ) %>%
    dplyr::arrange(Time_Num)

  list(
    detail = long_tbl,
    ready = ready_tbl
  )
}

table4_2_result <- make_outcome_time_results(
  patient_period_outcomes,
  "HbA1c"
)

table4_3_result <- make_outcome_time_results(
  patient_period_outcomes,
  "Glu_AC"
)

table4_4_result <- make_outcome_time_results(
  patient_period_outcomes,
  "LDL"
)

table4_5_result <- make_outcome_time_results(
  patient_period_outcomes,
  "TG"
)

appendix_A1_result <- make_outcome_time_results(
  patient_period_outcomes,
  "HDL"
)

appendix_A2_result <- make_outcome_time_results(
  patient_period_outcomes,
  "TC"
)

table4_2_HbA1c <- table4_2_result$ready
table4_3_Glu_AC <- table4_3_result$ready
table4_4_LDL <- table4_4_result$ready
table4_5_TG <- table4_5_result$ready

table4_2_HbA1c_detail <- table4_2_result$detail
table4_3_Glu_AC_detail <- table4_3_result$detail
table4_4_LDL_detail <- table4_4_result$detail
table4_5_TG_detail <- table4_5_result$detail

appendix_A1_HDL <- appendix_A1_result$ready
appendix_A2_TC <- appendix_A2_result$ready
appendix_A1_HDL_detail <- appendix_A1_result$detail
appendix_A2_TC_detail <- appendix_A2_result$detail

# ------------------------------------------------------------
# 7.4 表 4.1：每位病患第一筆合格觀測
# ------------------------------------------------------------

dt_baseline <- dt_model %>%
  dplyr::arrange(
    ID,
    報告日期,
    就診號
  ) %>%
  dplyr::group_by(ID) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

baseline_cat_vars <- c(
  "Sex",
  "CVD_History",
  "HTN_History",
  "CKD_History",
  "Liver_History"
)

build_baseline_results <- function(
    baseline_dat,
    bmi_var,
    bmi_label = "BMI"
) {
  continuous_vars <- c(
    "Age_baseline",
    bmi_var
  )

  continuous_long <- baseline_dat %>%
    dplyr::select(
      ID,
      Group,
      dplyr::all_of(continuous_vars)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(continuous_vars),
      names_to = "variable",
      values_to = "value"
    )

  cont_overall <- continuous_long %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      Group = "Overall",
      total_patient_n = dplyr::n(),
      observed_n = sum(!is.na(value)),
      missing_n = sum(is.na(value)),
      missing_rate_percent = round(
        missing_n / total_patient_n * 100,
        2
      ),
      mean = safe_mean(value),
      sd = safe_sd(value),
      median = safe_median(value),
      q1 = safe_quantile(value, 0.25),
      q3 = safe_quantile(value, 0.75),
      .groups = "drop"
    )

  cont_group <- continuous_long %>%
    dplyr::group_by(
      variable,
      Group
    ) %>%
    dplyr::summarise(
      total_patient_n = dplyr::n(),
      observed_n = sum(!is.na(value)),
      missing_n = sum(is.na(value)),
      missing_rate_percent = round(
        missing_n / total_patient_n * 100,
        2
      ),
      mean = safe_mean(value),
      sd = safe_sd(value),
      median = safe_median(value),
      q1 = safe_quantile(value, 0.25),
      q3 = safe_quantile(value, 0.75),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Group = as.character(Group)
    )

  cont_summary <- dplyr::bind_rows(
    cont_overall,
    cont_group
  ) %>%
    dplyr::mutate(
      Group = factor(
        Group,
        levels = c(
          "Overall",
          "Statin",
          "Control"
        )
      )
    ) %>%
    dplyr::arrange(
      match(
        variable,
        continuous_vars
      ),
      Group
    )

  cont_p <- purrr::map_dfr(
    continuous_vars,
    function(v) {
      test_info <- continuous_test_info(
        baseline_dat[[v]],
        baseline_dat$Group
      )

      tibble::tibble(
        variable = v
      ) %>%
        dplyr::bind_cols(test_info)
    }
  )

  categorical_long <- baseline_dat %>%
    dplyr::select(
      ID,
      Group,
      dplyr::all_of(baseline_cat_vars)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(baseline_cat_vars),
        as.character
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(baseline_cat_vars),
      names_to = "variable",
      values_to = "value"
    )

  cat_overall_denominator <- categorical_long %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      Group = "Overall",
      total_patient_n = dplyr::n(),
      missing_n = sum(is.na(value)),
      .groups = "drop"
    )

  cat_overall_count <- categorical_long %>%
    dplyr::filter(
      !is.na(value)
    ) %>%
    dplyr::count(
      variable,
      value,
      name = "n"
    ) %>%
    dplyr::mutate(
      Group = "Overall"
    )

  cat_group_denominator <- categorical_long %>%
    dplyr::group_by(
      variable,
      Group
    ) %>%
    dplyr::summarise(
      total_patient_n = dplyr::n(),
      missing_n = sum(is.na(value)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Group = as.character(Group)
    )

  cat_group_count <- categorical_long %>%
    dplyr::filter(
      !is.na(value)
    ) %>%
    dplyr::count(
      variable,
      Group,
      value,
      name = "n"
    ) %>%
    dplyr::mutate(
      Group = as.character(Group)
    )

  cat_summary <- dplyr::bind_rows(
    cat_overall_count,
    cat_group_count
  ) %>%
    dplyr::left_join(
      dplyr::bind_rows(
        cat_overall_denominator,
        cat_group_denominator
      ),
      by = c(
        "variable",
        "Group"
      ),
      relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      percent_of_total_group = round(
        n / total_patient_n * 100,
        2
      ),
      missing_rate_percent = round(
        missing_n / total_patient_n * 100,
        2
      ),
      Group = factor(
        Group,
        levels = c(
          "Overall",
          "Statin",
          "Control"
        )
      )
    ) %>%
    dplyr::arrange(
      match(
        variable,
        baseline_cat_vars
      ),
      value,
      Group
    )

  cat_p <- purrr::map_dfr(
    baseline_cat_vars,
    function(v) {
      test_info <- categorical_test_info(
        baseline_dat[[v]],
        baseline_dat$Group
      )

      tibble::tibble(
        variable = v
      ) %>%
        dplyr::bind_cols(test_info)
    }
  )

  cont_ready_main <- purrr::map_dfr(
    continuous_vars,
    function(v) {
      one <- cont_summary %>%
        dplyr::filter(
          variable == v
        ) %>%
        dplyr::mutate(
          display = purrr::map2_chr(
            mean,
            sd,
            ~ format_mean_sd(
              .x,
              .y,
              digits = 1L
            )
          )
        ) %>%
        dplyr::select(
          Group,
          display
        ) %>%
        tidyr::pivot_wider(
          names_from = Group,
          values_from = display
        )

      p_one <- cont_p %>%
        dplyr::filter(
          variable == v
        )

      tibble::tibble(
        row_type = "continuous",
        variable = v,
        variable_label = dplyr::case_when(
          v == "Age_baseline" ~
            "基準年齡（歲）",
          TRUE ~
            "身體質量指數（BMI）"
        ),
        category = "",
        Overall = one$Overall,
        Statin = one$Statin,
        Control = one$Control,
        p_value = p_one$p_value,
        p_value_formatted = format_p_for_thesis(
          p_one$p_value
        ),
        test = p_one$test
      )
    }
  )

  cont_ready_missing <- purrr::map_dfr(
    continuous_vars,
    function(v) {
      one <- cont_summary %>%
        dplyr::filter(
          variable == v
        ) %>%
        dplyr::mutate(
          display = purrr::map2_chr(
            missing_n,
            total_patient_n,
            format_missing_n_percent
          )
        ) %>%
        dplyr::select(
          Group,
          display,
          missing_n
        ) %>%
        tidyr::pivot_wider(
          names_from = Group,
          values_from = c(
            display,
            missing_n
          )
        )

      if (
        all(
          c(
            one$missing_n_Overall,
            one$missing_n_Statin,
            one$missing_n_Control
          ) == 0
        )
      ) {
        return(
          tibble::tibble(
            row_type = character(),
            variable = character(),
            variable_label = character(),
            category = character(),
            Overall = character(),
            Statin = character(),
            Control = character(),
            p_value = double(),
            p_value_formatted = character(),
            test = character()
          )
        )
      }

      tibble::tibble(
        row_type = "missing",
        variable = v,
        variable_label = "",
        category = "缺失，n (%)",
        Overall = one$display_Overall,
        Statin = one$display_Statin,
        Control = one$display_Control,
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        test = NA_character_
      )
    }
  )

  sex_p <- cat_p %>%
    dplyr::filter(
      variable == "Sex"
    )

  sex_ready <- cat_summary %>%
    dplyr::filter(
      variable == "Sex"
    ) %>%
    dplyr::mutate(
      category = dplyr::recode(
        value,
        Female = "女性",
        Male = "男性"
      ),
      display = purrr::map2_chr(
        n,
        total_patient_n,
        format_n_percent
      )
    ) %>%
    dplyr::select(
      category,
      Group,
      display
    ) %>%
    tidyr::pivot_wider(
      names_from = Group,
      values_from = display
    ) %>%
    dplyr::mutate(
      row_type = "categorical",
      variable = "Sex",
      variable_label = dplyr::if_else(
        dplyr::row_number() == 1L,
        "性別",
        ""
      ),
      p_value = dplyr::if_else(
        dplyr::row_number() == 1L,
        sex_p$p_value[1L],
        NA_real_
      ),
      p_value_formatted = dplyr::if_else(
        dplyr::row_number() == 1L,
        format_p_for_thesis(
          sex_p$p_value[1L]
        ),
        NA_character_
      ),
      test = dplyr::if_else(
        dplyr::row_number() == 1L,
        sex_p$test[1L],
        NA_character_
      )
    ) %>%
    dplyr::select(
      row_type,
      variable,
      variable_label,
      category,
      Overall,
      Statin,
      Control,
      p_value,
      p_value_formatted,
      test
    )

  comorbidity_labels <- c(
    CVD_History = "心血管疾病病史",
    HTN_History = "高血壓病史",
    CKD_History = "慢性腎臟疾病病史",
    Liver_History = "肝臟疾病病史"
  )

  comorbidity_ready <- purrr::map_dfr(
    names(comorbidity_labels),
    function(v) {
      one <- cat_summary %>%
        dplyr::filter(
          variable == v,
          value == "1"
        ) %>%
        dplyr::mutate(
          display = purrr::map2_chr(
            n,
            total_patient_n,
            format_n_percent
          )
        ) %>%
        dplyr::select(
          Group,
          display
        ) %>%
        tidyr::pivot_wider(
          names_from = Group,
          values_from = display
        )

      p_one <- cat_p %>%
        dplyr::filter(
          variable == v
        )

      tibble::tibble(
        row_type = "binary_yes",
        variable = v,
        variable_label = unname(
          comorbidity_labels[v]
        ),
        category = "有",
        Overall = one$Overall,
        Statin = one$Statin,
        Control = one$Control,
        p_value = p_one$p_value,
        p_value_formatted = format_p_for_thesis(
          p_one$p_value
        ),
        test = p_one$test
      )
    }
  )

  ready_tbl <- dplyr::bind_rows(
    cont_ready_main %>%
      dplyr::filter(
        variable == "Age_baseline"
      ),
    sex_ready,
    cont_ready_main %>%
      dplyr::filter(
        variable == bmi_var
      ),
    cont_ready_missing %>%
      dplyr::filter(
        variable == bmi_var
      ),
    comorbidity_ready
  ) %>%
    dplyr::mutate(
      row_order = dplyr::row_number()
    ) %>%
    dplyr::select(
      row_order,
      dplyr::everything()
    )

  list(
    continuous = cont_summary,
    continuous_p = cont_p,
    categorical = cat_summary,
    categorical_p = cat_p,
    ready = ready_tbl
  )
}

baseline_results_raw <- build_baseline_results(
  baseline_dat = dt_baseline,
  bmi_var = "BMI_baseline_raw"
)

baseline_cont_summary <- baseline_results_raw$continuous
baseline_cont_p <- baseline_results_raw$continuous_p
baseline_cat_summary <- baseline_results_raw$categorical
baseline_cat_p <- baseline_results_raw$categorical_p
table4_1_ready <- baseline_results_raw$ready

# ------------------------------------------------------------
# 7.5 血壓附錄缺失率
# ------------------------------------------------------------

bp_vars_available <- intersect(
  c(
    "Systolic_BP",
    "Diastolic_BP"
  ),
  names(dt_model)
)

make_bp_missing_or_note <- function(
    dat,
    group_vars = NULL,
    unit_label
) {
  if (
    length(bp_vars_available) == 0L
  ) {
    return(
      tibble::tibble(
        unit = unit_label,
        status = "not available",
        note = paste0(
          "最新版資料串檔主檔未提供血壓欄位；",
          "本次不計算血壓缺失率，也不將其解讀為 100% 缺失。"
        )
      )
    )
  }

  make_missing_summary(
    dat,
    bp_vars_available,
    group_vars = group_vars
  ) %>%
    dplyr::mutate(
      unit = unit_label,
      status = "available",
      note = ""
    )
}

bp_missing_event_overall <-
  make_bp_missing_or_note(
    dt_model,
    group_vars = NULL,
    unit_label = "event"
  )

bp_missing_event_by_group <-
  make_bp_missing_or_note(
    dt_model,
    group_vars = "Group",
    unit_label = "event"
  )

bp_missing_baseline_overall <-
  make_bp_missing_or_note(
    dt_baseline,
    group_vars = NULL,
    unit_label = "baseline patient"
  )

bp_missing_baseline_by_group <-
  make_bp_missing_or_note(
    dt_baseline,
    group_vars = "Group",
    unit_label = "baseline patient"
  )

# ============================================================
# 9. 圖形：中英文成對輸出
# ============================================================

plot_language_config <- list(
  zh = list(
    suffix = "zh",
    control = "控制組",
    statin = "史他汀組",
    follow_year = "追蹤時間（年）",
    follow_period = "追蹤時間區間",
    group_title = "組別",
    trajectory_title = "重複量測縱向軌跡圖",
    mean_title = "分組平均趨勢與 95% 信賴區間",
    combined_trajectory_title =
      "四項代謝指標之抽樣病患重複量測縱向軌跡圖",
    combined_mean_title =
      "四項代謝指標之分組平均趨勢與 95% 信賴區間",
    n_note = paste0(
      "n 為各組於該時間區間具有非缺失值之病患數；",
      "同一病患於同一區間有多筆紀錄時先取平均。"
    )
  ),
  en = list(
    suffix = "en",
    control = "Control group",
    statin = "Statin group",
    follow_year = "Follow-up time (years)",
    follow_period = "Follow-up interval",
    group_title = "Group",
    trajectory_title =
      "Longitudinal trajectories of repeated measurements",
    mean_title =
      "Group-specific mean trends with 95% confidence intervals",
    combined_trajectory_title =
      "Sampled-patient longitudinal trajectories for four metabolic indicators",
    combined_mean_title =
      "Mean trends with 95% confidence intervals for four metabolic indicators",
    n_note = paste0(
      "n denotes patients with an observed value in each group and interval; ",
      "multiple records from the same patient within an interval are averaged first."
    )
  )
)

recode_group_for_plot <- function(x, lang_cfg) {
  factor(
    dplyr::recode(
      as.character(x),
      Control = lang_cfg$control,
      Statin = lang_cfg$statin
    ),
    levels = c(
      lang_cfg$control,
      lang_cfg$statin
    )
  )
}

plot_mean_trend_ci <- function(
    stat_df,
    outcome_name,
    lang_cfg
) {
  outcome_label <- unname(
    outcome_labels[outcome_name]
  )

  if (is.na(outcome_label)) {
    outcome_label <- outcome_name
  }

  dat_p <- stat_df %>%
    dplyr::filter(
      outcome == outcome_name,
      Group %in% c(
        "Control",
        "Statin"
      ),
      !is.na(mean)
    ) %>%
    dplyr::mutate(
      Group_plot = recode_group_for_plot(
        Group,
        lang_cfg
      )
    )

  y_range <- range(
    c(
      dat_p$ci_lower,
      dat_p$ci_upper
    ),
    na.rm = TRUE
  )

  y_span <- diff(y_range)

  if (
    !is.finite(y_span) ||
      y_span <= 0
  ) {
    y_span <- max(
      abs(dat_p$mean),
      na.rm = TRUE
    ) * 0.05
  }

  if (
    !is.finite(y_span) ||
      y_span <= 0
  ) {
    y_span <- 1
  }

  dat_p <- dat_p %>%
    dplyr::mutate(
      n_label = dplyr::if_else(
        Time_Num %in% N_LABEL_TIMEPOINTS,
        paste0(
          "n=",
          patient_n
        ),
        NA_character_
      ),
      n_label_y = dplyr::if_else(
        as.character(Group) == "Control",
        ci_lower - 0.035 * y_span,
        ci_upper + 0.035 * y_span
      )
    )

  title_text <- if (
    lang_cfg$suffix == "zh"
  ) {
    paste0(
      outcome_label,
      " ",
      lang_cfg$mean_title
    )
  } else {
    paste0(
      outcome_label,
      ": ",
      lang_cfg$mean_title
    )
  }

  p <- ggplot2::ggplot(
    dat_p,
    ggplot2::aes(
      x = Time_Num,
      y = mean,
      color = Group_plot,
      group = Group_plot
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 1.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = ci_lower,
        ymax = ci_upper
      ),
      width = 0.15,
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = n_label_y,
        label = n_label
      ),
      size = 2.1,
      show.legend = FALSE,
      check_overlap = TRUE,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = 1:19,
      labels = paste0(
        "T",
        1:19
      )
    ) +
    ggplot2::labs(
      title = title_text,
      x = lang_cfg$follow_period,
      y = outcome_label,
      color = lang_cfg$group_title,
      caption = lang_cfg$n_note
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold"
      ),
      axis.title = ggplot2::element_text(
        size = 12
      ),
      legend.title = ggplot2::element_text(
        size = 11
      ),
      legend.text = ggplot2::element_text(
        size = 10
      ),
      plot.caption = ggplot2::element_text(
        size = 8,
        hjust = 0
      ),
      plot.margin = ggplot2::margin(
        8,
        8,
        8,
        8
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.10,
          0.10
        )
      )
    )

  p
}

plot_longitudinal_trajectory <- function(
    dat,
    outcome_name,
    lang_cfg
) {
  outcome_label <- unname(
    outcome_labels[outcome_name]
  )

  if (is.na(outcome_label)) {
    outcome_label <- outcome_name
  }

  dat_p <- dat %>%
    dplyr::filter(
      !is.na(.data[[outcome_name]]),
      !is.na(Time_in_Years),
      !is.na(Group),
      !is.na(ID)
    ) %>%
    dplyr::mutate(
      Group_plot = recode_group_for_plot(
        Group,
        lang_cfg
      )
    )

  title_text <- if (
    lang_cfg$suffix == "zh"
  ) {
    paste0(
      outcome_label,
      " ",
      lang_cfg$trajectory_title
    )
  } else {
    paste0(
      outcome_label,
      ": ",
      lang_cfg$trajectory_title
    )
  }

  p <- ggplot2::ggplot(
    dat_p,
    ggplot2::aes(
      x = Time_in_Years,
      y = .data[[outcome_name]],
      group = ID
    )
  ) +
    ggplot2::geom_line(
      alpha = 0.18,
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      alpha = 0.18,
      size = 0.45,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~ Group_plot
    ) +
    ggplot2::labs(
      title = title_text,
      x = lang_cfg$follow_year,
      y = outcome_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold"
      ),
      axis.title = ggplot2::element_text(
        size = 12
      ),
      strip.text = ggplot2::element_text(
        size = 12
      )
    )

  p
}

# 主要圖：分組平均趨勢與 95% CI。
mean_plots_by_lang <- list()

for (
  lang_name in names(
    plot_language_config
  )
) {
  lang_cfg <- plot_language_config[[lang_name]]

  mean_plots <- purrr::map(
    primary_outcomes,
    function(outcome_name) {
      plot_mean_trend_ci(
        stat_df = mean_sd_ci_figure4_1,
        outcome_name = outcome_name,
        lang_cfg = lang_cfg
      )
    }
  )

  names(mean_plots) <- primary_outcomes

  # 合併圖使用簡潔 panel title，並將 n 的說明集中放在整張圖下方，
  # 避免英文長標題被裁切及四個 panel caption 彼此重疊。
  mean_plots_compact <- purrr::imap(
    mean_plots[primary_outcomes],
    function(p, outcome_name) {
      p +
        ggplot2::labs(
          title = unname(outcome_labels[outcome_name]),
          caption = NULL
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 13,
            face = "bold",
            hjust = 0.5
          )
        )
    }
  )

  combined_mean <- patchwork::wrap_plots(
    mean_plots_compact,
    ncol = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = lang_cfg$combined_mean_title,
      caption = lang_cfg$n_note
    )

  combined_mean <- combined_mean &
    ggplot2::theme(
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(
        size = 8,
        hjust = 0
      )
    )

  ggplot2::ggsave(
    file.path(
      figure_folder,
      paste0(
        "combined_mean_trend_ci_primary_outcomes_",
        lang_cfg$suffix,
        ".png"
      )
    ),
    combined_mean,
    width = 14,
    height = 10,
    dpi = PLOT_DPI,
    bg = "white"
  )

  mean_plots_by_lang[[lang_name]] <- mean_plots
}

# 全體 ID 軌跡圖僅作 QC，預設不輸出。
if (isTRUE(make_all_trajectory)) {
  for (
    lang_name in names(
      plot_language_config
    )
  ) {
    lang_cfg <- plot_language_config[[lang_name]]

    all_plots <- purrr::map(
      primary_outcomes,
      function(outcome_name) {
        plot_longitudinal_trajectory(
          dat = dt_model,
          outcome_name = outcome_name,
          lang_cfg = lang_cfg
        )
      }
    )

    names(all_plots) <- primary_outcomes

    all_plots_compact <- purrr::imap(
      all_plots[primary_outcomes],
      function(p, outcome_name) {
        p +
          ggplot2::labs(
            title = unname(outcome_labels[outcome_name])
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = 13,
              face = "bold",
              hjust = 0.5
            )
          )
      }
    )

    combined_all <- patchwork::wrap_plots(
      all_plots_compact,
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = if (
          lang_cfg$suffix == "zh"
        ) {
          "四項代謝指標之全體病患重複量測縱向軌跡圖"
        } else {
          "Longitudinal trajectories of all patients for four metabolic indicators"
        }
      )

    ggplot2::ggsave(
      file.path(
        figure_folder,
        paste0(
          "combined_trajectory_all_primary_outcomes_",
          lang_cfg$suffix,
          ".png"
        )
      ),
      combined_all,
      width = 13,
      height = 9.5,
      dpi = PLOT_DPI,
      bg = "white"
    )
  }

  all_trajectory_summary <- tibble::tibble(
    status = "exported",
    note = "全體病患軌跡圖僅供 QC，不建議放入正文。"
  )
} else {
  # 避免舊版本留下的全體軌跡圖被誤認為本次輸出。
  stale_all_trajectory_files <- list.files(
    figure_folder,
    pattern = paste0(
      "^(trajectory_all_|",
      "combined_trajectory_all_).*[.]png$"
    ),
    full.names = TRUE
  )

  if (
    length(
      stale_all_trajectory_files
    ) > 0L
  ) {
    unlink(
      stale_all_trajectory_files
    )
  }

  all_trajectory_summary <- tibble::tibble(
    status = "not exported",
    removed_stale_file_n = length(
      stale_all_trajectory_files
    ),
    note = paste0(
      "全體病患軌跡圖因嚴重 overplotting 預設不輸出；",
      "正文使用分組平均趨勢、95% CI 與各區間 n。"
    )
  )
}

# 附錄／QC：各組固定隨機抽樣病患軌跡圖。
if (isTRUE(make_sampled_trajectory)) {
  set.seed(RANDOM_SEED)

  sampled_id_df <- dt_model %>%
    dplyr::distinct(
      ID,
      Group
    ) %>%
    dplyr::group_by(Group) %>%
    dplyr::group_modify(
      ~ dplyr::slice_sample(
        .x,
        n = min(
          sample_per_group_trajectory,
          nrow(.x)
        )
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      Group,
      ID
    )

  sampled_ids <- sampled_id_df$ID

  sampled_trajectory_id_summary <- sampled_id_df %>%
    dplyr::count(
      Group,
      name = "sampled_patient_n"
    )

  readr::write_csv(
    sampled_id_df,
    file.path(
      output_root,
      "sampled_ids_for_trajectory_plots.csv"
    )
  )

  for (
    lang_name in names(
      plot_language_config
    )
  ) {
    lang_cfg <- plot_language_config[[lang_name]]

    sampled_plots <- purrr::map(
      primary_outcomes,
      function(outcome_name) {
        plot_longitudinal_trajectory(
          dat = dt_model %>%
            dplyr::filter(
              ID %in% sampled_ids
            ),
          outcome_name = outcome_name,
          lang_cfg = lang_cfg
        )
      }
    )

    names(sampled_plots) <- primary_outcomes

    sampled_plots_compact <- purrr::imap(
      sampled_plots[primary_outcomes],
      function(p, outcome_name) {
        p +
          ggplot2::labs(
            title = unname(outcome_labels[outcome_name])
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = 13,
              face = "bold",
              hjust = 0.5
            )
          )
      }
    )

    combined_sampled <- patchwork::wrap_plots(
      sampled_plots_compact,
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = lang_cfg$combined_trajectory_title
      )

    ggplot2::ggsave(
      file.path(
        figure_folder,
        paste0(
          "combined_trajectory_sampled_primary_outcomes_",
          lang_cfg$suffix,
          ".png"
        )
      ),
      combined_sampled,
      width = 13,
      height = 9.5,
      dpi = PLOT_DPI,
      bg = "white"
    )
  }
} else {
  sampled_trajectory_id_summary <- tibble::tibble(
    note = "未輸出抽樣病患軌跡圖。"
  )
}

# ============================================================
# 10. 病患層級 train／validation／test 分割
# ============================================================

set.seed(RANDOM_SEED)

id_df <- dt_model %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    Group = dplyr::first(Group),
    n_visit = dplyr::n(),
    distinct_report_date_n = dplyr::n_distinct(
      報告日期
    ),
    .groups = "drop"
  )

split_manifest <- id_df %>%
  dplyr::group_by(Group) %>%
  dplyr::group_modify(
    ~ split_ids_by_group(
      .x,
      train_prop = 0.70,
      validation_prop = 0.15
    )
  ) %>%
  dplyr::ungroup()

dt_dl_raw_full <- dt_model %>%
  dplyr::left_join(
    split_manifest %>%
      dplyr::select(
        ID,
        split
      ),
    by = "ID",
    relationship = "many-to-one"
  )

split_check <- dt_dl_raw_full %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    n_split = dplyr::n_distinct(split),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_split != 1L
  )

split_count_patient <- split_manifest %>%
  dplyr::count(
    split,
    Group,
    name = "patient_n"
  )

split_count_observation <- dt_dl_raw_full %>%
  dplyr::count(
    split,
    Group,
    name = "observation_n"
  )

split_outcome_missing <- dt_dl_raw_full %>%
  dplyr::group_by(
    split,
    Group
  ) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    missing_summary = purrr::map(
      data,
      ~ make_missing_summary(
        .x,
        c(
          primary_outcomes,
          "BMI"
        )
      )
    )
  ) %>%
  dplyr::select(
    -data
  ) %>%
  tidyr::unnest(
    missing_summary
  ) %>%
  dplyr::ungroup()

# DL 標準輸出只保留建模需要欄位，避免把 89 個 QA／追溯欄位
# 與精確測量時間誤當成模型特徵。
dl_model_vars <- c(
  "row_id",
  "ID",
  "Group",
  "split",
  "報告日期",
  "Time_in_Days",
  "Time_in_Years",
  "Time_Category",
  "HbA1c",
  "Glu_AC",
  "LDL",
  "TG",
  "BMI",
  "Age_baseline",
  "Sex",
  "CVD_History",
  "HTN_History",
  "CKD_History",
  "Liver_History",
  "HbA1c_observed",
  "Glu_AC_observed",
  "LDL_observed",
  "TG_observed",
  "BMI_observed",
  "source_row_n",
  "source_visit_n",
  "source_visit_ids"
)

dt_dl_raw <- dt_dl_raw_full %>%
  dplyr::select(
    dplyr::all_of(dl_model_vars)
  )

readr::write_csv(
  split_manifest,
  file.path(
    output_root,
    "dl_id_split_manifest_type2.csv"
  )
)

# `dt_dl_raw` 與 `dt_dl_raw_full` 僅保留於本次 R session 作 split／QC 使用。
# 最終 DL 與 GEE/LME 預測統一讀取 prediction_long_with_split_type2.csv，
# 因此不再輸出重複的 dl_raw_* CSV。

# ============================================================
# 10.1 缺失機制與 MICE 適用性檢查
# ============================================================

# ------------------------------------------------------------
# 0. 套件與基本設定
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "geepack",
  "writexl"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "尚未安裝下列套件：",
      paste(missing_packages, collapse = ", ")
    )
  )
}

# 用於完整缺失情形檢查。
# HDL 與 TC 僅作為排除依據之資料品質確認，
# 不列入本研究後續 MICE、GEE、LME 或預測模型。
missingness_review_vars <- c(
  "Glu_AC",
  "LDL",
  "TG",
  "BMI",
  "HDL",
  "TC"
)

# 實際需要評估 MICE 適用性的候選變數。
mice_candidate_vars <- c(
  "Glu_AC",
  "LDL",
  "TG",
  "BMI"
)

# 確認必要欄位存在。
required_columns <- unique(
  c(
    "ID",
    "Group",
    "Time_in_Years",
    "HbA1c",
    "Sex",
    "Age_baseline",
    "CVD_History",
    "HTN_History",
    "CKD_History",
    "Liver_History",
    missingness_review_vars
  )
)

missing_columns <- setdiff(
  required_columns,
  names(dt_model)
)

if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "dt_model 缺少下列欄位：",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# 確認每位病患只有一個固定組別。
group_check <- dt_model %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    n_group = dplyr::n_distinct(Group),
    .groups = "drop"
  )

if (any(group_check$n_group > 1L)) {
  stop(
    "部分病患同時出現在多個 Group，請先確認 Statin／Control 分組。"
  )
}

# 因子層級設定。
# 若 Group 原本就是 Control、Statin，此設定可確保 Control 為參考組。
missing_test_data <- dt_model %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c("Control", "Statin")
    ),
    Sex = factor(Sex),
    miss_Glu_AC = as.integer(is.na(Glu_AC)),
    miss_LDL = as.integer(is.na(LDL)),
    miss_TG = as.integer(is.na(TG)),
    miss_BMI = as.integer(is.na(BMI))
  )

# ------------------------------------------------------------
# 1. 事件層級整體缺失率
# ------------------------------------------------------------

missing_rate_overall <- purrr::map_dfr(
  missingness_review_vars,
  function(v) {
    tibble::tibble(
      variable = v,
      total_observation_n = nrow(dt_model),
      observed_n = sum(
        !is.na(dt_model[[v]])
      ),
      missing_n = sum(
        is.na(dt_model[[v]])
      ),
      missing_rate_percent = round(
        mean(is.na(dt_model[[v]])) * 100,
        2
      )
    )
  }
)

# ------------------------------------------------------------
# 2. 事件層級分組缺失率
# ------------------------------------------------------------

missing_rate_by_group <- purrr::map_dfr(
  missingness_review_vars,
  function(v) {
    dt_model %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        variable = v,
        total_observation_n = dplyr::n(),
        observed_n = sum(
          !is.na(.data[[v]])
        ),
        missing_n = sum(
          is.na(.data[[v]])
        ),
        missing_rate_percent = round(
          mean(is.na(.data[[v]])) * 100,
          2
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        variable,
        Group,
        dplyr::everything()
      )
  }
)

# ------------------------------------------------------------
# 3. 病患層級：是否曾經量測
# ------------------------------------------------------------

never_measured_summary <- purrr::map_dfr(
  missingness_review_vars,
  function(v) {
    patient_level <- dt_model %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(
        ever_observed = any(
          !is.na(.data[[v]])
        ),
        .groups = "drop"
      )
    
    tibble::tibble(
      variable = v,
      total_patient_n = nrow(patient_level),
      ever_observed_patient_n = sum(
        patient_level$ever_observed
      ),
      never_observed_patient_n = sum(
        !patient_level$ever_observed
      ),
      never_observed_rate_percent = round(
        mean(!patient_level$ever_observed) * 100,
        2
      )
    )
  }
)

# ------------------------------------------------------------
# 4. 病患層級分組：是否曾經量測
# ------------------------------------------------------------

never_measured_by_group <- purrr::map_dfr(
  missingness_review_vars,
  function(v) {
    dt_model %>%
      dplyr::group_by(ID, Group) %>%
      dplyr::summarise(
        ever_observed = any(
          !is.na(.data[[v]])
        ),
        .groups = "drop"
      ) %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        variable = v,
        total_patient_n = dplyr::n(),
        ever_observed_patient_n = sum(
          ever_observed
        ),
        never_observed_patient_n = sum(
          !ever_observed
        ),
        never_observed_rate_percent = round(
          mean(!ever_observed) * 100,
          2
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        variable,
        Group,
        dplyr::everything()
      )
  }
)

# ------------------------------------------------------------
# 5. 缺失指標規格
# ------------------------------------------------------------

missingness_specs <- tibble::tribble(
  ~target,  ~missing_indicator,
  "Glu_AC", "miss_Glu_AC",
  "LDL",    "miss_LDL",
  "TG",     "miss_TG",
  "BMI",    "miss_BMI"
)

# 缺失模型使用完整度較高，且於缺失發生時仍可觀測之變數。
missingness_predictors <- c(
  "Group",
  "Time_in_Years",
  "HbA1c",
  "Sex",
  "Age_baseline",
  "CVD_History",
  "HTN_History",
  "CKD_History",
  "Liver_History"
)

# ------------------------------------------------------------
# 6. 一般 logistic regression
#    僅作為初步檢查
# ------------------------------------------------------------

fit_missingness_glm <- function(
    data,
    target,
    missing_indicator,
    predictors
) {
  
  predictors <- intersect(
    predictors,
    names(data)
  )
  
  required_vars <- c(
    missing_indicator,
    predictors
  )
  
  model_data <- data %>%
    dplyr::select(
      dplyr::all_of(required_vars)
    ) %>%
    tidyr::drop_na()
  
  if (
    nrow(model_data) == 0L ||
    dplyr::n_distinct(
      model_data[[missing_indicator]]
    ) < 2L
  ) {
    return(
      list(
        model = NULL,
        coefficient_result = tibble::tibble(
          target = target,
          analysis_n = nrow(model_data),
          note = paste0(
            "Insufficient data or no variation ",
            "in missingness indicator."
          )
        )
      )
    )
  }
  
  form <- stats::reformulate(
    termlabels = predictors,
    response = missing_indicator
  )
  
  fit <- stats::glm(
    formula = form,
    data = model_data,
    family = stats::binomial(
      link = "logit"
    ),
    model = TRUE,
    x = FALSE,
    y = TRUE
  )
  
  coef_table <- summary(fit)$coefficients
  
  coefficient_result <- tibble::as_tibble(
    coef_table,
    rownames = "term"
  ) %>%
    dplyr::rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      z_value = `z value`,
      p_value = `Pr(>|z|)`
    ) %>%
    dplyr::mutate(
      target = target,
      analysis_n = nrow(model_data),
      odds_ratio = exp(estimate),
      ci_lower = exp(
        estimate - 1.96 * std_error
      ),
      ci_upper = exp(
        estimate + 1.96 * std_error
      ),
      note = "",
      .before = 1
    )
  
  list(
    model = fit,
    coefficient_result = coefficient_result
  )
}

missingness_glm_models <- purrr::pmap(
  missingness_specs,
  function(target, missing_indicator) {
    fit_missingness_glm(
      data = missing_test_data,
      target = target,
      missing_indicator = missing_indicator,
      predictors = missingness_predictors
    )
  }
)

names(missingness_glm_models) <- missingness_specs$target

missingness_glm_results <- purrr::map_dfr(missingness_glm_models, "coefficient_result")

# ------------------------------------------------------------
# 7. GLM 整體 likelihood-ratio test
# ------------------------------------------------------------

missingness_glm_global_tests <-
  purrr::imap_dfr(
    missingness_glm_models,
    function(obj, target_name) {
      
      if (is.null(obj$model)) {
        return(
          tibble::tibble(
            target = target_name,
            analysis_n = NA_integer_,
            likelihood_ratio_statistic = NA_real_,
            degrees_of_freedom = NA_integer_,
            likelihood_ratio_p = NA_real_,
            note = "Model unavailable."
          )
        )
      }
      
      model_frame <- stats::model.frame(
        obj$model
      )
      
      response_vector <- stats::model.response(
        model_frame
      )
      
      null_fit <- stats::glm(
        response_vector ~ 1,
        family = stats::binomial(
          link = "logit"
        )
      )
      
      loglik_full <- stats::logLik(
        obj$model
      )
      
      loglik_null <- stats::logLik(
        null_fit
      )
      
      lr_statistic <- 2 * (
        as.numeric(loglik_full) -
          as.numeric(loglik_null)
      )
      
      df_difference <-
        attr(loglik_full, "df") -
        attr(loglik_null, "df")
      
      lr_p_value <- stats::pchisq(
        q = lr_statistic,
        df = df_difference,
        lower.tail = FALSE
      )
      
      tibble::tibble(
        target = target_name,
        analysis_n = stats::nobs(
          obj$model
        ),
        likelihood_ratio_statistic =
          lr_statistic,
        degrees_of_freedom =
          as.integer(df_difference),
        likelihood_ratio_p =
          lr_p_value,
        note = ""
      )
    }
  )

# ------------------------------------------------------------
# 8. GLM 模型診斷
# ------------------------------------------------------------

missingness_glm_diagnostics <-
  purrr::imap_dfr(
    missingness_glm_models,
    function(obj, target_name) {
      
      if (is.null(obj$model)) {
        return(
          tibble::tibble(
            target = target_name,
            converged = FALSE,
            analysis_n = NA_integer_,
            min_fitted = NA_real_,
            max_fitted = NA_real_,
            max_abs_coefficient = NA_real_,
            note = "Model unavailable."
          )
        )
      }
      
      fitted_values <- stats::fitted(
        obj$model
      )
      
      coefficient_values <- stats::coef(
        obj$model
      )
      
      tibble::tibble(
        target = target_name,
        converged = isTRUE(
          obj$model$converged
        ),
        analysis_n = stats::nobs(
          obj$model
        ),
        min_fitted = min(
          fitted_values,
          na.rm = TRUE
        ),
        max_fitted = max(
          fitted_values,
          na.rm = TRUE
        ),
        max_abs_coefficient = max(
          abs(coefficient_values),
          na.rm = TRUE
        ),
        note = dplyr::case_when(
          !isTRUE(obj$model$converged) ~
            "Model did not converge.",
          max(
            abs(coefficient_values),
            na.rm = TRUE
          ) > 10 ~
            paste0(
              "Large coefficient detected; ",
              "check possible separation."
            ),
          min(
            fitted_values,
            na.rm = TRUE
          ) < 1e-6 ||
            max(
              fitted_values,
              na.rm = TRUE
            ) > 1 - 1e-6 ~
            paste0(
              "Extreme fitted probability detected; ",
              "check possible separation."
            ),
          TRUE ~ ""
        )
      )
    }
  )

# ------------------------------------------------------------
# 9. Cluster-robust logistic regression
#    以病患 ID 為群集計算穩健標準誤
# ------------------------------------------------------------

fit_missingness_cluster_glm <- function(
    data,
    target,
    missing_indicator,
    predictors
) {
  
  predictors <- intersect(predictors, names(data))
  
  required_vars <- c(
    "ID",
    missing_indicator,
    predictors
  )
  
  model_data <- data %>%
    dplyr::select(
      dplyr::all_of(required_vars)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::arrange(
      ID,
      Time_in_Years
    )
  
  if (
    nrow(model_data) == 0L ||
    dplyr::n_distinct(
      model_data[[missing_indicator]]
    ) < 2L
  ) {
    return(
      list(
        model = NULL,
        robust_vcov = NULL,
        analysis_n = nrow(model_data),
        patient_n = dplyr::n_distinct(model_data$ID),
        coefficient_result = tibble::tibble(
          target = target,
          analysis_n = nrow(model_data),
          patient_n = dplyr::n_distinct(model_data$ID),
          term = NA_character_,
          estimate = NA_real_,
          robust_std_error = NA_real_,
          test_statistic = NA_real_,
          p_value = NA_real_,
          odds_ratio = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          note = paste0(
            "Insufficient data or no variation ",
            "in missingness indicator."
          )
        )
      )
    )
  }
  
  form <- stats::reformulate(
    termlabels = predictors,
    response = missing_indicator
  )
  
  fit <- stats::glm(
    formula = form,
    data = model_data,
    family = stats::binomial(
      link = "logit"
    ),
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  cluster_vcov <- sandwich::vcovCL(
    fit,
    cluster = model_data$ID,
    type = "HC0"
  )
  
  robust_test <- lmtest::coeftest(
    fit,
    vcov. = cluster_vcov
  )
  
  robust_matrix <- as.matrix(
    robust_test
  )
  
  if (ncol(robust_matrix) < 4L) {
    stop(
      paste0(
        "coeftest 輸出欄數不足。實際欄位為：",
        paste(
          colnames(robust_matrix),
          collapse = ", "
        )
      )
    )
  }
  
  coefficient_result <- tibble::tibble(
    target = target,
    analysis_n = nrow(model_data),
    patient_n = dplyr::n_distinct(
      model_data$ID
    ),
    term = rownames(robust_matrix),
    estimate = as.numeric(
      robust_matrix[, 1]
    ),
    robust_std_error = as.numeric(
      robust_matrix[, 2]
    ),
    test_statistic = as.numeric(
      robust_matrix[, 3]
    ),
    p_value = as.numeric(
      robust_matrix[, 4]
    )
  ) %>%
    dplyr::mutate(
      odds_ratio = exp(estimate),
      ci_lower = exp(
        estimate - 1.96 * robust_std_error
      ),
      ci_upper = exp(
        estimate + 1.96 * robust_std_error
      ),
      note = ""
    )
  
  list(
    model = fit,
    robust_vcov = cluster_vcov,
    analysis_n = nrow(model_data),
    patient_n = dplyr::n_distinct(model_data$ID),
    coefficient_result = coefficient_result
  )
}

missingness_cluster_models <- purrr::pmap(
  missingness_specs,
  function(target, missing_indicator) {
    
    message(
      format(Sys.time(), "%H:%M:%S"),
      " 開始：",
      target
    )
    
    result <- fit_missingness_cluster_glm(
      data = missing_test_data,
      target = target,
      missing_indicator = missing_indicator,
      predictors = missingness_predictors
    )
    
    message(
      format(Sys.time(), "%H:%M:%S"),
      " 完成：",
      target
    )
    
    result
  }
)

names(missingness_cluster_models) <- missingness_specs$target
missingness_cluster_results <- purrr::map_dfr(missingness_cluster_models, "coefficient_result")
print(missingness_cluster_results, n = Inf)

# ------------------------------------------------------------
# 10. Cluster-robust logistic regression 模型診斷摘要
# ------------------------------------------------------------

missingness_cluster_diagnostics <- purrr::imap_dfr(missingness_cluster_models, function(obj, target_name) {
      
      if (is.null(obj$model)) {
        return(
          tibble::tibble(
            target = target_name,
            converged = FALSE,
            analysis_n = NA_integer_,
            patient_n = NA_integer_,
            finite_coefficients = FALSE,
            finite_robust_variance = FALSE,
            min_fitted = NA_real_,
            max_fitted = NA_real_,
            max_abs_coefficient = NA_real_,
            note = "Model unavailable."
          )
        )
      }
      
      coefficient_values <- stats::coef(obj$model)
      fitted_values <- stats::fitted(obj$model)
      robust_variance_values <- as.numeric(obj$robust_vcov)
      finite_coefficients <- all(is.finite(coefficient_values))
      finite_robust_variance <- all(is.finite(robust_variance_values))
      max_abs_coefficient <- max(abs(coefficient_values), na.rm = TRUE)
      min_fitted <- min(fitted_values, na.rm = TRUE)
      max_fitted <- max(fitted_values, na.rm = TRUE)
      
      tibble::tibble(
        target = target_name,
        converged = isTRUE(
          obj$model$converged
        ),
        analysis_n = stats::nobs(
          obj$model
        ),
        patient_n = obj$patient_n,
        finite_coefficients =
          finite_coefficients,
        finite_robust_variance =
          finite_robust_variance,
        min_fitted =
          min_fitted,
        max_fitted =
          max_fitted,
        max_abs_coefficient =
          max_abs_coefficient,
        note = dplyr::case_when(
          !isTRUE(obj$model$converged) ~
            "Model did not converge.",
          !finite_coefficients ~
            "Non-finite coefficient detected.",
          !finite_robust_variance ~
            "Non-finite cluster-robust variance detected.",
          max_abs_coefficient > 10 ~
            paste0(
              "Large coefficient detected; ",
              "check possible separation."
            ),
          min_fitted < 1e-6 |
            max_fitted > 1 - 1e-6 ~
            paste0(
              "Extreme fitted probability detected; ",
              "check possible separation."
            ),
          TRUE ~ ""
        )
      )
    }
  )

# ------------------------------------------------------------
# 11. 依老師提出之 15% 原則建立初步判斷
# ------------------------------------------------------------

mice_eligibility_summary <- missing_rate_overall %>%
  dplyr::filter(
    variable %in% mice_candidate_vars
  ) %>%
  dplyr::left_join(
    never_measured_summary,
    by = "variable"
  ) %>%
  dplyr::left_join(
    missingness_glm_global_tests %>%
      dplyr::rename(
        glm_global_lr_statistic =
          likelihood_ratio_statistic,
        glm_global_df =
          degrees_of_freedom,
        glm_global_lr_p =
          likelihood_ratio_p
      ),
    by = c(
      "variable" = "target"
    )
  ) %>%
  dplyr::mutate(
    missing_rate_over_15 =
      missing_rate_percent > 15,
    never_observed_rate_over_15 =
      never_observed_rate_percent > 15,
    glm_missingness_association =
      !is.na(glm_global_lr_p) &
      glm_global_lr_p < 0.05,
    preliminary_note =
      dplyr::case_when(
        missing_rate_over_15 &
          never_observed_rate_over_15 ~
          paste0(
            "Overall missing rate exceeds 15%, ",
            "and more than 15% of patients were never observed. ",
            "Comprehensive imputation requires reconsideration."
          ),
        missing_rate_over_15 ~
          paste0(
            "Overall missing rate exceeds 15%; ",
            "review with advisor before imputation."
          ),
        glm_missingness_association ~
          paste0(
            "The conventional logistic model indicates ",
            "an association between missingness and observed variables. ",
            "Cluster-robust coefficient results should be used ",
            "for individual predictor interpretation."
          ),
        TRUE ~
          paste0(
            "Overall missing rate does not exceed 15%, ",
            "but imputation suitability still requires ",
            "clinical and statistical review."
          )
      )
  )

# ------------------------------------------------------------
# 12. HDL 與 TC 排除依據摘要
# ------------------------------------------------------------

excluded_variable_summary <-
  missing_rate_overall %>%
  dplyr::filter(
    variable %in% c(
      "HDL",
      "TC"
    )
  ) %>%
  dplyr::left_join(
    never_measured_summary,
    by = "variable"
  ) %>%
  dplyr::mutate(
    analysis_status =
      paste0(
        "Excluded from MICE and subsequent main analyses; ",
        "retained only for missingness quality review."
      )
  )

# ------------------------------------------------------------
# 13. 統一輸出 Excel
# ------------------------------------------------------------

output_file <- file.path(qc_folder, "missingness_mechanism_and_mice_review.xlsx")

writexl::write_xlsx(
  list(
    missing_rate_overall =
      missing_rate_overall,
    missing_rate_by_group =
      missing_rate_by_group,
    never_measured_overall =
      never_measured_summary,
    never_measured_by_group =
      never_measured_by_group,
    glm_coefficients =
      missingness_glm_results,
    glm_global_tests =
      missingness_glm_global_tests,
    glm_diagnostics =
      missingness_glm_diagnostics,
    cluster_coefficients =
      missingness_cluster_results,
    cluster_diagnostics =
      missingness_cluster_diagnostics,
    mice_eligibility =
      mice_eligibility_summary,
    excluded_HDL_TC =
      excluded_variable_summary
  ),
  output_file
)

message("缺失機制與 MICE 適用性檢查已完成：")
message(output_file)

# ============================================================
# 11. MICE：完整 cohort 推論用多重插補
# ============================================================

# 目前論文設定：
# - 使用 PMM 插補 Glu-AC、LDL、TG 與時變 BMI。
# - HbA1c 無缺失，只作預測因子。
# - Age、Sex、Group、時間、BP、HDL、TC 與共病不插補。
# - 插補在原始尺度執行；本程式不先取 log，也不先標準化。
# - GEE / LME 的 log 轉換與 DL 的 log + training-set z-score，均在後續模型程式執行。
# - 此為完整 cohort 推論資料；不將其視為 leakage-safe 的 DL 前處理。

mice_core_vars <- c(
  "row_id",
  "Group",
  "Time_in_Years",
  "HbA1c",
  "Glu_AC",
  "LDL",
  "TG",
  "BMI",
  "Sex",
  "Age_baseline",
  "CVD_History",
  "HTN_History",
  "CKD_History",
  "Liver_History"
)

dt_mice <- dt_model %>%
  dplyr::select(
    dplyr::all_of(mice_core_vars)
  ) %>%
  dplyr::mutate(
    row_id = as.integer(row_id),
    Group = factor(
      Group,
      levels = c(
        "Control",
        "Statin"
      )
    ),
    Sex = factor(
      Sex,
      levels = c(
        "Female",
        "Male"
      )
    ),
    dplyr::across(
      c(
        CVD_History,
        HTN_History,
        CKD_History,
        Liver_History
      ),
      factor
    )
  )

if (
  any(is.na(dt_mice$Age_baseline)) ||
    any(is.na(dt_mice$Sex)) ||
    any(is.na(dt_mice$Group))
) {
  stop(
    "MICE 前 Age_baseline、Sex 或 Group 仍有缺失。"
  )
}

ini <- mice::mice(
  dt_mice,
  maxit = 0,
  printFlag = FALSE
)

pred_matrix <- ini$predictorMatrix
meth_vec <- ini$method

pred_matrix[,] <- 0
meth_vec[] <- ""

mice_targets <- intersect(
  c(
    "Glu_AC",
    "LDL",
    "TG",
    "BMI"
  ),
  names(dt_mice)
)

mice_predictors <- intersect(
  c(
    "Group",
    "Time_in_Years",
    "HbA1c",
    "Glu_AC",
    "LDL",
    "TG",
    "BMI",
    "Sex",
    "Age_baseline",
    "CVD_History",
    "HTN_History",
    "CKD_History",
    "Liver_History"
  ),
  names(dt_mice)
)

for (target in mice_targets) {
  if (any(is.na(dt_mice[[target]]))) {
    meth_vec[target] <- "pmm"

    target_predictors <- setdiff(
      mice_predictors,
      target
    )

    pred_matrix[
      target,
      target_predictors
    ] <- 1
  }
}

meth_vec["row_id"] <- ""
pred_matrix["row_id", ] <- 0
pred_matrix[, "row_id"] <- 0

mice_setup_method <- tibble::tibble(
  variable = names(meth_vec),
  method = as.character(meth_vec),
  missing_n = purrr::map_int(
    names(meth_vec),
    ~ sum(
      is.na(
        dt_mice[[.x]]
      )
    )
  )
)

mice_predictor_matrix_export <- as.data.frame(
  pred_matrix
) %>%
  tibble::rownames_to_column(
    "target_variable"
  )

set.seed(RANDOM_SEED)

imputed_data_mids <- mice::mice(
  dt_mice,
  m = MICE_M,
  maxit = MICE_MAXIT,
  predictorMatrix = pred_matrix,
  method = meth_vec,
  seed = RANDOM_SEED,
  printFlag = TRUE
)

saveRDS(
  imputed_data_mids,
  file = file.path(
    rds_folder,
    "imputed_data_mids_type2_inference.rds"
  )
)

mice_metadata_cols <- setdiff(
  names(dt_model),
  names(dt_mice)
)

mice_metadata <- dt_model %>%
  dplyr::select(
    row_id,
    dplyr::all_of(
      mice_metadata_cols
    )
  )

imp_list <- lapply(
  seq_len(
    imputed_data_mids$m
  ),
  function(i) {
    completed_core <- mice::complete(
      imputed_data_mids,
      action = i
    )

    completed_core %>%
      dplyr::left_join(
        mice_metadata,
        by = "row_id",
        relationship = "one-to-one"
      ) %>%
      dplyr::arrange(row_id) %>%
      dplyr::mutate(
        .imp = i,
        ID = factor(ID),
        Group = factor(
          Group,
          levels = c(
            "Control",
            "Statin"
          )
        ),
        Sex = factor(
          Sex,
          levels = c(
            "Female",
            "Male"
          )
        ),
        Time_Category = factor(
          Time_Category,
          levels = paste0(
            "T",
            1:19
          ),
          ordered = TRUE
        )
      )
  }
)

saveRDS(
  imp_list,
  file = file.path(
    rds_folder,
    "imp_list_type2_inference.rds"
  )
)

mice_logged_events <- imputed_data_mids$loggedEvents

if (is.null(mice_logged_events)) {
  mice_logged_events <- tibble::tibble(
    note = "No logged events."
  )
}

# ------------------------------------------------------------
# 10.1 插補完整性與觀測值保留檢查
# ------------------------------------------------------------

imputation_completeness_qc <- purrr::map_dfr(
  seq_along(imp_list),
  function(i) {
    purrr::map_dfr(
      mice_targets,
      function(v) {
        tibble::tibble(
          .imp = i,
          variable = v,
          remaining_missing_n = sum(
            is.na(
              imp_list[[i]][[v]]
            )
          )
        )
      }
    )
  }
)

imputation_observed_preservation_qc <- purrr::map_dfr(
  seq_along(imp_list),
  function(i) {
    purrr::map_dfr(
      mice_targets,
      function(v) {
        original_value <- dt_model[[v]]
        completed_value <- imp_list[[i]][[v]]
        observed_index <- !is.na(original_value)

        changed_n <- sum(
          observed_index &
            (
              is.na(completed_value) |
                abs(
                  completed_value -
                    original_value
                ) > 1e-12
            ),
          na.rm = TRUE
        )

        tibble::tibble(
          .imp = i,
          variable = v,
          observed_value_changed_n = changed_n
        )
      }
    )
  }
)

# ------------------------------------------------------------
# 10.2 MICE 診斷圖
# ------------------------------------------------------------

mice_variable_labels <- c(
  Glu_AC = "Glu-AC",
  LDL = "LDL",
  TG = "TG",
  BMI = "BMI"
)

extract_mice_trace_df <- function(
    imp_object,
    vars
) {
  if (
    is.null(imp_object$chainMean) ||
      is.null(imp_object$chainVar)
  ) {
    return(tibble::tibble())
  }

  vars_available <- dimnames(
    imp_object$chainMean
  )[[1]]

  target_vars <- intersect(
    vars,
    vars_available
  )

  if (length(target_vars) == 0L) {
    return(tibble::tibble())
  }

  make_trace <- function(
      arr,
      stat_name,
      transform_fun
  ) {
    vars_arr <- dimnames(arr)[[1]]
    iterations <- seq_len(
      dim(arr)[2]
    )
    chains <- seq_len(
      dim(arr)[3]
    )

    df <- expand.grid(
      variable = vars_arr,
      iteration = iterations,
      chain = chains,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    df$value <- as.vector(arr)

    df %>%
      dplyr::filter(
        variable %in% target_vars
      ) %>%
      dplyr::mutate(
        value = transform_fun(value),
        stat = stat_name
      )
  }

  dplyr::bind_rows(
    make_trace(
      imp_object$chainMean,
      "mean",
      identity
    ),
    make_trace(
      imp_object$chainVar,
      "sd",
      function(x) {
        sqrt(
          pmax(
            x,
            0
          )
        )
      }
    )
  ) %>%
    dplyr::filter(
      !is.na(value),
      is.finite(value)
    ) %>%
    dplyr::mutate(
      variable_label = factor(
        mice_variable_labels[variable],
        levels = unname(
          mice_variable_labels[target_vars]
        )
      ),
      chain = factor(chain)
    )
}

trace_df <- extract_mice_trace_df(
  imputed_data_mids,
  mice_targets
)

make_trace_plot <- function(
    trace_data,
    lang
) {
  if (nrow(trace_data) == 0L) {
    return(NULL)
  }
  
  if (lang == "zh") {
    title_text <- "多重插補變數之收斂診斷圖"
    
    caption_text <- paste0(
      "圖中呈現 Glu-AC、LDL、TG 與 BMI 的 10 組插補鏈；",
      "HbA1c 因無缺失值，未進行插補。"
    )
    
    x_text <- "迭代次數"
    mean_text <- "平均值"
    sd_text <- "標準差"
    suffix <- "zh"
    
  } else if (lang == "en") {
    title_text <-
      "Convergence diagnostics for multiply imputed variables"
    
    caption_text <- paste0(
      "The 10 imputation chains are shown for Glu-AC, LDL, TG, and BMI. ",
      "HbA1c was not imputed because it had no missing values."
    )
    
    x_text <- "Iteration"
    mean_text <- "Mean"
    sd_text <- "Standard deviation"
    suffix <- "en"
    
  } else {
    stop("lang 僅可設定為 'zh' 或 'en'。")
  }
  
  panel_fun <- function(
    dat_panel,
    panel_title
  ) {
    ggplot2::ggplot(
      dat_panel,
      ggplot2::aes(
        x = iteration,
        y = value,
        group = chain,
        color = chain
      )
    ) +
      ggplot2::geom_line(
        linewidth = 0.6,
        alpha = 0.9,
        na.rm = TRUE
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(variable_label),
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = panel_title,
        x = x_text,
        y = NULL
      ) +
      ggplot2::theme_bw(
        base_size = 10
      ) +
      ggplot2::theme(
        legend.position = "none",
        
        plot.title = ggplot2::element_text(
          size = 12,
          face = "bold",
          hjust = 0.5
        ),
        
        strip.text = ggplot2::element_text(
          size = 11,
          face = "bold"
        ),
        
        axis.title.x = ggplot2::element_text(
          size = 10
        ),
        
        axis.text = ggplot2::element_text(
          size = 9
        ),
        
        panel.spacing = grid::unit(
          0.12,
          "cm"
        ),
        
        plot.margin = ggplot2::margin(
          t = 4,
          r = 4,
          b = 4,
          l = 4
        )
      )
  }
  
  mean_panel <- panel_fun(
    trace_data %>%
      dplyr::filter(stat == "mean"),
    mean_text
  )
  
  sd_panel <- panel_fun(
    trace_data %>%
      dplyr::filter(stat == "sd"),
    sd_text
  )
  
  panel_grid <- patchwork::wrap_plots(
    mean_panel,
    sd_panel,
    ncol = 2
  )
  
  title_grob <- grid::textGrob(
    label = title_text,
    x = 0.5,
    y = 0.5,
    hjust = 0.5,
    vjust = 0.5,
    gp = grid::gpar(
      fontsize = 15,
      fontface = "bold"
    )
  )
  
  caption_grob <- grid::textGrob(
    label = caption_text,
    x = 0,
    y = 0.5,
    hjust = 0,
    vjust = 0.5,
    gp = grid::gpar(
      fontsize = 8
    )
  )
  
  title_panel <- patchwork::wrap_elements(
    full = title_grob
  )
  
  caption_panel <- patchwork::wrap_elements(
    full = caption_grob
  )
  
  p <- (
    title_panel /
      panel_grid /
      caption_panel
  ) +
    patchwork::plot_layout(
      heights = c(
        0.075,
        1,
        0.055
      )
    )
  
  output_path <- file.path(
    figure_folder,
    paste0(
      "appendix_mice_convergence_",
      suffix,
      ".png"
    )
  )
  
  ggplot2::ggsave(
    filename = output_path,
    plot = p,
    width = 10,
    height = 7,
    units = "in",
    dpi = PLOT_DPI,
    limitsize = FALSE,
    bg = "white"
  )
  
  invisible(p)
}

mice_trace_plot_zh <- make_trace_plot(trace_df, "zh")
mice_trace_plot_en <- make_trace_plot(trace_df, "en")

make_density_df <- function(
    imp_object,
    vars
) {
  original_data <- imp_object$data
  vars <- intersect(
    vars,
    names(original_data)
  )

  purrr::map_dfr(
    vars,
    function(v) {
      observed_df <- original_data %>%
        dplyr::transmute(
          variable = v,
          value = .data[[v]],
          type = "observed",
          chain = "observed"
        ) %>%
        dplyr::filter(
          !is.na(value),
          is.finite(value)
        )

      if (
        !is.null(
          imp_object$imp[[v]]
        ) &&
          nrow(
            imp_object$imp[[v]]
          ) > 0L
      ) {
        imputed_df <- imp_object$imp[[v]] %>%
          as.data.frame() %>%
          tidyr::pivot_longer(
            cols = dplyr::everything(),
            names_to = "chain",
            values_to = "value"
          ) %>%
          dplyr::transmute(
            variable = v,
            value = value,
            type = "imputed",
            chain = chain
          ) %>%
          dplyr::filter(
            !is.na(value),
            is.finite(value)
          )

        dplyr::bind_rows(
          observed_df,
          imputed_df
        )
      } else {
        observed_df
      }
    }
  )
}

density_df <- make_density_df(
  imputed_data_mids,
  mice_targets
) %>%
  dplyr::mutate(
    variable_label = factor(
      dplyr::recode(
        variable,
        Glu_AC = "Glu-AC",
        LDL = "LDL",
        TG = "TG",
        BMI = "BMI"
      ),
      levels = c(
        "Glu-AC",
        "LDL",
        "TG",
        "BMI"
      )
    )
  )

make_density_plot <- function(
    density_data,
    lang
) {
  if (nrow(density_data) == 0L) {
    return(NULL)
  }
  
  if (lang == "zh") {
    dat_plot <- density_data %>%
      dplyr::mutate(
        type_plot = factor(
          dplyr::recode(
            type,
            observed = "觀測值",
            imputed = "插補值"
          ),
          levels = c(
            "觀測值",
            "插補值"
          )
        )
      )
    
    title_text <- "觀測值與插補值分布比較"
    
    caption_text <- paste0(
      "比較 Glu-AC、LDL、TG 與 BMI 的觀測值及 10 組插補值分布；",
      "HbA1c 因無缺失值，未納入本圖。"
    )
    
    y_text <- "密度"
    suffix <- "zh"
    
    legend_breaks <- c(
      "觀測值",
      "插補值"
    )
    
  } else if (lang == "en") {
    dat_plot <- density_data %>%
      dplyr::mutate(
        type_plot = factor(
          dplyr::recode(
            type,
            observed = "Observed values",
            imputed = "Imputed values"
          ),
          levels = c(
            "Observed values",
            "Imputed values"
          )
        )
      )
    
    title_text <-
      "Comparison of observed and imputed distributions"
    
    caption_text <- paste0(
      "Observed values are compared with the 10 imputation chains ",
      "for Glu-AC, LDL, TG, and BMI. ",
      "HbA1c is not shown because it had no missing values."
    )
    
    y_text <- "Density"
    suffix <- "en"
    
    legend_breaks <- c(
      "Observed values",
      "Imputed values"
    )
    
  } else {
    stop(
      "lang 僅可設定為 'zh' 或 'en'。"
    )
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = dat_plot %>%
        dplyr::filter(
          type == "imputed"
        ),
      ggplot2::aes(
        x = value,
        group = chain,
        color = type_plot
      ),
      linewidth = 0.45,
      alpha = 0.65,
      na.rm = TRUE
    ) +
    ggplot2::geom_density(
      data = dat_plot %>%
        dplyr::filter(
          type == "observed"
        ),
      ggplot2::aes(
        x = value,
        color = type_plot
      ),
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~ variable_label,
      ncol = 2,
      scales = "free"
    ) +
    ggplot2::scale_color_discrete(
      breaks = legend_breaks
    ) +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = y_text,
      color = NULL,
      caption = caption_text
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      
      strip.text = ggplot2::element_text(
        size = 11,
        face = "bold"
      ),
      
      plot.caption = ggplot2::element_text(
        size = 8,
        hjust = 0
      )
    )
  
  output_path <- file.path(
    figure_folder,
    paste0(
      "appendix_observed_vs_imputed_density_",
      suffix,
      ".png"
    )
  )
  
  ggplot2::ggsave(
    filename = output_path,
    plot = p,
    width = 9,
    height = 7,
    units = "in",
    dpi = PLOT_DPI,
    limitsize = FALSE,
    bg = "white"
  )
  
  invisible(p)
}

density_plot_zh <- make_density_plot(
  density_df,
  "zh"
)

density_plot_en <- make_density_plot(
  density_df,
  "en"
)

# 原始尺度與自然對數尺度分布：
# 使用多重插補前的非缺失觀測值，以呈現實際觀測分布並說明 log 轉換理由；
# 不使用第一組完成資料，避免將插補值誤標示為「原始」資料。
distribution_source_data <- dt_model

log_distribution_df <- distribution_source_data %>%
  dplyr::select(
    dplyr::all_of(
      primary_outcomes
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "variable",
    values_to = "value_original"
  ) %>%
  dplyr::filter(
    !is.na(value_original),
    is.finite(value_original),
    value_original > 0
  ) %>%
  dplyr::mutate(
    value_log = log(value_original),
    variable_label = factor(
      outcome_labels[variable],
      levels = unname(
        outcome_labels[
          primary_outcomes
        ]
      )
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      value_original,
      value_log
    ),
    names_to = "scale",
    values_to = "value"
  )

make_log_distribution_plot <- function(
    dat,
    lang
) {
  if (lang == "zh") {
    title_text <-
      "四項主要代謝指標觀測值之原始尺度與自然對數尺度分布"
    caption_text <- "本圖使用多重插補前之非缺失且大於 0 的觀測值。"
    x_text <- "數值"
    y_text <- "密度"
    suffix <- "zh"
    scale_levels <- c(
      "原始尺度",
      "自然對數尺度"
    )
  } else {
    title_text <-
      "Distributions of observed metabolic indicators on original and natural-log scales"
    caption_text <- "The figure uses nonmissing positive observations before multiple imputation."
    x_text <- "Value"
    y_text <- "Density"
    suffix <- "en"
    scale_levels <- c(
      "Original scale",
      "Natural-log scale"
    )
  }

  dat_plot <- dat %>%
    dplyr::mutate(
      scale_label = dplyr::case_when(
        scale == "value_original" ~
          scale_levels[1],
        TRUE ~
          scale_levels[2]
      ),
      panel_label = paste0(
        variable_label,
        if (
          lang == "zh"
        ) {
          "："
        } else {
          ": "
        },
        scale_label
      )
    )

  # 明確固定版面：
  # 第一列為 HbA1c、Glu-AC、LDL、TG 的原始尺度；
  # 第二列為 HbA1c、Glu-AC、LDL、TG 的自然對數尺度。
  separator_text <- if (lang == "zh") "：" else ": "

  desired_panel_levels <- c(
    paste0(
      unname(outcome_labels[primary_outcomes]),
      separator_text,
      scale_levels[1]
    ),
    paste0(
      unname(outcome_labels[primary_outcomes]),
      separator_text,
      scale_levels[2]
    )
  )

  dat_plot <- dat_plot %>%
    dplyr::mutate(
      panel_label = factor(
        panel_label,
        levels = desired_panel_levels
      )
    )

  p <- ggplot2::ggplot(
    dat_plot,
    ggplot2::aes(
      x = value
    )
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(
        y = ggplot2::after_stat(
          density
        )
      ),
      bins = 35
    ) +
    ggplot2::facet_wrap(
      ~ panel_label,
      ncol = 4,
      scales = "free"
    ) +
    ggplot2::labs(
      title = title_text,
      x = x_text,
      y = y_text,
      caption = caption_text
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      strip.text = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      plot.caption = ggplot2::element_text(
        size = 8,
        hjust = 0
      )
    )

  ggplot2::ggsave(
    file.path(
      figure_folder,
      paste0(
        "appendix_original_vs_log_distribution_",
        suffix,
        ".png"
      )
    ),
    p,
    width = 12,
    height = 6.5,
    dpi = PLOT_DPI,
    bg = "white"
  )

  p
}

log_distribution_plot_zh <- make_log_distribution_plot(
  log_distribution_df,
  "zh"
)

log_distribution_plot_en <- make_log_distribution_plot(
  log_distribution_df,
  "en"
)

# ============================================================
# 12. 10 組插補資料、imp_long 與附錄 A3－A4
# ============================================================

imp_list_split <- lapply(
  seq_along(imp_list),
  function(i) {
    imp_list[[i]] %>%
      dplyr::left_join(
        split_manifest %>%
          dplyr::select(
            ID,
            split
          ),
        by = "ID",
        relationship = "many-to-one"
      ) %>%
      dplyr::mutate(
        .imp = i,
        data_role = paste0(
          "inference-only 10-imputation completed data; ",
          "the full eligible cohort is used for GEE/LME statistical inference; ",
          "split labels are retained only for data alignment and QC; ",
          "this file is not used for empirical prediction evaluation"
        )
      )
  }
)

imp_long_full <- dplyr::bind_rows(
  imp_list_split
)

shared_model_vars <- c(
  ".imp",
  "row_id",
  "split",
  "data_role",
  "ID",
  "Group",
  "報告日期",
  "Time_in_Days",
  "Time_in_Years",
  "Time_Category",
  "HbA1c",
  "Glu_AC",
  "LDL",
  "TG",
  "BMI",
  "Age_baseline",
  "Sex",
  "CVD_History",
  "HTN_History",
  "CKD_History",
  "Liver_History",
  "HbA1c_observed",
  "Glu_AC_observed",
  "LDL_observed",
  "TG_observed",
  "BMI_observed"
)

imp_long <- imp_long_full %>%
  dplyr::select(
    dplyr::all_of(
      shared_model_vars
    )
  ) %>%
  dplyr::arrange(
    .imp,
    ID,
    Time_in_Years,
    row_id
  )

# canonical RDS：與 canonical CSV 同名邏輯。
saveRDS(
  imp_long,
  file = file.path(
    rds_folder,
    "imputed_long_with_split_type2.rds"
  )
)

# 完整 cohort MICE 完成資料：僅供 GEE/LME 統計推論與插補後描述性附錄使用。
readr::write_csv(
  imp_long,
  file.path(
    dl_folder,
    "imputed_long_with_split_type2.csv"
  )
)

# ------------------------------------------------------------
# 11.0A 無資料洩漏之預測用資料
# ------------------------------------------------------------
# 原則：
# 1. 先使用既定病患層級 split。
# 2. 四項 outcome 不做 MICE；模型序列僅使用 target 之前實際觀測值。
# 3. BMI 僅以 train set 估計之組別中位數補值；若組別中位數不可得，
#    再使用 train set overall 中位數。validation/test 不參與補值參數估計。
# 4. 預測資料只有一個 deterministic preprocessing 版本（.imp = 1），
#    不宣稱為多重插補資料。
prediction_bmi_group_median <- dt_model %>%
  dplyr::left_join(
    split_manifest %>% dplyr::select(ID, split),
    by = "ID",
    relationship = "many-to-one"
  ) %>%
  dplyr::filter(split == "train", !is.na(BMI)) %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    BMI_train_group_median = stats::median(BMI, na.rm = TRUE),
    .groups = "drop"
  )

prediction_bmi_overall_median <- dt_model %>%
  dplyr::left_join(
    split_manifest %>% dplyr::select(ID, split),
    by = "ID",
    relationship = "many-to-one"
  ) %>%
  dplyr::filter(split == "train", !is.na(BMI)) %>%
  dplyr::summarise(
    value = stats::median(BMI, na.rm = TRUE)
  ) %>%
  dplyr::pull(value)

if (length(prediction_bmi_overall_median) != 1L ||
    !is.finite(prediction_bmi_overall_median)) {
  stop("train set 無法估計 BMI overall 中位數，無法建立 prediction-safe 資料。")
}

prediction_long <- dt_model %>%
  dplyr::left_join(
    split_manifest %>% dplyr::select(ID, split),
    by = "ID",
    relationship = "many-to-one"
  ) %>%
  dplyr::left_join(
    prediction_bmi_group_median,
    by = "Group",
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    BMI_raw_before_prediction_fill = BMI,
    BMI = dplyr::coalesce(
      BMI,
      BMI_train_group_median,
      prediction_bmi_overall_median
    ),
    BMI_prediction_filled = is.na(BMI_raw_before_prediction_fill),
    .imp = 1L,
    data_role = paste0(
      "prediction-safe deterministic data; ",
      "same-patient same-date records were aggregated before model preprocessing; ",
      "patient split was fixed before prediction preprocessing; ",
      "outcomes were not imputed; ",
      "BMI missing values were filled using training-set group-specific medians ",
      "with the training-set overall median as fallback"
    )
  ) %>%
  dplyr::select(
    dplyr::all_of(shared_model_vars),
    BMI_prediction_filled
  ) %>%
  dplyr::arrange(ID, Time_in_Years, row_id)

prediction_safe_qc <- tibble::tibble(
  check = c(
    "prediction_imp_n",
    "duplicate_row_id",
    "patient_multiple_split",
    "remaining_BMI_missing",
    "outcome_values_changed_from_aggregated_dt_model"
  ),
  value = c(
    dplyr::n_distinct(prediction_long$.imp),
    nrow(prediction_long %>% dplyr::count(row_id) %>% dplyr::filter(n > 1L)),
    nrow(prediction_long %>% dplyr::distinct(ID, split) %>% dplyr::count(ID) %>% dplyr::filter(n != 1L)),
    sum(is.na(prediction_long$BMI)),
    sum(vapply(
      primary_outcomes,
      function(v) {
        raw_v <- dt_model[[v]][match(prediction_long$row_id, dt_model$row_id)]
        sum(
          xor(is.na(raw_v), is.na(prediction_long[[v]])) |
            (!is.na(raw_v) & !is.na(prediction_long[[v]]) & raw_v != prediction_long[[v]]),
          na.rm = TRUE
        )
      },
      numeric(1)
    ))
  ),
  expected = c(1L, 0L, 0L, 0L, 0L)
) %>%
  dplyr::mutate(passed = value == expected)

if (!all(prediction_safe_qc$passed)) {
  print(prediction_safe_qc)
  stop("prediction-safe 資料 QC 未通過。")
}

readr::write_csv(
  prediction_long,
  file.path(dl_folder, "prediction_long_with_split_type2.csv")
)

readr::write_csv(
  prediction_safe_qc,
  file.path(qc_folder, "prediction_safe_data_qc.csv")
)

readr::write_csv(
  prediction_bmi_group_median,
  file.path(qc_folder, "prediction_BMI_train_group_median.csv")
)

# canonical inference 檔為 imputed_long_with_split_type2.csv，
# split 可直接由 split 欄位篩選。

readr::write_csv(
  model_data_usage,
  file.path(
    dl_folder,
    "analysis_data_usage_manifest.csv"
  )
)

readr::write_csv(
  imputation_split_scope,
  file.path(
    dl_folder,
    "imputation_split_scope_manifest.csv"
  )
)

imp_long_split_check <- imp_long %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    n_split = dplyr::n_distinct(split),
    n_imp = dplyr::n_distinct(.imp),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_split != 1L |
      n_imp != MICE_M
  )

imp_long_count_qc <- tibble::tibble(
  item = c(
    "dt_model_row_n",
    "number_of_imputations",
    "expected_imp_long_row_n",
    "actual_imp_long_row_n",
    "unique_row_id_n_per_imputation"
  ),
  value = c(
    nrow(dt_model),
    MICE_M,
    nrow(dt_model) * MICE_M,
    nrow(imp_long),
    n_distinct(
      imp_long$row_id[
        imp_long$.imp == 1L
      ]
    )
  )
)

imputed_log_eligibility_qc <- purrr::map_dfr(
  primary_outcomes,
  function(v) {
    x <- imp_long[[v]]

    tibble::tibble(
      variable = v,
      completed_value_n = sum(!is.na(x)),
      remaining_missing_n = sum(is.na(x)),
      nonpositive_n = sum(x <= 0, na.rm = TRUE),
      min_completed_value = ifelse(
        any(!is.na(x)),
        min(x, na.rm = TRUE),
        NA_real_
      )
    )
  }
)

if (
  any(
    imputed_log_eligibility_qc$
      remaining_missing_n > 0L
  ) ||
    any(
      imputed_log_eligibility_qc$
        nonpositive_n > 0L
    )
) {
  stop(
    "完成插補資料中仍有缺失或 <= 0 的主要 outcome，",
    "不可進入後續 log 模型；請檢查 imputed_log_eligibility_qc。"
  )
}

# ------------------------------------------------------------
# 11.1 附錄 A3：第一組補值資料之基線特徵
# ------------------------------------------------------------

imp1_baseline <- imp_list_split[[1]] %>%
  dplyr::arrange(
    ID,
    報告日期,
    就診號
  ) %>%
  dplyr::group_by(ID) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    BMI_baseline_imp = BMI
  )

baseline_results_imp1 <- build_baseline_results(
  baseline_dat = imp1_baseline,
  bmi_var = "BMI_baseline_imp"
)

appendix_A3_ready <- baseline_results_imp1$ready
appendix_A3_cont <- baseline_results_imp1$continuous
appendix_A3_cont_p <- baseline_results_imp1$continuous_p
appendix_A3_cat <- baseline_results_imp1$categorical
appendix_A3_cat_p <- baseline_results_imp1$categorical_p

# ------------------------------------------------------------
# 11.2 附錄 A4：第一組補值資料之 T1－T19 敘述統計
# ------------------------------------------------------------

imp1_patient_period <- imp_list_split[[1]] %>%
  dplyr::select(
    ID,
    Group,
    Time_Category,
    dplyr::all_of(
      primary_outcomes
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      primary_outcomes
    ),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  dplyr::group_by(
    ID,
    Group,
    Time_Category,
    outcome
  ) %>%
  dplyr::summarise(
    value = mean(
      value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

appendix_A4_results <- purrr::map(
  primary_outcomes,
  ~ make_outcome_time_results(
    imp1_patient_period,
    .x
  )
)

names(appendix_A4_results) <- primary_outcomes

appendix_A4_imp1_primary_by_time <- purrr::map_dfr(
  appendix_A4_results,
  "ready"
) %>%
  dplyr::arrange(
    match(
      outcome,
      primary_outcomes
    ),
    Time_Num
  )

appendix_A4_imp1_detail <- purrr::map_dfr(
  appendix_A4_results,
  "detail"
) %>%
  dplyr::arrange(
    match(
      outcome,
      primary_outcomes
    ),
    Time_Num,
    Group
  )

# ============================================================
# 13. 分析規則與最終自動檢查
# ============================================================

analysis_rules <- tibble::tribble(
  ~item, ~rule,
  "Main input",
  paste0(
    "If CSV and XLSX both pass the current schema, analysis-critical columns must be identical; ",
    "otherwise execution stops. XLSX is the first-priority source."
  ),
  "Primary cohort",
  "Type 2 diabetes with baseline age >=18 years",
  "Age missing",
  "Excluded before adult cohort; not imputed by MICE",
  "Exposure definition",
  "Patient must remain fixed Control or fixed Statin across observed visits",
  "Time categories",
  paste0(
    "Raw data are checked through all observed categories. ",
    "T20 and T21 are reported in QC and excluded; main analysis retains T1-T19."
  ),
  "Longitudinal eligibility",
  "At least two distinct report dates within T1-T19",
  "Same-date aggregation",
  paste0(
    "Before MICE and all downstream models, multiple records from the same patient ",
    "on the same report date are aggregated into one time point. Continuous laboratory ",
    "values and BMI are averaged over nonmissing values; fixed patient characteristics ",
    "must be consistent within date."
  ),
  "Table 4.1 baseline",
  paste0(
    "First eligible observation; includes Age, Sex, baseline BMI and comorbidities. ",
    "BP and metabolic outcomes are not placed in Table 4.1."
  ),
  "Tables 4.2-4.5",
  paste0(
    "Patient-period unit; repeated records within the same 3-month interval ",
    "are averaged before mean/SD, missingness and Wilcoxon comparison."
  ),
  "Outcome cleaning",
  paste0(
    "Clinical-range, adjacent-jump and IQR rules are flags only; ",
    "values are not automatically deleted."
  ),
  "MICE targets",
  "Glu-AC, LDL, TG and time-varying BMI by PMM",
  "MICE exclusions",
  "Age, Sex, BP, HDL, TC, ID and time variables are not imputed",
  "Inference RDS objects",
  paste0(
    "imputed_data_mids_type2_inference.rds contains the mids object; ",
    "imp_list_type2_inference.rds contains 10 completed data sets; ",
    "imputed_long_with_split_type2.rds contains the canonical long-format completed data."
  ),
  "imputed_long_with_split_type2.csv",
  paste0(
    "Full-cohort MICE completed data used only for GEE/LME statistical inference ",
    "and imputation-based descriptive appendices."
  ),
  "prediction_long_with_split_type2.csv",
  paste0(
    "Prediction-safe deterministic data used by GEE/LME prediction and DL. ",
    "Patient split is fixed first; outcomes are not imputed; BMI is filled using train-only medians."
  ),
  "Prediction input",
  paste0(
    "prediction_long_with_split_type2.csv is the sole standard input for ",
    "GEE/LME prediction and the three deep-learning models. ",
    "Redundant dl_raw_* files are not exported."
  ),
  "BP role",
  paste0(
    "Optional input only. If present, BP is used for appendix data availability; ",
    "if absent in the latest preprocessing file, no BP missingness percentage is fabricated."
  ),
  "Figure languages",
  "Paired Chinese and English combined figures; individual outcome PNG files are not exported; all-ID trajectory is disabled by default due overplotting"
)

stopifnot(
  nrow(repeat_sequence_qc) == 0L,
  nrow(final_event_duplicate_qc) == 0L,
  nrow(same_patient_same_time_qc) == 0L,
  nrow(same_date_duplicate_after_qc) == 0L,
  nrow(age_baseline_constancy_qc) == 0L,
  nrow(split_check) == 0L,
  nrow(imp_long_split_check) == 0L,
  !any(is.na(dt_model$ID)),
  !any(is.na(dt_model$就診號)),
  !any(is.na(dt_model$報告日期)),
  !any(is.na(dt_model$Group)),
  !any(is.na(dt_model$Sex)),
  !any(is.na(dt_model$Age_baseline)),
  all(
    dt_model$Age_baseline >= 18
  ),
  all(
    as.integer(
      stringr::str_remove(
        as.character(
          dt_model$Time_Category
        ),
        "T"
      )
    ) >= 1L
  ),
  all(
    as.integer(
      stringr::str_remove(
        as.character(
          dt_model$Time_Category
        ),
        "T"
      )
    ) <= 19L
  ),
  all(
    dt_model %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(
        n_date = dplyr::n_distinct(
          報告日期
        ),
        .groups = "drop"
      ) %>%
      dplyr::pull(n_date) >= 2L
  ),
  rosis_exclusion$
    rosis_misclassified_as_statin_n[1L] == 0L,
  dplyr::n_distinct(
    dt_model$row_id
  ) == nrow(dt_model),
  nrow(imp_long) ==
    nrow(dt_model) * MICE_M,
  all(
    imputation_completeness_qc$
      remaining_missing_n == 0L
  ),
  all(
    imputation_observed_preservation_qc$
      observed_value_changed_n == 0L
  ),
  nrow(baseline_cat_summary) <
    100L,
  nrow(table4_1_ready) <
    30L,
  all(
    log_eligibility_qc$nonpositive_n == 0L
  ),
  all(
    imputed_log_eligibility_qc$
      remaining_missing_n == 0L
  ),
  all(
    imputed_log_eligibility_qc$
      nonpositive_n == 0L
  )
)

# ============================================================
# 14. 輸出
# ============================================================

analysis_run_manifest <- tibble::tibble(
  item = c(
    "script_version",
    "run_time_asia_taipei",
    "main_data_path",
    "final_patient_n",
    "final_observation_n",
    "control_patient_n",
    "statin_patient_n",
    "max_raw_time_category",
    "excluded_after_T19_observation_n",
    "mice_m",
    "mice_maxit"
  ),
  value = as.character(
    c(
      SCRIPT_VERSION,
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S",
        tz = LOCAL_TZ
      ),
      main_data_path,
      dplyr::n_distinct(
        dt_model$ID
      ),
      nrow(dt_model),
      group_count_final_patient$
        patient_n[
          as.character(
            group_count_final_patient$Group
          ) == "Control"
        ],
      group_count_final_patient$
        patient_n[
          as.character(
            group_count_final_patient$Group
          ) == "Statin"
        ],
      max(
        time_category_count_before$time_num,
        na.rm = TRUE
      ),
      sum(
        time_category_excluded_after_t19$
          observation_n
      ),
      MICE_M,
      MICE_MAXIT
    )
  )
)

summary_sheets <- list(
  run_manifest = analysis_run_manifest,
  main_source = main_source_summary,
  main_source_consistency =
    main_source_consistency_qc,
  optional_input_qc = main_optional_column_qc,
  patient_flow = patient_flow,
  type_patient = type_count_raw_patient,
  type_observation = type_count_raw_observation,
  adult_eligibility = adult_eligibility_summary,
  under18_detail = under18_patient_detail,
  age_missing_detail = age_missing_patient_detail,
  fixed_group_patient = fixed_exposure_group_patient,
  fixed_group_observation =
    fixed_exposure_group_observation,
  longitudinal_elig =
    longitudinal_eligibility_summary,
  time_range_qc = time_category_range_qc,
  time_before = time_category_count_before,
  time_by_group =
    time_category_count_by_group_before,
  time_decision = time_category_decision,
  time_excluded =
    time_category_excluded_after_t19,
  time_exclusion_total =
    time_category_exclusion_totals,
  time_final = time_category_count_final,
  statin_dictionary = statin_drug_dictionary,
  statin_drug_check = statin_drug_check,
  known_brand_qc = known_statin_brand_qc,
  unclassified_lipid_qc =
    unclassified_lipid_name_qc,
  rosis_exclusion = rosis_exclusion,
  ezetimibe_drugs = ezetimibe_drug_check,
  ezetimibe_exposure =
    ezetimibe_exposure_check,
  exposure_summary =
    excluded_exposure_summary,
  group_patient = group_count_final_patient,
  group_observation =
    group_count_final_observation,
  repeat_summary =
    repeat_patient_summary_final,
  same_date_aggregation =
    same_date_aggregation_summary,
  same_date_detail =
    same_date_duplicate_before_qc,
  same_date_fixed_qc =
    same_date_fixed_constancy_qc,
  final_event_duplicate =
    final_event_duplicate_qc,
  same_time_duplicate =
    same_patient_same_time_qc,
  same_date_multi_visit =
    same_patient_date_multi_visit_qc,
  quality_flags = quality_flag_summary,
  missing_outcome_all =
    missing_response_overall,
  missing_outcome_group =
    missing_response_by_group,
  missing_outcome_time =
    missing_response_by_group_time,
  missing_covariate_all =
    missing_covariate_overall,
  missing_covariate_group =
    missing_covariate_by_group,
  bp_missing_event_all =
    bp_missing_event_overall,
  bp_missing_event_group =
    bp_missing_event_by_group,
  bp_missing_base_all =
    bp_missing_baseline_overall,
  bp_missing_base_group =
    bp_missing_baseline_by_group,
  fig4_1_stats = mean_sd_ci_figure4_1,
  table4_1_ready = table4_1_ready,
  table4_1_cont = baseline_cont_summary,
  table4_1_cont_p = baseline_cont_p,
  table4_1_cat = baseline_cat_summary,
  table4_1_cat_p = baseline_cat_p,
  table4_2_HbA1c = table4_2_HbA1c,
  detail4_2_HbA1c =
    table4_2_HbA1c_detail,
  table4_3_Glu_AC = table4_3_Glu_AC,
  detail4_3_Glu_AC =
    table4_3_Glu_AC_detail,
  table4_4_LDL = table4_4_LDL,
  detail4_4_LDL =
    table4_4_LDL_detail,
  table4_5_TG = table4_5_TG,
  detail4_5_TG =
    table4_5_TG_detail,
  appendix_A1_HDL = appendix_A1_HDL,
  detail_A1_HDL = appendix_A1_HDL_detail,
  appendix_A2_TC = appendix_A2_TC,
  detail_A2_TC = appendix_A2_TC_detail,
  appendix_A3_ready = appendix_A3_ready,
  appendix_A3_cont = appendix_A3_cont,
  appendix_A3_cont_p =
    appendix_A3_cont_p,
  appendix_A3_cat = appendix_A3_cat,
  appendix_A3_cat_p =
    appendix_A3_cat_p,
  appendix_A4_imp1 =
    appendix_A4_imp1_primary_by_time,
  appendix_A4_detail =
    appendix_A4_imp1_detail,
  all_trajectory =
    all_trajectory_summary,
  sampled_trajectory =
    sampled_trajectory_id_summary,
  mice_method = mice_setup_method,
  mice_predictors =
    mice_predictor_matrix_export,
  mice_logged_events = mice_logged_events,
  mice_complete_qc =
    imputation_completeness_qc,
  mice_observed_qc =
    imputation_observed_preservation_qc,
  imp_long_count_qc = imp_long_count_qc,
  imputed_log_qc =
    imputed_log_eligibility_qc,
  split_patient = split_count_patient,
  split_observation =
    split_count_observation,
  split_missing =
    split_outcome_missing,
  split_check = split_check,
  imp_split_check =
    imp_long_split_check,
  comorbidity_def =
    comorbidity_definition,
  transformation_rules =
    transformation_rules,
  model_data_usage =
    model_data_usage,
  imputation_split_scope =
    imputation_split_scope,
  log_eligibility_qc =
    log_eligibility_qc,
  analysis_rules = analysis_rules
)

readr::write_csv(
  same_date_aggregation_summary,
  file.path(
    qc_folder,
    "same_date_aggregation_summary.csv"
  )
)

readr::write_csv(
  same_date_duplicate_before_qc,
  file.path(
    qc_folder,
    "same_date_aggregation_detail.csv"
  )
)

readr::write_csv(
  same_date_fixed_constancy_qc,
  file.path(
    qc_folder,
    "same_date_fixed_variable_qc.csv"
  )
)

writexl::write_xlsx(
  summary_sheets,
  file.path(
    table_folder,
    "table_outputs_type2_adult_main_analysis.xlsx"
  )
)

dt_model_export <- prepare_export_data(
  dt_model
)

readr::write_csv(
  dt_model_export,
  file.path(
    output_root,
    "dt_model_raw_type2_adult.csv"
  )
)

writexl::write_xlsx(
  list(
    dt_model_raw = dt_model_export
  ),
  file.path(
    output_root,
    "dt_model_raw_type2_adult.xlsx"
  )
)

message(
  "表格.R 已執行完成，且最終自動檢查通過。"
)

message(
  "主要輸出資料夾：",
  output_root
)

message("重要檔案：")

message(
  "1. ",
  file.path(
    table_folder,
    "table_outputs_type2_adult_main_analysis.xlsx"
  )
)

message(
  "2. ",
  file.path(
    dl_folder,
    "imputed_long_with_split_type2.csv"
  ),
  "（僅供推論）"
)

message(
  "2A. ",
  file.path(
    dl_folder,
    "prediction_long_with_split_type2.csv"
  ),
  "（供 GEE/LME 預測與 DL）"
)

message(
  "3. ",
  file.path(
    rds_folder,
    "imputed_data_mids_type2_inference.rds"
  )
)

message(
  "4. ",
  file.path(
    rds_folder,
    "imp_list_type2_inference.rds"
  )
)

message(
  "5. ",
  file.path(
    rds_folder,
    "imputed_long_with_split_type2.rds"
  )
)

message(
  "6. ",
  file.path(
    figure_folder,
    "combined_mean_trend_ci_primary_outcomes_zh.png"
  )
)

message(
  "7. ",
  file.path(
    figure_folder,
    "combined_mean_trend_ci_primary_outcomes_en.png"
  )
)

message(
  "8. ",
  file.path(
    figure_folder,
    "appendix_original_vs_log_distribution_zh.png"
  )
)

message(
  "9. ",
  file.path(
    figure_folder,
    "appendix_original_vs_log_distribution_en.png"
  )
)

message(
  "T20/T21 排除筆數：",
  sum(
    time_category_excluded_after_t19$
      observation_n
  )
)
# T20/T21 排除筆數：55
