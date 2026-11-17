# ============================================================
# 目的：
# 1. 建立糖尿病 cohort 與 Type 1 / Type 2 分類。
# 2. 以 HbA1c 檢驗事件作為縱向追蹤資料之 anchor。
# 3. 串接 Glu-AC、LDL、TG，以及描述性備用之 HDL、TC。
# 4. 避免逐項 full_join 造成 many-to-many 資料膨脹。
# 5. 同一事件同一檢驗出現多個不同數值時，固定以中位數彙整，
#    並將衝突事件完整輸出供 QA。
# 6. 性別以病人基本資料之病患層級眾數為主，檢驗資料僅作備援與 QC。
# 7. 保留基本資料完整日期時間；身高與體重必須取自同一筆真實紀錄，
#    不再分別取中位數後組合成不存在的 BMI。
# 8. 每一檢驗事件使用「當日或檢驗日前最近一次」有效身高體重配對。
#    不設定任意 look-back 截止日，但輸出 >365 日與 >730 日之 QA。
# 9. 主檔移除精確出生日期及多餘 QA 欄位，降低準識別資訊暴露。
# 10. 所有日期時間以 Asia/Taipei 處理，並以當地時間字串輸出，避免 UTC 位移。
# 11. 同時輸出 xlsx 與 csv；後續分析前應確認兩者由同一次執行產生。
#
# 重要原則：
# - 主要縱向 outcome：HbA1c、Glu-AC、LDL、TG。
# - HDL 與 TC 僅保留作描述性／附錄備用，不進入主要模型。
# - Creatinine、GPT、Uric acid 與其他未使用檢驗不讀取。
# - 本程式不使用相鄰差值、IQR 或三點式規則自動刪除 outcome。
# - 「<a」以 a/2 代入；「>a」以門檻 a 代入並保留 censor flag。
# - 寬鬆技術性範圍只用於辨識明顯輸入錯誤，不作臨床異常值判定。
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(writexl)
  library(lubridate)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ============================================================
# 0. 路徑與固定設定
# ============================================================

raw_root <- "raw_data"

lab_root <- file.path(
  raw_root,
  "郝醫師-糖尿病品質指標資料2020至2024年_114.1.20 _刪除姓名和身分證字號",
  "血液"
)

save_folder <- "./糖尿病患抽血資料整理"

# 所有日期時間均以台灣當地時間處理。
LOCAL_TZ <- "Asia/Taipei"
SCRIPT_VERSION <- "data-linkage-qc-enhanced-v2"

dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(lab_root)) {
  stop("找不到血液檢驗資料夾：", lab_root)
}

# 目前原始檢驗資料沒有可穩定使用的採檢時間，因此以日期定義事件。
EVENT_KEYS <- c("病歷號", "就診號", "報告日期")

MAIN_LABS <- c("HbA1C", "Glucose_AC", "LDL-C", "Triglycerides")

lab_map <- c(
  "91111100-HbA1C" = "HbA1C",
  "91111105-Glucose AC" = "Glucose_AC",
  "90111614-LDL-C" = "LDL-C",
  "90111200-三酸甘油脂" = "Triglycerides",
  "90111617-HDL-C" = "HDL-C",
  "90111200-總膽固醇" = "Total_cholesterol"
)

# 僅用於篩除明顯不可能之輸入值，不代表臨床正常範圍。
ABS_HEIGHT_RANGE_CM <- c(40, 250)
ABS_WEIGHT_RANGE_KG <- c(2, 500)
ABS_BMI_RANGE <- c(5, 100)

main_output_xlsx_path <- file.path(
  save_folder, "糖尿病患回診資料_with_type.xlsx"
)
main_output_csv_path <- file.path(
  save_folder, "糖尿病患回診資料_with_type.csv"
)
type_output_path <- file.path(
  save_folder, "patient_diabetes_type_summary.xlsx"
)
lab_conflict_output_path <- file.path(
  save_folder, "lab_value_conflicts.xlsx"
)
qc_output_path <- file.path(
  save_folder, "preprocess_qc_summary.xlsx"
)

# ============================================================
# 1. 工具函數
# ============================================================

read_xlsx_hospital <- function(path, drop_first_col = TRUE) {
  if (!file.exists(path)) {
    stop("找不到檔案：", path)
  }
  
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
    select(-matches("^\\.\\.\\.\\d+$"))
}

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

standardize_id <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

standardize_sex <- function(x) {
  x <- str_squish(as.character(x))
  
  case_when(
    x %in% c("M", "m", "男", "1", "male", "Male", "MALE") ~ "男",
    x %in% c("F", "f", "女", "0", "female", "Female", "FEMALE") ~ "女",
    TRUE ~ NA_character_
  )
}

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  
  if (inherits(x, c("POSIXct", "POSIXt"))) {
    return(as.Date(x, tz = LOCAL_TZ))
  }
  
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  
  parsed <- suppressWarnings(
    parse_date_time(
      x_chr,
      orders = c(
        "Ymd HMS", "Ymd HM", "Ymd",
        "ymd HMS", "ymd HM", "ymd",
        "mdy HMS", "mdy HM", "mdy"
      ),
      tz = LOCAL_TZ
    )
  )
  
  out <- as.Date(parsed, tz = LOCAL_TZ)
  
  # Excel serial date 備援。
  x_num <- suppressWarnings(as.numeric(x_chr))
  is_excel_serial <- (
    is.na(out) &
      !is.na(x_num) &
      x_num >= 20000 &
      x_num <= 80000 &
      !str_detect(x_chr, "[-/:]")
  )
  
  out[is_excel_serial] <- as.Date(
    floor(x_num[is_excel_serial]),
    origin = "1899-12-30"
  )
  
  out
}

parse_datetime_safe <- function(x) {
  if (inherits(x, c("POSIXct", "POSIXt"))) {
    return(as.POSIXct(x, tz = LOCAL_TZ))
  }
  
  if (inherits(x, "Date")) {
    return(as.POSIXct(x, tz = LOCAL_TZ))
  }
  
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  
  out <- as.POSIXct(
    rep(NA_character_, length(x_chr)),
    tz = LOCAL_TZ
  )
  
  # 醫院基本資料常見 14 碼格式：YYYYMMDDhhmmss。
  is_14_digit <- !is.na(x_chr) & str_detect(x_chr, "^\\d{14}$")
  
  if (any(is_14_digit)) {
    out[is_14_digit] <- suppressWarnings(
      as.POSIXct(
        x_chr[is_14_digit],
        format = "%Y%m%d%H%M%S",
        tz = LOCAL_TZ
      )
    )
  }
  
  need_fallback <- is.na(out) & !is.na(x_chr)
  
  if (any(need_fallback)) {
    out[need_fallback] <- suppressWarnings(
      as.POSIXct(
        parse_date_time(
          x_chr[need_fallback],
          orders = c(
            "Ymd HMS", "Ymd HM", "Ymd",
            "ymd HMS", "ymd HM", "ymd",
            "mdy HMS", "mdy HM", "mdy"
          ),
          tz = LOCAL_TZ
        ),
        tz = LOCAL_TZ
      )
    )
  }
  
  # Excel serial datetime 備援；使用台灣當地時間原點，避免 UTC 位移。
  x_num <- suppressWarnings(as.numeric(x_chr))
  is_excel_serial <- (
    is.na(out) &
      !is.na(x_num) &
      x_num >= 20000 &
      x_num <= 80000 &
      !str_detect(x_chr, "[-/:]")
  )
  
  if (any(is_excel_serial)) {
    local_origin <- as.POSIXct(
      "1899-12-30 00:00:00",
      format = "%Y-%m-%d %H:%M:%S",
      tz = LOCAL_TZ
    )
    
    out[is_excel_serial] <- (
      local_origin +
        x_num[is_excel_serial] * 86400
    )
  }
  
  out
}

as_local_date <- function(x) {
  as.Date(x, tz = LOCAL_TZ)
}

format_local_datetime <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  
  out[ok] <- format(
    x[ok],
    format = "%Y-%m-%d %H:%M:%S",
    tz = LOCAL_TZ
  )
  
  out
}

make_local_day_end <- function(x) {
  as.POSIXct(
    paste0(format(x, "%Y-%m-%d"), " 23:59:59"),
    format = "%Y-%m-%d %H:%M:%S",
    tz = LOCAL_TZ
  )
}

resolve_character_mode <- function(x) {
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  tab <- table(x)
  candidates <- names(tab)[tab == max(tab)]
  sort(candidates)[1L]
}

resolve_date_mode <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(as.Date(NA))
  }
  
  x_chr <- as.character(x)
  tab <- table(x_chr)
  candidates <- names(tab)[tab == max(tab)]
  as.Date(sort(candidates)[1L])
}

clean_numeric_value <- function(x) {
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("", "NA", "N/A", "NaN", "NULL", "-", "--")] <- NA_character_
  x_chr <- str_replace_all(x_chr, ",", "")
  
  out <- suppressWarnings(as.numeric(x_chr))
  
  is_lt <- (
    !is.na(x_chr) &
      str_detect(x_chr, "^<\\s*[-+]?[0-9]*\\.?[0-9]+$")
  )
  
  if (any(is_lt)) {
    lt_value <- suppressWarnings(
      as.numeric(str_remove(x_chr[is_lt], "^<\\s*"))
    )
    out[is_lt] <- lt_value / 2
  }
  
  is_gt <- (
    !is.na(x_chr) &
      str_detect(x_chr, "^>\\s*[-+]?[0-9]*\\.?[0-9]+$")
  )
  
  if (any(is_gt)) {
    gt_value <- suppressWarnings(
      as.numeric(str_remove(x_chr[is_gt], "^>\\s*"))
    )
    out[is_gt] <- gt_value
  }
  
  out
}

is_left_censored <- function(x) {
  x_chr <- str_squish(as.character(x))
  
  !is.na(x_chr) &
    str_detect(x_chr, "^<\\s*[-+]?[0-9]*\\.?[0-9]+$")
}

is_right_censored <- function(x) {
  x_chr <- str_squish(as.character(x))
  
  !is.na(x_chr) &
    str_detect(x_chr, "^>\\s*[-+]?[0-9]*\\.?[0-9]+$")
}

check_non_numeric_vector <- function(x, variable_name) {
  x_chr <- str_squish(as.character(x))
  valid_pattern <- "^\\s*[<>]?\\s*[-+]?[0-9]*\\.?[0-9]+\\s*$"
  
  tibble(
    variable = variable_name,
    value = x_chr
  ) %>%
    filter(
      !is.na(value),
      value != "",
      !str_detect(value, valid_pattern)
    ) %>%
    count(variable, value, name = "n") %>%
    arrange(desc(n), value)
}

aggregate_numeric_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  median(x)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  mean(x)
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    NA_real_
  } else {
    min(x)
  }
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    NA_real_
  } else {
    max(x)
  }
}

format_value_list <- function(x) {
  x <- sort(unique(as.numeric(x)))
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  paste(
    format(x, trim = TRUE, scientific = FALSE),
    collapse = " | "
  )
}

format_pair_list <- function(height, weight) {
  ok <- is.finite(height) & is.finite(weight)
  
  if (!any(ok)) {
    return(NA_character_)
  }
  
  vals <- paste0(
    "H=", format(height[ok], trim = TRUE, scientific = FALSE),
    "; W=", format(weight[ok], trim = TRUE, scientific = FALSE)
  )
  
  paste(sort(unique(vals)), collapse = " | ")
}

# ============================================================
# 2. 疾病診斷資料：cohort 與糖尿病類型
# ============================================================

diag_path <- file.path(raw_root, "病人疾病診斷明細.xlsx")
dat_diag <- read_xlsx_hospital(diag_path)

assert_has_cols(
  dat_diag,
  c("病歷號", "就診號", "慢性病碼"),
  "疾病診斷資料"
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
  names(dat_diag)
)[1L]

if (is.na(diag_code_col)) {
  stop("疾病診斷資料中找不到 ICD-10 診斷碼欄位。")
}

dat_diag <- dat_diag %>%
  mutate(
    病歷號 = standardize_id(病歷號),
    就診號 = standardize_id(就診號),
    慢性病碼 = str_squish(as.character(慢性病碼)),
    icd_code_raw = str_squish(as.character(.data[[diag_code_col]])),
    icd_code = str_to_upper(
      str_replace_all(icd_code_raw, "\\.", "")
    )
  )

# cohort 納入仍依醫院慢性病碼 01 的糖尿病相關就診。
dat_diabetes_visits <- dat_diag %>%
  filter(慢性病碼 %in% c("01", "1")) %>%
  select(病歷號, 就診號) %>%
  filter(!is.na(病歷號), !is.na(就診號)) %>%
  distinct()

cohort_patients <- dat_diabetes_visits %>%
  distinct(病歷號)

# 類型判定使用 cohort 病患所有可得 E10/E11 診斷紀錄。
# E08、E09、E13 僅計入其他糖尿病碼，不混入一般非糖尿病診斷。
dat_dm_type <- dat_diag %>%
  semi_join(cohort_patients, by = "病歷號") %>%
  mutate(
    dm_type_code = case_when(
      str_detect(icd_code, "^E10") ~ "Type 1",
      str_detect(icd_code, "^E11") ~ "Type 2",
      str_detect(icd_code, "^(E08|E09|E13)") ~ "Other diabetes code",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(病歷號) %>%
  summarise(
    n_E10 = sum(dm_type_code == "Type 1", na.rm = TRUE),
    n_E11 = sum(dm_type_code == "Type 2", na.rm = TRUE),
    n_other_diabetes_code = sum(
      dm_type_code == "Other diabetes code",
      na.rm = TRUE
    ),
    diabetes_type = case_when(
      n_E10 > 0 & n_E11 == 0 ~ "Type 1",
      n_E11 > 0 & n_E10 == 0 ~ "Type 2",
      n_E10 > 0 & n_E11 > 0 ~ "Mixed E10 and E11",
      TRUE ~ "Other or unclear"
    ),
    .groups = "drop"
  )

# ============================================================
# 3. 病人基本資料：性別與真實身高體重配對
# ============================================================

basic_path <- file.path(
  raw_root,
  "病人基本資料：年齡,性別,身高,體重.xlsx"
)

dat_basic_raw <- read_xlsx_hospital(basic_path)

assert_has_cols(
  dat_basic_raw,
  c("病歷號", "性別", "測量時間", "身高", "體重"),
  "病人基本資料"
)

dat_basic <- dat_basic_raw %>%
  mutate(
    .basic_row_id = row_number(),
    病歷號 = standardize_id(病歷號),
    sex_label_basic_raw = standardize_sex(性別),
    measurement_time_raw = str_squish(as.character(測量時間)),
    測量日期時間 = parse_datetime_safe(測量時間),
    身高_raw = clean_numeric_value(身高),
    體重_raw = clean_numeric_value(體重)
  ) %>%
  filter(!is.na(病歷號))

basic_time_parse_failure_qc <- dat_basic %>%
  filter(
    !is.na(measurement_time_raw),
    measurement_time_raw != "",
    is.na(測量日期時間)
  ) %>%
  select(
    病歷號,
    .basic_row_id,
    measurement_time_raw
  )

# 已確認之特定輸入錯誤保留明確修正紀錄。
dat_basic <- dat_basic %>%
  mutate(
    height_manual_correction_flag = (
      !is.na(身高_raw) &
        身高_raw == 1181.6
    ),
    weight_manual_correction_flag = (
      !is.na(體重_raw) &
        病歷號 == "000054812J" &
        體重_raw == 956
    ),
    身高_corrected = case_when(
      height_manual_correction_flag ~ 181.6,
      TRUE ~ 身高_raw
    ),
    體重_corrected = case_when(
      weight_manual_correction_flag ~ 95.6,
      TRUE ~ 體重_raw
    ),
    height_abs_invalid_flag = (
      !is.na(身高_corrected) &
        (
          身高_corrected < ABS_HEIGHT_RANGE_CM[1] |
            身高_corrected > ABS_HEIGHT_RANGE_CM[2]
        )
    ),
    weight_abs_invalid_flag = (
      !is.na(體重_corrected) &
        (
          體重_corrected < ABS_WEIGHT_RANGE_KG[1] |
            體重_corrected > ABS_WEIGHT_RANGE_KG[2]
        )
    ),
    身高 = if_else(
      height_abs_invalid_flag,
      NA_real_,
      身高_corrected
    ),
    體重 = if_else(
      weight_abs_invalid_flag,
      NA_real_,
      體重_corrected
    ),
    BMI_raw_pair = if_else(
      !is.na(身高) &
        身高 > 0 &
        !is.na(體重),
      體重 / (身高 / 100)^2,
      NA_real_
    ),
    bmi_abs_invalid_flag = (
      !is.na(BMI_raw_pair) &
        (
          BMI_raw_pair < ABS_BMI_RANGE[1] |
            BMI_raw_pair > ABS_BMI_RANGE[2]
        )
    ),
    BMI = if_else(
      bmi_abs_invalid_flag,
      NA_real_,
      BMI_raw_pair
    )
  )

basic_invalid_anthro_qc <- dat_basic %>%
  filter(
    height_abs_invalid_flag |
      weight_abs_invalid_flag |
      bmi_abs_invalid_flag |
      height_manual_correction_flag |
      weight_manual_correction_flag
  ) %>%
  select(
    病歷號,
    測量時間,
    測量日期時間,
    身高_raw,
    體重_raw,
    身高_corrected,
    體重_corrected,
    BMI_raw_pair,
    height_manual_correction_flag,
    weight_manual_correction_flag,
    height_abs_invalid_flag,
    weight_abs_invalid_flag,
    bmi_abs_invalid_flag
  ) %>%
  arrange(病歷號, 測量日期時間)

# 病患層級性別以基本資料眾數為主；若票數同分，採字典順序固定決定，
# 並保留衝突 QC，不依原始列順序決定。
basic_sex_reference <- dat_basic %>%
  group_by(病歷號) %>%
  summarise(
    sex_label_basic = resolve_character_mode(sex_label_basic_raw),
    basic_sex_n_distinct = n_distinct(
      sex_label_basic_raw[!is.na(sex_label_basic_raw)]
    ),
    basic_sex_values = paste(
      sort(unique(sex_label_basic_raw[!is.na(sex_label_basic_raw)])),
      collapse = " | "
    ),
    .groups = "drop"
  )

basic_sex_conflict_qc <- basic_sex_reference %>%
  filter(basic_sex_n_distinct > 1L)

# 身高與體重必須來自同一原始列，且該配對 BMI 通過寬鬆技術性檢查。
basic_pair_candidates <- dat_basic %>%
  filter(
    !is.na(測量日期時間),
    is.finite(身高),
    is.finite(體重),
    is.finite(BMI)
  ) %>%
  transmute(
    病歷號,
    測量日期時間,
    身高,
    體重,
    BMI,
    .basic_row_id
  )

# 同一病患、完全相同時間若出現不同身高體重配對，不應自動平均。
basic_timestamp_pair_qc <- basic_pair_candidates %>%
  group_by(病歷號, 測量日期時間) %>%
  summarise(
    n_source_records = n(),
    n_distinct_pair = n_distinct(
      paste(身高, 體重, sep = "|")
    ),
    pair_values = format_pair_list(身高, 體重),
    .groups = "drop"
  ) %>%
  filter(n_distinct_pair > 1L)

if (nrow(basic_timestamp_pair_qc) > 0L) {
  write_xlsx(
    list(
      timestamp_pair_conflict = basic_timestamp_pair_qc,
      invalid_anthropometry = basic_invalid_anthro_qc
    ),
    file.path(
      save_folder,
      "basic_pair_conflicts_REVIEW_REQUIRED.xlsx"
    )
  )
  
  stop(
    "同一病患、相同測量日期時間存在不同身高體重配對；",
    "已輸出 basic_pair_conflicts_REVIEW_REQUIRED.xlsx，",
    "請人工確認後再執行。"
  )
}

# 同時間完全相同之配對僅保留一次；不以日為單位彙整。
dat_basic_pairs <- basic_pair_candidates %>%
  arrange(病歷號, 測量日期時間, .basic_row_id) %>%
  distinct(
    病歷號,
    測量日期時間,
    身高,
    體重,
    BMI,
    .keep_all = TRUE
  ) %>%
  group_by(病歷號, 測量日期時間) %>%
  slice(1L) %>%
  ungroup() %>%
  select(-.basic_row_id)

# 同日不同時間出現不同配對只作 QA；不再合併成中位數。
basic_same_date_qc <- dat_basic_pairs %>%
  mutate(測量日期 = as_local_date(測量日期時間)) %>%
  group_by(病歷號, 測量日期) %>%
  summarise(
    n_measurement_time = n_distinct(測量日期時間),
    n_distinct_pair = n_distinct(
      paste(身高, 體重, sep = "|")
    ),
    pair_values = format_pair_list(身高, 體重),
    .groups = "drop"
  ) %>%
  filter(n_distinct_pair > 1L)

basic_pair_summary <- tibble(
  item = c(
    "basic_raw_row_n",
    "basic_patient_n",
    "valid_paired_measurement_n",
    "timestamp_pair_conflict_n",
    "same_date_different_pair_n",
    "invalid_or_corrected_row_n",
    "measurement_time_parse_failure_n"
  ),
  value = c(
    nrow(dat_basic),
    n_distinct(dat_basic$病歷號),
    nrow(dat_basic_pairs),
    nrow(basic_timestamp_pair_qc),
    nrow(basic_same_date_qc),
    nrow(basic_invalid_anthro_qc),
    nrow(basic_time_parse_failure_qc)
  )
)


# 將「無效或經修正」拆成可解釋的原因，避免將所有列誤解為人工修正。
basic_anthro_flag_summary <- tibble(
  item = c(
    "height_manual_correction_n",
    "weight_manual_correction_n",
    "height_abs_invalid_n",
    "weight_abs_invalid_n",
    "bmi_abs_invalid_n",
    "height_weight_both_zero_n",
    "any_invalid_or_corrected_row_n"
  ),
  value = c(
    sum(
      dat_basic$height_manual_correction_flag,
      na.rm = TRUE
    ),
    sum(
      dat_basic$weight_manual_correction_flag,
      na.rm = TRUE
    ),
    sum(
      dat_basic$height_abs_invalid_flag,
      na.rm = TRUE
    ),
    sum(
      dat_basic$weight_abs_invalid_flag,
      na.rm = TRUE
    ),
    sum(
      dat_basic$bmi_abs_invalid_flag,
      na.rm = TRUE
    ),
    sum(
      dat_basic$身高_raw == 0 &
        dat_basic$體重_raw == 0,
      na.rm = TRUE
    ),
    nrow(basic_invalid_anthro_qc)
  )
)

# ============================================================
# 4. 血液檢驗資料
# ============================================================

folders <- list.files(lab_root, full.names = TRUE)
folder_lookup <- setNames(folders, basename(folders))

missing_lab_folders <- setdiff(
  names(lab_map),
  names(folder_lookup)
)

if (length(missing_lab_folders) > 0L) {
  warning(
    "下列指定檢驗資料夾不存在：",
    paste(missing_lab_folders, collapse = ", ")
  )
}

standardize_lab_keys <- function(dat, data_name = "血液檢驗資料") {
  assert_has_cols(
    dat,
    c(
      "病歷號",
      "就診號",
      "報告日期",
      "出生年月日",
      "性別"
    ),
    data_name
  )
  
  dat %>%
    mutate(
      病歷號 = standardize_id(病歷號),
      就診號 = standardize_id(就診號),
      報告日期 = parse_date_safe(報告日期),
      出生年月日 = parse_date_safe(出生年月日),
      sex_label_lab_raw = standardize_sex(性別)
    )
}

lab_event_list <- list()
lab_conflict_list <- list()
lab_ingestion_qc_list <- list()
non_numeric_lab_list <- list()

for (folder_name in names(lab_map)) {
  mapped_name <- unname(lab_map[folder_name])
  folder_path <- unname(folder_lookup[folder_name])
  
  if (is.na(folder_path) || !dir.exists(folder_path)) {
    next
  }
  
  message("正在讀取：", folder_name, " -> ", mapped_name)
  
  files <- list.files(
    folder_path,
    full.names = TRUE,
    recursive = TRUE,
    pattern = "\\.xlsx?$"
  )
  
  if (length(files) == 0L) {
    warning("資料夾中沒有 xlsx/xls 檔案：", folder_path)
    next
  }
  
  dat_lab_raw <- map_dfr(
    files,
    ~ read_xlsx_hospital(.x) %>%
      mutate(.source_file = basename(.x))
  )
  
  assert_has_cols(
    dat_lab_raw,
    c("檢驗結果", EVENT_KEYS),
    paste0("檢驗資料：", folder_name)
  )
  
  dat_lab_raw <- standardize_lab_keys(
    dat_lab_raw,
    paste0("檢驗資料：", folder_name)
  ) %>%
    mutate(
      raw_value = str_squish(as.character(檢驗結果)),
      numeric_value = clean_numeric_value(raw_value),
      left_censored = is_left_censored(raw_value),
      right_censored = is_right_censored(raw_value)
    ) %>%
    filter(
      !is.na(病歷號),
      !is.na(就診號),
      !is.na(報告日期)
    )
  
  non_numeric_lab_list[[mapped_name]] <- check_non_numeric_vector(
    dat_lab_raw$raw_value,
    mapped_name
  )
  
  n_raw <- nrow(dat_lab_raw)
  
  # 相同事件、相同原始結果之重複列只保留一次。
  dat_lab_dedup <- dat_lab_raw %>%
    distinct(
      across(all_of(EVENT_KEYS)),
      raw_value,
      .keep_all = TRUE
    )
  
  conflict_tbl <- dat_lab_dedup %>%
    group_by(across(all_of(EVENT_KEYS))) %>%
    summarise(
      lab_variable = mapped_name,
      n_source_records = n(),
      n_distinct_numeric = n_distinct(
        numeric_value[is.finite(numeric_value)]
      ),
      min_value = safe_min(numeric_value),
      max_value = safe_max(numeric_value),
      value_range = max_value - min_value,
      values = format_value_list(numeric_value),
      any_left_censored = any(left_censored, na.rm = TRUE),
      any_right_censored = any(right_censored, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_distinct_numeric > 1L)
  
  lab_conflict_list[[mapped_name]] <- conflict_tbl
  
  event_tbl <- dat_lab_dedup %>%
    group_by(across(all_of(EVENT_KEYS))) %>%
    summarise(
      value = aggregate_numeric_median(numeric_value),
      n_source_records = n(),
      n_distinct_numeric = n_distinct(
        numeric_value[is.finite(numeric_value)]
      ),
      left_censored_flag = any(left_censored, na.rm = TRUE),
      right_censored_flag = any(right_censored, na.rm = TRUE),
      出生年月日 = resolve_date_mode(出生年月日),
      sex_label_lab = resolve_character_mode(sex_label_lab_raw),
      birthdate_n_distinct = n_distinct(
        出生年月日[!is.na(出生年月日)]
      ),
      lab_sex_n_distinct = n_distinct(
        sex_label_lab_raw[!is.na(sex_label_lab_raw)]
      ),
      .groups = "drop"
    ) %>%
    rename(
      !!mapped_name := value,
      !!paste0(mapped_name, "_n_source") := n_source_records,
      !!paste0(mapped_name, "_n_distinct") := n_distinct_numeric,
      !!paste0(mapped_name, "_left_censored") := left_censored_flag,
      !!paste0(mapped_name, "_right_censored") := right_censored_flag
    )
  
  if (mapped_name != "HbA1C") {
    event_tbl <- event_tbl %>%
      select(
        all_of(EVENT_KEYS),
        all_of(mapped_name),
        all_of(paste0(mapped_name, "_n_source")),
        all_of(paste0(mapped_name, "_n_distinct")),
        all_of(paste0(mapped_name, "_left_censored")),
        all_of(paste0(mapped_name, "_right_censored"))
      )
  }
  
  lab_event_list[[mapped_name]] <- event_tbl
  
  lab_ingestion_qc_list[[mapped_name]] <- tibble(
    lab_variable = mapped_name,
    raw_row_n = n_raw,
    after_exact_dedup_n = nrow(dat_lab_dedup),
    event_n = nrow(event_tbl),
    conflict_event_n = nrow(conflict_tbl),
    missing_numeric_n = sum(is.na(dat_lab_dedup$numeric_value)),
    left_censored_record_n = sum(
      dat_lab_dedup$left_censored,
      na.rm = TRUE
    ),
    right_censored_record_n = sum(
      dat_lab_dedup$right_censored,
      na.rm = TRUE
    )
  )
}

missing_main_labs <- setdiff(
  MAIN_LABS,
  names(lab_event_list)
)

if (length(missing_main_labs) > 0L) {
  stop(
    "下列主要檢驗未成功讀入：",
    paste(missing_main_labs, collapse = ", ")
  )
}

lab_ingestion_qc <- bind_rows(lab_ingestion_qc_list)
lab_value_conflicts <- bind_rows(lab_conflict_list)
non_numeric_before <- bind_rows(non_numeric_lab_list)

if (ncol(lab_value_conflicts) == 0L) {
  lab_value_conflicts <- tibble(
    病歷號 = character(),
    就診號 = character(),
    報告日期 = as.Date(character()),
    lab_variable = character(),
    n_source_records = integer(),
    n_distinct_numeric = integer(),
    min_value = double(),
    max_value = double(),
    value_range = double(),
    values = character(),
    any_left_censored = logical(),
    any_right_censored = logical()
  )
}

if (ncol(non_numeric_before) == 0L) {
  non_numeric_before <- tibble(
    variable = character(),
    value = character(),
    n = integer()
  )
}

# 以 HbA1c 建立一列一事件之 anchor，再串接其他檢驗。
dat_blood <- lab_event_list[["HbA1C"]]

other_labs <- setdiff(
  names(lab_event_list),
  "HbA1C"
)

for (lab_name in other_labs) {
  dat_blood <- dat_blood %>%
    left_join(
      lab_event_list[[lab_name]],
      by = EVENT_KEYS,
      relationship = "one-to-one"
    )
}

duplicate_anchor_event <- dat_blood %>%
  count(
    across(all_of(EVENT_KEYS)),
    name = "n"
  ) %>%
  filter(n > 1L)

if (nrow(duplicate_anchor_event) > 0L) {
  stop("HbA1c anchor 事件仍有重複。")
}

# 病患層級出生日期只用於計算年齡；精確日期不輸出至正式主檔。
birthdate_reference <- dat_blood %>%
  group_by(病歷號) %>%
  summarise(
    出生年月日_ref = resolve_date_mode(出生年月日),
    birthdate_patient_n_distinct = n_distinct(
      出生年月日[!is.na(出生年月日)]
    ),
    birthdate_values = paste(
      sort(unique(as.character(出生年月日[!is.na(出生年月日)]))),
      collapse = " | "
    ),
    .groups = "drop"
  )

birthdate_conflict_qc <- birthdate_reference %>%
  filter(birthdate_patient_n_distinct > 1L)

dat_blood <- dat_blood %>%
  select(-出生年月日) %>%
  left_join(
    birthdate_reference,
    by = "病歷號",
    relationship = "many-to-one"
  ) %>%
  mutate(
    Age = floor(
      time_length(
        interval(出生年月日_ref, 報告日期),
        unit = "years"
      )
    ),
    age_out_of_range_flag = (
      !is.na(Age) &
        (Age < 0 | Age > 110)
    ),
    Age = if_else(
      age_out_of_range_flag,
      NA_real_,
      as.numeric(Age)
    )
  )

# ============================================================
# 5. 合併糖尿病 cohort 與 HbA1c 事件
# ============================================================

stage_counts <- tibble(
  stage = c(
    "糖尿病 cohort 病患數",
    "糖尿病 cohort 就診數",
    "全部 HbA1c anchor 事件數"
  ),
  patient_n = c(
    n_distinct(dat_diabetes_visits$病歷號),
    n_distinct(dat_diabetes_visits$病歷號),
    n_distinct(dat_blood$病歷號)
  ),
  observation_n = c(
    NA_integer_,
    nrow(dat_diabetes_visits),
    nrow(dat_blood)
  )
)

dat_diabetes_2 <- dat_diabetes_visits %>%
  inner_join(
    dat_blood,
    by = c("病歷號", "就診號"),
    relationship = "one-to-many"
  ) %>%
  filter(!is.na(HbA1C)) %>%
  mutate(lab_row_id = row_number())

stage_counts <- bind_rows(
  stage_counts,
  tibble(
    stage = "糖尿病 cohort 內 HbA1c 追蹤事件",
    patient_n = n_distinct(dat_diabetes_2$病歷號),
    observation_n = nrow(dat_diabetes_2)
  )
)

duplicate_cohort_event <- dat_diabetes_2 %>%
  count(
    across(all_of(EVENT_KEYS)),
    name = "n"
  ) %>%
  filter(n > 1L)

if (nrow(duplicate_cohort_event) > 0L) {
  stop("cohort 合併後同一檢驗事件重複。")
}

# ============================================================
# 6. 性別與最近一次真實身高體重配對
# ============================================================

dat_diabetes_2 <- dat_diabetes_2 %>%
  left_join(
    basic_sex_reference,
    by = "病歷號",
    relationship = "many-to-one"
  ) %>%
  mutate(
    sex_conflict_flag = (
      !is.na(sex_label_basic) &
        !is.na(sex_label_lab) &
        sex_label_basic != sex_label_lab
    ),
    sex_label = coalesce(
      sex_label_basic,
      sex_label_lab
    ),
    sex_source = case_when(
      !is.na(sex_label_basic) ~ "basic_file",
      is.na(sex_label_basic) & !is.na(sex_label_lab) ~ "lab_fallback",
      TRUE ~ "missing"
    ),
    性別 = case_when(
      sex_label == "男" ~ 1,
      sex_label == "女" ~ 0,
      TRUE ~ NA_real_
    ),
    # 檢驗資料只有日期，故以台灣當地時間 23:59:59 作為當日截止點。
    # 不使用 Date + days(1) 再轉 POSIXct，避免時區造成次日清晨資料被誤納。
    lab_event_cutoff = make_local_day_end(報告日期)
  )

sex_conflict_qc <- dat_diabetes_2 %>%
  filter(sex_conflict_flag) %>%
  select(
    病歷號,
    就診號,
    報告日期,
    sex_label_basic,
    sex_label_lab,
    sex_source
  ) %>%
  distinct()

# 僅使用同一原始列中同時存在且通過技術性 QA 的身高、體重與 BMI。
dat_basic_for_join <- dat_basic_pairs %>%
  semi_join(
    dat_diabetes_2 %>% distinct(病歷號),
    by = "病歷號"
  )

dat_diabetes_3 <- dat_diabetes_2 %>%
  left_join(
    dat_basic_for_join,
    by = "病歷號",
    relationship = "many-to-many"
  ) %>%
  mutate(
    valid_basic_pair = (
      !is.na(測量日期時間) &
        測量日期時間 <= lab_event_cutoff
    ),
    basic_gap_seconds = if_else(
      valid_basic_pair,
      as.numeric(
        difftime(
          lab_event_cutoff,
          測量日期時間,
          units = "secs"
        )
      ),
      Inf
    )
  )

dat_diabetes_4 <- dat_diabetes_3 %>%
  group_by(lab_row_id) %>%
  arrange(
    desc(valid_basic_pair),
    basic_gap_seconds,
    desc(測量日期時間),
    .by_group = TRUE
  ) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(
    測量日期時間 = if_else(
      valid_basic_pair,
      測量日期時間,
      as.POSIXct(NA, tz = LOCAL_TZ)
    ),
    身高 = if_else(
      valid_basic_pair,
      as.numeric(身高),
      NA_real_
    ),
    體重 = if_else(
      valid_basic_pair,
      as.numeric(體重),
      NA_real_
    ),
    BMI = if_else(
      valid_basic_pair,
      as.numeric(BMI),
      NA_real_
    ),
    測量日期 = if_else(
      valid_basic_pair,
      as_local_date(測量日期時間),
      as.Date(NA)
    ),
    basic_time_gap_days = if_else(
      valid_basic_pair,
      as.numeric(
        報告日期 - 測量日期
      ),
      NA_real_
    ),
    basic_gap_gt365_flag = (
      !is.na(basic_time_gap_days) &
        basic_time_gap_days > 365
    ),
    basic_gap_gt730_flag = (
      !is.na(basic_time_gap_days) &
        basic_time_gap_days > 730
    ),
    basic_match_source = case_when(
      is.na(測量日期時間) ~ "no_prior_valid_pair",
      basic_time_gap_days == 0 ~ "same_day_pair",
      basic_time_gap_days > 0 ~ "prior_pair",
      TRUE ~ "unclassified"
    )
  ) %>%
  select(
    -any_of(
      c(
        "lab_event_cutoff",
        "valid_basic_pair",
        "basic_gap_seconds"
      )
    )
  )

if (nrow(dat_diabetes_4) != nrow(dat_diabetes_2)) {
  stop("基本資料串接後事件筆數改變。")
}

if (n_distinct(dat_diabetes_4$lab_row_id) != nrow(dat_diabetes_4)) {
  stop("基本資料串接後 lab_row_id 不唯一。")
}

selected_basic_after_report_qc <- dat_diabetes_4 %>%
  filter(
    !is.na(測量日期),
    測量日期 > 報告日期
  ) %>%
  select(
    病歷號,
    就診號,
    報告日期,
    測量日期,
    測量日期時間,
    身高,
    體重,
    BMI
  )

if (any(dat_diabetes_4$basic_time_gap_days < 0, na.rm = TRUE)) {
  stop("發現負的基本資料時間差。")
}

if (nrow(selected_basic_after_report_qc) > 0L) {
  stop("發現檢驗報告日之後的基本資料被誤選。")
}

# ============================================================
# 7. 糖尿病類型、事件唯一性與 Repeat
# ============================================================

dat_diabetes_4 <- dat_diabetes_4 %>%
  left_join(
    dat_dm_type,
    by = "病歷號",
    relationship = "many-to-one"
  )

duplicate_final_event <- dat_diabetes_4 %>%
  count(
    across(all_of(EVENT_KEYS)),
    name = "n"
  ) %>%
  filter(n > 1L)

if (nrow(duplicate_final_event) > 0L) {
  stop("最終資料仍有重複檢驗事件。")
}

dat_diabetes_4 <- dat_diabetes_4 %>%
  arrange(
    病歷號,
    報告日期,
    就診號,
    lab_row_id
  ) %>%
  group_by(病歷號) %>%
  mutate(Repeat = row_number()) %>%
  ungroup()

# ============================================================
# 8. 建立正式主檔與 QA
# ============================================================

# 正式主檔僅保留下游分析需要或需追溯之欄位。
# 精確出生日期、來源筆數與衝突細節留在 QC，不放入分析主檔。
analysis_main <- dat_diabetes_4 %>%
  select(
    病歷號,
    就診號,
    報告日期,
    Repeat,
    diabetes_type,
    n_E10,
    n_E11,
    n_other_diabetes_code,
    性別,
    sex_label,
    Age,
    測量日期,
    測量日期時間,
    basic_time_gap_days,
    basic_gap_gt365_flag,
    basic_gap_gt730_flag,
    身高,
    體重,
    BMI,
    HbA1C,
    Glucose_AC,
    `LDL-C`,
    Triglycerides,
    `HDL-C`,
    Total_cholesterol,
    HbA1C_left_censored,
    HbA1C_right_censored,
    Glucose_AC_left_censored,
    Glucose_AC_right_censored,
    `LDL-C_left_censored`,
    `LDL-C_right_censored`,
    Triglycerides_left_censored,
    Triglycerides_right_censored,
    `HDL-C_left_censored`,
    `HDL-C_right_censored`,
    Total_cholesterol_left_censored,
    Total_cholesterol_right_censored
  ) %>%
  mutate(
    # 正式主分析是否符合第二型糖尿病分類。
    # 後續仍須另依表格.R 的成人與其他納入排除條件篩選。
    analysis_type2_eligible = (
      diabetes_type == "Type 2"
    ),
    BMI_source_age_category = case_when(
      is.na(basic_time_gap_days) ~
        "No valid prior measurement",
      basic_time_gap_days == 0 ~
        "Same day",
      basic_time_gap_days <= 365 ~
        "1-365 days",
      basic_time_gap_days <= 730 ~
        "366-730 days",
      TRUE ~
        ">730 days"
    ),
    # 輸出為明確的台灣當地時間字串，避免 xlsx/csv 儲存時轉成 UTC。
    測量日期時間 = format_local_datetime(測量日期時間)
  )

unique_event_n <- analysis_main %>%
  distinct(across(all_of(EVENT_KEYS))) %>%
  nrow()

final_event_keys <- analysis_main %>%
  distinct(
    across(all_of(EVENT_KEYS))
  ) %>%
  mutate(
    in_final_cohort = TRUE
  )

lab_value_conflicts <- lab_value_conflicts %>%
  left_join(
    final_event_keys,
    by = EVENT_KEYS,
    relationship = "many-to-one"
  ) %>%
  mutate(
    in_final_cohort = coalesce(
      in_final_cohort,
      FALSE
    )
  ) %>%
  arrange(
    desc(in_final_cohort),
    lab_variable,
    病歷號,
    報告日期,
    就診號
  )

lab_conflict_final_cohort_summary <- lab_value_conflicts %>%
  count(
    lab_variable,
    in_final_cohort,
    name = "conflict_event_n"
  ) %>%
  arrange(
    lab_variable,
    desc(in_final_cohort)
  )

overview <- tibble(
  item = c(
    "final_patient_n",
    "final_observation_n",
    "unique_event_n",
    "sex_conflict_event_n",
    "basic_sex_conflict_patient_n",
    "birthdate_conflict_patient_n",
    "lab_conflict_event_n",
    "basic_timestamp_pair_conflict_n",
    "basic_same_date_different_pair_n",
    "basic_invalid_or_corrected_row_n",
    "measurement_time_parse_failure_n",
    "selected_basic_after_report_n"
  ),
  value = c(
    n_distinct(analysis_main$病歷號),
    nrow(analysis_main),
    unique_event_n,
    nrow(sex_conflict_qc),
    nrow(basic_sex_conflict_qc),
    nrow(birthdate_conflict_qc),
    nrow(lab_value_conflicts),
    nrow(basic_timestamp_pair_qc),
    nrow(basic_same_date_qc),
    nrow(basic_invalid_anthro_qc),
    nrow(basic_time_parse_failure_qc),
    nrow(selected_basic_after_report_qc)
  )
)

stage_counts <- bind_rows(
  stage_counts,
  tibble(
    stage = "最終串接資料",
    patient_n = n_distinct(analysis_main$病歷號),
    observation_n = nrow(analysis_main)
  )
)

repeat_patient_summary <- analysis_main %>%
  count(病歷號, name = "n_obs") %>%
  summarise(
    patient_n = n(),
    mean_n_obs = mean(n_obs),
    sd_n_obs = sd(n_obs),
    median_n_obs = median(n_obs),
    q1_n_obs = quantile(n_obs, 0.25),
    q3_n_obs = quantile(n_obs, 0.75),
    min_n_obs = min(n_obs),
    max_n_obs = max(n_obs),
    .groups = "drop"
  )

repeat_count_by_order <- analysis_main %>%
  count(Repeat, name = "n_records") %>%
  arrange(Repeat)

missing_vars <- c(
  "HbA1C",
  "Glucose_AC",
  "LDL-C",
  "Triglycerides",
  "HDL-C",
  "Total_cholesterol",
  "BMI",
  "Age",
  "性別",
  "diabetes_type"
)

missing_summary <- map_dfr(
  missing_vars,
  function(v) {
    tibble(
      variable = v,
      total_n = nrow(analysis_main),
      missing_n = sum(is.na(analysis_main[[v]])),
      observed_n = sum(!is.na(analysis_main[[v]])),
      missing_rate_percent = round(
        mean(is.na(analysis_main[[v]])) * 100,
        2
      )
    )
  }
)

patient_diabetes_type_summary <- dat_dm_type %>%
  semi_join(
    analysis_main %>% distinct(病歷號),
    by = "病歷號"
  ) %>%
  arrange(diabetes_type, 病歷號)

dm_type_patient_count <- patient_diabetes_type_summary %>%
  count(diabetes_type, name = "patient_n") %>%
  arrange(diabetes_type)

dm_type_observation_count <- analysis_main %>%
  count(diabetes_type, name = "observation_n") %>%
  arrange(diabetes_type)

basic_match_summary <- dat_diabetes_4 %>%
  count(basic_match_source, name = "observation_n") %>%
  mutate(
    observation_percent = round(
      observation_n / sum(observation_n) * 100,
      2
    )
  ) %>%
  arrange(basic_match_source)

basic_gap_summary <- dat_diabetes_4 %>%
  summarise(
    total_observation_n = n(),
    matched_observation_n = sum(
      !is.na(basic_time_gap_days)
    ),
    unmatched_observation_n = sum(
      is.na(basic_time_gap_days)
    ),
    matched_observation_percent = round(
      mean(!is.na(basic_time_gap_days)) * 100,
      2
    ),
    unmatched_observation_percent = round(
      mean(is.na(basic_time_gap_days)) * 100,
      2
    ),
    mean_gap_days = safe_mean(
      basic_time_gap_days
    ),
    sd_gap_days = sd(
      basic_time_gap_days,
      na.rm = TRUE
    ),
    median_gap_days = median(
      basic_time_gap_days,
      na.rm = TRUE
    ),
    q1_gap_days = as.numeric(
      quantile(
        basic_time_gap_days,
        0.25,
        na.rm = TRUE
      )
    ),
    q3_gap_days = as.numeric(
      quantile(
        basic_time_gap_days,
        0.75,
        na.rm = TRUE
      )
    ),
    max_gap_days = safe_max(
      basic_time_gap_days
    ),
    gap_gt365_n = sum(
      basic_gap_gt365_flag,
      na.rm = TRUE
    ),
    gap_gt365_percent_of_matched = round(
      if_else(
        matched_observation_n > 0,
        gap_gt365_n /
          matched_observation_n * 100,
        NA_real_
      ),
      2
    ),
    gap_gt730_n = sum(
      basic_gap_gt730_flag,
      na.rm = TRUE
    ),
    gap_gt730_percent_of_matched = round(
      if_else(
        matched_observation_n > 0,
        gap_gt730_n /
          matched_observation_n * 100,
        NA_real_
      ),
      2
    )
  )

basic_gap_patient_summary <- dat_diabetes_4 %>%
  group_by(病歷號) %>%
  summarise(
    any_valid_basic_pair = any(
      !is.na(basic_time_gap_days)
    ),
    any_gap_gt365 = any(
      basic_gap_gt365_flag,
      na.rm = TRUE
    ),
    any_gap_gt730 = any(
      basic_gap_gt730_flag,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  summarise(
    total_patient_n = n(),
    patient_with_valid_basic_pair_n = sum(
      any_valid_basic_pair
    ),
    patient_without_valid_basic_pair_n = sum(
      !any_valid_basic_pair
    ),
    gap_gt365_patient_n = sum(
      any_gap_gt365
    ),
    gap_gt730_patient_n = sum(
      any_gap_gt730
    ),
    gap_gt365_patient_percent = round(
      mean(any_gap_gt365) * 100,
      2
    ),
    gap_gt730_patient_percent = round(
      mean(any_gap_gt730) * 100,
      2
    )
  )

long_basic_gap_qc <- dat_diabetes_4 %>%
  filter(basic_gap_gt730_flag) %>%
  select(
    病歷號,
    就診號,
    報告日期,
    測量日期時間,
    basic_time_gap_days,
    身高,
    體重,
    BMI
  ) %>%
  arrange(desc(basic_time_gap_days), 病歷號, 報告日期)

same_patient_date_multiple_visit_keys <- analysis_main %>%
  count(
    病歷號,
    報告日期,
    name = "n_events"
  ) %>%
  filter(n_events > 1L)

same_patient_date_multiple_visit_qc <- analysis_main %>%
  semi_join(
    same_patient_date_multiple_visit_keys,
    by = c("病歷號", "報告日期")
  ) %>%
  arrange(
    病歷號,
    報告日期,
    就診號
  )

# QA 表內的 POSIXct 欄位也以台灣當地時間字串輸出。
basic_invalid_anthro_export <- basic_invalid_anthro_qc %>%
  mutate(
    測量日期時間 = format_local_datetime(測量日期時間)
  )

basic_timestamp_pair_export <- basic_timestamp_pair_qc %>%
  mutate(
    測量日期時間 = format_local_datetime(測量日期時間)
  )

long_basic_gap_export <- long_basic_gap_qc %>%
  mutate(
    測量日期時間 = format_local_datetime(測量日期時間)
  )

selected_basic_after_report_export <- selected_basic_after_report_qc %>%
  mutate(
    測量日期時間 = format_local_datetime(測量日期時間)
  )

sex_source_summary <- dat_diabetes_4 %>%
  count(sex_source, name = "observation_n") %>%
  arrange(sex_source)

event_source_summary <- dat_diabetes_4 %>%
  summarise(
    hbA1c_multi_source_event_n = sum(
      HbA1C_n_source > 1,
      na.rm = TRUE
    ),
    glu_ac_multi_source_event_n = sum(
      Glucose_AC_n_source > 1,
      na.rm = TRUE
    ),
    ldl_multi_source_event_n = sum(
      `LDL-C_n_source` > 1,
      na.rm = TRUE
    ),
    tg_multi_source_event_n = sum(
      Triglycerides_n_source > 1,
      na.rm = TRUE
    ),
    hbA1c_right_censored_event_n = sum(
      HbA1C_right_censored,
      na.rm = TRUE
    )
  )

repeat_sequence_qc <- analysis_main %>%
  group_by(病歷號) %>%
  summarise(
    n_obs = n(),
    repeat_min = min(Repeat),
    repeat_max = max(Repeat),
    repeat_unique_n = n_distinct(Repeat),
    repeat_sequence_ok = (
      repeat_min == 1L &
        repeat_max == n_obs &
        repeat_unique_n == n_obs
    ),
    .groups = "drop"
  ) %>%
  filter(!repeat_sequence_ok)

data_dictionary <- tribble(
  ~variable, ~description, ~role,
  "病歷號", "去識別化病患識別碼", "主鍵之一",
  "就診號", "去識別化就診識別碼", "事件鍵之一",
  "報告日期", "HbA1c 檢驗報告日期", "事件鍵與時間",
  "Repeat", "病患內依日期排序之觀測序號", "縱向索引",
  "diabetes_type", "Type 1、Type 2、Mixed 或 Other/unclear", "分層／納入條件",
  "analysis_type2_eligible", "糖尿病類型是否為 Type 2；不代表已完成成人及其他後續納入排除", "下游分析篩選輔助欄位",
  "性別", "女性=0、男性=1；基本資料優先", "共變數",
  "Age", "報告日期時之實歲", "共變數來源",
  "測量日期", "檢驗當日或之前最近一次有效身高體重配對之當地日期", "BMI 來源追溯",
  "測量日期時間", "檢驗當日或之前最近一次有效身高體重配對之台灣當地時間", "BMI 來源追溯",
  "basic_time_gap_days", "報告日期與身高體重測量日期之日數差", "QA",
  "BMI_source_age_category", "BMI 來源與檢驗日期之時間差類別", "QA／追溯",
  "BMI", "由同一原始列之身高與體重計算", "共變數",
  "HbA1C", "糖化血色素", "主要 outcome",
  "Glucose_AC", "空腹血糖", "主要 outcome",
  "LDL-C", "低密度脂蛋白膽固醇", "主要 outcome",
  "Triglycerides", "三酸甘油脂", "主要 outcome",
  "HDL-C", "高密度脂蛋白膽固醇", "描述性備用",
  "Total_cholesterol", "總膽固醇", "描述性備用"
)

# ============================================================
# 9. 輸出
# ============================================================

write_xlsx(
  analysis_main,
  main_output_xlsx_path
)

write_csv(
  analysis_main,
  main_output_csv_path,
  na = ""
)

write_xlsx(
  list(
    patient_detail =
      patient_diabetes_type_summary,
    patient_count =
      dm_type_patient_count,
    observation_count =
      dm_type_observation_count
  ),
  type_output_path
)

write_xlsx(
  list(
    lab_value_conflicts =
      lab_value_conflicts,
    final_cohort_summary =
      lab_conflict_final_cohort_summary,
    lab_ingestion_qc =
      lab_ingestion_qc
  ),
  lab_conflict_output_path
)

write_xlsx(
  list(
    overview = overview,
    stage_counts = stage_counts,
    repeat_patient_summary = repeat_patient_summary,
    repeat_count_by_order = repeat_count_by_order,
    missing_summary = missing_summary,
    dm_type_patient_count = dm_type_patient_count,
    dm_type_observation_count = dm_type_observation_count,
    basic_pair_summary = basic_pair_summary,
    basic_anthro_flags = basic_anthro_flag_summary,
    basic_match_summary = basic_match_summary,
    basic_gap_summary = basic_gap_summary,
    basic_gap_patient = basic_gap_patient_summary,
    long_basic_gap_qc = long_basic_gap_export,
    basic_invalid_anthro = basic_invalid_anthro_export,
    basic_timestamp_pair = basic_timestamp_pair_export,
    basic_same_date_qc = basic_same_date_qc,
    basic_time_parse_failure = basic_time_parse_failure_qc,
    selected_after_report = selected_basic_after_report_export,
    same_patient_date_multi_visit = same_patient_date_multiple_visit_qc,
    sex_source_summary = sex_source_summary,
    sex_conflict_qc = sex_conflict_qc,
    basic_sex_conflict = basic_sex_conflict_qc,
    birthdate_conflict = birthdate_conflict_qc,
    event_source_summary = event_source_summary,
    lab_ingestion_qc = lab_ingestion_qc,
    lab_conflict_final = lab_conflict_final_cohort_summary,
    non_numeric_before = non_numeric_before,
    duplicate_anchor = duplicate_anchor_event,
    duplicate_cohort = duplicate_cohort_event,
    duplicate_final = duplicate_final_event,
    repeat_sequence_qc = repeat_sequence_qc,
    data_dictionary = data_dictionary
  ),
  qc_output_path
)

# ============================================================
# 10. 最終自動檢查
# ============================================================

stopifnot(
  nrow(analysis_main) == unique_event_n,
  nrow(duplicate_anchor_event) == 0L,
  nrow(duplicate_cohort_event) == 0L,
  nrow(duplicate_final_event) == 0L,
  nrow(basic_timestamp_pair_qc) == 0L,
  nrow(basic_time_parse_failure_qc) == 0L,
  nrow(selected_basic_after_report_qc) == 0L,
  nrow(repeat_sequence_qc) == 0L,
  !any(is.na(analysis_main$病歷號)),
  !any(is.na(analysis_main$就診號)),
  !any(is.na(analysis_main$報告日期)),
  all(analysis_main$Repeat >= 1L),
  all(
    analysis_main$analysis_type2_eligible ==
      (analysis_main$diabetes_type == "Type 2")
  ),
  all(
    is.na(analysis_main$BMI) |
      abs(
        analysis_main$BMI -
          analysis_main$體重 / (analysis_main$身高 / 100)^2
      ) < 1e-10
  ),
  all(
    is.na(analysis_main$BMI) |
      (
        analysis_main$BMI >= ABS_BMI_RANGE[1] &
          analysis_main$BMI <= ABS_BMI_RANGE[2]
      )
  )
)

message("前處理完成，且最終自動檢查通過。")
message("程式版本：", SCRIPT_VERSION)
message("正式主檔（xlsx）：", main_output_xlsx_path)
message("正式主檔（csv）：", main_output_csv_path)
message("糖尿病類型摘要：", type_output_path)
message("檢驗衝突 QC：", lab_conflict_output_path)
message("前處理 QC：", qc_output_path)

print(overview)
# A tibble: 12 × 2
# item                             value
# <chr>                            <int>
# 1 final_patient_n                   3617
# 2 final_observation_n              29352
# 3 unique_event_n                   29352
# 4 sex_conflict_event_n                 0
# 5 basic_sex_conflict_patient_n         0
# 6 birthdate_conflict_patient_n         0
# 7 lab_conflict_event_n                15
# 8 basic_timestamp_pair_conflict_n      0
# 9 basic_same_date_different_pair_n  4157
# 10 basic_invalid_or_corrected_row_n  6888
# 11 measurement_time_parse_failure_n     0
# 12 selected_basic_after_report_n        0
print(stage_counts)
# A tibble: 5 × 3
# stage                           patient_n observation_n
# <chr>                               <int>         <int>
# 1 糖尿病 cohort 病患數                 3983            NA
# 2 糖尿病 cohort 就診數                 3983         32499
# 3 全部 HbA1c anchor 事件數            16054         78778
# 4 糖尿病 cohort 內 HbA1c 追蹤事件      3617         29352
# 5 最終串接資料                         3617         29352
print(repeat_patient_summary)
# A tibble: 1 × 8
#     patient_n mean_n_obs sd_n_obs median_n_obs q1_n_obs q3_n_obs min_n_obs max_n_obs
#       <int>      <dbl>    <dbl>        <int>    <dbl>    <dbl>     <int>     <int>
# 1      3617       8.12     7.17            5        2       14         1        29
print(dm_type_patient_count)
# A tibble: 4 × 2
# diabetes_type     patient_n
# <chr>                 <int>
# 1 Mixed E10 and E11         4
# 2 Other or unclear        125
# 3 Type 1                    5
# 4 Type 2                 3483
print(dm_type_observation_count)
# A tibble: 4 × 2
# diabetes_type     observation_n
# <chr>                     <int>
# 1 Mixed E10 and E11            36
# 2 Other or unclear           1607
# 3 Type 1                        6
# 4 Type 2                    27703
print(missing_summary)
# A tibble: 10 × 5
# variable          total_n missing_n observed_n missing_rate_percent
# <chr>               <int>     <int>      <int>                <dbl>
# 1 HbA1C               29352         0      29352                  0  
# 2 Glucose_AC          29352     13053      16299                 44.5
# 3 LDL-C               29352      3903      25449                 13.3
# 4 Triglycerides       29352      4066      25286                 13.8
# 5 HDL-C               29352     17982      11370                 61.3
# 6 Total_cholesterol   29352     19320      10032                 65.8
# 7 BMI                 29352      9316      20036                 31.7
# 8 Age                 29352        58      29294                  0.2
# 9 性別                29352         0      29352                  0  
# 10 diabetes_type      29352         0      29352                  0  
print(basic_pair_summary)
# A tibble: 7 × 2
# item                              value
# <chr>                             <int>
# 1 basic_raw_row_n                  131781
# 2 basic_patient_n                   26034
# 3 valid_paired_measurement_n        67831
# 4 timestamp_pair_conflict_n             0
# 5 same_date_different_pair_n         4157
# 6 invalid_or_corrected_row_n         6888
# 7 measurement_time_parse_failure_n      0
print(basic_match_summary)
# # A tibble: 3 × 3
# basic_match_source  observation_n observation_percent
# <chr>                       <int>               <dbl>
# 1 no_prior_valid_pair          9316               31.7 
# 2 prior_pair                  19203               65.4 
# 3 same_day_pair                 833                2.84
print(basic_gap_summary, n=15)
# A tibble: 1 × 15
#     total_observation_n matched_observation_n unmatched_observation_n matched_observation_percent unmatched_observation_per…¹ mean_gap_days sd_gap_days median_gap_days
#                 <int>                 <int>                   <int>                       <dbl>                       <dbl>         <dbl>       <dbl>           <dbl>
# 1               29352                 20036                    9316                        68.3                        31.7          151.        235.              77
# ℹ abbreviated name: ¹​unmatched_observation_percent
# ℹ 7 more variables: q1_gap_days <dbl>, q3_gap_days <dbl>, max_gap_days <dbl>, gap_gt365_n <int>, gap_gt365_percent_of_matched <dbl>, gap_gt730_n <int>,
#   gap_gt730_percent_of_matched <dbl>

# ============================================================
# 基準時點 Systolic / Diastolic 缺失率檢查
# 僅輸出 QA 結果，不加入正式分析主檔
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(writexl)
})

LOCAL_TZ <- "Asia/Taipei"

# ------------------------------------------------------------
# 1. 檔案路徑
# ------------------------------------------------------------

basic_path <- file.path(
  "raw_data",
  "病人基本資料：年齡,性別,身高,體重.xlsx"
)

final_analysis_path <- "analysis_outputs/dt_model_raw_type2_adult.xlsx"

output_path <- "糖尿病患抽血資料整理/baseline_bp_missingness_type2_adult.xlsx"

# ------------------------------------------------------------
# 2. 工具函數
# ------------------------------------------------------------

standardize_id <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

parse_datetime_safe <- function(x) {
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  
  suppressWarnings(
    parse_date_time(
      x_chr,
      orders = c(
        "Ymd HMS",
        "Ymd HM",
        "Ymd",
        "ymd HMS",
        "ymd HM",
        "ymd"
      ),
      tz = LOCAL_TZ
    )
  )
}

clean_numeric_value <- function(x) {
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c(
    "", "NA", "N/A", "NaN", "NULL", "-", "--"
  )] <- NA_character_
  
  x_chr <- str_replace_all(x_chr, ",", "")
  suppressWarnings(as.numeric(x_chr))
}

# ------------------------------------------------------------
# 3. 讀取目前論文使用的最終分析資料
# ------------------------------------------------------------

final_dat <- read_excel(
  final_analysis_path,
  sheet = "dt_model_raw"
) %>%
  mutate(
    ID = standardize_id(ID),
    報告日期 = as.Date(報告日期)
  )

# 每位病患最早的主要分析觀測，定義為基準時點 T1。
baseline_dat <- final_dat %>%
  filter(
    !is.na(ID),
    !is.na(報告日期)
  ) %>%
  arrange(ID, 報告日期, row_id) %>%
  group_by(ID) %>%
  slice(1L) %>%
  ungroup() %>%
  transmute(
    ID,
    baseline_report_date = 報告日期,
    baseline_cutoff = as.POSIXct(
      paste0(
        format(報告日期, "%Y-%m-%d"),
        " 23:59:59"
      ),
      format = "%Y-%m-%d %H:%M:%S",
      tz = LOCAL_TZ
    )
  )

stopifnot(
  n_distinct(baseline_dat$ID) == 1790L
)

# ------------------------------------------------------------
# 4. 讀取原始基本資料中的血壓欄位
# ------------------------------------------------------------

basic_raw <- read_excel(basic_path)

required_cols <- c(
  "病歷號",
  "測量時間",
  "收縮壓",
  "舒張壓"
)

missing_cols <- setdiff(required_cols, names(basic_raw))

if (length(missing_cols) > 0L) {
  stop(
    "基本資料缺少欄位：",
    paste(missing_cols, collapse = ", ")
  )
}

bp_raw <- basic_raw %>%
  transmute(
    ID = standardize_id(病歷號),
    measurement_datetime = parse_datetime_safe(測量時間),
    Systolic = clean_numeric_value(收縮壓),
    Diastolic = clean_numeric_value(舒張壓)
  ) %>%
  filter(
    !is.na(ID),
    !is.na(measurement_datetime)
  ) %>%
  mutate(
    # 0 或負值視為無效紀錄。
    Systolic = if_else(
      !is.na(Systolic) & Systolic > 0,
      Systolic,
      NA_real_
    ),
    Diastolic = if_else(
      !is.na(Diastolic) & Diastolic > 0,
      Diastolic,
      NA_real_
    )
  )

# ------------------------------------------------------------
# 5. 找每位病患在基準日期當日或之前最近一次 Systolic
# ------------------------------------------------------------

baseline_systolic <- baseline_dat %>%
  left_join(
    bp_raw %>%
      filter(!is.na(Systolic)),
    by = "ID",
    relationship = "many-to-many"
  ) %>%
  filter(
    is.na(measurement_datetime) |
      measurement_datetime <= baseline_cutoff
  ) %>%
  group_by(ID) %>%
  arrange(
    desc(measurement_datetime),
    .by_group = TRUE
  ) %>%
  slice(1L) %>%
  ungroup() %>%
  select(
    ID,
    baseline_report_date,
    systolic_measurement_datetime = measurement_datetime,
    Systolic
  )

# 將沒有任何符合條件血壓紀錄的病患補回。
baseline_systolic <- baseline_dat %>%
  select(ID, baseline_report_date) %>%
  left_join(
    baseline_systolic,
    by = c("ID", "baseline_report_date"),
    relationship = "one-to-one"
  )

# ------------------------------------------------------------
# 6. 找每位病患在基準日期當日或之前最近一次 Diastolic
# ------------------------------------------------------------

baseline_diastolic <- baseline_dat %>%
  left_join(
    bp_raw %>%
      filter(!is.na(Diastolic)),
    by = "ID",
    relationship = "many-to-many"
  ) %>%
  filter(
    is.na(measurement_datetime) |
      measurement_datetime <= baseline_cutoff
  ) %>%
  group_by(ID) %>%
  arrange(
    desc(measurement_datetime),
    .by_group = TRUE
  ) %>%
  slice(1L) %>%
  ungroup() %>%
  select(
    ID,
    baseline_report_date,
    diastolic_measurement_datetime = measurement_datetime,
    Diastolic
  )

baseline_diastolic <- baseline_dat %>%
  select(ID, baseline_report_date) %>%
  left_join(
    baseline_diastolic,
    by = c("ID", "baseline_report_date"),
    relationship = "one-to-one"
  )

# ------------------------------------------------------------
# 7. 合併並計算缺失率
# ------------------------------------------------------------

baseline_bp <- baseline_dat %>%
  select(ID, baseline_report_date) %>%
  left_join(
    baseline_systolic,
    by = c("ID", "baseline_report_date")
  ) %>%
  left_join(
    baseline_diastolic,
    by = c("ID", "baseline_report_date")
  ) %>%
  mutate(
    systolic_missing = is.na(Systolic),
    diastolic_missing = is.na(Diastolic),
    bp_pair_missing = is.na(Systolic) | is.na(Diastolic)
  )

bp_missing_summary <- tibble(
  variable = c(
    "Systolic",
    "Diastolic",
    "Systolic and Diastolic pair"
  ),
  total_patient_n = nrow(baseline_bp),
  observed_patient_n = c(
    sum(!baseline_bp$systolic_missing),
    sum(!baseline_bp$diastolic_missing),
    sum(!baseline_bp$bp_pair_missing)
  ),
  missing_patient_n = c(
    sum(baseline_bp$systolic_missing),
    sum(baseline_bp$diastolic_missing),
    sum(baseline_bp$bp_pair_missing)
  ),
  missing_rate_percent = round(
    missing_patient_n / total_patient_n * 100,
    2
  )
)

print(bp_missing_summary)
# A tibble: 3 × 5
# variable                 total_patient_n observed_patient_n missing_patient_n missing_rate_percent
# <chr>                              <int>              <int>             <int>                <dbl>
# 1 Systolic                            1790                902               888                 49.6
# 2 Diastolic                           1790                898               892                 49.8
# 3 Systolic and Diastolic …            1790                898               892                 49.8
