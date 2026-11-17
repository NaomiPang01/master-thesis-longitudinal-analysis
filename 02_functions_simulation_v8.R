# ============================================================
# 02_functions_simulation_updated.R
# 功能：
# 1. 建立模擬所需 empirical skeleton
# 2. 使用人工確認之 rounded LME-RI+RS 參數建立 S0
# 3. 依情境修改隨機效果與誤差變異
# 4. 產生四項主要代謝指標之模擬縱向資料
# ============================================================

# ------------------------------------------------------------
# 0. 全域名稱設定
# ------------------------------------------------------------

SIM_OUTCOMES <- c("HbA1c", "Glu_AC", "LDL", "TG")

SIM_OUTCOME_LABELS <- c(
  HbA1c = "HbA1c",
  Glu_AC = "Glu-AC",
  LDL = "LDL",
  TG = "TG"
)

SIM_GROUP_LABELS <- c(
  Overall = "整體",
  Control = "Control 組",
  Statin = "Statin 組"
)

# ------------------------------------------------------------
# 1. 小工具函數
# ------------------------------------------------------------

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

standardize_sex_to_binary <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    x_chr %in% c("1", "Male", "male", "M", "m", "男") ~ 1L,
    x_chr %in% c("0", "Female", "female", "F", "f", "女") ~ 0L,
    TRUE ~ NA_integer_
  )
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

# ------------------------------------------------------------
# 2. 讀取目前主分析資料作為 empirical skeleton
# ------------------------------------------------------------

load_empirical_skeleton <- function(
    data_path,
    imputation_id = 1,
    seed = 5014,
    create_split_if_missing = FALSE
) {
  
  if (is.na(data_path) || !file.exists(data_path)) {
    stop("找不到 empirical skeleton 資料檔：", data_path)
  }
  
  dat <- readr::read_csv(data_path, show_col_types = FALSE)
  
  # 若使用 imputed_long_with_split_type2.csv，先取第一組補值資料集。
  if (".imp" %in% names(dat)) {
    dat <- dat %>% dplyr::filter(.imp == imputation_id)
  }
  
  # 兼容舊欄位名稱。
  # 兼容舊欄位名稱。避免 dplyr::rename(any_of()) 在欄位不存在時造成版本相容問題。
  if ("HbA1C" %in% names(dat) && !("HbA1c" %in% names(dat))) {
    dat <- dat %>% dplyr::rename(HbA1c = HbA1C)
  }
  if ("Glucose_AC" %in% names(dat) && !("Glu_AC" %in% names(dat))) {
    dat <- dat %>% dplyr::rename(Glu_AC = Glucose_AC)
  }
  if ("Triglycerides" %in% names(dat) && !("TG" %in% names(dat))) {
    dat <- dat %>% dplyr::rename(TG = Triglycerides)
  }
  
  required_cols <- c("ID", "Group", "Time_in_Years", "Age_baseline", "Sex", "BMI")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("empirical skeleton 缺少必要欄位：", paste(missing_cols, collapse = ", "))
  }
  
  dat <- dat %>%
    dplyr::mutate(
      ID = as.character(ID),
      Group = as.character(Group),
      Group = dplyr::case_when(
        Group %in% c("Control", "Control 組") ~ "Control",
        Group %in% c("Statin", "Statin 組") ~ "Statin",
        TRUE ~ Group
      ),
      GroupStatin = ifelse(Group == "Statin", 1L, 0L),
      Sex1 = standardize_sex_to_binary(Sex),
      Time_in_Years = safe_num(Time_in_Years),
      Age_baseline = safe_num(Age_baseline),
      BMI = safe_num(BMI)
    ) %>%
    dplyr::filter(
      Group %in% c("Control", "Statin"),
      is.finite(Time_in_Years),
      is.finite(Age_baseline),
      is.finite(BMI),
      !is.na(Sex1)
    ) %>%
    dplyr::arrange(ID, Time_in_Years)
  
  if (!("split" %in% names(dat)) || all(is.na(dat$split))) {
    if (!isTRUE(create_split_if_missing)) {
      stop("資料中沒有 split 欄位。請先執行表格.R 產生 70/15/15 patient-level split。")
    }
    dat <- create_patient_split(dat, seed = seed)
  } else {
    dat <- dat %>%
      dplyr::mutate(split = as.character(split))
    
    valid_split <- c("train", "validation", "test")
    if (!all(unique(dat$split) %in% valid_split)) {
      stop("split 欄位只能包含 train、validation、test。")
    }
    
    split_check <- dat %>%
      dplyr::distinct(ID, split) %>%
      dplyr::count(ID, name = "n_split") %>%
      dplyr::filter(n_split > 1)
    
    if (nrow(split_check) > 0) {
      stop("同一病患出現在多個 split 中，請檢查表格.R 的資料切割。")
    }
  }
  
  dat
}

# ------------------------------------------------------------
# 3. 病患層級 70/15/15 split 與抽樣
# ------------------------------------------------------------

create_patient_split <- function(df, seed = 5014, train_prop = 0.70, val_prop = 0.15) {
  
  id_group <- df %>%
    dplyr::distinct(ID, Group)
  
  set.seed(seed)
  
  split_df <- id_group %>%
    dplyr::group_by(Group) %>%
    dplyr::group_modify(function(.x, .y) {
      ids <- sample(.x$ID)
      n <- length(ids)
      n_train <- floor(train_prop * n)
      n_val <- floor(val_prop * n)
      
      split_vec <- rep("test", n)
      if (n_train > 0) split_vec[seq_len(n_train)] <- "train"
      if (n_val > 0) {
        split_vec[seq.int(n_train + 1, n_train + n_val)] <- "validation"
      }
      
      tibble::tibble(ID = ids, split = split_vec)
    }) %>%
    dplyr::ungroup()
  
  df %>%
    dplyr::select(-dplyr::any_of("split")) %>%
    dplyr::left_join(split_df, by = "ID")
}

allocate_group_split_sample_size <- function(
    patient_manifest,
    target_n
) {
  allocation <- patient_manifest %>%
    dplyr::count(Group, split, name = "n_available") %>%
    dplyr::mutate(
      exact_n = target_n * n_available / sum(n_available),
      n_sample = floor(exact_n),
      remainder = exact_n - n_sample
    )
  
  n_remaining <- target_n - sum(allocation$n_sample)
  
  if (n_remaining > 0L) {
    add_index <- order(
      allocation$remainder,
      decreasing = TRUE
    )[seq_len(n_remaining)]
    
    allocation$n_sample[add_index] <-
      allocation$n_sample[add_index] + 1L
  }
  
  if (
    any(allocation$n_sample > allocation$n_available) ||
    sum(allocation$n_sample) != target_n
  ) {
    stop("Group × split 分層抽樣配置失敗。", call. = FALSE)
  }
  
  allocation
}

select_patients_by_group_split <- function(
    skeleton,
    target_n,
    seed = 5014
) {
  patient_manifest <- skeleton %>%
    dplyr::distinct(
      ID,
      Group,
      split
    ) %>%
    dplyr::mutate(
      ID = as.character(ID),
      Group = as.character(Group),
      split = as.character(split)
    )
  
  if (
    any(is.na(patient_manifest$ID)) ||
    any(is.na(patient_manifest$Group)) ||
    any(is.na(patient_manifest$split))
  ) {
    stop(
      "病患 manifest 的 ID、Group 或 split 含缺失值。",
      call. = FALSE
    )
  }
  
  patient_split_qc <- patient_manifest %>%
    dplyr::count(
      ID,
      name = "n_strata"
    ) %>%
    dplyr::filter(
      n_strata != 1L
    )
  
  if (nrow(patient_split_qc) > 0L) {
    stop(
      "至少一位病患出現在多個 Group × split 階層。",
      call. = FALSE
    )
  }
  
  available_patient_n <-
    dplyr::n_distinct(
      patient_manifest$ID
    )
  
  if (
    target_n > available_patient_n ||
    target_n <= 0L
  ) {
    stop(
      "target_n 必須介於 1 與可用病患總數之間。",
      call. = FALSE
    )
  }
  
  allocation <-
    allocate_group_split_sample_size(
      patient_manifest =
        patient_manifest,
      target_n = target_n
    ) %>%
    dplyr::mutate(
      Group = as.character(Group),
      split = as.character(split)
    )
  
  set.seed(seed)
  
  selected_patients <- purrr::pmap_dfr(
    allocation,
    function(
    Group,
    split,
    n_available,
    exact_n,
    n_sample,
    remainder
    ) {
      group_value <- as.character(Group)
      split_value <- as.character(split)
      n_take <- as.integer(n_sample)
      
      stratum_patients <-
        patient_manifest %>%
        dplyr::filter(
          .data$Group ==
            .env$group_value,
          .data$split ==
            .env$split_value
        )
      
      if (
        nrow(stratum_patients) !=
        as.integer(n_available)
      ) {
        stop(
          "階層可用病患數與 allocation 不一致：",
          group_value,
          " / ",
          split_value,
          call. = FALSE
        )
      }
      
      if (
        n_take < 0L ||
        n_take >
        nrow(stratum_patients)
      ) {
        stop(
          "階層抽樣數無效：",
          group_value,
          " / ",
          split_value,
          call. = FALSE
        )
      }
      
      if (n_take == 0L) {
        return(
          stratum_patients[
            0,
            ,
            drop = FALSE
          ]
        )
      }
      
      stratum_patients %>%
        dplyr::slice_sample(
          n = n_take,
          replace = FALSE
        )
    }
  ) %>%
    dplyr::distinct(
      ID,
      Group,
      split
    )
  
  selected_patient_n <-
    dplyr::n_distinct(
      selected_patients$ID
    )
  
  if (
    nrow(selected_patients) !=
    target_n ||
    selected_patient_n !=
    target_n
  ) {
    stop(
      paste0(
        "固定分層抽樣後的病患數不等於 target_n：",
        "selected rows = ",
        nrow(selected_patients),
        "；unique IDs = ",
        selected_patient_n,
        "；target_n = ",
        target_n,
        "。"
      ),
      call. = FALSE
    )
  }
  
  selected <- skeleton %>%
    dplyr::mutate(
      ID = as.character(ID),
      Group = as.character(Group),
      split = as.character(split)
    ) %>%
    dplyr::semi_join(
      selected_patients,
      by = c(
        "ID",
        "Group",
        "split"
      )
    ) %>%
    dplyr::arrange(
      ID,
      Time_in_Years,
      row_id
    )
  
  selected_manifest_qc <- selected %>%
    dplyr::distinct(
      ID,
      Group,
      split
    ) %>%
    dplyr::count(
      Group,
      split,
      name = "n_selected"
    ) %>%
    dplyr::full_join(
      allocation %>%
        dplyr::select(
          Group,
          split,
          n_available,
          exact_n,
          n_sample,
          remainder
        ),
      by = c(
        "Group",
        "split"
      )
    ) %>%
    dplyr::mutate(
      passed =
        !is.na(n_selected) &
        n_selected ==
        n_sample
    )
  
  if (
    any(
      !selected_manifest_qc$passed
    )
  ) {
    print(selected_manifest_qc)
    
    stop(
      "至少一個 Group × split 階層的實際抽樣數不符合配置數。",
      call. = FALSE
    )
  }
  
  if (
    dplyr::n_distinct(
      selected$ID
    ) != target_n
  ) {
    stop(
      "固定分層抽樣後的 skeleton 病患數不等於 target_n。",
      call. = FALSE
    )
  }
  
  attr(
    selected,
    "sampling_allocation"
  ) <- selected_manifest_qc
  
  selected
}

summarise_skeleton <- function(skeleton) {
  
  patient_by_group <- skeleton %>%
    dplyr::distinct(ID, Group, split) %>%
    dplyr::count(Group, split, name = "patient_n") %>%
    dplyr::arrange(Group, split)
  
  observation_by_group <- skeleton %>%
    dplyr::count(Group, split, name = "observation_n") %>%
    dplyr::arrange(Group, split)
  
  repeat_summary <- skeleton %>%
    dplyr::count(ID, name = "n_obs") %>%
    dplyr::summarise(
      patient_n = dplyr::n(),
      mean_n_obs = mean(n_obs),
      sd_n_obs = stats::sd(n_obs),
      median_n_obs = stats::median(n_obs),
      q1_n_obs = as.numeric(stats::quantile(n_obs, 0.25)),
      q3_n_obs = as.numeric(stats::quantile(n_obs, 0.75)),
      min_n_obs = min(n_obs),
      max_n_obs = max(n_obs)
    )
  
  list(
    patient_by_group = patient_by_group,
    observation_by_group = observation_by_group,
    repeat_summary = repeat_summary
  )
}

# ------------------------------------------------------------
# 4. 參數設定：預設為目前新版資料分析.R 之 LME-RI+RS 估計方向
# ------------------------------------------------------------


get_default_base_parameters <- function() {
  
  # 預設 S0 參數係依目前實證資料 LME-RI+RS 結果之量級設定，
  # 並採用四捨五入後之簡化數值，不逐字複製原始模型輸出。
  # 為使 S4 情境之「隨機效果相關方向改變」更清楚，
  # S0 中四項 outcome 的 random intercept 與 random slope 相關係數統一設為 -0.40；
  # S4 會將其反轉為 +0.40。
  common_rho <- -0.40
  
  list(
    HbA1c = list(
      fixed = c(
        beta0 = 2.0,
        beta_group = 0.010,
        beta_time = -0.007,
        beta_group_time = 0.006,
        beta_age = -0.002,
        beta_sex_male = -0.006,
        beta_bmi = 0.001
      ),
      sigma_b0 = 0.15,
      sigma_b1 = 0.04,
      rho = common_rho,
      sigma = 0.10
    ),
    
    Glu_AC = list(
      fixed = c(
        beta0 = 5.0,
        beta_group = -0.010,
        beta_time = -0.014,
        beta_group_time = 0.005,
        beta_age = -0.002,
        beta_sex_male = -0.026,
        beta_bmi = 0.002
      ),
      sigma_b0 = 0.20,
      sigma_b1 = 0.04,
      rho = common_rho,
      sigma = 0.20
    ),
    
    LDL = list(
      fixed = c(
        beta0 = 4.8,
        beta_group = -0.100,
        beta_time = -0.012,
        beta_group_time = 0.006,
        beta_age = -0.004,
        beta_sex_male = -0.032,
        beta_bmi = 0.002
      ),
      sigma_b0 = 0.25,
      sigma_b1 = 0.04,
      rho = common_rho,
      sigma = 0.20
    ),
    
    TG = list(
      fixed = c(
        beta0 = 5.2,
        beta_group = -0.040,
        beta_time = -0.018,
        beta_group_time = 0.012,
        beta_age = -0.006,
        beta_sex_male = -0.051,
        beta_bmi = 0.007
      ),
      sigma_b0 = 0.40,
      sigma_b1 = 0.07,
      rho = common_rho,
      sigma = 0.30
    )
  )
}

reorder_fixed_parameters <- function(fixed_vec) {
  # 論文附錄與輸出表一致使用：
  # beta0, beta_group, beta_time, beta_group_time, beta_age, beta_sex_male, beta_bmi
  wanted <- c(
    "beta0",
    "beta_group",
    "beta_time",
    "beta_group_time",
    "beta_age",
    "beta_sex_male",
    "beta_bmi"
  )
  fixed_vec[wanted]
}

round_fixed_parameter_for_simulation <- function(value, par_name) {
  # 與舊版附錄表 E2 的呈現邏輯一致：截距取較粗略尺度，
  # 其他效果取小數點後三位，避免把實證模型估計值逐字複製到模擬設定。
  if (!is.finite(value)) return(NA_real_)
  
  if (par_name == "beta0") {
    return(round(value, 2))
  }
  
  round(value, 3)
}

round_variance_parameter_for_simulation <- function(value, par_name) {
  # 隨機效果與殘差標準差以小數點後二位為主；
  # 相關係數保留二位，以保留 S4 情境之相關方向與大致強度。
  if (!is.finite(value)) return(NA_real_)
  
  if (par_name == "rho") {
    return(round(value, 2))
  }
  
  round(value, 2)
}

update_fixed_from_pooled <- function(base_params, pooled_path, model_to_use = "LME_RI_RS") {
  
  if (is.na(pooled_path) || !file.exists(pooled_path)) {
    message("未找到 pooled_parameter_estimates 檔案，固定效果使用預設 rounded 值。")
    return(base_params)
  }
  
  pooled <- readr::read_csv(pooled_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      outcome = dplyr::recode(
        as.character(outcome),
        HbA1C = "HbA1c",
        Glucose_AC = "Glu_AC",
        Triglycerides = "TG",
        .default = as.character(outcome)
      ),
      term = as.character(term),
      model = as.character(model)
    ) %>%
    dplyr::filter(model == model_to_use)
  
  term_map <- c(
    "(Intercept)" = "beta0",
    "GroupStatin" = "beta_group",
    "Time_in_Years" = "beta_time",
    "Age_baseline" = "beta_age",
    "SexMale" = "beta_sex_male",
    "Sex1" = "beta_sex_male",
    "BMI" = "beta_bmi",
    "GroupStatin:Time_in_Years" = "beta_group_time"
  )
  
  for (outcome_name in names(base_params)) {
    
    sub <- pooled %>%
      dplyr::filter(outcome == outcome_name, term %in% names(term_map))
    
    if (nrow(sub) == 0) {
      warning("pooled estimates 中找不到 outcome = ", outcome_name, " 的 ", model_to_use, " 固定效果，沿用預設 rounded 值。")
      next
    }
    
    fixed_new <- base_params[[outcome_name]]$fixed
    
    for (j in seq_len(nrow(sub))) {
      par_name <- unname(term_map[sub$term[j]])
      fixed_new[par_name] <- round_fixed_parameter_for_simulation(
        value = as.numeric(sub$estimate[j]),
        par_name = par_name
      )
    }
    
    base_params[[outcome_name]]$fixed <- reorder_fixed_parameters(fixed_new)
  }
  
  message("已使用 ", pooled_path, " 更新固定效果參數：", model_to_use, "；所有參數均已 rounded。")
  base_params
}

update_varcomp_from_file <- function(base_params, varcomp_path, model_to_use = "LME_RI_RS") {
  
  if (is.na(varcomp_path) || !file.exists(varcomp_path)) {
    message("未找到 lme_variance_components_summary_for_paper 檔案，隨機效果與殘差變異使用預設 rounded 值。")
    return(base_params)
  }
  
  vc <- readr::read_csv(varcomp_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      outcome = dplyr::recode(
        as.character(outcome),
        HbA1C = "HbA1c",
        Glucose_AC = "Glu_AC",
        Triglycerides = "TG",
        .default = as.character(outcome)
      ),
      model = as.character(model)
    ) %>%
    dplyr::filter(model == model_to_use)
  
  needed <- c("sigma_b0_sq_mean", "sigma_b1_sq_mean", "rho_b0b1_mean", "sigma_sq_mean")
  if (!all(needed %in% names(vc))) {
    warning("variance components 檔案缺少必要欄位，沿用預設 rounded 值。")
    return(base_params)
  }
  
  # 先根據各 outcome 最新 LME-RI+RS variance components 更新標準差。
  # rho 則不逐 outcome 直接複製；本研究為使 S4 情境定義清楚，
  # 採四項 outcome rho 估計值平均後，四捨五入為共同負相關。
  rho_pool <- c()
  
  for (outcome_name in names(base_params)) {
    
    sub <- vc %>% dplyr::filter(outcome == outcome_name)
    if (nrow(sub) == 0) {
      warning("variance components 中找不到 outcome = ", outcome_name, "，沿用預設 rounded 值。")
      next
    }
    
    base_params[[outcome_name]]$sigma_b0 <- round_variance_parameter_for_simulation(
      sqrt(max(as.numeric(sub$sigma_b0_sq_mean[1]), 0)),
      "sigma_b0"
    )
    base_params[[outcome_name]]$sigma_b1 <- round_variance_parameter_for_simulation(
      sqrt(max(as.numeric(sub$sigma_b1_sq_mean[1]), 0)),
      "sigma_b1"
    )
    base_params[[outcome_name]]$sigma <- round_variance_parameter_for_simulation(
      sqrt(max(as.numeric(sub$sigma_sq_mean[1]), 0)),
      "sigma"
    )
    
    rho_pool <- c(rho_pool, as.numeric(sub$rho_b0b1_mean[1]))
  }
  
  if (length(rho_pool) > 0 && any(is.finite(rho_pool))) {
    common_rho <- round(mean(rho_pool[is.finite(rho_pool)], na.rm = TRUE), 1)
    common_rho <- -abs(common_rho)
    for (outcome_name in names(base_params)) {
      base_params[[outcome_name]]$rho <- common_rho
    }
    message("四項 outcome 的 rho 估計值相近，已簡化為共同 S0 rho = ", common_rho,
            "；S4 會反轉為正相關。")
  }
  
  message("已使用 ", varcomp_path, " 更新 LME variance components：", model_to_use, "；所有標準差均已 rounded，rho 採共同 rounded 值。")
  base_params
}

get_base_parameters <- function(
    pooled_path = NA_character_,
    varcomp_path = NA_character_,
    use_current_fixed = TRUE,
    use_current_varcomp = TRUE
) {
  
  base_params <- get_default_base_parameters()
  
  if (isTRUE(use_current_fixed)) {
    base_params <- update_fixed_from_pooled(
      base_params = base_params,
      pooled_path = pooled_path,
      model_to_use = "LME_RI_RS"
    )
  }
  
  if (isTRUE(use_current_varcomp)) {
    base_params <- update_varcomp_from_file(
      base_params = base_params,
      varcomp_path = varcomp_path,
      model_to_use = "LME_RI_RS"
    )
  }
  
  base_params
}

base_parameters_to_table <- function(base_params) {
  
  purrr::map_dfr(names(base_params), function(outcome_name) {
    par <- base_params[[outcome_name]]
    
    tibble::tibble(
      outcome = outcome_name,
      outcome_label = SIM_OUTCOME_LABELS[[outcome_name]],
      beta0 = par$fixed["beta0"],
      beta_group = par$fixed["beta_group"],
      beta_time = par$fixed["beta_time"],
      beta_group_time = par$fixed["beta_group_time"],
      beta_age = par$fixed["beta_age"],
      beta_sex_male = par$fixed["beta_sex_male"],
      beta_bmi = par$fixed["beta_bmi"],
      sigma_b0 = par$sigma_b0,
      sigma_b1 = par$sigma_b1,
      rho_b0b1 = par$rho,
      sigma = par$sigma,
      parameter_source = "rounded_LME_RI_RS_estimates"
    )
  })
}

# ------------------------------------------------------------
# 5. 依情境修改參數
# ------------------------------------------------------------

modify_one_outcome_parameters <- function(par, scenario_setting) {
  
  sb_mult <- as.numeric(scenario_setting$sb_mult[1])
  sigma_mult <- as.numeric(scenario_setting$sigma_mult[1])
  rho_mode <- as.character(scenario_setting$rho_mode[1])
  
  par$sigma_b0 <- par$sigma_b0 * sb_mult
  par$sigma_b1 <- par$sigma_b1 * sb_mult
  par$sigma <- par$sigma * sigma_mult
  
  if (rho_mode == "positive") {
    par$rho <- abs(par$rho)
  } else if (rho_mode == "negative") {
    par$rho <- -abs(par$rho)
  } else if (rho_mode == "zero") {
    par$rho <- 0
  }
  
  par
}

apply_scenario_to_parameters <- function(base_params, scenario_setting) {
  
  # 將 S0--S4 的情境設定一次套用到四項主要代謝指標。
  # S1/S2：調整 random-effect SD；S3：調整 residual SD；S4：改變 rho 方向。
  # 這個函數名稱保留給主程式與論文附錄對照使用，避免情境修改邏輯散落在資料生成步驟中。
  
  purrr::map(
    base_params,
    ~ modify_one_outcome_parameters(.x, scenario_setting = scenario_setting)
  )
}

# 向下相容：若舊程式或使用者測試碼仍呼叫 modify_parameters()，
# 仍回傳單一 outcome 的情境修改後參數。
modify_parameters <- function(par, scenario_setting) {
  modify_one_outcome_parameters(par = par, scenario_setting = scenario_setting)
}

build_D <- function(sd_intercept, sd_slope, rho) {
  
  D <- matrix(
    c(
      sd_intercept^2,
      rho * sd_intercept * sd_slope,
      rho * sd_intercept * sd_slope,
      sd_slope^2
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  # 保護：若數值誤差導致非正定，微量調整對角線。
  eig <- eigen(D, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig <= 0)) {
    D <- D + diag(abs(min(eig)) + 1e-6, 2)
  }
  
  D
}

# ------------------------------------------------------------
# 6. 模擬資料生成
# ------------------------------------------------------------

simulate_one_outcome <- function(df, par, outcome_name) {
  
  ids <- unique(df$ID)
  
  D <- build_D(
    sd_intercept = par$sigma_b0,
    sd_slope = par$sigma_b1,
    rho = par$rho
  )
  
  b <- MASS::mvrnorm(
    n = length(ids),
    mu = c(0, 0),
    Sigma = D
  )
  
  bdf <- tibble::tibble(
    ID = ids,
    b0 = b[, 1],
    b1 = b[, 2]
  )
  
  df2 <- df %>%
    dplyr::left_join(bdf, by = "ID")
  
  fixed <- par$fixed
  
  eta <- with(
    df2,
    fixed["beta0"] +
      fixed["beta_group"] * GroupStatin +
      fixed["beta_time"] * Time_in_Years +
      fixed["beta_group_time"] * GroupStatin * Time_in_Years +
      fixed["beta_age"] * Age_baseline +
      fixed["beta_sex_male"] * Sex1 +
      fixed["beta_bmi"] * BMI +
      b0 +
      b1 * Time_in_Years
  )
  
  log_y <- eta + stats::rnorm(nrow(df2), mean = 0, sd = par$sigma)
  
  df2[[outcome_name]] <- exp(log_y)
  
  df2 %>%
    dplyr::select(-b0, -b1)
}

simulate_full_dataset <- function(skeleton, base_params, scenario_setting, seed) {
  
  set.seed(seed)
  
  skeleton_row_id <- as.integer(skeleton$row_id)
  skeleton_n <- nrow(skeleton)
  skeleton_patient_n <- dplyr::n_distinct(skeleton$ID)
  
  out <- skeleton
  scenario_params <- apply_scenario_to_parameters(
    base_params = base_params,
    scenario_setting = scenario_setting
  )
  
  for (outcome_name in names(scenario_params)) {
    out <- simulate_one_outcome(
      df = out,
      par = scenario_params[[outcome_name]],
      outcome_name = outcome_name
    )
  }
  
  if (
    nrow(out) != skeleton_n ||
    dplyr::n_distinct(out$ID) != skeleton_patient_n ||
    !identical(as.integer(out$row_id), skeleton_row_id)
  ) {
    stop("模擬生成後的資料列數、病患數或 row_id 與固定 skeleton 不一致。", call. = FALSE)
  }
  
  for (outcome_name in names(scenario_params)) {
    if (
      any(!is.finite(out[[outcome_name]])) ||
      any(out[[outcome_name]] <= 0)
    ) {
      stop("模擬 outcome 含非有限值或非正值：", outcome_name, call. = FALSE)
    }
  }
  
  out
}

# ------------------------------------------------------------
# 7. 參考平均軌跡
# ------------------------------------------------------------

create_reference_trajectories <- function(skeleton, base_params, tgrid = seq(0, 4.75, by = 0.25)) {
  
  all_ref <- list()
  
  age_mean <- mean(skeleton$Age_baseline, na.rm = TRUE)
  bmi_mean <- mean(skeleton$BMI, na.rm = TRUE)
  sex_mode <- ifelse(mean(skeleton$Sex1, na.rm = TRUE) >= 0.5, 1L, 0L)
  
  for (outcome_name in names(base_params)) {
    
    par <- base_params[[outcome_name]]
    fixed <- par$fixed
    
    newdat <- expand.grid(
      Group = c("Control", "Statin"),
      Time_in_Years = tgrid
    )
    
    newdat <- newdat %>%
      dplyr::mutate(
        GroupStatin = ifelse(Group == "Statin", 1L, 0L),
        Age_baseline = age_mean,
        Sex1 = sex_mode,
        BMI = bmi_mean
      )
    
    eta <- with(
      newdat,
      fixed["beta0"] +
        fixed["beta_group"] * GroupStatin +
        fixed["beta_time"] * Time_in_Years +
        fixed["beta_group_time"] * GroupStatin * Time_in_Years +
        fixed["beta_age"] * Age_baseline +
        fixed["beta_sex_male"] * Sex1 +
        fixed["beta_bmi"] * BMI
    )
    
    newdat$fit <- exp(eta)
    newdat$outcome <- outcome_name
    newdat$outcome_label <- SIM_OUTCOME_LABELS[[outcome_name]]
    newdat$Group_label <- SIM_GROUP_LABELS[newdat$Group]
    
    all_ref[[outcome_name]] <- newdat
  }
  
  dplyr::bind_rows(all_ref) %>%
    dplyr::mutate(
      outcome = factor(outcome, levels = names(base_params)),
      outcome_label = factor(outcome_label, levels = SIM_OUTCOME_LABELS[names(base_params)]),
      Group = factor(Group, levels = c("Control", "Statin")),
      Group_label = factor(Group_label, levels = unname(SIM_GROUP_LABELS[c("Control", "Statin")]))
    )
}
