# ESJ 风格脊柱外科统计 pipeline（per 1-point AAC effect）
# 请根据需要修改 DATA_PATH / OUT_DIR / FIG_DIR，然后运行整份脚本。

## 基础设置 -------------------------------------------------------------

# 数据与输出路径
DATA_PATH <- "data.csv"          # 请替换为实际数据路径
OUT_DIR   <- "outputs"           # 结果表格输出目录
FIG_DIR   <- file.path(OUT_DIR, "figs")  # 图像输出目录

# 统一加载所需包
lib_vec <- c(
  "tidyverse", "readr", "janitor", "gtsummary", "flextable", "officer",
  "broom", "sandwich", "lmtest", "rms", "irr", "clinfun", "logistf",
  "brglm2", "pwr", "patchwork"
)
invisible(lapply(lib_vec, require, character.only = TRUE))

# 创建结果/图像目录
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# 智能读数：UTF-8 失败则回退 GB18030
read_try <- function(path) {
  tryCatch(readr::read_csv(path, locale = locale(encoding = "UTF-8")),
           error = function(e) {
             message("UTF-8 读取失败，尝试 GB18030 …")
             readr::read_csv(path, locale = locale(encoding = "GB18030"))
           })
}

d <- read_try(DATA_PATH) %>% janitor::clean_names()

## 辅助函数 -------------------------------------------------------------

# 列名自动匹配
get_col <- function(df, candidates) {
  cand_lower <- tolower(candidates)
  nm_lower <- tolower(names(df))
  hit <- names(df)[match(cand_lower, nm_lower, nomatch = 0)][1]
  if (is.na(hit) || is.null(hit) || hit == "") return(NULL) else return(df[[hit]])
}

# 生成衍生变量
prep_data <- function(df) {
  df %>% mutate(
    Kp1 = get_col(df, c("Kp1", "AAC_Kp1", "Kp1_score")),
    Kp2 = get_col(df, c("Kp2", "AAC_Kp2")),
    Kauppila = get_col(df, c("Kauppila", "AAC", "AAC_total")),
    Fusion1 = get_col(df, c("Fusion1", "fusion_r1")),
    Fusion2 = get_col(df, c("Fusion2", "fusion_r2")),
    Fusion = coalesce(as.numeric(Fusion1), as.numeric(Fusion2), get_col(df, c("Fusion"))),
    JOA_pre = get_col(df, c("JOA_pre", "JOA0", "JOA_baseline")),
    JOA_1y = get_col(df, c("JOA_1y", "JOA1", "JOA_12m")),
    ODI_pre = get_col(df, c("ODI_pre", "ODI0", "ODI_baseline")),
    ODI_1y = get_col(df, c("ODI_1y", "ODI1", "ODI_12m")),
    YVAS = get_col(df, c("YVAS", "VAS_back", "Back_VAS")),
    TVAS = get_col(df, c("TVAS", "VAS_leg", "Leg_VAS")),
    VAS_total = get_col(df, c("VAS_total", "VAS_sum")),
    Gender = get_col(df, c("Gender", "Sex")),
    Age = get_col(df, c("Age")),
    BMI = get_col(df, c("BMI")),
    HT = get_col(df, c("HT", "Hypertension")),
    HP = get_col(df, c("HP", "Hyperlipidemia")),
    Diabetes = get_col(df, c("Diabetes")),
    AHUS = coalesce(get_col(df, c("AHUS", "HU", "AHU")), get_col(df, c("HU"))),
    Macnab = get_col(df, c("Macnab")),
    MF = get_col(df, c("MF", "MuscleFat")),
    OpTime = get_col(df, c("OpTime", "Operation_time")),
    BloodLoss = get_col(df, c("BloodLoss", "Blood_loss")),
    Tobacco = get_col(df, c("Tobacco", "Smoking", "Smoke"))
  ) %>%
    mutate(
      Fusion = as.numeric(Fusion),
      Fusion = case_when(Fusion1 %in% c(1, "1", "Yes") ~ 1,
                         Fusion2 %in% c(1, "1", "Yes") ~ 1,
                         TRUE ~ Fusion),
      Fusion = if_else(is.na(Fusion), NA_real_, if_else(Fusion >= 0.5, 1, 0)),
      Gender01 = case_when(tolower(as.character(Gender)) %in% c("m", "male", "1", "男") ~ 1,
                           tolower(as.character(Gender)) %in% c("f", "female", "0", "女") ~ 0,
                           TRUE ~ NA_real_),
      Smoking01 = case_when(is.na(Tobacco) ~ NA_character_,
                            Tobacco %in% c(1, "1", "Yes", "Y") ~ "1",
                            Tobacco %in% c(0, "0", "No", "N") ~ "0",
                            TRUE ~ as.character(Tobacco)) %>% factor(levels = c("0", "1")),
      AAC_score = coalesce(Kp1, Kauppila),
      AAC01 = if_else(AAC_score >= 1, 1, 0, missing = NA_real_),
      AAC_cat4 = case_when(AAC_score >= 9 ~ "≥9",
                           AAC_score >= 5 ~ "5–8",
                           AAC_score >= 1 ~ "1–4",
                           AAC_score == 0 ~ "0",
                           TRUE ~ NA_character_) %>% factor(levels = c("0", "1–4", "5–8", "≥9")),
      dJOA = JOA_1y - JOA_pre,
      dODI = ODI_1y - ODI_pre,
      MCID10 = if_else(dODI >= 10, 1, 0, missing = NA_real_),
      MCID30 = if_else(dODI / ODI_pre >= 0.30, 1, 0, missing = NA_real_),
      MacLab = case_when(Macnab == 1 ~ "Excellent",
                         Macnab == 2 ~ "Good",
                         Macnab == 3 ~ "Fair",
                         Macnab == 4 ~ "Poor",
                         TRUE ~ NA_character_) %>% factor(levels = c("Excellent", "Good", "Fair", "Poor"))
    )
}

d0 <- prep_data(d)

# AAC 审计与暴露选择

audit_aac <- function(df) {
  kp1_rng <- range(df$Kp1, na.rm = TRUE)
  kau_rng <- range(df$Kauppila, na.rm = TRUE)
  message(sprintf("Kp1 range: %s", paste(kp1_rng, collapse = "-")))
  message(sprintf("Kauppila range: %s", paste(kau_rng, collapse = "-")))
  if (is.finite(kau_rng[2]) && kau_rng[2] <= 4) message("提示：Kauppila 最大值 ≤4，疑似分档")
}

audit_aac(d0)

choose_exposure <- function(df) {
  if (sum(!is.na(df$Kp1)) >= 3) {
    message("使用 Kp1 作为 expo_var")
    return("Kp1")
  }
  message("使用 Kauppila 作为 expo_var")
  "Kauppila"
}

expo_var <- choose_exposure(d0)

## Table 1：基线特征 -----------------------------------------------------

# P 值角标
p_sup <- function(x) as_sup(" ", x)

# 统计设置
cont_stats_mean <- list(all_continuous() ~ "{mean} ± {sd}")
cont_stats_median <- list(all_continuous() ~ "{median} ({p25}, {p75})")

# 自定义检验函数
p_tests_t <- function(data, variable, by, ...) t.test(data[[variable]] ~ data[[by]])$p.value
p_tests_wilcox <- function(data, variable, by, ...) wilcox.test(data[[variable]] ~ data[[by]])$p.value
p_tests_chi <- function(data, variable, by, ...) chisq.test(table(data[[variable]], data[[by]]), simulate.p.value = TRUE)$p.value

# Table 1 变量
vars_t1 <- c("Age", "BMI", "Gender01", "HT", "HP", "Diabetes", "AHUS", "MF", "OpTime", "BloodLoss", "Smoking01")

create_table1 <- function(df) {
  df %>%
    select(all_of(c("AAC01", vars_t1))) %>%
    tbl_summary(
      by = AAC01,
      statistic = cont_stats_median,
      digits = everything() ~ 1,
      missing = "no"
    ) %>%
    modify_statistic(BMI ~ "{mean} ± {sd}") %>%
    add_p(
      test = list(
        BMI ~ "p_tests_t",
        all_continuous() ~ "p_tests_wilcox",
        all_categorical() ~ "p_tests_chi"
      )
    ) %>%
    modify_footnote(everything() ~ NA) %>%
    modify_header(label ~ "变量") %>%
    modify_caption("Table 1. Baseline characteristics by AAC01") %>%
    as_flex_table() %>%
    theme_vanilla() %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    hline_top(part = "header", border = fp_border(width = 1.5)) %>%
    hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(width = 0.8)) %>%
    set_caption("Table 1. Baseline characteristics by AAC01") %>%
    footnote(part = "footer", value = as_paragraph(
      as_sup("1"), " t 检验；", as_sup("2"), " Wilcoxon；",
      as_sup("3"), " 卡方 (模拟 p)"), ref_symbols = "1")
}

tbl1_ft <- create_table1(d0)

## Table 2：1 年结局 -----------------------------------------------------

vars_t2 <- c("JOA_pre", "JOA_1y", "ODI_pre", "ODI_1y", "YVAS", "TVAS", "VAS_total",
             "MacLab", "Fusion", "MCID10", "MCID30", "Smoking01")

create_table2 <- function(df) {
  df %>%
    select(all_of(c("AAC01", vars_t2))) %>%
    tbl_summary(
      by = AAC01,
      statistic = cont_stats_median,
      digits = everything() ~ 1,
      missing = "no"
    ) %>%
    modify_statistic(c(YVAS, TVAS, VAS_total) ~ "{median} ({p25}, {p75})") %>%
    add_p(test = list(all_continuous() ~ "p_tests_wilcox", all_categorical() ~ "p_tests_chi")) %>%
    modify_header(label ~ "变量") %>%
    modify_caption("Table 2. One-year outcomes by AAC01") %>%
    as_flex_table() %>%
    theme_vanilla() %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    hline_top(part = "header", border = fp_border(width = 1.5)) %>%
    hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(width = 0.8)) %>%
    set_caption("Table 2. One-year outcomes by AAC01") %>%
    footnote(part = "footer", value = as_paragraph(
      as_sup("1"), " Wilcoxon；", as_sup("3"), " 卡方 (模拟 p)"), ref_symbols = "1")
}

tbl2_ft <- create_table2(d0)

## Table 3：一致性 -------------------------------------------------------

calc_icc_kappa <- function(df) {
  res <- list()
  if (sum(complete.cases(df[, c("Kp1", "Kp2")])) >= 5) {
    icc_res <- irr::icc(df[, c("Kp1", "Kp2")], model = "twoway", type = "agreement", unit = "single")
    res$AAC <- tibble(
      Measure = "AAC (Kp1 vs Kp2)",
      Method = "ICC(A,1)",
      Estimate = icc_res$value,
      CI_low = icc_res$lbound,
      CI_high = icc_res$ubound,
      p = icc_res$p.value
    )
  }
  if (sum(complete.cases(df[, c("Fusion1", "Fusion2")])) >= 5) {
    kap <- irr::kappa2(df[, c("Fusion1", "Fusion2")])
    boot_fun <- function(i) irr::kappa2(df[i, c("Fusion1", "Fusion2")])$value
    boot_vals <- replicate(500, boot_fun(sample(seq_len(nrow(df)), replace = TRUE)), simplify = TRUE)
    ci <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
    res$Fusion <- tibble(
      Measure = "Fusion (R1 vs R2)",
      Method = "Cohen's kappa",
      Estimate = kap$value,
      CI_low = ci[1],
      CI_high = ci[2],
      p = kap$p.value
    )
  }
  bind_rows(res)
}

icc_tbl <- calc_icc_kappa(d0)

tbl3_ft <- icc_tbl %>%
  mutate(Estimate = round(Estimate, 3), CI = sprintf("%.3f–%.3f", CI_low, CI_high)) %>%
  select(Measure, Method, Estimate, CI, p) %>%
  flextable() %>%
  theme_vanilla() %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  set_caption("Table 3. Interrater agreement")

## Table 4：趋势检验（AAC 四分类） --------------------------------------

trend_test <- function(x, g) {
  ok <- complete.cases(x, g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 3 || length(unique(x)) < 3) return(list(stat = NA, p = NA, method = "Spearman"))
  jt <- try(clinfun::jonckheere.test(x, g), silent = TRUE)
  if (inherits(jt, "try-error")) {
    rho <- cor.test(x, as.numeric(g), method = "spearman")
    list(stat = rho$estimate, p = rho$p.value, method = "Spearman")
  } else {
    list(stat = jt$statistic, p = jt$p.value, method = "Jonckheere–Terpstra")
  }
}

tbl4_df <- bind_rows(
  {
    res <- trend_test(d0$JOA_1y, d0$AAC_cat4)
    tibble(Outcome = "JOA_1y", Statistic = res$stat, p_trend = res$p,
           n = sum(complete.cases(d0$JOA_1y, d0$AAC_cat4)), Method = res$method)
  },
  {
    res <- trend_test(d0$ODI_1y, d0$AAC_cat4)
    tibble(Outcome = "ODI_1y", Statistic = res$stat, p_trend = res$p,
           n = sum(complete.cases(d0$ODI_1y, d0$AAC_cat4)), Method = res$method)
  }
)

tbl4_ft <- tbl4_df %>%
  flextable() %>%
  theme_vanilla() %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  set_caption("Table 4. Trend test across AAC categories")

## ANCOVA：JOA/ODI (HC3) -------------------------------------------------

prune_vars <- function(df, vars_keep, outcome, threshold = 30) {
  cand <- setdiff(vars_keep, outcome)
  while (sum(complete.cases(df[, cand])) < threshold && length(cand) > 3) {
    na_counts <- sapply(df[, cand], function(x) sum(is.na(x)))
    drop_var <- names(sort(na_counts, decreasing = TRUE))[1]
    cand <- setdiff(cand, drop_var)
  }
  cand
}

hc3_glance <- function(fit) {
  co <- coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC3"))
  broom::tidy(co, conf.int = TRUE, conf.level = 0.95)
}

run_ancova <- function(df, outcome, exposure, adjust_vars) {
  df1 <- df %>% select(all_of(c(outcome, exposure, adjust_vars))) %>% drop_na()
  if (nrow(df1) < 10) return(tibble())
  formula_str <- sprintf("%s ~ %s + %s", outcome, exposure, paste(setdiff(adjust_vars, outcome), collapse = " + "))
  fit <- lm(as.formula(formula_str), data = df1)
  res <- hc3_glance(fit) %>% filter(term == exposure) %>%
    transmute(Outcome = outcome, term, estimate, conf.low, conf.high, p.value, n = nrow(df1))
  res
}

full_covars_JOA <- c("JOA_1y", expo_var, "JOA_pre", "Fusion", "Gender01", "Age", "BMI", "HT", "HP", "Diabetes", "AHUS", "MF", "OpTime", "BloodLoss", "Tobacco")
min_covars_JOA <- c("JOA_1y", expo_var, "JOA_pre", "Age", "Gender01", "BMI", "Tobacco")
full_covars_ODI <- c("ODI_1y", expo_var, "ODI_pre", "Fusion", "Gender01", "Age", "BMI", "HT", "HP", "Diabetes", "AHUS", "MF", "OpTime", "BloodLoss", "Tobacco")
min_covars_ODI <- c("ODI_1y", expo_var, "ODI_pre", "Age", "Gender01", "BMI", "Tobacco")

covars_JOA_full <- prune_vars(d0, full_covars_JOA, "JOA_1y")
covars_JOA_min <- prune_vars(d0, min_covars_JOA, "JOA_1y")
covars_ODI_full <- prune_vars(d0, full_covars_ODI, "ODI_1y")
covars_ODI_min <- prune_vars(d0, min_covars_ODI, "ODI_1y")

ancova_res <- bind_rows(
  run_ancova(d0, "JOA_1y", expo_var, covars_JOA_full),
  run_ancova(d0, "JOA_1y", expo_var, covars_JOA_min),
  run_ancova(d0, "ODI_1y", expo_var, covars_ODI_full),
  run_ancova(d0, "ODI_1y", expo_var, covars_ODI_min)
) %>%
  mutate(Model = c("JOA_full", "JOA_min", "ODI_full", "ODI_min")) %>%
  select(Model, everything())

ancova_ft <- ancova_res %>%
  mutate(beta = round(estimate, 3), CI = sprintf("%.3f–%.3f", conf.low, conf.high)) %>%
  select(Model, Outcome, beta, CI, p.value, n) %>%
  flextable() %>%
  theme_vanilla() %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  set_caption("ANCOVA results (per 1-point AAC)")

## RCS：非线性 -----------------------------------------------------------

rcs_knots <- function(x) {
  ux <- unique(na.omit(x))
  if (length(ux) >= 5) return(4)
  if (length(ux) == 3) return(3)
  return(0) # 线性
}

plot_rcs <- function(outcome, pre_var, ylab) {
  nk <- rcs_knots(d0$Kp1)
  if (nk == 0 || sum(complete.cases(d0[, c("Kp1", outcome, pre_var)])) < 20) return(NULL)
  dd <- datadist(d0); options(datadist = "dd")
  covars <- c(pre_var, "Age", "Gender01", "BMI", "Tobacco")
  df1 <- d0 %>% select(all_of(c(outcome, "Kp1", covars))) %>% drop_na()
  fml <- as.formula(sprintf("%s ~ rcs(Kp1, %s) + %s", outcome, nk, paste(covars, collapse = " + ")))
  fit <- rms::ols(fml, data = df1)
  pdat <- rms::Predict(fit, Kp1, fun = identity)
  ggplot(pdat, aes(x = Kp1, y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#6baed6") +
    geom_line(color = "#08519c", size = 1) +
    labs(x = "AAC (Kauppila 0–24)", y = ylab) +
    coord_cartesian(ylim = ifelse(outcome == "ODI_1y", c(0, 60), c(0, 30))) +
    theme_bw()
}

p_ODI <- plot_rcs("ODI_1y", "ODI_pre", "Adjusted 1-year ODI")
p_JOA <- plot_rcs("JOA_1y", "JOA_pre", "Adjusted 1-year JOA")

if (!is.null(p_ODI)) ggsave(file.path(FIG_DIR, "Fig_ODI_RCS.tiff"), p_ODI, width = 5, height = 4, dpi = 300)
if (!is.null(p_JOA)) ggsave(file.path(FIG_DIR, "Fig_JOA_RCS.tiff"), p_JOA, width = 5, height = 4, dpi = 300)
if (!is.null(p_ODI) && !is.null(p_JOA)) {
  p_combo <- p_ODI + p_JOA + plot_annotation(tag_levels = "A")
  ggsave(file.path(FIG_DIR, "Fig_RCS_ODI_JOA.tiff"), p_combo, width = 10, height = 4, dpi = 300)
}

## MCID 逻辑回归 + EPV ---------------------------------------------------

run_logit <- function(df, outcome, exposure, covars) {
  df1 <- df %>% select(all_of(c(outcome, exposure, covars))) %>% drop_na()
  if (nrow(df1) < 20) return(list(res = tibble(), epv = tibble()))
  fml <- as.formula(sprintf("%s ~ %s + %s", outcome, exposure, paste(setdiff(covars, outcome), collapse = " + ")))
  fit <- try(logistf(fml, data = df1), silent = TRUE)
  if (inherits(fit, "try-error")) fit <- brglm2::brglm(fml, data = df1, family = binomial())
  tidy_res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  res <- tidy_res %>% filter(term == exposure) %>%
    transmute(Outcome = outcome, OR = estimate, CI_low = conf.low, CI_high = conf.high, p.value, n = nrow(df1))
  events <- sum(df1[[outcome]] == 1, na.rm = TRUE)
  params <- length(coef(fit))
  epv <- tibble(Model = outcome, events = events, params = params, EPV = events / params)
  list(res = res, epv = epv)
}

logit_covars_full <- c("Fusion", "Gender01", "Age", "BMI", "HT", "HP", "Diabetes", "AHUS", "MF", "OpTime", "BloodLoss", "Tobacco")
logit_covars_min <- c("Age", "Gender01", "BMI", "Tobacco")

logit_res <- bind_rows(
  run_logit(d0, "MCID10", expo_var, logit_covars_full)$res %>% mutate(Model = "MCID10_full"),
  run_logit(d0, "MCID10", expo_var, logit_covars_min)$res %>% mutate(Model = "MCID10_min"),
  run_logit(d0, "MCID30", expo_var, logit_covars_full)$res %>% mutate(Model = "MCID30_full"),
  run_logit(d0, "MCID30", expo_var, logit_covars_min)$res %>% mutate(Model = "MCID30_min")
)

logit_epv <- bind_rows(
  run_logit(d0, "MCID10", expo_var, logit_covars_full)$epv %>% mutate(Model = "MCID10_full"),
  run_logit(d0, "MCID10", expo_var, logit_covars_min)$epv %>% mutate(Model = "MCID10_min"),
  run_logit(d0, "MCID30", expo_var, logit_covars_full)$epv %>% mutate(Model = "MCID30_full"),
  run_logit(d0, "MCID30", expo_var, logit_covars_min)$epv %>% mutate(Model = "MCID30_min")
)

logit_ft <- logit_res %>%
  mutate(OR = round(OR, 3), CI = sprintf("%.3f–%.3f", CI_low, CI_high)) %>%
  select(Model, Outcome, OR, CI, p.value, n) %>%
  flextable() %>%
  theme_vanilla() %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  set_caption("Logistic regression for MCID")

logit_epv_ft <- logit_epv %>% flextable() %>% theme_vanilla() %>% fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>% set_caption("Events per variable")

## 导出 Word / PPT -------------------------------------------------------

doc <- read_docx() %>%
  body_add_flextable(tbl1_ft) %>% body_add_par(" ") %>%
  body_add_flextable(tbl2_ft) %>% body_add_par(" ") %>%
  body_add_flextable(tbl3_ft) %>% body_add_par(" ") %>%
  body_add_flextable(tbl4_ft) %>% body_add_par(" ") %>%
  body_add_flextable(ancova_ft) %>% body_add_par(" ") %>%
  body_add_flextable(logit_ft) %>% body_add_par(" ") %>%
  body_add_flextable(logit_epv_ft)
print(doc, target = file.path(OUT_DIR, "ESJ_Tables_FULL_per1.docx"))

ppt <- read_pptx()
if (file.exists(file.path(FIG_DIR, "Fig_ODI_RCS.tiff"))) {
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme") %>%
    ph_with("RCS for ODI", location = ph_location_type(type = "title")) %>%
    ph_with(external_img(file.path(FIG_DIR, "Fig_ODI_RCS.tiff"), width = 7, height = 4), location = ph_location_type(type = "body"))
}
if (file.exists(file.path(FIG_DIR, "Fig_JOA_RCS.tiff"))) {
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme") %>%
    ph_with("RCS for JOA", location = ph_location_type(type = "title")) %>%
    ph_with(external_img(file.path(FIG_DIR, "Fig_JOA_RCS.tiff"), width = 7, height = 4), location = ph_location_type(type = "body"))
}
print(ppt, target = file.path(OUT_DIR, "ESJ_Figs.pptx"))

## 控制台输出 -----------------------------------------------------------

if (nrow(ancova_res) > 0) {
  ancova_res %>% mutate(beta_ci = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)) %>%
    select(Outcome, beta_ci, p.value) %>%
    {print("曝光对 1 年结局的回归系数 (per 1 point)"); print(.)}
}

cat(sprintf("ODI_1y range: %s\n", paste(range(d0$ODI_1y, na.rm = TRUE), collapse = "-")))
cat(sprintf("JOA_1y range: %s\n", paste(range(d0$JOA_1y, na.rm = TRUE), collapse = "-")))

## 脚本结束 -------------------------------------------------------------
