
# ESJ pipeline (simplified, modular)
# 功能：读数 → 基线/结局表 → 一致性 → 趋势 → ANCOVA(HC3) → RCS → MCID → 功效/EPV → 导出
# 关键定制：Table2 中 VAS 始终中位数(IQR)；RCS x 轴 = "AAC (Kauppila 0–24)"；生成 A/B 合并图

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(gtsummary)
  library(flextable)
  library(officer)
  library(broom)
  library(sandwich)
  library(lmtest)
  library(rms)
  library(irr)
  library(clinfun)
  library(logistf)
  library(brglm2)
  library(pwr)
  library(patchwork)
})

# ------------------------------ 基础配置 ------------------------------
DATA_PATH <- "C:/Users/wss29/Desktop/ESJ-2.csv"
OUT_DIR   <- "pub_results"
FIG_DIR   <- "pub_figs"
walk(list(OUT_DIR, FIG_DIR), dir.create, showWarnings = FALSE, recursive = TRUE)

read_try <- function(p){
  stopifnot(file.exists(p))
  tryCatch(
    read_csv(p, locale = locale(encoding = "UTF-8"), show_col_types = FALSE),
    error = function(e) read_csv(p, locale = locale(encoding = "GB18030"), show_col_types = FALSE)
  )
}
raw <- read_try(DATA_PATH)
cn  <- names(raw)

# ------------------------------ 列名匹配 ------------------------------
lookup_cols <- list(
  Kauppilla = c("^Kauppilla$", "^Kaupilla$", "aac.*24", "Kauppila"),
  Kp1       = c("^Kauppilla-?1$", "^Kaupilla-?1$", "Kaupilla-1", "Kp1"),
  Kp2       = c("^Kauppilla-?2$", "^Kaupilla-?2$", "Kaupilla-2"),
  Fusion1   = c("^Fusion[ _]?Status1$", "Fusion Status1", "^Fusion$", "Fusion1"),
  Fusion2   = c("^Fusion[ _]?Status2$", "Fusion Status2", "Fusion2"),
  JOA_pre   = c("^pre-?JOA$", "pre-JOA", "pre_JOA"),
  JOA_1y    = c("^JOA$", "JOA_1y", "JOA.*1"),
  ODI_pre   = c("^pre-?ODI$", "pre-ODI", "pre_ODI"),
  ODI_1y    = c("^ODI$", "ODI_1y", "ODI.*1"),
  YVAS      = c("^Yvas$", "^YVAS$", "Back.*VAS"),
  TVAS      = c("^Tvas$", "^TVAS$", "Leg.*VAS"),
  VAS       = c("^VAS$"),
  Gender    = c("^Gender$", "Sex"),
  Age       = c("^Age$"),
  BMI       = c("^BMI$"),
  HP        = c("^HP$", "hyper.*lip"),
  HT        = c("^Ht$", "^HT$", "hyper.*tens"),
  Diabetes  = c("^Diabetes$"),
  AHUS      = c("^AHUS$", "^HU.?OP", "Osteo.*HU"),
  Mac       = c("^Mac$", "Macnab"),
  MF        = c("^MF$", "multi.*(level|segment)"),
  OpTime    = c("^Surgery duration$", "op.*time", "Duration"),
  BloodLoss = c("^Blood loss$", "EBL"),
  Tobacco   = c("^Tobacco use$", "Tobacco$", "smok.*former")
)

pick_col <- function(patterns){
  for (pt in patterns){
    hit <- which(tolower(cn) == tolower(pt))
    if (length(hit)) return(cn[hit[1]])
    fuzzy <- grep(pt, cn, ignore.case = TRUE, value = TRUE)
    if (length(fuzzy)) return(fuzzy[1])
  }
  NA_character_
}
col_map <- purrr::map_chr(lookup_cols, pick_col)

numify <- function(x) suppressWarnings(as.numeric(x))

# ------------------------------ 数据整理 ------------------------------
d0 <- raw %>%
  transmute(
    Kauppilla_raw = if (!is.na(col_map["Kauppilla"])) .data[[col_map["Kauppilla"]]] else NA,
    Kp1           = if (!is.na(col_map["Kp1"])) numify(.data[[col_map["Kp1"]]]) else NA_real_,
    Kp2           = if (!is.na(col_map["Kp2"])) numify(.data[[col_map["Kp2"]]]) else NA_real_,
    Fusion1       = if (!is.na(col_map["Fusion1"])) as.integer(.data[[col_map["Fusion1"]]]) else NA_integer_,
    Fusion2       = if (!is.na(col_map["Fusion2"])) as.integer(.data[[col_map["Fusion2"]]]) else NA_integer_,
    JOA_pre       = numify(.data[[col_map["JOA_pre"]]]),
    JOA_1y        = numify(.data[[col_map["JOA_1y"]]]),
    ODI_pre       = numify(.data[[col_map["ODI_pre"]]]),
    ODI_1y        = numify(.data[[col_map["ODI_1y"]]]),
    YVAS          = if (!is.na(col_map["YVAS"])) numify(.data[[col_map["YVAS"]]]) else NA_real_,
    TVAS          = if (!is.na(col_map["TVAS"])) numify(.data[[col_map["TVAS"]]]) else NA_real_,
    VAS_total     = if (!is.na(col_map["VAS"])) numify(.data[[col_map["VAS"]]]) else (YVAS + TVAS),
    Gender        = .data[[col_map["Gender"]]],
    Age           = numify(.data[[col_map["Age"]]]),
    BMI           = numify(.data[[col_map["BMI"]]]),
    HP            = as.integer(.data[[col_map["HP"]]]),
    HT            = as.integer(.data[[col_map["HT"]]]),
    Diabetes      = as.integer(.data[[col_map["Diabetes"]]]),
    AHUS          = as.integer(.data[[col_map["AHUS"]]]),
    Mac           = .data[[col_map["Mac"]]],
    MF            = as.integer(.data[[col_map["MF"]]]),
    OpTime        = numify(.data[[col_map["OpTime"]]]),
    BloodLoss     = numify(.data[[col_map["BloodLoss"]]]),
    Tobacco01     = as.integer(.data[[col_map["Tobacco"]]])
  ) %>%
  mutate(
    Fusion      = if_else(Fusion1 %in% c(1, "1"), 1L, if_else(is.na(Fusion1), NA_integer_, 0L)),
    Gender01    = if_else(Gender %in% c(1, "1", "Male", "male", "M", "m"), 1L, if_else(is.na(Gender), NA_integer_, 0L)),
    Smoking01   = factor(Tobacco01, levels = c(0,1), labels = c("0","1")),
    Kauppilla   = numify(Kauppilla_raw),
    AAC_score   = dplyr::coalesce(Kp1, Kauppilla),
    AAC01       = factor(if_else(AAC_score >= 1, 1, 0), levels = c(0,1)),
    AAC_cat4    = cut(AAC_score, breaks = c(-Inf,0,4,8,Inf), labels = c("0","1–4","5–8","≥9"), ordered_result = TRUE),
    dJOA        = JOA_1y - JOA_pre,
    dODI        = ODI_pre - ODI_1y,
    MacLab      = recode(as.character(Mac), `1`="Excellent", `2`="Good", `3`="Fair", `4`="Poor", .default = NA_character_) %>% factor(levels = c("Excellent","Good","Fair","Poor")),
    MCID10      = as.integer(dODI >= 10),
    MCID30      = as.integer(dODI >= 0.30 * ODI_pre)
  )

range_report <- function(x) paste(range(x, na.rm = TRUE), collapse = " ~ ")
cat("\n[AUDIT] Kp1 range:", range_report(d0$Kp1), "| Kauppilla range:", range_report(d0$Kauppilla), "\n")
if (is.finite(max(d0$Kauppilla, na.rm = TRUE)) && max(d0$Kauppilla, na.rm = TRUE) <= 4)
  warning("Kauppilla 上界 ≤4，疑似分类分档")

expo_var <- if (is.finite(max(d0$Kp1, na.rm = TRUE)) && max(d0$Kp1, na.rm = TRUE) > 4) "Kp1"
else if (is.finite(max(d0$Kauppilla, na.rm = TRUE)) && max(d0$Kauppilla, na.rm = TRUE) > 4) "Kauppilla"
else stop("缺少可用 AAC-24 总分")
cat("[INFO] ANCOVA 暴露:", expo_var, "(per 1 point)\n")

# ------------------------------ 常用工具 ------------------------------
pfun3 <- function(x) gtsummary::style_pvalue(x, digits = 3)
`%||%` <- function(x, y) if (!is.null(x) && length(x)) x else y
fmt_median_iqr <- function(x){
  x <- x[is.finite(x)]; if (!length(x)) return("NA")
  sprintf("%.2f (%.2f, %.2f)", median(x), quantile(x, 0.25, type = 2), quantile(x, 0.75, type = 2))
}
three_line <- function(ft){
  std <- fp_border(color = "black", width = 1)
  ft <- theme_vanilla(ft) %>% border_remove()
  nh <- nrow_part(ft, "header"); nb <- nrow_part(ft, "body"); nc <- length(ft$col_keys)
  if (!is.na(nh) && nh>0){
    ft <- border(ft, i=1, j=1:nc, border.top=std, part="header")
    ft <- border(ft, i=nh, j=1:nc, border.bottom=std, part="header")
  }
  if (!is.na(nb) && nb>0) ft <- border(ft, i=nb, j=1:nc, border.bottom=std, part="body")
  font(ft, fontname = "Times New Roman", part = "all") %>% fontsize(size = 9, part = "all")
}
compose_p_with_sup <- function(ft, tbl, sup_vec){
  p_col <- intersect(ft$col_keys, "p.value"); if (!length(p_col)) return(ft)
  p_fmt <- gtsummary::style_pvalue(tbl$table_body$p.value, digits = 3)
  FONT <- "Times New Roman"; SIZE <- 9
  fp_norm <- officer::fp_text(font.family = FONT, font.size = SIZE)
  fp_sup  <- officer::fp_text(font.family = FONT, font.size = SIZE, vertical.align = "superscript")
  rows <- which(!is.na(tbl$table_body$p.value))
  for (i in rows){
    base <- p_fmt[i]; s <- sup_vec[i]
    ft <- flextable::compose(ft, i = i, j = p_col,
      value = flextable::as_paragraph(
        flextable::as_chunk(base %||% "", props = fp_norm),
        if (!is.na(s) && nzchar(s)) flextable::as_chunk(s, props = fp_sup)
      ))
  }
  ft
}

# ------------------------------ 正态性检验 ------------------------------
is_sw_ok <- function(x){
  x <- x[is.finite(x)]; if (length(x) < 3) return(NA_real_)
  if (length(x) > 5000) x <- sample(x, 5000)
  suppressWarnings(shapiro.test(x)$p.value)
}
write_csv(
  tibble(var = c("Age","BMI","OpTime","BloodLoss","JOA_1y","ODI_1y")) %>%
    mutate(shapiro_p = map_dbl(var, ~ is_sw_ok(d0[[.x]]))),
  file.path(OUT_DIR,"shapiro_wilk_tests.csv")
)

# ------------------------------ Table 1 ------------------------------
cont_t1 <- c("Age","BMI","OpTime","BloodLoss")
cat_t1  <- c("Gender01","HT","HP","Diabetes","AHUS","MF","Smoking01")

sup1 <- "¹"; sup2 <- "²"; sup3 <- "³"

tbl1 <- d0 %>%
  select(AAC01, Age, BMI, Gender01, HT, HP, Diabetes, AHUS, MF, OpTime, BloodLoss, Smoking01) %>%
  tbl_summary(
    by = AAC01,
    statistic = list(BMI ~ "{mean} ± {sd}", all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} / {N} ({p}%)"),
    digits = list(BMI ~ 2, all_continuous() ~ 2),
    label = list(BMI ~ "BMI (kg/m^2)", Gender01 ~ "Male (1 vs 0)", HT ~ "Hypertension", HP ~ "Hyperlipidemia", AHUS ~ "Osteoporosis (HU-based)", MF ~ "Multilevel fusion (≥2)", Smoking01 ~ "Tobacco use (former vs never)")
  ) %>%
  add_n() %>%
  add_p(test = list(BMI ~ "t.test", all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"), test.args = list(all_categorical() ~ list(simulate.p.value = TRUE)), pvalue_fun = pfun3) %>%
  modify_fmt_fun(columns = dplyr::all_of("p.value"), fun = identity)

sup_vec_t1 <- tbl1$table_body %>% mutate(.sup = case_when(variable == "BMI" ~ sup1, variable %in% cont_t1 ~ sup2, variable %in% cat_t1 ~ sup3, TRUE ~ "")) %>% pull(.sup)

ft1 <- as_flex_table(tbl1) %>% three_line() %>% compose_p_with_sup(tbl1, sup_vec_t1) %>%
  add_footer_lines(values = c(
    "Notes: ¹ Student’s t-test for BMI (mean ± SD); ² Wilcoxon rank–sum for other continuous variables (median [Q1,Q3]); ³ Pearson’s chi-square (Monte Carlo when sparse).",
    "Tobacco use: 0 = never, 1 = former. AAC = abdominal aortic calcification; HU = Hounsfield unit."
  )) %>% fontsize(size = 8, part = "footer")

# ------------------------------ Table 2 ------------------------------
cont_t2 <- c("JOA_pre","JOA_1y","ODI_pre","ODI_1y","YVAS","TVAS","VAS_total")
cat_t2  <- c("MacLab","Fusion","MCID10","MCID30","Smoking01")

force_median_iqr <- function(tbl, data, vars){
  tb <- tbl$table_body; stat_cols <- grep("^stat_", names(tb), value = TRUE); lev_by <- levels(data$AAC01)
  for (v in vars){
    idx <- which(tb$variable == v & tb$row_type == "label"); if (!length(idx)) next
    repl <- data %>% filter(!is.na(AAC01)) %>% group_by(AAC01) %>% summarise(txt = fmt_median_iqr(.data[[v]]), .groups = "drop")
    vals <- repl$txt[match(lev_by, repl$AAC01)]
    for (k in seq_along(stat_cols)) if (!is.na(vals[k])) tb[[stat_cols[k]]][idx] <- vals[k]
  }
  tbl$table_body <- tb; tbl
}

tbl2 <- d0 %>% mutate(across(c(YVAS,TVAS,VAS_total), numify)) %>%
  select(AAC01, JOA_pre, JOA_1y, ODI_pre, ODI_1y, YVAS, TVAS, VAS_total, MacLab, Fusion, MCID10, MCID30, Smoking01) %>%
  tbl_summary(
    by = AAC01,
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} / {N} ({p}%)"),
    digits = all_continuous() ~ 2,
    label = list(JOA_pre ~ "Preoperative JOA", JOA_1y ~ "JOA at 1 year", ODI_pre ~ "Preoperative ODI", ODI_1y ~ "ODI at 1 year", YVAS ~ "Back pain VAS", TVAS ~ "Leg pain VAS", VAS_total ~ "Total VAS (back + leg)", MacLab ~ "Modified MacNab", Fusion ~ "Fusion (1 vs 0)", MCID10 ~ "MCID ≥10 points", MCID30 ~ "MCID ≥30%", Smoking01 ~ "Tobacco use (former vs never)")
  ) %>%
  add_n() %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"), test.args = list(all_categorical() ~ list(simulate.p.value = TRUE)), pvalue_fun = pfun3) %>%
  modify_fmt_fun(columns = dplyr::all_of("p.value"), fun = identity) %>%
  force_median_iqr(data = d0, vars = c("YVAS","TVAS","VAS_total"))

sup_vec_t2 <- tbl2$table_body %>% mutate(.sup = case_when(variable %in% cont_t2 ~ sup2, variable %in% cat_t2 ~ sup3, TRUE ~ "")) %>% pull(.sup)
ft2 <- as_flex_table(tbl2) %>% three_line() %>% compose_p_with_sup(tbl2, sup_vec_t2)

# ------------------------------ Table 3：一致性 ------------------------------
try_icc <- function(dat) tryCatch(irr::icc(dat, model="twoway", type="agreement", unit="single", conf.level=0.95), error = function(e) NULL)
try_kappa <- function(dat) tryCatch(irr::kappa2(dat), error = function(e) NULL)

icc_res <- if (sum(complete.cases(d0[,c("Kp1","Kp2")])) >= 5) try_icc(as.matrix(d0[,c("Kp1","Kp2")])) else NULL
kappa_obj <- if (sum(complete.cases(d0[,c("Fusion1","Fusion2")])) >= 5) try_kappa(d0[,c("Fusion1","Fusion2")]) else NULL
kap_ci <- if (!is.null(kappa_obj)) {set.seed(2025); df <- d0 %>% select(Fusion1, Fusion2) %>% drop_na(); n <- nrow(df); boots <- replicate(2000, irr::kappa2(df[sample.int(n, n, replace = TRUE), ])$value); quantile(boots, c(0.025,0.975), na.rm = TRUE)} else c(NA_real_, NA_real_)

tab_reliab <- tibble(
  Measure = c("AAC-24 (Kp1 vs Kp2)", "Fusion (Rater1 vs Rater2)"),
  Method  = c("ICC(A,1)", "Cohen's κ"),
  Estimate= c(ifelse(is.null(icc_res), NA, icc_res$value), ifelse(is.null(kappa_obj), NA, kappa_obj$value)),
  `95% CI`= c(ifelse(is.null(icc_res), NA_character_, sprintf("(%.3f, %.3f)", icc_res$lbound, icc_res$ubound)), sprintf("(%.3f, %.3f)", kap_ci[1], kap_ci[2])),
  p       = c(ifelse(is.null(icc_res), NA, icc_res$p.value), ifelse(is.null(kappa_obj), NA, kappa_obj$p.value))
) %>% mutate(Estimate = ifelse(is.na(Estimate), NA, sprintf("%.3f", Estimate)), p = case_when(is.na(p) ~ NA_character_, p < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p)))
ft3 <- flextable(tab_reliab) %>% three_line()

# ------------------------------ Table 4：趋势检验 ------------------------------
trend_test <- function(y, g){
  ok <- is.finite(y) & is.finite(as.numeric(g)); if (sum(ok) < 3) return(list(stat=NA, p=NA, n=sum(ok)))
  jt <- tryCatch(clinfun::jonckheere.test(y[ok], as.numeric(g[ok]), alternative="two.sided"), error = function(e) NULL)
  if (is.null(jt)) return(list(stat=NA, p=NA, n=sum(ok)))
  list(stat=as.numeric(jt$statistic), p=as.numeric(jt$p.value), n=sum(ok))
}

tt1 <- trend_test(d0$JOA_1y, d0$AAC_cat4); tt2 <- trend_test(d0$ODI_1y, d0$AAC_cat4)
ft4 <- flextable(tibble(Outcome=c("JOA (1-year)", "ODI (1-year)"), Statistic=c(tt1$stat, tt2$stat), `p (trend)`=c(formatC(tt1$p, format="f", digits=4), formatC(tt2$p, format="f", digits=4)), n=c(tt1$n, tt2$n))) %>% three_line()

# ------------------------------ ANCOVA (HC3) ------------------------------
complete_n <- function(outcome, preds, data){
  frm <- as.formula(paste(outcome, "~", paste(preds, collapse=" + ")))
  sum(complete.cases(model.frame(frm, data = data, na.action = na.pass)))
}
prune_vars <- function(outcome, preds, data, min_n = 10){
  while (complete_n(outcome, preds, data) < min_n && length(preds) > 2){
    miss_rate <- sapply(preds, function(v) mean(is.na(data[[v]])))
    keep <- c(expo_var, if (outcome=="ODI_1y") "ODI_pre" else "JOA_pre")
    cand <- setdiff(preds, keep); if (!length(cand)) break
    preds <- setdiff(preds, cand[which.max(miss_rate[cand])])
  }
  preds
}

cand_full <- c(expo_var,"ODI_pre","Fusion","Gender01","Age","BMI","HT","HP","Diabetes","AHUS","MF","OpTime","BloodLoss","Tobacco01")
cand_min  <- c(expo_var,"ODI_pre","Age","Gender01","BMI","Tobacco01")

preds_odi_full <- prune_vars("ODI_1y", cand_full, d0)
preds_joa_full <- prune_vars("JOA_1y", setdiff(cand_full,"ODI_pre") %>% c("JOA_pre"), d0)
preds_odi_min  <- prune_vars("ODI_1y", cand_min, d0)
preds_joa_min  <- prune_vars("JOA_1y", setdiff(cand_min,"ODI_pre") %>% c("JOA_pre"), d0)

mk_form <- function(outcome, preds) as.formula(paste(outcome, "~", paste(preds, collapse = " + ")))
fit_safe <- function(formula, data){
  stopifnot(complete_n(all.vars(formula)[1], attr(terms(formula), "term.labels"), data) > 0)
  lm(formula, data = data, na.action = na.omit)
}

fit_odi_full <- fit_safe(mk_form("ODI_1y", preds_odi_full), d0)
fit_joa_full <- fit_safe(mk_form("JOA_1y", preds_joa_full), d0)
fit_odi_min  <- fit_safe(mk_form("ODI_1y", preds_odi_min), d0)
fit_joa_min  <- fit_safe(mk_form("JOA_1y", preds_joa_min), d0)

tidy_hc3 <- function(fit){
  vc <- sandwich::vcovHC(fit, type="HC3"); se <- sqrt(diag(vc)); co <- coef(fit)
  tibble(term = names(co), estimate = co, std.error = se, conf.low = co - qnorm(0.975)*se, conf.high = co + qnorm(0.975)*se, p.value = 2*pnorm(abs(co/se), lower.tail = FALSE))
}
var_label_lm <- c("(Intercept)"="Intercept", "Kauppilla"="AAC-24 (per 1 point)", "Kp1"="AAC-24 (per 1 point)", "ODI_pre"="Preoperative ODI", "JOA_pre"="Preoperative JOA", "Fusion"="Fusion (1 vs 0)", "Gender01"="Male (1 vs 0)", "Age"="Age", "BMI"="BMI", "HT"="Hypertension", "HP"="Hyperlipidemia", "Diabetes"="Diabetes", "AHUS"="Osteoporosis", "MF"="Multilevel fusion", "OpTime"="Surgery duration", "BloodLoss"="Blood loss", "Tobacco01"="Tobacco use")
mk_lm_table <- function(fit, title){
  tt <- tidy_hc3(fit) %>% filter(term != "(Intercept)")
  lab <- var_label_lm[tt$term]; lab[is.na(lab)] <- tt$term
  tb <- tt %>% transmute(Variable = lab, Beta = sprintf("%.2f", estimate), `Std. Error`=sprintf("%.2f", std.error), `95% CI` = paste0("(", sprintf("%.2f", conf.low), ", ", sprintf("%.2f", conf.high), ")"), p = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))
  flextable(tb) %>% three_line() %>% set_caption(title)
}

ft5_full <- mk_lm_table(fit_joa_full, "Table 5. Linear regression for 1-year JOA (HC3)")
ft6_full <- mk_lm_table(fit_odi_full, "Table 6. Linear regression for 1-year ODI (HC3)")
ft5_min  <- mk_lm_table(fit_joa_min,  "Table S1. Minimal model for 1-year JOA")
ft6_min  <- mk_lm_table(fit_odi_min,  "Table S2. Minimal model for 1-year ODI")

# ------------------------------ RCS ------------------------------
dd <- datadist(d0); options(datadist = "dd")
if (!is.finite(max(d0$Kp1, na.rm = TRUE))) stop("RCS 需要 Kp1")

pick_nk <- function(nu) dplyr::case_when(nu >= 4 ~ 4L, nu == 3 ~ 3L, TRUE ~ NA_integer_)
nk_odi <- pick_nk(dplyr::n_distinct(na.omit(d0$Kp1[d0$ODI_1y %>% is.finite()])))
nk_joa <- pick_nk(dplyr::n_distinct(na.omit(d0$Kp1[d0$JOA_1y %>% is.finite()])))
mk_rcs_term <- function(nk) if (is.na(nk)) "Kp1" else sprintf("rcs(Kp1,%d)", nk)

f_rcs_odi <- as.formula(sprintf("ODI_1y ~ %s + ODI_pre + Fusion + Gender01 + Age + BMI + HT + HP + Diabetes + AHUS + MF + OpTime + BloodLoss + Tobacco01", mk_rcs_term(nk_odi)))
f_rcs_joa <- as.formula(sprintf("JOA_1y ~ %s + JOA_pre + Fusion + Gender01 + Age + BMI + HT + HP + Diabetes + AHUS + MF + OpTime + BloodLoss + Tobacco01", mk_rcs_term(nk_joa)))
ols_odi <- rms::ols(f_rcs_odi, data = d0, x = TRUE, y = TRUE, na.action = na.omit)
ols_joa <- rms::ols(f_rcs_joa, data = d0, x = TRUE, y = TRUE, na.action = na.omit)

build_effect_plot <- function(fit, data, ylab){
  xr <- range(data$Kp1, na.rm = TRUE); grid <- list(Kp1 = seq(xr[1], xr[2], length.out = 200))
  pr <- as.data.frame(do.call(rms::Predict, c(list(fit), grid, list(conf.int = 0.95))))
  pr <- pr %>% filter(is.finite(Kp1), is.finite(yhat), is.finite(lower), is.finite(upper))
  ggplot(pr, aes(Kp1, yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey90") +
    geom_line(linewidth = 0.7) +
    labs(x = "AAC (Kauppila 0–24)", y = ylab) +
    theme_classic(base_size = 9) +
    theme(text = element_text(family = "Arial"), legend.position = "none") +
    coord_cartesian(clip = "off")
}

p_odi <- build_effect_plot(ols_odi, d0, "Adjusted 1-year ODI")
p_joa <- build_effect_plot(ols_joa, d0, "Adjusted 1-year JOA")
ggsave(file.path(FIG_DIR,"Fig_ODI_RCS.tiff"), plot=p_odi, width=84, height=60, units="mm", dpi=600, device="tiff", compression="lzw", limitsize=FALSE)
ggsave(file.path(FIG_DIR,"Fig_JOA_RCS.tiff"), plot=p_joa, width=84, height=60, units="mm", dpi=600, device="tiff", compression="lzw", limitsize=FALSE)
p_combined <- (p_odi + labs(tag = "A")) + (p_joa + labs(tag = "B")) + plot_layout(ncol = 2)
ggsave(filename = file.path(FIG_DIR, "Fig_RCS_ODI_JOA.tiff"), plot = p_combined, width = 170, height = 80, units = "mm", dpi = 600, device = "tiff", compression = "lzw", limitsize = FALSE)

# ------------------------------ 功效与 MCID 逻辑回归 ------------------------------
power_note <- function(f_full, term_to_drop, label){
  mf_cc <- model.frame(f_full); f_full0 <- lm(formula(f_full), data = mf_cc)
  f_min <- lm(update(formula(f_full), paste(". ~ . -", term_to_drop)), data = mf_cc)
  r2_diff <- summary(f_full0)$r.squared - summary(f_min)$r.squared
  f2 <- max(r2_diff / (1 - summary(f_full0)$r.squared), 0)
  pw <- tryCatch(pwr::pwr.f2.test(u = 1, v = df.residual(f_full0), f2 = f2)$power, error = function(e) NA_real_)
  paste0(label, ": ΔR2=", sprintf("%.4f", r2_diff), ", f^2=", sprintf("%.4f", f2), if (is.finite(pw)) paste0(", power≈", sprintf("%.3f", pw)) else "")
}

pw_txt <- c(power_note(f_full = fit_odi_full, term_to_drop = expo_var, label = "ODI ANCOVA (AAC effect)"), power_note(f_full = fit_joa_full, term_to_drop = expo_var, label = "JOA ANCOVA (AAC effect)"))
writeLines(pw_txt, file.path(OUT_DIR,"posthoc_power_linear.txt"))

mk_logit_forms <- function(outcome){
  list(full = as.formula(paste(outcome, "~", paste(c(expo_var,"ODI_pre","Fusion","Gender01","Age","BMI","HT","HP","Diabetes","AHUS","MF","OpTime","BloodLoss","Tobacco01"), collapse = " + ")))),
       min  = as.formula(paste(outcome, "~", paste(c(expo_var,"ODI_pre","Age","Gender01","BMI","Tobacco01"), collapse = " + ")))))
}
form_mcid10 <- mk_logit_forms("MCID10"); form_mcid30 <- mk_logit_forms("MCID30")

fit_firth_safe <- function(formula, data){
  dat <- model.frame(formula, data = data) %>% drop_na()
  tryCatch(logistf::logistf(formula, data = dat), error = function(e) NULL) %||%
    glm(formula, data = dat, family = binomial(), method = brglm2::brglmFit, type = "AS_mixed")
}

f_mcid10_full <- fit_firth_safe(form_mcid10$full, d0); f_mcid10_min <- fit_firth_safe(form_mcid10$min, d0)
f_mcid30_full <- fit_firth_safe(form_mcid30$full, d0); f_mcid30_min <- fit_firth_safe(form_mcid30$min, d0)

epv_calc <- function(y, p) {ev <- sum(y==1, na.rm = TRUE); nev <- sum(y==0, na.rm = TRUE); c(events = ev, nonevents = nev, params = p, EPV = ifelse(p>0, min(ev,nev)/p, NA))}
epv_out <- rbind(MCID10_full = epv_calc(d0$MCID10, length(attr(terms(form_mcid10$full), "term.labels"))), MCID10_min = epv_calc(d0$MCID10, length(attr(terms(form_mcid10$min), "term.labels"))), MCID30_full = epv_calc(d0$MCID30, length(attr(terms(form_mcid30$full), "term.labels"))), MCID30_min = epv_calc(d0$MCID30, length(attr(terms(form_mcid30$min), "term.labels")))) %>% as.data.frame()
write.csv(epv_out, file.path(OUT_DIR, "EPV_logistic_models.csv"))

tidy_logistf <- function(fit){
  if (inherits(fit, "logistf")) return(broom::tidy(fit, conf.int = TRUE))
  broom::tidy(fit, conf.int = TRUE)
}

show_or <- function(fit, title, epv_label){
  tt <- tidy_logistf(fit) %>% filter(term != "(Intercept)") %>%
    mutate(OR = exp(estimate), lo = exp(conf.low), hi = exp(conf.high), p = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value))) %>%
    transmute(Variable = recode(term, "Kauppilla"="AAC-24 (per 1 point)", "Kp1"="AAC-24 (per 1 point)", "ODI_pre"="Preoperative ODI", "Fusion"="Fusion (1 vs 0)", "Gender01"="Male (1 vs 0)", "Age"="Age", "BMI"="BMI", "HT"="Hypertension", "HP"="Hyperlipidemia", "Diabetes"="Diabetes", "AHUS"="Osteoporosis", "MF"="Multilevel fusion", "OpTime"="Surgery duration", "BloodLoss"="Blood loss", "Tobacco01"="Tobacco use", .default = term), OR = sprintf("%.2f", OR), `95% CI` = paste0("(", sprintf("%.2f", lo), ", ", sprintf("%.2f", hi), ")"), p)
  flextable(tt) %>% set_caption(title) %>% three_line() %>% add_footer_lines(values = paste0("EPV: ", sprintf("%.2f", epv_label)))
}

ft_mcid10_or_full <- show_or(f_mcid10_full, "Table S5. Logistic regression for MCID ≥10 (full)", epv_out["MCID10_full","EPV"])
ft_mcid10_or_min  <- show_or(f_mcid10_min,  "Table S6. Logistic regression for MCID ≥10 (minimal)", epv_out["MCID10_min","EPV"])
ft_mcid30_or_full <- show_or(f_mcid30_full, "Table S7. Logistic regression for MCID ≥30% (full)", epv_out["MCID30_full","EPV"])
ft_mcid30_or_min  <- show_or(f_mcid30_min,  "Table S8. Logistic regression for MCID ≥30% (minimal)", epv_out["MCID30_min","EPV"])

# ------------------------------ 导出 Word & PPT ------------------------------
doc <- read_docx() %>%
  body_add_par("Tables – ESJ style (per 1-point AAC effect)", style = "heading 1") %>%
  body_add_par(paste0("Generated on: ", Sys.Date()), style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("Table 1. Baseline characteristics by AAC presence", style = "heading 2") %>% body_add_flextable(ft1) %>%
  body_add_par("Table 2. One-year outcomes by AAC presence", style = "heading 2") %>% body_add_flextable(ft2) %>%
  body_add_par("Table 3. Reliability (ICC/κ)", style = "heading 2") %>% body_add_flextable(ft3) %>%
  body_add_par("Table 4. Trend across AAC categories", style = "heading 2") %>% body_add_flextable(ft4) %>%
  body_add_par("Table 5. Linear regression for 1-year JOA", style = "heading 2") %>% body_add_flextable(ft5_full) %>%
  body_add_par("Table 6. Linear regression for 1-year ODI", style = "heading 2") %>% body_add_flextable(ft6_full) %>%
  body_add_par("Table S1–S2. Minimal ANCOVA models", style = "heading 2") %>% body_add_flextable(ft5_min) %>% body_add_flextable(ft6_min) %>%
  body_add_par("Table S5–S8. Logistic models for MCID", style = "heading 2") %>% body_add_flextable(ft_mcid10_or_full) %>% body_add_flextable(ft_mcid10_or_min) %>% body_add_flextable(ft_mcid30_or_full) %>% body_add_flextable(ft_mcid30_or_min) %>%
  body_add_par("Figure Legends", style = "heading 2") %>% body_add_par("Figure legends to be added manually.", style = "Normal")

print(doc, target = file.path(OUT_DIR,"ESJ_Tables_FULL_per1.docx"))
ppt <- read_pptx() %>%
  add_slide(layout="Blank", master="Office Theme") %>% ph_with(external_img(file.path(FIG_DIR,"Fig_ODI_RCS.tiff")), location = ph_location_fullsize()) %>%
  add_slide(layout="Blank", master="Office Theme") %>% ph_with(external_img(file.path(FIG_DIR,"Fig_JOA_RCS.tiff")), location = ph_location_fullsize())
print(ppt, target = file.path(OUT_DIR, "ESJ_Figs.pptx"))

# ------------------------------ 控制台核对 ------------------------------
summ_expo <- function(fit, expo){
  tb <- tidy_hc3(fit) %>% filter(term == expo)
  if (nrow(tb)) cat(sprintf("\n[%s] β=%.2f (95%%CI %.2f, %.2f) | p=%s\n",
                          expo, tb$estimate, tb$conf.low, tb$conf.high,
                          ifelse(tb$p.value < 0.001, "<0.001", sprintf("%.3f", tb$p.value))))
}

summ_expo(fit_odi_full, expo_var)
summ_expo(fit_joa_full, expo_var)
cat(sprintf("\n[Range] ODI_1y: %s ~ %s ; JOA_1y: %s ~ %s\n",
            signif(range(d0$ODI_1y, na.rm=TRUE)),
            signif(range(d0$JOA_1y, na.rm=TRUE))))
