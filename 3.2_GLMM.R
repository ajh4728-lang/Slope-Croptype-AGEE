# ============================================================
# Gamma(log) GLMM: crop main effect adjusted for rainfall and slope
# 작물유형별 TN/TP 부하량 차이 분석
#
# 핵심 목적:
#  - 상호작용 없이, 강우량과 경사도를 통제한 상태에서
#    작물유형별 TN/TP 부하량 차이가 유의한지 검정
#
# 최종 해석 모델:
#  - TN ~ crop + ns(rainfall, df = 3) + slope + (1 | study)
#  - TP ~ crop + ns(rainfall, df = 3) + slope + (1 | study)
#
# 이 코드는 상호작용을 검정하지 않음.
# 따라서 해석은 다음으로 제한:
#  - "강우량과 경사도를 보정한 후에도 작물유형별 TN/TP 부하량 차이가 있는가?"
#
# 말하면 안 되는 해석:
#  - "작물유형에 따라 경사도 효과가 다르다"
#  - "작물유형과 경사도의 상호작용이 유의하다"
#
# 반영 사항:
#  - Rainfall <= 200 mm/event
#  - SlopeGradient < 30%
#  - 엑셀의 "-" 결측값 제거
#  - TN, TP, Rainfall, SlopeGradient 숫자형 변환
#  - Gamma(log) GLMM
#  - 상호작용항 없음
#  - study random effect: (1 | study)
#  - 작물 주효과 LRT
#  - 작물별 보정 추정평균 EMM
#  - Tukey 사후검정
#  - CLD 문자: 추정평균이 높은 작물부터 a, b, c...
#  - 과분산 대응용 dispformula = ~ crop 모델 함께 적합
#  - DHARMa 진단
#  - 엑셀 및 그림 저장
#  - 그림 색상: crop 값을 Cc/Fc/Lv/Lg/Oc/Pc/Rv 코드로 변환 후 수동 색상 적용
#  - 그림 제목 제거
#  - 원자료 점 제거
#  - 막대 검은 테두리
#  - 그래프 panel border로 상단/우측 포함 네모 테두리 표시
#  - Y축 숫자 소수점 1자리 표시
#  - TN 그래프 Y축 0.0–8.0 고정
# ============================================================


# ------------------------------------------------------------
# 0. 패키지 불러오기
# ------------------------------------------------------------

# 필요시 설치:
# install.packages(c("glmmTMB", "emmeans", "multcompView", "DHARMa",
#                    "performance", "openxlsx", "ggplot2", "dplyr",
#                    "readxl", "splines"))

library(glmmTMB)
library(emmeans)
library(multcompView)
library(DHARMa)
library(performance)
library(ggplot2)
library(dplyr)
library(readxl)
library(openxlsx)
library(splines)


# ------------------------------------------------------------
# 1. 저장 폴더 설정
# ------------------------------------------------------------

out_dir <- "C:/Users/ajh/OneDrive/00_수생태환경연구실/01_연구논문관련/1_논문작성/1_양분유출에 취약한 밭 환경조건 분석/Data/GLMM분석"


# ------------------------------------------------------------
# 2. 데이터 불러오기
# ------------------------------------------------------------

df_raw <- read_excel(
  file.path(out_dir, "선행연구 데이터 - 플롯 나누기(최종).xlsx"),
  sheet = "선행연구 데이터5"
)


# ------------------------------------------------------------
# 3. 필수 변수 확인 및 변수명 정리
# ------------------------------------------------------------

required_cols <- c("CropType", "Rainfall", "Study", "SlopeGradient", "TN", "TP")
missing_cols <- setdiff(required_cols, names(df_raw))

if (length(missing_cols) > 0) {
  stop("다음 필수 변수가 데이터에 없습니다: ",
       paste(missing_cols, collapse = ", "))
}

df <- df_raw %>%
  rename(
    crop = CropType,
    rainfall = Rainfall,
    study = Study,
    slope = SlopeGradient
  )


# ------------------------------------------------------------
# 4. 결측값 처리 및 숫자형 변환
# ------------------------------------------------------------

to_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("-", "", "NA", "NaN")] <- NA
  as.numeric(x)
}

to_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("-", "", "NA", "NaN")] <- NA
  x
}

df$TN <- to_num(df$TN)
df$TP <- to_num(df$TP)
df$rainfall <- to_num(df$rainfall)
df$slope <- to_num(df$slope)

df$crop <- to_chr(df$crop)
df$study <- to_chr(df$study)


# ------------------------------------------------------------
# 5. 분석 조건 필터링
# ------------------------------------------------------------

n_before <- nrow(df)

df <- df %>%
  filter(
    !is.na(TN),
    !is.na(TP),
    !is.na(rainfall),
    !is.na(slope),
    !is.na(crop),
    !is.na(study),
    TN > 0,
    TP > 0,
    rainfall <= 200,
    slope < 30
  )

n_after <- nrow(df)

df$crop <- factor(df$crop)
df$study <- factor(df$study)

cat("\n============================================================\n")
cat("Data summary\n")
cat("============================================================\n")
cat("분석 전 n =", n_before, "\n")
cat("분석 후 n =", n_after, "\n")
cat("제외된 행 n =", n_before - n_after, "\n\n")

cat("작물별 표본수:\n")
print(table(df$crop))

cat("\nStudy 수:", nlevels(df$study), "\n")
cat("TN=0:", sum(df$TN == 0), " TP=0:", sum(df$TP == 0), "\n")
cat("Rainfall 범위:", min(df$rainfall), "-", max(df$rainfall), "\n")
cat("Slope 범위:", min(df$slope), "-", max(df$slope), "\n")


# ------------------------------------------------------------
# 6. Random effect 구조
# ------------------------------------------------------------

re_term <- "(1 | study)"

cat("\nRandom effect structure:", re_term, "\n")


# ------------------------------------------------------------
# 7. 모델 적합 함수
# ------------------------------------------------------------

fit_model <- function(yvar, fixed_rhs, disp_rhs = "~ 1") {
  
  f <- as.formula(
    paste0(yvar, " ~ ", fixed_rhs, " + ", re_term)
  )
  
  glmmTMB(
    formula = f,
    dispformula = as.formula(disp_rhs),
    data = df,
    family = Gamma(link = "log"),
    REML = FALSE
  )
}


# ------------------------------------------------------------
# 8. 모델식 설정
# ------------------------------------------------------------

fixed_null <- "ns(rainfall, df = 3) + slope"
fixed_main <- "crop + ns(rainfall, df = 3) + slope"


# ------------------------------------------------------------
# 9. Gamma(log) GLMM 적합
# ------------------------------------------------------------

# 9-1. 기본 dispersion 모델
m_TN_null_base <- fit_model("TN", fixed_null, disp_rhs = "~ 1")
m_TN_main_base <- fit_model("TN", fixed_main, disp_rhs = "~ 1")

m_TP_null_base <- fit_model("TP", fixed_null, disp_rhs = "~ 1")
m_TP_main_base <- fit_model("TP", fixed_main, disp_rhs = "~ 1")


# 9-2. 작물별 dispersion 허용 모델
m_TN_null_disp <- fit_model("TN", fixed_null, disp_rhs = "~ crop")
m_TN_main_disp <- fit_model("TN", fixed_main, disp_rhs = "~ crop")

m_TP_null_disp <- fit_model("TP", fixed_null, disp_rhs = "~ crop")
m_TP_main_disp <- fit_model("TP", fixed_main, disp_rhs = "~ crop")


# ------------------------------------------------------------
# 10. 모델 요약 출력
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Model summaries: base dispersion\n")
cat("============================================================\n")

cat("\n--- TN main-effect model: base dispersion ---\n")
print(summary(m_TN_main_base))

cat("\n--- TP main-effect model: base dispersion ---\n")
print(summary(m_TP_main_base))


cat("\n============================================================\n")
cat("Model summaries: crop-specific dispersion\n")
cat("============================================================\n")

cat("\n--- TN main-effect model: crop-specific dispersion ---\n")
print(summary(m_TN_main_disp))

cat("\n--- TP main-effect model: crop-specific dispersion ---\n")
print(summary(m_TP_main_disp))


# ------------------------------------------------------------
# 11. 작물 주효과 LRT
# ------------------------------------------------------------

lrt_TN_base <- anova(m_TN_null_base, m_TN_main_base)
lrt_TP_base <- anova(m_TP_null_base, m_TP_main_base)

lrt_TN_disp <- anova(m_TN_null_disp, m_TN_main_disp)
lrt_TP_disp <- anova(m_TP_null_disp, m_TP_main_disp)

cat("\n============================================================\n")
cat("Crop main effect LRT\n")
cat("============================================================\n")

cat("\n--- TN crop main effect LRT: base dispersion ---\n")
print(lrt_TN_base)

cat("\n--- TP crop main effect LRT: base dispersion ---\n")
print(lrt_TP_base)

cat("\n--- TN crop main effect LRT: crop-specific dispersion ---\n")
print(lrt_TN_disp)

cat("\n--- TP crop main effect LRT: crop-specific dispersion ---\n")
print(lrt_TP_disp)


# ------------------------------------------------------------
# 12. AIC 비교
# ------------------------------------------------------------

aic_TN <- AIC(
  m_TN_main_base,
  m_TN_main_disp
)

aic_TP <- AIC(
  m_TP_main_base,
  m_TP_main_disp
)

cat("\n============================================================\n")
cat("AIC comparison\n")
cat("============================================================\n")

cat("\n--- TN AIC comparison ---\n")
print(aic_TN)

cat("\n--- TP AIC comparison ---\n")
print(aic_TP)


# ------------------------------------------------------------
# 13. 최종 모델 선택
# ------------------------------------------------------------

use_crop_specific_dispersion <- TRUE

if (use_crop_specific_dispersion) {
  
  m_TN_final <- m_TN_main_disp
  m_TP_final <- m_TP_main_disp
  
  m_TN_null_final <- m_TN_null_disp
  m_TP_null_final <- m_TP_null_disp
  
  lrt_TN_final <- lrt_TN_disp
  lrt_TP_final <- lrt_TP_disp
  
  final_model_type <- "crop-specific dispersion"
  
} else {
  
  m_TN_final <- m_TN_main_base
  m_TP_final <- m_TP_main_base
  
  m_TN_null_final <- m_TN_null_base
  m_TP_null_final <- m_TP_null_base
  
  lrt_TN_final <- lrt_TN_base
  lrt_TP_final <- lrt_TP_base
  
  final_model_type <- "base dispersion"
}

cat("\n============================================================\n")
cat("Final model type:", final_model_type, "\n")
cat("Final model formula:\n")
cat("Response ~ crop + ns(rainfall, df = 3) + slope + (1 | study)\n")
cat("No interaction terms included.\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 14. CLD 생성 함수
# ------------------------------------------------------------

make_cld_from_pairwise <- function(pw_df, emm_df) {
  
  pvals <- pw_df$p.value
  
  pair_names <- as.character(pw_df$contrast)
  pair_names <- gsub(" / ", "-", pair_names)
  pair_names <- gsub(" - ", "-", pair_names)
  names(pvals) <- pair_names
  
  raw_letters <- multcompView::multcompLetters(
    pvals,
    threshold = 0.05
  )$Letters
  
  cld_df <- data.frame(
    crop = names(raw_letters),
    raw_CLD = unname(raw_letters),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  est_col <- if ("response" %in% names(emm_df)) {
    "response"
  } else if ("emmean" %in% names(emm_df)) {
    "emmean"
  } else {
    stop("EMM 결과에서 response 또는 emmean 컬럼을 찾지 못했습니다.")
  }
  
  emm_order <- emm_df %>%
    dplyr::select(crop, estimate = all_of(est_col)) %>%
    dplyr::mutate(crop = as.character(crop)) %>%
    dplyr::arrange(dplyr::desc(estimate))
  
  cld_df <- cld_df %>%
    dplyr::mutate(crop = as.character(crop)) %>%
    dplyr::left_join(emm_order, by = "crop") %>%
    dplyr::arrange(dplyr::desc(estimate))
  
  old_letters <- unique(unlist(strsplit(cld_df$raw_CLD, "")))
  old_letters <- old_letters[old_letters != ""]
  
  letter_rank <- sapply(old_letters, function(ltr) {
    max(cld_df$estimate[grepl(ltr, cld_df$raw_CLD)], na.rm = TRUE)
  })
  
  old_letters_ordered <- old_letters[order(letter_rank, decreasing = TRUE)]
  
  new_letter_pool <- c(letters, LETTERS)
  
  if (length(old_letters_ordered) > length(new_letter_pool)) {
    stop("CLD 문자 수가 너무 많습니다. 문자 풀을 확장해야 합니다.")
  }
  
  new_letters <- new_letter_pool[seq_along(old_letters_ordered)]
  names(new_letters) <- old_letters_ordered
  
  cld_df$CLD <- sapply(cld_df$raw_CLD, function(x) {
    x_split <- unlist(strsplit(x, ""))
    x_split <- x_split[x_split != ""]
    mapped <- new_letters[x_split]
    mapped <- mapped[order(match(mapped, new_letter_pool))]
    paste0(mapped, collapse = "")
  })
  
  cld_df %>%
    dplyr::select(crop, CLD, estimate) %>%
    dplyr::arrange(dplyr::desc(estimate))
}


# ------------------------------------------------------------
# 15. emmeans + Tukey + CLD
# ------------------------------------------------------------

get_emm_results <- function(m) {
  
  emm <- emmeans(
    m,
    ~ crop,
    type = "response"
  )
  
  emm_df <- as.data.frame(emm)
  
  pw <- pairs(
    emm,
    adjust = "tukey"
  )
  
  pw_df <- as.data.frame(summary(pw))
  
  cld_df <- make_cld_from_pairwise(pw_df, emm_df)
  
  list(
    emm = emm_df,
    pairwise = pw_df,
    cld = cld_df
  )
}

res_TN <- get_emm_results(m_TN_final)
res_TP <- get_emm_results(m_TP_final)

cat("\n============================================================\n")
cat("Rainfall- and slope-adjusted EMM\n")
cat("============================================================\n")

cat("\n--- TN adjusted EMM by crop ---\n")
print(res_TN$emm)

cat("\n--- TN CLD: high-to-low letters ---\n")
print(res_TN$cld)

cat("\n--- TP adjusted EMM by crop ---\n")
print(res_TP$emm)

cat("\n--- TP CLD: high-to-low letters ---\n")
print(res_TP$cld)


# ------------------------------------------------------------
# 16. DHARMa 진단 함수
# ------------------------------------------------------------

diagnose <- function(m, label) {
  
  cat("\n================= ", label, " 진단 =================\n")
  
  sim <- simulateResiduals(m, n = 1000)
  plot(sim)
  
  uniform_p <- testUniformity(sim)$p.value
  dispersion_p <- testDispersion(sim)$p.value
  zero_p <- testZeroInflation(sim)$p.value
  singularity <- performance::check_singularity(m)
  r2v <- tryCatch(performance::r2(m), error = function(e) NULL)
  
  cat("[DHARMa 균일성] p =", round(uniform_p, 4),
      " / p > 0.05이면 잔차분포 양호\n")
  cat("[DHARMa 과분산] p =", round(dispersion_p, 4),
      " / p > 0.05이면 과분산 문제 없음\n")
  cat("[영과잉] p =", round(zero_p, 4),
      " / p > 0.05이면 영과잉 문제 없음\n")
  cat("[특이성] :", singularity,
      " / FALSE이면 임의효과 붕괴 문제 없음\n")
  
  if (!is.null(r2v)) {
    print(r2v)
  }
  
  invisible(sim)
}


# ------------------------------------------------------------
# 17. 모델 진단 수행
# ------------------------------------------------------------

sim_TN_base <- diagnose(
  m_TN_main_base,
  "TN main-effect model: base dispersion"
)

sim_TP_base <- diagnose(
  m_TP_main_base,
  "TP main-effect model: base dispersion"
)

sim_TN_disp <- diagnose(
  m_TN_main_disp,
  "TN main-effect model: crop-specific dispersion"
)

sim_TP_disp <- diagnose(
  m_TP_main_disp,
  "TP main-effect model: crop-specific dispersion"
)


# ------------------------------------------------------------
# 18. 공선성 확인
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Collinearity check\n")
cat("============================================================\n")

cat("\n[TN 공선성: final model]\n")
print(performance::check_collinearity(m_TN_final))

cat("\n[TP 공선성: final model]\n")
print(performance::check_collinearity(m_TP_final))

cat("\n[TN 공선성: base model]\n")
print(performance::check_collinearity(m_TN_main_base))

cat("\n[TP 공선성: base model]\n")
print(performance::check_collinearity(m_TP_main_base))


# ------------------------------------------------------------
# 19. 진단 결과 정리 함수
# ------------------------------------------------------------

make_diag_df <- function(m, sim, label) {
  
  r2v <- tryCatch(performance::r2(m), error = function(e) NULL)
  
  data.frame(
    Response = label,
    Uniformity_KS_p = round(testUniformity(sim)$p.value, 4),
    Dispersion_p = round(testDispersion(sim)$p.value, 4),
    ZeroInflation_p = round(testZeroInflation(sim)$p.value, 4),
    Singularity = performance::check_singularity(m),
    R2_marginal = if (!is.null(r2v)) {
      round(as.numeric(r2v$R2_marginal), 3)
    } else {
      NA
    },
    R2_conditional = if (!is.null(r2v)) {
      round(as.numeric(r2v$R2_conditional), 3)
    } else {
      NA
    }
  )
}

diag_df <- rbind(
  make_diag_df(m_TN_main_base, sim_TN_base, "TN_base"),
  make_diag_df(m_TP_main_base, sim_TP_base, "TP_base"),
  make_diag_df(m_TN_main_disp, sim_TN_disp, "TN_disp_crop"),
  make_diag_df(m_TP_main_disp, sim_TP_disp, "TP_disp_crop")
)


# ------------------------------------------------------------
# 20. LRT 결과 정리
# ------------------------------------------------------------

extract_lrt <- function(lrt_obj, response, model_type) {
  
  data.frame(
    Response = response,
    Model_type = model_type,
    Chisq = lrt_obj$Chisq[2],
    df = lrt_obj$`Chi Df`[2],
    p_value = lrt_obj$`Pr(>Chisq)`[2]
  )
}

lrt_all <- rbind(
  extract_lrt(lrt_TN_base, "TN", "base dispersion"),
  extract_lrt(lrt_TP_base, "TP", "base dispersion"),
  extract_lrt(lrt_TN_disp, "TN", "crop-specific dispersion"),
  extract_lrt(lrt_TP_disp, "TP", "crop-specific dispersion")
)

lrt_final <- rbind(
  extract_lrt(lrt_TN_final, "TN", final_model_type),
  extract_lrt(lrt_TP_final, "TP", final_model_type)
)


# ------------------------------------------------------------
# 21. 모델 정보 정리
# ------------------------------------------------------------

model_info <- data.frame(
  Item = c(
    "Analysis objective",
    "Final conditional model",
    "Interaction terms",
    "Rainfall treatment",
    "Slope treatment",
    "Random effect",
    "Distribution",
    "Link function",
    "Crop effect test",
    "Post-hoc test",
    "Interpretation allowed",
    "Interpretation not allowed"
  ),
  Description = c(
    "Compare TN/TP loads among crop types after adjusting for rainfall and slope",
    "Response ~ crop + ns(rainfall, df = 3) + slope + (1 | study)",
    "Not included",
    "Natural spline covariate, df = 3",
    "Linear covariate",
    "Study-level random intercept",
    "Gamma",
    "log",
    "Likelihood ratio test comparing models with and without crop",
    "Tukey-adjusted pairwise comparisons using emmeans",
    "Crop-type differences after adjustment for rainfall and slope",
    "Crop-by-slope or crop-by-rainfall interaction effects"
  )
)

final_model_choice <- data.frame(
  Response = c("TN", "TP"),
  Final_model_type = final_model_type,
  Conditional_formula = "Response ~ crop + ns(rainfall, df = 3) + slope + (1 | study)",
  Dispersion_formula = ifelse(use_crop_specific_dispersion, "~ crop", "~ 1")
)


# ------------------------------------------------------------
# 22. 엑셀 저장
# ------------------------------------------------------------

wb <- createWorkbook()

addWorksheet(wb, "Model_info")
writeData(wb, "Model_info", model_info)

addWorksheet(wb, "Final_model_choice")
writeData(wb, "Final_model_choice", final_model_choice)

addWorksheet(wb, "Data_summary")
writeData(
  wb,
  "Data_summary",
  data.frame(
    n_before = n_before,
    n_after = n_after,
    n_removed = n_before - n_after,
    n_study = nlevels(df$study),
    rainfall_min = min(df$rainfall),
    rainfall_max = max(df$rainfall),
    slope_min = min(df$slope),
    slope_max = max(df$slope)
  )
)

addWorksheet(wb, "Crop_n")
writeData(
  wb,
  "Crop_n",
  as.data.frame(table(df$crop)) %>%
    rename(crop = Var1, n = Freq)
)

addWorksheet(wb, "LRT_all")
writeData(wb, "LRT_all", lrt_all)

addWorksheet(wb, "LRT_final")
writeData(wb, "LRT_final", lrt_final)

addWorksheet(wb, "AIC_TN")
writeData(wb, "AIC_TN", as.data.frame(aic_TN))

addWorksheet(wb, "AIC_TP")
writeData(wb, "AIC_TP", as.data.frame(aic_TP))

addWorksheet(wb, "EMM_TN")
writeData(wb, "EMM_TN", res_TN$emm)

addWorksheet(wb, "EMM_TP")
writeData(wb, "EMM_TP", res_TP$emm)

addWorksheet(wb, "CLD_TN")
writeData(wb, "CLD_TN", res_TN$cld)

addWorksheet(wb, "CLD_TP")
writeData(wb, "CLD_TP", res_TP$cld)

addWorksheet(wb, "Pairwise_TN")
writeData(wb, "Pairwise_TN", res_TN$pairwise)

addWorksheet(wb, "Pairwise_TP")
writeData(wb, "Pairwise_TP", res_TP$pairwise)

addWorksheet(wb, "Diagnostics")
writeData(wb, "Diagnostics", diag_df)

saveWorkbook(
  wb,
  file.path(out_dir, "GLMM_crop_main_effect_rainfall_slope_adjusted_results.xlsx"),
  overwrite = TRUE
)

cat("\n============================================================\n")
cat("엑셀 저장 완료:\n")
cat(file.path(out_dir, "GLMM_crop_main_effect_rainfall_slope_adjusted_results.xlsx"), "\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 23. 그림 함수
# ------------------------------------------------------------
# 그림 의미:
#  - 막대: 강우량과 경사도를 보정한 작물별 추정평균
#  - 에러바: 95% CI
#  - 문자: Tukey CLD
#  - 오른쪽 위 p-value: 전체 crop main effect LRT
#
# 주의:
#  - 원자료 점은 표시하지 않음
#  - 그림 제목은 표시하지 않음
#  - 이 그림은 상호작용 그림이 아님
#  - 경사도별 반응 곡선을 보여주는 그림도 아님
#  - 보정평균 차이 그림임
#  - panel.border로 상단/우측 포함 네모 테두리 표시
# ------------------------------------------------------------

# 작물유형별 막대 색상
crop_colors <- c(
  "Cc" = "#FFD700",  # Cereal crops
  "Fc" = "#FF6347",  # Fruit crops
  "Lv" = "#32CD32",  # Leafy vegetables
  "Lg" = "#6A5ACD",  # Legumes
  "Oc" = "#FFA500",  # Oil crops
  "Pc" = "#20B2AA",  # Perennial crops
  "Rv" = "#8B4513"   # Root vegetables
)

# crop 값이 풀네임 또는 약어로 들어온 경우 모두 코드명으로 변환
crop_name_to_code <- c(
  "Cc" = "Cc",
  "Fc" = "Fc",
  "Lv" = "Lv",
  "Lg" = "Lg",
  "Oc" = "Oc",
  "Pc" = "Pc",
  "Rv" = "Rv",
  "Cereal crops" = "Cc",
  "Fruit crops" = "Fc",
  "Leafy vegetables" = "Lv",
  "Legumes" = "Lg",
  "Oil crops" = "Oc",
  "Perennial crops" = "Pc",
  "Root vegetables" = "Rv",
  "Cereal crop" = "Cc",
  "Fruit crop" = "Fc",
  "Leafy vegetable" = "Lv",
  "Oil crop" = "Oc",
  "Perennial crop" = "Pc",
  "Root vegetable" = "Rv"
)

standardize_crop_code <- function(x) {
  
  x <- trimws(as.character(x))
  
  out <- ifelse(
    x %in% names(crop_name_to_code),
    crop_name_to_code[x],
    x
  )
  
  out
}

# p-value 표시 형식
format_p <- function(p) {
  
  if (length(p) == 0 || is.na(p)) {
    return("p = NA")
  }
  
  if (p < 0.001) {
    return("p < 0.001")
  }
  
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

plot_bar <- function(res, yvar_raw, ylab, global_p = NULL,
                     y_digits = 1, y_max = NULL, y_break_by = NULL) {
  
  emm <- res$emm
  
  est_col <- if ("response" %in% names(emm)) {
    "response"
  } else {
    "emmean"
  }
  
  lo_col <- if ("asymp.LCL" %in% names(emm)) {
    "asymp.LCL"
  } else {
    "lower.CL"
  }
  
  hi_col <- if ("asymp.UCL" %in% names(emm)) {
    "asymp.UCL"
  } else {
    "upper.CL"
  }
  
  emm$est <- emm[[est_col]]
  emm$lo <- emm[[lo_col]]
  emm$hi <- emm[[hi_col]]
  
  emm <- emm %>%
    dplyr::mutate(
      crop_original = trimws(as.character(crop))
    ) %>%
    dplyr::left_join(
      res$cld %>%
        dplyr::mutate(
          crop_original = trimws(as.character(crop))
        ) %>%
        dplyr::select(crop_original, CLD),
      by = "crop_original"
    ) %>%
    dplyr::arrange(dplyr::desc(est)) %>%
    dplyr::mutate(
      crop_code = standardize_crop_code(crop_original)
    )
  
  emm$crop_code <- factor(emm$crop_code, levels = emm$crop_code)
  
  # y_max가 지정되면 고정 축 범위 사용, 아니면 자동 범위 사용
  ytop <- if (!is.null(y_max)) {
    y_max
  } else {
    max(emm$hi, na.rm = TRUE) * 1.25
  }
  
  y_breaks <- if (!is.null(y_break_by)) {
    seq(0, ytop, by = y_break_by)
  } else {
    waiver()
  }
  
  ggplot() +
    geom_col(
      data = emm,
      aes(x = crop_code, y = est, fill = crop_code),
      width = 0.65,
      alpha = 0.85,
      color = "black",
      linewidth = 0.4
    ) +
    geom_errorbar(
      data = emm,
      aes(x = crop_code, ymin = lo, ymax = hi),
      width = 0.18,
      linewidth = 0.7,
      color = "grey20"
    ) +
    geom_text(
      data = emm,
      aes(x = crop_code, y = hi, label = CLD),
      vjust = -0.6,
      fontface = "bold",
      size = 6,
      color = "black"
    ) +
    annotate(
      "text",
      x = Inf,
      y = ytop * 0.98,
      label = paste0("Crop type effect (LRT): ", format_p(global_p)),
      hjust = 1.05,
      vjust = 1,
      size = 5.0,
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(
      values = crop_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) formatC(
        x,
        format = "f",
        digits = y_digits
      )
    ) +
    coord_cartesian(ylim = c(0, ytop)) +
    labs(
      x = NULL,
      y = ylab,
      title = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      
      # 그래프 네모 테두리: 상단/우측 포함
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      # 왼쪽/아래 축선과 panel.border 중복 방지
      axis.line = element_blank(),
      
      # X축 작물 코드 글씨 크기
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        size = 16,
        color = "black"
      ),
      
      # Y축 눈금 숫자 글씨 크기
      axis.text.y = element_text(
        size = 16,
        color = "black"
      ),
      
      # Y축 라벨 글씨 크기
      axis.title.y = element_text(
        size = 18,
        color = "black"
      ),
      
      # X축 라벨은 x = NULL이라 표시되지 않음
      axis.title.x = element_blank(),
      
      # 그림 제목 제거
      plot.title = element_blank()
    )
}


# ------------------------------------------------------------
# 24. TN/TP 그림 생성
# ------------------------------------------------------------

p_TN <- plot_bar(
  res_TN,
  "TN",
  "TN Load (kg/ha/event)",
  global_p = lrt_final$p_value[lrt_final$Response == "TN"],
  y_digits = 1,
  y_max = 8.0,
  y_break_by = 2.0
)

p_TP <- plot_bar(
  res_TP,
  "TP",
  "TP Load (kg/ha/event)",
  global_p = lrt_final$p_value[lrt_final$Response == "TP"],
  y_digits = 1
)

print(p_TN)
print(p_TP)


# ------------------------------------------------------------
# 25. 그림 저장
# ------------------------------------------------------------

ggsave(
  file.path(out_dir, "Fig_crop_Gamma_TN_main_effect_rainfall_slope_adjusted.png"),
  p_TN,
  width = 5,
  height = 5,
  dpi = 500
)

ggsave(
  file.path(out_dir, "Fig_crop_Gamma_TP_main_effect_rainfall_slope_adjusted.png"),
  p_TP,
  width = 5,
  height = 5,
  dpi = 500
)

cat("\n============================================================\n")
cat("그림 저장 완료. 저장 폴더:\n")
cat(out_dir, "\n")
cat("============================================================\n")


# ============================================================
# 해석 메모
# ============================================================
#
# 1. 이 분석의 핵심 질문
#    강우량과 경사도를 보정한 후에도 작물유형별 TN/TP 부하량 차이가 있는가?
#
# 2. 최종 모델
#    TN ~ crop + ns(rainfall, df = 3) + slope + (1 | study)
#    TP ~ crop + ns(rainfall, df = 3) + slope + (1 | study)
#
# 3. 상호작용항
#    포함하지 않음.
#    따라서 crop:slope, crop:rainfall 효과는 검정하지 않음.
#
# 4. LRT_final
#    crop이 없는 모델과 crop이 있는 모델을 비교.
#    p < 0.05이면:
#    강우량, 경사도, study-level dependence를 고려한 후에도
#    작물유형이 TN/TP 부하량에 유의한 영향을 미친다고 해석.
#
# 5. EMM_TN / EMM_TP
#    강우량과 경사도를 보정한 작물별 추정평균.
#    원자료 평균이 아니라 모형 기반 보정평균.
#
# 6. Pairwise_TN / Pairwise_TP
#    Tukey 보정 사후검정.
#    p < 0.05이면 해당 작물 쌍 간 보정평균 차이가 유의함.
#
# 7. CLD_TN / CLD_TP
#    추정평균이 높은 작물부터 a, b, c...로 배정.
#    같은 문자를 공유하면 유의한 차이가 없고,
#    공유 문자가 없으면 Tukey 기준 유의한 차이가 있음.
#
# 8. DHARMa diagnostics
#    효과 검정과 달리 진단검정은 p > 0.05가 양호.
#    p < 0.05이면 잔차분포, 과분산 등 문제가 남아 있을 수 있음.
#
# 9. 논문식 해석 예시
#    A Gamma GLMM with a log link was used to compare TN and TP loads
#    among crop types after adjusting for rainfall and slope, with study
#    included as a random intercept. The significance of crop type was
#    assessed using likelihood ratio tests comparing models with and
#    without crop type. Adjusted marginal means were compared using
#    Tukey-adjusted pairwise tests.
#
# 10. 한국어 해석 예시
#    강우량과 경사도를 보정하고 연구 단위의 비독립성을 고려한
#    Gamma(log) GLMM 분석 결과, 작물유형이 TN/TP 부하량에 미치는
#    주효과를 검정하였다. 작물유형의 유의성은 작물유형을 포함한
#    모형과 제외한 모형 간 우도비 검정을 통해 평가하였으며,
#    작물유형 간 보정평균 차이는 Tukey 사후검정으로 비교하였다.
#
# ============================================================