# ---- 0. Packages ----
# install.packages(c("readxl","mgcv","dplyr","writexl","rstudioapi"))  # first run only
library(readxl)
library(mgcv)
library(dplyr)
library(writexl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

# ---- 1. Read data ----
# Place the data file in the same folder as this script, or edit the path below.
# 'sheet' must match the worksheet name inside the Excel file.
path <- "metadata.xlsx"

raw <- read_excel(path, sheet = "rawdata")

if (!"Sites" %in% names(raw)) {
  stop("The source data has no 'Sites' column. Check the actual site-name column and edit rename(site = Sites).")
}

# ---- 2. Preprocessing ----
keep_crops <- c("Leafy vegetables", "Root vegetables")  # the first level (leafy) is the reference

dat <- raw %>%
  mutate(across(c(SlopeGradient, Rainfall, TN, TP), as.numeric)) %>%
  filter(!is.na(SlopeGradient), !is.na(Rainfall), !is.na(TN), !is.na(TP)) %>%
  filter(SlopeGradient <= 30, Rainfall <= 200) %>%
  filter(CropType %in% keep_crops) %>%
  rename(
    slope = SlopeGradient,
    rainfall = Rainfall,
    crop = CropType,
    site = Sites
  ) %>%
  mutate(
    crop     = factor(crop, levels = keep_crops),
    crop_ord = as.ordered(crop),     # ordered factor for the difference curve (Leafy < Root)
    site     = factor(site),
    TN_s     = sqrt(pmax(TN, 0)),
    TP_s     = sqrt(pmax(TP, 0))
  )

cat("Analysis n =", nrow(dat),
    "| crops =", nlevels(dat$crop),
    "| sites =", nlevels(dat$site),
    "| reference crop =", levels(dat$crop_ord)[1], "\n")
print(table(dat$crop))

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# ---- 3. Fit the GI model (REML) ----
fit_GI <- function(resp) {
  f <- as.formula(
    paste0(
      resp, " ~ crop_ord + ",
      "s(slope, k = 4) + ",                 # reference (leafy) slope curve
      "s(slope, by = crop_ord, k = 4) + ",  # root difference curve = interaction
      "s(rainfall, k = 5) + ",
      "s(site, bs = 're')"
    )
  )
  gam(f, data = dat, method = "REML")
}

m_TN <- fit_GI("TN_s")
m_TP <- fit_GI("TP_s")

cat("\n===== TN summary =====\n"); print(summary(m_TN))
cat("\n===== TP summary =====\n"); print(summary(m_TP))

# ---- 4. Smooth-term table + interaction (difference-curve) p-value extraction ----
# Row names in s.table:
#   "s(slope)"                          -> reference (leafy) slope curve
#   "s(slope):crop_ordRoot vegetables"  -> root difference curve = crop x slope interaction
ref_label  <- "s(slope)"
get_diff_row <- function(model) {
  rn <- rownames(summary(model)$s.table)
  grep("slope.*crop_ord", rn, value = TRUE)[1]   # difference-curve row
}

smooth_tab <- function(model, resp) {
  s <- summary(model)$s.table
  data.frame(
    Response = resp,
    Term = rownames(s),
    edf  = round(s[, "edf"], 2),
    F    = round(s[, "F"], 2),
    p    = signif(s[, "p-value"], 3),
    p_label = fmt_p(s[, "p-value"]),
    row.names = NULL
  )
}

tab_smooth <- rbind(
  smooth_tab(m_TN, "sqrt(TN)"),
  smooth_tab(m_TP, "sqrt(TP)")
)
cat("\n===== Smooth-term significance (the difference-curve row is the interaction test) =====\n")
print(tab_smooth)

interaction_p <- function(model, resp) {
  s  <- summary(model)$s.table
  dr <- get_diff_row(model)
  data.frame(
    Response       = resp,
    Diff_term      = dr,
    Interaction_edf = round(s[dr, "edf"], 2),
    Interaction_F   = round(s[dr, "F"], 2),
    Interaction_p   = signif(s[dr, "p-value"], 3),
    Interaction_p_label = fmt_p(s[dr, "p-value"])
  )
}

inter_TN <- interaction_p(m_TN, "sqrt(TN)")
inter_TP <- interaction_p(m_TP, "sqrt(TP)")
interaction_summary <- rbind(inter_TN, inter_TP)

cat("\n===== crop x slope interaction (single difference-curve test) =====\n")
print(interaction_summary)

# ---- 5. Contribution: deviance-explained increment when the difference curve is removed (gam.hp substitute) ----
fit_reduced <- function(resp) {
  f <- as.formula(
    paste0(resp, " ~ crop_ord + s(slope, k = 4) + ",
           "s(rainfall, k = 5) + s(site, bs = 're')")
  )
  gam(f, data = dat, method = "REML")
}

dev_increment <- function(full, resp) {
  red <- fit_reduced(resp)
  data.frame(
    Response = resp,
    DevExpl_no_interaction = round(summary(red)$dev.expl  * 100, 2),
    DevExpl_with_interaction = round(summary(full)$dev.expl * 100, 2),
    Delta_DevExpl_percent = round(
      (summary(full)$dev.expl - summary(red)$dev.expl) * 100, 2
    )
  )
}

contrib_summary <- rbind(
  dev_increment(m_TN, "sqrt(TN)"),
  dev_increment(m_TP, "sqrt(TP)")
)
cat("\n===== Deviance-explained increment (relative contribution of the interaction, effect size) =====\n")
print(contrib_summary)

# Reporting labels
make_label <- function(ip) {
  paste0("Crop \u00d7 slope: p = ", ip$Interaction_p_label,
         ", edf = ", sprintf("%.1f", ip$Interaction_edf))
}
lab_TN <- make_label(inter_TN)
lab_TP <- make_label(inter_TP)
cat("\nTN:", lab_TN, "\nTP:", lab_TP, "\n")

# ---- 6. Figure: two-panel crop-specific slope curves ----
# predict() automatically composes (reference curve + difference curve),
# so only the crop level needs to be specified.
cols  <- c("Leafy vegetables" = "#2E7D32", "Root vegetables" = "#C62828")
short <- c("Leafy vegetables" = "Leafy vegetables",   "Root vegetables" = "Root vegetables")

curve_panel <- function(model, resp_raw, ylab, title, inter_label) {
  yo <- dat[[resp_raw]]
  plot(NA,
       xlim = range(dat$slope),
       ylim = c(0, quantile(yo, 0.99, na.rm = TRUE) * 1.10),
       xlab = "Slope gradient (%)", ylab = ylab, las = 1)
  mtext(title, adj = 0, line = 1.2, font = 2, cex = 1.0)
  
  usr <- par("usr")
  text(usr[1] + 0.03 * (usr[2] - usr[1]),
       usr[4] - 0.06 * (usr[4] - usr[3]),
       labels = inter_label, adj = c(0, 1), cex = 0.80)
  
  # 95% CI
  for (cg in levels(dat$crop_ord)) {
    s <- dat[dat$crop_ord == cg, ]; if (nrow(s) < 2) next
    xs <- seq(min(s$slope), max(s$slope), length = 100)
    nd <- data.frame(slope = xs, rainfall = mean(dat$rainfall),
                     crop_ord = ordered(cg, levels = levels(dat$crop_ord)),
                     site = factor(levels(dat$site)[1], levels = levels(dat$site)))
    pr <- predict(model, nd, se.fit = TRUE, exclude = "s(site)")
    upper <- (pr$fit + 1.96 * pr$se.fit)^2
    lower <- pmax(pr$fit - 1.96 * pr$se.fit, 0)^2
    polygon(c(xs, rev(xs)), c(upper, rev(lower)),
            col = adjustcolor(cols[cg], 0.12),
            border = adjustcolor(cols[cg], 0.6), lty = 3, lwd = 0.5)
  }
  # Observed values
  for (cg in levels(dat$crop_ord)) {
    s <- dat[dat$crop_ord == cg, ]
    points(s$slope, s[[resp_raw]], col = adjustcolor(cols[cg], 0.35),
           pch = 16, cex = 0.55)
  }
  # Predicted lines
  for (cg in levels(dat$crop_ord)) {
    s <- dat[dat$crop_ord == cg, ]; if (nrow(s) < 2) next
    xs <- seq(min(s$slope), max(s$slope), length = 100)
    nd <- data.frame(slope = xs, rainfall = mean(dat$rainfall),
                     crop_ord = ordered(cg, levels = levels(dat$crop_ord)),
                     site = factor(levels(dat$site)[1], levels = levels(dat$site)))
    pr <- predict(model, nd, exclude = "s(site)")
    lines(xs, pr^2, col = cols[cg], lwd = 2.6)
  }
  legend("topright",
         legend = c(short[levels(dat$crop_ord)], "95% CI"),
         col = c(cols[levels(dat$crop_ord)], "grey50"),
         lwd = c(rep(2.4, nlevels(dat$crop_ord)), 0.5),
         lty = c(rep(1, nlevels(dat$crop_ord)), 3),
         bty = "n", cex = 0.72)
}

draw_all <- function() {
  par(mfrow = c(1, 2), mar = c(4.4, 4.8, 2.8, 1.2), mgp = c(2.8, 0.8, 0))
  curve_panel(m_TN, "TN",
              expression("TN load (kg/ha/event)"),
              "(a) TN - crop-specific slope response", lab_TN)
  curve_panel(m_TP, "TP",
              expression("TP load (kg/ha/event)"),
              "(b) TP - crop-specific slope response", lab_TP)
}

draw_all()
png("Fig5_GAMM_diffsmooth_2crop.png", width = 2600, height = 1200, res = 215)
draw_all()
dev.off()
cat("\nFigure saved:", file.path(getwd(), "Fig5_GAMM_diffsmooth_2crop.png"), "\n")

# ---- 7. Check k adequacy (recommended) ----
cat("\n===== gam.check (k adequacy: TN) =====\n"); gam.check(m_TN)
cat("\n===== gam.check (k adequacy: TP) =====\n"); gam.check(m_TP)

# ---- 8. Export results to Excel ----
out <- list(
  data_summary = data.frame(
    n = nrow(dat), n_crop = nlevels(dat$crop), n_site = nlevels(dat$site),
    reference_crop = levels(dat$crop_ord)[1],
    slope_min = min(dat$slope), slope_max = max(dat$slope),
    rainfall_min = min(dat$rainfall), rainfall_max = max(dat$rainfall)
  ),
  crop_count          = as.data.frame(table(dat$crop)) %>% rename(Crop = Var1, n = Freq),
  smooth_significance = tab_smooth,
  interaction_test    = interaction_summary,   # single difference-curve p-value = interaction
  contribution_devexpl = contrib_summary        # gam.hp substitute: deviance-explained increment
)
write_xlsx(out, "GAMM_results_diffsmooth_2crop.xlsx")
cat("Tables saved:", file.path(getwd(), "GAMM_results_diffsmooth_2crop.xlsx"), "\n")