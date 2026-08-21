# =============================================================================
# SPATIAL CAUSAL MACHINE LEARNING: THEORY AND AN APPLICATION USING CAUSAL FORESTS
# Outcome:   Building permit density (BP_POST_DEN, cumulative 2015–2019)
# Treatment: 606/Bloomingdale Trail access points — 3 specifications
# Estimator: Causal Forest (grf), standard and spatially cross-validated
# =============================================================================

# =============================================================================
# SECTION 0: SETUP
# =============================================================================

rm(list = ls())
MASTER_SEED <- 12051984

library(sf)
library(spdep)
library(grf)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(leaflet)
library(RColorBrewer)
library(conflicted)

has_leafsync <- requireNamespace("leafsync", quietly = TRUE)
if (has_leafsync) library(leafsync)

conflict_prefer("filter",    "dplyr")
conflict_prefer("select",    "dplyr")
conflict_prefer("lag",       "dplyr")
conflict_prefer("mutate",    "dplyr")
conflict_prefer("summarise", "dplyr")

PROJECT_DIR <- "/Users/kevincredit/Library/CloudStorage/Dropbox/Packages/SArf paper/Book chapter/Spatial_Causal_ML/Spatial Causal ML"
DATA_DIR    <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR  <- file.path(PROJECT_DIR, "outputs")
setwd(PROJECT_DIR)

cat("=== Setup complete. MASTER_SEED =", MASTER_SEED, "===\n\n")


# =============================================================================
# SECTION 1: DATA LOADING & PREPARATION
# =============================================================================

chicago_sf <- st_read(file.path(DATA_DIR, "Chicago_BGs_Covariates.shp"), quiet = TRUE)

chicago_sf <- chicago_sf %>%
  filter(!is.na(MEDAGE10), !is.na(PCIN10), !is.na(MHHIN10),
         !is.na(MYRBLT10), !is.na(MYRMOV10))

# Neighbourhood control area
C_CAs <- chicago_sf %>% filter(TREAT == 1 | CONTR_1 == 1)
C_CA  <- as.data.frame(C_CAs)

cat("N =", nrow(C_CA), "| Treated:", sum(C_CA$TREAT),
    "| Control:", sum(1 - C_CA$TREAT), "\n\n")


# =============================================================================
# SECTION 2: VARIABLE CONSTRUCTION
# =============================================================================

# --- Spatial weights & CV folds (built first — needed for lag covariate) ------
coords <- st_coordinates(st_centroid(C_CAs))
knn5   <- knn2nb(knearneigh(coords, k = 5))
lw     <- nb2listw(knn5, style = "W")

K_FOLDS <- 5
set.seed(MASTER_SEED)
km_fit   <- kmeans(coords, centers = K_FOLDS, nstart = 25)
sp_folds <- km_fit$cluster

cat("Spatial CV fold sizes:", table(sp_folds), "\n")

# Map of spatial CV folds
p_folds <- ggplot() +
  geom_sf(data = C_CAs %>% mutate(fold = factor(sp_folds)),
          aes(fill = fold), colour = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = "Fold") +
  labs(title = "Spatial Cross-Validation Folds (K = 5)",
       subtitle = "K-means geographic blocks") +
  theme_void(base_size = 11) +
  theme(legend.position = "right")
print(p_folds)
ggsave(file.path(OUTPUT_DIR, "p_folds.png"), plot = p_folds,
       width = 7, height = 6, dpi = 300, bg = "white")

# --- Outcome ------------------------------------------------------------------
# Building permit densities (per km²)
# Transform to UTM zone 16N (EPSG:32616, metres) for correct metric area
bp_area_km2  <- as.numeric(st_area(st_transform(C_CAs, 32616))) / 1e6

# Pre-intervention: cumulative permits 2010–2014 (before 606 construction)
# Post-intervention: cumulative permits 2015–2019 (construction + early effects)
C_CA$BP_PRE_DEN  <- rowSums(C_CA[, paste0("BP", 2010:2014)]) / bp_area_km2
C_CA$BP_POST_DEN <- rowSums(C_CA[, paste0("BP", 2015:2019)]) / bp_area_km2

Y_bp_raw     <- C_CA$BP_POST_DEN         # density, for mapping
Y_bp         <- scale(Y_bp_raw)[, 1]     # standardised, for model
Y_bp_density <- Y_bp_raw                 # alias for map section

# --- Covariates ---------------------------------------------------------------
# All spatial lags computed in-script using knn5 weights matrix (k = 5 nearest
# neighbours, row-standardised). Overrides any pre-computed lag columns from the
# shapefile so that lag specification is consistent and transparent.
# NOTE: spatial lag of post-intervention outcome intentionally excluded — bad control.

cov_bp_base <- c("BP_PRE_DEN","BIZ_ZONEP","MEDAGE10","BLKP10","HSPP10",
                 "BACHP10","UNEMP10","MBSAP10","MHHIN10","OWNP10","MYRMOV10")

cov_bp_lag  <- c("BP_PRE_DEN_LAG","BIZZ_LAG","MEDA_LAG","BLKP_LAG","HSPP_LAG",
                 "BACH_LAG","UNEM_LAG","MBSA_LAG","MHHI_LAG","OWN_LAG","YRMV_LAG")

for (i in seq_along(cov_bp_base)) {
  C_CA[[cov_bp_lag[i]]] <- lag.listw(lw, C_CA[[cov_bp_base[i]]])
}

X_bp_base <- as.matrix(scale(C_CA[, cov_bp_base]))
X_bp_slx  <- as.matrix(scale(C_CA[, c(cov_bp_base, cov_bp_lag)]))
colnames(X_bp_base) <- cov_bp_base
colnames(X_bp_slx)  <- c(cov_bp_base, cov_bp_lag)

cat("BP base covariates:", length(cov_bp_base),
    "| SLX covariates:", length(c(cov_bp_base, cov_bp_lag)), "\n")

# --- Treatment specifications -------------------------------------------------
# (1) Binary T — do NOT scale for grf
T_bin <- C_CA$TREAT

# (2) T + WT combined — scaled for causal forest
T_twt <- scale(C_CA$TREAT + C_CA$TREAT_LAG)[, 1]

# (3) Logistic distance decay — scaled for causal forest
T_dec <- scale(
  1 - (1 / (exp((800/180) - (0.48/60) * ((C_CA$distance * 3600) / 5000)) + 1))
)[, 1]

cat("\n")


# --- Treatment specification shape plot (actual data) -------------------------------------------------
T_bin_raw <- C_CA$TREAT
T_twt_raw <- (C_CA$TREAT + C_CA$TREAT_LAG) /
              max(C_CA$TREAT + C_CA$TREAT_LAG, na.rm = TRUE)
T_dec_raw <- 1 - (1 / (exp((800/180) -
              (0.48/60) * ((C_CA$distance * 3600) / 5000)) + 1))

# Smooth logistic line over actual distance range
dist_seq   <- seq(min(C_CA$distance), max(C_CA$distance), length.out = 500)
logit_line <- data.frame(
  distance = dist_seq,
  value    = 1 - (1 / (exp((800/180) -
               (0.48/60) * ((dist_seq * 3600) / 5000)) + 1)),
  spec     = "(3) Logistic decay"
)

shape_df <- data.frame(
  distance = rep(C_CA$distance, 3),
  value    = c(T_bin_raw, T_twt_raw, T_dec_raw),
  spec     = factor(
    rep(c("(1) Binary", "(2) T + WT (norm.)", "(3) Logistic decay"),
        each = nrow(C_CA)),
    levels = c("(1) Binary", "(2) T + WT (norm.)", "(3) Logistic decay"))
)

spec_cols <- c("(1) Binary"         = "#08519C",
               "(2) T + WT (norm.)" = "#6BAED6",
               "(3) Logistic decay" = "#238B45")

p_treat_shape <- ggplot(shape_df,
                        aes(x = distance, y = value, colour = spec)) +
  geom_jitter(height = 0.015, width = 0, alpha = 0.5, size = 1.2) +
  # Binary: sharp step line
  geom_step(data = shape_df %>% dplyr::filter(spec == "(1) Binary") %>%
                   dplyr::arrange(distance),
            aes(x = distance, y = value),
            colour = "#08519C", linewidth = 1.0, inherit.aes = FALSE) +
  # T + WT: loess smooth
  geom_smooth(data = shape_df %>% dplyr::filter(spec == "(2) T + WT (norm.)"),
              aes(x = distance, y = value),
              method = "loess", se = FALSE,
              colour = "#6BAED6", linewidth = 1.0, inherit.aes = FALSE) +
  # Logistic decay: exact formula curve
  geom_line(data = logit_line, aes(x = distance, y = value),
            colour = "#238B45", linewidth = 1.1, inherit.aes = FALSE) +
  scale_colour_manual(values = spec_cols, name = NULL) +
  scale_y_continuous(limits = c(-0.05, 1.1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Distance from access points (m)",
       y = "Treatment intensity (0–1)",
       title = "Treatment specification shapes") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "bottom",
        plot.background  = element_rect(fill = "white", colour = NA))

print(p_treat_shape)
ggsave(file.path(OUTPUT_DIR, "p_treat_shape.png"), plot = p_treat_shape,
       width = 6, height = 4, dpi = 300, bg = "white")


# =============================================================================
# SECTION 3: SPATIAL AUTOCORRELATION DIAGNOSTIC (MORAN'S I ON RAW Y)
# =============================================================================
# Pre-modelling question: does BP_POST_DEN exhibit spatial autocorrelation?
# Significant → SLX specification.

mi_bp <- moran.test(C_CA$BP_POST_DEN, lw)

cat(sprintf("Moran's I on BP_POST_DEN (2015-2019): statistic = %.4f,  p = %.4f\n",
            mi_bp$statistic, mi_bp$p.value))

if (mi_bp$p.value < 0.05) {
  cat("=> Significant spatial autocorrelation. SLX specification suggested. \n\n")
} else {
  cat("=> No significant spatial autocorrelation. SLX recommended as robustness.\n\n")
}


# =============================================================================
# SECTION 4: STANDARD CAUSAL FOREST — BUILDING PERMITS (3 SPECIFICATIONS)
# =============================================================================
# Doubly-robust causal forest via grf (Athey et al. 2019).
# Nuisance functions Y.hat = E[Y|X] and W.hat = E[W|X] are estimated via
# regression forests using grf's internal OOB cross-fitting, then passed
# to causal_forest() which estimates tau(x) via the R-learner objective:
#   (Y - Y.hat) = tau(X) * (W - W.hat) + epsilon

cat("--- Section 4: Standard causal forest (building permits) ---\n")

fit_causal_forest <- function(X, Y, W, label, seed) {
  set.seed(seed)
  nuisance_Y <- regression_forest(X, Y, seed = seed)
  set.seed(seed)
  nuisance_W <- regression_forest(X, W, seed = seed)

  Y_hat <- predict(nuisance_Y)$predictions
  W_hat <- predict(nuisance_W)$predictions

  set.seed(seed)
  cf  <- causal_forest(X, Y, W, Y.hat = Y_hat, W.hat = W_hat, seed = seed)
  ate <- average_treatment_effect(cf, target.sample = "all")
  tau <- predict(cf)$predictions

  cat(sprintf("  [%s]: ATE = %.4f  SE = %.4f\n", label, ate[1], ate[2]))
  list(model = cf, ate = ate, tau = tau, Y_hat = Y_hat, W_hat = W_hat, label = label)
}

bp_cf1 <- fit_causal_forest(X_bp_slx, Y_bp, T_bin, "Binary T",       seed = MASTER_SEED)
bp_cf2 <- fit_causal_forest(X_bp_slx,  Y_bp, T_twt, "T+WT combined",  seed = MASTER_SEED)
bp_cf3 <- fit_causal_forest(X_bp_slx,  Y_bp, T_dec, "Distance decay", seed = MASTER_SEED)

cat("\n")


# =============================================================================
# SECTION 5: SPATIALLY CROSS-VALIDATED CAUSAL FOREST — BUILDING PERMITS
# =============================================================================
# Standard grf OOB cross-fitting holds out individual observations, not
# geographic blocks. Under spatial autocorrelation, a unit's spatial neighbours
# remain in its training set, allowing nuisance forests to exploit spatial
# proximity — producing over-optimistic residuals (Y - Y.hat) and (W - W.hat).
#
# The spatial CV version uses K=5 geographic blocks (k-means on coordinates).
# For each fold k, nuisance forests are trained on all blocks EXCEPT k, then
# predict for held-out block k. This ensures no spatial neighbours of any
# test unit appear in its nuisance training set.
#
# The spatially honest Y.hat and W.hat are then passed to causal_forest()
# which is trained on the full data as usual.

cat("--- Section 5: Spatially CV causal forest (building permits) ---\n")

fit_causal_forest_spcv <- function(X, Y, W, folds, label, seed) {
  n     <- nrow(X)
  K     <- max(folds)
  Y_hat <- numeric(n)
  W_hat <- numeric(n)

  for (k in seq_len(K)) {
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)

    set.seed(seed)
    rf_Y_k <- regression_forest(X[train_idx, , drop = FALSE], Y[train_idx], seed = seed)
    set.seed(seed)
    rf_W_k <- regression_forest(X[train_idx, , drop = FALSE], W[train_idx], seed = seed)

    Y_hat[test_idx] <- predict(rf_Y_k, newdata = X[test_idx, , drop = FALSE])$predictions
    W_hat[test_idx] <- predict(rf_W_k, newdata = X[test_idx, , drop = FALSE])$predictions
  }

  set.seed(seed)
  cf  <- causal_forest(X, Y, W, Y.hat = Y_hat, W.hat = W_hat, seed = seed)
  ate <- average_treatment_effect(cf, target.sample = "all")
  tau <- predict(cf)$predictions

  cat(sprintf("  [%s]: ATE = %.4f  SE = %.4f\n", label, ate[1], ate[2]))
  list(model = cf, ate = ate, tau = tau, label = paste0(label, " [Spatial CV]"))
}

bp_scf1 <- fit_causal_forest_spcv(X_bp_slx, Y_bp, T_bin, sp_folds, "Binary T",       seed = MASTER_SEED)
bp_scf2 <- fit_causal_forest_spcv(X_bp_slx,  Y_bp, T_twt, sp_folds, "T+WT combined",  seed = MASTER_SEED)
bp_scf3 <- fit_causal_forest_spcv(X_bp_slx,  Y_bp, T_dec, sp_folds, "Distance decay", seed = MASTER_SEED)

cat("\n")


# =============================================================================
# SECTION 6: ATE COMPARISON PLOT — BUILDING PERMITS
# =============================================================================

cat("--- Section 6: ATE comparison plot ---\n")

spec_levels <- c("Binary T", "T+WT\ncombined", "Distance\ndecay")

build_ate_row <- function(obj, cv, spec_idx) {
  data.frame(
    label    = obj$label,
    estimate = obj$ate[["estimate"]],
    std.err  = obj$ate[["std.err"]],
    cv       = cv,
    spec     = factor(spec_levels[spec_idx], levels = spec_levels),
    stringsAsFactors = FALSE
  )
}

ate_bp_df <- rbind(
  build_ate_row(bp_cf1,  "Standard",   1),
  build_ate_row(bp_cf2,  "Standard",   2),
  build_ate_row(bp_cf3,  "Standard",   3),
  build_ate_row(bp_scf1, "Spatial CV", 1),
  build_ate_row(bp_scf2, "Spatial CV", 2),
  build_ate_row(bp_scf3, "Spatial CV", 3)
)

p_ate_bp <- ggplot(ate_bp_df, aes(x = spec, y = estimate, colour = cv, shape = cv)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.err,
                    ymax = estimate + 1.96 * std.err),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Standard" = "#D6604D", "Spatial CV" = "#2166AC"),
                      name = "Nuisance estimation") +
  scale_shape_manual(values = c("Standard" = 16, "Spatial CV" = 17),
                     name = "Nuisance estimation") +
  labs(title = "Causal Forest ATE by Treatment Specification",
       subtitle = "Standard (OOB) vs Spatially Cross-Validated Nuisance Estimation",
       x = "Treatment specification", y = "ATE (standardised building permit density)",
       caption = "Error bars = 95% CI") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())

print(p_ate_bp)
ggsave(file.path(OUTPUT_DIR, "p_ate_bp.png"), plot = p_ate_bp,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("ATE plot printed.\n\n")


# =============================================================================
# SECTION 7: DEEP DIVE — SPATIAL CV DISTANCE DECAY MODEL
# =============================================================================
# Primary model: spatial CV causal forest with distance decay treatment (bp_scf3).
# Spatial CV nuisance estimation addresses spatial leakage in the OOB stage.

cat("--- Section 7: Spatial CV distance decay model deep dive ---\n")

best_bp_cf  <- bp_scf3
tau_bp_best <- best_bp_cf$tau
cat("Model:", best_bp_cf$label, "\n\n")


# ---- 7a. Permutation Variable Importance ------------------------------------
cat("--- 7a. Variable importance ---\n")

vi_bp_df <- data.frame(
  variable   = colnames(best_bp_cf$model$X.orig),
  importance = as.numeric(variable_importance(best_bp_cf$model))
) %>% arrange(desc(importance))

cat("Top 10 variables:\n")
print(head(vi_bp_df, 10), row.names = FALSE)

p_vi_bp <- ggplot(vi_bp_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "#4DAC26") +
  coord_flip() +
  labs(title = paste("Variable Importance —", best_bp_cf$label),
       x = NULL, y = "Permutation importance") +
  theme_minimal(base_size = 11)
print(p_vi_bp)
ggsave(file.path(OUTPUT_DIR, "p_vi_bp.png"), plot = p_vi_bp,
       width = 7, height = 5, dpi = 300, bg = "white")


# ---- 7b. CATE Heterogeneity Plots -------------------------------------------
cat("\n--- 7b. CATE heterogeneity ---\n")

# Top 7 non-lag variables by importance + distance = 8 panels.
top7_vars   <- vi_bp_df %>%
  filter(!grepl("_LAG$", variable)) %>%
  head(7) %>%
  pull(variable)
het_bp_vars   <- c(top7_vars, "distance")
het_bp_labels <- ifelse(het_bp_vars == "distance", "DIST (distance to trail access points)", het_bp_vars)

tau_ylim <- range(tau_bp_best, na.rm = TRUE) +
  c(-1, 1) * diff(range(tau_bp_best, na.rm = TRUE)) * 0.05

het_bp_plots <- mapply(function(v, lbl) {
  ggplot(data.frame(x = C_CA[[v]], tau = tau_bp_best), aes(x = x, y = tau)) +
    geom_point(alpha = 0.35, size = 0.8, colour = "#666666") +
    geom_smooth(method = "lm", se = TRUE, colour = "#4DAC26", fill = "#B8E186") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_cartesian(ylim = tau_ylim) +
    labs(title = lbl, x = lbl, y = expression(hat(tau))) +
    theme_minimal(base_size = 9)
}, het_bp_vars, het_bp_labels, SIMPLIFY = FALSE)

p_het_bp <- cowplot::plot_grid(plotlist = het_bp_plots, ncol = 4)
print(p_het_bp)
ggsave(file.path(OUTPUT_DIR, "p_het_bp.png"), plot = p_het_bp,
       width = 14, height = 6, dpi = 300, bg = "white")


# ---- 7c. Best Linear Projection ---------------------------------------------
cat("\n--- 7c. Best linear projection ---\n")
print(best_linear_projection(best_bp_cf$model, X_bp_slx))


# ---- 7d. Calibration Test ---------------------------------------------------
cat("\n--- 7d. Calibration test ---\n")
print(test_calibration(best_bp_cf$model))


# ---- 7e. Leaflet Maps (CATE and building permit density) --------------------
cat("\n--- 7e. Leaflet maps ---\n")

C_CAs_wgs                <- st_transform(C_CAs, crs = 4326)
C_CAs_wgs$tau_bp_best    <- tau_bp_best
C_CAs_wgs$Y_bp_density   <- Y_bp_density

# SD-based breaks: mean ± 0.5, 1, 1.5 SD, clipped to actual data range
make_sd_breaks <- function(x, n_sd = 3, step = 0.5) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  brks <- m + seq(-n_sd, n_sd, by = step) * s
  brks <- brks[brks >= min(x, na.rm = TRUE) & brks <= max(x, na.rm = TRUE)]
  unique(c(min(x, na.rm = TRUE), brks, max(x, na.rm = TRUE)))
}

# CATE: sequential blues, SD breaks
tau_breaks  <- make_sd_breaks(tau_bp_best)
pal_bp_cate <- colorBin("YlOrRd", domain = tau_bp_best, bins = tau_breaks)

# Density: YlGn, SD breaks (floor negative breaks at 0)
dens_breaks <- pmax(0, make_sd_breaks(Y_bp_density))
dens_breaks <- unique(dens_breaks)
pal_bp_dens <- colorBin("YlGn", domain = Y_bp_density, bins = dens_breaks)

map_bp_cate <- leaflet(C_CAs_wgs) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(fillColor = ~pal_bp_cate(tau_bp_best), fillOpacity = 0.5,
              weight = 0.5, color = "#FFFFFF",
              popup = ~paste0("CATE: ", round(tau_bp_best, 3))) %>%
  addLegend("bottomright", pal = pal_bp_cate, values = ~tau_bp_best,
            labFormat = labelFormat(digits = 2),
            title = "CATE (std.)", opacity = 0.9)

map_bp_dens <- leaflet(C_CAs_wgs) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(fillColor = ~pal_bp_dens(Y_bp_density), fillOpacity = 0.5,
              weight = 0.5, color = "#FFFFFF",
              popup = ~paste0("BP density (per km²): ", round(Y_bp_density, 1))) %>%
  addLegend("bottomright", pal = pal_bp_dens, values = ~Y_bp_density,
            labFormat = labelFormat(digits = 1),
            title = "Building permits\nper km² (2015–2019)", opacity = 0.9)

if (has_leafsync) leafsync::sync(map_bp_cate, map_bp_dens) else {
  print(map_bp_cate)
  print(map_bp_dens)
}


# =============================================================================
# SECTION 8: SAVE OUTPUTS
# =============================================================================

output_sf <- C_CAs %>%
  mutate(
    tau_bp_cf_bin  = bp_cf1$tau,
    tau_bp_cf_twt  = bp_cf2$tau,
    tau_bp_cf_dec  = bp_cf3$tau,
    tau_bp_scf_bin = bp_scf1$tau,
    tau_bp_scf_twt = bp_scf2$tau,
    tau_bp_scf_dec = bp_scf3$tau,
    tau_bp_best    = tau_bp_best,
    bp_density     = Y_bp_density
  )

st_write(output_sf, file.path(OUTPUT_DIR, "Chicago_CausalML_BP_Results.shp"), delete_dsn = TRUE)

cat("\n=====================================================\n")
cat("Script complete — Building permits application.\n")
cat("  Standard CF:    bp_cf1 (Binary T), bp_cf2 (T+WT), bp_cf3 (Distance decay)\n")
cat("  Spatial CV CF:  bp_scf1, bp_scf2, bp_scf3\n")
cat("  Deep dive:      Spatial CV distance decay (bp_scf3)\n")
cat("=====================================================\n")
