library(dagitty)
library(ggdag)
library(ggplot2)
library(patchwork)

# ── Helper: parse DAGitty string, flip y so layout matches the online viewer ──
parse_dag <- function(dag_str) {
  dag <- dagitty(dag_str)
  coords <- coordinates(dag)   # reads pos= attributes; DAGitty y is y-down
  coords$y <- -coords$y        # flip to ggplot2 y-up so layout matches viewer
  coordinates(dag) <- coords
  dag
}

# ── Time-indexed display labels (plotmath for subscripts) ─────────────────────
time_label <- function(name) {
  dplyr::case_when(
    name %in% c("Ti", "Ti_t")    ~ "T[it]",
    name %in% c("Tj", "Tj_t")    ~ "T[jt]",
    name %in% c("Xi", "Xi_t-1")  ~ "X[it-1]",
    name %in% c("Xj", "Xj_t-1")  ~ "X[jt-1]",
    name %in% c("Yi", "Yi_t+1")  ~ "Y[it+1]",
    name %in% c("Yj", "Yj_t+1")  ~ "Y[jt+1]",
    name == "Ui"                  ~ "U[i]",
    name == "Uj"                  ~ "U[j]",
    TRUE                          ~ name
  )
}

# ── Node colour lookup (by name) ──────────────────────────────────────────────
node_fill <- function(name) {
  dplyr::case_when(
    name %in% c("Ti", "Tj", "Ti_t", "Tj_t")       ~ "#7BBD5E",
    name %in% c("Yi", "Yi_t+1")                    ~ "#3A9EC2",
    name %in% c("Xi", "Xj", "Xi_t-1", "Xj_t-1")  ~ "white",
    name %in% c("Ui", "Uj")                        ~ "#CCCCCC",
    TRUE                                            ~ "#CCCCCC"
  )
}

# ── Helper: styled ggdag plot ─────────────────────────────────────────────────
# parse_labels = TRUE  → plotmath expressions (subscripts) for typology DAGs
# parse_labels = FALSE → plain wrapped text for application DAGs with long names
plot_dag <- function(dag, title = NULL, extra_fills = NULL,
                     parse_labels = TRUE, node_size = 16, label_size = 2.8,
                     wrap_width = 13) {
  td <- tidy_dagitty(dag)

  # Build fill lookup from distinct nodes only
  node_df <- td$data %>%
    dplyr::distinct(name, x, y) %>%
    dplyr::filter(!is.na(name)) %>%
    dplyr::mutate(
      fill_col = node_fill(name),
      disp_label = if (parse_labels) {
        time_label(name)
      } else {
        vapply(name, function(n)
          paste(strwrap(n, width = wrap_width), collapse = "\n"),
          character(1), USE.NAMES = FALSE)
      }
    )

  fills <- setNames(node_df$fill_col, node_df$name)
  if (!is.null(extra_fills)) fills[names(extra_fills)] <- extra_fills

  ggplot(td, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(
      arrow_directed   = grid::arrow(length = grid::unit(0.18, "cm"),
                                     type = "closed"),
      arrow_bidirected = grid::arrow(length = grid::unit(0.15, "cm"),
                                     type = "open", ends = "both")
    ) +
    geom_point(data = node_df,
               aes(x = x, y = y, fill = name, xend = NULL, yend = NULL),
               shape = 21, size = node_size, colour = "black", stroke = 1.2) +
    geom_text(data = node_df,
              aes(x = x, y = y, label = disp_label,
                  xend = NULL, yend = NULL),
              size = label_size, colour = "black", parse = parse_labels,
              lineheight = 0.9) +
    scale_fill_manual(values = fills, guide = "none") +
    labs(title = title) +
    theme_dag(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# ═════════════════════════════════════════════════════════════════════════════
# CELL 1 — Non-SAR, No interference:  y = τT + Xβ + ε
# ═════════════════════════════════════════════════════════════════════════════
dag_c1 <- parse_dag('dag {
  Ti_t    [exposure, pos="-1.5,0.0"]
  "Xi_t-1" [adjusted, pos="-0.5,-1.0"]
  "Yi_t+1" [outcome,  pos="1.0,0.0"]
  Ti_t -> "Yi_t+1"
  "Xi_t-1" -> Ti_t
  "Xi_t-1" -> "Yi_t+1"
}')

p1 <- plot_dag(dag_c1, title = "Standard: Non-SAR, No interference\ny_(t+1) = \u03c4T_t + X_(t-1)\u03b2 + \u03b5")

# ═════════════════════════════════════════════════════════════════════════════
# CELL 2 — SAR, No interference:  y = (I - ρW)⁻¹(τT + Xβ + ε)
# ═════════════════════════════════════════════════════════════════════════════
dag_c2 <- parse_dag('dag {
  Ti_t     [exposure, pos="-1.952,0.211"]
  Tj_t     [exposure, pos="-1.459,0.664"]
  "Xi_t-1" [adjusted, pos="-0.638,-0.762"]
  "Xj_t-1" [adjusted, pos="-0.631,-0.234"]
  "Yi_t+1" [outcome,  pos="1.035,0.206"]
  "Yj_t+1" [pos="0.253,0.656"]
  Ti_t -> "Yi_t+1"
  Ti_t -> "Yj_t+1"
  Ti_t <-> Tj_t
  Tj_t -> "Yi_t+1"
  Tj_t -> "Yj_t+1"
  "Xi_t-1" -> Ti_t
  "Xi_t-1" -> "Yi_t+1"
  "Xi_t-1" -> "Yj_t+1"
  "Xi_t-1" <-> "Xj_t-1"
  "Xj_t-1" -> Tj_t
  "Xj_t-1" -> "Yi_t+1"
  "Xj_t-1" -> "Yj_t+1"
  "Yi_t+1" <-> "Yj_t+1"
}')

p2 <- plot_dag(dag_c2,
               title = "SAR\ny_(t+1) = (I \u2212 \u03c1W)\u207b\u00b9(\u03c4T_t + X_(t-1)\u03b2 + \u03b5)\ny_(t+1) = (I \u2212 \u03c1W)\u207b\u00b9(\u03b4WT_t + \u03c4T_t + X_(t-1)\u03b2 + \u03b5)")

# ═════════════════════════════════════════════════════════════════════════════
# CELL 3 — Non-SAR, Interference:  y = δWT + τT + Xβ + ε
# ═════════════════════════════════════════════════════════════════════════════
dag_c3 <- parse_dag('dag {
  "Xi_t-1" [adjusted, pos="-0.638,-0.762"]
  "Xj_t-1" [adjusted, pos="-0.614,-0.559"]
  "Yi_t+1" [outcome,  pos="0.980,0.215"]
  "Yj_t+1" [pos="0.822,0.452"]
  Ti_t [exposure, pos="-1.952,0.211"]
  Tj_t [exposure, pos="-1.502,0.386"]
  "Xi_t-1" -> "Yi_t+1"
  "Xi_t-1" -> Ti_t
  "Xj_t-1" -> "Yj_t+1"
  "Xj_t-1" -> Tj_t
  Ti_t -> "Yi_t+1"
  Ti_t -> "Yj_t+1"
  Ti_t <-> Tj_t
  Tj_t -> "Yi_t+1"
  Tj_t -> "Yj_t+1"
}')

p3 <- plot_dag(dag_c3,
               title = "Interference: Non-SAR, Interference\ny_(t+1) = \u03b4WT_t + \u03c4T_t + X_(t-1)\u03b2 + \u03b5",
               extra_fills = c("Yj_t+1" = "white"))

# ═════════════════════════════════════════════════════════════════════════════
# SPATIAL ERROR — Non-SAR, No interference, spatially autocorrelated errors
# y_(t+1) = τT_t + X_(t-1)β + (I−λW)⁻¹ε
#
# Key DAG features:
#   - Same structural equation as Cell 1 (T→Y, X→T, X→Y)
#   - Unobserved U_i, U_j drive the spatially correlated errors
#   - U_i ↔ U_j: bidirectional association (spatial error correlation)
#   - NO Y_j → Y_i (no SAR feedback — distinguishes from SAR)
#   - NO T_j → Y_i (no interference)
#   - WY_{t-1} is a legitimate proxy because it absorbs W·U without
#     introducing collider bias (T has no path to WY through SAR)
# ═════════════════════════════════════════════════════════════════════════════
dag_sem <- parse_dag('dag {
  Ti_t     [exposure, pos="-1.5,-0.5"]
  "Xi_t-1" [adjusted, pos="-0.5,-1.2"]
  "Yi_t+1" [outcome,  pos="1.2,-0.5"]
  "Yj_t+1" [pos="1.2,0.5"]
  Tj_t     [exposure, pos="-1.5,0.5"]
  "Xj_t-1" [adjusted, pos="-0.5,1.2"]
  Ui       [latent,   pos="0.3,0.0"]
  Uj       [latent,   pos="0.3,0.8"]
  Ti_t     -> "Yi_t+1"
  "Xi_t-1" -> Ti_t
  "Xi_t-1" -> "Yi_t+1"
  Tj_t     -> "Yj_t+1"
  "Xj_t-1" -> Tj_t
  "Xj_t-1" -> "Yj_t+1"
  Ui       -> "Yi_t+1"
  Uj       -> "Yj_t+1"
  Ui       <-> Uj
}')

p_sem <- plot_dag(dag_sem,
                  title = "Spatial error structure: Non-SAR, No interference\ny_(t+1) = \u03c4T_t + X_(t-1)\u03b2 + (I\u2212\u03bbW)\u207b\u00b9\u03b5",
                  extra_fills = c("Yj_t+1" = "white"))

# ═════════════════════════════════════════════════════════════════════════════
# SPATIAL CONFOUNDER — Ui causes Xi, Ti, Yi; Ui <-> Uj (bidirectional)
# Same layout as dag_sem, with additional Ui→Xi, Ui→Ti (and Uj→Xj, Uj→Tj)
#
# Key DAG features:
#   - Ui (pink) is unobserved confounder: Ui → Xi, Ti, Yi
#   - Ui <-> Uj (bidirectional): spatially correlated unobserved features
#   - Uj (grey) mirrors Ui structure for j-unit: Uj → Xj, Tj, Yj
#   - Open backdoor paths Ti ← Ui → Yi
#   - WY is an imperfect proxy for Ui through Uj → Yj, but does NOT formally
#     block the backdoor and risks collider bias if SAR is also present
# ═════════════════════════════════════════════════════════════════════════════
dag_conf <- parse_dag('dag {
  Ti_t     [exposure, pos="-1.5,-0.5"]
  "Xi_t-1" [adjusted, pos="-0.5,-1.2"]
  "Yi_t+1" [outcome,  pos="1.2,-0.5"]
  "Yj_t+1" [pos="1.2,0.5"]
  Tj_t     [exposure, pos="-1.5,0.5"]
  "Xj_t-1" [adjusted, pos="-0.5,1.2"]
  Ui       [latent,   pos="0.3,0.0"]
  Uj       [latent,   pos="0.3,0.8"]
  Ti_t     -> "Yi_t+1"
  "Xi_t-1" -> Ti_t
  "Xi_t-1" -> "Yi_t+1"
  Tj_t     -> "Yj_t+1"
  "Xj_t-1" -> Tj_t
  "Xj_t-1" -> "Yj_t+1"
  Ui       -> "Yi_t+1"
  Ui       -> "Xi_t-1"
  Ui       -> Ti_t
  Ui       <-> Uj
  Uj       -> "Yj_t+1"
  Uj       -> "Xj_t-1"
  Uj       -> Tj_t
}')

p_conf <- plot_dag(dag_conf,
                   title = "Unobserved spatial confounding\ny_(t+1) = \u03c4T_t + X_(t-1)\u03b2 + \u03b3U_i + \u03b5",
                   extra_fills = c("Ui"     = "#F08080",
                                   "Yj_t+1" = "white"))

# ── Save individual panels ────────────────────────────────────────────────────
ggsave("dag_conf_spatial_confounder.png", plot = p_conf,
       width = 6, height = 5, dpi = 300, bg = "white")
cat("Saved: dag_conf_spatial_confounder.png\n")

ggsave("dag_sem_spatial_error.png", plot = p_sem,
       width = 6, height = 5, dpi = 300, bg = "white")
cat("Saved: dag_sem_spatial_error.png\n")

ggsave("dag_panel1_standard.png", plot = p1,
       width = 5, height = 5, dpi = 300, bg = "white")

ggsave("dag_panel2_sar.png", plot = p2,
       width = 6, height = 5, dpi = 300, bg = "white")

ggsave("dag_panel3_interference.png", plot = p3,
       width = 6, height = 5, dpi = 300, bg = "white")

cat("Saved: dag_panel1_standard.png\n")
cat("Saved: dag_panel2_sar.png\n")
cat("Saved: dag_panel3_interference.png\n")

# ── 3-panel vertical typology figure (Figure X.1) ────────────────────────────
# Top    = Panel 1: Standard (Non-SAR, No interference)
# Middle = Panel 2: SAR (= SAR+I in reduced form)
# Bottom = Panel 3: Interference (Non-SAR, Interference)

combined <- p1 / p2 / p3

ggsave("dag_3panel_typology.png", plot = combined,
       width = 7, height = 15, dpi = 300, bg = "white")

cat("Saved: dag_3panel_typology.png\n")

# ── Figure X.2: 2-panel spatial error structure (side by side) ───────────────
# Left  = classical SEM (U only causes Y, Ui <-> Uj correlated errors)
# Right = spatial confounder (Ui causes X, T, Y; Ui → Uj directed)

fig_x2 <- p_sem + p_conf +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "1",
                  tag_prefix = "(",
                  tag_suffix = ")")

ggsave("dag_fig_x2_spatial_error.png", plot = fig_x2,
       width = 12, height = 5, dpi = 300, bg = "white")

cat("Saved: dag_fig_x2_spatial_error.png\n")

cat("Saved: dag_3panel_typology.png\n")

# ═════════════════════════════════════════════════════════════════════════════
# APPLICATION DAG — 606 / Bloomingdale Trail: building permits
# Exposure: 606 Construction (2015)
# Outcome:  Building permit density (2017)
# Adjusted: socioeconomic, demographic, land-use covariates (2010/2012)
# ═════════════════════════════════════════════════════════════════════════════
dag_app <- parse_dag('dag {
  "606 Construction (2015)"        [exposure, pos="-0.023,1.616"]
  "Building permit density (2015-2019)"  [outcome,  pos="1.110,1.594"]
  "Building permit density (2010-2014)" [adjusted, pos="1.402,1.061"]
  "Business zoning (2012)"         [adjusted, pos="1.204,0.378"]
  "Education (2010)"               [adjusted, pos="-1.805,0.811"]
  "Employment and commuting (2010)"  [adjusted, pos="-1.685,1.281"]
  "Income (2010)"                  [adjusted, pos="-2.030,1.563"]
  "Median age (2010)"              [adjusted, pos="-1.227,-0.048"]
  "Race & ethnicity (2010)"        [adjusted, pos="-1.767,0.222"]
  "Occupancy and tenure (2010)"                  [adjusted, pos="0.405,-0.028"]
  "606 Construction (2015)"        -> "Building permit density (2015-2019)"
  "Building permit density (2010-2014)" -> "606 Construction (2015)"
  "Building permit density (2010-2014)" -> "Building permit density (2015-2019)"
  "Building permit density (2010-2014)" -> "Business zoning (2012)"
  "Business zoning (2012)"         -> "606 Construction (2015)"
  "Business zoning (2012)"         -> "Building permit density (2015-2019)"
  "Education (2010)"               -> "606 Construction (2015)"
  "Education (2010)"               -> "Building permit density (2015-2019)"
  "Education (2010)"               -> "Employment and commuting (2010)"
  "Education (2010)"               -> "Income (2010)"
  "Race & ethnicity (2010)"        -> "Education (2010)"
  "Employment and commuting (2010)"  -> "606 Construction (2015)"
  "Employment and commuting (2010)"  -> "Building permit density (2015-2019)"
  "Employment and commuting (2010)"  -> "Business zoning (2012)"
  "Employment and commuting (2010)"  -> "Income (2010)"
  "Income (2010)"                  -> "606 Construction (2015)"
  "Income (2010)"                  -> "Building permit density (2015-2019)"
  "Income (2010)"                  -> "Occupancy and tenure (2010)"
  "Median age (2010)"              -> "606 Construction (2015)"
  "Median age (2010)"              -> "Building permit density (2015-2019)"
  "Median age (2010)"              -> "Business zoning (2012)"
  "Median age (2010)"              -> "Education (2010)"
  "Median age (2010)"              -> "Employment and commuting (2010)"
  "Median age (2010)"              -> "Income (2010)"
  "Median age (2010)"              -> "Occupancy and tenure (2010)"
  "Race & ethnicity (2010)"        -> "606 Construction (2015)"
  "Race & ethnicity (2010)"        -> "Building permit density (2015-2019)"
  "Race & ethnicity (2010)"        -> "Median age (2010)"
  "Occupancy and tenure (2010)"                  -> "606 Construction (2015)"
  "Occupancy and tenure (2010)"                  -> "Building permit density (2015-2019)"
  "Occupancy and tenure (2010)"                  -> "Business zoning (2012)"
}')

p_app <- plot_dag(
  dag_app,
  title        = "Application DAG: 606 / Bloomingdale Trail",
  parse_labels = FALSE,
  node_size    = 20,
  label_size   = 1.9,
  wrap_width   = 13,
  extra_fills  = c(
    "606 Construction (2015)"        = "#7BBD5E",
    "Building permit density (2015-2019)"  = "#3A9EC2",
    "Building permit density (2010-2014)" = "white",
    "Business zoning (2012)"         = "white",
    "Education (2010)"               = "white",
    "Employment and commuting (2010)"  = "white",
    "Income (2010)"                  = "white",
    "Median age (2010)"              = "white",
    "Race & ethnicity (2010)"        = "white",
    "Occupancy and tenure (2010)"                  = "white"
  )
)

ggsave("dag_app_606.png", plot = p_app,
       width = 10, height = 7, dpi = 300, bg = "white")

cat("Saved: dag_app_606.png\n")
