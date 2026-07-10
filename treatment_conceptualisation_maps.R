library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(ggspatial)

# ── Grid and treatment source ─────────────────────────────────────────────────
grid <- expand.grid(x = 1:10, y = 1:10)
tx <- 5; ty <- 5

grid <- grid %>%
  mutate(
    dist      = sqrt((x - tx)^2 + (y - ty)^2),
    chebyshev = pmax(abs(x - tx), abs(y - ty))
  )

# ── Treatment buffer: T=1 within queen-contiguity (3×3) ──────────────────────
grid <- grid %>% mutate(T_buf = as.numeric(chebyshev <= 1))

# ── Row-standardised WT (queen contiguity) ────────────────────────────────────
compute_WT <- function(T_vec, grid) {
  WT <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    nb <- which(
      abs(grid$x - grid$x[i]) <= 1 &
      abs(grid$y - grid$y[i]) <= 1 &
      !(grid$x == grid$x[i] & grid$y == grid$y[i])
    )
    if (length(nb) > 0) WT[i] <- mean(T_vec[nb])
  }
  WT
}

WT_buf <- compute_WT(grid$T_buf, grid)
grid$WT_buf <- WT_buf

# ── Variables for each panel ──────────────────────────────────────────────────

# (1) Binary
grid$T_binary <- grid$T_buf

# (2a) T and WT as separate categorical rings
grid <- grid %>%
  mutate(
    T_WT_cat = case_when(
      T_buf == 1               ~ "Direct (T)",
      WT_buf > 0 & T_buf == 0  ~ "Spillover (WT)",
      TRUE                     ~ "Untreated"
    ),
    T_WT_cat = factor(T_WT_cat, levels = c("Direct (T)", "Spillover (WT)", "Untreated"))
  )

# (2b) Distance rings — 2 rings only
grid <- grid %>%
  mutate(
    ring2 = case_when(
      T_buf == 1  ~ "T (buffer)",
      dist <= 3.5 ~ "Ring 1",
      dist <= 5.5 ~ "Ring 2",
      TRUE        ~ "Untreated"
    ),
    ring2 = factor(ring2, levels = c("T (buffer)", "Ring 1", "Ring 2", "Untreated"))
  )

# (3) Continuous total effect: T + WT
grid$T_total <- grid$T_buf + WT_buf

# (4) Continuous distance decay
alpha <- 1.2; beta <- 3.5
grid <- grid %>%
  mutate(T_decay = 1 / (1 + exp(alpha * (dist - beta))))

# ── Themes ────────────────────────────────────────────────────────────────────
map_theme <- function(title_size = 10) {
  theme_minimal(base_size = 10) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = title_size, hjust = 0.5),
    plot.subtitle    = element_text(size = 8.5, hjust = 0.5, colour = "grey40"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7.5),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
}

source_point <- geom_point(
  data = data.frame(x = tx, y = ty),
  aes(x = x, y = y), inherit.aes = FALSE,
  shape = 4, size = 3.5, colour = "black", stroke = 1.8
)

treatment_blue <- "#2171B5"

# ── (1) Binary ────────────────────────────────────────────────────────────────
p1 <- ggplot(grid, aes(x = x, y = y, fill = T_binary)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  source_point +
  scale_fill_gradient(low = "white", high = treatment_blue,
                      name = expression(T[t]),
                      breaks = c(0, 1), labels = c("0", "1")) +
  coord_equal() +
  labs(title = "(1) Binary",
       subtitle = expression(T[t] == 1 ~ "within buffer")) +
  map_theme(title_size = 11)

# ── (2a) T direct (solid) / WT spillover (continuous) ────────────────────────
# WT has continuous values depending on how many buffer neighbours each cell has;
# visualise as gradient (white→light blue) while T buffer stays solid dark blue.
p2a <- ggplot(grid, aes(x = x, y = y)) +
  # Non-buffer cells: continuous WT value
  geom_tile(data = dplyr::filter(grid, T_buf == 0),
            aes(fill = WT_buf), colour = "grey85", linewidth = 0.3) +
  # Buffer cells: solid treatment blue (T = 1, separate variable)
  geom_tile(data = dplyr::filter(grid, T_buf == 1),
            fill = treatment_blue, colour = "grey85", linewidth = 0.3) +
  source_point +
  scale_fill_gradient(low = "white", high = treatment_blue,
                      name = expression(WT[t]), limits = c(0, 1)) +
  coord_equal() +
  labs(subtitle = expression(T[t] ~ "direct (solid);" ~ WT[t] ~ "continuous")) +
  map_theme(title_size = 9)

# ── (2b) Distance rings (2 rings) ────────────────────────────────────────────
ring_colours <- c(
  "T (buffer)" = treatment_blue,
  "Ring 1"     = "#6BAED6",
  "Ring 2"     = "#C6DBEF",
  "Untreated"  = "white"
)

p2b <- ggplot(grid, aes(x = x, y = y, fill = ring2)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  source_point +
  scale_fill_manual(values = ring_colours, name = NULL,
                    labels = c(expression(T[t]),
                               expression(T[td1]),
                               expression(T[td2]),
                               "Untreated")) +
  coord_equal() +
  labs(subtitle = expression(T[t] * "," ~ T[td1] * "," ~ T[td2])) +
  map_theme(title_size = 9) +
  guides(fill = guide_legend(nrow = 1))

# ── Nest (2a) and (2b) under shared Distance rings heading ───────────────────
p2_title <- ggplot() +
  labs(title = '(2) Distance "rings"') +
  theme_void() +
  theme(
    plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.background = element_rect(fill = "white", colour = NA)
  )

p2_group <- p2_title / (p2a | p2b) + plot_layout(heights = c(0.12, 1))

# ── (3) Continuous total effect ───────────────────────────────────────────────
p3 <- ggplot(grid, aes(x = x, y = y, fill = T_total)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  source_point +
  scale_fill_gradientn(
    colours = c("white", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits  = c(0, 2),
    name    = expression(T[t] + WT[t])
  ) +
  coord_equal() +
  labs(title = "(3) Continuous total effect",
       subtitle = expression(T[t] + WT[t] ~ "(combined)")) +
  map_theme(title_size = 11)

# ── (4) Continuous distance decay ────────────────────────────────────────────
p4 <- ggplot(grid, aes(x = x, y = y, fill = T_decay)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  source_point +
  scale_fill_gradient(low = "white", high = treatment_blue,
                      name = expression(kappa(d[i]*";"~theta)),
                      limits = c(0, 1)) +
  coord_equal() +
  labs(title = "(4) Continuous distance decay",
       subtitle = expression(T[t] == kappa(d[i]*";"~theta))) +
  map_theme(title_size = 11)

# ── Combine: (1) row / (2) row / (3)+(4) row ─────────────────────────────────
combined <- p1 / p2_group / (p3 | p4) +
  plot_layout(heights = c(1, 1.1, 1)) +
  plot_annotation(
    caption = "× marks treatment source; buffer = queen-contiguity neighbourhood (3×3 cells)",
    theme   = theme(
      plot.caption    = element_text(size = 8, colour = "grey50", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave("treatment_conceptualisation_maps.png", plot = combined,
       width = 8, height = 14, dpi = 300, bg = "white")

cat("Saved: treatment_conceptualisation_maps.png\n")


# ═════════════════════════════════════════════════════════════════════════════
# STUDY AREA OVERVIEW MAP — The 606 / Bloomingdale Trail
# Block groups (light gray border, transparent), community areas (thick black
# border with labels), access points, 500m buffer (dashed yellow)
# ═════════════════════════════════════════════════════════════════════════════

BG_SHP     <- "/Users/kevincredit/Library/CloudStorage/Dropbox/Packages/SArf paper/Book chapter/Spatial_Causal_ML/Spatial Causal ML/data/Chicago_BGs_Covariates.shp"
CA_SHP     <- "/Users/kevincredit/Library/CloudStorage/Dropbox/Courses/Casual ML/606 causal analysis/Data/CA_606.shp"
POINTS_SHP <- "/Users/kevincredit/Library/CloudStorage/Dropbox/Courses/Casual ML/606 causal analysis/Data/Access_Points_606_2163.shp"
BUFFER_SHP <- "/Users/kevincredit/Library/CloudStorage/Dropbox/Courses/Casual ML/606 causal analysis/Data/Access_Points_606_2163_500m_Buffer2.shp"

bg_sf  <- st_read(BG_SHP,     quiet = TRUE) %>% st_transform(3857)
ca_sf  <- st_read(CA_SHP,     quiet = TRUE) %>% st_transform(3857)
pts_sf <- st_read(POINTS_SHP, quiet = TRUE) %>% st_transform(3857)
buf_sf <- st_read(BUFFER_SHP, quiet = TRUE) %>% st_transform(3857)

ca_cents <- suppressWarnings(st_centroid(ca_sf))

ca_bbox <- st_bbox(ca_sf)
x_pad   <- (ca_bbox["xmax"] - ca_bbox["xmin"]) * 0.06
y_pad   <- (ca_bbox["ymax"] - ca_bbox["ymin"]) * 0.06

p_study <- ggplot() +
  annotation_map_tile(type = "osm", zoom = 14, alpha = 0.55) +
  geom_sf(data = bg_sf,  fill = NA, colour = "gray70",  linewidth = 0.2) +
  geom_sf(data = ca_sf,  fill = NA, colour = "black",   linewidth = 1.1) +
  geom_sf(data = buf_sf, fill = NA, colour = "yellow",  linewidth = 0.9,
          linetype = "dashed") +
  geom_sf(data = pts_sf, colour = "black", fill = "#2171B5",
          shape = 21, size = 2.8, stroke = 1.0) +
  geom_sf_text(data = ca_cents,
               aes(label = community),
               size = 2.8, fontface = "bold", colour = "black") +
  coord_sf(
    xlim = c(ca_bbox["xmin"] - x_pad, ca_bbox["xmax"] + x_pad),
    ylim = c(ca_bbox["ymin"] - y_pad, ca_bbox["ymax"] + y_pad)
  ) +
  labs(
    title    = "The 606 / Bloomingdale Trail: Study Area",
    subtitle = "Community areas (bold border) · Access points · 500m buffer (dashed)"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40", hjust = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave("study_area_map.png", plot = p_study,
       width = 7, height = 6, dpi = 300, bg = "white")

cat("Saved: study_area_map.png\n")
