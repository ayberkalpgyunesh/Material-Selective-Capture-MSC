# ==============================================================================
# MSC Figure 6 (main) and Supplementary Figure S9 — three-dataset panels
# Version 6
# ==============================================================================
# Input:  the 10_r_input folders produced by the three harmonised notebooks
#         (Middleton, Perez, Tushev). They may sit in one folder or in three.
#
# Usage:  Rscript MSC_Figure6_and_FigureS9_panels_v6.r /path/to/inputs
#
# This script does not simulate and does not re-fit any model. It reads the
# canonical plotting tables, harmonises their column names, assembles the main
# and supplementary panels, and exports every panel as PDF, PNG, SVG and TIFF
# together with the exact plotted data.
#
# Changes in version 6, after visual review of the version 5 output
# -----------------------------------------------------------------
#   1. Facet strip labels on the right-hand side of S9A and S9B were clipped
#      ("liddleton (mouse, gen(" instead of "Middleton (mouse, gene)").
#      Short strip labels, an explicit strip margin and a wider canvas fix it.
#   2. Panel 6C listed the analysis arms in the opposite order to 6D, 6E, S9A
#      and S9B. All preference and shift panels now put the primary score on
#      top.
#   3. Panels 6C, S9A and S9B show the same quantity but used three different
#      x-ranges. They now share one fixed range so that docking and productive
#      capture can be compared directly.
#
# Panels produced
# ---------------
#   Figure 6B   paired sequence-score differences, nine scores, three datasets
#   Figure 6C   absolute preference in each PSD state, primary observable
#   Figure 6D   WT-like minus ME-like shift in preference fraction
#   Figure 6E   mapping contrasts at the primary observable
#
#   Figure S9A  relative preference at physical docking, ME-like and WT-like,
#               all three PSD observables and all three analysis arms
#   Figure S9B  the same at productive capture
#   Figure S9C  paired sequence-score differences with significance stars
#   Figure S9D  WT-like minus ME-like shift in the log localized/soma rate ratio
#   Figure S9E  PSD-coordinate scaling sweep. Optional: drawn only if a scaling
#               table is present, because that sweep is not part of the
#               standard 18-condition grid.
#
# Panels carry no explanatory text. Captions are added during layout.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Install and load packages
# ------------------------------------------------------------------------------
# data.table, ggplot2 and scales are required. ragg and svglite improve raster
# and vector export but the script runs without them.

required_packages <- c("data.table", "ggplot2", "scales")

optional_packages <- c("ragg", "svglite")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

for (package_name in optional_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    try(
      install.packages(package_name, repos = "https://cloud.r-project.org"),
      silent = TRUE
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

have_ragg <- requireNamespace("ragg", quietly = TRUE)
have_svglite <- requireNamespace("svglite", quietly = TRUE)


# ------------------------------------------------------------------------------
# 2. Paths and output folders
# ------------------------------------------------------------------------------
# The first command-line argument is the search root. Every required CSV is
# located by recursive search, so the three 10_r_input folders may live side by
# side or inside a common parent.

arguments <- commandArgs(trailingOnly = TRUE)

INPUT_DIR <- if (length(arguments) >= 1 && nzchar(arguments[1])) {
  arguments[1]
} else {
  getwd()
}

OUTPUT_DIR <- file.path(INPUT_DIR, "R_output_figure6_S9_v6")

for (folder in c("pdf", "png", "svg", "tiff", "tables")) {
  dir.create(file.path(OUTPUT_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}

# Files written by an earlier run must never be picked up as input.

all_files <- list.files(
  INPUT_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

all_files <- all_files[
  !startsWith(
    normalizePath(all_files, mustWork = FALSE),
    normalizePath(OUTPUT_DIR, mustWork = FALSE)
  )
]


find_csv <- function(required_name, optional = FALSE) {

  matches <- all_files[basename(all_files) == required_name]

  if (length(matches) == 0) {
    if (optional) return(NA_character_)
    stop("Could not find required CSV under ", INPUT_DIR, ":\n  ",
         required_name, call. = FALSE)
  }

  if (length(matches) > 1) {
    matches <- matches[which.max(file.info(matches)$mtime)]
  }

  message("Reading: ", matches[[1]])
  matches[[1]]
}


# ------------------------------------------------------------------------------
# 3. Design definitions
# ------------------------------------------------------------------------------
# All naming, ordering, colour and axis decisions live here so that the
# plotting sections below contain layout only.

DPI <- 600

alpha_level <- 0.05

# Fix 3: one shared x-range for every panel that shows a preference fraction.
preference_limits <- c(0.40, 0.56)
preference_breaks <- seq(0.40, 0.55, 0.05)

dataset_prefixes <- c(
  middleton = "Middleton_2019_GSE115480",
  perez     = "Perez_2021_GSE157204",
  tushev    = "Tushev_2018_CA1_isoform"
)

dataset_labels <- c(
  Middleton_2019_GSE115480 = "Middleton (mouse, gene)",
  Perez_2021_GSE157204     = "Perez (rat, gene)",
  Tushev_2018_CA1_isoform  = "Tushev (rat, isoform)"
)

# Fix 1: short labels are used only inside vertical facet strips, where the
# long form does not fit and was previously clipped.
dataset_strip_labels <- c(
  "Middleton (mouse, gene)" = "Middleton",
  "Perez (rat, gene)"       = "Perez",
  "Tushev (rat, isoform)"   = "Tushev"
)

# Okabe-Ito palette: blue, green and vermillion remain distinguishable in
# greyscale and under the common colour-vision deficiencies.

dataset_colours <- c(
  "Middleton (mouse, gene)" = "#0072B2",
  "Perez (rat, gene)"       = "#009E73",
  "Tushev (rat, isoform)"   = "#D55E00"
)

dataset_shapes <- c(
  "Middleton (mouse, gene)" = 21,
  "Perez (rat, gene)"       = 22,
  "Tushev (rat, isoform)"   = 24
)

# Panels S9A and S9B separate the two PSD states by colour rather than by
# dataset, because there the state contrast is the message.

state_colours <- c("ME-like" = "#8C8C8C", "WT-like" = "#2166AC")

state_shapes <- c("ME-like" = 21, "WT-like" = 22)

score_order <- c(
  "primary_g_rg4",
  "g_only",
  "rg4_only",
  "rg4_residual",
  "relaxed_rg4",
  "gc_control",
  "length_control",
  "shuffled_primary",
  "random_score"
)

score_labels <- c(
  primary_g_rg4    = "G + rG4 (primary)",
  g_only           = "G content",
  rg4_only         = "rG4 density",
  rg4_residual     = "rG4 | G, length",
  relaxed_rg4      = "rG4 (relaxed)",
  gc_control       = "GC content",
  length_control   = "3'UTR length",
  shuffled_primary = "Shuffled score",
  random_score     = "Random score"
)

analysis_order <- c("primary_g_rg4", "g_only", "primary_mapping_shuffled")

analysis_labels <- c(
  primary_g_rg4            = "G + rG4 (primary)",
  g_only                   = "G content",
  primary_mapping_shuffled = "Shuffled mapping"
)

contrast_order <- c(
  "primary_minus_shuffled",
  "g_only_minus_shuffled",
  "g_only_minus_primary"
)

contrast_labels <- c(
  primary_minus_shuffled = "Primary \u2212 shuffled mapping",
  g_only_minus_shuffled  = "G content \u2212 shuffled mapping",
  g_only_minus_primary   = "G content \u2212 primary"
)

outcome_labels <- c(
  physical_docking   = "Physical docking",
  productive_capture = "Productive capture"
)

state_labels <- c(ME_like = "ME-like", WT_like = "WT-like")

# PSD anchors follow from the published wild-type versus oligomerisation-
# deficient contrasts with gain 0.55 and are identical in all three datasets.

observable_order <- c("tau2_ms", "E_modulus_kPa", "mobile_fraction")

observable_meta <- data.table(
  observable = observable_order,
  label = c(
    "Relaxation time (\u03c4\u2082)\nanchor \u00b10.83, 100 seeds",
    "Elastic modulus\nanchor \u00b10.44, 30 seeds",
    "Mobile fraction\nanchor \u00b10.38, 30 seeds"
  )
)


# Significance stars follow the convention stated in the S9 legend:
# * q < 0.05, ** q < 0.01, *** q < 0.001.

significance_stars <- function(q_value) {

  fifelse(
    is.na(q_value), "",
    fifelse(q_value < 0.001, "***",
            fifelse(q_value < 0.01, "**",
                    fifelse(q_value < alpha_level, "*", "")))
  )
}


# ------------------------------------------------------------------------------
# 4. Load the three datasets and harmonise their column names
# ------------------------------------------------------------------------------
# The notebooks differ in two harmless ways that must be absorbed here:
#
#   (a) the seed count is called 'n_seeds' in the Middleton and Perez state
#       tables but 'n' in the Tushev tables;
#   (b) Panel B stores the localized-group mean as 'mean_dendritic' for the
#       gene-level datasets and 'mean_distal' for the isoform-level dataset.
#
# Resolving both once here means that no plotting code needs to know which
# dataset it is looking at.

read_panel <- function(filename) {

  table <- fread(find_csv(filename))

  if (!"n_seeds" %in% names(table) && "n" %in% names(table)) {
    table[, n_seeds := n]
  }

  if (!"mean_localized" %in% names(table)) {
    if ("mean_distal" %in% names(table)) {
      table[, mean_localized := mean_distal]
    } else if ("mean_dendritic" %in% names(table)) {
      table[, mean_localized := mean_dendritic]
    }
  }

  table
}


read_dataset <- function(prefix) {

  list(
    sequence_effects = read_panel(
      paste0("panel_B_", prefix, "_sequence_effects.csv")
    ),
    state_preferences = read_panel(
      paste0("panel_C_", prefix, "_state_preferences.csv")
    ),
    wt_me_shifts = read_panel(
      paste0("panel_D_", prefix, "_wt_me_shifts.csv")
    ),
    mapping_contrasts = read_panel(
      paste0("panel_D_", prefix, "_mapping_contrasts.csv")
    )
  )
}


datasets <- lapply(names(dataset_prefixes), read_dataset)
names(datasets) <- names(dataset_prefixes)

sequence_effects <- rbindlist(lapply(datasets, `[[`, "sequence_effects"), fill = TRUE)
state_preferences <- rbindlist(lapply(datasets, `[[`, "state_preferences"), fill = TRUE)
wt_me_shifts <- rbindlist(lapply(datasets, `[[`, "wt_me_shifts"), fill = TRUE)
mapping_contrasts <- rbindlist(lapply(datasets, `[[`, "mapping_contrasts"), fill = TRUE)


# ------------------------------------------------------------------------------
# 5. Validate the combined input before anything is plotted
# ------------------------------------------------------------------------------
# These checks are deliberately strict. A silently missing dataset or a wrong
# seed count would otherwise produce a plausible-looking but wrong figure.

expected_datasets <- unname(dataset_prefixes)

for (table_name in c("sequence_effects", "state_preferences",
                     "wt_me_shifts", "mapping_contrasts")) {

  table <- get(table_name)

  found <- sort(unique(table$dataset))

  if (!setequal(found, expected_datasets)) {
    stop(table_name, " does not contain all three datasets.\n  found: ",
         paste(found, collapse = ", "), call. = FALSE)
  }
}

if (nrow(sequence_effects) != 3L * length(score_order)) {
  stop("Panel B must contain ", length(score_order),
       " scores per dataset; found ", nrow(sequence_effects), " rows.")
}

if (!setequal(unique(sequence_effects$score), score_order)) {
  stop("Panel B score names do not match the expected set.")
}

if (!setequal(sort(unique(wt_me_shifts$n_seeds)), c(30L, 100L))) {
  stop("Expected seed counts 30 and 100 in the shift tables.")
}

if (!setequal(unique(wt_me_shifts$observable), observable_order)) {
  stop("Shift tables do not contain the three expected PSD observables.")
}

message("Input checks passed for all three datasets.")


# ------------------------------------------------------------------------------
# 6. Shared reshaping helpers
# ------------------------------------------------------------------------------
# Every panel needs the same factor construction, so it is written once.
#
# Note on ordering: ggplot draws the first factor level at the bottom of a
# discrete axis. Passing reverse = TRUE therefore puts the primary score at the
# top, which is the convention used by every preference and shift panel.

add_dataset_label <- function(table) {

  table[, dataset_label := factor(
    dataset_labels[dataset],
    levels = unname(dataset_labels)
  )]

  # Fix 1: short form for vertical facet strips.
  table[, dataset_strip := factor(
    dataset_strip_labels[as.character(dataset_label)],
    levels = unname(dataset_strip_labels)
  )]

  table
}

add_analysis_label <- function(table, reverse = TRUE) {

  levels_used <- analysis_labels[analysis_order]
  if (reverse) levels_used <- rev(levels_used)

  table[, analysis_label := factor(analysis_labels[analysis_name],
                                   levels = levels_used)]
  table
}

add_outcome_label <- function(table) {
  table[, outcome_label := factor(outcome_labels[outcome],
                                  levels = unname(outcome_labels))]
  table
}

add_observable_label <- function(table) {
  table <- merge(table, observable_meta, by = "observable", all.x = TRUE)
  table[, observable_label := factor(label, levels = observable_meta$label)]
  table
}


# ------------------------------------------------------------------------------
# 7. Data for the main Figure 6 panels
# ------------------------------------------------------------------------------

# Figure 6B. One point per score and dataset. The estimate is the mean paired
# difference between the localized transcript and its matched somatic partner,
# in units of the frozen-panel z-scale. Zero means no difference.

panel_6b <- copy(sequence_effects)

panel_6b <- add_dataset_label(panel_6b)

panel_6b[, score_label := factor(score_labels[score],
                                 levels = rev(score_labels[score_order]))]

panel_6b[, stars := significance_stars(q_bh_empirical)]

panel_6b[, significant := q_bh_empirical < alpha_level]


# Figure 6C. Absolute preference level in each PSD state at the primary
# observable. The dashed line at 0.5 marks the point of no preference.
# Fix 2: reverse = TRUE so that the primary score sits on top, as in 6D and 6E.

panel_6c <- state_preferences[
  observable == "tau2_ms" &
    metric == "preference_fraction" &
    outcome == "productive_capture"
]

panel_6c <- add_dataset_label(panel_6c)
panel_6c <- add_analysis_label(panel_6c, reverse = TRUE)

panel_6c[, state_label := factor(state_labels[psd_state],
                                 levels = unname(state_labels))]


# Figure 6D. The main result: WT-like minus ME-like shift in preference
# fraction, across all three observables and both outcomes.

panel_6d <- wt_me_shifts[metric == "shift_preference_fraction"]

panel_6d <- add_dataset_label(panel_6d)
panel_6d <- add_analysis_label(panel_6d, reverse = TRUE)
panel_6d <- add_outcome_label(panel_6d)
panel_6d <- add_observable_label(panel_6d)

panel_6d[, stars := significance_stars(q_bh_global)]
panel_6d[, significant := q_bh_global < alpha_level]


# Figure 6E. Mapping contrasts at the primary observable. The
# primary-minus-shuffled row is the anti-tautology test.

panel_6e <- mapping_contrasts[
  metric == "shift_preference_fraction" & observable == "tau2_ms"
]

panel_6e <- add_dataset_label(panel_6e)
panel_6e <- add_outcome_label(panel_6e)

panel_6e[, contrast_label := factor(contrast_labels[contrast],
                                    levels = rev(contrast_labels[contrast_order]))]

panel_6e[, stars := significance_stars(q_bh_global)]
panel_6e[, significant := q_bh_global < alpha_level]


# ------------------------------------------------------------------------------
# 8. Data for the Supplementary Figure S9 panels
# ------------------------------------------------------------------------------

# S9A and S9B. Relative preference for localized transcripts under each PSD
# state, at physical docking and at productive capture respectively. All three
# observables and all three analysis arms are shown.

build_state_panel <- function(target_outcome) {

  table <- state_preferences[
    metric == "preference_fraction" & outcome == target_outcome
  ]

  table <- add_dataset_label(table)
  table <- add_analysis_label(table, reverse = TRUE)
  table <- add_observable_label(table)

  table[, state_label := factor(state_labels[psd_state],
                                levels = unname(state_labels))]
  table
}

panel_s9a <- build_state_panel("physical_docking")
panel_s9b <- build_state_panel("productive_capture")


# S9C. The supplementary version of the paired sequence-score differences,
# with significance stars. Kept as a separate object so that the main and
# supplementary versions can be styled independently.

panel_s9c <- copy(panel_6b)


# S9D. WT-like minus ME-like shift in the log localized-to-soma rate ratio.
# This is the same contrast as Figure 6D on a multiplicative scale, which is
# easier to read as a relative rate change.

panel_s9d <- wt_me_shifts[metric == "shift_log_rate_ratio"]

panel_s9d <- add_dataset_label(panel_s9d)
panel_s9d <- add_analysis_label(panel_s9d, reverse = TRUE)
panel_s9d <- add_outcome_label(panel_s9d)
panel_s9d <- add_observable_label(panel_s9d)

panel_s9d[, stars := significance_stars(q_bh_global)]
panel_s9d[, significant := q_bh_global < alpha_level]

# A readable secondary column: the shift expressed as a percentage change in
# the localized-to-soma rate ratio.
panel_s9d[, rate_ratio_change_percent := 100 * (exp(estimate) - 1)]


# S9E. PSD-coordinate scaling sweep. This analysis varies the WT-minus-ME
# separation from zero to twice the reference scaling and is not part of the
# standard 18-condition grid, so its table is optional. If it is absent the
# panel is skipped and a message is printed.

scaling_file <- find_csv("panel_S9E_psd_scaling_sweep.csv", optional = TRUE)

has_scaling <- !is.na(scaling_file)

if (has_scaling) {

  panel_s9e <- fread(scaling_file)

  if (!"n_seeds" %in% names(panel_s9e) && "n" %in% names(panel_s9e)) {
    panel_s9e[, n_seeds := n]
  }

  panel_s9e <- add_dataset_label(panel_s9e)

  panel_s9e[, mapping_label := factor(
    fifelse(mode == "shuffled", "Shuffled mapping", "Authentic mapping"),
    levels = c("Authentic mapping", "Shuffled mapping")
  )]

} else {

  message(
    "Optional table panel_S9E_psd_scaling_sweep.csv not found. ",
    "Panel S9E is skipped."
  )
}


# ------------------------------------------------------------------------------
# 9. Plotting theme
# ------------------------------------------------------------------------------
# Deliberately close to theme_msc in the Figure 5 script, so that Figure 5,
# Figure 6 and Figure S9 look like one manuscript.
#
# Fix 1: strip.text.y is rotated and given explicit margins, and the right plot
# margin is widened, so that vertical facet labels are never clipped.

theme_msc <- function(base_size = 9) {

  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.40),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = element_text(size = base_size + 0.2),
      strip.background = element_blank(),
      strip.text.x = element_text(size = base_size, face = "plain",
                                  margin = margin(b = 2.5)),
      strip.text.y = element_text(size = base_size, face = "plain",
                                  angle = -90,
                                  margin = margin(l = 2.5, r = 2.5)),
      strip.clip = "off",
      panel.spacing = grid::unit(3.2, "mm"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      plot.margin = margin(5, 11, 5, 5)
    )
}

dodge_width <- 0.62

fill_scale <- scale_fill_manual(
  values = c(`TRUE` = "black", `FALSE` = "white"),
  guide = "none"
)


# ------------------------------------------------------------------------------
# 10. Main Figure 6 panels
# ------------------------------------------------------------------------------

panel_B <- ggplot(
  panel_6b,
  aes(x = estimate, y = score_label, colour = dataset_label,
      shape = dataset_label, fill = significant)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0,
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    size = 2.6,
    stroke = 0.80
  ) +
  scale_colour_manual(values = dataset_colours) +
  scale_shape_manual(values = dataset_shapes) +
  fill_scale +
  labs(x = "Paired localized \u2212 soma score difference (z units)", y = NULL) +
  theme_msc()


panel_C <- ggplot(
  panel_6c,
  aes(x = estimate, y = analysis_label, colour = dataset_label,
      shape = dataset_label)
) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0,
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    size = 2.6,
    stroke = 0.80,
    fill = "white"
  ) +
  facet_wrap(~ state_label, nrow = 1) +
  scale_colour_manual(values = dataset_colours) +
  scale_shape_manual(values = dataset_shapes) +
  scale_x_continuous(limits = preference_limits, breaks = preference_breaks) +
  labs(x = "Localized preference fraction, productive capture", y = NULL) +
  theme_msc()


panel_D <- ggplot(
  panel_6d,
  aes(x = estimate, y = analysis_label, colour = dataset_label,
      shape = dataset_label, fill = significant)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0,
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    size = 2.6,
    stroke = 0.80
  ) +
  facet_grid(outcome_label ~ observable_label) +
  scale_colour_manual(values = dataset_colours) +
  scale_shape_manual(values = dataset_shapes) +
  fill_scale +
  labs(x = "WT-like \u2212 ME-like shift in localized preference fraction",
       y = NULL) +
  theme_msc()


panel_E <- ggplot(
  panel_6e,
  aes(x = estimate, y = contrast_label, colour = dataset_label,
      shape = dataset_label, fill = significant)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0,
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    size = 2.6,
    stroke = 0.80
  ) +
  facet_wrap(~ outcome_label, nrow = 1) +
  scale_colour_manual(values = dataset_colours) +
  scale_shape_manual(values = dataset_shapes) +
  fill_scale +
  labs(x = "Contrast in WT-like \u2212 ME-like shift", y = NULL) +
  theme_msc()


# ------------------------------------------------------------------------------
# 11. Supplementary Figure S9 panels
# ------------------------------------------------------------------------------
# S9A and S9B share one builder because they differ only in the outcome stage.
# The vertical facet now uses the short dataset label (fix 1) and both panels
# share the common preference range (fix 3).

build_state_plot <- function(table, axis_label) {

  ggplot(
    table,
    aes(x = estimate, y = analysis_label, colour = state_label,
        shape = state_label)
  ) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.50) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      position = position_dodge(width = dodge_width),
      width = 0,
      linewidth = 0.70
    ) +
    geom_point(
      position = position_dodge(width = dodge_width),
      size = 2.6,
      stroke = 0.80,
      fill = "white"
    ) +
    facet_grid(dataset_strip ~ observable_label) +
    scale_colour_manual(values = state_colours) +
    scale_shape_manual(values = state_shapes) +
    scale_x_continuous(limits = preference_limits, breaks = preference_breaks) +
    labs(x = axis_label, y = NULL) +
    theme_msc()
}

panel_S9A <- build_state_plot(
  panel_s9a,
  "Localized preference fraction, physical docking"
)

panel_S9B <- build_state_plot(
  panel_s9b,
  "Localized preference fraction, productive capture"
)


panel_S9C <- ggplot(
  panel_s9c,
  aes(x = estimate, y = score_label, colour = dataset_label,
      shape = dataset_label, fill = significant)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0,
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    size = 2.6,
    stroke = 0.80
  ) +
  geom_text(
    aes(x = ci_high, label = stars),
    position = position_dodge(width = dodge_width),
    hjust = -0.35,
    size = 2.6,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = dataset_colours) +
  scale_shape_manual(values = dataset_shapes) +
  fill_scale +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = "Paired localized \u2212 soma score difference (z units)", y = NULL) +
  theme_msc()


panel_S9D <- ggplot(
  panel_s9d,
  aes(x = estimate, y = analysis_label, colour = dataset_label,
      shape = dataset_label, fill = significant)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0,
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    size = 2.6,
    stroke = 0.80
  ) +
  geom_text(
    aes(x = ci_high, label = stars),
    position = position_dodge(width = dodge_width),
    hjust = -0.35,
    size = 2.4,
    show.legend = FALSE
  ) +
  facet_grid(outcome_label ~ observable_label) +
  scale_colour_manual(values = dataset_colours) +
  scale_shape_manual(values = dataset_shapes) +
  fill_scale +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = "WT-like \u2212 ME-like shift in log localized/soma rate ratio",
       y = NULL) +
  theme_msc()


if (has_scaling) {

  panel_S9E <- ggplot(
    panel_s9e,
    aes(x = scaling, y = estimate, colour = dataset_label,
        linetype = mapping_label, fill = dataset_label)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.50) +
    geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.50) +
    geom_ribbon(
      aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 0.70) +
    scale_colour_manual(values = dataset_colours) +
    scale_fill_manual(values = dataset_colours, guide = "none") +
    labs(
      x = "PSD separation relative to the reference scaling",
      y = "WT-like \u2212 ME-like shift, productive capture"
    ) +
    theme_msc()
}


# ------------------------------------------------------------------------------
# 12. Export helper: every panel in PDF, PNG, SVG and TIFF
# ------------------------------------------------------------------------------
# A failure in one format must not abort the run. Every attempt is recorded in
# an export log written next to the figures.

export_log <- data.table(
  panel = character(),
  format = character(),
  status = character(),
  message = character()
)


save_one <- function(plot_object, file_stem, format, width_mm, height_mm) {

  filename <- file.path(OUTPUT_DIR, format, paste0(file_stem, ".", format))

  result <- tryCatch({

    device <- switch(
      format,
      pdf  = grDevices::pdf,
      png  = if (have_ragg) ragg::agg_png else grDevices::png,
      svg  = if (have_svglite) svglite::svglite else grDevices::svg,
      tiff = if (have_ragg) ragg::agg_tiff else grDevices::tiff
    )

    save_arguments <- list(
      filename = filename,
      plot = plot_object,
      device = device,
      width = width_mm,
      height = height_mm,
      units = "mm",
      bg = "white",
      limitsize = FALSE
    )

    if (format %in% c("png", "tiff")) save_arguments$dpi <- DPI
    if (format == "tiff") save_arguments$compression <- "lzw"

    do.call(ggsave, save_arguments)

    if (!file.exists(filename) || file.info(filename)$size <= 0) {
      stop("No non-empty output file was written.")
    }

    list(ok = TRUE, message = "")

  }, error = function(error) {

    if (file.exists(filename)) try(file.remove(filename), silent = TRUE)
    list(ok = FALSE, message = conditionMessage(error))
  })

  export_log <<- rbindlist(list(
    export_log,
    data.table(
      panel = file_stem,
      format = format,
      status = if (result$ok) "ok" else "failed",
      message = result$message
    )
  ))

  invisible(result$ok)
}


save_panel <- function(plot_object, file_stem, width_mm, height_mm) {
  for (format in c("pdf", "png", "svg", "tiff")) {
    save_one(plot_object, file_stem, format, width_mm, height_mm)
  }
}


# Main figure.

save_panel(panel_B, "Figure6B_sequence_effects",  width_mm = 120, height_mm = 110)
save_panel(panel_C, "Figure6C_state_preferences", width_mm = 155, height_mm = 75)
save_panel(panel_D, "Figure6D_wt_me_shifts",      width_mm = 195, height_mm = 110)
save_panel(panel_E, "Figure6E_mapping_contrasts", width_mm = 165, height_mm = 80)

# Supplementary figure. Fix 1: S9A and S9B are wider so that the vertical
# dataset strips have room.

save_panel(panel_S9A, "FigureS9A_preference_docking",    width_mm = 205, height_mm = 130)
save_panel(panel_S9B, "FigureS9B_preference_productive", width_mm = 205, height_mm = 130)
save_panel(panel_S9C, "FigureS9C_sequence_effects",      width_mm = 130, height_mm = 115)
save_panel(panel_S9D, "FigureS9D_log_rate_ratio_shifts", width_mm = 195, height_mm = 115)

if (has_scaling) {
  save_panel(panel_S9E, "FigureS9E_psd_scaling_sweep", width_mm = 150, height_mm = 95)
}


# ------------------------------------------------------------------------------
# 13. Canonical tables and session information
# ------------------------------------------------------------------------------
# The exact plotted data is written out so that every number in both figures can
# be traced without re-running the notebooks.

table_dir <- file.path(OUTPUT_DIR, "tables")

fwrite(panel_6b,  file.path(table_dir, "figure6B_plotted_data.csv"))
fwrite(panel_6c,  file.path(table_dir, "figure6C_plotted_data.csv"))
fwrite(panel_6d,  file.path(table_dir, "figure6D_plotted_data.csv"))
fwrite(panel_6e,  file.path(table_dir, "figure6E_plotted_data.csv"))

fwrite(panel_s9a, file.path(table_dir, "figureS9A_plotted_data.csv"))
fwrite(panel_s9b, file.path(table_dir, "figureS9B_plotted_data.csv"))
fwrite(panel_s9c, file.path(table_dir, "figureS9C_plotted_data.csv"))
fwrite(panel_s9d, file.path(table_dir, "figureS9D_plotted_data.csv"))

if (has_scaling) {
  fwrite(panel_s9e, file.path(table_dir, "figureS9E_plotted_data.csv"))
}

fwrite(export_log, file.path(table_dir, "figure6_S9_export_log.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "figure6_S9_session_info.txt")
)


# ------------------------------------------------------------------------------
# 14. Layout and range diagnostics
# ------------------------------------------------------------------------------
# Two checks that correspond to the version 6 fixes, so that a regression would
# be visible in the console rather than only in the rendered figure.

preference_values <- rbindlist(list(
  panel_6c[, .(panel = "6C", value = c(ci_low, ci_high))],
  panel_s9a[, .(panel = "S9A", value = c(ci_low, ci_high))],
  panel_s9b[, .(panel = "S9B", value = c(ci_low, ci_high))]
))

outside_range <- preference_values[
  value < preference_limits[1] | value > preference_limits[2]
]

if (nrow(outside_range) > 0) {
  cat("\nWARNING: some preference intervals fall outside the fixed x-range ",
      preference_limits[1], " to ", preference_limits[2], ".\n", sep = "")
  print(outside_range[, .(min = min(value), max = max(value)), by = panel])
} else {
  cat("\nAll preference intervals fit inside the shared x-range.\n")
}

analysis_top_row <- levels(panel_6c$analysis_label)[nlevels(panel_6c$analysis_label)]

cat("Top row of the preference and shift panels: ", analysis_top_row, "\n", sep = "")

if (analysis_top_row != analysis_labels[["primary_g_rg4"]]) {
  cat("WARNING: the primary score is not on the top row.\n")
}


# ------------------------------------------------------------------------------
# 15. Console summary
# ------------------------------------------------------------------------------
# A short readable report so that problems are visible without opening files.

cat("\n============ MSC Figure 6 and Figure S9 summary ============\n")

cat("\nPanel 6B / S9C, primary score by dataset:\n")
print(panel_6b[
  score == "primary_g_rg4",
  .(dataset_label,
    estimate = round(estimate, 4),
    ci_low = round(ci_low, 4),
    ci_high = round(ci_high, 4),
    q = signif(q_bh_empirical, 3),
    stars)
])

cat("\nPanel 6D, tau2 productive capture:\n")
print(panel_6d[
  observable == "tau2_ms" & outcome == "productive_capture",
  .(dataset_label, analysis_label,
    estimate = round(estimate, 4),
    ci_low = round(ci_low, 4),
    ci_high = round(ci_high, 4),
    q = signif(q_bh_global, 3),
    stars)
][order(dataset_label, analysis_label)])

cat("\nPanel S9D, tau2 productive capture as a relative rate change:\n")
print(panel_s9d[
  observable == "tau2_ms" & outcome == "productive_capture",
  .(dataset_label, analysis_label,
    log_shift = round(estimate, 4),
    percent_change = round(rate_ratio_change_percent, 1),
    q = signif(q_bh_global, 3),
    stars)
][order(dataset_label, analysis_label)])

cat("\nPanel 6E, anti-tautology contrast at tau2:\n")
print(panel_6e[
  contrast == "primary_minus_shuffled",
  .(dataset_label, outcome_label,
    estimate = round(estimate, 4),
    q = signif(q_bh_global, 3),
    stars)
][order(dataset_label, outcome_label)])

cat("\nExport results:\n")
print(export_log[, .N, by = .(format, status)])

failed_exports <- export_log[status == "failed"]

if (nrow(failed_exports) > 0) {
  cat("\nSome exports failed:\n")
  print(failed_exports)
} else {
  cat("\nAll panels exported in all four formats.\n")
}

if (!has_scaling) {
  cat("\nNote: panel S9E was skipped because the PSD scaling sweep table\n")
  cat("      panel_S9E_psd_scaling_sweep.csv was not found.\n")
}

cat("\nOutput folder: ", OUTPUT_DIR, "\n", sep = "")
cat("============================================================\n")
