# ==============================================================================
# MSC PREPRINT - FIGURE 3
# Null-model analysis: organized material identity and feedback support distinct
# components of synaptic specificity.
#
# Companion to MSC_Figure2_reproducible.ipynb. Figure 3 needs NO separate
# simulation: msc_feedback_nullmodel_metrics.csv and
# msc_feedback_paired_effects.csv are written by the same Pass 1 run that
# produces Figure 2. Make sure the notebook stages the per-seed metrics file
# into 05_r_input as well (see note at the bottom of this header).
#
# CHANGES RELATIVE TO THE EARLIER figure3.R
#   1. No setwd() and no absolute desktop path. The input folder is the first
#      command-line argument.
#   2. Four output formats per panel: pdf, png, svg, tiff. Each is attempted
#      independently with device fallbacks, so a missing cairo backend can no
#      longer prevent the other formats from being written.
#   3. New distribution panel for strong-granule docking enrichment, which
#      previously appeared only in the paired-effects grid.
#   4. New crosstalk panels. Full MSC minus activity-only gives -0.0475 for
#      consolidation crosstalk and -0.0825 for material crosstalk; these are the
#      comparisons where material selectivity is unambiguous, and they were not
#      plotted anywhere before.
#   5. Fixed axis limits are checked against the data. If a value or CI falls
#      outside, the limits expand and a warning is logged instead of silently
#      clipping the point.
#   6. Guard log, session info and per-file export log are written for the
#      reproducibility archive.
#
# STATISTICAL READING OF THE 400-SEED DATASET (for the legend and Results text)
#   Material organization dominates cargo composition:
#     productive-capture enrichment, Full MSC minus shuffled  +0.410
#     productive-capture enrichment, Full MSC minus random    +0.598
#     strong-docking enrichment,     Full MSC minus shuffled  +0.365
#   Feedback dominates the persistent state:
#     material-memory specificity,   Full MSC minus capture-only +0.323
#   Activity alone reproduces consolidation magnitude:
#     consolidation specificity,     Full MSC minus activity-only +0.019,
#     95% CI -0.005 to +0.044, i.e. NOT distinguishable from zero.
#   But activity alone does not reproduce selectivity or crosstalk suppression:
#     strong-docking enrichment,     Full MSC minus activity-only +0.381
#     material crosstalk rate,       Full MSC minus activity-only -0.083
#   Dimensionality is irrelevant:
#     projected-1D differs from Full MSC on none of the enrichment metrics.
#
#   These four statements are the honest summary of Figure 3 and should replace
#   any wording implying that Full MSC produces more consolidation than
#   activity-dependent capture alone.
#
# INPUT FILES
#   msc_feedback_nullmodel_metrics.csv   per-seed, per-variant metrics (required)
#   msc_feedback_paired_effects.csv      Full MSC minus each control (required)
#   MANIFEST_figure2.json                expected seed count (optional)
#
# USAGE
#   Rscript figure3_final.R /path/to/05_r_input
#   Rscript figure3_final.R /path/to/05_r_input strict
#
# NOTE FOR THE NOTEBOOK
#   Add this line to the staging cell of MSC_Figure2_reproducible.ipynb:
#       (D1, "msc_feedback_nullmodel_metrics.csv"),
#   Without it, 05_r_input lacks the per-seed file this script needs.
#
# Ubuntu-compatible. No absolute paths, no network access at plot time.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "tibble",
                       "scales", "patchwork")

optional_packages <- c("ragg", "svglite", "jsonlite")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE,
          FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop(paste0(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))"
  ))
}

for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    try(install.packages(pkg, repos = "https://cloud.r-project.org"),
        silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(scales)
  library(patchwork)
})

have_ragg     <- requireNamespace("ragg", quietly = TRUE)
have_svglite  <- requireNamespace("svglite", quietly = TRUE)
have_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)

has_cairo <- isTRUE(unname(capabilities("cairo")))


# ------------------------------------------------------------------------------
# 1. Paths and run mode
# ------------------------------------------------------------------------------

cli_args <- commandArgs(trailingOnly = TRUE)

INPUT_DIR <- if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  cli_args[1]
} else {
  getwd()
}

GUARD_MODE <- if (length(cli_args) >= 2 && tolower(cli_args[2]) == "strict") {
  "strict"
} else {
  "warn"
}

output_dir <- file.path(INPUT_DIR, "R_output_figure3")
for (sub in c("pdf", "png", "svg", "tiff", "plotted_data")) {
  dir.create(file.path(output_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

metrics_path <- file.path(INPUT_DIR, "msc_feedback_nullmodel_metrics.csv")
paired_path  <- file.path(INPUT_DIR, "msc_feedback_paired_effects.csv")
manifest_path <- file.path(INPUT_DIR, "MANIFEST_figure2.json")

if (!file.exists(metrics_path)) {
  stop("Per-seed metrics file not found:\n", metrics_path,
       "\nAdd msc_feedback_nullmodel_metrics.csv to the notebook staging cell.")
}
if (!file.exists(paired_path)) {
  stop("Paired-effects file not found:\n", paired_path)
}

# Expected number of parameter-generating seeds. Read from the manifest when
# available so the check follows the run that produced the data.
EXPECTED_SEEDS <- 400
if (have_jsonlite && file.exists(manifest_path)) {
  mf <- jsonlite::fromJSON(manifest_path)
  if (!is.null(mf$seeds_run)) EXPECTED_SEEDS <- as.integer(mf$seeds_run)
}

guard_log <- character(0)
log_guard <- function(label, ok, detail = "") {
  guard_log[[length(guard_log) + 1]] <<-
    sprintf("%-40s %-5s %s", label, if (ok) "PASS" else "FAIL", detail)
  if (!ok) {
    msg <- paste0(label, " failed. ", detail)
    if (GUARD_MODE == "strict") stop(msg) else warning(msg, call. = FALSE)
  }
  invisible(ok)
}

cat("\nGRAPHICS BACKEND\n")
cat("  capabilities(\"cairo\"): ", has_cairo, "\n", sep = "")
cat("  ragg installed:        ", have_ragg, "\n", sep = "")
cat("  svglite installed:     ", have_svglite, "\n", sep = "")
cat("\nFIGURE 3 PATHS\n")
cat("  Input folder:   ", normalizePath(INPUT_DIR), "\n", sep = "")
cat("  Guard mode:     ", GUARD_MODE, "\n", sep = "")
cat("  Expected seeds: ", EXPECTED_SEEDS, "\n", sep = "")
cat("  Output:         ", normalizePath(output_dir, mustWork = FALSE), "\n\n", sep = "")


# ------------------------------------------------------------------------------
# 2. Import
# ------------------------------------------------------------------------------

metrics <- read_csv(metrics_path, show_col_types = FALSE)
paired  <- read_csv(paired_path,  show_col_types = FALSE)


# ------------------------------------------------------------------------------
# 3. Validate columns
#
# crosstalk_rate and material_crosstalk_rate are new requirements relative to
# the earlier script; they drive the crosstalk panels.
# ------------------------------------------------------------------------------

required_metrics <- c("seed", "variant",
                      "continuous_specificity", "material_memory_specificity",
                      "active_productive_enrichment", "active_strong_enrichment",
                      "crosstalk_rate", "material_crosstalk_rate")

required_paired <- c("metric", "comparison", "n_pairs",
                     "mean_difference", "sem_difference")

missing_metrics <- setdiff(required_metrics, names(metrics))
missing_paired  <- setdiff(required_paired,  names(paired))

if (length(missing_metrics) > 0) {
  stop("Missing columns in metrics file: ",
       paste(missing_metrics, collapse = ", "))
}
if (length(missing_paired) > 0) {
  stop("Missing columns in paired-effects file: ",
       paste(missing_paired, collapse = ", "))
}


# ------------------------------------------------------------------------------
# 4. Validate seed counts
#
# Enrichment metrics are undefined for seeds in which a variant produced no
# docking events, so n_pairs can legitimately fall slightly below the seed
# count. The exact numbers are logged rather than hidden.
# ------------------------------------------------------------------------------

seed_check <- metrics %>% count(variant, name = "n_seeds")
print(seed_check)

log_guard(
  "metrics rows per variant",
  all(seed_check$n_seeds == EXPECTED_SEEDS),
  sprintf("observed %s, expected %d",
          paste(unique(seed_check$n_seeds), collapse = "/"), EXPECTED_SEEDS)
)

pair_counts <- paired %>%
  group_by(metric) %>%
  summarise(min_pairs = min(n_pairs), max_pairs = max(n_pairs), .groups = "drop")
print(pair_counts)

log_guard(
  "paired n_pairs within 2% of seeds",
  all(pair_counts$min_pairs >= 0.98 * EXPECTED_SEEDS),
  sprintf("minimum n_pairs = %d", min(pair_counts$min_pairs))
)


# ------------------------------------------------------------------------------
# 5. Model labels and colours  (unchanged from the earlier script)
# ------------------------------------------------------------------------------

variant_order <- c("full_msc_feedback", "msc_capture_only",
                   "projected1d_feedback", "shuffled_material_feedback",
                   "random_material_feedback", "tag_only_feedback")

variant_labels <- c(
  full_msc_feedback          = "Full MSC\nfeedback",
  msc_capture_only           = "MSC capture\nonly",
  projected1d_feedback       = "Projected 1D\nfeedback",
  shuffled_material_feedback = "Shuffled\nmaterial",
  random_material_feedback   = "Random\nmaterial",
  tag_only_feedback          = "Activity-only\nfeedback"
)

variant_colours <- c(
  full_msc_feedback          = "#3A923A",
  msc_capture_only           = "#3C88C9",
  projected1d_feedback       = "#8752B5",
  shuffled_material_feedback = "#F28E2B",
  random_material_feedback   = "#A7A9AC",
  tag_only_feedback          = "#D94C4C"
)

comparison_order <- c("msc_capture_only", "projected1d_feedback",
                      "shuffled_material_feedback",
                      "random_material_feedback", "tag_only_feedback")

comparison_labels <- c(
  msc_capture_only           = "Capture only",
  projected1d_feedback       = "Projected 1D",
  shuffled_material_feedback = "Shuffled material",
  random_material_feedback   = "Random material",
  tag_only_feedback          = "Activity only"
)

metrics <- metrics %>% mutate(variant = factor(variant, levels = variant_order))


# ------------------------------------------------------------------------------
# 6. Distribution-panel theme  (unchanged)
# ------------------------------------------------------------------------------

theme_nature_distribution <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.45),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.6, "mm"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1,
                                 margin = margin(t = 4)),
      axis.title = element_text(colour = "black", size = base_size + 0.5),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.title = element_text(face = "bold", colour = "black",
                                size = base_size + 1.4, lineheight = 1.05,
                                hjust = 0, margin = margin(t = 2, b = 8)),
      panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 12, r = 16, b = 12, l = 14)
    )
}


# ------------------------------------------------------------------------------
# 7. Paired-effect theme  (unchanged)
# ------------------------------------------------------------------------------

theme_nature_effect <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.45),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.4, "mm"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(margin = margin(r = 2)),
      axis.title = element_text(colour = "black", size = base_size + 0.3),
      plot.title = element_text(face = "bold", colour = "black",
                                size = base_size + 0.9, lineheight = 1.0,
                                hjust = 0.5, margin = margin(t = 2, b = 7)),
      panel.grid.major.y = element_line(colour = "#EAEAEA", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 8, r = 10, b = 8, l = 10)
    )
}

theme_nature_effect_legend <- function(base_size = 8.5) {
  theme_nature_effect(base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(4.6, "mm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(t = -2)
    )
}


# ------------------------------------------------------------------------------
# 8. Helpers
# ------------------------------------------------------------------------------

summarise_metric <- function(data, variable) {
  data %>%
    filter(is.finite(.data[[variable]])) %>%
    group_by(variant) %>%
    summarise(
      n = n(),
      mean = mean(.data[[variable]], na.rm = TRUE),
      sd = sd(.data[[variable]], na.rm = TRUE),
      sem = sd / sqrt(n),
      lower = mean - 1.96 * sem,
      upper = mean + 1.96 * sem,
      .groups = "drop"
    )
}

# Fixed limits keep the published look, but silent clipping is dangerous when a
# dataset is regenerated. safe_limits() expands the range if needed and logs it.
safe_limits <- function(fixed, values, label) {
  values <- values[is.finite(values)]
  if (length(values) == 0 || is.null(fixed)) return(fixed)
  lo <- min(values); hi <- max(values)
  if (lo >= fixed[1] && hi <= fixed[2]) {
    log_guard(paste0("axis limits: ", label), TRUE,
              sprintf("data [%.4g, %.4g] inside [%.4g, %.4g]",
                      lo, hi, fixed[1], fixed[2]))
    return(fixed)
  }
  pad <- 0.05 * (hi - lo)
  expanded <- c(min(fixed[1], lo - pad), max(fixed[2], hi + pad))
  log_guard(paste0("axis limits: ", label), FALSE,
            sprintf("data [%.4g, %.4g] exceeded [%.4g, %.4g]; expanded to [%.4g, %.4g]",
                    lo, hi, fixed[1], fixed[2], expanded[1], expanded[2]))
  expanded
}

# Wilson interval for the binary crosstalk rates.
wilson <- function(k, n, z = 1.959964) {
  p <- k / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / denom
  tibble(rate = p,
         lower = pmax(0, centre - half),
         upper = pmin(1, centre + half))
}


# ------------------------------------------------------------------------------
# 9. Distribution-panel function  (geometry unchanged)
# ------------------------------------------------------------------------------

distribution_panel <- function(data, variable, title, y_label,
                               reference_line, reference_label = NULL,
                               y_limits = NULL, y_breaks = waiver(),
                               limit_label = variable) {

  plot_data <- data %>% filter(is.finite(.data[[variable]]))
  summary_data <- summarise_metric(data = data, variable = variable)

  y_limits <- safe_limits(
    y_limits,
    c(plot_data[[variable]], summary_data$lower, summary_data$upper),
    limit_label
  )

  p <- ggplot(plot_data,
              aes(x = variant, y = .data[[variable]],
                  colour = variant, fill = variant)) +
    geom_hline(yintercept = reference_line, linetype = "dashed",
               linewidth = 0.45, colour = "#777777") +
    geom_violin(width = 0.78, trim = TRUE, scale = "width",
                linewidth = 0.35, alpha = 0.16) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 31),
               size = 0.75, alpha = 0.28, stroke = 0) +
    geom_errorbar(data = summary_data,
                  aes(x = variant, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, width = 0.13, linewidth = 0.65,
                  colour = "black") +
    geom_point(data = summary_data, aes(x = variant, y = mean),
               inherit.aes = FALSE, shape = 21, size = 2.7, stroke = 0.75,
               fill = "white", colour = "black") +
    scale_colour_manual(values = variant_colours, drop = FALSE) +
    scale_fill_manual(values = variant_colours, drop = FALSE) +
    scale_x_discrete(labels = variant_labels, drop = FALSE) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks,
                       expand = expansion(mult = c(0.04, 0.08))) +
    labs(title = title, x = NULL, y = y_label) +
    theme_nature_distribution()

  if (!is.null(reference_label)) {
    p <- p + annotate(geom = "text", x = 0.58, y = reference_line,
                      label = reference_label, hjust = 0, vjust = -0.65,
                      size = 2.45, colour = "#666666")
  }
  p
}


# ------------------------------------------------------------------------------
# 10. Distribution panels
# ------------------------------------------------------------------------------

panel_specificity <- distribution_panel(
  data = metrics, variable = "continuous_specificity",
  title = paste0("Organized material identity and feedback\n",
                 "support distinct components of synaptic specificity"),
  y_label = paste0("Consolidation-state specificity\n",
                   "active spine \u2212 maximum control spine"),
  reference_line = 0, reference_label = "no specificity",
  y_limits = c(-0.5, 1.05), y_breaks = seq(-0.5, 1.0, by = 0.5),
  limit_label = "consolidation specificity"
)

panel_memory <- distribution_panel(
  data = metrics, variable = "material_memory_specificity",
  title = paste0("Feedback establishes local\n", "material memory"),
  y_label = paste0("Material-memory specificity\n",
                   "active spine \u2212 maximum control spine"),
  reference_line = 0, reference_label = "no specificity",
  y_limits = c(-1.0, 1.05), y_breaks = seq(-1.0, 1.0, by = 0.5),
  limit_label = "material-memory specificity"
)

panel_productive <- distribution_panel(
  data = metrics, variable = "active_productive_enrichment",
  title = paste0("Material structure supports\n", "productive capture"),
  y_label = paste0("Productive-capture enrichment\n", "active spine / pool"),
  reference_line = 1, reference_label = "no enrichment",
  y_limits = c(0, 2.5), y_breaks = seq(0, 2.5, by = 0.5),
  limit_label = "productive-capture enrichment"
)

# NEW: strong-granule docking enrichment as a distribution. Previously this
# metric appeared only in the paired-effects grid, although it carries the
# largest Full-MSC-minus-activity-only difference (+0.381).
panel_strong <- distribution_panel(
  data = metrics, variable = "active_strong_enrichment",
  title = paste0("Material identity biases\n", "which granules dock"),
  y_label = paste0("Strong-granule docking enrichment\n", "active spine / pool"),
  reference_line = 1, reference_label = "no enrichment",
  y_limits = c(0, 3.0), y_breaks = seq(0, 3.0, by = 0.5),
  limit_label = "strong-docking enrichment"
)


# ------------------------------------------------------------------------------
# 11. NEW: crosstalk panel
#
# crosstalk_rate and material_crosstalk_rate are binary per seed, so a violin is
# inappropriate. Mean rates with Wilson intervals are shown instead. This is the
# comparison in which activity-only feedback fails most clearly: it crosses the
# consolidation threshold at a control synapse in 4.75 percentage points more
# seeds, and the material threshold in 8.25 more.
# ------------------------------------------------------------------------------

crosstalk_summary <- metrics %>%
  select(variant, crosstalk_rate, material_crosstalk_rate) %>%
  pivot_longer(-variant, names_to = "metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  group_by(variant, metric) %>%
  summarise(k = sum(value > 0), n = n(), .groups = "drop") %>%
  bind_cols(., wilson(.$k, .$n)) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("crosstalk_rate", "material_crosstalk_rate"),
      labels = c("Consolidation crosstalk", "Material-memory crosstalk")
    )
  )

panel_crosstalk <- ggplot(
  crosstalk_summary,
  aes(x = variant, y = rate, fill = metric, colour = metric)
) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           alpha = 0.85, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.72),
                width = 0.16, linewidth = 0.55, colour = "black") +
  scale_fill_manual(values = c("Consolidation crosstalk" = "#2E7D32",
                               "Material-memory crosstalk" = "#7B1FA2")) +
  scale_colour_manual(values = c("Consolidation crosstalk" = "#2E7D32",
                                 "Material-memory crosstalk" = "#7B1FA2")) +
  scale_x_discrete(labels = variant_labels, drop = FALSE) +
  scale_y_continuous(labels = label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = paste0("Material selectivity suppresses\n", "off-target crosstalk"),
    x = NULL,
    y = paste0("Fraction of seeds with a control spine\n",
               "above the readout threshold")
  ) +
  theme_nature_distribution() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8, colour = "black"),
    legend.key.height = unit(3.4, "mm"),
    legend.key.width = unit(4.6, "mm"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(b = -2)
  )


# ------------------------------------------------------------------------------
# 12. Prepare paired-effect data
# ------------------------------------------------------------------------------

effect_metrics <- c("continuous_specificity", "material_memory_specificity",
                    "active_productive_enrichment", "active_strong_enrichment",
                    "crosstalk_rate", "material_crosstalk_rate")

paired_effects <- paired %>%
  filter(metric %in% effect_metrics, comparison %in% comparison_order) %>%
  mutate(
    comparison = factor(comparison, levels = rev(comparison_order)),
    ci_lower = mean_difference - 1.96 * sem_difference,
    ci_upper = mean_difference + 1.96 * sem_difference,
    # Flagging which intervals exclude zero makes the honest reading explicit
    # in the exported table: consolidation specificity vs activity-only does
    # not, whereas strong-docking enrichment vs activity-only clearly does.
    excludes_zero = (ci_lower > 0) | (ci_upper < 0)
  )

paired_check <- paired_effects %>% count(metric, name = "n_comparisons")
print(paired_check)

log_guard("five controls per paired metric",
          all(paired_check$n_comparisons == 5),
          sprintf("observed %s",
                  paste(unique(paired_check$n_comparisons), collapse = "/")))


# ------------------------------------------------------------------------------
# 13. Paired-effect plotting function  (geometry unchanged)
# ------------------------------------------------------------------------------

effect_panel <- function(data, selected_metric, title, colour,
                         x_limits = NULL, x_breaks = waiver(),
                         show_y_labels = TRUE, base_size = 8.5) {

  plot_data <- data %>% filter(metric == selected_metric)

  x_limits <- safe_limits(x_limits,
                          c(plot_data$ci_lower, plot_data$ci_upper),
                          paste0("effect ", selected_metric))

  p <- ggplot(plot_data, aes(x = mean_difference, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "#777777", linewidth = 0.45) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                  width = 0.16, linewidth = 0.7, colour = colour) +
    geom_point(shape = 21, size = 2.8, stroke = 0.5,
               fill = colour, colour = colour) +
    scale_y_discrete(labels = comparison_labels, drop = FALSE) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       expand = expansion(mult = c(0.05, 0.08))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_nature_effect(base_size = base_size)

  if (!show_y_labels) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}


# ------------------------------------------------------------------------------
# 14. Individual paired-effect panels for the combined grid
# ------------------------------------------------------------------------------

effect_specificity <- effect_panel(
  paired_effects, "continuous_specificity", "Consolidation\nspecificity",
  "#2E7D32", x_limits = c(-0.04, 0.22), x_breaks = seq(0, 0.2, by = 0.1),
  show_y_labels = TRUE, base_size = 8.0
)

effect_memory <- effect_panel(
  paired_effects, "material_memory_specificity", "Material\nmemory",
  "#7B1FA2", x_limits = c(-0.04, 0.38), x_breaks = seq(0, 0.3, by = 0.1),
  show_y_labels = FALSE, base_size = 8.0
)

effect_productive <- effect_panel(
  paired_effects, "active_productive_enrichment", "Productive\ncapture",
  "#008C83", x_limits = c(-0.08, 0.68), x_breaks = seq(0, 0.6, by = 0.2),
  show_y_labels = TRUE, base_size = 8.0
)

effect_strong <- effect_panel(
  paired_effects, "active_strong_enrichment", "Strong\ndocking",
  "#F39C12", x_limits = c(-0.10, 0.46), x_breaks = seq(0, 0.4, by = 0.2),
  show_y_labels = FALSE, base_size = 8.0
)

# NEW: crosstalk effects. Negative values mean Full MSC has LOWER crosstalk
# than the control, which is the direction of interest for activity-only.
effect_crosstalk <- paired_effects %>%
  filter(metric %in% c("crosstalk_rate", "material_crosstalk_rate")) %>%
  mutate(metric = factor(
    metric,
    levels = c("crosstalk_rate", "material_crosstalk_rate"),
    labels = c("Consolidation", "Material memory")
  )) %>%
  ggplot(aes(x = mean_difference, y = comparison,
             colour = metric, fill = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#777777", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                position = position_dodge(width = 0.62),
                width = 0.14, linewidth = 0.65) +
  geom_point(shape = 21, size = 2.4, stroke = 0.5,
             position = position_dodge(width = 0.62)) +
  scale_colour_manual(values = c("Consolidation" = "#2E7D32",
                                 "Material memory" = "#7B1FA2")) +
  scale_fill_manual(values = c("Consolidation" = "#2E7D32",
                               "Material memory" = "#7B1FA2")) +
  scale_y_discrete(labels = comparison_labels, drop = FALSE) +
  scale_x_continuous(breaks = pretty_breaks(n = 5),
                     labels = label_number(accuracy = 0.02),
                     expand = expansion(mult = c(0.08, 0.08))) +
  labs(title = "Crosstalk rate\n(negative = less crosstalk in Full MSC)",
       x = NULL, y = NULL) +
  theme_nature_effect_legend(base_size = 8.0)


# ------------------------------------------------------------------------------
# 15. Combined paired-effects panel
#
# Caption to add below the panel during assembly:
#   Paired effect: Full MSC minus null model; points show means and bars show
#   approximate 95% confidence intervals across matched parameter seeds.
# ------------------------------------------------------------------------------

panel_paired_effects <- (
  (effect_specificity | effect_memory) /
  (effect_productive | effect_strong)
) +
  plot_layout(widths = c(1.06, 0.94), heights = c(1, 1)) +
  plot_annotation(
    title = "Full MSC effects relative to null models",
    theme = theme(
      plot.title = element_text(family = "sans", face = "bold", size = 11,
                                colour = "black", hjust = 0,
                                margin = margin(t = 2, b = 8)),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
  )


# ------------------------------------------------------------------------------
# 16. Individual-effect versions with full labels, for separate exports
# ------------------------------------------------------------------------------

effect_specificity_individual <- effect_panel(
  paired_effects, "continuous_specificity", "Consolidation specificity",
  "#2E7D32", x_limits = c(-0.04, 0.22), x_breaks = seq(0, 0.2, by = 0.1),
  show_y_labels = TRUE, base_size = 9.0
)

effect_memory_individual <- effect_panel(
  paired_effects, "material_memory_specificity", "Material memory",
  "#7B1FA2", x_limits = c(-0.04, 0.38), x_breaks = seq(0, 0.3, by = 0.1),
  show_y_labels = TRUE, base_size = 9.0
)

effect_productive_individual <- effect_panel(
  paired_effects, "active_productive_enrichment", "Productive capture",
  "#008C83", x_limits = c(-0.08, 0.68), x_breaks = seq(0, 0.6, by = 0.2),
  show_y_labels = TRUE, base_size = 9.0
)

effect_strong_individual <- effect_panel(
  paired_effects, "active_strong_enrichment", "Strong docking",
  "#F39C12", x_limits = c(-0.10, 0.46), x_breaks = seq(0, 0.4, by = 0.2),
  show_y_labels = TRUE, base_size = 9.0
)


# ------------------------------------------------------------------------------
# 17. Combined previews
# ------------------------------------------------------------------------------

figure_3_preview <- (
  (panel_specificity | panel_memory) /
  (panel_productive | panel_paired_effects)
) +
  plot_layout(widths = c(1.02, 0.98), heights = c(1, 1.08)) &
  theme(plot.margin = margin(t = 8, r = 10, b = 8, l = 10))

# Second preview for the new panels, kept separate so the published layout is
# untouched.
figure_3_supp_preview <- (
  (panel_strong | panel_crosstalk) / effect_crosstalk
) +
  plot_layout(heights = c(1, 0.85)) &
  theme(plot.margin = margin(t = 8, r = 10, b = 8, l = 10))


# ------------------------------------------------------------------------------
# 18. Robust export helper: pdf, png, svg, tiff
#
# Every format is attempted independently with device fallbacks, and each file
# is verified to exist and be non-empty. Failures are logged, never fatal.
# ------------------------------------------------------------------------------

export_log <- data.frame(panel = character(0), format = character(0),
                         device = character(0), status = character(0),
                         bytes = numeric(0), message = character(0),
                         stringsAsFactors = FALSE)

try_devices <- function(panel_name, fmt, file_path, plot_object,
                        width_mm, height_mm, attempts) {
  for (att in attempts) {
    res <- tryCatch({
      att$fun(file_path, plot_object, width_mm, height_mm)
      if (!file.exists(file_path)) stop("device reported success but no file was written")
      sz <- file.info(file_path)$size
      if (is.na(sz) || sz <= 0) stop("file written but empty (0 bytes)")
      list(ok = TRUE, device = att$label, bytes = sz, message = "")
    }, error = function(e) {
      if (file.exists(file_path)) try(file.remove(file_path), silent = TRUE)
      list(ok = FALSE, device = att$label, bytes = 0,
           message = conditionMessage(e))
    })
    export_log <<- rbind(export_log, data.frame(
      panel = panel_name, format = fmt, device = res$device,
      status = if (isTRUE(res$ok)) "OK" else "FAILED",
      bytes = res$bytes, message = res$message, stringsAsFactors = FALSE
    ))
    if (isTRUE(res$ok)) return(invisible(TRUE))
  }
  invisible(FALSE)
}

pdf_attempts <- list(
  list(label = "cairo_pdf", fun = function(f, p, w, h)
    ggsave(f, p, device = grDevices::cairo_pdf, width = w, height = h,
           units = "mm", bg = "white", limitsize = FALSE)),
  list(label = "grDevices::pdf", fun = function(f, p, w, h)
    ggsave(f, p, device = grDevices::pdf, width = w, height = h,
           units = "mm", bg = "white", limitsize = FALSE))
)

png_attempts <- list()
if (have_ragg) {
  png_attempts <- c(png_attempts, list(list(
    label = "ragg::agg_png", fun = function(f, p, w, h)
      ggsave(f, p, device = ragg::agg_png, width = w, height = h,
             units = "mm", dpi = 600, bg = "white", limitsize = FALSE))))
}
if (has_cairo) {
  png_attempts <- c(png_attempts, list(list(
    label = "png(type=cairo)", fun = function(f, p, w, h)
      ggsave(f, p, device = function(filename, ...)
        grDevices::png(filename, ..., type = "cairo"),
        width = w, height = h, units = "mm", dpi = 600,
        bg = "white", limitsize = FALSE))))
}
png_attempts <- c(png_attempts, list(list(
  label = "grDevices::png", fun = function(f, p, w, h)
    ggsave(f, p, device = grDevices::png, width = w, height = h,
           units = "mm", dpi = 600, bg = "white", limitsize = FALSE))))

svg_attempts <- list()
if (have_svglite) {
  svg_attempts <- c(svg_attempts, list(list(
    label = "svglite::svglite", fun = function(f, p, w, h)
      ggsave(f, p, device = svglite::svglite, width = w, height = h,
             units = "mm", bg = "white", limitsize = FALSE))))
}
svg_attempts <- c(svg_attempts, list(list(
  label = "grDevices::svg", fun = function(f, p, w, h)
    ggsave(f, p, device = grDevices::svg, width = w, height = h,
           units = "mm", bg = "white", limitsize = FALSE))))

tiff_attempts <- list()
if (have_ragg) {
  tiff_attempts <- c(tiff_attempts, list(list(
    label = "ragg::agg_tiff", fun = function(f, p, w, h)
      ggsave(f, p, device = ragg::agg_tiff, width = w, height = h,
             units = "mm", dpi = 600, bg = "white", limitsize = FALSE,
             compression = "lzw"))))
}
tiff_attempts <- c(tiff_attempts, list(list(
  label = "grDevices::tiff(lzw)", fun = function(f, p, w, h)
    ggsave(f, p, device = grDevices::tiff, width = w, height = h,
           units = "mm", dpi = 600, bg = "white", limitsize = FALSE,
           compression = "lzw"))))
tiff_attempts <- c(tiff_attempts, list(list(
  label = "grDevices::tiff", fun = function(f, p, w, h)
    ggsave(f, p, device = grDevices::tiff, width = w, height = h,
           units = "mm", dpi = 600, bg = "white", limitsize = FALSE))))

export_panel <- function(plot_object, file_stem, width_mm, height_mm) {
  if (is.null(plot_object)) {
    cat(sprintf("skipped  %-50s (panel not built)\n", file_stem))
    return(invisible(NULL))
  }
  ok <- c(
    pdf  = try_devices(file_stem, "pdf",
                       file.path(output_dir, "pdf",  paste0(file_stem, ".pdf")),
                       plot_object, width_mm, height_mm, pdf_attempts),
    png  = try_devices(file_stem, "png",
                       file.path(output_dir, "png",  paste0(file_stem, ".png")),
                       plot_object, width_mm, height_mm, png_attempts),
    svg  = try_devices(file_stem, "svg",
                       file.path(output_dir, "svg",  paste0(file_stem, ".svg")),
                       plot_object, width_mm, height_mm, svg_attempts),
    tiff = try_devices(file_stem, "tiff",
                       file.path(output_dir, "tiff", paste0(file_stem, ".tiff")),
                       plot_object, width_mm, height_mm, tiff_attempts)
  )
  cat(sprintf("saved    %-50s pdf:%-3s png:%-3s svg:%-3s tiff:%-3s (%g x %g mm)\n",
              file_stem,
              if (ok["pdf"])  "ok" else "NO",
              if (ok["png"])  "ok" else "NO",
              if (ok["svg"])  "ok" else "NO",
              if (ok["tiff"]) "ok" else "NO",
              width_mm, height_mm))
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# 19. Export
# ------------------------------------------------------------------------------

cat("\nEXPORTING PANELS\n")

export_panel(panel_specificity, "Figure3_specificity",        110, 95)
export_panel(panel_memory,      "Figure3_material_memory",    110, 95)
export_panel(panel_productive,  "Figure3_productive_capture", 110, 95)
export_panel(panel_strong,      "FigS3_strong_docking",       110, 95)
export_panel(panel_crosstalk,   "FigS3_crosstalk_rates",      120, 100)

export_panel(panel_paired_effects, "Figure3_paired_effects", 120, 100)

export_panel(effect_specificity_individual,
             "Figure3_effect_consolidation_specificity", 90, 65)
export_panel(effect_memory_individual,
             "Figure3_effect_material_memory",           90, 65)
export_panel(effect_productive_individual,
             "Figure3_effect_productive_capture",        90, 65)
export_panel(effect_strong_individual,
             "Figure3_effect_strong_docking",            90, 65)
export_panel(effect_crosstalk,
             "FigS3_effect_crosstalk",                  100, 78)

export_panel(figure_3_preview,      "Figure3_combined_preview_no_letters", 210, 180)
export_panel(figure_3_supp_preview, "FigS3_supplementary_preview",         220, 170)


# ------------------------------------------------------------------------------
# 20. Supporting tables
# ------------------------------------------------------------------------------

data_dir <- file.path(output_dir, "plotted_data")

write_csv(paired_effects %>% arrange(metric, comparison),
          file.path(data_dir, "Figure3_paired_effects_with_95CI.csv"))
write_csv(seed_check, file.path(data_dir, "Figure3_seed_count_check.csv"))
write_csv(crosstalk_summary, file.path(data_dir, "FigS3_crosstalk_summary.csv"))

distribution_summary <- bind_rows(
  summarise_metric(metrics, "continuous_specificity") %>%
    mutate(metric = "continuous_specificity"),
  summarise_metric(metrics, "material_memory_specificity") %>%
    mutate(metric = "material_memory_specificity"),
  summarise_metric(metrics, "active_productive_enrichment") %>%
    mutate(metric = "active_productive_enrichment"),
  summarise_metric(metrics, "active_strong_enrichment") %>%
    mutate(metric = "active_strong_enrichment")
) %>%
  select(metric, everything())

write_csv(distribution_summary,
          file.path(data_dir, "Figure3_distribution_summary.csv"))

write_csv(export_log, file.path(output_dir, "Figure3_R_export_log.csv"))
writeLines(guard_log, con = file.path(output_dir, "Figure3_R_guard_log.txt"))
writeLines(capture.output(sessionInfo()),
           con = file.path(output_dir, "Figure3_R_sessionInfo.txt"))


# ------------------------------------------------------------------------------
# 21. Interpretation printout
#
# Prints which paired comparisons exclude zero, so the Results text can be
# written from the data rather than from memory.
# ------------------------------------------------------------------------------

cat("\nPAIRED EFFECTS: WHICH 95% INTERVALS EXCLUDE ZERO\n")
interp <- paired_effects %>%
  arrange(metric, desc(abs(mean_difference))) %>%
  mutate(
    comparison = as.character(comparison),
    effect = sprintf("%+.4f [%+.4f, %+.4f]",
                     mean_difference, ci_lower, ci_upper),
    verdict = ifelse(excludes_zero, "excludes 0", "INCLUDES 0")
  ) %>%
  select(metric, comparison, n_pairs, effect, verdict)
print(interp, n = Inf, width = Inf)

cat("\nNOTE FOR THE RESULTS TEXT\n")
ns <- paired_effects %>%
  filter(!excludes_zero) %>%
  mutate(comparison = as.character(comparison))
if (nrow(ns) > 0) {
  cat("  The following comparisons are NOT distinguishable from zero and must\n")
  cat("  not be described as differences:\n")
  for (i in seq_len(nrow(ns))) {
    cat(sprintf("    %-30s vs %-28s %+.4f [%+.4f, %+.4f]\n",
                ns$metric[i], ns$comparison[i], ns$mean_difference[i],
                ns$ci_lower[i], ns$ci_upper[i]))
  }
} else {
  cat("  All plotted paired comparisons exclude zero.\n")
}


# ------------------------------------------------------------------------------
# 22. Export status
# ------------------------------------------------------------------------------

cat("\nEXPORT STATUS BY FORMAT\n")
for (fmt in c("pdf", "png", "svg", "tiff")) {
  sub <- export_log[export_log$format == fmt, , drop = FALSE]
  n_ok <- sum(sub$status == "OK")
  devs <- unique(sub$device[sub$status == "OK"])
  cat(sprintf("  %-5s %2d file(s)  device: %-22s %6.2f MB\n",
              toupper(fmt), n_ok,
              if (length(devs)) paste(devs, collapse = ", ") else "-",
              sum(sub$bytes[sub$status == "OK"]) / 1e6))
  if (n_ok == 0) {
    fails <- sub[sub$status == "FAILED", , drop = FALSE]
    cat("        NO FILES WRITTEN. Failed attempts:\n")
    if (nrow(fails) == 0) {
      cat("        (no device available for this format)\n")
    } else {
      for (i in seq_len(nrow(fails))) {
        cat(sprintf("        %-22s %s\n", fails$device[i], fails$message[i]))
      }
    }
    if (fmt %in% c("png", "tiff")) {
      cat("        Hint: install.packages(\"ragg\") provides cairo-independent devices.\n")
    }
  }
}

cat(sprintf("\nGraphics files on disk: %d\n",
            length(list.files(output_dir, recursive = TRUE,
                              pattern = "\\.(pdf|png|svg|tiff)$"))))

cat("\nGUARD LOG\n")
cat(paste(guard_log, collapse = "\n"), "\n")

cat("\nFIGURE 3 EXPORT COMPLETED\n")
cat("Output directory:", normalizePath(output_dir, mustWork = FALSE), "\n")
cat("No panel letters were generated; add them during assembly.\n")
