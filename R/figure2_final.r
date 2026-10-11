# ==============================================================================
# MSC PREPRINT - FIGURE 2 (+ SUPPLEMENTARY PANELS)
#
# Companion to MSC_Figure2_reproducible.ipynb. The notebook writes every input
# consumed here into its 05_r_input/ folder; this script turns them into the
# publication panels. Nothing is computed from scratch: R only plots.
#
# EVERY PANEL IS EXPORTED IN THREE FORMATS, INDEPENDENTLY
#   .pdf   vector, for LaTeX / journal submission
#   .png   raster, 600 dpi, for slides and quick review
#   .svg   vector with editable text, for Inkscape / Illustrator / BioRender
#
# Each format is written inside its own tryCatch with device fallbacks, so a
# missing cairo backend can no longer prevent the other formats from being
# written. Section 17 prints a per-format status table; if any format has zero
# files, the reason is printed there.
#
#   PDF fallback chain:  cairo_pdf  ->  grDevices::pdf
#   PNG fallback chain:  ragg::agg_png  ->  png(type="cairo")  ->  png()
#   SVG fallback chain:  svglite::svglite  ->  grDevices::svg
#
# MAIN PANELS
#   A   Docking probability across 12 quantile bins of material mismatch d^2
#   B   Logistic observation-model coefficients for T, M and d^2
#
# SUPPLEMENTARY PANELS
#   B2  The same observation model INCLUDING the active-synapse indicator.
#       Panel B shows only the three terms that enter the docking hazard; B2
#       exists so that the fourth fitted predictor is never silently dropped.
#   C   Observed / kernel-predicted docking probability per bin. The docking
#       hazard contains exp(-0.5 d^2 / ell_h^2), so log p_dock must be linear
#       in d^2. This panel is an implementation validation, not a result.
#   D   Panel A stratified by synapse state. The lowest-mismatch bin of the
#       pooled curve sits above the kernel expectation because the active
#       synapse carries its target vector at the origin: it contributes ~21%
#       of that bin's observations but ~72% of its docking events.
#   E   Binned d^2 -> docking curves for all six model variants, replacing the
#       earlier single-variant "tautology test" figure.
#
# INPUT FILES (all in one flat folder, produced by the notebook)
#   mismatch_binned_docking_probability.csv     Panel A          (required)
#   msc_feedback_binary_observation_model.csv   Panels B and B2  (required)
#   FigS_panelA_kernel_consistency.csv          Panel C          (optional)
#   FigS_panelA_stratified_active.csv           Panel D          (optional)
#   tautology_binned_by_variant.csv             Panel E          (optional)
#   figure2_paper_numbers.csv                   console summary  (optional)
#   MANIFEST_figure2.json                       expected values  (optional)
#
# USAGE
#   Rscript figure2_final.R /path/to/05_r_input
#   Rscript figure2_final.R /path/to/05_r_input strict
#
# Second argument controls the version guards:
#   warn   (default) mismatches warn and plotting continues
#   strict mismatches abort; use this to prove exact reproduction
#
# Ubuntu-compatible. No absolute paths, no network access at plot time.
# ==============================================================================


# ==============================================================================
# 1. PACKAGE SETUP
# ==============================================================================

required_packages <- c("ggplot2", "dplyr", "readr", "scales", "patchwork")

# ragg     better PNG device, independent of the cairo build
# svglite  SVG with text kept as text
# jsonlite read expected values from the notebook manifest
optional_packages <- c("ragg", "svglite", "jsonlite")

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

for (pkg in optional_packages) {
  if (!pkg %in% rownames(installed.packages())) {
    try(install.packages(pkg, repos = "https://cloud.r-project.org"),
        silent = TRUE)
  }
}

library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(patchwork)

have_ragg     <- requireNamespace("ragg", quietly = TRUE)
have_svglite  <- requireNamespace("svglite", quietly = TRUE)
have_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)

# capabilities() tells us up front whether the cairo backend exists. On a build
# without cairo, ggsave(device = png, type = "cairo") throws, which is the usual
# reason a run ends up with only some of the formats.
has_cairo <- isTRUE(unname(capabilities("cairo")))
has_png   <- isTRUE(unname(capabilities("png")))


# ==============================================================================
# 2. FILE PATHS
# ==============================================================================

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

OUTPUT_DIR <- file.path(INPUT_DIR, "R_output_figure2")

for (sub in c("pdf", "png", "svg", "plotted_data")) {
  dir.create(file.path(OUTPUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}

input_path <- function(name) file.path(INPUT_DIR, name)

mismatch_file    <- input_path("mismatch_binned_docking_probability.csv")
observation_file <- input_path("msc_feedback_binary_observation_model.csv")
kernel_file      <- input_path("FigS_panelA_kernel_consistency.csv")
stratified_file  <- input_path("FigS_panelA_stratified_active.csv")
tautology_file   <- input_path("tautology_binned_by_variant.csv")
numbers_file     <- input_path("figure2_paper_numbers.csv")
manifest_file    <- input_path("MANIFEST_figure2.json")

if (!file.exists(mismatch_file)) {
  stop("Panel A input file not found:\n", mismatch_file,
       "\nPass the notebook's 05_r_input folder as the first argument.")
}
if (!file.exists(observation_file)) {
  stop("Panel B input file not found:\n", observation_file)
}

cat("\nGRAPHICS BACKEND\n")
cat("  capabilities(\"cairo\"): ", has_cairo, "\n", sep = "")
cat("  capabilities(\"png\"):   ", has_png, "\n", sep = "")
cat("  ragg installed:        ", have_ragg, "\n", sep = "")
cat("  svglite installed:     ", have_svglite, "\n", sep = "")

cat("\nFILES USED BY R\n")
cat("  Input folder: ", normalizePath(INPUT_DIR), "\n", sep = "")
cat("  Guard mode:   ", GUARD_MODE, "\n", sep = "")
cat("  Output:       ", normalizePath(OUTPUT_DIR, mustWork = FALSE), "\n\n", sep = "")


# ==============================================================================
# 3. VERIFICATION HELPERS
# ==============================================================================

EXPECTED <- list(
  encounters   = 2000000,
  docks        = 6079,
  first_p_dock = 0.00971998056003888,
  last_p_dock  = 0.00048001344037633056,
  beta_T       = 3.42239207834849,
  beta_M       = 2.915147227029582,
  beta_d2      = -0.19045122928569386
)

if (have_jsonlite && file.exists(manifest_file)) {
  manifest <- jsonlite::fromJSON(manifest_file)
  ref <- manifest$canonical_reference
  if (!is.null(ref)) {
    if (!is.null(ref$first_bin_p_dock)) EXPECTED$first_p_dock <- ref$first_bin_p_dock
    if (!is.null(ref$last_bin_p_dock))  EXPECTED$last_p_dock  <- ref$last_bin_p_dock
    if (!is.null(ref$beta_T))           EXPECTED$beta_T       <- ref$beta_T
    if (!is.null(ref$beta_M))           EXPECTED$beta_M       <- ref$beta_M
    if (!is.null(ref$beta_d2))          EXPECTED$beta_d2      <- ref$beta_d2
  }
  if (!is.null(manifest$run_config$CANONICAL_ENCOUNTERS)) {
    EXPECTED$encounters <- manifest$run_config$CANONICAL_ENCOUNTERS
  }
  if (!is.null(manifest$run_config$CANONICAL_DOCKS)) {
    EXPECTED$docks <- manifest$run_config$CANONICAL_DOCKS
  }
  cat("Expected values loaded from MANIFEST_figure2.json\n\n")
} else {
  cat("Using built-in expected values (no manifest / jsonlite available)\n\n")
}

guard_log <- character(0)

check_value <- function(label, observed, expected, tolerance = 1e-10) {
  ok <- isTRUE(all.equal(unname(observed), unname(expected), tolerance = tolerance))
  guard_log[[length(guard_log) + 1]] <<- sprintf(
    "%-34s %-5s observed %s | expected %s",
    label,
    if (ok) "PASS" else "FAIL",
    paste(format(observed, digits = 12), collapse = ", "),
    paste(format(expected, digits = 12), collapse = ", ")
  )
  if (!ok) {
    msg <- sprintf("%s does not match the reference dataset (observed %s, expected %s).",
                   label,
                   paste(format(observed, digits = 12), collapse = ", "),
                   paste(format(expected, digits = 12), collapse = ", "))
    if (GUARD_MODE == "strict") stop(msg) else warning(msg, call. = FALSE)
  }
  invisible(ok)
}


# ==============================================================================
# 4. STYLE CONSTANTS  (unchanged from the original figure2.R)
# ==============================================================================

font_family <- "sans"

colour_teal_dark  <- "#007C83"
colour_teal_light <- "#72C7C2"

colour_activity <- "#D55E00"
colour_memory   <- "#7651A8"
colour_mismatch <- "#4F565A"

colour_text <- "#222222"
colour_axis <- "#2D2D2D"
colour_grid <- "#E5E8EA"
colour_zero <- "#73797D"

colour_inactive <- "#4F565A"
colour_active   <- "#D55E00"

variant_colours <- c(
  full_msc_feedback          = "#007C83",
  msc_capture_only           = "#72C7C2",
  projected1d_feedback       = "#0072B2",
  shuffled_material_feedback = "#E69F00",
  random_material_feedback   = "#CC79A7",
  tag_only_feedback          = "#4F565A"
)

variant_labels <- c(
  full_msc_feedback          = "Full MSC + feedback",
  msc_capture_only           = "MSC capture only",
  projected1d_feedback       = "Projected 1D + feedback",
  shuffled_material_feedback = "Shuffled material",
  random_material_feedback   = "Random material",
  tag_only_feedback          = "Activity tag only"
)

base_size       <- 9.5
title_size      <- 11.5
axis_title_size <- 10.5
axis_text_size  <- 9.5


# ==============================================================================
# 5. COMMON FIGURE THEME
# ==============================================================================

theme_msc <- function() {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = title_size, face = "bold", colour = colour_text,
        hjust = 0, lineheight = 1.0, margin = margin(b = 8)
      ),
      axis.title.x = element_text(
        size = axis_title_size, face = "plain", colour = colour_text,
        margin = margin(t = 7)
      ),
      axis.title.y = element_text(
        size = axis_title_size, face = "plain", colour = colour_text,
        lineheight = 0.95, margin = margin(r = 8)
      ),
      axis.text = element_text(size = axis_text_size, colour = colour_text),
      axis.line = element_line(
        colour = colour_axis, linewidth = 0.55, lineend = "square"
      ),
      axis.ticks = element_line(colour = colour_axis, linewidth = 0.45),
      axis.ticks.length = grid::unit(2, "mm"),
      panel.grid.major.y = element_line(colour = colour_grid, linewidth = 0.45),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 7, r = 10, b = 8, l = 9)
    )
}

theme_msc_legend <- function() {
  theme_msc() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.5, colour = colour_text),
      legend.key.height = grid::unit(4.2, "mm"),
      legend.key.width = grid::unit(5.5, "mm"),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}


# ==============================================================================
# 6. ROBUST EXPORT HELPER
#
# One panel, three formats, each attempted independently. A failure in one
# format is recorded and reported but never stops the others, and never stops
# the script. Every write is verified to exist and be non-empty afterwards.
# ==============================================================================

export_log <- data.frame(
  panel = character(0), format = character(0), device = character(0),
  status = character(0), bytes = numeric(0), message = character(0),
  stringsAsFactors = FALSE
)

# Attempt a list of (device label, ggsave call) pairs until one succeeds.
try_devices <- function(panel_name, fmt, file_path, plot_object,
                        width_mm, height_mm, attempts) {

  for (att in attempts) {
    result <- tryCatch({
      att$fun(file_path, plot_object, width_mm, height_mm)
      # Verify the file really landed and is not a zero-byte stub.
      if (!file.exists(file_path)) stop("device reported success but no file was written")
      sz <- file.info(file_path)$size
      if (is.na(sz) || sz <= 0) stop("file written but empty (0 bytes)")
      list(ok = TRUE, device = att$label, bytes = sz, message = "")
    }, error = function(e) {
      if (file.exists(file_path)) try(file.remove(file_path), silent = TRUE)
      list(ok = FALSE, device = att$label, bytes = 0,
           message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      export_log <<- rbind(export_log, data.frame(
        panel = panel_name, format = fmt, device = result$device,
        status = "OK", bytes = result$bytes, message = "",
        stringsAsFactors = FALSE
      ))
      return(invisible(TRUE))
    } else {
      # Record the failed attempt, then fall through to the next device.
      export_log <<- rbind(export_log, data.frame(
        panel = panel_name, format = fmt, device = result$device,
        status = "FAILED", bytes = 0, message = result$message,
        stringsAsFactors = FALSE
      ))
    }
  }
  invisible(FALSE)
}

# ---- PDF attempts -----------------------------------------------------------

pdf_attempts <- list(
  list(
    label = "cairo_pdf",
    fun = function(f, p, w, h) {
      ggsave(filename = f, plot = p, device = grDevices::cairo_pdf,
             width = w, height = h, units = "mm", bg = "white",
             limitsize = FALSE)
    }
  ),
  list(
    label = "grDevices::pdf",
    fun = function(f, p, w, h) {
      ggsave(filename = f, plot = p, device = grDevices::pdf,
             width = w, height = h, units = "mm", bg = "white",
             limitsize = FALSE)
    }
  )
)

# ---- PNG attempts -----------------------------------------------------------
# ragg first: it does not depend on the cairo build and gives the cleanest text.
# Then cairo png. Then the plain png device, which exists on every R build.

png_attempts <- list()

if (have_ragg) {
  png_attempts <- c(png_attempts, list(list(
    label = "ragg::agg_png",
    fun = function(f, p, w, h) {
      ggsave(filename = f, plot = p, device = ragg::agg_png,
             width = w, height = h, units = "mm", dpi = 600,
             bg = "white", limitsize = FALSE)
    }
  )))
}

if (has_cairo) {
  png_attempts <- c(png_attempts, list(list(
    label = "png(type=cairo)",
    fun = function(f, p, w, h) {
      ggsave(filename = f, plot = p,
             device = function(filename, ...) {
               grDevices::png(filename, ..., type = "cairo")
             },
             width = w, height = h, units = "mm", dpi = 600,
             bg = "white", limitsize = FALSE)
    }
  )))
}

png_attempts <- c(png_attempts, list(list(
  label = "grDevices::png",
  fun = function(f, p, w, h) {
    ggsave(filename = f, plot = p, device = grDevices::png,
           width = w, height = h, units = "mm", dpi = 600,
           bg = "white", limitsize = FALSE)
  }
)))

# ---- SVG attempts -----------------------------------------------------------

svg_attempts <- list()

if (have_svglite) {
  svg_attempts <- c(svg_attempts, list(list(
    label = "svglite::svglite",
    fun = function(f, p, w, h) {
      ggsave(filename = f, plot = p, device = svglite::svglite,
             width = w, height = h, units = "mm", bg = "white",
             limitsize = FALSE)
    }
  )))
}

svg_attempts <- c(svg_attempts, list(list(
  label = "grDevices::svg",
  fun = function(f, p, w, h) {
    ggsave(filename = f, plot = p, device = grDevices::svg,
           width = w, height = h, units = "mm", bg = "white",
           limitsize = FALSE)
  }
)))

# ---- one call, three formats ------------------------------------------------

save_panel <- function(plot_object, filename, width_mm, height_mm) {

  if (is.null(plot_object)) {
    cat(sprintf("skipped  %-34s (panel not built)\n", filename))
    return(invisible(NULL))
  }

  ok_pdf <- try_devices(filename, "pdf",
                        file.path(OUTPUT_DIR, "pdf", paste0(filename, ".pdf")),
                        plot_object, width_mm, height_mm, pdf_attempts)
  ok_png <- try_devices(filename, "png",
                        file.path(OUTPUT_DIR, "png", paste0(filename, ".png")),
                        plot_object, width_mm, height_mm, png_attempts)
  ok_svg <- try_devices(filename, "svg",
                        file.path(OUTPUT_DIR, "svg", paste0(filename, ".svg")),
                        plot_object, width_mm, height_mm, svg_attempts)

  cat(sprintf("saved    %-34s pdf:%-3s png:%-3s svg:%-3s (%g x %g mm)\n",
              filename,
              if (ok_pdf) "ok" else "NO",
              if (ok_png) "ok" else "NO",
              if (ok_svg) "ok" else "NO",
              width_mm, height_mm))
  invisible(NULL)
}


# ==============================================================================
# 7. PANEL A - read, validate, add Wilson intervals
# ==============================================================================

mismatch_data <- read_csv(mismatch_file, show_col_types = FALSE)

required_A_columns <- c("bin", "d2_mean", "p_dock", "n")
missing_A_columns <- setdiff(required_A_columns, names(mismatch_data))
if (length(missing_A_columns) > 0) {
  stop("Panel A is missing required columns: ",
       paste(missing_A_columns, collapse = ", "))
}

mismatch_data <- mismatch_data %>%
  transmute(
    bin = as.character(bin),
    d2_mean = as.numeric(d2_mean),
    p_dock = as.numeric(p_dock),
    n = as.numeric(n)
  ) %>%
  arrange(d2_mean)

if (nrow(mismatch_data) == 0) stop("Panel A contains no rows.")
if (nrow(mismatch_data) != 12) {
  stop("Panel A must contain 12 mismatch bins, found ", nrow(mismatch_data), ".")
}
if (any(!is.finite(mismatch_data$d2_mean)) ||
    any(!is.finite(mismatch_data$p_dock)) ||
    any(!is.finite(mismatch_data$n))) {
  stop("Panel A contains missing or non-finite values.")
}
if (any(mismatch_data$p_dock < 0) || any(mismatch_data$p_dock > 1)) {
  stop("Panel A docking probabilities must lie between 0 and 1.")
}
if (any(mismatch_data$n <= 0)) stop("Panel A encounter counts must be positive.")

total_encounters_A <- sum(mismatch_data$n)
check_value("Panel A total encounters", total_encounters_A,
            EXPECTED$encounters, tolerance = 0)
check_value("Panel A first-bin p_dock", mismatch_data$p_dock[1],
            EXPECTED$first_p_dock)
check_value("Panel A last-bin p_dock",
            mismatch_data$p_dock[nrow(mismatch_data)], EXPECTED$last_p_dock)

# Approximate Wilson binomial intervals, for descriptive visualisation only.
# These are event-level intervals over simulated encounter opportunities and do
# NOT represent variation across independent biological replicates.
z_95 <- qnorm(0.975)

mismatch_data <- mismatch_data %>%
  mutate(
    estimated_docks = p_dock * n,
    wilson_denominator = 1 + z_95^2 / n,
    wilson_centre = (p_dock + z_95^2 / (2 * n)) / wilson_denominator,
    wilson_half_width = (
      z_95 * sqrt(p_dock * (1 - p_dock) / n + z_95^2 / (4 * n^2))
    ) / wilson_denominator,
    ci_low = pmax(0, wilson_centre - wilson_half_width),
    ci_high = pmin(1, wilson_centre + wilson_half_width)
  )

x_upper_A <- ceiling(max(mismatch_data$d2_mean) / 3) * 3
if (!is.finite(x_upper_A) || x_upper_A <= 0) x_upper_A <- 12

y_upper_A <- ceiling(max(mismatch_data$ci_high) * 1000) / 1000
y_upper_A <- max(y_upper_A * 1.08, 0.004)


# ==============================================================================
# 8. PANEL A - plot
# ==============================================================================

panel_A <- ggplot(mismatch_data, aes(x = d2_mean, y = p_dock)) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0, linewidth = 0.45, colour = colour_teal_dark, alpha = 0.80
  ) +
  geom_line(
    colour = colour_teal_dark, linewidth = 0.90,
    lineend = "round", linejoin = "round"
  ) +
  geom_point(
    shape = 21, size = 2.70, stroke = 0.65,
    colour = colour_teal_dark, fill = colour_teal_light
  ) +
  scale_x_continuous(
    name = expression("Material mismatch, " * d[ig]^2),
    breaks = pretty_breaks(n = 5),
    limits = c(0, x_upper_A),
    expand = expansion(mult = c(0.015, 0.02))
  ) +
  scale_y_continuous(
    name = "Docking probability\nper encounter",
    labels = label_percent(accuracy = 0.1),
    breaks = pretty_breaks(n = 5),
    limits = c(0, y_upper_A),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(title = "Docking probability declines with material mismatch") +
  theme_msc() +
  theme(plot.margin = margin(t = 8, r = 11, b = 8, l = 12))


# ==============================================================================
# 9. PANELS B AND B2 - observation model
# ==============================================================================

observation_data <- read_csv(observation_file, show_col_types = FALSE)

required_B_columns <- c("term", "estimate", "std_error", "z", "pvalue",
                        "n_rows", "n_docks")
missing_B_columns <- setdiff(required_B_columns, names(observation_data))
if (length(missing_B_columns) > 0) {
  stop("Panel B is missing required columns: ",
       paste(missing_B_columns, collapse = ", "))
}

observation_data <- observation_data %>%
  mutate(
    term = as.character(term),
    estimate = as.numeric(estimate),
    std_error = as.numeric(std_error),
    z = as.numeric(z),
    pvalue = as.numeric(pvalue),
    n_rows = as.numeric(n_rows),
    n_docks = as.numeric(n_docks)
  )

n_rows_B  <- unique(observation_data$n_rows)
n_docks_B <- unique(observation_data$n_docks)
if (length(n_rows_B) != 1 || length(n_docks_B) != 1) {
  stop("Panel B contains inconsistent n_rows or n_docks values.")
}

check_value("Panel B observations", n_rows_B, EXPECTED$encounters, tolerance = 0)
check_value("Panel B docking events", n_docks_B, EXPECTED$docks, tolerance = 0)

label_terms <- function(df) {
  df %>%
    mutate(
      ci_low = estimate - 1.96 * std_error,
      ci_high = estimate + 1.96 * std_error,
      term_label = case_when(
        term == "T"              ~ "Activity tag, T",
        term == "M"              ~ "Material-memory tag, M",
        term == "d2"             ~ "Material mismatch, d\u00b2",
        term == "active_synapse" ~ "Active-synapse indicator",
        TRUE                     ~ term
      )
    )
}

coefficient_axis <- function(df) {
  x_min <- min(df$ci_low, 0)
  x_max <- max(df$ci_high, 0)
  x_range <- x_max - x_min
  if (!is.finite(x_range) || x_range <= 0) x_range <- 1
  c(min(x_min - 0.10 * x_range, -0.35), x_max + 0.10 * x_range)
}

observation_plot_data <- observation_data %>%
  filter(term %in% c("T", "M", "d2")) %>%
  label_terms() %>%
  mutate(
    term_label = factor(
      term_label,
      levels = c("Material mismatch, d\u00b2",
                 "Material-memory tag, M",
                 "Activity tag, T")
    ),
    colour_group = factor(term, levels = c("T", "M", "d2"))
  )

if (nrow(observation_plot_data) != 3) {
  stop("Panel B must contain exactly the three mechanistic predictors T, M and d2.")
}
if (any(!is.finite(observation_plot_data$estimate)) ||
    any(!is.finite(observation_plot_data$std_error))) {
  stop("Panel B contains missing or non-finite model values.")
}
if (any(observation_plot_data$std_error <= 0)) {
  stop("Panel B contains a non-positive standard error.")
}

expected_B <- c(T = EXPECTED$beta_T, M = EXPECTED$beta_M, d2 = EXPECTED$beta_d2)
observed_B <- observation_plot_data$estimate[
  match(names(expected_B), observation_plot_data$term)
]
check_value("Panel B coefficients (T, M, d2)", observed_B, expected_B)

lim_B <- coefficient_axis(observation_plot_data)

panel_B <- ggplot(observation_plot_data, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = colour_zero, linewidth = 0.60) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, y = term_label, yend = term_label,
        colour = colour_group),
    linewidth = 1.15, lineend = "round"
  ) +
  geom_point(aes(colour = colour_group), size = 3.60) +
  scale_colour_manual(
    values = c(T = colour_activity, M = colour_memory, d2 = colour_mismatch),
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Logistic coefficient estimate",
    breaks = pretty_breaks(n = 5), limits = lim_B,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.55, 0.55))) +
  labs(title = "Observation model recovers encoded effect directions") +
  theme_msc() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9.5, colour = colour_text,
                               hjust = 1, margin = margin(r = 7)),
    panel.grid.major.y = element_line(colour = colour_grid, linewidth = 0.45),
    plot.margin = margin(t = 8, r = 14, b = 8, l = 12)
  )

observation_full_data <- observation_data %>%
  filter(term != "Intercept") %>%
  label_terms() %>%
  mutate(
    term_label = factor(
      term_label,
      levels = c("Active-synapse indicator",
                 "Material mismatch, d\u00b2",
                 "Material-memory tag, M",
                 "Activity tag, T")
    ),
    colour_group = factor(term, levels = c("T", "M", "d2", "active_synapse"))
  )

lim_B2 <- coefficient_axis(observation_full_data)

panel_B2 <- ggplot(observation_full_data, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = colour_zero, linewidth = 0.60) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, y = term_label, yend = term_label,
        colour = colour_group),
    linewidth = 1.15, lineend = "round"
  ) +
  geom_point(aes(colour = colour_group), size = 3.60) +
  scale_colour_manual(
    values = c(T = colour_activity, M = colour_memory,
               d2 = colour_mismatch, active_synapse = "#0072B2"),
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Logistic coefficient estimate",
    breaks = pretty_breaks(n = 5), limits = lim_B2,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.55, 0.55))) +
  labs(title = "Complete fitted model, including the active-synapse indicator") +
  theme_msc() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9.5, colour = colour_text,
                               hjust = 1, margin = margin(r = 7)),
    panel.grid.major.y = element_line(colour = colour_grid, linewidth = 0.45),
    plot.margin = margin(t = 8, r = 14, b = 8, l = 12)
  )


# ==============================================================================
# 10. PANEL C - kernel consistency
# ==============================================================================

panel_C <- NULL
kernel_data <- NULL

if (file.exists(kernel_file)) {
  kernel_data <- read_csv(kernel_file, show_col_types = FALSE) %>%
    arrange(d2_mean)
  if (all(c("d2_mean", "obs_over_pred") %in% names(kernel_data))) {
    panel_C <- ggplot(kernel_data, aes(x = d2_mean, y = obs_over_pred)) +
      geom_hline(yintercept = 1, linetype = "dashed",
                 colour = colour_zero, linewidth = 0.60) +
      geom_line(colour = colour_teal_dark, linewidth = 0.80) +
      geom_point(shape = 21, size = 2.70, stroke = 0.65,
                 colour = colour_teal_dark, fill = colour_teal_light) +
      scale_x_continuous(
        name = expression("Material mismatch, " * d[ig]^2),
        breaks = pretty_breaks(n = 5),
        expand = expansion(mult = c(0.02, 0.03))
      ) +
      scale_y_continuous(
        name = "Observed / kernel-predicted\ndocking probability",
        breaks = pretty_breaks(n = 5),
        expand = expansion(mult = c(0.05, 0.08))
      ) +
      labs(title = "Simulated docking follows the encoded compatibility kernel") +
      theme_msc() +
      theme(plot.margin = margin(t = 8, r = 11, b = 8, l = 12))
  } else {
    warning("FigS_panelA_kernel_consistency.csv lacks expected columns; skipping panel C.",
            call. = FALSE)
    kernel_data <- NULL
  }
} else {
  cat("Panel C input not found; skipping.\n")
}


# ==============================================================================
# 11. PANEL D - Panel A stratified by synapse state
# ==============================================================================

panel_D <- NULL
stratified_data <- NULL

if (file.exists(stratified_file)) {
  stratified_data <- read_csv(stratified_file, show_col_types = FALSE)
  if (all(c("active_synapse", "d2_mean", "p_dock", "n") %in% names(stratified_data))) {
    stratified_data <- stratified_data %>%
      filter(is.finite(p_dock), p_dock > 0) %>%
      mutate(
        stratum = factor(
          ifelse(active_synapse == 1, "Active synapse", "Inactive synapses"),
          levels = c("Inactive synapses", "Active synapse")
        )
      ) %>%
      arrange(stratum, d2_mean)

    panel_D <- ggplot(
      stratified_data,
      aes(x = d2_mean, y = p_dock, colour = stratum, fill = stratum)
    ) +
      geom_line(linewidth = 0.85, lineend = "round") +
      geom_point(shape = 21, size = 2.40, stroke = 0.60, colour = "white") +
      scale_colour_manual(values = c("Inactive synapses" = colour_inactive,
                                     "Active synapse" = colour_active)) +
      scale_fill_manual(values = c("Inactive synapses" = colour_inactive,
                                   "Active synapse" = colour_active)) +
      scale_x_continuous(
        name = expression("Material mismatch, " * d[ig]^2),
        breaks = pretty_breaks(n = 5),
        expand = expansion(mult = c(0.02, 0.03))
      ) +
      scale_y_log10(
        name = "Docking probability\nper encounter",
        labels = label_percent(accuracy = 0.01),
        expand = expansion(mult = c(0.06, 0.06))
      ) +
      labs(title = "Stratification removes the lowest-bin composition effect") +
      theme_msc_legend() +
      theme(plot.margin = margin(t = 8, r = 11, b = 8, l = 12))
  } else {
    warning("FigS_panelA_stratified_active.csv lacks expected columns; skipping panel D.",
            call. = FALSE)
    stratified_data <- NULL
  }
} else {
  cat("Panel D input not found; skipping.\n")
}


# ==============================================================================
# 12. PANEL E - d^2 dependence across all six model variants
# ==============================================================================

panel_E <- NULL
tautology_data <- NULL

if (file.exists(tautology_file)) {
  tautology_data <- read_csv(tautology_file, show_col_types = FALSE)
  if (all(c("variant", "d2_mean", "p_dock") %in% names(tautology_data))) {
    tautology_data <- tautology_data %>%
      filter(is.finite(p_dock), p_dock > 0) %>%
      mutate(variant = factor(variant, levels = names(variant_colours))) %>%
      filter(!is.na(variant)) %>%
      arrange(variant, d2_mean)

    panel_E <- ggplot(tautology_data, aes(x = d2_mean, y = p_dock, colour = variant)) +
      geom_line(linewidth = 0.80, lineend = "round") +
      geom_point(size = 1.60) +
      scale_colour_manual(values = variant_colours,
                          labels = variant_labels, drop = FALSE) +
      scale_x_continuous(
        name = expression("Material mismatch, " * d[ig]^2),
        breaks = pretty_breaks(n = 5),
        expand = expansion(mult = c(0.02, 0.03))
      ) +
      scale_y_log10(
        name = "Docking probability\nper encounter",
        labels = label_percent(accuracy = 0.01),
        expand = expansion(mult = c(0.06, 0.06))
      ) +
      labs(title = "Mismatch dependence is confined to material-aware variants") +
      theme_msc_legend() +
      theme(plot.margin = margin(t = 8, r = 11, b = 8, l = 12))
  } else {
    warning("tautology_binned_by_variant.csv lacks expected columns; skipping panel E.",
            call. = FALSE)
    tautology_data <- NULL
  }
} else {
  cat("Panel E input not found; skipping.\n")
}


# ==============================================================================
# 13. COMBINED PREVIEWS
# ==============================================================================

main_preview <- panel_A + panel_B +
  plot_layout(widths = c(1.00, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(family = font_family, size = 12,
                            face = "bold", colour = colour_text),
    plot.tag.position = c(0.005, 0.995)
  )

supp_panels <- Filter(Negate(is.null), list(panel_B2, panel_C, panel_D, panel_E))

supp_preview <- if (length(supp_panels) > 0) {
  wrap_plots(supp_panels, ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(family = font_family, size = 12,
                              face = "bold", colour = colour_text),
      plot.tag.position = c(0.005, 0.995)
    )
} else {
  NULL
}


# ==============================================================================
# 14. EXPORT ALL PANELS AND PREVIEWS
# ==============================================================================

cat("\nEXPORTING PANELS\n")

save_panel(panel_A,  "Figure2A_mismatch_docking",       125, 82)
save_panel(panel_B,  "Figure2B_observation_model",      150, 82)
save_panel(panel_B2, "FigS2B2_observation_model_full",  165, 88)
save_panel(panel_C,  "FigS2C_kernel_consistency",       125, 82)
save_panel(panel_D,  "FigS2D_stratified_by_synapse",    150, 82)
save_panel(panel_E,  "FigS2E_variant_mismatch_curves",  165, 88)
save_panel(main_preview, "Figure2_AB_preview", 275, 88)

if (!is.null(supp_preview)) {
  save_panel(supp_preview, "FigS2_supplementary_preview", 330, 185)
}


# ==============================================================================
# 15. PLOTTED VALUES, GUARD LOG, EXPORT LOG AND SESSION INFO
# ==============================================================================

data_dir <- file.path(OUTPUT_DIR, "plotted_data")

write_csv(
  mismatch_data %>%
    select(bin, d2_mean, p_dock, ci_low, ci_high, n, estimated_docks),
  file.path(data_dir, "Figure2A_plotted_data.csv")
)
write_csv(
  observation_plot_data %>%
    select(term, estimate, std_error, ci_low, ci_high, z, pvalue),
  file.path(data_dir, "Figure2B_plotted_data.csv")
)
write_csv(
  observation_full_data %>%
    select(term, estimate, std_error, ci_low, ci_high, z, pvalue),
  file.path(data_dir, "FigS2B2_plotted_data.csv")
)
if (!is.null(kernel_data)) {
  write_csv(kernel_data, file.path(data_dir, "FigS2C_plotted_data.csv"))
}
if (!is.null(stratified_data)) {
  write_csv(stratified_data, file.path(data_dir, "FigS2D_plotted_data.csv"))
}
if (!is.null(tautology_data)) {
  write_csv(tautology_data, file.path(data_dir, "FigS2E_plotted_data.csv"))
}

write_csv(export_log, file.path(OUTPUT_DIR, "Figure2_R_export_log.csv"))
writeLines(guard_log, con = file.path(OUTPUT_DIR, "Figure2_R_guard_log.txt"))
writeLines(capture.output(sessionInfo()),
           con = file.path(OUTPUT_DIR, "Figure2_R_sessionInfo.txt"))


# ==============================================================================
# 16. NUMERIC SUMMARY
# ==============================================================================

calculated_docking_frequency <- n_docks_B / n_rows_B

get_coef <- function(df, term_name) df$estimate[df$term == term_name]

cat("\n")
cat("FIGURE 2 NUMERIC SUMMARY\n")
cat("------------------------\n")
cat("Panel A bins:                 ", nrow(mismatch_data), "\n")
cat("Panel A encounters:           ",
    format(total_encounters_A, big.mark = ",", scientific = FALSE), "\n")
cat("Panel B observations:         ",
    format(n_rows_B, big.mark = ",", scientific = FALSE), "\n")
cat("Panel B docking events:       ",
    format(n_docks_B, big.mark = ",", scientific = FALSE), "\n")
cat("Overall docking frequency:    ",
    sprintf("%.4f %%", 100 * calculated_docking_frequency), "\n")
cat("First-bin docking probability:",
    sprintf("%.4f %%", 100 * mismatch_data$p_dock[1]), "\n")
cat("Last-bin docking probability: ",
    sprintf("%.4f %%", 100 * mismatch_data$p_dock[nrow(mismatch_data)]), "\n")
cat("T coefficient:                ",
    sprintf("%.6f", get_coef(observation_plot_data, "T")), "\n")
cat("M coefficient:                ",
    sprintf("%.6f", get_coef(observation_plot_data, "M")), "\n")
cat("d2 coefficient:               ",
    sprintf("%.7f", get_coef(observation_plot_data, "d2")), "\n")

if ("active_synapse" %in% observation_full_data$term) {
  cat("active-synapse indicator:     ",
      sprintf("%.6f", get_coef(observation_full_data, "active_synapse")),
      " (shown in FigS2B2; not part of the docking hazard)\n")
}

cat("\nVERSION GUARDS (mode: ", GUARD_MODE, ")\n", sep = "")
cat(paste(guard_log, collapse = "\n"), "\n")

if (file.exists(numbers_file)) {
  cat("\nPAPER-FACING NUMBERS FROM THE NOTEBOOK\n")
  print(read_csv(numbers_file, show_col_types = FALSE), n = Inf, width = Inf)
}


# ==============================================================================
# 17. EXPORT STATUS TABLE
#
# If a format shows zero files, the failing device and its error message are
# printed here and stored in Figure2_R_export_log.csv.
# ==============================================================================

cat("\nEXPORT STATUS BY FORMAT\n")
for (fmt in c("pdf", "png", "svg")) {
  sub <- export_log[export_log$format == fmt, , drop = FALSE]
  n_ok <- sum(sub$status == "OK")
  devs <- unique(sub$device[sub$status == "OK"])
  total_mb <- sum(sub$bytes[sub$status == "OK"]) / 1e6
  cat(sprintf("  %-4s %2d file(s)  device: %-18s %6.2f MB\n",
              toupper(fmt), n_ok,
              if (length(devs)) paste(devs, collapse = ", ") else "-",
              total_mb))
  if (n_ok == 0) {
    cat("       NO FILES WRITTEN. Failed attempts:\n")
    fails <- sub[sub$status == "FAILED", , drop = FALSE]
    if (nrow(fails) == 0) {
      cat("       (no device was available for this format)\n")
    } else {
      for (i in seq_len(nrow(fails))) {
        cat(sprintf("       %-18s %s\n", fails$device[i], fails$message[i]))
      }
    }
    if (fmt == "png") {
      cat("       Hint: install.packages(\"ragg\") gives a cairo-independent PNG device.\n")
    }
    if (fmt == "pdf") {
      cat("       Hint: cairo_pdf needs a cairo-enabled R build; the plain pdf device is the fallback.\n")
    }
  }
}

n_actual <- length(list.files(OUTPUT_DIR, recursive = TRUE,
                              pattern = "\\.(pdf|png|svg)$"))
cat(sprintf("\nGraphics files on disk: %d\n", n_actual))

cat("\nFILES WRITTEN\n")
cat("  pdf/                      vector, for submission\n")
cat("  png/                      600 dpi raster\n")
cat("  svg/                      vector, editable text\n")
cat("  plotted_data/             the exact numbers drawn in each panel\n")
cat("  Figure2_R_export_log.csv  per-file device and status\n")
cat("  Figure2_R_guard_log.txt   reproduction checks\n")
cat("  Figure2_R_sessionInfo.txt R and package versions\n")

cat("\nEXPORT COMPLETED\n")
cat("Output folder:",
    normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("Individual panel exports contain no panel letters.\n")
cat("The previews carry temporary tags for layout checking only.\n")
