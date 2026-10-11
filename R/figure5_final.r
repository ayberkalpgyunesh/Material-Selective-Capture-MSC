# ==============================================================================
# MSC Figure 5 — Canonical anti-tautology analysis
# ==============================================================================
# Input: 04_r_input from MSC_Figure5_reproducible.ipynb
# Usage: Rscript figure5_final_v3.r /path/to/04_r_input
#
# This script does not simulate. It reads the 6-million-row parquet dataset,
# fits the pre-specified anti-tautology models, writes canonical analysis tables,
# and exports Figure 5B-D as PDF, PNG, SVG and TIFF.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Install and load packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "arrow", "data.table", "fixest", "ggplot2", "scales", "jsonlite"
)

optional_packages <- c("ragg", "svglite")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

for (package_name in optional_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    try(install.packages(package_name, repos = "https://cloud.r-project.org"), silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(scales)
  library(jsonlite)
})

have_ragg <- requireNamespace("ragg", quietly = TRUE)
have_svglite <- requireNamespace("svglite", quietly = TRUE)

# ------------------------------------------------------------------------------
# 2. Paths and output folders
# ------------------------------------------------------------------------------

arguments <- commandArgs(trailingOnly = TRUE)

INPUT_DIR <- if (length(arguments) >= 1 && nzchar(arguments[1])) {
  arguments[1]
} else {
  getwd()
}

OUTPUT_DIR <- file.path(INPUT_DIR, "R_output_figure5_v3")

for (folder in c("pdf", "png", "svg", "tiff", "tables")) {
  dir.create(file.path(OUTPUT_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}

input_file <- function(filename) file.path(INPUT_DIR, filename)

raw_file <- input_file("encounter_observations_all_variants.parquet")
plan_file <- input_file("analysis_plan.json")
diagnostic_file <- input_file("Figure5_true_used_separation_diagnostics.csv")

stopifnot(file.exists(raw_file))
stopifnot(file.exists(plan_file))
stopifnot(file.exists(diagnostic_file))

# ------------------------------------------------------------------------------
# 3. Design definitions
# ------------------------------------------------------------------------------

plan <- jsonlite::fromJSON(plan_file)

variant_order <- plan$variants
equivalence_bound <- as.numeric(plan$equivalence_bound_logit)
alpha_level <- 0.05

material_preserving <- c(
  "full_msc_feedback",
  "msc_capture_only",
  "projected1d_feedback"
)

joint_model_variants <- c(
  "shuffled_material_feedback",
  "random_material_feedback"
)

equivalence_controls <- c(
  "random_material_feedback",
  "tag_only_feedback"
)

variant_labels <- c(
  full_msc_feedback = "Full MSC + feedback",
  msc_capture_only = "MSC capture only",
  projected1d_feedback = "Projected 1D + feedback",
  shuffled_material_feedback = "Shuffled material",
  random_material_feedback = "Random material",
  tag_only_feedback = "Activity only"
)

variant_group <- function(variant) {
  fifelse(
    variant %in% material_preserving,
    "Material preserving",
    fifelse(
      variant == "shuffled_material_feedback",
      "Intermediate",
      "Material disrupted"
    )
  )
}

group_colours <- c(
  "Material preserving" = "#2C7FB8",
  "Intermediate" = "#E66101",
  "Material disrupted" = "#777777"
)

mismatch_colours <- c(
  "True mismatch" = "#2166AC",
  "Model-used mismatch" = "#B2182B"
)

# ------------------------------------------------------------------------------
# 4. Load raw encounter observations and validate the dataset
# ------------------------------------------------------------------------------

cat("Reading parquet dataset...\n")
encounters <- as.data.table(arrow::read_parquet(raw_file, as_data_frame = FALSE))

required_columns <- c(
  "seed", "variant", "T", "M", "d2", "d2_true", "d2_used",
  "active_synapse", "dock"
)

missing_columns <- setdiff(required_columns, names(encounters))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

encounters[, variant := factor(as.character(variant), levels = variant_order)]
encounters[, seed := factor(seed)]

# Only these columns are numeric. Factors seed and variant must not be tested
# with is.finite(), because that would generate a false non-finite error.
numeric_columns <- c(
  "T", "M", "d2", "d2_true", "d2_used", "active_synapse", "dock"
)

nonfinite_counts <- encounters[
  ,
  lapply(.SD, function(x) sum(!is.finite(x))),
  .SDcols = numeric_columns
]

if (sum(unlist(nonfinite_counts)) > 0) {
  print(nonfinite_counts)
  stop("Non-finite values detected in numeric encounter columns.")
}

input_counts <- encounters[
  ,
  .(
    n_encounters = .N,
    n_seeds = uniqueN(seed),
    n_docks = sum(dock),
    docking_rate = mean(dock)
  ),
  by = variant
]

if (nrow(encounters) != 6000000L) {
  stop("Expected 6,000,000 total observations; found ", nrow(encounters))
}

if (any(input_counts$n_encounters != 1000000L)) {
  stop("Every variant must contain exactly 1,000,000 observations.")
}

if (any(input_counts$n_seeds != 100L)) {
  stop("Every variant must contain exactly 100 seeds.")
}

if (any(input_counts$n_docks == 0L)) {
  stop("At least one variant contains no docking events.")
}

if (any(abs(encounters$d2 - encounters$d2_used) > 1e-12)) {
  stop("The compatibility alias d2 must be identical to d2_used.")
}

cat("Input checks passed.\n")
print(input_counts)

# ------------------------------------------------------------------------------
# 5. Seed-fixed-effect logistic models
# ------------------------------------------------------------------------------

extract_coefficient <- function(model, variant_name, mismatch_name) {

  coefficients <- coef(model)
  standard_errors <- se(model)

  if (!mismatch_name %in% names(coefficients)) {
    stop("Mismatch term was not estimable: ", variant_name, " / ", mismatch_name)
  }

  beta <- unname(coefficients[mismatch_name])
  standard_error <- unname(standard_errors[mismatch_name])
  z_value <- beta / standard_error

  data.table(
    variant = variant_name,
    mismatch = mismatch_name,
    beta = beta,
    se = standard_error,
    z = z_value,
    p = 2 * pnorm(-abs(z_value)),
    ci_low = beta - 1.96 * standard_error,
    ci_high = beta + 1.96 * standard_error,
    n_model = nobs(model),
    omitted_terms = paste(model$collin.var, collapse = ";")
  )
}

fit_variant_model <- function(variant_name, mismatch_variable) {

  variant_data <- encounters[variant == variant_name]

  formula <- as.formula(
    paste0("dock ~ ", mismatch_variable, " + T + M | seed")
  )

  model <- fixest::feglm(
    formula,
    data = variant_data,
    family = "binomial",
    vcov = ~seed,
    fixef.rm = "perfect_fit",
    notes = FALSE
  )

  extract_coefficient(model, variant_name, mismatch_variable)
}

cat("Fitting true-mismatch models...\n")
true_results <- rbindlist(
  lapply(variant_order, fit_variant_model, mismatch_variable = "d2_true")
)

cat("Fitting model-used mismatch models...\n")
used_results <- rbindlist(
  lapply(variant_order, fit_variant_model, mismatch_variable = "d2_used")
)

true_results[, group := variant_group(as.character(variant))]
used_results[, group := variant_group(as.character(variant))]

true_results[, p_holm := p.adjust(p, method = "holm")]
used_results[, p_holm := p.adjust(p, method = "holm")]

# ------------------------------------------------------------------------------
# 6. Holm-adjusted TOST for pre-specified destroyed controls
# ------------------------------------------------------------------------------

run_tost <- function(beta, standard_error, bound) {

  p_lower <- 1 - pnorm((beta + bound) / standard_error)
  p_upper <- pnorm((beta - bound) / standard_error)

  data.table(
    tost_p_lower = p_lower,
    tost_p_upper = p_upper,
    p_tost = max(p_lower, p_upper),
    ci90_low = beta - qnorm(0.95) * standard_error,
    ci90_high = beta + qnorm(0.95) * standard_error
  )
}

true_results[
  ,
  c("tost_p_lower", "tost_p_upper", "p_tost", "p_tost_holm",
    "ci90_low", "ci90_high", "equivalent_to_zero") :=
    list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, FALSE)
]

for (variant_name in equivalence_controls) {

  row_id <- which(as.character(true_results$variant) == variant_name)

  tost_result <- run_tost(
    true_results$beta[row_id],
    true_results$se[row_id],
    equivalence_bound
  )

  true_results[row_id,
    c("tost_p_lower", "tost_p_upper", "p_tost", "ci90_low", "ci90_high") := tost_result
  ]
}

control_rows <- which(as.character(true_results$variant) %in% equivalence_controls)

true_results[control_rows,
  p_tost_holm := p.adjust(p_tost, method = "holm")
]

true_results[control_rows,
  equivalent_to_zero :=
    p_tost_holm < alpha_level &
    ci90_low > -equivalence_bound &
    ci90_high < equivalence_bound
]

if (!all(true_results$equivalent_to_zero[control_rows])) {
  stop("Pre-specified equivalence controls did not pass Holm-adjusted TOST.")
}

# ------------------------------------------------------------------------------
# 7. Joint true-versus-used models after material disruption
# ------------------------------------------------------------------------------

fit_joint_model <- function(variant_name) {

  variant_data <- encounters[variant == variant_name]

  model <- fixest::feglm(
    dock ~ d2_true + d2_used + T + M | seed,
    data = variant_data,
    family = "binomial",
    vcov = ~seed,
    fixef.rm = "perfect_fit",
    notes = FALSE
  )

  output <- rbindlist(list(
    extract_coefficient(model, variant_name, "d2_true"),
    extract_coefficient(model, variant_name, "d2_used")
  ))

  output[, mismatch := fifelse(
    mismatch == "d2_true",
    "True mismatch",
    "Model-used mismatch"
  )]

  output
}

cat("Fitting joint true-versus-used models...\n")
joint_results <- rbindlist(lapply(joint_model_variants, fit_joint_model))
joint_results[, p_holm := p.adjust(p, method = "holm")]

# ------------------------------------------------------------------------------
# 8. IQR odds ratios
# ------------------------------------------------------------------------------

true_quantiles <- quantile(encounters$d2_true, c(0.25, 0.75), type = 7)
true_q25 <- unname(true_quantiles[1])
true_q75 <- unname(true_quantiles[2])
true_iqr <- true_q75 - true_q25

odds_results <- copy(true_results)

# The calculations are deliberately separated. In data.table, a column created
# inside one := call cannot safely be reused by another expression in that same
# call; this avoids the odds_ratio-not-found error from the earlier script.
odds_results[, q25_true := true_q25]
odds_results[, q75_true := true_q75]
odds_results[, true_iqr := true_iqr]
odds_results[, odds_ratio := exp(beta * true_iqr)]
odds_results[, or_low := exp(ci_low * true_iqr)]
odds_results[, or_high := exp(ci_high * true_iqr)]
odds_results[, percent_change_in_odds := 100 * (odds_ratio - 1)]

# ------------------------------------------------------------------------------
# 9. Plotting theme and panels B-D
# ------------------------------------------------------------------------------

theme_msc <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.40),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = element_text(size = base_size + 0.2),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      plot.margin = margin(5, 7, 5, 5)
    )
}

plot_B <- copy(true_results)
plot_B[, label := factor(variant_labels[as.character(variant)],
                         levels = rev(variant_labels[variant_order]))]
plot_B[, group := factor(group, levels = names(group_colours))]

panel_B <- ggplot(plot_B, aes(x = beta, y = label, colour = group, shape = group)) +
  annotate("rect", xmin = -equivalence_bound, xmax = equivalence_bound,
           ymin = -Inf, ymax = Inf, fill = "#EBEBEB", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0, linewidth = 0.80) +
  geom_point(size = 3.2, stroke = 0.85, fill = "white") +
  scale_colour_manual(values = group_colours) +
  scale_shape_manual(values = c("Material preserving" = 21,
                                "Intermediate" = 22,
                                "Material disrupted" = 24)) +
  scale_x_continuous(limits = c(-0.215, 0.04), breaks = seq(-0.20, 0, 0.05)) +
  labs(x = "True-mismatch coefficient", y = NULL) +
  theme_msc() +
  theme(legend.position = "none")

plot_C <- copy(joint_results)
plot_C[, label := factor(variant_labels[as.character(variant)],
                         levels = rev(variant_labels[joint_model_variants]))]
plot_C[, mismatch := factor(mismatch, levels = c("True mismatch", "Model-used mismatch"))]

panel_C <- ggplot(plot_C, aes(x = beta, y = label, colour = mismatch, shape = mismatch)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                position = position_dodge(width = 0.42), width = 0, linewidth = 0.80) +
  geom_point(position = position_dodge(width = 0.42), size = 3.2,
             stroke = 0.85, fill = "white") +
  scale_colour_manual(values = mismatch_colours) +
  scale_shape_manual(values = c("True mismatch" = 21, "Model-used mismatch" = 22)) +
  scale_x_continuous(limits = c(-0.215, 0.04), breaks = seq(-0.20, 0, 0.05)) +
  labs(x = "Joint-model mismatch coefficient", y = NULL) +
  theme_msc()

plot_D <- copy(odds_results)
plot_D[, label := factor(variant_labels[as.character(variant)],
                         levels = rev(variant_labels[variant_order]))]
plot_D[, group := factor(group, levels = names(group_colours))]

panel_D <- ggplot(plot_D, aes(x = odds_ratio, y = label, colour = group, shape = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.50) +
  geom_errorbar(aes(xmin = or_low, xmax = or_high), orientation = "y", width = 0, linewidth = 0.80) +
  geom_point(size = 3.2, stroke = 0.85, fill = "white") +
  scale_colour_manual(values = group_colours) +
  scale_shape_manual(values = c("Material preserving" = 21,
                                "Intermediate" = 22,
                                "Material disrupted" = 24)) +
  scale_x_continuous(limits = c(0.40, 1.04), breaks = seq(0.4, 1.0, 0.1)) +
  labs(x = "Docking odds ratio: Q75 versus Q25 true mismatch", y = NULL) +
  theme_msc() +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# 10. Export helper: each panel in PDF, PNG, SVG and TIFF
# ------------------------------------------------------------------------------

export_log <- data.table(panel = character(), format = character(),
                         status = character(), message = character())

save_one <- function(plot_object, file_stem, format, width_mm, height_mm) {

  filename <- file.path(OUTPUT_DIR, format, paste0(file_stem, ".", format))

  result <- tryCatch({

    device <- switch(
      format,
      pdf = grDevices::pdf,
      png = if (have_ragg) ragg::agg_png else grDevices::png,
      svg = if (have_svglite) svglite::svglite else grDevices::svg,
      tiff = if (have_ragg) ragg::agg_tiff else grDevices::tiff
    )

    arguments <- list(
      filename = filename,
      plot = plot_object,
      device = device,
      width = width_mm,
      height = height_mm,
      units = "mm",
      bg = "white",
      limitsize = FALSE
    )

    if (format %in% c("png", "tiff")) arguments$dpi <- 600
    if (format == "tiff") arguments$compression <- "lzw"

    do.call(ggsave, arguments)

    if (!file.exists(filename) || file.info(filename)$size <= 0) {
      stop("No non-empty output file was written.")
    }

    list(ok = TRUE, message = "")

  }, error = function(error) {
    if (file.exists(filename)) try(file.remove(filename), silent = TRUE)
    list(ok = FALSE, message = conditionMessage(error))
  })

  export_log <<- rbind(
    export_log,
    data.table(
      panel = file_stem,
      format = format,
      status = ifelse(result$ok, "OK", "FAILED"),
      message = result$message
    )
  )

  result$ok
}

export_panel <- function(plot_object, file_stem, width_mm, height_mm) {

  formats <- c("pdf", "png", "svg", "tiff")
  result <- sapply(formats, function(format) {
    save_one(plot_object, file_stem, format, width_mm, height_mm)
  })

  cat(
    sprintf(
      "%-42s PDF:%-4s PNG:%-4s SVG:%-4s TIFF:%-4s\n",
      file_stem,
      ifelse(result["pdf"], "OK", "FAIL"),
      ifelse(result["png"], "OK", "FAIL"),
      ifelse(result["svg"], "OK", "FAIL"),
      ifelse(result["tiff"], "OK", "FAIL")
    )
  )
}

export_panel(panel_B, "Figure5B_true_mismatch_coefficient", 112, 68)
export_panel(panel_C, "Figure5C_true_vs_used_joint_model", 125, 72)
export_panel(panel_D, "Figure5D_true_mismatch_odds_ratio", 112, 68)

# ------------------------------------------------------------------------------
# 11. Canonical tables and session information
# ------------------------------------------------------------------------------

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

fwrite(true_results, file.path(TABLE_DIR, "Figure5_true_mismatch_seedFE_TOST.csv"))
fwrite(used_results, file.path(TABLE_DIR, "Figure5_used_mismatch_seedFE.csv"))
fwrite(joint_results, file.path(TABLE_DIR, "Figure5_joint_true_used.csv"))
fwrite(odds_results, file.path(TABLE_DIR, "Figure5_odds_ratios.csv"))
fwrite(input_counts, file.path(TABLE_DIR, "Figure5_input_counts.csv"))
fwrite(export_log, file.path(OUTPUT_DIR, "Figure5_R_export_log.csv"))

writeLines(capture.output(sessionInfo()),
           file.path(OUTPUT_DIR, "Figure5_R_sessionInfo.txt"))

# ------------------------------------------------------------------------------
# 12. Console summary
# ------------------------------------------------------------------------------

cat("\nCANONICAL FIGURE 5 RESULTS\n")
cat("True-mismatch seed-fixed-effect models:\n")
print(true_results[, .(variant, beta, se, ci_low, ci_high, p_holm, equivalent_to_zero)])

cat("\nJoint true-versus-used models:\n")
print(joint_results[, .(variant, mismatch, beta, se, ci_low, ci_high, p_holm)])

cat("\nTrue-mismatch IQR odds ratios:\n")
print(odds_results[, .(variant, odds_ratio, or_low, or_high, percent_change_in_odds)])

cat("\nEXPORT STATUS\n")
print(export_log[, .N, by = .(format, status)])

cat("\nOutput directory: ", normalizePath(OUTPUT_DIR, mustWork = FALSE), "\n", sep = "")
