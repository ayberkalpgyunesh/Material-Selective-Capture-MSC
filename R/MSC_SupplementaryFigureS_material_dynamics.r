# ==============================================================================
# MSC Supplementary Figure S1 — Time-varying material identity
# Version 3: corrected plotting, no embedded panel letters
# ==============================================================================
# Input:  10_r_input from the reproducible material-dynamics notebook
#
# Usage:  Rscript MSC_SupplementaryFigureS1_material_dynamics_v3_corrected.r \
#                  /path/to/10_r_input
#
# Required files:
#   panel_A_granule_material_drift.csv
#   panel_C_timestep_control.csv
#   panel_D_copy_number_enrichment.csv
#   panel_F_psd_target_drift.csv
#   msc_drift_reference_controls_v7.csv
#   msc_drift_paired_v7.csv
#   summary_material_dynamics_v7.json
#   tau_enc_per_seed.json
#
# This script does not simulate. It validates and summarizes the completed
# seed-level outputs, exports panels A–F and an assembled figure, and writes the
# exact plotted data. No A–F letters are embedded in the plots; add them later
# in BioRender.
#
# Corrections in version 3:
#   1. sigma_rel labels use as.character(), preventing 0.5 and 1.0 from becoming
#      NA after factor conversion.
#   2. Panel A suppresses the duplicate fill legend.
#   3. Panel E explicitly fixes the copy-number order to 1, 3, 10.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Install and load packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table", "ggplot2", "scales", "jsonlite", "patchwork"
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
  library(jsonlite)
  library(patchwork)
})

have_ragg <- requireNamespace("ragg", quietly = TRUE)
have_svglite <- requireNamespace("svglite", quietly = TRUE)


# ------------------------------------------------------------------------------
# 2. Paths and input discovery
# ------------------------------------------------------------------------------

arguments <- commandArgs(trailingOnly = TRUE)

INPUT_DIR <- if (length(arguments) >= 1 && nzchar(arguments[1])) {
  arguments[1]
} else {
  getwd()
}

OUTPUT_DIR <- file.path(INPUT_DIR, "R_output_SupplementaryFigureS1_v3")

for (folder in c("pdf", "png", "svg", "tiff", "tables")) {
  dir.create(file.path(OUTPUT_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}

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


find_input <- function(filename) {

  matches <- all_files[basename(all_files) == filename]

  if (length(matches) == 0) {
    stop("Could not find required input under ", INPUT_DIR, ":\n  ",
         filename, call. = FALSE)
  }

  if (length(matches) > 1) {
    matches <- matches[which.max(file.info(matches)$mtime)]
  }

  message("Reading: ", matches[[1]])
  matches[[1]]
}


drift_file <- find_input("panel_A_granule_material_drift.csv")
dt_file <- find_input("panel_C_timestep_control.csv")
copy_file <- find_input("panel_D_copy_number_enrichment.csv")
psd_file <- find_input("panel_F_psd_target_drift.csv")
reference_file <- find_input("msc_drift_reference_controls_v7.csv")
paired_file <- find_input("msc_drift_paired_v7.csv")
summary_file <- find_input("summary_material_dynamics_v7.json")
tau_enc_file <- find_input("tau_enc_per_seed.json")


# ------------------------------------------------------------------------------
# 3. Frozen design and visual definitions
# ------------------------------------------------------------------------------

expected_rho <- c(
  0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 30, 100, 300
)

expected_sigma_rel <- c(0.25, 0.50, 1.00)
expected_dt_rho <- c(0.03, 0.10, 0.30, 1.00)
expected_copies <- c(1L, 3L, 10L)
expected_kc <- c(0, 1e-5, 1e-4, 1e-3, 1e-2)

DPI <- 600
N_BOOT <- 10000L
BOOTSTRAP_SEED <- 20260911L

sigma_colours <- c(
  "0.25" = "#3C78A8",
  "0.5"  = "#3FA35C",
  "1"    = "#D55E00"
)

metric_colours <- c(
  "Docking selectivity" = "#1B9E77",
  "Material-memory specificity" = "#7570B3"
)

dt_colours <- c(
  "Original timestep" = "#4C78A8",
  "Half timestep" = "#E45756"
)

copy_colour <- "#3FAE5A"
psd_colour <- "#7570B3"
static_colour <- "#666666"
random_floor_colour <- "#B2182B"


# ------------------------------------------------------------------------------
# 4. Load and validate completed outputs
# ------------------------------------------------------------------------------

drift <- fread(drift_file)
dt_control <- fread(dt_file)
copy_number <- fread(copy_file)
psd_drift <- fread(psd_file)
reference_controls <- fread(reference_file)
paired_python <- fread(paired_file)
summary_json <- jsonlite::fromJSON(summary_file, simplifyVector = TRUE)
tau_enc <- jsonlite::fromJSON(tau_enc_file, simplifyVector = TRUE)

static_drift <- drift[is.infinite(rho)]
dynamic_drift <- drift[is.finite(rho)]

if (!isTRUE(all.equal(sort(unique(dynamic_drift$rho)), expected_rho))) {
  stop("The finite rho grid does not match the frozen v7 design.")
}

if (!isTRUE(all.equal(sort(unique(dynamic_drift$sigma_rel)),
                      expected_sigma_rel))) {
  stop("The sigma_rel grid does not match the frozen v7 design.")
}

if (nrow(drift) != 4300L || nrow(static_drift) != 100L) {
  stop("Expected 4,300 drift rows, including 100 static references.")
}

if (!isTRUE(all.equal(sort(unique(copy_number$n_copies)), expected_copies))) {
  stop("Copy-number grid does not match 1, 3 and 10.")
}

if (!isTRUE(all.equal(sort(unique(psd_drift$k_c)), expected_kc))) {
  stop("PSD-target-drift grid does not match the frozen design.")
}

if (dynamic_drift[, uniqueN(seed), by = .(rho, sigma_rel)][, any(V1 != 100L)]) {
  stop("At least one drift condition does not contain 100 seeds.")
}

if (copy_number[, uniqueN(seed), by = n_copies][, any(V1 != 100L)]) {
  stop("Every copy-number condition must contain 100 seeds.")
}

if (psd_drift[, uniqueN(seed), by = k_c][, any(V1 != 50L)]) {
  stop("Every PSD-target-drift condition must contain 50 seeds.")
}

if (anyDuplicated(drift) || anyDuplicated(dt_control) ||
    anyDuplicated(copy_number) || anyDuplicated(psd_drift)) {
  stop("Duplicated rows detected in an input table.")
}

if (!all(is.finite(dynamic_drift$continuous_specificity)) ||
    !all(is.finite(dynamic_drift$material_memory_specificity))) {
  stop("Non-finite primary specificity values detected.")
}

n_undefined_enrichment <- sum(!is.finite(drift$active_strong_enrichment))
undefined_enrichment_seeds <- unique(
  drift[!is.finite(active_strong_enrichment), seed]
)

message("Input checks passed.")
message("Undefined enrichment rows: ", n_undefined_enrichment,
        "; affected seed(s): ",
        paste(undefined_enrichment_seeds, collapse = ", "))


# ------------------------------------------------------------------------------
# 5. Seed-bootstrap helper
# ------------------------------------------------------------------------------

bootstrap_mean_ci <- function(values, seed = BOOTSTRAP_SEED) {

  values <- values[is.finite(values)]

  if (length(values) < 2) {
    return(list(
      estimate = if (length(values)) mean(values) else NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      n = length(values)
    ))
  }

  set.seed(seed)

  bootstrap_means <- replicate(
    N_BOOT,
    mean(sample(values, size = length(values), replace = TRUE))
  )

  list(
    estimate = mean(values),
    ci_low = unname(quantile(bootstrap_means, 0.025)),
    ci_high = unname(quantile(bootstrap_means, 0.975)),
    n = length(values)
  )
}


summarise_bootstrap <- function(table, value_column, group_columns,
                                seed_offset = 0L) {

  table[
    ,
    {
      result <- bootstrap_mean_ci(
        get(value_column),
        seed = BOOTSTRAP_SEED + seed_offset + .GRP
      )

      .(
        estimate = result$estimate,
        ci_low = result$ci_low,
        ci_high = result$ci_high,
        n = result$n
      )
    },
    by = group_columns
  ]
}


# ------------------------------------------------------------------------------
# 6. Panel A data
# ------------------------------------------------------------------------------
# Correction: use as.character() rather than format(). Vector-wide formatting
# had produced "0.50" and "1.00", which did not match factor levels "0.5" and
# "1" and were therefore converted to NA.

panel_A_data <- copy(paired_python)

panel_A_data[, sigma_label := factor(
  as.character(sigma_rel),
  levels = c("0.25", "0.5", "1")
)]

if (anyNA(panel_A_data$sigma_label)) {
  stop("sigma_rel factor conversion produced NA values.")
}

panel_A_data[, unresolved := !as.logical(dt_resolved)]

random_floor <- as.numeric(summary_json$floor_relative)

rho_dip_high <- if (!is.null(summary_json$dip[["1.0"]])) {
  as.numeric(summary_json$dip[["1.0"]][["rho_dip"]])
} else {
  as.numeric(summary_json$dip[["1"]][["rho_dip"]])
}


# ------------------------------------------------------------------------------
# 7. Panel B data
# ------------------------------------------------------------------------------

finite_enrichment_counts <- dynamic_drift[
  ,
  .(n = sum(is.finite(active_strong_enrichment))),
  by = .(sigma_rel, rho)
]

panel_B_docking <- merge(
  panel_A_data[
    ,
    .(
      sigma_rel,
      sigma_label,
      rho,
      metric = "Docking selectivity",
      estimate = strongenr_rel
    )
  ],
  finite_enrichment_counts,
  by = c("sigma_rel", "rho"),
  all.x = TRUE
)

panel_B_memory <- panel_A_data[
  ,
  .(
    sigma_rel,
    sigma_label,
    rho,
    metric = "Material-memory specificity",
    estimate = Mspec_rel,
    n = n
  )
]

panel_B_data <- rbindlist(
  list(panel_B_docking, panel_B_memory),
  use.names = TRUE,
  fill = TRUE
)

panel_B_data[, metric := factor(
  metric,
  levels = c("Docking selectivity", "Material-memory specificity")
)]

if (anyNA(panel_B_data$sigma_label)) {
  stop("Panel B contains an undefined sigma_rel facet.")
}


# ------------------------------------------------------------------------------
# 8. Panel C data
# ------------------------------------------------------------------------------

static_dt <- dt_control[
  is.infinite(rho),
  .(seed, dt_factor, static_specificity = continuous_specificity)
]

dynamic_dt <- merge(
  dt_control[is.finite(rho)],
  static_dt,
  by = c("seed", "dt_factor"),
  all.x = TRUE
)

dynamic_dt[, specificity_relative :=
             continuous_specificity / static_specificity]

dynamic_dt[, timestep := factor(
  fifelse(dt_factor == 1, "Original timestep", "Half timestep"),
  levels = c("Half timestep", "Original timestep")
)]

panel_C_data <- dynamic_dt[
  ,
  .(
    estimate = mean(specificity_relative, na.rm = TRUE),
    sem = sd(specificity_relative, na.rm = TRUE) /
      sqrt(sum(is.finite(specificity_relative))),
    n = sum(is.finite(specificity_relative))
  ),
  by = .(rho, dt_factor, timestep)
]

panel_C_data[, ci_low := estimate - 1.96 * sem]
panel_C_data[, ci_high := estimate + 1.96 * sem]

# Paired numerical audit: half timestep minus original timestep.

dt_wide <- dcast(
  dynamic_dt,
  seed + rho ~ dt_factor,
  value.var = "specificity_relative"
)

setnames(dt_wide, c("0.5", "1"), c("half_dt", "original_dt"))
dt_wide[, half_minus_original := half_dt - original_dt]

panel_C_paired <- summarise_bootstrap(
  dt_wide,
  value_column = "half_minus_original",
  group_columns = "rho",
  seed_offset = 1000L
)


# ------------------------------------------------------------------------------
# 9. Panels D and E data
# ------------------------------------------------------------------------------

panel_D_data <- summarise_bootstrap(
  copy_number,
  value_column = "active_strong_enrichment",
  group_columns = "n_copies",
  seed_offset = 2000L
)

copy_specificity <- dcast(
  copy_number,
  seed ~ n_copies,
  value.var = "continuous_specificity"
)

setnames(
  copy_specificity,
  c("1", "3", "10"),
  c("specificity_N1", "specificity_N3", "specificity_N10")
)

panel_E_data <- rbindlist(list(
  copy_specificity[
    , .(seed, n_copies = 1L,
        specificity_relative = specificity_N1 / specificity_N1)
  ],
  copy_specificity[
    , .(seed, n_copies = 3L,
        specificity_relative = specificity_N3 / specificity_N1)
  ],
  copy_specificity[
    , .(seed, n_copies = 10L,
        specificity_relative = specificity_N10 / specificity_N1)
  ]
))

panel_E_data[, copy_label := factor(
  as.character(n_copies),
  levels = c("1", "3", "10")
)]

panel_E_summary <- summarise_bootstrap(
  panel_E_data,
  value_column = "specificity_relative",
  group_columns = "n_copies",
  seed_offset = 3000L
)


# ------------------------------------------------------------------------------
# 10. Panel F data
# ------------------------------------------------------------------------------

psd_wide <- dcast(
  psd_drift,
  seed ~ k_c,
  value.var = "continuous_specificity"
)

kc_names <- setdiff(names(psd_wide), "seed")

panel_F_data <- melt(
  psd_wide,
  id.vars = "seed",
  measure.vars = kc_names,
  variable.name = "k_c_character",
  value.name = "continuous_specificity"
)

panel_F_data[, k_c := as.numeric(as.character(k_c_character))]
panel_F_data[, static_specificity :=
               continuous_specificity[k_c == 0], by = seed]
panel_F_data[, specificity_relative :=
               continuous_specificity / static_specificity]

panel_F_summary <- summarise_bootstrap(
  panel_F_data,
  value_column = "specificity_relative",
  group_columns = "k_c",
  seed_offset = 4000L
)

panel_F_summary[, x_position := match(k_c, expected_kc) - 1L]

kc_labels <- c(
  "0",
  expression(10^-5),
  expression(10^-4),
  expression(10^-3),
  expression(10^-2)
)


# ------------------------------------------------------------------------------
# 11. Plotting theme
# ------------------------------------------------------------------------------
# No panel letters are embedded. Panel identity remains in filenames only.

theme_msc <- function(base_size = 9) {

  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.40),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = element_text(size = base_size + 0.2),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "plain"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      plot.margin = margin(6, 7, 5, 6)
    )
}

rho_breaks <- c(0.003, 0.01, 0.1, 1, 10, 100, 300)
rho_labels <- c("0.003", "0.01", "0.1", "1", "10", "100", "300")


# ------------------------------------------------------------------------------
# 12. Corrected panels A–F
# ------------------------------------------------------------------------------

panel_A <- ggplot(
  panel_A_data,
  aes(x = rho, y = spec_rel_mean, colour = sigma_label,
      fill = sigma_label, group = sigma_label)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.50,
    colour = static_colour
  ) +
  geom_hline(
    yintercept = random_floor,
    linetype = "dotted",
    linewidth = 0.60,
    colour = random_floor_colour
  ) +
  geom_vline(
    xintercept = rho_dip_high,
    linetype = "dotdash",
    linewidth = 0.50,
    colour = sigma_colours[["1"]]
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    alpha = 0.16,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.70) +
  geom_point(size = 2.0) +
  geom_point(
    data = panel_A_data[unresolved == TRUE],
    aes(x = rho, y = spec_rel_mean),
    inherit.aes = FALSE,
    shape = 4,
    size = 2.5,
    stroke = 0.8,
    colour = "black"
  ) +
  scale_x_log10(
    breaks = rho_breaks,
    labels = rho_labels
  ) +
  scale_colour_manual(values = sigma_colours) +
  scale_fill_manual(values = sigma_colours) +
  guides(fill = "none") +
  coord_cartesian(ylim = c(0.52, 1.32)) +
  labs(
    x = expression(rho == tau[m] / tau[enc]),
    y = "Relative consolidation specificity",
    colour = expression(sigma[rel])
  ) +
  theme_msc()


panel_B <- ggplot(
  panel_B_data,
  aes(x = rho, y = estimate, colour = metric,
      shape = metric, group = metric)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.50,
    colour = static_colour
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ sigma_label,
    nrow = 1,
    drop = FALSE,
    labeller = labeller(
      sigma_label = c(
        "0.25" = "sigma_rel = 0.25",
        "0.5" = "sigma_rel = 0.5",
        "1" = "sigma_rel = 1"
      )
    )
  ) +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  scale_colour_manual(values = metric_colours) +
  scale_shape_manual(values = c(
    "Docking selectivity" = 15,
    "Material-memory specificity" = 17
  )) +
  coord_cartesian(ylim = c(0.81, 1.015)) +
  labs(
    x = expression(rho == tau[m] / tau[enc]),
    y = "Relative to static MSC"
  ) +
  theme_msc(base_size = 8.5) +
  theme(
    axis.text.x = element_text(size = 7.3),
    panel.spacing = grid::unit(2.5, "mm")
  )


panel_C <- ggplot(
  panel_C_data,
  aes(x = rho, y = estimate, colour = timestep,
      shape = timestep, group = timestep)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.50,
    colour = static_colour
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.025,
    linewidth = 0.55
  ) +
  geom_line(linewidth = 0.70) +
  geom_point(size = 2.4) +
  scale_x_log10(
    breaks = expected_dt_rho,
    labels = c("0.03", "0.1", "0.3", "1")
  ) +
  scale_colour_manual(values = dt_colours) +
  scale_shape_manual(values = c(
    "Half timestep" = 15,
    "Original timestep" = 16
  )) +
  labs(
    x = expression(rho == tau[m] / tau[enc]),
    y = "Relative consolidation specificity"
  ) +
  theme_msc()


panel_D <- ggplot(
  panel_D_data,
  aes(x = factor(n_copies, levels = c(1, 3, 10)), y = estimate)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    linewidth = 0.55,
    colour = random_floor_colour
  ) +
  geom_col(
    width = 0.65,
    fill = copy_colour,
    colour = "black",
    linewidth = 0.45
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.18,
    linewidth = 0.60
  ) +
  scale_x_discrete(limits = c("1", "3", "10"), drop = FALSE) +
  coord_cartesian(ylim = c(0, 1.48)) +
  labs(
    x = "Transcript copies per granule",
    y = "Strong-granule enrichment"
  ) +
  theme_msc() +
  theme(legend.position = "none")


panel_E <- ggplot(
  panel_E_data,
  aes(x = copy_label, y = specificity_relative)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    linewidth = 0.55,
    colour = random_floor_colour
  ) +
  geom_violin(
    data = panel_E_data[n_copies != 1],
    aes(group = copy_label),
    width = 0.82,
    fill = "#99CC99",
    colour = "#3C5B45",
    linewidth = 0.45,
    trim = FALSE
  ) +
  geom_boxplot(
    data = panel_E_data[n_copies != 1],
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.45
  ) +
  geom_point(
    data = panel_E_summary,
    aes(x = factor(as.character(n_copies), levels = c("1", "3", "10")),
        y = estimate),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.8,
    fill = "black",
    colour = "white",
    stroke = 0.5
  ) +
  scale_x_discrete(
    limits = c("1", "3", "10"),
    drop = FALSE
  ) +
  labs(
    x = "Transcript copies per granule",
    y = "Specificity relative to N = 1"
  ) +
  theme_msc() +
  theme(legend.position = "none")


panel_F <- ggplot(
  panel_F_summary,
  aes(x = x_position, y = estimate, group = 1)
) +
  annotate(
    "rect",
    xmin = -0.35,
    xmax = 2.35,
    ymin = -Inf,
    ymax = Inf,
    fill = "#DFF2E5",
    alpha = 0.55
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.50,
    colour = static_colour
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.10,
    linewidth = 0.60,
    colour = psd_colour
  ) +
  geom_line(
    linewidth = 0.70,
    colour = psd_colour
  ) +
  geom_point(
    size = 2.5,
    colour = psd_colour
  ) +
  scale_x_continuous(
    breaks = 0:4,
    labels = kc_labels,
    limits = c(-0.4, 4.4)
  ) +
  labs(
    x = expression(k[c]~(s^{-1})),
    y = expression("Specificity relative to "*k[c]*" = 0")
  ) +
  theme_msc() +
  theme(legend.position = "none")


# ------------------------------------------------------------------------------
# 13. Assemble without panel letters
# ------------------------------------------------------------------------------

combined_figure <- (
  (panel_A | panel_B) /
  (panel_D | panel_E) /
  (panel_C | panel_F)
) +
  plot_layout(heights = c(1, 1, 1))


# ------------------------------------------------------------------------------
# 14. Export helper
# ------------------------------------------------------------------------------

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

save_panel(panel_A, "FigureS1A_granule_material_drift_corrected", 95, 78)
save_panel(panel_B, "FigureS1B_docking_memory_dissociation_corrected", 145, 78)
save_panel(panel_C, "FigureS1C_timestep_control", 95, 78)
save_panel(panel_D, "FigureS1D_copy_number_enrichment", 80, 78)
save_panel(panel_E, "FigureS1E_copy_number_specificity_corrected", 85, 78)
save_panel(panel_F, "FigureS1F_psd_target_drift", 90, 78)
save_panel(
  combined_figure,
  "SupplementaryFigureS1_material_dynamics_corrected_no_panel_labels",
  190,
  225
)


# ------------------------------------------------------------------------------
# 15. Export exact plotted data and audit information
# ------------------------------------------------------------------------------

table_dir <- file.path(OUTPUT_DIR, "tables")

fwrite(panel_A_data, file.path(table_dir, "FigureS1A_plotted_data.csv"))
fwrite(panel_B_data, file.path(table_dir, "FigureS1B_plotted_data.csv"))
fwrite(panel_C_data, file.path(table_dir, "FigureS1C_plotted_data.csv"))
fwrite(panel_C_paired,
       file.path(table_dir, "FigureS1C_paired_dt_differences.csv"))
fwrite(panel_D_data, file.path(table_dir, "FigureS1D_plotted_data.csv"))
fwrite(panel_E_data, file.path(table_dir, "FigureS1E_seed_level_data.csv"))
fwrite(panel_E_summary, file.path(table_dir, "FigureS1E_plotted_summary.csv"))
fwrite(panel_F_data, file.path(table_dir, "FigureS1F_seed_level_data.csv"))
fwrite(panel_F_summary, file.path(table_dir, "FigureS1F_plotted_summary.csv"))
fwrite(export_log, file.path(table_dir, "FigureS1_export_log.csv"))

undefined_audit <- data.table(
  n_undefined_enrichment_rows = n_undefined_enrichment,
  affected_seeds = paste(undefined_enrichment_seeds, collapse = ";"),
  interpretation = paste(
    "active_strong_enrichment is undefined when a simulation contains",
    "zero active-synapse docking events"
  )
)

fwrite(
  undefined_audit,
  file.path(table_dir, "FigureS1_undefined_enrichment_audit.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "FigureS1_session_info.txt")
)


# ------------------------------------------------------------------------------
# 16. Final visual and numerical checks
# ------------------------------------------------------------------------------

expected_sigma_labels <- c("0.25", "0.5", "1")
observed_sigma_labels <- levels(droplevels(panel_A_data$sigma_label))

if (!identical(observed_sigma_labels, expected_sigma_labels)) {
  stop("Corrected sigma_rel labels were not recovered.")
}

observed_copy_order <- levels(panel_E_data$copy_label)

if (!identical(observed_copy_order, c("1", "3", "10"))) {
  stop("Panel E copy-number order is not 1, 3, 10.")
}

cat("\n================ MSC Supplementary Figure S1 ================\n")
cat("Correct sigma_rel facets:", paste(observed_sigma_labels, collapse = ", "), "\n")
cat("Correct copy-number order:", paste(observed_copy_order, collapse = ", "), "\n")
cat("No panel letters were embedded in the plots.\n")
cat("Undefined enrichment rows:", n_undefined_enrichment,
    "| affected seed(s):", paste(undefined_enrichment_seeds, collapse = ", "), "\n")

cat("\nExport results:\n")
print(export_log[, .N, by = .(format, status)])

failed_exports <- export_log[status == "failed"]

if (nrow(failed_exports) > 0) {
  cat("\nSome exports failed:\n")
  print(failed_exports)
} else {
  cat("\nAll corrected panels exported in all four formats.\n")
}

cat("\nOutput folder: ", OUTPUT_DIR, "\n", sep = "")
cat("==============================================================\n")
