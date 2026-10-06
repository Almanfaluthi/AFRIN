#' @title Epidemiological Choropleth Mapping and Empirical Bayes Smoothing
#' @description Generates high-resolution spatial disease maps. Supports raw disease
#' burden mapping for logistical planning, and Empirical Bayes Smoothing to correct
#' extreme variance (crude rate bias) in small-area populations. Automatically
#' coercing character data types to numeric to prevent shapefile DBF reading errors.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries.
#' @param cases_col Character. The column name representing the number of cases.
#' @param pop_col Character. The column name representing the population. Required if smoothing = TRUE.
#' @param smoothing Logical. If TRUE, applies Global Empirical Bayes smoothing and returns 4 diagnostic maps.
#' @param rate_multiplier Numeric. Multiplier for the rates (e.g., 100000 for per 100,000 population).
#' @param color_palette Character. Name of the RColorBrewer palette (default: "YlOrRd").
#' @param reverse_palette Logical. If TRUE, reverses the color palette vector.
#' @param legend Logical. If TRUE, displays the legend.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" for legend layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, saves the plot(s) as high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return The original 'sf' object appended with calculated epidemiological metrics
#' (Crude_Rate, Expected_Cases, SMR, EB_Smoothed_Rate) invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Mapping Malaria with Auto-Coercion on character columns
#' mapped_smooth <- A2_mapping(
#'   spatial_data = df_jateng,
#'   cases_col = "MALARIA",
#'   pop_col = "AREA",
#'   smoothing = TRUE,
#'   legend = TRUE,
#'   legend_position = "bottomleft",
#'   legend_orientation = "vertical",
#'   save_plot = TRUE,
#'   save_prefix = "Malaria_EB_Smoothed"
#' )
#' }
A2_mapping <- function(spatial_data,
                       cases_col,
                       pop_col = NULL,
                       smoothing = TRUE,
                       rate_multiplier = 100000,
                       color_palette = "YlOrRd",
                       reverse_palette = FALSE,
                       legend = TRUE,
                       legend_position = "bottomleft",
                       legend_orientation = "vertical",
                       footer = TRUE,
                       save_plot = FALSE,
                       save_prefix = "A2_mapping",
                       res = 300) {

  # 1. Strict Input Validation & Engine Check
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("AFRIN Error: Package 'sf' is required. Please install it.", call. = FALSE)
  }
  if (!inherits(spatial_data, "sf")) {
    stop("AFRIN Error: 'spatial_data' must be an 'sf' object.", call. = FALSE)
  }
  if (!(cases_col %in% names(spatial_data))) {
    stop(paste("AFRIN Error: Column", cases_col, "not found in spatial_data."), call. = FALSE)
  }

  # 2. Auto-Coercion Engine for Cases Column
  cases_raw <- spatial_data[[cases_col]]
  if (!is.numeric(cases_raw)) {
    warning(paste("AFRIN Warning: 'cases_col' (", cases_col, ") is not numeric. Auto-coercing to numeric..."), call. = FALSE)
    cases_num <- suppressWarnings(as.numeric(as.character(cases_raw)))

    na_diff <- sum(is.na(cases_num) & !is.na(cases_raw))
    if (na_diff > 0) {
      warning(paste("AFRIN Warning: Coercion introduced NAs in", na_diff, "rows (possibly string text). Treating them as 0 for epidemiological conservatism."), call. = FALSE)
    }
    cases_num[is.na(cases_num)] <- 0
    spatial_data[[cases_col]] <- cases_num
  } else {
    cases_raw[is.na(cases_raw)] <- 0
    spatial_data[[cases_col]] <- cases_raw
  }

  # 3. Epidemiological Math & Empirical Bayes Calculation
  if (smoothing) {
    if (is.null(pop_col) || !(pop_col %in% names(spatial_data))) {
      stop("AFRIN Error: 'pop_col' is required and must exist in spatial_data when smoothing = TRUE.", call. = FALSE)
    }

    # Auto-Coercion Engine for Population Column
    pop_raw <- spatial_data[[pop_col]]
    if (!is.numeric(pop_raw)) {
      warning(paste("AFRIN Warning: 'pop_col' (", pop_col, ") is not numeric. Auto-coercing to numeric..."), call. = FALSE)
      pop_num <- suppressWarnings(as.numeric(as.character(pop_raw)))

      na_diff_pop <- sum(is.na(pop_num) & !is.na(pop_raw))
      if (na_diff_pop > 0) {
        warning(paste("AFRIN Warning: Coercion introduced NAs in 'pop_col' in", na_diff_pop, "rows. Treating them as 1 to prevent division by zero."), call. = FALSE)
      }
      pop_num[is.na(pop_num) | pop_num <= 0] <- 1
      spatial_data[[pop_col]] <- pop_num
    } else {
      pop_raw[is.na(pop_raw) | pop_raw <= 0] <- 1
      spatial_data[[pop_col]] <- pop_raw
    }

    cases <- spatial_data[[cases_col]]
    pop <- spatial_data[[pop_col]]

    # Global metrics
    total_cases <- sum(cases, na.rm = TRUE)
    total_pop <- sum(pop, na.rm = TRUE)
    global_rate <- total_cases / total_pop

    # Crude Rate & Expected Cases
    crude_rate <- cases / pop
    expected_cases <- pop * global_rate
    smr <- ifelse(expected_cases == 0, 0, cases / expected_cases)

    # Global Empirical Bayes (Marshall 1991) using Base R stats
    variance_crude <- stats::var(crude_rate, na.rm = TRUE)
    mean_pop <- mean(pop, na.rm = TRUE)

    # Calculate prior variance (a)
    a <- variance_crude - (global_rate / mean_pop)
    if (is.na(a) || a < 0) a <- 0

    # Calculate shrinkage weight (w) and EB Smoothed Rate
    w <- a / (a + (global_rate / pop))
    eb_rate <- (w * crude_rate) + ((1 - w) * global_rate)

    # Append to sf object
    spatial_data[["Crude_Rate"]] <- crude_rate * rate_multiplier
    spatial_data[["Expected_Cases"]] <- expected_cases
    spatial_data[["SMR"]] <- smr
    spatial_data[["EB_Smoothed_Rate"]] <- eb_rate * rate_multiplier
  }

  # 4. High-Res Visualization Setup
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  # Color Palette Setup
  pal_colors <- NULL
  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    n_col <- RColorBrewer::brewer.pal.info[color_palette, "maxcolors"]
    pal_colors <- RColorBrewer::brewer.pal(n_col, color_palette)
  } else {
    pal_colors <- grDevices::heat.colors(7)
  }
  if (reverse_palette) {
    pal_colors <- rev(pal_colors)
  }
  pal_generator <- grDevices::colorRampPalette(pal_colors)

  # Internal function to generate choropleth map cleanly
  plot_choropleth <- function(var_name, title, suffix) {

    if (save_plot) {
      filename <- paste0(save_prefix, "_A2_mapping_", timestamp, "_", suffix, ".png")
      grDevices::png(filename = filename, width = 3840, height = 2160, res = res)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    # Outer margins strategy to prevent overlap with footer and legends
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = c(2, 2, 2, 2))

    vals <- spatial_data[[var_name]]

    # Quantile binning for choropleth
    if (length(unique(vals[!is.na(vals)])) > 1) {
      breaks <- unique(stats::quantile(vals, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
      if (length(breaks) > 1) {
        cats <- cut(vals, breaks = breaks, include.lowest = TRUE)
        colors_mapped <- pal_generator(length(levels(cats)))[as.numeric(cats)]
        legend_labels <- levels(cats)
        legend_colors <- pal_generator(length(levels(cats)))
      } else {
        colors_mapped <- rep(pal_colors[length(pal_colors)], length(vals))
        legend_labels <- "Homogeneous"
        legend_colors <- pal_colors[length(pal_colors)]
      }
    } else {
      colors_mapped <- rep(pal_colors[length(pal_colors)], length(vals))
      legend_labels <- "Homogeneous/Zero"
      legend_colors <- pal_colors[length(pal_colors)]
    }

    # Plot Base Geometry
    graphics::plot(sf::st_geometry(spatial_data), col = colors_mapped,
                   border = "#555555", lwd = 0.5)
    graphics::title(main = title, outer = TRUE, line = 0, cex.main = 1.5)

    # Legend Placement if TRUE
    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = legend_labels, fill = legend_colors,
                       title = var_name, bty = "n", horiz = horiz, cex = 0.9, xpd = NA)
    }

    # Strict Footer Placement (Bottom Center)
    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 4, adj = 0.5,
                      outer = TRUE, col = "darkgray", cex = 0.8)
    }

    if (save_plot) {
      grDevices::dev.off()
      message(paste("AFRIN Note: 4K Plot saved ->", paste0(save_prefix, "_A2_mapping_", timestamp, "_", suffix, ".png")))
    }

    graphics::par(old_par)
  }

  # 5. Generate 1 or 4 4K Plots based on smoothing parameter
  if (!smoothing) {
    # 1 Plot Only: Raw Volume
    plot_choropleth(cases_col, "Raw Disease Burden / Case Volume", "1_RawCases")
  } else {
    # 4 Panel Profiling for Spatial Epi
    plot_choropleth(cases_col, "Raw Disease Burden (Count)", "1_RawCases")

    rate_title <- paste0("Crude Rate (per ", format(rate_multiplier, scientific = FALSE), ")")
    plot_choropleth("Crude_Rate", rate_title, "2_CrudeRate")

    eb_title <- paste0("Empirical Bayes Smoothed Rate (per ", format(rate_multiplier, scientific = FALSE), ")")
    plot_choropleth("EB_Smoothed_Rate", eb_title, "3_EBSmoothedRate")

    plot_choropleth("SMR", "Standardized Morbidity/Mortality Ratio (SMR > 1.0 indicates excess risk)", "4_SMR")
  }

  message("AFRIN Note: Mapping complete. Variables auto-coerced (if needed) and appended successfully.")
  return(invisible(spatial_data))
}
