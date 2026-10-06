#' @title Spatial Choropleth Mapping with Empirical Bayes Smoothing and Unified Legends
#'
#' @description
#' Generates high-resolution 4K spatial choropleth maps for disease incidence.
#' Implements a dual-panel layout comparing Raw Incidence Rates against
#' Global Empirical Bayes (EB) Smoothed Rates. Features RColorBrewer integration,
#' dynamic legend orientations, and unified color scaling for epidemiological accuracy.
#'
#' @param map_data An \code{sf} spatial data frame containing polygons.
#' @param cases_col Character string. The column name containing disease event counts.
#' @param pop_col Character string. The column name containing population base.
#' @param smoothing Logical. If \code{TRUE}, applies Empirical Bayes smoothing. Default is \code{TRUE}.
#' @param auto_colour Character string. Name of an RColorBrewer palette (e.g., "YlOrRd"). If \code{NULL}, uses manual color_palette. Default is \code{NULL}.
#' @param color_palette Character vector of manual colors. Used if auto_colour is \code{NULL}. Default is a custom clinical red palette.
#' @param reverse_palette Logical. If \code{TRUE}, reverses the color palette array. Default is \code{FALSE}.
#' @param show_legend Logical. If \code{TRUE}, displays a legend on the maps. Default is \code{TRUE}.
#' @param legend_pos Character string. Position of the legend (e.g., "bottomright", "topleft", "bottom"). Default is "bottomright".
#' @param legend_orientation Character string. "vertical" or "horizontal". Default is "vertical".
#' @param common_legend Logical. If \code{TRUE} and smoothing is \code{TRUE}, applies a unified scale and single legend. Default is \code{FALSE}.
#' @param footer Logical. If \code{TRUE}, adds a reproducible timestamp footer. Default is \code{TRUE}.
#' @param save_plot Logical. If \code{TRUE}, exports the plot in locked 4K resolution (3840x2160). Default is \code{FALSE}.
#' @param res Numeric. Resolution (DPI) for the exported plot. Default is 300. Can be increased for journal submission.
#' @param save_prefix Character string. Prefix for the exported file name. Default is "AFRIN".
#'
#' @return A base R plot is rendered. If \code{save_plot = TRUE}, a 4K PNG file is written to the working directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   A2_mapping(
#'     map_data = nc,
#'     cases_col = "SID74",
#'     pop_col = "BIR74",
#'     smoothing = TRUE,
#'     auto_colour = "YlOrRd",
#'     show_legend = TRUE,
#'     legend_pos = "bottom",
#'     legend_orientation = "horizontal",
#'     common_legend = TRUE,
#'     save_plot = FALSE
#'   )
#' }
A2_mapping <- function(map_data,
                             cases_col,
                             pop_col,
                             smoothing = TRUE,
                             auto_colour = NULL,
                             color_palette = c("#FFF5F0", "#FB6A4A", "#CB181D", "#67000D"),
                             reverse_palette = FALSE,
                             show_legend = TRUE,
                             legend_pos = "bottomright",
                             legend_orientation = "vertical",
                             common_legend = FALSE,
                             footer = TRUE,
                             save_plot = FALSE,
                             res = 300,
                             save_prefix = "AFRIN") {

  # 1. Dependency & Input Validation
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("CRITICAL ERROR: Package 'sf' is required to plot spatial polygons.")
  }
  stopifnot("cases_col must be a character string" = is.character(cases_col))
  stopifnot("pop_col must be a character string" = is.character(pop_col))

  # Legend orientation logic
  is_horiz <- if (legend_orientation == "horizontal") TRUE else FALSE

  # 2. Extract Data
  O_i <- as.numeric(map_data[[cases_col]])
  N_i <- as.numeric(map_data[[pop_col]])

  if (any(N_i == 0, na.rm = TRUE)) {
    warning("Clinical Warning: Regions with zero population detected. Replaced with 1 to avoid Inf.")
    N_i[N_i == 0 & !is.na(N_i)] <- 1
  }

  # 3. Calculate Rates
  raw_rate <- (O_i / N_i) * 1000

  if (smoothing) {
    m <- sum(O_i, na.rm = TRUE) / sum(N_i, na.rm = TRUE)
    s2 <- stats::var(raw_rate / 1000, na.rm = TRUE)
    mean_N <- mean(N_i, na.rm = TRUE)

    a <- s2 - (m / mean_N)
    if (is.na(a) || a < 0) a <- 0

    w_i <- a / (a + (m / N_i))
    eb_rate <- (w_i * (raw_rate / 1000) + (1 - w_i) * m) * 1000
  } else {
    eb_rate <- raw_rate
  }

  # 4. Automated Color Palette Logic
  n_colors <- 5
  if (!is.null(auto_colour)) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      warning("Package 'RColorBrewer' is not installed. Falling back to manual color_palette.")
      colors_mapped <- grDevices::colorRampPalette(color_palette)(n_colors)
    } else {
      colors_mapped <- RColorBrewer::brewer.pal(n = n_colors, name = auto_colour)
    }
  } else {
    colors_mapped <- grDevices::colorRampPalette(color_palette)(n_colors)
  }

  if (reverse_palette) {
    colors_mapped <- rev(colors_mapped)
  }

  # 5. Break Intervals & Unified Scale Logic
  if (common_legend && smoothing) {
    # Combine data to create a unified global scale
    combined_rates <- c(raw_rate, eb_rate)
    global_breaks <- stats::quantile(combined_rates, probs = seq(0, 1, length.out = n_colors + 1), na.rm = TRUE)

    raw_cut <- cut(raw_rate, breaks = global_breaks, include.lowest = TRUE)
    eb_cut <- cut(eb_rate, breaks = global_breaks, include.lowest = TRUE)

    color_idx_raw <- as.integer(raw_cut)
    color_idx_eb <- as.integer(eb_cut)
    legend_labels <- levels(raw_cut)

  } else {
    # Individual scales
    raw_breaks <- stats::quantile(raw_rate, probs = seq(0, 1, length.out = n_colors + 1), na.rm = TRUE)
    raw_cut <- cut(raw_rate, breaks = raw_breaks, include.lowest = TRUE)
    color_idx_raw <- as.integer(raw_cut)
    raw_legend_labels <- levels(raw_cut)

    if (smoothing) {
      eb_breaks <- stats::quantile(eb_rate, probs = seq(0, 1, length.out = n_colors + 1), na.rm = TRUE)
      eb_cut <- cut(eb_rate, breaks = eb_breaks, include.lowest = TRUE)
      color_idx_eb <- as.integer(eb_cut)
      eb_legend_labels <- levels(eb_cut)
    }
  }

  # 6. Plotting Engine (High-Res 4K Ready with strict margins)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (smoothing) {
    # mar = margins of individual plots, oma = outer margins for the footer/legends
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 4, 2) + 0.1, oma = c(5, 0, 0, 0), xpd = TRUE)
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4, 4, 4, 2) + 0.1, oma = c(5, 0, 0, 0), xpd = TRUE)
  }

  # Plot Left: Raw Data
  graphics::plot(sf::st_geometry(map_data),
                 col = colors_mapped[color_idx_raw],
                 main = "Raw Incidence Rate\n(per 1,000 population)",
                 border = "darkgrey")

  # Draw separate legend for raw data if NOT using common legend
  if (show_legend && (!common_legend || !smoothing)) {
    graphics::legend(legend_pos, legend = raw_legend_labels, fill = colors_mapped,
                     title = "Raw Rate", bty = "n", cex = 0.8, inset = 0.05, horiz = is_horiz)
  }

  # Plot Right: Smoothed Data
  if (smoothing) {
    graphics::plot(sf::st_geometry(map_data),
                   col = colors_mapped[color_idx_eb],
                   main = "Empirical Bayes Smoothed Rate\n(Adjusted Risk)",
                   border = "darkgrey")

    # Draw either the common legend or the specific EB legend on the right plot
    if (show_legend) {
      if (common_legend) {
        graphics::legend(legend_pos, legend = legend_labels, fill = colors_mapped,
                         title = "Unified Rate (per 1,000)", bty = "n", cex = 0.8, inset = 0.05, horiz = is_horiz)
      } else {
        graphics::legend(legend_pos, legend = eb_legend_labels, fill = colors_mapped,
                         title = "EB Rate", bty = "n", cex = 0.8, inset = 0.05, horiz = is_horiz)
      }
    }
  }

  # 7. Reproducibility Footer (Strict Center Align)
  if (footer) {
    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
    graphics::mtext(footer_text, side = 1, line = 4, outer = TRUE, adj = 0.5, cex = 0.8, col = "black")
  }

  # 8. Locked 4K High-Resolution Export Module with Adjustable DPI
  if (save_plot) {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    rand_code <- paste0(sample(LETTERS, 4, replace = TRUE), collapse = "")
    file_name <- paste0(save_prefix, "_A2_mapping_", stamp, "_", rand_code, ".png")

    # Force 3840x2160, but inject the user-defined `res`
    grDevices::dev.copy(grDevices::png,
                        filename = file_name,
                        width = 3840,
                        height = 2160,
                        res = res)
    grDevices::dev.off()

    message(paste("SUCCESS:", file_name, "has been exported in 4K resolution at", res, "DPI."))
  }

  invisible(data.frame(Raw_Rate = raw_rate, EB_Smoothed_Rate = eb_rate))
}
