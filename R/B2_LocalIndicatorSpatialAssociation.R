#' @title Local Indicator of Spatial Association (LISA) Analysis
#'
#' @description
#' Conducts Local Moran's I spatial autocorrelation to identify exact locations
#' of spatial clusters and outliers. Generates two 1x2 4K high-resolution panels:
#' Image 1: Outbreak Mitigation (LISA Cluster Map & Significance Map).
#' Image 2: Statistical Diagnostic (Local Moran Scatterplot & Raw Magnitude Map).
#' Incorporates strict layout matrices for unified legends and robust outer margins.
#'
#' @param map_data An \code{sf} spatial data frame containing polygons.
#' @param var_col Character string. The column name containing the clinical variable.
#' @param sig_level Numeric. The p-value threshold for significance. Default is 0.05.
#' @param auto_colour Character string. Name of an RColorBrewer palette. Default is "YlOrRd".
#' @param reverse_palette Logical. If \code{TRUE}, reverses the color palette. Default is \code{FALSE}.
#' @param show_legend Logical. If \code{TRUE}, displays legends on the maps. Default is \code{TRUE}.
#' @param legend_orientation Character string. "vertical" or "horizontal". Default is "vertical".
#' @param common_legend Logical. If \code{TRUE}, uses a layout matrix to build a single unified legend. Default is \code{FALSE}.
#' @param footer Logical. If \code{TRUE}, adds a reproducible timestamp footer. Default is \code{TRUE}.
#' @param save_plot Logical. If \code{TRUE}, exports TWO plots in locked 4K resolution (3840x2160). Default is \code{FALSE}.
#' @param res Numeric. Resolution (DPI) for the exported plots. Default is 300.
#' @param save_prefix Character string. Prefix for the exported file names. Default is "B2".
#'
#' @return A list containing LISA statistics, Quadrant classifications, and p-values.
#' Generates two 4K visualizations in the active graphics device.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Execute LISA analysis with horizontal common legend using strict layout
#'   results <- B2_LocalIndicatorSpatialAssociation(
#'     map_data = nc,
#'     var_col = "SID74",
#'     sig_level = 0.05,
#'     legend_orientation = "horizontal",
#'     common_legend = TRUE,
#'     save_plot = FALSE
#'   )
#' }
B2_LocalIndicatorSpatialAssociation <- function(map_data,
                                                var_col,
                                                sig_level = 0.05,
                                                auto_colour = "YlOrRd",
                                                reverse_palette = FALSE,
                                                show_legend = TRUE,
                                                legend_orientation = "vertical",
                                                common_legend = FALSE,
                                                footer = TRUE,
                                                save_plot = FALSE,
                                                res = 300,
                                                save_prefix = "B2") {

  # 1. Strict Dependency & Input Validation
  if (!requireNamespace("sf", quietly = TRUE)) stop("CRITICAL ERROR: Package 'sf' is required.")
  if (!requireNamespace("spdep", quietly = TRUE)) stop("CRITICAL ERROR: Package 'spdep' is required for LISA.")
  stopifnot("var_col must be a character string" = is.character(var_col))
  if (!(var_col %in% names(map_data))) stop("var_col not found in map_data.")

  is_horiz <- if (legend_orientation == "horizontal") TRUE else FALSE

  # 2. Extract Data & NAs Subsetting
  var_val <- as.numeric(map_data[[var_col]])
  valid_idx <- stats::complete.cases(var_val)
  if (any(!valid_idx)) {
    warning("NAs detected. Subsetting spatial data to complete cases to protect Queen Contiguity.")
    map_data <- map_data[valid_idx, ]
    var_val <- var_val[valid_idx]
  }

  # 3. Spatial Weights & LISA Engine
  message("Building Connectivity Matrix & Calculating Local Moran's I...")
  nb <- spdep::poly2nb(map_data, queen = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  locm <- spdep::localmoran(var_val, lw, zero.policy = TRUE)
  Ii_val <- locm[, 1]  # Local Moran's I statistic
  p_val <- locm[, 5]   # Pr(z != 0) Two-sided p-value

  # 4. Standardized Lag & Quadrant Classification
  z_var <- as.numeric(scale(var_val))
  z_lag <- as.numeric(scale(spdep::lag.listw(lw, var_val, zero.policy = TRUE)))

  # Strict Clinical Boundary Check
  quadrant <- rep("Not Significant", length(var_val))
  quadrant[z_var > 0 & z_lag > 0 & p_val < sig_level] <- "High-High"
  quadrant[z_var < 0 & z_lag < 0 & p_val < sig_level] <- "Low-Low"
  quadrant[z_var > 0 & z_lag < 0 & p_val < sig_level] <- "High-Low"
  quadrant[z_var < 0 & z_lag > 0 & p_val < sig_level] <- "Low-High"

  quad_levels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant")
  quadrant_f <- factor(quadrant, levels = quad_levels)

  # Categorical Colors
  cat_cols <- c("High-High" = "#CB181D", "Low-Low" = "#2171B5", "High-Low" = "#FDBB84",
                "Low-High" = "#9ECAE1", "Not Significant" = "#F0F0F0")

  # 5. Significance Mapping Logic
  sig_cat <- rep("ns", length(p_val))
  sig_cat[p_val < 0.05] <- "p < 0.05"
  sig_cat[p_val < 0.01] <- "p < 0.01"
  sig_cat[p_val < 0.001] <- "p < 0.001"
  sig_levels <- c("p < 0.001", "p < 0.01", "p < 0.05", "ns")
  sig_cat_f <- factor(sig_cat, levels = sig_levels)
  sig_cols <- c("p < 0.001" = "#005A32", "p < 0.01" = "#41AB5D", "p < 0.05" = "#A1D99B", "ns" = "#F0F0F0")

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  rand_code <- paste0(sample(LETTERS, 4, replace = TRUE), collapse = "")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # =========================================================================
  # IMAGE 1: OUTBREAK MITIGATION (CLUSTER & SIGNIFICANCE)
  # =========================================================================
  if (common_legend) {
    graphics::layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4, 1))
    graphics::par(oma = c(7, 0, 3, 0), mar = c(2, 4, 3, 2) + 0.1, xpd = NA)
  } else {
    graphics::par(mfrow = c(1, 2), oma = c(7, 0, 3, 0), mar = c(7, 4, 3, 2) + 0.1, xpd = NA)
  }

  # Plot 1.1: LISA Cluster Map
  graphics::plot(sf::st_geometry(map_data), col = cat_cols[as.character(quadrant_f)],
                 border = "darkgrey", main = "LISA Cluster Map\n(Actionable Targeting)")
  if (show_legend && !common_legend) {
    graphics::legend("bottom", legend = names(cat_cols), fill = cat_cols,
                     title = "Spatial Quadrant", bty = "n", cex = 0.8, inset = c(0, -0.4), xpd = NA, horiz = is_horiz)
  }

  # Plot 1.2: LISA Significance Map
  graphics::plot(sf::st_geometry(map_data), col = sig_cols[as.character(sig_cat_f)],
                 border = "darkgrey", main = "LISA Significance Map\n(Statistical Priority)")
  if (show_legend && !common_legend) {
    graphics::legend("bottom", legend = names(sig_cols), fill = sig_cols,
                     title = "Significance Level", bty = "n", cex = 0.8, inset = c(0, -0.4), xpd = NA, horiz = is_horiz)
  }

  # Plot 1.3: Unified Legend Space (Only if common_legend = TRUE)
  if (common_legend) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    if (show_legend) {
      graphics::legend("center", legend = names(cat_cols), fill = cat_cols,
                       title = paste("LISA Cluster Zones (Threshold: p <", sig_level, ")"),
                       bty = "n", cex = 1.1, horiz = is_horiz, xpd = NA)
    }
  }

  # Top & Bottom Annotations (Outer Margins)
  graphics::mtext("EPIDEMIOLOGICAL MITIGATION TARGETS", side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  if (footer) {
    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
    graphics::mtext(footer_text, side = 1, line = 4, outer = TRUE, adj = 0.5, cex = 0.8, col = "black")
  }

  if (save_plot) {
    f1 <- paste0(save_prefix, "_OutbreakMitigation_", stamp, "_", rand_code, ".png")
    grDevices::dev.copy(grDevices::png, filename = f1, width = 3840, height = 2160, res = res)
    grDevices::dev.off()
    message(paste("SUCCESS: Image 1 (Outbreak Mitigation) exported ->", f1))
  }

  Sys.sleep(1) # Flush buffer

  # =========================================================================
  # IMAGE 2: STATISTICAL DIAGNOSTIC (SCATTERPLOT & MAGNITUDE MAP)
  # =========================================================================
  # Reset Par for Image 2
  if (common_legend) {
    graphics::layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4, 1))
    graphics::par(oma = c(7, 0, 3, 0), mar = c(4, 4, 3, 2) + 0.1, xpd = NA)
  } else {
    graphics::par(mfrow = c(1, 2), oma = c(7, 0, 3, 0), mar = c(7, 4, 3, 2) + 0.1, xpd = NA)
  }

  # Plot 2.1: Local Moran Scatterplot
  point_cols <- cat_cols[as.character(quadrant_f)]
  graphics::plot(z_var, z_lag, pch = 21, bg = point_cols, col = "darkgrey", cex = 1.5,
                 xlab = "Standardized Z-Score", ylab = "Spatial Lag Z-Score",
                 main = "Local Moran Scatterplot\n(Quadrant Distribution)")
  graphics::abline(h = 0, v = 0, lty = 2, col = "grey50")
  # (No legend usually needed for scatterplot as it matches the map)

  # Plot 2.2: Local Moran's I Magnitude Map
  n_colors <- 5
  if (!is.null(auto_colour)) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      mag_mapped <- grDevices::colorRampPalette(c("#FFF5F0", "#FB6A4A", "#CB181D", "#67000D"))(n_colors)
    } else {
      mag_mapped <- RColorBrewer::brewer.pal(n = n_colors, name = auto_colour)
    }
  } else {
    mag_mapped <- grDevices::colorRampPalette(c("#FFF5F0", "#FB6A4A", "#CB181D", "#67000D"))(n_colors)
  }
  if (reverse_palette) mag_mapped <- rev(mag_mapped)

  mag_breaks <- stats::quantile(Ii_val, probs = seq(0, 1, length.out = n_colors + 1), na.rm = TRUE)
  mag_cut <- cut(Ii_val, breaks = mag_breaks, include.lowest = TRUE)

  graphics::plot(sf::st_geometry(map_data), col = mag_mapped[as.integer(mag_cut)], border = "darkgrey",
                 main = "Local Moran's I Magnitude Map\n(Raw Autocorrelation Intensity)")

  if (show_legend && !common_legend) {
    graphics::legend("bottom", legend = levels(mag_cut), fill = mag_mapped,
                     title = "I_i Statistic", bty = "n", cex = 0.8, inset = c(0, -0.4), xpd = NA, horiz = is_horiz)
  }

  # Plot 2.3: Unified Legend Space (Only if common_legend = TRUE)
  if (common_legend) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    if (show_legend) {
      graphics::legend("center", legend = levels(mag_cut), fill = mag_mapped,
                       title = "Local Moran's I (I_i) Intensity Scale",
                       bty = "n", cex = 1.1, horiz = is_horiz, xpd = NA)
    }
  }

  # Top & Bottom Annotations (Outer Margins)
  graphics::mtext("LOCAL SPATIAL DIAGNOSTICS & MAGNITUDE", side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  if (footer) {
    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
    graphics::mtext(footer_text, side = 1, line = 4, outer = TRUE, adj = 0.5, cex = 0.8, col = "black")
  }

  if (save_plot) {
    f2 <- paste0(save_prefix, "_StatisticalDiagnostic_", stamp, "_", rand_code, ".png")
    grDevices::dev.copy(grDevices::png, filename = f2, width = 3840, height = 2160, res = res)
    grDevices::dev.off()
    message(paste("SUCCESS: Image 2 (Statistical Diagnostic) exported ->", f2))
  }

  # Silent Data Return
  invisible(data.frame(
    Local_I = Ii_val,
    p_value = p_val,
    LISA_Quadrant = quadrant_f
  ))
}
