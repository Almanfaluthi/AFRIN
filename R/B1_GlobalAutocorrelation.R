#' @title Global Spatial Autocorrelation Analysis (Moran's I)
#'
#' @description
#' Conducts rigorous Global Spatial Autocorrelation testing using Moran's I
#' and Monte Carlo simulations. Generates two highly detailed 1x2 4K panels:
#' Plot 1: Moran Scatterplot and Monte Carlo Permutation Density Plot.
#' Plot 2: Spatial Connectivity (Weights) Network and Spatial Lag Choropleth Map.
#' Essential for identifying infectious disease clustering, poverty contagion,
#' or healthcare fraud (e.g., BPJS moral hazard clustering).
#'
#' @param map_data An \code{sf} spatial data frame containing polygons.
#' @param var_col Character string. The column name containing the clinical/variable of interest.
#' @param nsim Integer. Number of Monte Carlo permutations. Default is 999.
#' @param auto_colour Character string. Name of an RColorBrewer palette. Default is "YlOrRd". If \code{NULL}, uses manual color_palette.
#' @param color_palette Character vector of manual colors. Used if auto_colour is \code{NULL}.
#' @param reverse_palette Logical. If \code{TRUE}, reverses the color palette array. Default is \code{FALSE}.
#' @param show_legend Logical. If \code{TRUE}, displays a legend on the maps. Default is \code{TRUE}.
#' @param legend_pos Character string. Position of the legend. Default is "bottomright".
#' @param legend_orientation Character string. "vertical" or "horizontal". Default is "vertical".
#' @param common_legend Logical. If \code{TRUE}, centers a unified legend across the layout. Default is \code{FALSE}.
#' @param footer Logical. If \code{TRUE}, adds a reproducible timestamp footer. Default is \code{TRUE}.
#' @param save_plot Logical. If \code{TRUE}, exports TWO plots in locked 4K resolution (3840x2160). Default is \code{FALSE}.
#' @param res Numeric. Resolution (DPI) for the exported plots. Default is 300.
#' @param save_prefix Character string. Prefix for the exported file names. Default is "AFRIN".
#'
#' @return A list containing Moran's I statistics, p-values, and spatial lag values.
#' Generates two 4K visualizations in the active graphics device.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Analyzing spatial clustering with a horizontal legend
#'   results <- B1_GlobalAutocorrelation(
#'     map_data = nc,
#'     var_col = "SID74",
#'     nsim = 999,
#'     legend_orientation = "horizontal",
#'     common_legend = TRUE,
#'     save_plot = FALSE
#'   )
#' }
B1_GlobalAutocorrelation <- function(map_data,
                                           var_col,
                                           nsim = 999,
                                           auto_colour = "YlOrRd",
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

  # 1. Strict Dependency & Input Validation
  if (!requireNamespace("sf", quietly = TRUE)) stop("CRITICAL ERROR: Package 'sf' is required.")
  if (!requireNamespace("spdep", quietly = TRUE)) stop("CRITICAL ERROR: Package 'spdep' is required for Moran's I.")
  stopifnot("var_col must be a character string" = is.character(var_col))
  if (!(var_col %in% names(map_data))) stop("var_col not found in map_data.")

  is_horiz <- if (legend_orientation == "horizontal") TRUE else FALSE

  # 2. Extract Data and Remove NAs (Spatial weights fail with NAs)
  var_val <- as.numeric(map_data[[var_col]])
  valid_idx <- stats::complete.cases(var_val)
  if (any(!valid_idx)) {
    warning("NAs detected in var_col. Subsetting spatial data to complete cases only.")
    map_data <- map_data[valid_idx, ]
    var_val <- var_val[valid_idx]
  }

  # 3. Spatial Weights & Autocorrelation Engine
  message("Building Spatial Connectivity Matrix (Queen's Contiguity)...")
  nb <- spdep::poly2nb(map_data, queen = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  message("Executing Monte Carlo Moran's I Simulation...")
  set.seed(42) # Reproducibility for permutation
  mc_test <- spdep::moran.mc(var_val, lw, nsim = nsim, zero.policy = TRUE)

  moran_I <- mc_test$statistic
  p_val <- mc_test$p.value

  # 4. Strict Clinical Boundary Check & Impression Engine
  is_sig <- p_val < 0.05
  if (!is_sig) {
    impression <- "Impression: (ns) Random Spatial Distribution. No true clustering detected."
    col_sig <- "darkgrey"
  } else if (moran_I > 0) {
    impression <- "Impression: Significant Spatial Clustering (High-High / Low-Low) detected."
    col_sig <- "#CB181D" # Clinical Red
  } else {
    impression <- "Impression: Significant Spatial Dispersion (Checkerboard Pattern) detected."
    col_sig <- "#2171B5" # Clinical Blue
  }

  stat_title <- sprintf("Moran's I: %.3f | p-value: %.3f", moran_I, p_val)

  # 5. Calculate Spatial Lags for Scatterplot & Maps
  lag_val <- spdep::lag.listw(lw, var_val, zero.policy = TRUE)
  z_var <- scale(var_val)
  z_lag <- scale(lag_val)

  # Base Setup
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  rand_code <- paste0(sample(LETTERS, 4, replace = TRUE), collapse = "")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # =========================================================================
  # PLOT 1: MORAN SCATTERPLOT + MONTE CARLO DENSITY
  # =========================================================================
  graphics::par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1, oma = c(6, 0, 3, 0), xpd = TRUE)

  # Left: Moran Scatterplot
  graphics::plot(z_var, z_lag, pch = 21, bg = "lightgrey", col = "darkgrey",
                 xlab = "Standardized Value (Z-score)", ylab = "Spatial Lag (Z-score)",
                 main = "Moran Scatterplot\n(Neighborhood Influence)",
                 cex.main = 1.1, cex.lab = 1)
  graphics::abline(h = 0, v = 0, lty = 2, col = "grey50")
  fit <- stats::lm(z_lag ~ z_var)
  graphics::abline(fit, col = col_sig, lwd = 2)
  graphics::text(max(z_var, na.rm=TRUE), max(z_lag, na.rm=TRUE), "High-High", pos=2, col="red", cex=0.8)
  graphics::text(min(z_var, na.rm=TRUE), min(z_lag, na.rm=TRUE), "Low-Low", pos=4, col="blue", cex=0.8)

  # Right: Monte Carlo Density Plot
  dens <- stats::density(mc_test$res)
  graphics::plot(dens, main = "Monte Carlo Permutation\n(Null Hypothesis Distribution)",
                 xlab = "Simulated Moran's I Values", ylab = "Density", col = "black", lwd = 2)
  graphics::polygon(dens, col = "#F0F0F0", border = "black")
  graphics::abline(v = moran_I, col = col_sig, lwd = 3, lty = 1)
  graphics::text(moran_I, max(dens$y) * 0.9, "Observed\nStatistic", col = col_sig, pos = if(moran_I > 0) 2 else 4, cex=0.8)

  # Annotations & Strict Footer (Plot 1)
  graphics::mtext(paste("GLOBAL SPATIAL AUTOCORRELATION:", stat_title), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  graphics::mtext(impression, side = 3, line = -0.5, outer = TRUE, col = col_sig, font = 3, cex = 1)
  if (footer) {
    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
    graphics::mtext(footer_text, side = 1, line = 4, outer = TRUE, adj = 0.5, cex = 0.8, col = "black")
  }

  # Export Plot 1
  if (save_plot) {
    file_name1 <- paste0(save_prefix, "_B1_StatValidation_", stamp, "_", rand_code, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name1, width = 3840, height = 2160, res = res)
    grDevices::dev.off()
    message(paste("SUCCESS: Plot 1 (Statistics) exported ->", file_name1))
  }

  Sys.sleep(0.5) # Flush graphics buffer

  # =========================================================================
  # PLOT 2: CONNECTIVITY GRAPH + SPATIAL LAG MAP
  # =========================================================================
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 4, 2) + 0.1, oma = c(6, 0, 3, 0), xpd = TRUE)

  # Left: Connectivity Network
  graphics::plot(sf::st_geometry(map_data), border = "lightgrey", main = "Spatial Weights Network\n(Queen Contiguity)")
  coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(map_data)))
  spdep::plot.nb(nb, coords, add = TRUE, col = "#525252", lwd = 1.2, pch = 20, cex = 1.5)

  # Right: Spatial Lag Map Logic & RColorBrewer
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

  if (reverse_palette) colors_mapped <- rev(colors_mapped)

  lag_breaks <- stats::quantile(lag_val, probs = seq(0, 1, length.out = n_colors + 1), na.rm = TRUE)
  lag_cut <- cut(lag_val, breaks = lag_breaks, include.lowest = TRUE)
  color_idx_lag <- as.integer(lag_cut)

  graphics::plot(sf::st_geometry(map_data), col = colors_mapped[color_idx_lag], border = "darkgrey",
                 main = "Spatial Lag Map\n(Average of Neighbors)")

  # Dynamic Legend Engine
  if (show_legend) {
    if (common_legend) {
      # Unified legend placed centrally at the bottom of the entire 1x2 panel
      graphics::legend("bottom", legend = levels(lag_cut), fill = colors_mapped,
                       title = "Spatial Lag Rate (Quintiles)", bty = "n", cex = 1, horiz = is_horiz, xpd = NA, inset = -0.15)
    } else {
      # Specific legend placed purely inside the right-hand Lag map
      graphics::legend(legend_pos, legend = levels(lag_cut), fill = colors_mapped,
                       title = "Lag Variable", bty = "n", cex = 0.8, inset = 0.05, horiz = is_horiz)
    }
  }

  # Annotations & Strict Footer (Plot 2)
  graphics::mtext(paste("SPATIAL DISTRIBUTION TOPOLOGY:", stat_title), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  graphics::mtext(impression, side = 3, line = -0.5, outer = TRUE, col = col_sig, font = 3, cex = 1)
  if (footer) {
    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
    graphics::mtext(footer_text, side = 1, line = 4, outer = TRUE, adj = 0.5, cex = 0.8, col = "black")
  }

  # Export Plot 2
  if (save_plot) {
    file_name2 <- paste0(save_prefix, "_B1_TopologyMap_", stamp, "_", rand_code, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name2, width = 3840, height = 2160, res = res)
    grDevices::dev.off()
    message(paste("SUCCESS: Plot 2 (Topology) exported ->", file_name2))
  }

  # 6. Silent Return
  invisible(list(
    Moran_I = moran_I,
    p_value = p_val,
    Interpretation = impression,
    Raw_Values = var_val,
    Spatial_Lag = lag_val
  ))
}
