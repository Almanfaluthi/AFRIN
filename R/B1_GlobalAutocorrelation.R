#' @title Global Spatial Autocorrelation (Moran's I) and Diagnostic Mapping
#' @description Computes Global Moran's I using pure Base R matrix algebra to detect
#' spatial clustering (High-High/Low-Low) or dispersion. It generates a 4-stage
#' diagnostic profile including Spatial Connectivity Network, Moran Scatterplot,
#' Monte Carlo Permutation Density, and Spatial Lag Map. It uses Queen contiguity
#' weights directly derived from 'sf' geometries.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries.
#' @param var_col Character. The column name of the variable to analyze.
#' @param nsim Numeric. Number of Monte Carlo permutations for p-value (default: 999).
#' @param color_palette Character. RColorBrewer palette for the Spatial Lag Map (default: "YlOrRd").
#' @param reverse_palette Logical. If TRUE, reverses the color palette vector.
#' @param legend Logical. If TRUE, displays legends on applicable plots.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" for legend layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, saves the 4 diagnostic plots as high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return A list containing the Moran's I statistic, p-value, and the original
#' 'sf' object appended with Spatial_Lag, Z_Score, and Moran_Quadrant invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Analyzing spatial clustering of BPJS Claims / Stunting
#' moran_results <- B1_GlobalAutocorrelationMoranI(
#'   spatial_data = df_jateng,
#'   var_col = "MALARIA",
#'   nsim = 999,
#'   save_plot = TRUE,
#'   save_prefix = "B1_Morans'I"
#' )
#'
#' # Print Empirical Moran's I and p-value
#' print(moran_results$Morans_I)
#' print(moran_results$p_value)
#' }
B1_GlobalAutocorrelation <- function(spatial_data,
                                           var_col,
                                           nsim = 999,
                                           color_palette = "YlOrRd",
                                           reverse_palette = FALSE,
                                           legend = TRUE,
                                           legend_position = "bottomleft",
                                           legend_orientation = "vertical",
                                           footer = TRUE,
                                           save_plot = FALSE,
                                           save_prefix = "B1_Morans'I",
                                           res = 300) {

  # 1. Strict Input Validation & Engine Check
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("AFRIN Error: Package 'sf' is required for geometry operations.", call. = FALSE)
  }
  if (!inherits(spatial_data, "sf")) {
    stop("AFRIN Error: 'spatial_data' must be an 'sf' object.", call. = FALSE)
  }
  if (!(var_col %in% names(spatial_data))) {
    stop(paste("AFRIN Error: Column", var_col, "not found in spatial_data."), call. = FALSE)
  }

  # 2. Auto-Coercion Engine (Character to Numeric Handling)
  raw_var <- spatial_data[[var_col]]
  if (!is.numeric(raw_var)) {
    warning(paste("AFRIN Warning: 'var_col' (", var_col, ") is not numeric. Auto-coercing to numeric..."), call. = FALSE)
    num_var <- suppressWarnings(as.numeric(as.character(raw_var)))

    na_diff <- sum(is.na(num_var) & !is.na(raw_var))
    if (na_diff > 0) {
      warning(paste("AFRIN Warning: Coercion introduced NAs in", na_diff, "rows. Treating them as 0 to maintain spatial matrix dimensions."), call. = FALSE)
    }
    num_var[is.na(num_var)] <- 0
    spatial_data[[var_col]] <- num_var
  } else {
    raw_var[is.na(raw_var)] <- 0
    spatial_data[[var_col]] <- raw_var
  }

  x <- spatial_data[[var_col]]
  n <- length(x)

  # 3. Pure Base R Spatial Weights Matrix (Queen Contiguity)
  message("AFRIN Note: Building Spatial Weights Matrix (Queen) without external dependencies...")
  touch_list <- sf::st_touches(spatial_data, sparse = FALSE)
  W_bin <- matrix(as.numeric(touch_list), nrow = n, ncol = n)

  # Row standardization
  row_sums <- rowSums(W_bin)
  # Prevent division by zero for spatial islands
  row_sums_safe <- ifelse(row_sums == 0, 1, row_sums)
  W_std <- sweep(W_bin, 1, row_sums_safe, FUN = "/")

  # 4. Global Moran's I Calculation (Matrix Algebra)
  x_bar <- mean(x, na.rm = TRUE)
  z <- x - x_bar
  # Spatial Lag of standardized values
  lag_z <- as.numeric(W_std %*% z)
  # Spatial Lag of raw values (for mapping)
  lag_x <- as.numeric(W_std %*% x)

  S0 <- sum(W_std)
  var_z <- sum(z^2)

  # Empirical Moran's I formula: (n/S0) * (z' W z) / (z' z)
  W_zz <- sum(W_std * (z %*% t(z)))
  I_emp <- (n / S0) * (W_zz / var_z)

  # 5. Monte Carlo Permutations (Simulated p-value)
  message(paste("AFRIN Note: Running", nsim, "Monte Carlo permutations for statistical inference..."))
  I_sim <- numeric(nsim)
  const_multiplier <- (n / S0) / var_z

  for (i in seq_len(nsim)) {
    z_rand <- sample(z)
    W_zz_rand <- sum(W_std * (z_rand %*% t(z_rand)))
    I_sim[i] <- const_multiplier * W_zz_rand
  }

  # Calculate pseudo p-value
  if (I_emp > 0) {
    p_value <- (sum(I_sim >= I_emp) + 1) / (nsim + 1)
  } else {
    p_value <- (sum(I_sim <= I_emp) + 1) / (nsim + 1)
  }

  # 6. Spatial Quadrants Definition (Clinical Risk Grouping)
  quadrants <- rep("Not Significant / Island", n)
  quadrants[z > 0 & lag_z > 0] <- "High-High (Hotspot)"
  quadrants[z < 0 & lag_z < 0] <- "Low-Low (Coldspot)"
  quadrants[z > 0 & lag_z < 0] <- "High-Low (Outlier)"
  quadrants[z < 0 & lag_z > 0] <- "Low-High (Outlier)"
  quadrants[row_sums == 0] <- "Island (No Neighbors)"

  # Append to sf object
  spatial_data[["Z_Score"]] <- z
  spatial_data[["Spatial_Lag_Z"]] <- lag_z
  spatial_data[["Spatial_Lag_Raw"]] <- lag_x
  spatial_data[["Moran_Quadrant"]] <- quadrants

  # 7. Visualization Engine (4K 4-Panel Profiling)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  generate_plot <- function(plot_expr, suffix, title) {
    if (save_plot) {
      filename <- paste0(save_prefix, "_", timestamp, "_", suffix, ".png")
      grDevices::png(filename = filename, width = 3840, height = 2160, res = res)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = c(5, 5, 4, 2))

    eval(plot_expr)
    graphics::title(main = title, outer = TRUE, line = 0, cex.main = 1.5)

    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 4, adj = 0.5, outer = TRUE, col = "darkgray", cex = 0.8)
    }

    if (save_plot) {
      grDevices::dev.off()
      message(paste("AFRIN Note: Saved", filename))
    }
    graphics::par(old_par)
  }

  # PLOT 1: Spatial Connectivity Network Map
  generate_plot(expression({
    graphics::plot(sf::st_geometry(spatial_data), border = "#E0E0E0", col = "#FAFAFA")
    cents <- suppressWarnings(sf::st_centroid(spatial_data))
    coords <- sf::st_coordinates(cents)

    # Draw network edges based on Binary Weights Matrix
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (W_bin[i, j] > 0) {
          graphics::segments(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2],
                             col = grDevices::rgb(0.2, 0.5, 0.8, 0.5), lwd = 0.5)
        }
      }
    }
    graphics::plot(sf::st_geometry(cents), col = "#2C3E50", pch = 16, cex = 0.8, add = TRUE)
    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = c("Centroids", "Contiguity Edges"),
                       col = c("#2C3E50", grDevices::rgb(0.2, 0.5, 0.8, 0.5)),
                       pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n", horiz = horiz)
    }
  }), "1_Connectivity", "Spatial Connectivity Matrix (Queen Contiguity)")

  # PLOT 2: Moran Scatterplot (Z vs Spatial Lag Z)
  generate_plot(expression({
    # Colors for quadrants
    pt_cols <- rep("gray", n)
    pt_cols[quadrants == "High-High (Hotspot)"] <- "red"
    pt_cols[quadrants == "Low-Low (Coldspot)"] <- "blue"
    pt_cols[quadrants == "High-Low (Outlier)"] <- "pink"
    pt_cols[quadrants == "Low-High (Outlier)"] <- "lightblue"

    graphics::plot(z, lag_z, pch = 21, bg = pt_cols, col = "darkgray", cex = 1.5,
                   xlab = "Standardized Observation (Z)", ylab = "Spatial Lag of Z",
                   main = paste("Global Moran's I:", round(I_emp, 4), "| p-value:", round(p_value, 4)))

    # Grid lines & Regression line
    graphics::abline(h = 0, v = 0, lty = 2, col = "darkgray")
    if (stats::var(z) > 0) {
      reg <- stats::lm(lag_z ~ z)
      graphics::abline(reg, col = "black", lwd = 2)
    }

    if (legend) {
      graphics::legend("topleft", legend = c("High-High", "Low-Low", "High-Low", "Low-High"),
                       pt.bg = c("red", "blue", "pink", "lightblue"), pch = 21, bty = "n", cex = 0.9)
    }
  }), "2_Scatterplot", "Moran's I Scatterplot (Spatial Quadrants)")

  # PLOT 3: Monte Carlo Permutation Density
  generate_plot(expression({
    dens <- stats::density(I_sim)
    graphics::plot(dens, main = paste("Monte Carlo Permutations (nsim =", nsim, ")"),
                   xlab = "Simulated Moran's I Values", ylab = "Density", lwd = 2, col = "darkblue")
    graphics::polygon(dens, col = grDevices::rgb(0.2, 0.5, 0.8, 0.3), border = NA)

    # Empirical Line
    graphics::abline(v = I_emp, col = "red", lwd = 3, lty = 1)

    if (legend) {
      graphics::legend("topright", legend = c("Simulated Distribution", "Empirical Observation"),
                       fill = c(grDevices::rgb(0.2, 0.5, 0.8, 0.3), "red"), bty = "n", cex = 0.9)
    }
  }), "3_MonteCarlo", "Monte Carlo Permutation Distribution")

  # PLOT 4: Spatial Lag Reference Map (Choropleth)
  generate_plot(expression({
    # Color Setup
    pal_colors <- NULL
    if (requireNamespace("RColorBrewer", quietly = TRUE) &&
        color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      n_col <- RColorBrewer::brewer.pal.info[color_palette, "maxcolors"]
      pal_colors <- RColorBrewer::brewer.pal(n_col, color_palette)
    } else {
      pal_colors <- grDevices::heat.colors(7)
    }
    if (reverse_palette) pal_colors <- rev(pal_colors)
    pal_generator <- grDevices::colorRampPalette(pal_colors)

    # Mapping logic
    breaks <- unique(stats::quantile(lag_x, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
    if (length(breaks) > 1) {
      cats <- cut(lag_x, breaks = breaks, include.lowest = TRUE)
      colors_mapped <- pal_generator(length(levels(cats)))[as.numeric(cats)]
      lgd_labels <- levels(cats)
      lgd_cols <- pal_generator(length(levels(cats)))
    } else {
      colors_mapped <- rep(pal_colors[length(pal_colors)], n)
      lgd_labels <- "Homogeneous"
      lgd_cols <- pal_colors[length(pal_colors)]
    }

    graphics::plot(sf::st_geometry(spatial_data), col = colors_mapped, border = "#555555", lwd = 0.5)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = lgd_labels, fill = lgd_cols,
                       title = "Avg Neighbor Risk (Lag)", bty = "n", horiz = horiz, cex = 0.9, xpd = NA)
    }
  }), "4_SpatialLagMap", paste("Spatial Lag Reference Map -", var_col))

  # 8. Compile and Return Results
  results <- list(
    Morans_I = I_emp,
    p_value = p_value,
    permutations = nsim,
    spatial_data_appended = spatial_data
  )

  message("AFRIN Note: Global Spatial Autocorrelation completed successfully.")
  return(invisible(results))
}
