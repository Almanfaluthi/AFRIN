#' @title Spatially Constrained Regionalization and Syndemic Zoning
#' @description Performs spatially constrained multivariate clustering to group contiguous
#' regions based on similar disease burdens. Incorporates a Zone Core vs. Borderline
#' Fit Map to evaluate zone stability and syndemic homogeneity.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries.
#' @param var_cols Character vector. Columns representing the variables for clustering.
#' @param k_zones Numeric. The number of regions/zones to create.
#' @param color_palette Character. RColorBrewer palette name (default: "Set3").
#' @param reverse_palette Logical. If TRUE, reverses the color palette vector.
#' @param legend Logical. If TRUE, displays legends on applicable plots.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" for legend layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, saves the 4 diagnostic plots as high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return A list containing the zone assignments, fit scores, and the appended 'sf' object invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Grouping districts into 3 Syndemic Quarantine Zones
#' zoning_results <- C3_Regionalization(
#'   spatial_data = df_jateng,
#'   var_cols = c("MALARIA", "DHF"),
#'   k_zones = 3,
#'   save_plot = TRUE,
#'   save_prefix = "Syndemic_Zone"
#' )
#' }
C3_Regionalization <- function(spatial_data,
                               var_cols,
                               k_zones = 3,
                               color_palette = "Set3",
                               reverse_palette = FALSE,
                               legend = TRUE,
                               legend_position = "bottom",
                               legend_orientation = "horizontal",
                               footer = TRUE,
                               save_plot = FALSE,
                               save_prefix = "C3_Regionalization",
                               res = 300) {

  # 1. Strict Validation
  if (!requireNamespace("sf", quietly = TRUE)) stop("AFRIN Error: Package 'sf' is required.", call. = FALSE)
  if (!inherits(spatial_data, "sf")) stop("AFRIN Error: 'spatial_data' must be an 'sf' object.", call. = FALSE)
  if (k_zones < 2) stop("AFRIN Error: 'k_zones' must be at least 2.", call. = FALSE)

  n <- nrow(spatial_data)
  if (k_zones >= n) stop("AFRIN Error: 'k_zones' cannot exceed the number of polygons.", call. = FALSE)

  # 2. Intelligent Auto-Coercion Engine
  data_mat <- matrix(NA, nrow = n, ncol = length(var_cols))
  colnames(data_mat) <- var_cols

  for (i in seq_along(var_cols)) {
    vc <- var_cols[i]
    if (!(vc %in% names(spatial_data))) stop(paste("AFRIN Error: Column", vc, "not found."), call. = FALSE)

    raw_val <- spatial_data[[vc]]
    if (!is.numeric(raw_val)) {
      warning(paste("AFRIN Warning:", vc, "is not numeric. Auto-coercing..."), call. = FALSE)
      num_val <- suppressWarnings(as.numeric(as.character(raw_val)))
      num_val[is.na(num_val)] <- 0
      spatial_data[[vc]] <- num_val
      data_mat[, i] <- num_val
    } else {
      raw_val[is.na(raw_val)] <- 0
      data_mat[, i] <- raw_val
    }
  }

  # 3. Multivariate Scaling & Modified Spatial Penalty Engine
  message("AFRIN Note: Computing Spatially Constrained Distance Matrix...")
  scaled_data <- scale(data_mat)
  dist_euclid <- as.matrix(stats::dist(scaled_data))

  # Queen Contiguity Matrix
  touch_list <- sf::st_touches(spatial_data, sparse = FALSE)
  W_bin <- matrix(as.numeric(touch_list), nrow = n, ncol = n)

  # Apply massive penalty to non-adjacent polygons to force contiguous clustering
  max_dist <- max(dist_euclid, na.rm = TRUE)
  penalty <- max_dist * 1000

  dist_constrained <- dist_euclid
  dist_constrained[W_bin == 0] <- dist_constrained[W_bin == 0] + penalty
  diag(dist_constrained) <- 0 # Distance to self is 0

  dist_obj <- stats::as.dist(dist_constrained)

  # 4. Hierarchical Clustering (Ward's Method)
  hc <- stats::hclust(dist_obj, method = "ward.D2")
  zones <- stats::cutree(hc, k = k_zones)

  zone_labels <- paste("Zone", zones)
  spatial_data[["Assigned_Zone"]] <- factor(zone_labels, levels = paste("Zone", 1:k_zones))

  # 5. Zone Profiles & Fit Score Calculation
  message("AFRIN Note: Calculating Multivariate Fit Scores for Core vs Borderline Map...")
  zone_means <- stats::aggregate(scaled_data, by = list(Zone = zones), FUN = mean)

  dist_to_center <- numeric(n)
  for(i in 1:n) {
    z_id <- zones[i]
    # Extract the centroid vector for this zone
    center_vec <- as.numeric(zone_means[zone_means$Zone == z_id, -1])
    # Euclidean distance in multivariate space
    dist_to_center[i] <- sqrt(sum((scaled_data[i, ] - center_vec)^2))
  }

  # Normalize Fit Score (0 to 1) -> 1 is best fit (Core), 0 is worst fit (Borderline)
  max_dev <- max(dist_to_center, na.rm = TRUE)
  fit_score <- 1 - (dist_to_center / (max_dev + 1e-9))
  spatial_data[["Zone_Fit_Score"]] <- fit_score

  # 6. High-Res Visualization Engine
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  # Dynamic Categorical Palette
  pal_colors <- NULL
  if (requireNamespace("RColorBrewer", quietly = TRUE) && color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal_colors <- RColorBrewer::brewer.pal(max(3, k_zones), color_palette)
    if (k_zones > length(pal_colors)) {
      pal_colors <- grDevices::colorRampPalette(pal_colors)(k_zones)
    } else {
      pal_colors <- pal_colors[1:k_zones]
    }
  } else {
    pal_colors <- grDevices::rainbow(k_zones)
  }
  if (reverse_palette) pal_colors <- rev(pal_colors)

  generate_plot <- function(plot_expr, suffix, title, custom_mar = c(4, 4, 3, 2)) {
    if (save_plot) {
      filename <- paste0(save_prefix, "_", timestamp, "_", suffix, ".png")
      grDevices::png(filename = filename, width = 3840, height = 2160, res = res)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = custom_mar)

    eval(plot_expr)
    graphics::title(main = paste("AFRIN:", title), outer = TRUE, line = 0, cex.main = 1.5)

    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 4, adj = 0.5, outer = TRUE, col = "darkgray", cex = 0.8)
    }

    if (save_plot) {
      grDevices::dev.off()
      message(paste("AFRIN Note: Saved", filename))
    }
    graphics::par(old_par)
  }

  lgd_horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)

  # PLOT 1: Spatially Constrained Zoning Map
  generate_plot(expression({
    geom_cols <- pal_colors[zones]
    graphics::plot(sf::st_geometry(spatial_data), col = geom_cols, border = "#555555", lwd = 0.5)

    if (legend) {
      graphics::legend(legend_position, legend = levels(spatial_data[["Assigned_Zone"]]),
                       fill = pal_colors, bty = "n", horiz = lgd_horiz, cex = 0.9, xpd = NA)
    }
  }), "1_ZoningMap", "Spatially Constrained Syndemic Zones")

  # PLOT 2: Zone Core vs. Borderline Fit Map (NEW!)
  generate_plot(expression({
    # Convert base colors to RGB and apply Fit Score as Alpha (Transparency)
    # Range alpha from 0.15 (very transparent) to 0.95 (solid core)
    alpha_scaled <- 0.15 + (fit_score * 0.8)
    rgb_mat <- grDevices::col2rgb(pal_colors[zones]) / 255
    fit_cols <- grDevices::rgb(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ], alpha = alpha_scaled)

    graphics::plot(sf::st_geometry(spatial_data), col = fit_cols, border = "#555555", lwd = 0.5)

    if (legend) {
      # Legend demonstrating opacity gradient
      graphics::legend(legend_position, legend = c("Zone Core (High Fit Score)", "Zone Borderline (Low Fit Score)"),
                       fill = c(grDevices::rgb(0.5, 0.5, 0.5, 0.95), grDevices::rgb(0.5, 0.5, 0.5, 0.3)),
                       border = c("black", "black"),
                       bty = "n", horiz = lgd_horiz, cex = 0.9, xpd = NA)
    }
  }), "2_FitScoreMap", "Zone Core vs. Borderline (Multivariate Fit Score)")

  # PLOT 3: Connectivity & Edge Cut Network
  generate_plot(expression({
    graphics::plot(sf::st_geometry(spatial_data), col = "#FAFAFA", border = "#D0D0D0")
    cents <- suppressWarnings(sf::st_centroid(spatial_data))
    coords <- sf::st_coordinates(cents)

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (W_bin[i, j] > 0) {
          if (zones[i] == zones[j]) {
            graphics::segments(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2],
                               col = pal_colors[zones[i]], lwd = 2)
          } else {
            graphics::segments(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2],
                               col = "#E74C3C", lwd = 1, lty = 3)
          }
        }
      }
    }

    graphics::plot(sf::st_geometry(cents), pch = 21, bg = pal_colors[zones], col = "#333333", cex = 1.2, add = TRUE)

    if (legend) {
      graphics::legend(legend_position, legend = c("Within-Zone Link", "Zone Border Cut (Edge Cut)"),
                       col = c("darkgray", "#E74C3C"), lty = c(1, 3), lwd = c(2, 1),
                       bty = "n", horiz = lgd_horiz, xpd = NA, cex = 0.9)
    }
  }), "3_ConnectivityNetwork", "Spatial Spanning Tree & Boundary Edge Cuts")

  # PLOT 4: Multivariate Zone Profile (Z-Scores)
  generate_plot(expression({
    plot_data <- t(zone_means[, -1, drop = FALSE])

    graphics::matplot(plot_data, type = "b", pch = 15:(14+k_zones), lty = 1, lwd = 2,
                      col = pal_colors, xaxt = "n",
                      ylab = "Standardized Disease Burden (Z-Score Mean)",
                      xlab = "Syndemic Variables",
                      main = "")

    graphics::axis(1, at = 1:nrow(plot_data), labels = rownames(plot_data), cex.axis = 0.8)
    graphics::abline(h = 0, lty = 2, col = "darkgray")

    if (legend) {
      graphics::legend(legend_position, legend = levels(spatial_data[["Assigned_Zone"]]),
                       col = pal_colors, pch = 15:(14+k_zones), lty = 1, lwd = 2,
                       bty = "n", horiz = lgd_horiz, xpd = NA, cex = 0.9)
    }
  }), "4_MultivariateProfile", "Multivariate Disease Profile per Zone", custom_mar = c(5, 5, 3, 2))

  # 7. Compile Results
  res_list <- list(
    Zone_Assignments = sf::st_drop_geometry(spatial_data)[, c(var_cols, "Assigned_Zone", "Zone_Fit_Score")],
    Zone_Profiles = zone_means,
    spatial_data_appended = spatial_data
  )

  message("AFRIN Note: Spatially Constrained Regionalization with Fit Score completed successfully.")
  return(invisible(res_list))
}
