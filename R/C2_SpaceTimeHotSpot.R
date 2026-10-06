#' @title Space-Time Hotspot and Epidemic Trajectory Analysis
#' @description Analyzes the spatiotemporal movement of disease epicenters across multiple
#' time periods. Calculates the Center of Gravity (CoG) trajectory, Emerging Hotspot trends,
#' and generates a Hovmoller Space-Time diagram to track disease propagation waves.
#' Includes an integrated data simulator for longitudinal testing.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries.
#' @param time_cols Character vector. Column names representing sequential time periods.
#' @param simulate_data Logical. If TRUE, injects simulated longitudinal data (Month_1 to Month_4) mimicking a West-to-East epidemic spread.
#' @param color_palette Character. RColorBrewer palette name (default: "YlOrRd").
#' @param reverse_palette Logical. If TRUE, reverses the color palette vector.
#' @param legend Logical. If TRUE, displays legends on applicable plots.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, saves the 4 diagnostic plots as high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return A list containing CoG coordinates, classification data, and the appended 'sf' object invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Simulating and tracking a West-to-East outbreak across 4 months
#' st_results <- C2_SpaceTimeHotSpot(
#'   spatial_data = df_jateng,
#'   time_cols = c("Month_1", "Month_2", "Month_3", "Month_4"),
#'   simulate_data = TRUE,
#'   save_plot = TRUE,
#'   save_prefix = "Trajectory"
#' )
#' }
C2_SpaceTimeHotSpot <- function(spatial_data,
                                time_cols = c("Month_1", "Month_2", "Month_3", "Month_4"),
                                simulate_data = FALSE,
                                color_palette = "YlOrRd",
                                reverse_palette = FALSE,
                                legend = TRUE,
                                legend_position = "bottom",
                                legend_orientation = "vertical",
                                footer = TRUE,
                                save_plot = FALSE,
                                save_prefix = "C2_SpaceTime",
                                res = 300) {

  # 1. Strict Validation
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("AFRIN Error: Package 'sf' is required.", call. = FALSE)
  }
  if (!inherits(spatial_data, "sf")) stop("AFRIN Error: 'spatial_data' must be an 'sf' object.", call. = FALSE)

  cents <- suppressWarnings(sf::st_centroid(spatial_data))
  coords <- sf::st_coordinates(cents)
  n <- nrow(spatial_data)

  # 2. Epidemic Simulator Engine (West-to-East Propagation)
  if (simulate_data) {
    message("AFRIN Note: Simulating West-to-East longitudinal disease propagation...")
    min_x <- min(coords[, 1])
    max_x <- max(coords[, 1])
    x_range <- max_x - min_x

    # Generate 4 time points based on longitude proximity
    for (t in 1:4) {
      target_x <- min_x + (x_range * (t / 5)) # Shifting peak to the east
      dist_to_peak <- abs(coords[, 1] - target_x)

      # Inverse distance with some noise
      cases_sim <- 1000 / (dist_to_peak + 0.1) + stats::runif(n, 0, 50)
      spatial_data[[paste0("Month_", t)]] <- round(cases_sim)
    }
    time_cols <- c("Month_1", "Month_2", "Month_3", "Month_4")
  }

  # Ensure time columns exist and coerce to numeric
  for (tc in time_cols) {
    if (!(tc %in% names(spatial_data))) {
      stop(paste("AFRIN Error: Time column", tc, "not found."), call. = FALSE)
    }
    if (!is.numeric(spatial_data[[tc]])) {
      spatial_data[[tc]] <- suppressWarnings(as.numeric(as.character(spatial_data[[tc]])))
      spatial_data[[tc]][is.na(spatial_data[[tc]])] <- 0
    }
  }

  # 3. Center of Gravity (CoG) Trajectory Calculation
  message("AFRIN Note: Computing Center of Gravity Trajectory...")
  num_times <- length(time_cols)
  cog_coords <- matrix(NA, nrow = num_times, ncol = 2)

  for (i in seq_len(num_times)) {
    cases <- spatial_data[[time_cols[i]]]
    total_cases <- sum(cases, na.rm = TRUE)
    if (total_cases > 0) {
      cog_coords[i, 1] <- sum(coords[, 1] * cases, na.rm = TRUE) / total_cases
      cog_coords[i, 2] <- sum(coords[, 2] * cases, na.rm = TRUE) / total_cases
    } else {
      cog_coords[i, ] <- c(NA, NA) # Fallback if 0 cases
    }
  }

  # 4. Emerging Hotspot Trend Classification (Early vs Late phase)
  message("AFRIN Note: Classifying Spatiotemporal Trends...")
  mid_point <- floor(num_times / 2)
  early_avg <- rowMeans(sf::st_drop_geometry(spatial_data)[, time_cols[1:mid_point], drop = FALSE], na.rm = TRUE)
  late_avg <- rowMeans(sf::st_drop_geometry(spatial_data)[, time_cols[(mid_point+1):num_times], drop = FALSE], na.rm = TRUE)

  global_avg <- mean(c(early_avg, late_avg), na.rm = TRUE)

  trend_class <- rep("Coldspot / Quiet", n)
  trend_class[early_avg < global_avg & late_avg >= global_avg] <- "Emerging Hotspot"
  trend_class[early_avg >= global_avg & late_avg >= global_avg] <- "Persistent Hotspot"
  trend_class[early_avg >= global_avg & late_avg < global_avg] <- "Diminishing (Cooling)"

  spatial_data[["ST_Trend_Class"]] <- trend_class

  # 5. Hovmoller Data Prep (Distance from Initial Epicenter)
  init_cog <- cog_coords[1, ]
  dists_to_init <- sqrt((coords[, 1] - init_cog[1])^2 + (coords[, 2] - init_cog[2])^2)

  dist_breaks <- seq(0, max(dists_to_init, na.rm = TRUE), length.out = 10)
  dist_bins <- cut(dists_to_init, breaks = dist_breaks, include.lowest = TRUE)
  hov_matrix <- matrix(NA, nrow = num_times, ncol = length(levels(dist_bins)))

  for (i in seq_len(num_times)) {
    cases <- spatial_data[[time_cols[i]]]
    hov_matrix[i, ] <- tapply(cases, dist_bins, mean, na.rm = TRUE)
  }
  hov_matrix[is.na(hov_matrix)] <- 0

  # 6. High-Res Visualization Engine
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  pal_colors <- NULL
  if (requireNamespace("RColorBrewer", quietly = TRUE) && color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal_colors <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[color_palette, "maxcolors"], color_palette)
  } else {
    pal_colors <- grDevices::heat.colors(7)
  }
  if (reverse_palette) pal_colors <- rev(pal_colors)
  pal_func <- grDevices::colorRampPalette(pal_colors)

  generate_plot <- function(plot_expr, suffix, title, custom_mar = c(4, 4, 3, 2)) {
    if (save_plot) {
      filename <- paste0(save_prefix, "_", timestamp, "_", suffix, ".png")
      grDevices::png(filename = filename, width = 3840, height = 2160, res = res)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = custom_mar)

    eval(plot_expr)

    if (suffix != "1_TimeLapse") { # Prevent double title in mfrow
      graphics::title(main = paste("AFRIN:", title), outer = TRUE, line = 0, cex.main = 1.5)
    }

    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 4, adj = 0.5, outer = TRUE, col = "darkgray", cex = 0.8)
    }

    if (save_plot) {
      grDevices::dev.off()
      message(paste("AFRIN Note: Saved", filename))
    }
    graphics::par(old_par)
  }

  # PLOT 1: Time-Lapse Facet Map (Small Multiples)
  generate_plot(expression({
    plot_rows <- ceiling(sqrt(num_times))
    plot_cols <- ceiling(num_times / plot_rows)
    graphics::par(mfrow = c(plot_rows, plot_cols), mar = c(1, 1, 2, 1))

    all_vals <- unlist(sf::st_drop_geometry(spatial_data)[, time_cols])
    brks <- unique(stats::quantile(all_vals, probs = seq(0, 1, length.out = 6), na.rm = TRUE))

    for (i in seq_len(num_times)) {
      t_vals <- spatial_data[[time_cols[i]]]

      if (length(brks) > 1) {
        cats <- cut(t_vals, breaks = brks, include.lowest = TRUE)
        cols_mapped <- pal_func(length(levels(cats)))[as.numeric(cats)]
      } else {
        cols_mapped <- rep(pal_colors[3], n)
      }

      graphics::plot(sf::st_geometry(spatial_data), col = cols_mapped, border = "#7F8C8D", lwd = 0.3)
      graphics::title(main = time_cols[i], line = -1)
    }
    graphics::title(main = "AFRIN: Time-Lapse Facet Map (Disease Propagation)", outer = TRUE, line = 0, cex.main = 1.5)
  }), "1_TimeLapse", "Time-Lapse Facet Map", custom_mar = c(2, 2, 2, 2))

  # PLOT 2: Center of Gravity (CoG) Trajectory Map
  generate_plot(expression({
    graphics::plot(sf::st_geometry(spatial_data), col = "#F5F5F5", border = "#BDC3C7")

    # Plot complete trajectory arrows
    if (num_times > 1) {
      for (i in 1:(num_times - 1)) {
        if (!is.na(cog_coords[i, 1]) && !is.na(cog_coords[i+1, 1])) {
          graphics::arrows(cog_coords[i, 1], cog_coords[i, 2],
                           cog_coords[i+1, 1], cog_coords[i+1, 2],
                           col = "#C0392B", lwd = 3, length = 0.15, angle = 20)
        }
      }
    }

    # Plot CoG Nodes
    graphics::points(cog_coords[, 1], cog_coords[, 2], pch = 21, bg = "#F1C40F", col = "#C0392B", cex = 2)
    graphics::text(cog_coords[, 1], cog_coords[, 2], labels = 1:num_times, pos = 3, col = "black", font = 2)

    if (legend) {
      graphics::legend(legend_position, legend = c("CoG Node (Time Order)", "Trajectory Vector"),
                       pt.bg = c("#F1C40F", NA), col = c("#C0392B", "#C0392B"),
                       pch = c(21, NA), lty = c(NA, 1), lwd = c(NA, 3), bty = "n", xpd = NA)
    }
  }), "2_CoG_Trajectory", "Center of Gravity (CoG) Epidemic Trajectory")

  # PLOT 3: Emerging Hotspot Trend Classification Map
  generate_plot(expression({
    class_cols <- c("Emerging Hotspot" = "#E74C3C",
                    "Persistent Hotspot" = "#8E44AD",
                    "Diminishing (Cooling)" = "#3498DB",
                    "Coldspot / Quiet" = "#ECF0F1")

    geom_cols <- class_cols[trend_class]
    graphics::plot(sf::st_geometry(spatial_data), col = geom_cols, border = "#555555", lwd = 0.5)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      active_cats <- intersect(names(class_cols), unique(trend_class))
      graphics::legend(legend_position, legend = active_cats, fill = class_cols[active_cats],
                       title = "Spatiotemporal Phase", bty = "n", horiz = horiz, cex = 0.9, xpd = NA)
    }
  }), "3_EmergingHotspots", "Spatiotemporal Trend & Phase Classification")

  # PLOT 4: Hovmoller Space-Time Heatmap
  generate_plot(expression({
    # Prepare matrix for image()
    z_mat <- hov_matrix
    x_val <- 1:num_times
    y_val <- 1:ncol(hov_matrix)

    # Draw heatmap
    graphics::image(x = x_val, y = y_val, z = z_mat, col = pal_func(50),
                    xlab = "Time Period", ylab = "Distance Zone from Initial Epicenter",
                    main = "", xaxt = "n", yaxt = "n")

    graphics::axis(1, at = x_val, labels = time_cols)
    graphics::axis(2, at = y_val, labels = paste("Zone", y_val), las = 1)

    # Contour overlay for peak intensity
    graphics::contour(x = x_val, y = y_val, z = z_mat, add = TRUE, col = "#2C3E50", lwd = 1.5)

  }), "4_HovmollerDiagram", "Hovmoller Space-Time Diagram (Propagation Waves)", custom_mar = c(5, 6, 3, 2))

  # 7. Compile Results
  res_list <- list(
    CoG_Coordinates = data.frame(Time = time_cols, X = cog_coords[, 1], Y = cog_coords[, 2]),
    Classification = table(trend_class),
    spatial_data_appended = spatial_data
  )

  message("AFRIN Note: Space-Time Hotspot Analysis completed successfully.")
  return(invisible(res_list))
}
