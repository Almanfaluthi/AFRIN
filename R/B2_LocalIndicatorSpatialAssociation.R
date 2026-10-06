#' @title Local Indicator of Spatial Association (LISA) / Local Moran's I
#' @description Computes Local Moran's I to detect precise spatial epicenters
#' (Hotspots/Coldspots) and Spatial Outliers. Includes a strict clinical boundary
#' check that masks non-significant clusters (p > 0.05) to prevent false-alarm
#' interventions. Automatically generates a 4-panel high-resolution diagnostic profile.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries.
#' @param var_col Character. The column name of the variable to analyze.
#' @param nsim Numeric. Number of Monte Carlo permutations for pseudo p-value (default: 999).
#' @param p_threshold Numeric. The significance threshold for clinical masking (default: 0.05).
#' @param color_palette Character. Not heavily used here as LISA uses fixed clinical colors, but kept for compatibility.
#' @param reverse_palette Logical. If TRUE, reverses the generic palette.
#' @param legend Logical. If TRUE, displays legends on applicable plots.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" for legend layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, copies the viewer plots to high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return A list containing the Local Moran statistics and the original 'sf'
#' object appended with Local_I, p_value, and LISA_Cluster invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Detecting Malaria Epicenters
#' lisa_results <- B2_LocalIndicatorSpatialAssociation(
#'   spatial_data = df_jateng,
#'   var_col = "MALARIA",
#'   nsim = 999,
#'   p_threshold = 0.05,
#'   save_plot = TRUE,
#'   save_prefix = "Malaria_LISA"
#' )
#' }
B2_LocalIndicatorSpatialAssociation <- function(spatial_data,
                                                var_col,
                                                nsim = 999,
                                                p_threshold = 0.05,
                                                color_palette = "YlOrRd",
                                                reverse_palette = FALSE,
                                                legend = TRUE,
                                                legend_position = "bottom",
                                                legend_orientation = "vertical",
                                                footer = TRUE,
                                                save_plot = FALSE,
                                                save_prefix = "B2_LISA",
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
    warning(paste("AFRIN Warning: 'var_col' (", var_col, ") is not numeric. Auto-coercing..."), call. = FALSE)
    num_var <- suppressWarnings(as.numeric(as.character(raw_var)))

    na_diff <- sum(is.na(num_var) & !is.na(raw_var))
    if (na_diff > 0) {
      warning(paste("AFRIN Warning: Coercion introduced NAs in", na_diff, "rows. Treated as 0."), call. = FALSE)
    }
    num_var[is.na(num_var)] <- 0
    spatial_data[[var_col]] <- num_var
  } else {
    raw_var[is.na(raw_var)] <- 0
    spatial_data[[var_col]] <- raw_var
  }

  x <- spatial_data[[var_col]]
  n <- length(x)

  # 3. Base R Spatial Weights Matrix (Queen Contiguity)
  message("AFRIN Note: Building Spatial Weights Matrix & Computing LISA...")
  touch_list <- sf::st_touches(spatial_data, sparse = FALSE)
  W_bin <- matrix(as.numeric(touch_list), nrow = n, ncol = n)

  row_sums <- rowSums(W_bin)
  row_sums_safe <- ifelse(row_sums == 0, 1, row_sums)
  W_std <- sweep(W_bin, 1, row_sums_safe, FUN = "/")

  # 4. Local Moran's I Calculation
  x_bar <- mean(x, na.rm = TRUE)
  z <- x - x_bar
  lag_z <- as.numeric(W_std %*% z)
  m2 <- sum(z^2) / n

  # Empirical Local Moran's I formula
  I_local <- (z / m2) * lag_z

  # 5. Local Monte Carlo Permutations (Simulated p-value matrix)
  message(paste("AFRIN Note: Running", nsim, "Local Monte Carlo permutations..."))
  I_sim_mat <- matrix(0, nrow = n, ncol = nsim)

  for (s in seq_len(nsim)) {
    z_rand <- sample(z)
    lag_z_sim <- as.numeric(W_std %*% z_rand)
    I_sim_mat[, s] <- (z / m2) * lag_z_sim
  }

  # Pseudo p-value calculation (two-tailed logic approximation for local I)
  p_val <- numeric(n)
  for (i in seq_len(n)) {
    if (row_sums[i] == 0) {
      p_val[i] <- NA # Island geometries have no neighbors
    } else {
      if (I_local[i] > 0) {
        p_val[i] <- (sum(I_sim_mat[i, ] >= I_local[i]) + 1) / (nsim + 1)
      } else {
        p_val[i] <- (sum(I_sim_mat[i, ] <= I_local[i]) + 1) / (nsim + 1)
      }
    }
  }

  # 6. Clinical Boundary Filter & Quadrant Classification
  quadrants <- rep("Not Significant", n)
  quadrants[row_sums == 0] <- "Island"

  # STRICT MASK: Only classify if p-value <= threshold
  sig_mask <- !is.na(p_val) & (p_val <= p_threshold)

  quadrants[z > 0 & lag_z > 0 & sig_mask] <- "High-High (Hotspot)"
  quadrants[z < 0 & lag_z < 0 & sig_mask] <- "Low-Low (Coldspot)"
  quadrants[z > 0 & lag_z < 0 & sig_mask] <- "High-Low (Outlier)"
  quadrants[z < 0 & lag_z > 0 & sig_mask] <- "Low-High (Outlier)"

  # Append to sf object
  spatial_data[["Z_Score"]] <- z
  spatial_data[["Spatial_Lag_Z"]] <- lag_z
  spatial_data[["Local_Morans_I"]] <- I_local
  spatial_data[["LISA_P_Value"]] <- p_val
  spatial_data[["LISA_Cluster"]] <- quadrants

  # 7. High-Res Visualization Engine
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  # LISA Colors Strategy
  lisa_cols <- c("High-High (Hotspot)" = "#E31A1C",  # Strong Red
                 "Low-Low (Coldspot)"  = "#1F78B4",  # Strong Blue
                 "High-Low (Outlier)"  = "#FB9A99",  # Pink
                 "Low-High (Outlier)"  = "#A6CEE3",  # Light Blue
                 "Not Significant"     = "#E0E0E0",  # Grey
                 "Island"              = "#FFFFFF")  # White

  # Internal plotting wrapper
  generate_lisa_plot <- function(plot_expr, suffix, title) {
    old_par <- graphics::par(no.readonly = TRUE)
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = c(4, 4, 3, 2))

    eval(plot_expr)
    graphics::title(main = paste("AFRIN LISA:", title), outer = TRUE, line = 0, cex.main = 1.5)

    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 4, adj = 0.5, outer = TRUE, col = "darkgray", cex = 0.8)
    }

    if (save_plot) {
      filename <- paste0(save_prefix, "_", timestamp, "_", suffix, ".png")
      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = res)
      grDevices::dev.off()
      message(paste("AFRIN Note: Saved", filename))
    }
    graphics::par(old_par)
  }

  # PLOT 1: LISA Cluster Map (The Epicenter Map)
  generate_lisa_plot(expression({
    geom_cols <- lisa_cols[quadrants]
    graphics::plot(sf::st_geometry(spatial_data), col = geom_cols, border = "#555555", lwd = 0.5)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      active_cats <- intersect(names(lisa_cols), unique(quadrants))
      graphics::legend(legend_position, legend = active_cats, fill = lisa_cols[active_cats],
                       title = paste("Cluster (p <", p_threshold, ")"),
                       bty = "n", horiz = horiz, cex = 0.9, xpd = NA)
    }
  }), "1_LISA_Epicenters", "Epicenter & Cluster Distribution")

  # PLOT 2: LISA Significance Map
  generate_lisa_plot(expression({
    sig_cats <- rep("Not Significant", n)
    sig_cats[p_val <= 0.05] <- "p <= 0.05"
    sig_cats[p_val <= 0.01] <- "p <= 0.01"
    sig_cats[p_val <= 0.001] <- "p <= 0.001"
    sig_cats[row_sums == 0] <- "Island"

    sig_cols <- c("p <= 0.001" = "#005a32", "p <= 0.01" = "#238b45",
                  "p <= 0.05" = "#74c476", "Not Significant" = "#E0E0E0", "Island" = "#FFFFFF")

    graphics::plot(sf::st_geometry(spatial_data), col = sig_cols[sig_cats], border = "#555555", lwd = 0.5)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      active_sig <- intersect(names(sig_cols), unique(sig_cats))
      graphics::legend(legend_position, legend = active_sig, fill = sig_cols[active_sig],
                       title = "Monte Carlo Significance", bty = "n", horiz = horiz, cex = 0.9, xpd = NA)
    }
  }), "2_LISA_Significance", "Statistical Significance Map")

  # PLOT 3: Spatial Outlier Isolation Map (Tactical Map)
  generate_lisa_plot(expression({
    # Darken background
    graphics::par(bg = "#222222", fg = "white", col.axis = "white", col.main = "white")

    outlier_cols <- rep("#333333", n) # Dark grey for non-outliers
    outlier_cols[quadrants == "High-Low (Outlier)"] <- "#FB9A99"
    outlier_cols[quadrants == "Low-High (Outlier)"] <- "#A6CEE3"

    graphics::plot(sf::st_geometry(spatial_data), col = outlier_cols, border = "#555555", lwd = 0.5)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = c("High-Low (Isolated Risk)", "Low-High (Surrounded Risk)"),
                       fill = c("#FB9A99", "#A6CEE3"), title = "Spatial Outliers Only",
                       bty = "n", horiz = horiz, text.col = "white", cex = 0.9, xpd = NA)
    }
  }), "3_LISA_Outliers", "Outlier Isolation & Containment Strategy")

  # PLOT 4: Highlighted Moran Scatterplot
  generate_lisa_plot(expression({
    # Reset bg to white for scatterplot
    graphics::par(bg = "white", fg = "black", col.axis = "black", col.main = "black")

    pt_cols <- lisa_cols[quadrants]
    graphics::plot(z, lag_z, pch = 21, bg = pt_cols, col = "#555555", cex = 1.5,
                   xlab = "Standardized Observation (Z)", ylab = "Spatial Lag of Z",
                   main = paste("Local Moran Scatterplot (Masked at p <", p_threshold, ")"))

    graphics::abline(h = 0, v = 0, lty = 2, col = "darkgray")
    if (stats::var(z) > 0) {
      reg <- stats::lm(lag_z ~ z)
      graphics::abline(reg, col = "black", lwd = 2)
    }

    if (legend) {
      active_cats <- intersect(names(lisa_cols), unique(quadrants))
      graphics::legend("topleft", legend = active_cats, pt.bg = lisa_cols[active_cats],
                       pch = 21, bty = "n", cex = 0.8)
    }
  }), "4_LISA_Scatterplot", "Clinical Risk Quadrants Projection")

  # 8. Compile Results
  results <- list(
    LISA_Stats = data.frame(Local_I = I_local, p_value = p_val, Cluster = quadrants),
    spatial_data_appended = spatial_data
  )

  message("AFRIN Note: Local Indicator of Spatial Association completed.")
  return(invisible(results))
}
