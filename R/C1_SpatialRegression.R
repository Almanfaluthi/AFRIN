#' @title Spatial Lag Regression Model (Maximum Likelihood Estimation)
#' @description Fits a Spatial Lag Model (SLM) using pure Base R Maximum Likelihood
#' Estimation. It handles spatial autocorrelation by introducing a spatial lag term
#' (rho * Wy). The function auto-coerces character columns to numeric and generates
#' a 4-panel diagnostic profile including a Strict-Boundary Forest Plot, Predicted
#' Risk Map, Spillover Impact Map, and Spatial Residuals Map.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries.
#' @param formula An object of class 'formula' (e.g., MALARIA ~ DHF + AREA).
#' @param color_palette Character. RColorBrewer palette name (default: "YlOrRd").
#' @param reverse_palette Logical. If TRUE, reverses the color palette vector.
#' @param legend Logical. If TRUE, displays legends on applicable plots.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" for legend layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, saves the 4 diagnostic plots as high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return A list containing the optimal spatial lag parameter (rho), coefficients,
#' standard errors, and the appended 'sf' object invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Modeling Malaria driven by DHF and Spillover Effects
#' slm_results <- C1_SpatialRegression(
#'   spatial_data = df_jateng,
#'   formula = MALARIA ~ DHF,
#'   save_plot = TRUE,
#'   save_prefix = "Malaria_DHF_SpatialReg"
#' )
#' }
C1_SpatialRegression <- function(spatial_data,
                                 formula,
                                 color_palette = "YlOrRd",
                                 reverse_palette = FALSE,
                                 legend = TRUE,
                                 legend_position = "bottom",
                                 legend_orientation = "vertical",
                                 footer = TRUE,
                                 save_plot = FALSE,
                                 save_prefix = "C1_SpatialReg",
                                 res = 300) {

  # 1. Strict Validation & sf Check
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("AFRIN Error: Package 'sf' is required.", call. = FALSE)
  }
  if (!inherits(spatial_data, "sf")) stop("AFRIN Error: 'spatial_data' must be an 'sf' object.", call. = FALSE)
  if (!inherits(formula, "formula")) stop("AFRIN Error: 'formula' must be a valid R formula.", call. = FALSE)

  # 2. Intelligent Auto-Coercion Engine based on Formula
  vars <- all.vars(formula)
  for (v in vars) {
    if (v %in% names(spatial_data)) {
      raw_val <- spatial_data[[v]]
      if (!is.numeric(raw_val)) {
        warning(paste("AFRIN Warning: Variable", v, "is not numeric. Auto-coercing..."), call. = FALSE)
        num_val <- suppressWarnings(as.numeric(as.character(raw_val)))
        num_val[is.na(num_val)] <- 0
        spatial_data[[v]] <- num_val
      }
    } else {
      stop(paste("AFRIN Error: Variable", v, "not found in spatial_data."), call. = FALSE)
    }
  }

  # 3. Data Preparation & Matrix Extraction
  mf <- stats::model.frame(formula, data = spatial_data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  n <- length(y)

  # 4. Pure Base R Spatial Weights Matrix (Row-Standardized Queen)
  message("AFRIN Note: Building Spatial Weights Matrix...")
  touch_list <- sf::st_touches(spatial_data, sparse = FALSE)
  W_bin <- matrix(as.numeric(touch_list), nrow = n, ncol = n)
  row_sums <- rowSums(W_bin)
  row_sums_safe <- ifelse(row_sums == 0, 1, row_sums)
  W_std <- sweep(W_bin, 1, row_sums_safe, FUN = "/")

  Wy <- as.numeric(W_std %*% y)

  # 5. Maximum Likelihood Estimation (MLE) for Spatial Lag Model (rho)
  message("AFRIN Note: Optimizing Spatial Lag Parameter (MLE)...")
  eigW <- Re(eigen(W_std)$values)

  # Concentrated Log-Likelihood Function
  log_lik_rho <- function(rho) {
    A <- diag(n) - rho * W_std
    Ay <- as.numeric(A %*% y)
    beta_tmp <- solve(crossprod(X)) %*% crossprod(X, Ay)
    res_tmp <- Ay - as.numeric(X %*% beta_tmp)
    sig2_tmp <- sum(res_tmp^2) / n
    # MLE equation
    ll <- sum(log(1 - rho * eigW)) - (n / 2) * log(sig2_tmp)
    return(ll)
  }

  # Optimize rho bounded strictly between -0.99 and 0.99 to prevent explosion
  opt <- stats::optimize(log_lik_rho, interval = c(-0.99, 0.99), maximum = TRUE)
  rho_opt <- opt$maximum

  # 6. Final Parameter Estimation
  A_opt <- diag(n) - rho_opt * W_std
  Ay_opt <- as.numeric(A_opt %*% y)
  beta_opt <- solve(crossprod(X)) %*% crossprod(X, Ay_opt)

  fitted_direct <- as.numeric(X %*% beta_opt)
  spillover_effect <- rho_opt * Wy
  fitted_total <- fitted_direct + spillover_effect
  residuals_opt <- y - fitted_total

  # Variance & Standard Errors (Beta)
  sig2_opt <- sum(residuals_opt^2) / n
  var_beta <- sig2_opt * solve(crossprod(X))
  se_beta <- sqrt(diag(var_beta))

  # Calculate 95% Confidence Intervals
  ci_lower <- beta_opt[, 1] - 1.96 * se_beta
  ci_upper <- beta_opt[, 1] + 1.96 * se_beta

  spatial_data[["Predicted_Risk"]] <- fitted_total
  spatial_data[["Spillover_Impact"]] <- spillover_effect
  spatial_data[["Spatial_Residuals"]] <- residuals_opt

  # 7. High-Res Visualization Engine
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

  generate_plot <- function(plot_expr, suffix, title) {
    if (save_plot) {
      filename <- paste0(save_prefix, "_", timestamp, "_", suffix, ".png")
      grDevices::png(filename = filename, width = 3840, height = 2160, res = res)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = c(4, 4, 3, 2))

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

  # PLOT 1: Strict-Boundary Coefficient Forest Plot
  generate_plot(expression({
    n_coef <- length(beta_opt)
    term_names <- rownames(beta_opt)

    # Strict Clinical Logic: If CI crosses 0, it is Not Significant
    is_sig <- (ci_lower * ci_upper) > 0
    plot_cols <- ifelse(is_sig, "#27AE60", "#BDC3C7") # Green if sig, Grey if ns

    graphics::plot(beta_opt[, 1], 1:n_coef, xlim = range(c(ci_lower, ci_upper, 0)),
                   ylim = c(0.5, n_coef + 0.5), yaxt = "n", xlab = "Estimate (95% CI)",
                   ylab = "", pch = 16, col = plot_cols, cex = 1.5)

    graphics::axis(2, at = 1:n_coef, labels = term_names, las = 1, cex.axis = 0.8)
    graphics::abline(v = 0, lty = 2, col = "#E74C3C", lwd = 2) # Null line

    # Draw Confidence Intervals
    graphics::arrows(ci_lower, 1:n_coef, ci_upper, 1:n_coef, angle = 90, code = 3,
                     length = 0.05, col = plot_cols, lwd = 2)

    # Add Spatial Lag (Rho) text context
    graphics::text(x = min(ci_lower), y = n_coef + 0.4,
                   labels = paste("Spatial Lag (rho):", round(rho_opt, 4)),
                   pos = 4, font = 2, col = "#2980B9", cex = 0.9)

    if (legend) {
      graphics::legend("topright", legend = c("Significant (CI excludes 0)", "Not Significant"),
                       col = c("#27AE60", "#BDC3C7"), pch = 16, lty = 1, lwd = 2, bty = "n")
    }
  }), "1_ForestPlot", "Spatial Regression Coefficient Forest Plot")

  # Helper to plot choropleths cleanly
  plot_choro <- function(var_data, legend_title, div_palette = FALSE) {
    brks <- unique(stats::quantile(var_data, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
    if (length(brks) > 1) {
      cats <- cut(var_data, breaks = brks, include.lowest = TRUE)

      if (div_palette) {
        # Diverging palette for residuals (Blue to Red)
        div_cols <- grDevices::colorRampPalette(c("#313695", "#E0F3F8", "#A50026"))(length(levels(cats)))
        cols_mapped <- div_cols[as.numeric(cats)]
        lgd_cols <- div_cols
      } else {
        cols_mapped <- pal_func(length(levels(cats)))[as.numeric(cats)]
        lgd_cols <- pal_func(length(levels(cats)))
      }
      lgd_labels <- levels(cats)
    } else {
      cols_mapped <- rep(pal_colors[3], n)
      lgd_labels <- "Homogeneous"
      lgd_cols <- pal_colors[3]
    }

    graphics::plot(sf::st_geometry(spatial_data), col = cols_mapped, border = "#555555", lwd = 0.5)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = lgd_labels, fill = lgd_cols,
                       title = legend_title, bty = "n", horiz = horiz, cex = 0.8, xpd = NA)
    }
  }

  # PLOT 2: Predicted Risk Map
  generate_plot(expression({
    plot_choro(fitted_total, "Fitted Values")
  }), "2_PredictedRisk", "Model Predicted Risk Distribution")

  # PLOT 3: Spillover Impact Map
  generate_plot(expression({
    plot_choro(spillover_effect, "Indirect / Spillover Effect")
  }), "3_SpilloverImpact", "Spatial Spillover Impact (Rho * Wy)")

  # PLOT 4: Spatial Residuals Map
  generate_plot(expression({
    plot_choro(residuals_opt, "Residual Error", div_palette = TRUE)
  }), "4_Residuals", "Spatial Residuals Distribution (Randomness Check)")

  # 8. Compile Results
  res_list <- list(
    rho = rho_opt,
    coefficients = beta_opt,
    standard_errors = se_beta,
    spatial_data_appended = spatial_data
  )

  message("AFRIN Note: Spatial Regression (SLM) via MLE completed successfully.")
  return(invisible(res_list))
}
