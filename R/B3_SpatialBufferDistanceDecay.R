#' @title Spatial Buffer, Blank Spot Analysis, and Distance Decay
#' @description Computes spatial buffers around healthcare nodes (e.g., hospitals)
#' to evaluate accessibility, identify population blank spots, and measure distance
#' decay. Generates a 4-panel diagnostic profile to support clinical planning
#' (Golden Hour safety) and strategic hospital expansion.
#'
#' @param spatial_data An object of class 'sf' containing polygon geometries (e.g., districts).
#' @param nodes An object of class 'sf' containing point geometries (e.g., hospitals).
#' @param events Optional. An object of class 'sf' (points) representing patient cases.
#'        If NULL, the function uses spatial_data centroids as a proxy.
#' @param buffer_dist Numeric. The radius for the primary buffer (in the units of the CRS, usually meters).
#' @param color_palette Character. RColorBrewer palette name for thematic mapping (default: "YlOrRd").
#' @param reverse_palette Logical. If TRUE, reverses the color palette vector.
#' @param legend Logical. If TRUE, displays legends on applicable plots.
#' @param legend_position Character. Legend position (default: "bottom").
#' @param legend_orientation Character. "vertical" or "horizontal" for legend layout.
#' @param footer Logical. If TRUE, adds a reproducible footer at the bottom center.
#' @param save_plot Logical. If TRUE, saves the 4 diagnostic plots as high-resolution 4K PNGs.
#' @param save_prefix Character. Prefix for output files.
#' @param res Numeric. Resolution in DPI for the saved plot (default is 300).
#'
#' @return A list containing distance metrics and appended 'sf' objects indicating
#' coverage status for polygons and events invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Evaluating Malaria patient access to PHC (5km radius)
#' access_results <- B3_SpatialBufferDistanceDecay(
#'   spatial_data = df_jateng,
#'   nodes = WBC,
#'   events = MALARIA,
#'   buffer_dist = 5000,
#'   save_plot = TRUE,
#'   save_prefix = "Malaria_GoldenHour"
#' )
#' }
B3_SpatialBufferDistanceDecay <- function(spatial_data,
                                          nodes,
                                          events = NULL,
                                          buffer_dist,
                                          color_palette = "YlOrRd",
                                          reverse_palette = FALSE,
                                          legend = TRUE,
                                          legend_position = "bottom",
                                          legend_orientation = "vertical",
                                          footer = TRUE,
                                          save_plot = FALSE,
                                          save_prefix = "B3_BufferDecay",
                                          res = 300) {

  # 1. Strict Input Validation & Engine Check
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("AFRIN Error: Package 'sf' is required for geometry operations.", call. = FALSE)
  }
  if (!inherits(spatial_data, "sf")) stop("AFRIN Error: 'spatial_data' must be an 'sf' object.", call. = FALSE)
  if (!inherits(nodes, "sf")) stop("AFRIN Error: 'nodes' must be an 'sf' object.", call. = FALSE)

  # Ensure CRS matches
  crs_spatial <- sf::st_crs(spatial_data)
  crs_nodes <- sf::st_crs(nodes)
  if (crs_spatial != crs_nodes) {
    warning("AFRIN Warning: CRS mismatch detected. Automatically projecting nodes to match spatial_data.", call. = FALSE)
    nodes <- sf::st_transform(nodes, crs_spatial)
  }

  # Process events or use centroids as proxy
  has_events <- !is.null(events)
  if (has_events) {
    if (!inherits(events, "sf")) stop("AFRIN Error: 'events' must be an 'sf' object.", call. = FALSE)
    if (sf::st_crs(events) != crs_spatial) {
      events <- sf::st_transform(events, crs_spatial)
    }
    eval_pts <- events
  } else {
    message("AFRIN Note: No 'events' provided. Using spatial_data centroids as population proxies.")
    eval_pts <- suppressWarnings(sf::st_centroid(spatial_data))
  }

  # 2. Spatial Geoprocessing (Base sf Engine)
  message("AFRIN Note: Computing Spatial Buffers and Distance Matrices...")

  # Generate Concentric Buffers
  buf_primary <- sf::st_buffer(nodes, dist = buffer_dist)
  buf_half <- sf::st_buffer(nodes, dist = buffer_dist / 2)
  buf_double <- sf::st_buffer(nodes, dist = buffer_dist * 2)

  # Dissolved Primary Buffer (for coverage calculation)
  buf_dissolved <- sf::st_union(buf_primary)

  # 3. Clinical & Managerial Coverage Checks
  # Polygon Coverage (Is the polygon intersecting the buffer?)
  intersect_poly <- suppressWarnings(sf::st_intersects(spatial_data, buf_dissolved))
  covered_poly <- lengths(intersect_poly) > 0
  spatial_data[["Coverage_Status"]] <- ifelse(covered_poly, "Covered", "Blank Spot")

  # Event Coverage & Distance
  dist_matrix <- sf::st_distance(eval_pts, nodes)
  min_dists <- apply(dist_matrix, 1, min)

  eval_pts[["Dist_Nearest_Node"]] <- min_dists
  eval_pts[["Golden_Hour_Status"]] <- ifelse(min_dists <= buffer_dist, "Safe (Within Buffer)", "At Risk (Out of Reach)")

  # 4. Visualization Engine (4-Panel 4K Profiling)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  generate_plot <- function(plot_expr, suffix, title) {
    if (save_plot) {
      filename <- paste0(save_prefix, "_", timestamp, "_", suffix, ".png")
      grDevices::png(filename = filename, width = 3840, height = 2160, res = res)
    }

    old_par <- graphics::par(no.readonly = TRUE)
    oma_bottom <- ifelse(footer, 6, 2)
    graphics::par(oma = c(oma_bottom, 2, 3, 2), mar = c(4, 4, 3, 2))

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

  # PLOT 1: Concentric Buffer Zone Map
  generate_plot(expression({
    graphics::plot(sf::st_geometry(spatial_data), col = "#F5F5F5", border = "#D0D0D0")

    # Plot from largest to smallest to avoid masking
    graphics::plot(sf::st_geometry(buf_double), col = grDevices::rgb(0.9, 0.9, 0.6, 0.3), border = NA, add = TRUE)
    graphics::plot(sf::st_geometry(buf_primary), col = grDevices::rgb(0.9, 0.6, 0.3, 0.4), border = NA, add = TRUE)
    graphics::plot(sf::st_geometry(buf_half), col = grDevices::rgb(0.8, 0.2, 0.2, 0.5), border = NA, add = TRUE)
    graphics::plot(sf::st_geometry(nodes), pch = 23, col = "black", bg = "white", cex = 1.2, add = TRUE)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = c(paste("Node (Hospital)"),
                                                   paste("Radius 0.5x"),
                                                   paste("Radius 1x (", buffer_dist, ")"),
                                                   paste("Radius 2x")),
                       pt.bg = c("white", grDevices::rgb(0.8, 0.2, 0.2, 0.5),
                                 grDevices::rgb(0.9, 0.6, 0.3, 0.4), grDevices::rgb(0.9, 0.9, 0.6, 0.3)),
                       pch = c(23, 22, 22, 22), bty = "n", horiz = horiz, xpd = NA)
    }
  }), "1_ConcentricBuffers", "Concentric Healthcare Accessibility Zones")

  # PLOT 2: Event-Buffer Intersection Map (Golden Hour Risk)
  generate_plot(expression({
    graphics::plot(sf::st_geometry(spatial_data), col = "#E0E0E0", border = "white")
    graphics::plot(sf::st_geometry(buf_dissolved), col = grDevices::rgb(0.2, 0.6, 0.8, 0.3), border = "darkblue", lwd = 1.5, add = TRUE)

    # Event Colors
    pt_cols <- ifelse(eval_pts[["Golden_Hour_Status"]] == "Safe (Within Buffer)", "#27AE60", "#C0392B")
    graphics::plot(sf::st_geometry(eval_pts), pch = 20, col = pt_cols, cex = 1.2, add = TRUE)
    graphics::plot(sf::st_geometry(nodes), pch = 4, col = "black", cex = 1.5, lwd = 2, add = TRUE)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = c("Safe (Within Golden Hour)", "At Risk (Out of Reach)", "Healthcare Nodes"),
                       col = c("#27AE60", "#C0392B", "black"), pch = c(20, 20, 4),
                       bty = "n", horiz = horiz, xpd = NA)
    }
  }), "2_EventIntersection", "Patient Golden Hour Safety & Coverage Status")

  # PLOT 3: Coverage vs. Blank Spot Map (Strategic Expansion)
  generate_plot(expression({
    poly_cols <- ifelse(covered_poly, "#BDC3C7", "#34495E") # Light Grey vs Dark Blue/Grey for Blank Spots
    graphics::plot(sf::st_geometry(spatial_data), col = poly_cols, border = "white")
    graphics::plot(sf::st_geometry(nodes), pch = 21, col = "white", bg = "#E74C3C", cex = 1.5, add = TRUE)

    # Overlay an outline of the buffer to prove the logic
    graphics::plot(sf::st_geometry(buf_dissolved), border = "#E74C3C", col = NA, lty = 2, lwd = 2, add = TRUE)

    if (legend) {
      horiz <- ifelse(legend_orientation == "horizontal", TRUE, FALSE)
      graphics::legend(legend_position, legend = c("Covered Area", "Blank Spot (Target Expansion)", "Buffer Limit"),
                       fill = c("#BDC3C7", "#34495E", NA), border = c("black", "black", NA),
                       lty = c(NA, NA, 2), col = c(NA, NA, "#E74C3C"), lwd = c(NA, NA, 2),
                       bty = "n", horiz = horiz, xpd = NA)
    }
  }), "3_BlankSpotMap", "Healthcare Blank Spot & Network Expansion Zones")

  # PLOT 4: Distance Decay Curve (Statistical Proof)
  generate_plot(expression({
    # Calculate Empirical Cumulative Distribution Function (ECDF)
    dists <- as.numeric(eval_pts[["Dist_Nearest_Node"]])
    ecdf_fun <- stats::ecdf(dists)

    graphics::plot(ecdf_fun, main = "", xlab = "Distance to Nearest Healthcare Node (Units)",
                   ylab = "Cumulative Proportion of Cases/Population",
                   col = "darkblue", lwd = 2, do.points = FALSE)

    # Golden Hour Cutoff Line
    graphics::abline(v = buffer_dist, col = "red", lty = 2, lwd = 2)

    # Calculate exact % at threshold
    pct_covered <- ecdf_fun(buffer_dist) * 100
    graphics::text(x = buffer_dist, y = 0.1, labels = paste0(round(pct_covered, 1), "% Covered"),
                   col = "red", pos = 4, font = 2)

    if (legend) {
      graphics::legend("bottomright", legend = c("Distance Decay Curve", paste("Clinical Threshold (", buffer_dist, ")")),
                       col = c("darkblue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
    }
  }), "4_DistanceDecayCurve", "Distance Decay & Clinical Golden Hour Function")

  # 5. Result Compilation
  results <- list(
    spatial_coverage = spatial_data,
    event_status = eval_pts,
    buffer_polygons = buf_primary
  )

  message("AFRIN Note: Spatial Buffer and Distance Decay analysis completed successfully.")
  return(invisible(results))
}
