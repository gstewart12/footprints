
plot_footprint <- function(x, type = c("weights", "cumulative"), 
                           colors = c("greyscale")) {
  
  # good example: calc_footprint_kormann(data[8000, ])
  
  data <- x %>%
    accumulate_weights(max = 1, zero_as_max = TRUE) %>% 
    pivot_matrix() 
  
  data %>% 
    ggplot2::ggplot(ggplot2::aes(x, y, z = z)) + 
    ggplot2::geom_contour_filled(
      na.rm = TRUE, breaks = seq(0, 0.9, by = 0.1), alpha = 0.6
    ) + 
    ggplot2::coord_fixed() + 
    #ggplot2::scale_fill_brewer(palette = "Greys", direction = 1) +
    ggplot2::scale_fill_grey(start = 0.95, end = 0.05) + 
    ggplot2::theme_minimal() + 
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    #ggplot2::guides(fill = guide_bins(show.limits = TRUE)) +
    ggplot2::guides(fill = guide_colorsteps(show.limits = TRUE))
}


#' Fast plotting of a generic matrix
#'
#' @param x A matrix.
#' 
#' @importFrom raster "raster" "plot"
#' 
#' @export
plot_matrix <- function(x) {
  
  if (storage.mode(x) == "character") {
    x <- with_matrix(x, ~ unclass(as.factor(.x)))
  } 
  
  raster::plot(raster::raster(x))
}


#' Plot a generic matrix using ggplot
#'
#' @param x A matrix.
#' @param trans Name of a transformation object as defined in 
#'   ggplot2::continuous_scale. "sqrt" and "log" are particularly useful.
#'
#' @return A ggplot object.
#' 
#' @export
plot_tidy_matrix <- function(x, trans = NA) {
  
  # trans = "sqrt" best for individual footprints
  # trans = "log" best for footprint topology
  if (is.na(trans)) trans <- "identity"
  
  data <- pivot_matrix(x)
  
  data %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_distiller(palette = "Spectral", trans = trans) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}


plot_windrose <- function(data, ws = ws, wd = wd, ws_res = 0.5, wd_res = 30,
                          ws_min = 0, ws_max = 4, ws_seq = NULL,
                          palette = "Spectral", labels = c("range", "center"), 
                          countmax = NA) {
  
  ws <- data %>% pull(!!rlang::enquo(ws))
  wd <- data %>% pull(!!rlang::enquo(wd))
  labels <- rlang::arg_match(labels)
  
  # wd_bin <- ggplot2::cut_number(wd, 360 / wd_res)
  # ws_bin <- ggplot2::cut_interval(ws, length = ws_res)
  
  # Tidy up input data ----
  n.in <- NROW(data)
  dnu <- (is.na(ws) | is.na(wd))
  ws[dnu] <- NA
  wd[dnu] <- NA
  
  # figure out the wind speed bins
  if (missing(ws_seq)) {
    ws_seq <- seq(ws_min, ws_max, ws_res)
  }
  
  # get some information about the number of bins, etc.
  n_ws_seq <- length(ws_seq)
  n_colors <- n_ws_seq - 1
  
  # create the color map
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
    min(max(3, n_colors), min(9, n_colors)), palette
  )))(n_colors)
  
  if (max(ws, na.rm = TRUE) > ws_max) {
    ws_breaks <- c(ws_seq, max(ws, na.rm = TRUE))
    ws_labs <- c(
      paste(c(ws_seq[1:n_ws_seq - 1]), "-", c(ws_seq[2:n_ws_seq])),
      paste(ws_max, "-", round(max(ws, na.rm = TRUE), 1))
    )
    colors <- c(colors, "grey50")
  } else {
    ws_breaks <- ws_seq
    ws_labs <- paste(c(ws_seq[1:n_ws_seq - 1]), "-", c(ws_seq[2:n_ws_seq]))
  }
  
  ws_bin <- cut(
    x = ws,
    breaks = ws_breaks,
    labels = ws_labs,
    ordered_result = TRUE
  )
  
  # figure out the wind direction bins
  wd_breaks <- c(
    -wd_res / 2, seq(wd_res / 2, 360 - wd_res / 2, by = wd_res),
    360 + wd_res / 2
  )
  
  if (labels == "range") {
    wd_labs <- c(
      paste(360 - wd_res / 2, "-", wd_res / 2),
      paste(
        seq(wd_res / 2, 360 - 3 * wd_res / 2, by = wd_res), "-",
        seq(3 * wd_res / 2, 360 - wd_res / 2, by = wd_res)
      ),
      paste(360 - wd_res / 2, "-", wd_res / 2)
    )
  } else if (labels == "center") {
    wd_labs <- c(seq(0, 360 - wd_res, by = wd_res), 0)
  }
  
  # assign each wind direction to a bin
  wd_bin <- cut(wd, breaks = wd_breaks, ordered_result = TRUE)
  levels(wd_bin) <- wd_labs
  
  # deal with change in ordering introduced somewhere around version 2.2
  if (packageVersion("ggplot2") > "2.2") {
    ws_bin <- factor(ws_bin, levels = rev(levels(ws_bin)))
    colors <- rev(colors)
  }
  
  # Put everything in a data frame
  df <- data.frame(wd_bin = wd_bin, ws_bin = ws_bin)
  
  # create the plot
  out <- df %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(wd_bin, fill = ws_bin)) +
    ggplot2::geom_bar() +
    ggplot2::scale_x_discrete(drop = FALSE, labels = waiver()) +
    ggplot2::coord_polar(start = -((wd_res / 2) / 360) * 2 * pi) +
    ggplot2::scale_fill_brewer(
      name = "Wind speed (m/s)", palette = "Spectral", drop = FALSE
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90")
    )
  
  # adjust axes if required
  if (!is.na(countmax)) {
    out <- out + ggplot2::ylim(c(0, countmax))
  }
  
  out
}