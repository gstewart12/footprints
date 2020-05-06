
check_percents <- function(pct) {
  
  # Convert to numeric, drop any non-numeric
  out <- na.omit(purrr::quietly(as.numeric)(pct)$result)
  attributes(out) <- NULL
  
  # Throw error if any out of range
  if (any(out > 100)) {
    stop("'pct' cannot be higher than 100.", call. = FALSE)
  }
  
  # Truncate excessive precision
  out <- signif(out, 2)
  
  # Ensure decimal format
  out <- pmax(out %% 1, out / 100)
  
  # Return sorted decimal percents
  sort(out)
}


grid_template_rst <- function(grid, coords) {
  
  # lon range, lat range (not positive about the order)
  grid_range <- c(range(grid$y + coords[1]), range(grid$x + coords[2]))
  extent <- raster::extent(grid_range)
  
  # Return raster template
  raster::raster(extent, nrows = nrow(grid$x), ncols = ncol(grid$x))
}