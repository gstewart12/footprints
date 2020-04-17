
## =============================================================================
#' Construct matrix template for calculating footprint
#'
#' @param extent A numeric vector of length four
#' @param fetch An integer, maximum length in each direction to calculate
#'   footprint probabilities. Standard is 100 * (z - zd).
#' @param res An integer, width of grid cells to be used in calculating
#'   footprint probabilities. Defaults to one (i.e. a 1-m resolution).
#'
#' @export
#'
fp_grid <- function(extent, fetch, res = 1) {
  
  dx <- res
  dy <- res
  
  if (missing(extent) & missing(fetch)) {
    stop("Must provide either extent limits or fetch.", call. = FALSE)
  } else if (!missing(extent)) {
    nx <- (extent[2] - extent[1]) / dx
    ny <- (extent[4] - extent[3]) / dy
  } else if (!missing(fetch)) {
    fetch <- floor(fetch)
    nx <- fetch * 2
    ny <- fetch * 2
    extent <- c(-nx * dx / 2, nx * dx / 2, -ny * dy / 2, ny * dy / 2)
  }
  
  # Put extent into more convenient vars
  xmin <- extent[1]
  xmax <- extent[2]
  ymin <- extent[3]
  ymax <- extent[4]
  
  # Define physical domain in cartesian coordinates
  x <- seq(xmin, xmax, length.out = nx + 1)
  y <- -seq(ymin, ymax, length.out = ny + 1)
  
  x_2d <- matrix(rep(x, each = ny + 1), nrow = ny + 1)
  y_2d <- matrix(rep(y, nx + 1), nrow = ny + 1)
  
  out <- list(x = x_2d, y = y_2d)
  attributes(out) <- list("names" = names(out), "fetch" = fetch, "res" = res)
  out
}


aoi_to_grid <- function(aoi, grid, coords) {
  
  # Add geographical coordinates to grid
  grid$y <- grid$y + coords[1] # lat
  grid$x <- grid$x + coords[2] # lon
  
  # Convert list of matrices to XYZ data frame
  temp_df <- data.frame(
    # Necessary order: lat (y), lon (x), values
    y = as.vector(grid$y),
    x = as.vector(grid$x)
  )
  temp_df$z <- 1L
  
  # Get CRS from AOI
  aoi_crs <- raster::crs(aoi)
  
  # Rasterize
  temp_rst <- raster::rasterFromXYZ(as.matrix(temp_df), crs = aoi_crs)
  
  # Mask AOI
  aoi_mask <- raster::mask(temp_rst, aoi, updatevalue = 0L)
  
  # Return matrix
  as.matrix(aoi_mask)
}


get_trim_extent <- function(x) {
  
  keep_x <- which(colSums(x) != 0)
  keep_y <- which(rowSums(x) != 0)
  
  # xmin, xmax, ymin, ymax
  c(min(keep_x), max(keep_x), min(keep_y), max(keep_y))
}


trim_matrix <- function(x, extent) {
  
  x[extent[1]:extent[2], extent[3]:extent[4]]
}


#' Rasterize a footprint matrix
#'
#' @param x A footprint object
#' @param grid A list of length two containing matrices of equal dimensions,
#'   indicating x and y coordinates. Template returned by \link{fp_grid}.
#' @param coords A numeric vector of length two
#' @param crs A character or object of class CRS.
#'
#' @export
#'
fp_rasterize <- function(x, grid, coords, crs) {
  
  # Add geographical coordinates to grid
  grid$y <- grid$y + coords[1] # lat
  grid$x <- grid$x + coords[2] # lon
  
  # Convert list of matrices to XYZ data frame
  fp_grid <- as.matrix(data.frame(
    # Necessary order: lat (y), lon (x), values
    y = as.vector(grid$y),
    x = as.vector(grid$x),
    z = as.vector(x)
  ))
  
  # Rasterize
  raster::rasterFromXYZ(fp_grid, crs = crs)
  
}


fp_polygonize <- function(x) {
  class <- class(x)
  #if (fp_is_empty(x)) return(NA)
  if (class == "RasterLayer") {
    # Method for RasterLayer class footprint
    # Set all raster values to 0 - they don't matter anymore
    homog <- raster::calc(x, function(x) x * 0)
    # Return polygon of the total area
    raster::rasterToPolygons(homog, dissolve = TRUE)
  } else if (class == "SpatialLinesDataFrame") {
    # Method for SpatialLinesDataFrame class footprint
    # Tidy ftp for creating polygon object
    fort <- broom::tidy(x)
    # Create spatial polygons object
    poly <- sp::Polygon(fort[, 1:2])
    polys <- sp::Polygons(list(sp::Polygon(poly)), "poly")
    sp::SpatialPolygons(
      list(polys), proj4string = sp::CRS("+proj=utm +zone=18 +ellps=GRS80")
    )
  }
}