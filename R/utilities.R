
# Performs vector operations on matrices while retaining matrix dimensions
with_matrix <- function(x, .f) {
  
  dims <- dim(x)
  mat_fun <- rlang::as_function(.f)
  
  vec <- mat_fun(as.vector(x))
  
  matrix(vec, nrow = dims[1], ncol = dims[2])
}


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
grid_init <- function(extent, fetch, res = 1) {
  
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


snap_to_grid <- function(x, grid, coords, resample_method = "ngb") {
  
  template <- grid_template_rst(grid, coords)
  
  # Assign CRS from original object
  raster::crs(template) <- raster::crs(x)
  
  # Resample according to grid parameters
  resampled <- raster::resample(x, template, method = resample_method)
  
  out <- raster::as.matrix(resampled)
  
  # Mask AOI if exists in grid list (maybe not, easier if there are no NAs)
  if ("aoi" %in% names(grid)) {
    aoi_mask <- grid$aoi
    aoi_mask[which(aoi_mask == 0)] <- NA
    #out <- out * aoi_mask
  } 
  
  # Return matrix object fit to grid
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
  raster::as.matrix(aoi_mask)
}


get_trim_extent <- function(x, values = 0) {
  
  # Only support for one trim value at the moment
  keep_y <- which(rowSums(x) != values)
  keep_x <- which(colSums(x) != values)
  
  # ymin, ymax, xmin, xmax
  c(min(keep_y), max(keep_y), min(keep_x), max(keep_x))
}


trim_matrix <- function(x, extent = get_trim_extent(x)) {
  
  x[extent[1]:extent[2], extent[3]:extent[4]]
}


plot_matrix <- function(x) {
  
  if (storage.mode(x) == "character") {
    x <- with_matrix(x, ~ unclass(as.factor(.x)))
  } 
  
  raster::plot(raster::raster(x))
}


#' @importFrom tidyr "pivot_longer" tibble "as_tibble"
pivot_matrix <- function(x) {
  
  if (storage.mode(x) == "character") {
    x <- with_matrix(x, ~ unclass(as.factor(.x)))
  }
  
  dims <- dim(x)
  if (is.null(dimnames(x))) dimnames(x) <- list(rev(1:dims[1]), 1:dims[2])
  
  x %>%
    tibble::as_tibble(rownames = "y") %>%
    tidyr::pivot_longer(-y, names_to = "x", values_to = "z") %>%
    dplyr::mutate(dplyr::across(is.character, as.numeric))
}


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


#' @importFrom readr "write_delim"
#' @importFrom tibble "as_tibble"
write_matrix <- function(x, file, trunc = 9, compress = TRUE) {
  
  # Convert data type to integer for more efficient storage
  # - this means weights below (1 / 'trunc') are effectively zero
  if (!is.na(trunc)) {
    x <- trunc(x * 10^trunc)
    storage.mode(x) <- "integer"
  }
  
  # Set path
  file_path <- paste0(file, ".txt")
  if (compress) file_path <- paste0(file_path, ".gz")
  
  # Convert to tbl for writing
  tbl <- tibble::as_tibble(x, .name_repair = "minimal")
  
  # Write to file
  readr::write_delim(tbl, file_path, col_names = FALSE)
}


#' @importFrom readr "read_table2" "cols" "col_integer" "col_double"
read_matrix <- function(file, trunc = 9) {
  
  if (!is.na(trunc)) {
    scale <- 10^trunc
    data_type <- readr::col_integer()
  } else {
    scale <- 1
    data_type <- readr::col_double()
  }
  
  data <- readr::read_table2(
    file, col_names = FALSE, col_types = readr::cols(.default = data_type),
    progress = FALSE
  )
  #data <- read.table(file, sep = " ")
  data_mat <- as.matrix(data)
  dimnames(data_mat) <- NULL
  
  # Scale based on how data was stored
  data_mat <- data_mat / scale
  
  data_mat
}


#' Read grid coordinates from two .txt files
#' 
#' @param path A character vector containing a single path
#' @param names A character vector of length two, the file names of grid 
#'   elements (omitting the extension)
#' 
#' @importFrom magrittr "%>%"
#' 
#' @export
read_grid <- function(path, names = c("x", "y")) {
  
  file_names <- paste0(names, ".txt")
  grid_files <- file.path(path, file_names)
  
  grid_files %>% 
    purrr::map(read_matrix, trunc = NA) %>%
    rlang::set_names(names)
}


pluck_cell <- function(x, grid, coords) {
  
  cell_loc <- which(grid$x == coords[1] & grid$y == coords[2])
  x[cell_loc]
}


#' Rasterize a footprint matrix
#'
#' @param x A footprint object
#' @param grid A list of length two containing matrices of equal dimensions,
#'   indicating x and y coordinates. Template returned by \link{grid_init}.
#' @param coords A numeric vector of length two
#' @param crs A character or object of class CRS.
#'
#' @export
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


rotate_grid <- function(grid, dir) {
  #browser()
  # Calculate wind direction angle
  #theta <- ((360 - dir) %% 360) * (pi / 180)
  
  #out <- grid
  
  # Rotate coordinates toward wind direction
  #out$x <- grid$x * cos(theta) - grid$y * sin(theta)
  #out$y <- grid$x * sin(theta) + grid$y * cos(theta)
  
  out <- grid
  
  # Distance and direction of each grid cell
  fpd <- sqrt(grid$y^2 + grid$x^2)
  fpa <- atan2(grid$x, grid$y) * 180 / pi - dir
  
  # Rotate original coordinates using grid vectors
  out$x <- cos(fpa * pi / 180) * fpd
  out$y <- -1 * sin(fpa * pi / 180) * fpd
  
  out
}


#' Mask footprint to its analytical source area
#'
#' As the analytical footprint approaches 100%, the coverage increases rapidly
#' over area that contributes little to the total measured flux (Matthes et al.,
#' 2014). Thus, it may be desirable to consider a percantage of the calculated
#' footprint as the analytical footprint. This is acheived by masking the
#' footprint to a given maximum cumulative weight.
#'
#' @param x The footprint as a numeric matrix or RasterLayer. 
#' @param p A numeric value between 0 and 1, recommended by Kljun et al. (2015)
#'   to be less than 0.9. Defaults to 0.85.
#' @param mask_value A value to assign cells beyond the maximum percent weight.
#'   Defaults to NA, but 0 may be useful in some situations.
#'
#' @return An object of the same class as x, either a numeric matrix or 
#'   RasterLayer.
#' 
#' @export
mask_source_area <- function(x, p = 0.85, mask_value = NA) {
  
  # Check inputs
  
  # Assume that p is given as a percentage if >1
  if (p > 1) p <- p / 100
  
  # Sort matrix values from high to low
  rs <- sort(as.vector(x), decreasing = TRUE)
  
  # Calculate cumulative sums of footprint weights
  cs <- cumsum(rs)
  
  # Get the minimum value of cumulative sum less than p
  rp <- min(rs[cs < p])
  out <- x
  
  # Set all other cells to fill value
  out[out < rp] <- mask_value
  
  out
  
}


#' Convert footprint to a polygon representing a given source area
#'
#' @param x The footprint as a numeric matrix or RasterLayer. 
#' @param p A numeric value between 0 and 1, recommended by Kljun et al. (2015)
#'   to be less than 0.9. Defaults to 0.85.
#'
#' @export
source_area_polygon <- function(x, p = 0.85) {
  
  masked <- mask_source_area(x, p = p, mask_value = NA)
  
  class <- class(x)
  
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