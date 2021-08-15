#' Two-dimensional footprint calculation
#'
#' @param grid A list of matrices as returned by [grid_init()].
#' @param wd A numeric vector, wind direction in degrees from North.
#' @param ustar A numeric vector, friction velocity in m+1s-1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param v_sigma A numeric vector, standard deviation of the crosswind 
#'   component.
#' @param blh A numeric vector, boundary layer height in m.
#' @param z A numeric vector, measurement height in m.
#' @param zd A numeric vector, zero-plane displacement height in m.
#' @param z0 A numeric vector, aerodynamic roughness length in m.
#' @param ws A numeric vector, wind speed in m+1s-1. If NULL (the default),
#'   roughness length will be used as provided.
#' @param ... To be left empty, for flexible calls to a footprint function.
#'
#' @return A matrix with same dimensions as each item of grid.
#' 
#' @details
#' The authors assert the following criteria must be met for their model to be 
#' valid: 
#' * ustar > 0.1
#' * (zm - zd) / mo_length >= -15.5
#' * (20)zo < zm < blh
#' 
#' @references
#' Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple 
#' two-dimensional parameterisation for Flux Footprint Prediction (FFP). 
#' Geoscientific Model Development, 8(11), 3695–3713. 
#' https://doi.org/10.5194/gmd-8-3695-2015
#' 
#' @family footprint
#' @export
calc_footprint_kljun <- function(grid, wd, ustar, mo_length, v_sigma, blh, 
                                 z, zd, z0, ws = NULL, ...) {
  
  # TODO give choice of whether to choose z0 or calculate it (Eqs. 6-9)
  
  # Get grid dimensions
  dims <- dim(grid$x)
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(wd, ustar, mo_length, v_sigma, blh))) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  dw <- grid_rot$x > 0
  
  args <- list(
    x = grid_rot$x[dw], y = grid_rot$y[dw], 
    ustar = ustar, mo_length = mo_length, v_sigma = v_sigma, blh = blh, 
    z = z, zd = zd, z0 = z0
  )
  
  # If WS is given, add it to list of parameters
  if (!is.null(ws)) args <- append(args, list(ws = ws))
  
  # Subset downwind area of grid for calculating contributions
  footprint <- rlang::exec(kljun_2d, !!!args)
  
  # Check model results, return empty matrix if unexpected
  if (all(is.na(footprint)) | length(footprint) != sum(dw)) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  phi[dw] <- footprint
  
  phi
}

#' Two-dimensional footprint calculation
#'
#' @param grid A list of matrices as returned by [grid_init()].
#' @param wd A numeric vector, wind direction in degrees from North.
#' @param ustar A numeric vector, friction velocity in m+1s-1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param v_sigma A numeric vector, standard deviation of the crosswind 
#'   component.
#' @param z A numeric vector, measurement height in m,
#' @param zd A numeric vector, zero-plane displacement height in m.
#' @param z0 A numeric vector, aerodynamic roughness length in m.
#' @param ... To be left empty, for flexible calls to a footprint function.
#'
#' @return A matrix with same dimensions as each item of grid.
#' 
#' @references
#' Detto, M., Montaldo, N., Albertson, J. D., Mancini, M., & Katul, G. (2006). 
#' Soil moisture and vegetation controls on evapotranspiration in a 
#' heterogeneous Mediterranean ecosystem on Sardinia, Italy. Water Resources 
#' Research, 42(8), 1–16. https://doi.org/10.1029/2005WR004693
#' 
#' Hsieh, C.-I., Katul, G., & Chi, T. (2000). An approximate analytical model 
#' for footprint estimation of scalar fluxes in thermally stratified atmospheric 
#' flows. Advances in Water Resources, 23(7), 765–772. 
#' https://doi.org/10.1016/S0309-1708(99)00042-1
#' 
#' @family footprint
#' @export
calc_footprint_hsieh <- function(grid, wd, ustar, mo_length, v_sigma, 
                                 z, zd, z0, ...) {
  
  # Get grid dimensions
  dims <- dim(grid$x)
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(wd, ustar, mo_length, v_sigma))) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  dw <- grid_rot$x > 0
  
  # Execute model
  footprint <- hsieh_2d(
    grid_rot$x[dw], grid_rot$y[dw], ustar, mo_length, v_sigma, z, zd, z0
  )
  
  # Check model results, return empty matrix if unexpected
  if (all(is.na(footprint)) | length(footprint) != sum(dw)) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  phi[dw] <- footprint
  
  phi
}

#' Two-dimensional footprint calculation
#'
#' @param grid A list of matrices as returned by [grid_init()].
#' @param wd A numeric vector, wind direction in degrees from North.
#' @param ws A numeric vector, wind speed in m+1s-1.
#' @param ustar A numeric vector, friction velocity in m+1s-1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param v_sigma A numeric vector, standard deviation of the crosswind 
#'   component.
#' @param z A numeric vector, measurement height in m.
#' @param zd A numeric vector, zero-plane displacement height in m.
#' @param approach Either "analytical" (the default) or "numerical". Corresponds
#'   to the two different approaches for relating the power law to Monin-Obukhov
#'   similarity as described in Kormann & Meixner (2001). The analytical
#'   approach is more computationally efficient, while according to the authors
#'   the numerical approach produces better results.
#' @param ... To be left empty, for flexible calls to a footprint function.
#'
#' @return A matrix with same dimensions as each item of grid
#' 
#' @references
#' Kormann, R., & Meixner, F. X. (2001). An Analytical Footprint Model For 
#' Non-Neutral Stratification. Boundary-Layer Meteorology, 99(2), 207–224. 
#' <https://doi.org/10.1023/A:1018991015119>
#' 
#' @family footprint
#' @export
calc_footprint_kormann <- function(grid, wd, ws, ustar, mo_length, v_sigma, z,
                                   zd, approach = c("analytical", "numerical"), 
                                   ...) {
  
  approach <- match.arg(approach)
  
  # Get grid dimensions
  dims <- dim(grid$x)
  
  # Check required variables/conditions
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(ws, wd, ustar, mo_length, v_sigma))) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  dw <- grid_rot$x > 0
  
  # Execute model
  footprint <- kormann_2d(
    grid_rot$x[dw], grid_rot$y[dw], ws, ustar, mo_length, v_sigma, z, zd, 
    approach
  )
  
  # Check model results, return empty matrix if unexpected
  if (all(is.na(footprint)) | length(footprint) != sum(dw)) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  phi[dw] <- footprint
  
  phi
}

#' One-dimensional fetch calculation
#'
#' @param pct Vector of percent contributions for which to calculate fetch 
#'   distances. All values must be < 1.
#' @param ustar A numeric vector, friction velocity in m+1s-1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param blh A numeric vector, boundary layer height in m.
#' @param z A numeric vector, measurement height in m.
#' @param zd A numeric vector, zero-plane displacement height in m.
#' @param z0 A numeric vector, aerodynamic roughness length in m.
#' @param ws A numeric vector, wind speed in m+1s-1.
#' @param max_x A numeric value of length 1, maximum distance to consider
#' @param dx A numeric value of length 1, horizontal resolution
#' @param include_offset Logical, should offset distance be included in output?
#' @param include_peak Logical, should peak distance be included in output?
#' @param ... To be left empty, for flexible calls to a footprint function.
#'
#' @return
#' 
#' @family fetch
#' @export
calc_fetch_kljun <- function(pct, ustar, mo_length, blh, z, zd, z0, ws = NULL, 
                             max_x = 1000, dx = 1, include_offset = TRUE, 
                             include_peak = TRUE, ...) {
  
  if (any(pct >= 1)) {
    stop("All 'pct' values must be < 1.", call. = FALSE)
  }
  
  x <- seq(dx, max_x, by = dx)
  
  args <- list(
    x = x, ustar = ustar, mo_length = mo_length, blh = blh, 
    z = z, zd = zd, z0 = z0
  )
  
  # If WS is given, add it to list of parameters
  if (!is.null(ws)) args <- append(args, list(ws = ws))
  
  # Calculate crosswind-integrated footprint along x values
  ci <- rlang::exec(kljun_1d, !!!args)
  
  # Calculate cumulative footprint weight (i.e. integrate the curve)
  f_cume <- cumsum(ci$f) * dx
  
  #x_pct <- lapply(pct, function(.x) x[f_cume > .x][1])
  #x_pct <- lapply(x_pct, function(.x) if (.x == 100) NA else .x)
  
  # Find first x value where f_cume is greater than pct (fetch distance)
  x_pct <- purrr::map_dbl(pct, ~ x[f_cume > .x][1])
  # Remove fetch distances at upper integration limit
  x_pct[x_pct == max_x] <- NA
  
  x_pct <- as.list(x_pct)
  names(x_pct) <- paste0("fetch_", pct * 100, "_k15")
  
  if (include_offset) {
    # Add offset distance (location of 1% contribution)
    x_pct[["fetch_offset_k15"]] <- x[f_cume > 0.01][1]
  }
  
  if (include_peak) {
    # Add peak distance
    x_pct[["fetch_peak_k15"]] <- ci$x_max
  }
  
  x_pct
}

#' One-dimensional fetch calculation
#' 
#' @param pct Vector of percent contributions for which to calculate fetch 
#'   distances. All values must be < 1.
#' @param ws A numeric vector, wind speed in m+1s-1.
#' @param ustar A numeric vector, friction velocity in m+1s-1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param z A numeric vector, measurement height in m.
#' @param zd A numeric vector, zero-plane displacement height in m.
#' @param z0 A numeric vector, aerodynamic roughness length in m.
#' @param approach Either "analytical" (the default) or "numerical". Corresponds
#'   to the two different approaches for relating the power law to Monin-Obukhov
#'   similarity as described in Kormann & Meixner (2001). According to the 
#'   authors, the analytical approach is more computationally efficient while
#'   the numerical approach produces better results.
#' @param max_x A numeric value of length 1, maximum distance to consider
#' @param dx A numeric value of length 1, horizontal resolution
#' @param include_offset Logical, should offset distance be included in output?
#' @param include_peak Logical, should peak distance be included in output?
#' @param ... To be left empty, for flexible calls to a footprint function.
#'
#' @return
#' 
#' @family fetch
#' @export
calc_fetch_kormann <- function(pct, ws, ustar, mo_length, z, zd, z0,
                               approach = c("analytical", "numerical"),
                               max_x = 1000, dx = 1, include_offset = TRUE, 
                               include_peak = TRUE, ...) {
  
  approach <- match.arg(approach)
  
  if (any(pct >= 1)) {
    stop("All 'pct' values must be < 1.", call. = FALSE)
  }
  
  x <- seq(dx, max_x, by = dx)
  
  args <- list(
    x = x, ws = ws, ustar = ustar, mo_length = mo_length,
    z = z, zd = zd, z0 = z0
  )
  
  # Calculate crosswind-integrated footprint along x values
  ci_fun <- if (approach == "numerical") kormann_1d_num else kormann_1d
  ci <- rlang::exec(ci_fun, !!!args)
  
  # Calculate cumulative footprint weight (i.e. integrate the curve)
  f_cume <- cumsum(ci$f) * dx
  
  #x_pct <- lapply(pct, function(.x) x[f_cume > .x][1])
  #x_pct <- lapply(x_pct, function(.x) if (.x == 100) NA else .x)
  
  # Find first x value where f_cume is greater than pct (fetch distance)
  x_pct <- purrr::map_dbl(pct, ~ x[f_cume > .x][1])
  # Remove fetch distances at upper integration limit
  x_pct[x_pct == max_x] <- NA
  
  x_pct <- as.list(x_pct)
  names(x_pct) <- paste0("fetch_", pct * 100, "_km01")
  
  if (include_offset) {
    # Add offset distance (location of 1% contribution)
    x_pct[["fetch_offset_km01"]] <- x[f_cume > 0.01][1]
  }
  
  if (include_peak) {
    # Add peak distance
    x_pct[["fetch_peak_km01"]] <- ci$x_max
  }
  
  x_pct
}

#' One-dimensional fetch calculation
#' 
#' @param pct Vector of percent contributions for which to calculate fetch 
#'   distances. All values must be < 1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param z A numeric vector, measurement height in m.
#' @param zd A numeric vector, zero-plane displacement height in m.
#' @param z0 A numeric vector, aerodynamic roughness length in m.
#' @param max_x A numeric value of length 1, maximum distance to consider
#' @param dx A numeric value of length 1, horizontal resolution
#' @param include_offset Logical, should offset distance be included in output?
#' @param include_peak Logical, should peak distance be included in output?
#' @param ... To be left empty, for flexible calls to a footprint function.
#'
#' @return
#' 
#' @family fetch
#' @export
calc_fetch_hsieh <- function(pct, mo_length, z, zd, z0, max_x = 1000, dx = 1, 
                             include_offset = TRUE, include_peak = TRUE, ...) {
  
  if (any(pct >= 1)) {
    stop("All 'pct' values must be < 1.", call. = FALSE)
  }
  
  x <- seq(dx, max_x, by = dx)
  
  args <- list(
    x = x, mo_length = mo_length, z = z, zd = zd, z0 = z0
  )
  
  # Calculate crosswind-integrated footprint along x values
  ci <- rlang::exec(hsieh_1d, !!!args)
  
  # Calculate cumulative footprint weight (i.e. integrate the curve)
  f_cume <- cumsum(ci$f) * dx
  
  #x_pct <- lapply(pct, function(.x) x[f_cume > .x][1])
  #x_pct <- lapply(x_pct, function(.x) if (.x == 100) NA else .x)
  
  # Find first x value where f_cume is greater than pct (fetch distance)
  x_pct <- purrr::map_dbl(pct, ~ x[f_cume > .x][1])
  # Remove fetch distances at upper integration limit
  x_pct[x_pct == max_x] <- NA
  
  x_pct <- as.list(x_pct)
  names(x_pct) <- paste0("fetch_", pct * 100, "_h00")
  
  if (include_offset) {
    # Add offset distance (location of 1% contribution)
    x_pct[["fetch_offset_h00"]] <- x[f_cume > 0.01][1]
  }
  
  if (include_peak) {
    # Add peak distance
    x_pct[["fetch_peak_h00"]] <- ci$x_max
  }
  
  x_pct
}
