
calc_fp_kljun <- function(grid, wd, ustar, mo_length, v_sigma, blh, z, zd, zo) {
  
  # Get grid dimensions
  n <- dim(grid$x)[1]
  
  # Check required variables/conditions
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(wd, ustar, mo_length, v_sigma, blh))) {
    return(matrix(NA, nrow = n, ncol = n))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = n, ncol = n)
  
  # Subset downwind area of grid for calculating contributions
  phi[grid_rot$x > 0] <- kljun_model(
    grid_rot$x, grid_rot$y, ustar, mo_length, v_sigma, blh, z, zd, zo, n
  )
  
  phi
}


## =============================================================================
#' Calculate two-dimensional footprints
#'
#' Create a two-dimensional matrix of footprint probabilities.
#'
#' @param wd A numeric vector, wind direction in degrees from North.
#' @param ustar A numeric vector, friction velocity in ms+1s-1.
#' @param mo_length A numeric vector, Monin-Obukhov length in m.
#' @param v_sigma A numeric vector, standard deviation of the crosswind 
#'   component.
#' @param ws A numeric vector, wind speed in ms+1s-1.
#' @param blh A numeric value, boundary layer height in m.
#' @param z A numeric value, tower height in m.
#' @param zd A numeric value, zero-plane displacement height in m.
#' @param zo The aerodynamic roughness length in m.
#' @param grid A list of length two containing matrices of equal dimensions,
#'   indicating x and y coordinates. Template returned by \link{fp_grid}.
#' @param model A character string naming the model to be used in footprint
#'   calculations. Can be "KM01" for the Kormann & Meixner (2001) model or "H00"
#'   for the Hsieh et al. (2000) model.
#'
#' @details Support only for square grids at the moment.
#'
#' @export
#'
calc_fp <- function(grid, wd, ustar, mo_length, v_sigma, z, zd, zo, ws = NULL, 
                    blh = NULL, model = c("H00", "K15", "KM01")) {
  
  model <- match.arg(model)
  
  # Get grid dimensions
  n <- dim(grid$x)[1]
  
  # Check for equal dimensions
  #if (n != dim(grid$x)[2]) stop("Only square grids supported.", call. = FALSE)
  
  # Check required variables/conditions
  vars <- c(wd, ustar, mo_length, v_sigma)
  invalid <- FALSE
  if (model == "KM01") {
    vars <- c(ws, vars)
    # Kormann model assumptions must be met: abs(zeta) > 3
    invalid <- !dplyr::between((z - zd) / mo_length, -3, 3)
  } else invalid <- FALSE
  if (model == "K15") {
    vars <- c(blh, vars)
    # Kljun model assumptions must be met: ustar > 0.1, zm / L >= -15.5
    #invalid <- ustar >= 0.1 & ((z - zd) / mo_length) >= -15.5
  }
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(vars) | invalid) return(matrix(NA, nrow = n, ncol = n))
  
  # Set initial output grid
  phi <- matrix(0, nrow = n, ncol = n)
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Subset downwind area of grid for calculating contributions
  dw <- grid_rot$x > 0
  
  # Execute chosen model
  if (model == "KM01") {
    # Kormann model
    phi[dw] <- kormann_model(
      grid_rot$x, grid_rot$y, ws, ustar, mo_length, v_sigma, z, zd, zo, n
    )
  } else if (model == "H00") {
    # Hsieh model
    phi[dw] <- hsieh_model(
      grid_rot$x, grid_rot$y, ustar, mo_length, v_sigma, z, zd, zo, n
    )
  } else if (model == "K15") {
    # Kljun model
    phi[dw] <- kljun_model(
      grid_rot$x, grid_rot$y, ustar, mo_length, v_sigma, blh, z, zd, zo, n
    )
  }
  
  # In case of invalid output within footprint function, return entire grid NA
  #if (all(is.na(phi[dw]))) phi <- matrix(NA, nrow = n, ncol = n)
  
  # Output matrix holding footprint probabilities
  phi
}


rotate_grid <- function(grid, dir) {
  
  # Calculate wind direction angle
  theta <- ((360 - dir) %% 360) * (pi / 180)
  
  out <- grid
  
  # Rotate coordinates toward wind direction
  out$x <- grid$x * cos(theta) - grid$y * sin(theta)
  out$y <- grid$x * sin(theta) + grid$y * cos(theta)
  
  out
}

