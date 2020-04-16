
## =============================================================================
#' Calculate two-dimensional footprints
#'
#' Create a two-dimensional matrix of footprint probabilities.
#'
#' @param wd A numeric vector, wind direction in degrees from North.
#' @param ustar A numeric vector, friction velocity in ms+1s-1.
#' @param L A numeric vector, Monin-Obukhov length in m.
#' @param v_sd A numeric vector, standard deviation of the crosswind component.
#' @param ws A numeric vector, wind speed in ms+1s-1.
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
#' @return
#' @export
#'
#' @examples
fp_calculate <- function(wd, ustar, mo_length, v_sigma, ws, blh, z, zd, zo,
                         grid, model = c("H00", "K15", "KM01")) {
  
  model <- match.arg(model)
  
  # Get grid dimensions
  n <- dim(grid$x)[1]
  
  # Check for equal dimensions
  if (n != dim(grid$x)[2]) stop("Only square grids supported.", call. = FALSE)
  
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
  if (anyNA(vars) | invalid) {
    phi <- matrix(NA, nrow = n, ncol = n)
    return(phi)
  }
  
  # Set initial output grid
  phi <- matrix(0, nrow = n, ncol = n)
  
  # Calculate wind direction angle
  theta <- ((360 - wd) %% 360) * (pi / 180)
  
  # Rotate coordinates toward wind direction
  x_rot <- grid$x * cos(theta) - grid$y * sin(theta)
  y_rot <- grid$x * sin(theta) + grid$y * cos(theta)
  
  # Subset downwind area of grid for calculating contributions
  dw <- x_rot > 0
  
  # Execute chosen model
  if (model == "KM01") {
    # Kormann model
    phi[dw] <- kormann_model(
      x_rot, y_rot, ws, ustar, mo_length, v_sigma, z, zd, zo, n
    )
  } else if (model == "H00") {
    # Hsieh model
    phi[dw] <- hsieh_model(
      x_rot, y_rot, ustar, mo_length, v_sigma, z, zd, zo, n
    )
  } else if (model == "K15") {
    # Hsieh model
    phi[dw] <- kljun_model(
      x_rot, y_rot, ustar, mo_length, v_sigma, blh, z, zd, zo, n
    )
  }
  
  # In case of invalid output within footprint function, return entire grid NA
  #if (all(is.na(phi[dw]))) phi <- matrix(NA, nrow = n, ncol = n)
  
  # Output matrix holding footprint probabilities
  phi
  
}