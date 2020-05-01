
calc_footprint_kljun <- function(grid, wd, ustar, mo_length, v_sigma, blh, 
                                 z, zd, zo, ...) {
  
  # Get grid dimensions
  dims <- dim(grid$x)
  
  # Check required variables/conditions
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(wd, ustar, mo_length, v_sigma, blh))) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  footprint <- kljun_model(
    grid_rot$x, grid_rot$y, ustar, mo_length, v_sigma, blh, z, zd, zo, dims
  )
  
  # Check model results, return empty matrix if unexpected
  if (length(footprint) != prod(dims)) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  phi[grid_rot$x > 0] <- footprint
  
  phi
}

calc_footprint_hsieh <- function(grid, wd, ustar, mo_length, v_sigma, 
                                 z, zd, zo, ...) {
  
  # Get grid dimensions
  dims <- dim(grid$x)
  
  # Check required variables/conditions
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(wd, ustar, mo_length, v_sigma))) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  footprint <- hsieh_model(
    grid_rot$x, grid_rot$y, ustar, mo_length, v_sigma, z, zd, zo, dims
  )
  
  # Check model results, return empty matrix if unexpected
  if (length(footprint) != prod(dims)) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  phi[grid_rot$x > 0] <- footprint
  
  phi
}

calc_footprint_kormann <- function(grid, ws, wd, ustar, mo_length, v_sigma, 
                                   z, zd, zo, ...) {
  
  # Get grid dimensions
  dims <- dim(grid$x)
  
  # Check required variables/conditions
  
  # Exit gracefully for invalid cases (return empty matrix)
  if (anyNA(c(wd, ustar, mo_length, v_sigma))) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  # Rotate grid toward wind direction
  grid_rot <- rotate_grid(grid, wd)
  
  # Initialize output grid
  phi <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  footprint <- kormann_model(
    grid_rot$x, grid_rot$y, ws, ustar, mo_length, v_sigma, z, zd, zo, dims
  )
  
  # Check model results, return empty matrix if unexpected
  if (length(footprint) != prod(dims)) {
    return(matrix(NA, nrow = dims[1], ncol = dims[2]))
  } 
  
  phi[grid_rot$x > 0] <- footprint
  
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
#'   indicating x and y coordinates. Template returned by \link{grid_init}.
#' @param model A character string naming the model to be used in footprint
#'   calculations. Can be "KM01" for the Kormann & Meixner (2001) model or "H00"
#'   for the Hsieh et al. (2000) model.
#'
#' @details Support only for square grids at the moment.
#'
#' @export
#'
calc_footprint <- function(grid, wd, ustar, mo_length, v_sigma, z, zd, zo, 
                           ws = NULL, blh = NULL, 
                           model = c("H00", "K15", "KM01")) {
  
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
  if (anyNA(vars) | invalid) {
    return(matrix(NA, nrow = nrow(grid$x), ncol = ncol(grid$x)))
  } 
  
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


## =============================================================================
#' Calculate Fetch Distances for a Dataset
#'
#' @description
#'
#' @param data
#' @param z
#' @param zd
#' @param roughness_length
#' @param percent
#' @param method
#'
#' @return
#' @export
#'
#' @examples
fp_fetch <- function(data, ws = "ws", ustar = "ustar", zeta = "zeta",
                     mo_length = "mo_length", z, zd, zo, percent,
                     method = c("KM01", "H00", "K04")) {
  
  # Simplify variables names
  x <- rep(NA, nrow(data))
  p <- sort(percent)
  method <- match.arg(method)
  
  pbar <- dplyr::progress_estimated(length(x))
  if (method == "KM01") {
    vars <- c(ws, ustar, zeta)
    if (!all(vars %in% names(data))) {
      stop("Missing ", paste0(
        vars[!(vars %in% names(data))], collapse = ", "
      ), ".", call. = FALSE)
    }
    out <- data[, 0]
    names <- dplyr::if_else(
      !suppressWarnings(is.na(as.numeric(p))),
      paste0("x_", p, "perc"), paste0("x_", p)
    )
    out[, names] <- NA
    for (i in seq_along(x)) {
      # Set model inputs - Kormann & Meixner 2001
      ws_ <- data[i, ws]
      ustar_ <- data[i, ustar]
      zeta_ <- data[i, zeta]
      # Execute model function
      out[i, ] <- kormann_ftp_dist(ws_, ustar_, zeta_, z, zd, p)
      pbar$tick()$print()
    }
  } else if (method == "H00") {
    if (!mo_length %in% names(data)) stop("Missing mo_length.", call. = FALSE)
    out <- data[, 0]
    names <- dplyr::if_else(
      !suppressWarnings(is.na(as.numeric(p))),
      paste0("x_", p, "perc"), paste0("x_", p)
    )
    names <- c("xoffset", names)
    out[, names] <- NA
    for (i in seq_along(x)) {
      # Set model inputs - Hsieh 2000
      l_ <- data[i, mo_length]
      # Execute model function
      out[i, ] <- hsieh_ftp_dist(z, zd, zo, l_, p)
      pbar$tick()$print()
    }
  }
  # Return vector if only one p was requested
  if (ncol(out) < 2) out <- out[, 1]
  out
}


#' Time-averaged footprint
#'
#' @param x A list of footprint matrices.
#' @param grid A list of length two containing matrices of equal dimensions,
#'   indicating x and y coordinates. Template returned by \link{grid_init}.
#' @param weights A vector of same length as x
#'
#' @export
calc_climatology <- function(x, grid, weights = NULL) {
  
  # Get grid dimensions
  n <- dim(grid$x)[1]
  len <- length(x)
  
  if (!is.null(weights)) {
    avg_ftp <- list(
      z = matrix(0, nrow = n, ncol = n),
      fpw = matrix(0, nrow = n, ncol = n)
    )
  } else avg_ftp <- list(z = matrix(0, nrow = n, ncol = n))
  
  count <- 0
  #pbar <- dplyr::progress_estimated(len)
  
  for (i in 1:len) {
    
    fp <- x[[i]]
    
    if (any(!is.na(fp$z))) {
      if (!is.null(weights)) {
        # If using weights, only include record if both fp and weight exist
        if (!is.na(weights[i])) {
          avg_ftp$z <- avg_ftp$z + fp$z
          avg_ftp$fpw <- avg_ftp$fpw + fp$z * weights[i]
          count <- count + 1
        }
      } else if (is.null(weights)) {
        # Only need fp if not using weights
        avg_ftp$z <- avg_ftp$z + fp$z
        count <- count + 1
      }
    }
    
    #pbar$tick()$print()
  }
  
  avg_ftp$z <- avg_ftp$z / count
  
  if (!is.null(weights)) {
    # Similar weighting routine to Budishchev et al. (2014)
    wsum <- sum(avg_ftp$fpw, na.rm = TRUE) # equals sum of weights
    avg_ftp$z <- avg_ftp$fpw / wsum
    #avg_ftp$fpw <- avg_ftp$fpw / count # this preserves the weights units
  }
  
  avg_ftp$z
}


summarize_cover <- function(x, y, type = c("factor", "numeric"), levels, 
                            ignore_levels = 0) {
  
  type <- rlang::arg_match(type)
  
  if (type == "numeric") {
    
    # Convolve matrices, normalize by total footprint weight
    #out <- sum((x / sum(x)) * y)
    out <- sum(x * y) / sum(x)
    
    return(out)
  }
  
  # Detect levels if necessary
  if (missing(levels)) {
    levels <- sort(na.omit(unique(as.vector(y))))
    if (length(levels) > length(x) / 2) {
      stop("Too many levels detected. Is type = 'numeric' more appropriate?")
    }
  } 
  
  # Ignore level(s) if indicated
  levels <- levels[!levels %in% ignore_levels]
  
  # Add names if not given
  if (!rlang::is_named(levels)) levels <- rlang::set_names(levels)
  
  n <- length(levels)
  out <- rlang::set_names(rep(NA, n), names(levels))
  
  # Add up weights for each cover type
  for (i in 1:n) {
    cells <- which(y == levels[i])
    if (length(cells) == 0) out[i] <- 0 else out[i] <- sum(x[cells])
  }
  
  # Return named vector equal in length to number of level
  out
}
