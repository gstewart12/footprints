
## =============================================================================
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
#' @family footprints
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


## =============================================================================
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
#' @family footprints
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


## =============================================================================
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
#' @family footprints
#' @export
calc_footprint_kormann <- function(grid, wd, ws, ustar, mo_length, v_sigma, z,
                                   zd, approach = c("analytical", "numerical"), 
                                   ...) {
  
  approach <- rlang::arg_match(approach)
  
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


calc_fetch_kljun <- function(pct, ustar, mo_length, blh, z, zd, z0, ws = NULL, 
                             max_x = 1000, dx = 1, offset = 0.01, peak = TRUE, 
                             ...) {
  
  if (any(pct >= 1)) {
    stop("All 'pct' values must be < 1.0.", call. = FALSE)
  }
  #browser()
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
  
  if (!is.null(offset)) {
    # Add offset distance (location of 1% contribution)
    x_pct[["fetch_offset_k15"]] <- x[f_cume > offset][1]
  }
  
  if (peak) {
    # Add peak distance
    x_pct[["fetch_peak_k15"]] <- ci$x_max
  }
  
  x_pct
}


calc_fetch_kormann <- function(pct, ws, ustar, mo_length, zm, zd, ...) {
  
  # Calculate x if flux percents are given
  where <- pct
  pct <- na.omit(purrr::quietly(as.numeric)(pct)$result)
  pct <- sort(pct)
  
  if (length(pct) > 0) {
    # Define footprint model as a function of x (distance from tower)
    f <- function(x) xi^mu * exp(-xi / x) / (x^(1 + mu) * gamma(mu))
    p <- sort(p)
    fx <- rep(NA, length(p))
    x <- 1
    while (x < 10000) {
      fp <- try(integrate(f, 0, x)$value, silent = TRUE)
      if (class(fp) == "try-error") {
        # Certain (uncommon) meteorological conditions result in a non-finite
        # integral of the footprint function: set these records to NA
        fp <- NA
        break
      }
      fx[which(((fp > (p / 100) & is.na(fx))) == TRUE)[1]] <- x
      if (!anyNA(fx)) break # finished if all fetch lengths are found
      x <- x + 1
    }
  }
  if ("peak" %in% all_p) {
    peak <- xi / (1 + mu)
    fx <- c(fx, peak)
  }
}


calc_fetch_hsieh <- function(pct = c("offset", "peak", 30, 50, 70, 80, 90), 
                             mo_length, zm, zd, zo, ...) {
  
  # Fail gracefully if missing values
  if (anyNA(c(mo_length, zm, zd, zo))) {
    return(rep(NA, length(pct)))
  } 
  
  # Prepare decimal percents at which distances should be found
  pos <- pct
  pct <- check_percents(pct)
  if ("offset" %in% pos) pct <- c(0.01, pct)
  
  # Calculate x if flux percents are given
  if (length(pct) > 0) {
    
    min_x <- (mu / 100) * zo
    max_x <- mu * zm
    x_bin <- (max_x - min_x) / 1000
    
    pct <- sort(pct)
    fx <- rep(NA, length(pct))
    
    x <- min_x
    
    while (x < max_x) {
      
      fp <- hsieh_1d(x, mo_length, zm, zd, zo)$Fy
      fx[which(((fp > pct & is.na(fx))) == TRUE)[1]] <- x
      
      if (!anyNA(fx)) break # finished if all fetch lengths are found
      x <- x + x_bin
    }
  }
  
  # Calculate x if peak distance is requested
  if ("peak" %in% pct) {
    peak <- DP / (2 * k^2)
    fx <- c(fx, peak)
  }
}


fp_fetch <- function(data, ws = "ws", ustar = "ustar", zeta = "zeta",
                     mo_length = "mo_length", z, zd, zo, percent,
                     method = c("KM01", "H00", "K04")) {
  
  # Simplify variables names
  x <- rep(NA, nrow(data))
  p <- sort(percent)
  method <- rlang::arg_match(method)
  
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
#'   indicating x and y coordinates. Template returned by grid_init.
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


#' Footprint-weighted land cover
#' 
#' @param x A numeric matrix
#' @param y A numeric matrix
#' @param type The type of convolution to take place, either "factor" or 
#'   "numeric"
#' @param levels A vector
#' @param ignore_levels A vector
#'
#' @importFrom stats "na.omit"
#' @export
summarize_cover <- function(x, y, type = c("factor", "numeric"), levels, 
                            ignore_levels = 0) {
  
  type <- rlang::arg_match(type)
  
  if (type == "numeric") {
    
    # Convolve matrices, normalize by total footprint weight
    #out <- sum((x / sum(x)) * y)
    #out <- sum(x * y) / sum(x)
    out <- sum(x * y)
    
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
