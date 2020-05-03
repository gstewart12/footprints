
## =============================================================================
kljun_model <- function(x, y, ustar, mo_length, v_sigma, blh, z, zd, zo, dims) {
  
  # KLJUN MODEL - via Kljun et al. (2015)
  
  # Input variables
  ol <- mo_length
  h <- blh
  
  # Construct helper matrices
  xstar <- matrix(0, nrow = dims[1], ncol = dims[2])
  fstar <- matrix(0, nrow = dims[1], ncol = dims[2])
  sigystar <- matrix(0, nrow = dims[1], ncol = dims[2])
  f <- matrix(0, nrow = dims[1], ncol = dims[2])
  sigy <- matrix(0, nrow = dims[1], ncol = dims[2])
  dy <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Set model constants
  
  # Crosswind-integrated footprint parameters (Eq. 17)
  a <- 1.4524
  b <- -1.9914
  c <- 1.4622
  d <- 0.1359
  
  # Deviation of crosswind distance parameters (Eq. 19)
  ac <- 2.17
  bc <- 1.66
  cc <- 20.0
  
  ol_n <- 5000 # Limit for neutral scaling
  k <- 0.4  # Von Karman constant
  
  zm <- z - zd  # aerodynamic measurement height
  
  # Non-dimensional wind shear Psi (Eq. 5)
  if (ol <= 0 | ol >= ol_n) {
    xx <- (1 - 19.0 * zm / ol)^0.25
    psi <- log((1 + xx^2) / 2) + 2 * log((1 + xx) / 2) - 2 * atan(xx) + pi / 2
  } else if (ol > 0 & ol < ol_n) {
    psi <- -5.3 * zm / ol
  }
  
  # Subset downwind area of grid for calculating contributions
  dw <- x > 0
  
  # Fail gracefully under invalid conditions
  if ((log(zm / zo) - psi) <= 0) {
    f[dw] <- NA
    return(f)
  }
  
  # 1. Non-dimensional upwind distance X* for given x (Eq. 7)
  xstar[dw] <- (x[dw] / zm) * (1 - zm / h) / (log(zm / zo) - psi)
  
  # 2. Parameterized crosswind-integrated footprint F* (Eq. 14)
  fstar[dw] <- a * (xstar[dw] - d)^b * exp(-c / (xstar[dw] - d))
  fstar[which(xstar - d < 0)] <- 0
  
  # 3. Parameterized deviation of crosswind distance sigma_y* (Eq. 18)
  sigystar[dw] <- ac * (bc * xstar[dw]^2 / (1 + cc * xstar[dw]))^0.5
  
  # 4. Crosswind-integrated footprint F (Eq. 9)
  f[dw] <- fstar[dw] / zm * (1 - (zm / h)) / (log(zm / zo) - psi)
  
  # 5. Deviation of crosswind distance sigma_y (Eq. 13)
  if (abs(ol) > ol_n) ol <- -1000000
  if (ol <= 0 ) p <- 0.8 else if (ol > 0) p <- 0.55
  ps1 <- 1e-5 * abs(zm / ol)^(-1) + p
  if (ps1 > 1) ps1 <- 1
  
  sigy[dw] <- sigystar[dw] / ps1 * zm * v_sigma / ustar
  
  # 6. Crosswind dispersion (Eq. 10)
  dy[dw] <- 1 / (sqrt(2 * pi) * sigy[dw]) * exp(-y[dw]^2 / (2 * sigy[dw]^2))
  
  # Return 2D footprint adjusted by lateral diffusion
  dy[dw] * f[dw]
}


## =============================================================================
hsieh_model <- function(x, y, ustar, mo_length, v_sigma, z, zd, zo, dims) {
  
  ## HSIEH MODEL - via Hsieh et al. (2000)
  
  # Construct helper matrices
  f <- matrix(0, nrow = dims, ncol = dims)
  sig <- matrix(0, nrow = dims, ncol = dims)
  dy <- matrix(0, nrow = dims, ncol = dims)
  
  # Set model constants
  zm <- z - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  
  # Define and calculate footprint parameters
  zu <- zm * (log(zm / zo) - 1 + (zo / zm))
  zl <- zu / mo_length
  if (abs(zl) < 0.04) {
    d <- 0.97
    p <- 1
  } else if (zl < 0) {
    d <- 0.28
    p <- 0.59
  } else if (zl > 0) {
    d <- 2.44
    p <- 1.33
  }
  dp <- d * zu^p * abs(mo_length)^(1 - p)
  
  # Subset downwind area of grid for calculating contributions
  dw <- x > 0
  
  # Calculate cross-wind integrated footprint
  f[dw] <- (dp / (k^2 * x[dw]^2)) * (exp(-dp / (k^2 * x[dw])))
  
  # Define and calculate 2D extension parameters
  a1 <- 0.3
  p1 <- 0.86
  sig[dw] <- a1 * zo * (v_sigma / ustar) * (x[dw] / zo) ^ p1
  dy[dw] <- 1 / (sqrt(2 * pi) * sig[dw]) * exp(-0.5 * (y[dw] / sig[dw]) ^ 2)
  
  # Return 2D footprint adjusted by lateral diffusion
  dy[dw] * f[dw]
  
}


## =============================================================================
kormann_model <- function(x, y, ws, ustar, mo_length, v_sigma, zm, zd, dims) {
  
  ## KORMANN MODEL - via Kormann & Meixner (2001)
  
  # Construct helper matrices
  ubar <- matrix(0, nrow = dims[1], ncol = dims[2])
  f <- matrix(0, nrow = dims[1], ncol = dims[2])
  sigmay <- matrix(0, nrow = dims[1], ncol = dims[2])
  Dy <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Subset downwind area of grid for calculating contributions
  dw <- x > 0
  
  # Set model constants
  z <- zm - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  zl <- z / mo_length
  
  # Dimensionless gradient functions of wind and temp profiles (Dyer 1974)
  # Non-dimensional wind shear (Eq. 33)
  phim <- if (mo_length > 0) { 
    1 + 5 * zl
  } else if (mo_length < 0) {
    (1 - 16 * zl)^(-1 / 4)
  }
  
  # Heat (Eq. 34)
  phic <- if (mo_length > 0) { 
    1 + 5 * zl
  } else if (mo_length < 0) {
    (1 - 16 * zl)^(-1 / 2)
  }
  
  # Eddy diffusivity (Eq. 32)
  K <- (k * ustar * z) / (phic * zl)
  
  # Power law exponents (Eqs. 36)
  m <- (ustar / k) * (phim / ws) # wind velocity
  n <- if (mo_length > 0) { # eddy diffusivity
    1 / (1 + 5 * zl)
  } else if (mo_length < 0) {
    (1 - 24 * zl) / (1 - 16 * zl)
  } 
  
  # Constants in power-law profiles (Eqs. 11)
  U <- ws / z^m # wind velocity
  kappa <- K / z^n # eddy diffusivity
  
  # Intermediate parameters
  r <- 2 + m - n # shape factor
  mu <- (1 + m) / r # constant
  
  # Flux length scale (Eq. 19)
  xi <- (U * z^r) / (r^2 * kappa)
  
  # Effective plume velocity (Eq. 18)
  ubar[dw] <- (
    (gamma(mu) / gamma(1 / r)) * (r^2 * kappa / U)^(m / r) * U * x[dw]^(m / r)
  )
  
  # Crosswind integrated flux (Eq. 21)
  f[dw] <- (1 / gamma(mu)) * (xi^mu / x[dw]^(1 + mu)) * exp(-xi / x[dw])
  
  # Crosswind distribution function (Eq. 9)
  sigmay[dw] <- v_sigma * x[dw] / ubar[dw]
  Dy[dw] <- (1 / (sqrt(2 * pi) * sigmay[dw])) * exp(-y[dw]^2/(2 * sigmay[dw]^2))
  
  # Crosswind distributed flux (Eq. 1)
  f[dw] * Dy[dw]
}

## =============================================================================
#' Fetch Distance Using Hsieh (2000) Model
#'
#' @param z
#' @param zd
#' @param zo
#' @param mo_length
#' @param p
#'
#' @export
hsieh_ftp_dist <- function(z, zd, zo, mo_length, p) {
  
  # Recommended to gap-fill L prior to using this model
  if (!is.na(mo_length)) {
    # Define and calculate footprint parameters
    k <- 0.4
    dx <- 5
    zm <- z - zd
    zu <- zm * (log(zm / zo) - 1 + (zo / zm))
    zL <- zu / mo_length
    if (abs(zL) < 0.04) {
      # Near neutral and neutral conditions
      D <- 0.97
      P <- 1
      mu <- 500
    } else if (zL < 0) {
      # Unstable conditions
      D <- 0.28
      P <- 0.59
      mu <- 100
    } else if (zL > 0) {
      # Stable conditions
      D <- 2.44
      P <- 1.33
      mu <- 2000
    }
    DP <- D * zu^P * abs(mo_length)^(1 - P)
    # Calculate x if flux percents are given
    all_p <- p
    p <- as.numeric(all_p[!suppressWarnings(is.na(as.numeric(all_p)))])
    if (length(p) > 0) {
      # Define footprint model
      #f_int <- function(x) exp(-(D * zu^P * abs(L)^(1 - P) / (k^2 * x)))
      #f_int <- function(x) {
      #  fact <- D * zu^P * abs(L)^(1 - P) / (k^2 * (x * dx))
      #  exp(-fact)
      #}
      f_int <- function(x) {
        c1 <- (-1 / k^2) * (D * zu^P * abs(mo_length)^(1 - P)) / x
        exp(c1)
        #fc <- exp(c1) # Cumulative source contribution with distance
        #-(c1 / x) * fc # Source-weight function
      }
      
      # Peak distance
      #xp <- (1 / 2 / k^2) * (D * zu^P * abs(L)^(1 - P))
      # Fetch to height ratio
      #f2h <- (D / 0.105 / k^2) * (zm^(-1) * abs(L)^(1 - P) * zu^(P))
      
      min_x <- (mu / 100) * zo
      max_x <- mu * zm
      x_bin <- (max_x - min_x) / 1000
      
      p <- c(1, sort(p))
      fx <- rep(NA, length(p))
      #x <- 1
      x <- min_x
      #browser()
      while (x < max_x) {
        fp <- f_int(x)
        fx[which(((fp > (p / 100) & is.na(fx))) == TRUE)[1]] <- x
        #fx[which(((fp > (p / 100) & is.na(fx))) == TRUE)[1]] <- x * dx
        if (!anyNA(fx)) break # finished if all fetch lengths are found
        x <- x + x_bin
      }
    }
    # Calculate x if flux peak is requested
    if ("peak" %in% all_p) {
      peak <- DP / (2 * k^2)
      fx <- c(fx, peak)
    }
  } else {
    # Cannot get fetch lengths if meteorological parameters are absent
    fx <- rep(NA, length(p) + 1)
  }
  
  fx
}

## =============================================================================
#' Fetch Distance Using Kormann & Meixner (2001) Model
#'
#' @param ws
#' @param ustar
#' @param zeta
#' @param z
#' @param zd
#' @param p
#'
#' @export
kormann_ftp_dist <- function(ws, ustar, zeta, z, zd, p) {
  
  # Simplify variable names
  if (!is.na(ws) & !is.na(ustar) & !is.na(zeta)) {
    # Calculate model parameters
    zm <- z - zd  # measurement height
    k <- 0.4  # Von Karman constant
    if (zeta > 0) {
      # similarity relations
      phi_m <- 1 + 5 * zeta
      phi_c <- phi_m
      # exponent of diffusivity power law
      n <- 1 / phi_m
    } else if (zeta <= 0) {
      # similarity relations
      phi_m <- (1 - 16 * zeta)^-0.25
      phi_c <- (1 - 16 * zeta)^-0.5
      eta <- (1 - 16 * zeta)^0.25
      # exponent of the diffusivity power law
      n <- (1 - 24 * zeta) / (1 - 16 * zeta)
    }
    # Proportionality constant of the diffusivity power law (Eqs. 11 & 32)
    key <- k * ustar * zm / (phi_c * zm^n)
    # Exponent of the wind speed power law
    m <- ustar * phi_m / (k * ws)
    # Proportionality constant of the wind speed power law (Eqs. 11 & 31)
    U <- ws / (zm^m)
    # Intermediate parameters
    r <- 2 + m - n  # shape factor
    mu <- (1 + m) / r  # constant
    xi <- U * zm^r / (r^2 * key)  # flux length scale
    # Calculate x if flux percents are given
    all_p <- p
    p <- as.numeric(all_p[!suppressWarnings(is.na(as.numeric(all_p)))])
    if (length(p) > 0) {
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
  } else {
    # Cannot get fetch lengths if meteorological parameters are absent
    fx <- rep(NA, length(p))
  }
  fx
}