
## =============================================================================
kljun_1d <- function(x, ustar, mo_length, blh, zm, zd, zo) {
  
  # Set model constants
  n_lim <- 5000 # Limit for neutral scaling
  k <- 0.4  # Von Karman constant
  z <- zm - zd  # aerodynamic measurement height
  
  # Crosswind-integrated footprint parameters (Eq. 17)
  a <- 1.4524
  b <- -1.9914
  c <- 1.4622
  d <- 0.1359
  
  # Non-dimensional wind shear Psi (Eq. 5)
  if (mo_length <= 0 | mo_length >= n_lim) {
    xx <- (1 - 19.0 * z / mo_length)^0.25
    psi <- log((1 + xx^2) / 2) + 2 * log((1 + xx) / 2) - 2 * atan(xx) + pi / 2
  } else if (mo_length > 0 & mo_length < n_lim) {
    psi <- -5.3 * z / mo_length
  }
  
  # Fail gracefully under invalid conditions
  if ((log(z / zo) - psi) <= 0) {
    return(list(z = z, xstar = rep(NA, length(x)), f = rep(NA, length(x))))
  }
  
  # Non-dimensional upwind distance X* for given x (Eq. 7)
  xstar <- (x / z) * (1 - z / blh) / (log(z / zo) - psi)
  
  # Parameterized crosswind-integrated footprint F* (Eq. 14)
  fstar <- a * (xstar - d)^b * exp(-c / (xstar - d))
  fstar[which(xstar - d < 0)] <- 0
  
  # Crosswind-integrated footprint F (Eq. 9)
  f <- fstar / z * (1 - (z / blh)) / (log(z / zo) - psi)
  
  # Return list including intermediates
  list(z = z, xstar = xstar, f = f)
}


## =============================================================================
kljun_2d <- function(x, y, ustar, mo_length, v_sigma, blh, zm, zd, zo) {
  
  # Set model constants
  n_lim <- 5000 # Limit for neutral scaling
  # Deviation of crosswind distance parameters (Eq. 19)
  ac <- 2.17
  bc <- 1.66
  cc <- 20.0
  
  # Crosswind integrated footprint
  ci <- kljun_1d(x, ustar, mo_length, blh, zm, zd, zo)
  
  # Fail gracefully under invalid conditions
  if (all(is.na(ci$f))) {
    return(rep(NA, length(x)))
  }
  
  # Parameterized deviation of crosswind distance sigma_y* (Eq. 18)
  sigystar <- ac * (bc * ci$xstar^2 / (1 + cc * ci$xstar))^0.5
  
  # Deviation of crosswind distance sigma_y (Eq. 13)
  if (abs(mo_length) > n_lim) mo_length <- -1000000
  p <- if (mo_length <= 0) 0.8 else 0.55
  
  ps1 <- 1e-5 * abs(ci$z / mo_length)^(-1) + p
  if (ps1 > 1) ps1 <- 1
  
  sigy <- sigystar / ps1 * ci$z * v_sigma / ustar
  
  # 6. Crosswind dispersion (Eq. 10)
  dy <- 1 / (sqrt(2 * pi) * sigy) * exp(-y^2 / (2 * sigy^2))
  
  # Return 2D footprint adjusted by lateral diffusion
  ci$f * dy 
}


## =============================================================================
kormann_1d <- function(x, ws, ustar, mo_length, zm, zd) {
  
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
  K <- (k * ustar * z) / phic
  
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
  
  # Crosswind integrated flux (Eq. 21)
  f <- (1 / gamma(mu)) * (xi^mu / x^(1 + mu)) * exp(-xi / x)
  
  # Return list including intermediates
  list(m = m, U = U, kappa = kappa, r = r, mu = mu, f = f)
}


## =============================================================================
kormann_2d <- function(x, y, ws, ustar, mo_length, v_sigma, zm, zd) {
  
  # Crosswind integrated flux
  ci <- kormann_1d(x, ws, ustar, mo_length, zm, zd)
  
  # Effective plume velocity (Eq. 18)
  ubar <- (gamma(ci$mu) / gamma(1 / ci$r)) * (ci$r^2 * ci$kappa / ci$U)^(ci$m / ci$r) * ci$U * x^(ci$m / ci$r)
  
  # Crosswind distribution function (Eq. 9)
  sigmay <- v_sigma * x / ubar
  Dy <- (1 / (sqrt(2 * pi) * sigmay)) * exp(-y^2 / (2 * sigmay^2))
  
  # Crosswind distributed flux (Eq. 1)
  ci$f * Dy
}


hsieh_1d <- function(x, mo_length, zm, zd, zo) {
  
  ## HSIEH MODEL - via Hsieh et al. (2000)
  
  # Set model constants
  z <- zm - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  
  # Length scale
  zu <- z * (log(z / zo) - 1 + (zo / z))
  zl <- zu / mo_length
  
  # Similarity constants
  if (abs(zl) < 0.04) {
    D <- 0.97
    P <- 1
    mu <- 500
  } else if (zl < 0) {
    D <- 0.28
    P <- 0.59
    mu <- 100
  } else if (zl > 0) {
    D <- 2.44
    P <- 1.33
    mu <- 2000
  }
  
  dp <- D * zu^P * abs(mo_length)^(1 - P)
  
  # Surface flux (rearranged Eq. 15)
  So <- (1 / (k^2 * x^2)) * D * zu^P * abs(mo_length)^(1 - P)
  
  # Cross-wind integrated flux (Eq. 16)
  Fy <- exp((-1 / (k^2 * x)) * D * zu^P * abs(mo_length)^(1 - P))
  
  # Cross-wind integrated footprint (Eq. 17)
  f <- So * Fy
  
  # Return list including intermediates
  list(mu = mu, Fy = Fy, f = f)
}


hsieh_2d <- function(x, y, ustar, mo_length, v_sigma, zm, zd, zo) {
  
  # Crosswind integrated flux
  ci <- hsieh_1d(x, mo_length, zm, zd, zo)
  
  ## HSIEH 2D EXTENSION - via Detto et al. (2006)
  
  # Define constants (similarity parameters)
  a1 <- 0.3
  p1 <- 0.86
  
  # Standard deviation (Eq. B4)
  sigmay <- a1 * zo * (v_sigma / ustar) * (x / zo) ^ p1
  
  # Lateral diffusion (Eq. B3)
  Dy <- 1 / (sqrt(2 * pi) * sigmay) * exp(-0.5 * (y / sigmay) ^ 2)
  
  # Return 2D footprint adjusted by lateral diffusion
  ci$f * Dy
}


## =============================================================================
hsieh_model <- function(x, y, ustar, mo_length, v_sigma, zm, zd, zo, dims) {
  
  ## HSIEH MODEL - via Hsieh et al. (2000)
  
  # Construct helper matrices
  f <- matrix(0, nrow = dims[1], ncol = dims[2])
  sig <- matrix(0, nrow = dims[1], ncol = dims[2])
  dy <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Set model constants
  z <- zm - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  
  # Define and calculate footprint parameters
  zu <- z * (log(z / zo) - 1 + (zo / z))
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


kormann_ftp_dist <- function(ws, ustar, mo_length, zm, zd, p) {
  
  # Simplify variable names
  if (!is.na(ws) & !is.na(ustar) & !is.na(mo_length)) {
    # Calculate model parameters
    z <- zm - zd  # measurement height
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
    K <- (k * ustar * z) / phic
    
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


lateral_diffusion <- function(y, sigma_y) {
  
  # Common lateral diffusion (Dy) estimation as in Schmid (1994)
  # - also referred to as crosswind dispersion/distribution
  
  1 / (sqrt(2 * pi) * sigma_y) * exp(-y^2 / (2 * sigma_y^2))
}

