
## =============================================================================
kljun_model <- function(x, y, ustar, mo_length, v_sigma, blh, z, zd, zo, dims) {
  
  # KLJUN MODEL - via Kljun et al. (2015)
  
  # Input variables
  ol <- mo_length
  sigmav <- v_sigma
  ustar <- ustar
  h <- blh
  
  # Construct helper matrices
  xstar <- matrix(0, nrow = dims, ncol = dims)
  fstar <- matrix(0, nrow = dims, ncol = dims)
  sigystar <- matrix(0, nrow = dims, ncol = dims)
  f <- matrix(0, nrow = dims, ncol = dims)
  sigy <- matrix(0, nrow = dims, ncol = dims)
  dy <- matrix(0, nrow = dims, ncol = dims)
  
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
  
  sigy[dw] <- sigystar[dw] / ps1 * zm * sigmav / ustar
  
  # 6. Crosswind dispersion (Eq. 10)
  dy[dw] <- 1 / (sqrt(2 * pi) * sigy[dw]) * exp(-y[dw]^2 / (2 * sigy[dw]^2))
  
  # Return 2D footprint adjusted by lateral diffusion
  dy[dw] * f[dw]
}


## =============================================================================
hsieh_model <- function(x, y, ustar, mo_length, v_sigma, z, zd, zo, dims) {
  
  ## HSIEH MODEL - via Hsieh et al. (2000)
  
  # Input variables
  ol <- mo_length
  sigmav <- v_sigma
  ustar <- ustar
  
  # Construct helper matrices
  f <- matrix(0, nrow = dims, ncol = dims)
  sig <- matrix(0, nrow = dims, ncol = dims)
  dy <- matrix(0, nrow = dims, ncol = dims)
  
  # Set model constants
  zm <- z - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  
  # Define and calculate footprint parameters
  zu <- zm * (log(zm / zo) - 1 + (zo / zm))
  zl <- zu / ol
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
  dp <- d * zu^p * abs(ol)^(1 - p)
  
  # Subset downwind area of grid for calculating contributions
  dw <- x > 0
  
  # Calculate cross-wind integrated footprint
  f[dw] <- (dp / (k^2 * x[dw]^2)) * (exp(-dp / (k^2 * x[dw])))
  
  # Define and calculate 2D extension parameters
  a1 <- 0.3
  p1 <- 0.86
  sig[dw] <- a1 * zo * (sigmav / ustar) * (x[dw] / zo) ^ p1
  dy[dw] <- 1 / (sqrt(2 * pi) * sig[dw]) * exp(-0.5 * (y[dw] / sig[dw]) ^ 2)
  
  # Return 2D footprint adjusted by lateral diffusion
  dy[dw] * f[dw]
  
}


## =============================================================================
kormann_model <- function(x, y, ws, ustar, mo_length, v_sigma, z, zd, zo,
                          dims) {
  
  ## KORMANN MODEL - via Forbrich et al. (2011)
  
  # Input variables
  ol <- mo_length
  sigmav <- v_sigma
  ustar <- ustar
  
  # Construct helper matrices
  h0 <- matrix(0, nrow = dims, ncol = dims)
  h1 <- matrix(0, nrow = dims, ncol = dims)
  h2 <- matrix(0, nrow = dims, ncol = dims)
  h3 <- matrix(0, nrow = dims, ncol = dims)
  
  # Set model constants
  zm <- z - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  zmol <- zm / ol
  
  # Dimensionless gradient functions of wind and temp profiles (Dyer 1974)
  if (l > 0) {
    phi <- 1 + 5 * zmol
    phim <- 1 + 5 * zmol
    phic <- 1 + 5 * zmol
    m <- (1 + 5 * (zmol)) / (log(zm / zo) + 5 * zmol)
    n <- 1 / (1 + 5 * zmol)
  } else if (l < 0) {
    zeta <- (1 - 16 * (zmol))^(-1 / 4)
    phi <- 1 / (zeta^2)
    phim <- (1 - 16 * (zmol))^(-1 / 4)
    phic <- (1 - 16 * (zmol))^(-1 / 2)
    m <- ustar / k / zeta / ws
    n <- (1 - 24 * zmol) / (1 - 16 * zmol)
  }
  
  # Coefficients of diffusivity power law (Eqs. A.4-A.5)
  #kappa <- (k * ustar * zm) / (phic * zm^n)
  kappa <- (k * ustar * zm) / (phi * zm^n)
  # Many different ways to calculate U
  #U <- (ustar * log((zm / zo) + (phim))) / (k * zm^m)
  #U <- (ustar * phim) / (k * m * zm^m) # neftel
  U <- ws / zm^n
  
  # Intermediate parameters (Eqs. A.6-A.7; p = m)
  r <- 2 + m - n
  mu <- (1 + m) / r
  
  # Flux length scale (Eq. A.8)
  xi <- (U * zm^r) / (r^2 * kappa)
  
  # Neftel parameters
  A <- 1 + mu
  B <- (U * (zm^r)) / (r^2 * kappa)
  C <- B^mu / gamma(mu)
  D <- sigmav / U * gamma(1 / r) / gamma(mu) * (U / (r^2 * kappa))^(m / r)
  E <- (r - m) / r
  
  # Subset downwind area of grid for calculating contributions
  dw <- x > 0
  
  # Splitting Eq 2.1 into a few helpers as per OzFlux
  h0[dw] <- x[dw]^E
  h1[dw] <- (sqrt(2 * pi) * D * (h0[dw]))
  h2[dw] <- exp(-(y[dw]^2) / (2 * (D * (h0[dw]))^2))
  h3[dw] <- C * exp(-B / x[dw]) * (x[dw]^(-A))
  
  # Combine helpers - two-dimensional footprint function (Eq. 2.1)
  h2[dw] * h3[dw] / h1[dw]
  
}