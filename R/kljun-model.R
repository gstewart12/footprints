#' Kljun et al. (2015) footprint model
#'
#' @name Kljun
#' @param x 
#' @param ustar 
#' @param mo_length 
#' @param blh 
#' @param z 
#' @param zd 
#' @param z0 
#' @param ws 
#'
#' @keywords internal
kljun_1d <- function(x, ustar, mo_length, blh, z, zd, z0, ws = NULL) {
  
  # Set model constants
  n_lim <- 5000 # Limit for neutral scaling
  k <- 0.4  # Von Karman constant
  zm <- z - zd  # aerodynamic measurement height
  L <- mo_length
  
  # Crosswind-integrated footprint parameters (Eq. 17)
  a <- 1.4524
  b <- -1.9914
  c <- 1.4622
  d <- 0.1359
  
  # Non-dimensional wind shear (Eq. 5)
  if (L <= 0 | L >= n_lim) {
    chi <- (1 - 19 * zm/L)^(1/4)
    psim <- log((1 + chi^2) / 2) + 2 * log((1 + chi) / 2) - 2 * atan(chi) + pi/2
  } else if (L > 0 & L < n_lim) {
    psim <- -5.3 * zm/L
  }
  
  # Dimensionless Pi groups (Eqs. 4)
  pi2 <- x/zm
  pi3 <- 1 - zm/blh
  pi4 <- if (!is.null(ws)) {
    ws/ustar * k
  } else {
    log(zm/z0) - psim
  }
  
  # Fail gracefully under invalid conditions
  if (pi4 <= 0) {
    return(list(xstar = rep(NA, length(x)), f = rep(NA, length(x))))
  }
  
  # Non-dimensional upwind distance X* for given x (Eqs. 7-8)
  xstar <- pi2 * pi3/pi4
  
  # Peak location of influence (Eqs. 20-22)
  xstar_max <- -c/b + d
  x_max <- xstar_max * zm/pi3 * pi4
  
  # Parameterized crosswind-integrated footprint F* (Eq. 14)
  fstar <- a * (xstar - d)^b * exp(-c / (xstar - d))
  fstar[which(xstar - d < 0)] <- 0
  
  # Crosswind-integrated footprint F (Eqs. 9-10)
  f <- fstar/zm * pi3/pi4
  
  # Return list including intermediates
  list(xstar = xstar, x_max = x_max, f = f)
}


#' Kljun et al. (2015) footprint model
#' 
#' @rdname Kljun
#' @keywords internal
kljun_2d <- function(x, y, ustar, mo_length, v_sigma, blh, z, zd, z0, 
                     ws = NULL) {
  
  # Set model constants
  n_lim <- 5000 # Limit for neutral scaling
  zm <- z - zd  # aerodynamic measurement height
  L <- mo_length
  # Deviation of crosswind distance parameters (Eq. 19)
  ac <- 2.17
  bc <- 1.66
  cc <- 20.0
  
  # Crosswind integrated footprint
  ci <- kljun_1d(x, ustar, L, blh, zm, zd, z0, ws)
  
  # Fail gracefully under invalid conditions
  if (all(is.na(ci$f))) {
    return(rep(NA, length(x)))
  }
  
  # Parameterized deviation of crosswind distance sigma_y* (Eq. 18)
  sigystar <- ac * (bc * ci$xstar^2 / (1 + cc * ci$xstar))^(1/2)
  
  # Deviation of crosswind distance sigma_y (Eq. 13)
  if (abs(L) > n_lim) L <- -1000000
  p <- if (L <= 0) 0.8 else 0.55
  
  ps1 <- 1e-5 * abs(zm/L)^(-1) + p
  if (ps1 > 1) ps1 <- 1
  
  sigy <- sigystar/ps1 * zm * v_sigma/ustar
  
  # 6. Crosswind dispersion (Eq. 10)
  dy <- 1 / (sqrt(2 * pi) * sigy) * exp(-y^2 / (2 * sigy^2))
  
  # Return 2D footprint adjusted by lateral diffusion
  ci$f * dy 
}