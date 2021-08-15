#' Kormann & Meixner (2001) footprint model
#' 
#' @name Kormann-Meixner
#' @param x 
#' @param ws 
#' @param ustar 
#' @param mo_length 
#' @param z 
#' @param zd 
#'
#' @keywords internal
kormann_1d <- function(x, ws, ustar, mo_length, z, zd) {
  
  # https://github.com/KIT-HYD/bridget/blob/3d0c8a7047b84714282d2f033c99a92d2372bbd4/bridget/eddy/footprint.py
  
  # Set model constants
  zm <- z - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  L <- mo_length
  
  # Dimensionless gradient functions of wind and temp profiles (Dyer 1974)
  # Non-dimensional wind shear (Eq. 33)
  phi_m <- if (L > 0) { 
    1 + 5 * zm/L
  } else if (L < 0) {
    (1 - 16 * zm/L)^(-1/4)
  }
  
  # Heat (Eq. 34)
  phi_c <- f_phi_c(zm, L)
  
  # Eddy diffusivity (Eq. 32)
  K <- (k * ustar * zm) / phi_c
  
  # Power law exponents (Eqs. 36)
  m <- (ustar/k) * (phi_m/ws) # wind velocity
  n <- if (L > 0) { # eddy diffusivity
    1 / (1 + 5 * zm/L)
  } else if (L < 0) {
    (1 - 24 * zm/L) / (1 - 16 * zm/L)
  } 
  
  # Constants in power-law profiles (Eqs. 11)
  U <- ws / zm^m # wind velocity
  kappa <- K / zm^n # eddy diffusivity
  
  # Intermediate parameters
  r <- 2 + m - n # shape factor
  mu <- (1 + m) / r # constant
  
  # Flux length scale (Eq. 19)
  xi <- (U * zm^r) / (r^2 * kappa)
  
  # Crosswind integrated flux (Eq. 21)
  f <- (1/gamma(mu)) * (xi^mu / x^(1 + mu)) * exp(-xi/x)
  
  # Return list including intermediates
  list(m = m, U = U, kappa = kappa, r = r, mu = mu, f = f)
}


#' @rdname Kormann-Meixner
#' @keywords internal
kormann_1d_num <- function(x, ws, ustar, mo_length, z, zd) {
  
  # Set model constants
  zm <- z - zd  # aerodynamic measurement height
  k <- 0.4  # Von Karman constant
  L <- mo_length
  
  # Roughness length from rearranging Eq. 31
  z0 <- z / exp(ws/ustar * k + f_psi_m(z, L))
  
  # Integration bounds for Eqs. 42-46
  z1 <- 3 * z0
  z2 <- (1 + k) * zm
  
  # Find power law exponent roots
  m_root <- stats::uniroot(
    f_m, c(0, 3), z0 = z0, zm = zm, L = L, z1 = z1, z2 = z2, extendInt = "upX", 
    maxiter = 1000
  )
  m <- m_root$root
  n_root <- stats::uniroot(
    f_n, c(0, 3), z0 = z0, zm = zm, L = L, z1 = z1, z2 = z2, extendInt = "upX",
    maxiter = 1000
  )
  n <- n_root$root
  
  # Constants in power-law profiles (Eqs. 41)
  # Wind velocity
  U <- ustar/k * (I2(m, z0/zm, z1, z2, zm) + J1(m, z1, z2, zm, L, f_psi_m)) / 
    (I1(2 * m, z1, z2, zm) * zm^m)
  # Eddy diffusivity
  kappa <- k * ustar * J1(n, z1, z2, zm, L, f_z_phi_c) / 
    (I1(2 * n, z1, z2, zm) * zm^(n - 1))
  
  # Intermediate parameters
  r <- 2 + m - n # shape factor
  mu <- (1 + m) / r # constant
  
  # Flux length scale (Eq. 19)
  xi <- U * zm^r / (r^2 * kappa)
  
  # Crosswind integrated flux (Eq. 21)
  f <- 1/gamma(mu) * xi^mu/x^(1 + mu) * exp(-xi/x)
  
  # Return list including intermediates
  list(m = m, U = U, kappa = kappa, r = r, mu = mu, f = f)
}


#' @rdname Kormann-Meixner
#' @keywords internal
kormann_2d <- function(x, y, ws, ustar, mo_length, v_sigma, z, zd, 
                       approach = c("analytical", "numerical")) {
  
  approach <- match.arg(approach)
  ci_fun <- if (approach == "numerical") kormann_1d_num else kormann_1d
  L <- mo_length
  
  # Crosswind integrated flux
  #ci <- kormann_1d(x, ws, ustar, mo_length, z, zd)
  ci <- rlang::exec(
    ci_fun, x = x, ws = ws, ustar = ustar, mo_length = L, z = z, zd = zd
  )
  m <- ci$m
  U <- ci$U
  kappa <- ci$kappa
  r <- ci$r
  mu <- ci$mu
  f <- ci$f
  
  # Effective plume velocity (eq. 18)
  ubar <- gamma(mu) / gamma(1/r) * (r^2 * kappa/U)^(m/r) * U * x^(m/r)
  
  # Crosswind distribution function (eq. 9)
  sigmay <- v_sigma * x/ubar
  # Dy <- 1/(sqrt(2 * pi) * sigmay) * exp(-y^2 / (2 * sigmay^2))
  Dy <- 1/(sqrt(2 * pi) * sigmay) * exp(-(y^2 / (2 * sigmay^2)))
  
  # Crosswind distributed flux (eq. 1)
  f * Dy
}

#' Helper functions for Kormann & Meixner (2001) model
#' @name KM-helpers
#' @keywords internal
NULL
#> NULL

#' Heat (Eq. 34)
#' @rdname KM-helpers
#' @keywords internal
f_phi_c <- function(zm, L, ...) {
  if (L > 0) { 
    phi_c <- 1 + 5 * zm/L
  } else if (L < 0) {
    phi_c <- (1 - 16 * zm/L)^(-1/2)
  }
  phi_c
}

#' Diabatic integration of the wind profile (Eq. 35)
#' @rdname KM-helpers
#' @keywords internal
f_psi_m <- function(zm, L, ...) {
  zeta <- (1 - 16 * zm/L)^(1/4)
  if (L > 0) { 
    psi_m <- 5 * zm/L
  } else if (L < 0) {
    psi_m <- -2 * log((1 + zeta) / 2) - log((1 + zeta^2) / 2) + 
      2 * atan(zeta) - pi/2
  }
  psi_m
}

#' Helper for Eqs. 40 & 41
#' @rdname KM-helpers
#' @keywords internal
f_z_phi_c <- function(zeta, zm, L) {
  zeta / f_phi_c(zeta, L) * zm
}

#' Power law exponent - wind velocity (Eq. 39)
#' @rdname KM-helpers
#' @keywords internal
f_m <- function(m, z0, zm, L, z1, z2) {
  a <- I1(2 * m, z1, z2, zm) * 
    (I3(m, z0/zm, z1, z2, zm) + J2(m, z1, z2, zm, L, f_psi_m))
  b <- I2(2 * m, 1, z1, z2, zm) * 
    (I2(m, z0/zm, z1, z2, zm) + J1(m, z1, z2, zm, L, f_psi_m))
  b - a
}

#' Power law exponent - eddy diffusivity (Eq. 40)
#' @rdname KM-helpers
#' @keywords internal
f_n <- function(n, z0, zm, L, z1, z2) {
  a <- I1(2 * n, z1, z2, zm) * J2(n, z1, z2, zm, L, f_z_phi_c)
  b <- I2(2 * n, 1, z1, z2, zm) * J1(n, z1, z2, zm, L, f_z_phi_c)
  b - a
}

#' Integration functions (Eqs. 42-46)
#' @rdname KM-helpers
#' @keywords internal
I1 <- function(p, z1, z2, zm) {
  f <- function(zeta, p) zeta^p
  i <- stats::integrate(
    f = f, lower = z1/zm, upper = z2/zm, p = p, subdivisions = 1000
  )
  i$value
}

#' Integration functions (Eqs. 42-46)
#' @rdname KM-helpers
#' @keywords internal
I2 <- function(p, zeta0, z1, z2, zm) {
  f <- function(zeta, p, zeta0) zeta^p * log(zeta/zeta0)
  i <- stats::integrate(
    f = f, lower = z1/zm, upper = z2/zm, p = p, zeta0 = zeta0,
    subdivisions = 1000
  )
  i$value
}

#' Integration functions (Eqs. 42-46)
#' @rdname KM-helpers
#' @keywords internal
I3 <- function(p, zeta0, z1, z2, zm) {
  f <- function(zeta, p, zeta0) zeta^p * log(zeta) * log(zeta/zeta0)
  i <- stats::integrate(
    f = f, lower = z1/zm, upper = z2/zm, p = p, zeta0 = zeta0,
    subdivisions = 1000
  )
  i$value
}

#' Integration functions (Eqs. 42-46)
#' @rdname KM-helpers
#' @keywords internal
J1 <- function(p, z1, z2, zm, L, .f) {
  f <- function(zeta, p, zm, L, .f) zeta^p * rlang::exec(.f, zeta * zm, zm, L)
  i <- stats::integrate(
    f = f, lower = z1/zm, upper = z2/zm, p = p, zm = zm, L = L, .f = .f,
    subdivisions = 1000
  )
  i$value
}

#' Integration functions (Eqs. 42-46)
#' @rdname KM-helpers
#' @keywords internal
J2 <- function(p, z1, z2, zm, L, .f) {
  f <- function(zeta, p, zm, L, .f) {
    zeta^p * rlang::exec(.f, zeta * zm, zm, L) * log(zeta)
  }
  i <- stats::integrate(
    f = f, lower = z1/zm, upper = z2/zm, p = p, zm = zm, L = L, .f = .f,
    subdivisions = 1000
  )
  i$value
}