#' Hsieh et al. (2000) footprint model
#'
#' @param x 
#' @param mo_length 
#' @param zm 
#' @param zd 
#' @param zo 
#'
#' @return
#' @export
#'
#' @examples
hsieh_1d <- function(x, mo_length, zm, zd, zo) {
  
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