
#' Chebyshev Interpolation Function
#' 
#' This function takes in a minimum and a maximum range, a function, and a
#' tolerance, and returns a chebyshev approximation to a function
#' @keywords chebyshev approximation
#' @param fxn The function we wish to approximate
#' @param lower The lower bound of the approximation of the range
#' @param upper The upper bound of the approximation of the range
#' @param tol Tolerance of the functional approximation
#' @export
#' @examples
#' f <- function(x) {
#'    x^5
#' }
#' chebyshev_compute(f, lower = -2, upper = 2, tol = 1E-10)

chebyshev_compute <- function(fxn, lower, upper, tol = 1E-15, n_subdivisions = 100, maxN = 3000, maxEfficiency = F) {
  # Start with three points - a quadratic approximation
  N = 2
  
  while (TRUE) {
    cat(sprintf("Iteration: N = %*.0f + 1 | ", 4, N))
    weights = c(1/2, rep(1, N-1), 1/2)
    k = 0:N
    f_points = -cos(k * pi / N)
    f_points_scaled = (f_points + 1) / 2 * (upper - lower) + lower
    # Calculate the maximum difference
    diff = laply(.data = seq(from = lower + 1E-14, to = upper - 1E-14, length = n_subdivisions), .fun = function(x){
      # Move the point from [lower,upper]->[-1,1]
      x <- 2 * (x - lower) / (upper - lower) - 1
      # Ensure that y_j points in [lower,upper] are usable
      p_x = sum(weights * (-1)^k * fxn(f_points_scaled) / (x - f_points)) / sum(weights * (-1)^k / (x - f_points))
       # Return the absolute deviation
      browser()
      abs(p_x - fxn(x))
    })
    # Break if tolerance is good, else evaluate updating conditions
    cat(sprintf("Max Diff in Digits: %f\n", abs(log10(max(diff, na.rm = T)))))
    if (max(diff, na.rm = T) < tol) {
      break
    } else if (N > maxN) {
      stop("The function needs a lot of precision, function evaluation will stop past 3000 maxN")
    } else if (maxEfficiency) {
      # in maxEfficiency, we attempt to find the optimum number of points
      N = N + 1 
    } else {
      # Do an expotential growth of n
      N = floor(N * 1.5)
    }
  }
  cat("Choice of n:", N, "\n")
  # Now we return a function that interpolates between [lower,upper] for our function
  function(val) {
    # Calculate the weights and the x_j
    weights = c(1/2, rep(1, N-1), 1/2)
    k = 0:N
    f_points = -cos(k * pi / N)
    
    # Ensure we are not extrapolating
    stopifnot(val >= lower & val <= upper)
    
    # For each point to eval, compute the polynomial approximation
    result = sapply(val, function(x) {
      # Scale x to the [-1,1] interval
      x <- 2 * (x - lower) / (upper - lower) - 1
      # Ensure that y_j points in [lower,upper] are usable
      f_points_scaled = (f_points + 1) / 2 * (upper - lower) + lower
      p_x = sum(weights * (-1)^k * fxn(f_points_scaled) / (x - f_points)) / sum(weights * (-1)^k / (x - f_points))
    })
    return(result)
  }
  
  
}

f <- function(x) {
  cos(x)
}

l = -2
u = 2
# laply(seq(1/2,8), function(x) {
  myfunc = chebyshev_compute(f, lower = l, upper = u, tol = 10^(-8), maxEfficiency = T, n_subdivisions = 1000)
  curve(f, l,u, lty = 2, col = "red")
  curve(myfunc, l,u, add = T, col = "blue", lty = 4, pch = 2)
  p = seq(l,u, length = 1000)
  plot(log10(abs(f(p) - myfunc(p))), type = "p", ylim = c(-15,0))
# })


