
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

chebyshev_compute <- function(fxn, lower, upper, tol = 1E-10, n_subdivisions = 100, maxN = 3000, maxEfficiency = F) {
  # Start with three points - a quadratic approximation
  N = 2
  
  while (TRUE) {
    cat(sprintf("Iteration: N = %*.0f + 1 | ", 4, N))
    weights = c(1/2, rep(1, N-1), 1/2)
    k = 0:N
    f_points = -cos(k * pi / N)
    f_points_scaled = (f_points + 1) / 2 * (upper - lower) + lower
    # Calculate the maximum difference
    diff = sapply(seq(from = lower, to = upper, length = n_subdivisions), function(x) {
      # Move the point from [lower,upper]->[-1,1]
      x <- 2 * (x - lower) / (upper - lower) - 1
      # Ensure that y_j points in [lower,upper] are usable
      p_x = sum(weights * (-1)^k * fxn(f_points_scaled) / (x - f_points)) / sum(weights * (-1)^k / (x - f_points))
       # Return the absolute deviation
      abs(p_x - fxn((x + 1) / 2 * (upper - lower) + lower))
    })
    # Break if tolerance is good, else evaluate updating conditions
    cat(sprintf("Max Diff Digits: %*.5f | ", 8, abs(log10(max(diff, na.rm = T)))))
    cat(sprintf("ActualTol: %.5E\n", max(diff,na.rm=T)))
    if (max(diff, na.rm = T) < tol) {
      break
    } else if (N > maxN) {
      cat("The function needs a lot of precision, function evaluation stops past N > maxN")
      break
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
  # N = 50
  function(val) {
    # Ensure we are not extrapolating
    stopifnot(all(val >= lower) & all(val <= upper))
    
    # Calculate the weights and the x_j
    weights = c(1/2, rep(1, N-1), 1/2)
    k = 0:N
    f_points = -cos(k * pi / N)
    # Ensure that y_j points in [lower,upper] are usable
    f_points_scaled = (f_points + 1) / 2 * (upper - lower) + lower
  
    # For each point to eval, compute the polynomial approximation
    result = sapply(val, function(x) {
      
      # Scale x to the [-1,1] interval
      x <- 2 * (x - lower) / (upper - lower) - 1

      p_x = sum(weights * (-1)^k * fxn(f_points_scaled) / (x - f_points)) / sum(weights * (-1)^k / (x - f_points))
    })
    return(result)
  }
  
  
}

f <- function(x) {
  exp(x) * sin(x) * gamma(x) * exp(abs(x))
}

l = 1
u = 3
# laply(seq(1/2,8), function(x) {
  myfunc = chebyshev_compute(f, lower = l, upper = u, maxEfficiency = T, n_subdivisions = 1000)
  curve(f, l,u, lty = 2, col = "red", n = 300)
  curve(myfunc, l,u, add = T, col = "blue", lty = 4, pch = 2, n = 300)
  p = seq(l,u, length = 1000)
  plot(log10(abs(f(p) - myfunc(p))), type = "s", ylim = c(-15,0))
  abline(h = log10(1E-10))
# })


