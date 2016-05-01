

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


#' Compute the Chebyshev Coefficients for a given function
#'
#' This function takes in a minimum and a maximum range, a function, and a
#' tolerance, and returns a the chebyshev coefficients of a function - useful
#' for determining the proper N to use. Sources: SIAM Chebyshev Expansions
#' EQ3.62 #' @keywords chebyshev approximation
#' @param fxn The function we wish to approximate
#' @param lower The lower bound of the approximation of the range
#' @param upper The upper bound of the approximation of the range
#' @param tol Any coefficient smaller than this will be disreguarded
#' @param maxN The maximum number of coefficients to calculate before the
#'   program errors and says there is no resolution to this approximation
#' @export
#' @examples
#' f <- function(x) {
#'    x^5
#' }
#' chebyshev_coef(f, lower = -2, upper = 2, tol = 1E-10)
chebyshev_coef <- function(fxn, lower, upper, tol = 1E-10, maxN = 10000) {
  N = 5
  # cat(sprintf("The first value of N is: %i\n", N))
  while (TRUE) {
    cx = cos(0:N * pi / N) # Get the chebyshev nodes
    shifted_cx = (cx + 1) / 2 * (upper - lower) + lower
    weights = c(1/2, rep(1,N-1), 1/2)
    # Apply the formula
    coef <- 2/N * sapply(0:N, function(k) {
      sum(weights * fxn(cx) * cos(k * 0:N * pi / N))
    })
    # Get the average of the last five coefficients
    tolCheck = mean(coef[seq(from = N-5, to = N)])
    # If the average of the last five coefficeints is less than our tol and we
    # have not exceeded our maxN
    if (tolCheck > tol) {
      if(N > maxN) {
        stop("The chebyshev coefficients are not bounded under tol at maxN iterations")
      }
      N = N * 3
    } else {
      return(coef)
    }
  }
}


#' Compute the value of N that resolves a Chebyshev Approximation
#'
#' This function takes in a minimum and a maximum range, a function, and a
#' tolerance, and returns a the chebyshev coefficients of a function - useful
#' for determining the proper N to use. Sources: SIAM Chebyshev Expansions EQ3.62
#' @keywords chebyshev approximation
#' @param fxn The function we wish to approximate
#' @param lower The lower bound of the approximation of the range
#' @param upper The upper bound of the approximation of the range
#' @param tol Any coefficient smaller than this will be disreguarded
#' @param maxN The maximum number of coefficients to calculate before the
#'   program errors and says there is no resolution to this approximation
#' @export
#' @examples
#' f <- function(x) {
#'    x^5
#' }
#' chebyshev_bestN(f, N = 1, lower = l, upper = u, tol = 1E-15)

chebyshev_bestN <- function(fxn, lower, upper, tol = 1E-10, maxN = 10000) {
  # Get the coefs from the chebyshev expansion and find where they drop off
  coefs = chebyshev_coef(fxn = fxn, lower = lower, upper = upper, tol = tol, maxN = maxN)
  # Use this code which does the following:
  # 1) Find out which of the coefficients are less than the tolerance
  # 2) Find which elements of 1:length(coefs) are not less than the tolerance
  # 3) Use the values from (2) to determine the index of the coefficient that
  #    is not
  bestN = max(which(!(1:length(coefs) %in% which(coefs < tol)))) + 1
  return(bestN)
}


f <- function(x) {
  exp(x)
}

l = -1
u = 1
  myfunc = chebyshev_coef(f,lower = l, upper = u)
  plot(log10(abs(myfunc)), type = "s")
  abline(v = chebyshev_bestN(f, lower = l, upper = u, tol = 1E-15))





