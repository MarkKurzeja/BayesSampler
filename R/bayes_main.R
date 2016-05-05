

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
      cat("Can't be resolved in maxN points, function evaluation stops past N > maxN")
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
#' @param verbose logical; if true, the chebyshev_coef function will print a
#'   verbose message output indicating the coefficient sizes and the amount of
#'   digits (base 10) that each coefficient represents
#' @export
#' @examples
#' f <- function(x) {
#'    x^5
#' }
#' chebyshev_coef(f, lower = -2, upper = 2, tol = 1E-10)
chebyshev_coef <- function(fxn, lower, upper, tol = 1E-10, maxN = 10000, verbose = F) {
  N = 5
  while (TRUE) {
    if (verbose) cat(sprintf("Trying N = %*i\n", 5, N))
    cx = cos(0:N * pi / N) # Get the chebyshev nodes
    shifted_cx = (cx + 1) / 2 * (upper - lower) + lower
    weights = c(1/2, rep(1,N-1), 1/2)
    # Apply the formula
    if (verbose) cat(sprintf("   Computing the coefficients..."))
    coef <- 2/N * sapply(0:N, function(k) {
      sum(weights * fxn(shifted_cx) * cos(k * 0:N * pi / N))
    })
    if (verbose) cat(sprintf("done\n"))
    # Get the average of the last five coefficients
    tolCheck = mean(coef[seq(from = N-1, to = N)])
    if (verbose) cat(sprintf("   Tolerance: %.20f\n", abs(tolCheck)))
    if (verbose) cat(sprintf("   Tolerance Digits: %.3f\n", -1 * log10(abs(tolCheck))))
    # If the average of the last five coefficeints is less than our tol and we
    # have not exceeded our maxN
    if (abs(tolCheck) > tol) {
      N = N * 3
       if(N > maxN) {
        stop("The chebyshev coefficients are not bounded under tol at maxN iterations")
      }
    } else {
      return(coef)
    }
  }
}


#' Compute the value of N that resolves a Chebyshev Approximation
#' 
#' This function takes in a minimum and a maximum range, a function, and a 
#' tolerance, and returns a the chebyshev coefficients of a function - useful 
#' for determining the proper N to use. Sources: SIAM Chebyshev Expansions
#' EQ3.62. It returns the best value of N and the coefficients that produce it
#' @keywords chebyshev approximation
#' @param fxn The function we wish to approximate
#' @param lower The lower bound of the approximation of the range
#' @param upper The upper bound of the approximation of the range
#' @param tol Any coefficient smaller than this will be disreguarded
#' @param maxN The maximum number of coefficients to calculate before the 
#'   program errors and says there is no resolution to this approximation
#' @param verbose logical; if true, the chebyshev_coef function will print a
#'   verbose message output indicating the coefficient sizes and the amount of
#'   digits (base 10) that each coefficient represents
#' @export
#' @examples
#' f <- function(x) {
#'    x^5
#' }
#' chebyshev_bestN(f, N = 1, lower = l, upper = u, tol = 1E-15)

chebyshev_bestN <- function(fxn, lower, upper, tol = 1E-10, maxN = 10000, verbose = F) {
  # Get the coefs from the chebyshev expansion and find where they drop off
  coefs = chebyshev_coef(fxn = fxn, lower = lower, upper = upper, tol = tol, maxN = maxN, verbose = verbose)
  # Use this code which does the following:
  # 1) Find out which of the coefficients are less than the tolerance
  # 2) Find which elements of 1:length(coefs) are not less than the tolerance
  # 3) Use the values from (2) to determine the index of the coefficient that
  #    is not
  bestN = max(which(!(1:length(coefs) %in% which(abs(coefs) < tol)))) + 1
  coefs = coefs[1:bestN]
  return(list(bestN = bestN, coefs = coefs))
}

#' Approximate a function with a Chebyshev Series
#' 
#' This function takes in a minimum and a maximum range, a function, and a 
#' tolerance, and returns a chebyshev series evalued based on the coefficients
#' of the series. This function is actually faster than the interpolation
#' function for a simular value of N
#' Sources: SIAM Chebyshev Expansions: EQ3.62
#' for coefficients Algorithm 3.2 for evaluating the sum (Modified Clenshaws
#' algorithm) EQ3.70 for error bounds
#' @keywords chebyshev approximation
#' @param fxn The function we wish to approximate
#' @param lower The lower bound of the approximation of the range
#' @param upper The upper bound of the approximation of the range
#' @param tol Any coefficient smaller than this will be disreguarded
#' @param maxN The maximum number of coefficients to calculate before the 
#'   program errors and says there is no resolution to this approximation
#' @param verbose logical; if true, the chebyshev_coef function will print a 
#'   verbose message output indicating the coefficient sizes and the amount of 
#'   digits (base 10) that each coefficient represents
#' @export
#' @examples
#' f <- function(x) {
#'    x^5
#' }
#' chebfun(f, N = 1, lower = l, upper = u, tol = 1E-15)


chebfun <- function(fxn, lower, upper, tol = 1E-10, maxN = 10000, verbose = F) {
  # Get the coefficients and truncate them to the largest N
  coef_list = chebyshev_bestN(fxn = fxn, lower = lower, upper = upper, tol = tol, maxN = maxN)
  N = coef_list$bestN - 1
  coef = coef_list$coefs
  
  result <- function(x) {
    # Ensure x is in [lower, upper]
    if (any(x < lower) | any(x > upper)) {
      stop("x is not in [lower,upper] - extrapolation is not stable")
    }
    
    # Scale the x such that they are in [-1,1]
    x <- (x - lower) / (upper - lower) * 2 - 1
    
    y_val <- sapply(x, function(loopx) {
      weights = c(1/2, rep(1,N - 1), 1/2) 
      # browser()
      sum(weights * coef * cos(0:N * acos(loopx)))
    })
    return(y_val)
  }
    return(result)
}
  
#' Meromorphic function for function weighting 
#' 
#' Meromorphic Sinosoid function that makes this possible Very cool properties - bounded, parameterized, monotonic, C^Inf ' Sources: Ulrich Multz Documentation on using this function for a parabolic interpolation program 
#' @keywords parabola meromorphic weighting weight
#' @param x The input value to the function x \in [0,1] '
#' @export 
#' @examples  
#' meroSFunc(0:100 / 100)
meroSFunc <- function(x) {
  # Set a value of p to 0.810747 which minimizes the L2 norm of the second
  # derivative
  p = 0.810747

  # Vectorize the function
  sapply(x, function(loopx) {
    if (loopx <= 0)
    { 
      return(0)
    } else if(loopx >= 1){
      return(1)
    }
    
    return(1 / (1 + exp(p * (1 / loopx + 1 / (loopx - 1)))))
  })
};

curve(meroSFunc,-1/2,3/2)


#' Interpolate a function with a parabola 
#' 
#' Take in a vector of x values and a vector of y values for a parabola, three in total, and return the parabola approximation at the value of x 
#' @keywords parabola interpolation interpolate
#' @param x The input value to the function
#' @param xvec A vector of length three of the x values
#' @param yvec A vector of length three of the y values
#' @export 
#' @examples  
#' parabolaInterpolate(0.5, c(1,0,4), c(5,6,2))
parabolaInterpolate <- function(x, xvec, yvec) {
  sapply(x, function(loopx) {
    yvec[1] * (loopx - xvec[2])*(loopx - xvec[3]) / ((xvec[1] - xvec[2])*(xvec[1] - xvec[3])) +
    yvec[2] * (loopx - xvec[1])*(loopx - xvec[3]) / ((xvec[2] - xvec[1])*(xvec[2] - xvec[3])) +
    yvec[3] * (loopx - xvec[1])*(loopx - xvec[2]) / ((xvec[3] - xvec[1])*(xvec[3] - xvec[2]))
  })
}

#' C^Inf Parabola Smoother - Ulrich Multze
#' 
#' Take in a vector of x and y values and return the meromorphic weighting of parabola splines as suggested by Ulrich Multz. This function is resulting C^Inf meaning all derivatives are continuious!
#' This funcion returns a function that allows you to compute this for an arbitrary x
#' @keywords parabola interpolation interpolate
#' @param x The input value to the function
#' @param xvec A vector of the x values
#' @param yvec A vector of the y values
#' @export 
#' @examples  
#' This example aims to show the boundary polynomials that are used in this 
#' method of interpolation. The black line is the meromorphic weight of the two
#' interpolating polynomials that are to the left and to the right of the 
#' given value of x
#' xvec = 1:10
#' yvec = runif(10) * x
#' Get the interpolating function
#' myfunc <- CInfParabolaSmoother(xvec,yvec)
#' Plot the curve and the points
#' curve(myfunc, 1,10)
#' points(cbind(xvec,yvec), pch = 3)
#' Draw in the boundary polynomials
#' for (i in 2:9) {
#' curve(parabolaInterpolate(x, xvec[seq(i-1,i+1)],yvec[seq(i-1,i+1)]), xvec[i-1], xvec[i + 1], add = T, col = "blue", lty = 2)
#' }

CInfParabolaSmoother <- function(xvec, yvec) {
  # Check the lengths to ensure that they are the same
  stopifnot(length(xvec) == length(yvec))
  # Check to see that we have enough points to make this work
  stopifnot(length(xvec) >= 4)
  
  # Sort the yvec and the xvec to ensure stable ordering
  myorder = order(xvec)
  xvec <- xvec[myorder]
  yvec <- yvec[myorder]
  N = length(xvec)
  
  # Create the vectorized function
  result <- function(x) {
    sapply(x, function(loopx) {
      # Check the boundary conditions
      stopifnot(loopx >= min(xvec))
      stopifnot(loopx <= max(xvec))
      
      # If the value is equal to any of the inputs, then just return the input
      if (any(loopx == xvec)) {
        return(yvec[loopx == xvec])
      }
      
      # Find the breaking point
      breakP <- max(which(xvec <= loopx))
      cat(sprintf(">> Loopx at:      %.4f\n", loopx))
      cat(sprintf("   Breakpoint at: %i -> [%0.0f,%0.0f]\n", breakP, xvec[breakP], xvec[breakP + 1]))
      if (breakP == 1) {
        # We only are using the first three points to approximate it
        return(parabolaInterpolate(loopx, xvec[1:3], yvec[1:3]))
      } else if (breakP == N - 1) {
        # We are only using the last three points to approximate it
        indicies = seq(N-2,N)
        return(parabolaInterpolate(loopx, xvec[indicies], yvec[indicies]))
      }
      
      # Else, we are using two approximations to weight these
      # browser()
      weight = (loopx - xvec[breakP]) / (xvec[breakP + 1] - xvec[breakP])
      if (is.na(weight)) {
        browser()
      }
      w.first = meroSFunc(weight)
      w.second = 1 - w.first
      # Get the left and the right parabola approximations centered on the
      # points breakP and breakP + 1
      leftIndicies = seq(breakP-1, breakP + 1)
      rightIndicies = seq(breakP, breakP + 2)
      leftApprox = parabolaInterpolate(loopx, xvec[leftIndicies], yvec[leftIndicies])
      rightApprox = parabolaInterpolate(loopx, xvec[rightIndicies], yvec[rightIndicies])
      
      # Return the approximation
      return(w.second * leftApprox + w.first * rightApprox)
    })
  }
}








