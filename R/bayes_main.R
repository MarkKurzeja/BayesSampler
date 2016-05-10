#' For this file, the following is useful: library(roxygen2) library(devtools) 
#' From the package parent directory: Create a package with:
#' devtools::create("packagename") You can document it with
#' devtools::document("packagename") and documenting it allows the function to
#' be visible Install it with devtools::install("packagename") which will make
#' it usable to all. 
#' Download from github with: devtools::install_github("username/projectName")

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
        error("The chebyshev coefficients are not bounded under tol at maxN iterations")
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
      if (breakP == 1) {
        # We only are using the first three points to approximate it
        return(parabolaInterpolate(loopx, xvec[1:3], yvec[1:3]))
      } else if (breakP == N - 1) {
        # We are only using the last three points to approximate it
        indicies = seq(N-2,N)
        return(parabolaInterpolate(loopx, xvec[indicies], yvec[indicies]))
      }
      
      # Else, we are using two approximations to weight this point value
      weight = (loopx - xvec[breakP]) / (xvec[breakP + 1] - xvec[breakP])
      stopifnot(!(is.na(weight)))
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





#' TODO:
#' Add a function for BarycentricInterpolation(xvals, yvals, weightfunc)
#' Add a function for BarycentricInterpolation.Polynomial(xvals,yvals)
#' Add a function for BarycentricInterpolation.Chebyshev(xvals,yvals)
#' Use the new barycentric interpolation functions to simplify previous code for approximations
#' Debug the barycentric interpolation function
#' Add a function for FAST interpolation with BarycentricInterpolation.Berrut(xvals,yvals)
#' 
#' Bayes Stuff:
#' Add a function that computes, given only value of x's as input:
#'  gamma posterior for 
#'  lambda posterior for number of arrivals
#'  mean & var posterior
#'  var posterior for population
#'  pareto posterior for uniform
#'  dirlecht distribution for categorical variance
#'  beta posterior for probability 
#' 


#' Compute the Barycentric Interpolant for a given weight function
#' 
#' Take in a vector of x,y pairs and a weight function, and return a function that provides the barycentric interpolant through those points. This function is guaranteed to be rational 
#' @keywords rational interpolation interpolate
#' @param xvec A vector of the x values
#' @param yvec A vector of the y values
#' @param weights The weights to interpolate on
#' @export 
#' @examples  
#' f <- function(x) {
#'   sin(x)
#' }
#' 
#' x <- seq(0,1,length = 100)
#' # Caution, this interpolatnt is highly unstable!!!
#' curve(barycentricInterpolation(x, f(x), function(xvals, yvals) {return(-1)^length(x)})

BarycentricInterpolation.main <- function(xvals, yvals, weights){
  # Check the following:
  # The length of the xvals and yvals are the same
  # The length of the points is >= the degree of the rational approx
  stopifnot(length(xvals) == length(yvals))
  mymin = min(xvals)
  mymax = max(xvals)
  
  # Scale the points to [-2,2] for the weights for numerical stabality
  xvals <- (xvals - mymin) / (mymax - mymin) * 4 - 2
  
  # The weights are already computed

  # Build a function that returns the barycentric interpolant of the points
  # according to the rational function
  result <- function(x) {
    sapply(x, function(loopx) {
      # Ensure points that we are evaluating are not being extrapolated
      stopifnot(loopx >= mymin & x <= mymax)
      # Scale the x into [-2,2] interval for numerical stability. This trick was
      # borrowed from Barycentric interpolation of Chebyshev points, and it 
      # reduces the chance of overflow or underflow when evalutaing the 
      # polynomial for high degrees of N
      loopx <- (loopx - mymin) / (mymax - mymin) * 4 - 2
      # Check and see if the points is in xvals - if it is, return its y value
      if (any(loopx == xvals)) {
        return(yvals[loopx == xvals])
      }
      # Compute the Barycentric quotient to get our result
      sum(weights / (loopx - xvals) * yvals) / sum(weights / (loopx - xvals))
    })
  }
  return(result)
}


#' Floater-Hormann Rational Interpolation
#' 
#' Floater-Hormann interpolation, as suggested in "Barycentric rational 
#' interpolation with no poles and high rates of approximation Michael S. 
#' Floater & Kai Hormann (Equations 11 and 18)" is a method of interpolation 
#' that is far more stable than polynomial interpolation in equidistant or 
#' irregular points. While this function does benefit from points distributed 
#' near the end points such as Chebyshev points, it works well in almost any set
#' of points. It does this by moving out of the realm of polynomial functions, 
#' and instead, represents the function as a rational function (i.e. a 
#' polynomial of degree d <= N divided by another polynomial of degree f <= N). 
#' Using the Barycentric formula, this set of interpolates can be computed 
#' quickly once the weights have been computed. This function is guaranteed to 
#' have no poles on the real numbers, and thus wont suffer from issues with 
#' discontinuous functions.
#' 
#' The parameter d, the degree of the numerator and denominator, determines the 
#' convergence property of the interpolating rational function. The order of 
#' convergence is O(h^(d + 1)) if the function is smooth, and thus, for 
#' computational savings, select a small d (this is relative to N, the number of
#' points), and for very precise approximations, choose a value of d that is 
#' large relative to N. A caveat when choosing d is that for points distributed 
#' equidistant, only low values of d can be used else we will experience Runge 
#' phenomenon, and for points closer to Chebyshev, lower order approximations 
#' tend to produce inferior results. As such, when points are equidistant, keep 
#' the number of points <= 25 as a rule of thumb, and for points distributed 
#' closer to Chebyshev, take d = N for better results. d in [7,15] is a good
#' choice overall and generally yields good results for equidistant points
#' 
#' Because this function is an interpolate, it will be exactly equal to f(x) N+1
#' times, and thus, the precision is a function of the number of N+1 points 
#' primary where greater orders of approximation are honed with different values
#' of d. Choosing d = 3 uses third degree polynomials for the numerator and 
#' denominator, and thus the order of the error is O(h^4) - the same as cubic 
#' splines.
#' @param xvals The x values where f(x) is evaluated
#' @param yvals The values of f(x) at their respective xvals
#' @param d The degree of the polynomial. d <= N should be set to a low value 
#'   (<=10) for a general use function and set to d \\approx N for a high 
#'   precision evaluation
#' @export
#' @examples
#'  f <- function(x) {
#'    sin(x * 5) 
#'  }
#'  # Plot the function
#'  curve(f, -3, -1)
#'  
#'  # Set the graphical interface to a 2x1 grid
#'  par(mfrow = c(1,2))
#'  
#'  # Plot the {2,...,15} degree interpolants for d = N (precision oriented)
#'  # Do so through the equidistant points or through the chebyshev points
#'  # depending on which myxvals you comment out
#'  sapply(seq(2, 15, by = 1), function(len) {
#'    myxvals = seq(-3, -1,length = len)
#' #    myxvals = ((-cos(seq(0, len) * pi / len) + 1) / 2) * 2 - 3
#'    myyvals = f(myxvals)
#'    myfunc <- FloaterHormannInterpolation(myxvals, myyvals, d = len - 1)
#'    curve(f(x), -3, -1, col = "blue",  
#'          main = paste("Nodes:", len, "Interp. Result" ))
#'    abline(v = myxvals, col = "grey")
#'    curve(myfunc(x), -3, -1, add = T, col = "red", lty = 2)
#'    curve(log10(abs(myfunc(x) - f(x))), -3, -1, 
#'          col = "red", 
#'          n = 300, 
#'          main = paste("Error for", len, "nodes"), 
#'          ylim = c(-20,0))
#'  })

FloaterHormannInterpolation <- function(xvals, yvals, d) {
  # Checks - all others done in the barycentric formula
  stopifnot(length(xvals) == length(yvals))
  stopifnot(length(xvals) >= d)
  
  # Scale the values for the weight calculation to [-2,2] This prevents massive
  # over/underflow from happening when we are taking products later on
  bounds = range(xvals)
  xvals <- (xvals - bounds[1]) / (bounds[2]-bounds[1]) * 4 - 2
  # Compute the weights wk determined by equations 11 and 18, and we use an
  # alternative of 18 on the bottom of page 7 that allows us to compute the
  # product in a numerically stable way
  N = length(xvals)
  wk = sapply(seq(0, N - 1), function(k) {
    littleN = N - 1
    # Get the sequence of J_k which is a set from eq 11 and set it equal to iseq
    iseq = seq(0, littleN - d) %>% {.[k-d <= . & . <= k]}
    # For each i, compute the sumproduct
    (-1)^(k-d) * sum(sapply(iseq, function(i) {
        jseq = seq(i, i + d) %>% {.[. != k]}
        # We compute this product by taking the sum of the logs of the values
        # for numerical stability
        exp(sum(log(abs(1 / (xvals[k + 1] - xvals[jseq + 1])))))
     }))
  })
  
  # Transfrom back the xvals onto the desired range original range now that the
  # computation of the coefficients is over!
  xvals <- (xvals + 2) / (4) * (bounds[2] - bounds[1]) + bounds[1]
  # Return the barycentric interpolant with these weights!
  return(BarycentricInterpolation.main(xvals = xvals, yvals = yvals, weights = wk))
}

# Play code for testing out various forms of d and n for floater-hormann interpolation
# 
# par(mfrow = c(1,1))
# 
# plot(2,2,xlim = c(-1,1), ylim = c(-20,0), type = "n")
# f<- function(x) {
#   sin(50 * x)
# }
# 
# # curve(f, -1,1)
# # points(seq(-1,1,length = 100), f(seq(-1,1,length = 100)))
# 
# for (i in seq(30,500,by=50)) {
#   cat("i:", i, "\n")
# 
#   l = i
#   xvals <- seq(-1,1,length = l)
#   # xvals <- cos(seq(0,pi,length= l))
# 
#   yvals <- f(xvals)
# 
# 
#   d = 15#i - 1
#   cat("  d:", d, "\n")
#   myfunc <- FloaterHormannInterpolation(xvals, yvals, d)
#   # plot(xvals,yvals, main = paste("d:", d, "with", l, "points"))
#   # curve(myfunc, from = -1, to = 1, add = T)
# 
#   diff <- function(xvals) {
#     sapply(xvals, function(x) {
#       max(log10(abs(myfunc(x) - f(x))), -15)
#     })
#   }
# 
#   # curve(myfunc,-1,1,add = T, col = i, n = 300)
#   curve(diff, from = -1, to = 1, add = T, n = 399, col = i, main = paste("d:", d, "with", i, "points"))
#   # plotseq <- seq(-1,1,length = 100)
#   # wm <- plotseq[which.max(diff(plotseq))]
# 
#   # points(wm, diff(wm))
#   # cat("   Maxerror: ", diff(wm))
#   readline(prompt = "======")
# }
























