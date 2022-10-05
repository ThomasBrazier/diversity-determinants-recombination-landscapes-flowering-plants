#============================================================================#
# The PDF of the Gamma function ----
#============================================================================#
# PDF of the Gamma function
g = function(x, alpha, beta) {
  (x^(alpha - 1)*exp(-x/beta))/(gamma(alpha)*beta^alpha)
}


#============================================================================#
# The PDF of the Gamma model ----
#============================================================================#
Gamma_Telomere_PDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, dx = 0.01) {
  # dx = 0.01
  x = seq(0, L, by = dx)
  # Return the PDF of the Gamma Telomere Model for given parameters
  # Equation of the PDF
  pdf = data.frame(x = x, pdf = (p * g(x, Alpha1, Beta1) + (1 - p) * g(L - x, Alpha2, Beta2)))
  return(pdf)
}

# plot(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L))

# Gamma_Telomere_PDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, dx = 0.01) {
#   # dx = 0.01
#   x = seq(0, L, by = dx)
#   # Return the PDF of the Gamma Telomere Model for given parameters
#   # Equation of the PDF
#   pdf = data.frame(x = x, pdf = (p * dgamma(x, Alpha1, Beta1) + (1 - p) * dgamma(L - x, Alpha2, Beta2)))
#   return(pdf)
# }
# 
# plot(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L))


#============================================================================#
# The CDF of the Gamma model ----
#============================================================================#
# So to get CDF from Probability Density Function(PDF), you need to integrate on PDF:
# fx <- Vectorize(fx)
# dx <- 0.01
# x <- seq(0, 10, by = dx)
# plot(x, cumsum(fx(x) * dx), type = "l", ylab = "cummulative probability", main = "My CDF")

# Another solution
# cdf
# Fx <- function(x, dx) {
#   cumsum(fx(x)*dx)
# }

# Numerical integration of the PDF to compute the CDF across the chromosome
# Take the same arguments as the PDF in argument and return the CDF
# Gamma_Telomere_CDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, dx = 0.01) {
#   # dx = 0.01
#   x = seq(0, L, by = dx)
#   # CDF function, integrate the PDF
#   # Scaled to be constrained between 0 and 1
#   CDF = function(dx) {
#     cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf*dx)/max(cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf*dx))
#   }
#   # Apply on values of x
#   cdf = data.frame(x = x , cdf = CDF(dx))
#   return(cdf)
# }


Gamma_Telomere_CDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, dx = 0.01) {
  library(pracma)
  x = seq(0, L, by = dx)
  # CDF function, integrate the PDF
  # Scaled to be constrained between 0 and 1
  CDF = function() {
    # cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf)/max(cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf), na.rm = TRUE)
    cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf*dx)
    # cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf*dx)/max(cumsum(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L, dx)$pdf*dx))
    
    # Use the formulae given by Mathematica
    # p - ((p * gammainc(Alpha1, (x/Beta1)))/gamma(Alpha1)) + (((-1 + p) * (gammainc(Alpha2, L/Beta2) - gammainc(Alpha2, (L-x)/Beta2)))/(gamma(Alpha2)))
  }
  # Apply on values of x
  cdf = data.frame(x = x , cdf = CDF())
  return(cdf)
}

#============================================================================#
# The Empirical PDF drawn from a vector X of data ----
#============================================================================#
# TODOne Implement the empirical CDF to fit on data
# Take the list X as argument to predict Y instead of creating a sequence from 0 to L by dx
Gamma_Telomere_EPDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, X) {
  # Return the EPDF of the Gamma Telomere Model for given parameters
  # Equation of the PDF
  pdf = data.frame(x = X, pdf = (p * g(X, Alpha1, Beta1) + (1 - p) * g(L - X, Alpha2, Beta2)))
  return(pdf)
}

# Gamma_Telomere_EPDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, X) {
#   # Return the EPDF of the Gamma Telomere Model for given parameters
#   # Equation of the PDF
#   pdf = data.frame(x = X, pdf = (p * dgamma(X, Alpha1, Beta1) + (1 - p) * dgamma(L - X, Alpha2, Beta2)))
#   return(pdf)
# }
#  DEBUG
# pdf = data.frame(x = X, pdf = (p * dgamma(X, Alpha1, Beta1) + (1 - p) * dgamma(L - X, Alpha2, Beta2)))
# plot(pdf)

#============================================================================#
# The Empirical CDF drawn from a vector X of data ----
#============================================================================#
# Gamma_Telomere_ECDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, X) {
#   # A vector of dx values
#   # As many dx than X values
#   dx = diff(X)
#   dx = c(0, dx) 
#   # CDF function, integrate the PDF
#   CDF = function(X) {
#     # A scaled CDF from 0 to 1
#     # Scaled to be constrained between 0 and 1
#     cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf*dx)/max(cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf*dx))
#   }
#   # Apply on values of x
#   cdf = data.frame(x = X, cdf = CDF(X))
#   return(cdf)
# }

Gamma_Telomere_ECDF = function(Alpha1, Alpha2, Beta1, Beta2, p, L, X) {
  # A vector of dx values
  # As many dx than X values
  dx = diff(X)
  dx = c(0, dx)
  # CDF function, integrate the PDF
  CDF = function(X) {
    # A scaled CDF from 0 to 1
    # Scaled to be constrained between 0 and 1
    # cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf)/max(cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf), na.rm = TRUE)
    cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf*dx)
    # cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf)/max(cumsum(Gamma_Telomere_EPDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$pdf))
  }
  # Apply on values of x
  cdf = data.frame(x = X, cdf = CDF(X))
  return(cdf)
}


#============================================================================#
# The objective function to minimize to fit the model ----
#============================================================================#
# Fit a CDF to data with a Sum of Squared Errors criterion
fitcdf = function(x) {
  # parameters to estimate
  Alpha1 = x[1]
  Alpha2 = x[2]
  Beta1 = x[3]
  Beta2 = x[4]
  p = x[5]
  # L is fixed to the empirical chromosome size in Mb
  L = max(X)
  
  # The fitted CDF
  fitted = Gamma_Telomere_ECDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X)$cdf
  # scaledfitted = fitted/max(fitted, na.rm = TRUE)
  # The goodness-of-fit criterion is the sum of squared differences
  goodness_fit = sum((Ys-fitted)^2, na.rm = TRUE) # Criterion: sum of squared differences
  
  return(goodness_fit)
}


#============================================================================#
# The iterative seach of the best set of parameters to fit the data ----
#============================================================================#
# Multiple starts from randomized initializations, take the best-of-N results.
# Look for the posterior mode with a random draw of multiple starting values?
# N (default = 1,000) iterations and selection of the set of parameters minimizing the value
# Explore the entire parameter space to search for the best local minima

iterative_fitcdf = function(N = 1000) {
  library(R.utils)
  # Search parameters among N random set of starting values
  # Very relaxed constraints on priors and boundaries of the function
  # Alpha diverge in O, so keep alpha >=1
  lower = c(1, 1, 0.2, 0.2, 0.1)
  upper = c(10, 10, 20, 20, 0.9)
  priorlower = c(1, 1, 0.5, 0.5, 0.2)
  priorupper = c(2, 2, 5, 5, 0.8)
  
  # Optim function
  optimpar = function() {
    parameters = c(runif(5, min = priorlower, max = priorupper))
    fit.optim = optimx(par = parameters, fitcdf, method = "L-BFGS-B", hessian = TRUE, lower = lower,
                       upper = upper)
    # Return results
    return(fit.optim)
  }
  
  try_optim = function(i) {
    # Try the optim function and do not abort if errors returned
    fit.optim = try(optimpar(), silent = TRUE)
    
    # Get results if try was successfull
    if (grepl("Error", fit.optim)[1] != 1) {
      fitoptim = unlist(fit.optim[1:6])
      # iterated_parameters$a1[i] = fit.optim$p1
      # iterated_parameters$a2[i] = fit.optim$p2
      # iterated_parameters$b1[i] = fit.optim$p3
      # iterated_parameters$b2[i] = fit.optim$p4
      # iterated_parameters$p[i] = fit.optim$p5
      # iterated_parameters$crit[i] = fit.optim$value
    } else {
      # Init results at NA if try failed
      fitoptim = rep(NA, 6)
      # iterated_parameters$a1[i] = NA
      # iterated_parameters$a2[i] = NA
      # iterated_parameters$b1[i] = NA
      # iterated_parameters$b2[i] = NA
      # iterated_parameters$p[i] = NA
      # iterated_parameters$crit[i] = NA
    }
    return(fitoptim)
  }
  
  library(pbmcapply)
  results = data.frame(matrix(unlist(pbmclapply(X = 1:N, function(x) try_optim(x))), nrow=N, byrow=TRUE),stringsAsFactors=FALSE)
  colnames(results) = c("a1", "a2", "b1", "b2", "p", "crit")

  # Select the best out of N runs
  best_parameter =  results[which.min(results$crit),]
  return(best_parameter)
}


#============================================================================#
# Fit a Marey map ----
#============================================================================#
# Fit a Marey map (species and chromosome) using iterative search of the best parameter
# and a Sum of Squares criterion of goodness-of-fit
# Returns the set of parameters

fit_Mareymap = function(set, chromosome, niter = 100, X, Ys, L) {
  # # Import data
  # data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
  # data.final$set = as.character(data.final$set)
  # data.final$map = as.character(data.final$map)
  # 
  # X = as.numeric(as.character(data.final$phys[data.final$set == set & data.final$map == chromosome]))/1000000
  # Y = data.final$gen[data.final$set == set & data.final$map == chromosome] 
  # # Sort X/Y by X (ascending) to have a continuous ascending function from 0 to L
  # Y = Y[order(X)]
  # X = X[order(X)]
  # # Work only on a scaled dataset for the recombination
  # Ys = Y/max(Y, na.rm = TRUE)
  # # Standard starting parameter
  # L = max(X)
  # Iterative fitting to find the best local minima out of N random sets of starting parameters
  param = iterative_fitcdf(N = niter)
  res = list(L = L, param = param)
  
  # Save the graph
  png(filename = paste("output/model_gamma_telomere/figures_fit/Gammatelomere_",
                       set,
                       "_chr", chromosome,
                       "_fit.png", sep = ""), width = 800, height = 600)
  # Plot of the observed ECDF with the fitted CDF (parameters fitted on the ECDF)
  plot(X, Ys, main = paste(set, "chromosome", chromosome, sep = " "),
       xlab = "Physical position (Mb)", ylab = "Genetic position (cM)")
  lines(Gamma_Telomere_ECDF(res$param$a1,  res$param$a2, res$param$b1, res$param$b2,
                            res$param$p, res$L, X), col = "Red", lwd = 2)
  lines(Gamma_Telomere_ECDF(res$param$a1,  res$param$a2, res$param$b1, res$param$b2,
                            res$param$p, res$L, seq(0, res$L, by = 0.1)),
        col = "Green", lty = "dashed", lwd = 2)
  # Add the putative centromere position, p
  abline(v = res$param$p*res$L)
  dev.off()
  
  return(res)
}




#============================================================================#
# END ----
#============================================================================#