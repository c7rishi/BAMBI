#' `BAMBI`: An R package for Bivariate Angular Mixture Models
#'
#' `BAMBI` is an R package that provides functions for fitting
#' (using Bayesian methods) and simulating mixtures of univariate
#' and bivariate angular distributions. Please see the reference for a
#' detailed description of the functionalities of `BAMBI`.
#'
#' @references
#'
#' Chakraborty, S., & Wong, S. W. (2021). BAMBI: An R package for
#' fitting bivariate angular mixture models. *Journal of Statistical Software*.
#' \doi{10.18637/jss.v099.i11}.
#'
#' @name BAMBI
#' @md
NULL
#> NULL


#' @useDynLib BAMBI, .registration = TRUE
#' @import stats
#' @import graphics
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom utils tail txtProgressBar
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp evalCpp

NULL





#' @export
print.angmcmc <- function(x, ...) {

  # browser()
  output <-  paste("Dataset consists of", x$n.data, "observations.")


  if(grepl(x$method, "hmc")) {
    output[2] <- paste(x$ncomp, "component", x$model, "mixture fitted via HMC for model parameters.",
                       ifelse(x$n.chains == 1, "Number of chain = ", "Number of chains = "),
                       paste0(x$n.chains, "."))

    if(x$epsilon.random) {
      output[3] <- paste("epsilon chosen randomly at each iteration with average epsilon =",
                         paste0(format(x$epsilon, scientific=TRUE, digits = 2), collapse = ", "))
    } else {
        output[3] <- paste("epsilon fixed at", x$epsilon )
    }
    if(x$L.random){
      output[4] <- paste("L chosen randomly at each iteration with average L =",
                         paste(round(x$L, 2), collapse = ", "),
                         "across the ", x$n.chains, ifelse(x$n.chains == 1, "chain.", "chains."))
    } else {
      output[4] <- paste("L fixed at", paste(x$L, collapse = ", "),
                         "across the", x$n.chains, ifelse(x$n.chains == 1, "chain.", "chains."))
    }
    output[5] <- paste("acceptance rate for model parameters = ",
                       round(100*mean(x$accpt.modelpar), 2), "%.")
  }


   else if(grepl(x$method, "rwmh")) {

    output[2] <- paste(x$ncomp, "component", x$model, "mixture fitted via RWMH for model parameters.",
                       ifelse(x$n.chains == 1, "Number of chain = ", "Number of chains = "),
                       paste0(x$n.chains, "."))

    output[3] <- paste("proposals are independent normal with variances",
                       paste(format(x$propscale.final, scientific=TRUE, digits = 2), sep = "", collapse = ", "),
                       "for", paste(x$par.name[-1], sep = "", collapse=", "))

    output[4] <- paste("acceptance rate for model parameters = ",
                       round(100*mean(x$accpt.modelpar), 2), "%.")
  }

  output[5] <- paste("Number of iterations =", x$n.iter)
  cat(output, sep = "\n")
}


.onUnload <- function (libpath) {
  library.dynam.unload("BAMBI", libpath)
}

#' @export
print.stepfit <- function(x, ...)
{
  if(x$check_min) {
    output <- paste("Optimum component size =", x$ncomp.best)
    output[2] <- paste("Extract the best fit using the function \'bestmodel()\'")
    cat("\n")
    cat(output, sep = "\n")
  } else {
    warning(paste(toupper(x$crit), "did not attend a first minimum. Probably more clusters are needed."))
    cat("\n")
    cat("Extract the best fit using the function \'bestmodel()\'")
  }
}


#' Angular MCMC (\code{angmcmc}) Object
#' @description Checking for and creating an angmcmc object
#' @param object any R object
#' @param ... arguments required to make an angmcmc object. See details
#' @return logical. Is the input an angmcmc object?
#' @details
#' \code{angmcmc} objects are classified lists that are created when any of the five mixture model fitting
#' functions, viz., \code{fit_vmmix}, \code{fit_wnormmix}, \code{fit_vmsinmix}, \code{fit_vmcosmix} and
#' \code{fit_wnorm2mix} is used. An \code{angmcmc} object contains a number of elements, including the dataset, the
#' model being fitted on the dataset and dimension of the model (univariate or bivariate), the tuning parameters
#' used, MCMC samples for the mixture model parameters, the (hidden) component or cluster indicators for  data
#' points in each iteration and the (iteration-wise) log likelihood and log posterior density values (both calculated
#' upto some normalizing constants). When printed, an angmcmc object returns a brief summary of the function
#' arguments used to produce the object and the average acceptance rate of the proposals (in HMC and RWMH) used
#' over iterations. An \code{angmcmc} object can be used as an argument for the diagnostic and post-processing
#' functions available in \code{BAMBI} for making further inferences.
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' is.angmcmc(fit.vmsin.20)
#' @export

is.angmcmc <- function(object) {
  inherits(object, "angmcmc")
}

#' @rdname is.angmcmc
#' @export
angmcmc <- function(...) {
 ell <- list(...)
 class(ell) <- "angmcmc"
 ell
}
