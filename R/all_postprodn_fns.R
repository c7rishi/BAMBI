#' Add (extra) burnin and thin to angmcmc object after original run
#' @param object angmcmc object
#' @inheritParams fit_angmix
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' lpdtrace(fit.vmsin.20)
#' # Now add extra burn-in
#' fit.vmsin.20.burn <- add_burnin_thin(fit.vmsin.20, 0.3)
#' lpdtrace(fit.vmsin.20.burn)
#' @export
add_burnin_thin <- function(object, burnin.prop=0, thin=1)
{
  if (!is.angmcmc(object))
    stop("object must be an \'angmcmc\' object")
  if(burnin.prop < 0 || burnin.prop >= 1)
    stop("\"burnin.prop\" must be in [0, 1)")
  if(thin < 1)
    stop("\"thin\" must be a positive integer")


  # browser()
  final_iter.orig <- object$final_iter
  final_iter.burn <- setdiff(final_iter.orig,
                             final_iter.orig[seq_len(ceiling(length(final_iter.orig)*burnin.prop))])

  thin <- round(thin)
  thin.final <- object$thin*thin
  final_iter.thin <- final_iter.burn[c(TRUE, rep(FALSE, thin.final-1))]

  object$n.burnin <- final_iter.thin[1] - 1
  object$thin <- thin.final
  object$n.iter.final <- length(final_iter.thin)
  object$final_iter <- final_iter.thin

  object
}

#' Select chains from angmcmc objects
#' @inheritParams pointest
#' @param chain.no labels of chains to be retained in the final sample. If missing,
#' all chains are used.
#' @param ... unused
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              L = c(10, 12), chains_parallel = FALSE,
#'                              n.chains = 2)
#' fit.vmsin.20
#' fit.vmsin.20.1 <- select_chains(fit.vmsin.20, 1)
#' fit.vmsin.20.1
#'
#' @return Returns another angmcmc object with only selected chains passed through \code{chain.no}
#' @export

select_chains <- function(object, chain.no, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  if (missing(chain.no)) {
    chain.no <- 1:object$n.chains
  } else if (any(!chain.no %in% 1:object$n.chains)) {
    stop("invalid chain number")
  }

  object$par.value <- object$par.value[, , , chain.no, drop=FALSE]
  object$clus.ind <- object$clus.ind[, , chain.no, drop=FALSE]
  object$llik.contri <- object$llik.contri[, , chain.no, drop=FALSE]
  object$mem.prob <- object$mem.prob[, , , chain.no, drop=FALSE]
  object$llik <- object$llik[,  chain.no, drop=FALSE]
  object$lpd <- object$lpd[, chain.no, drop=FALSE]
  object$lprior <- object$lprior[, chain.no, drop=FALSE]
  object$accpt.modelpar <- object$accpt.modelpar[, , chain.no, drop=FALSE]
  if (object$return_tune_param)
    object$tune_param <- object$tune_param[, , chain.no, drop=FALSE]
  object$epsilon <- object$epsilon[, chain.no, drop=FALSE]
  object$L <- object$L[chain.no]
  object$propscale <- object$propscale[, chain.no, drop=FALSE]
  object$n.chains <- length(chain.no)

  object
}

#' Create an mcmc.list object from an angmcmc object
#' @param x angmcmc object
#' @param ... unused
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#'
#' # now convert it to mcmc.list
#' library(coda)
#' fit.vmsin.20.mcmc <- as.mcmc.list(fit.vmsin.20)
#' @importFrom coda as.mcmc.list
#' @method as.mcmc.list angmcmc
#'
#' @export
as.mcmc.list.angmcmc <- function(x, ...)
{
  if (!is.angmcmc(x))
    stop("\'x\' must be an angmcmc object")

  object <- x
  npars <- length(object$par.name)
  ncomp <- object$ncomp
  n.iter <- object$n.iter
  n.data <- object$n.data
  n.chains <- object$n.chains

  if (object$ncomp == 1) {
    par_names_long <- object$par.name
  } else {
    par_names_long <- apply(as.matrix(expand.grid(object$par.name, 1:object$ncomp)),
                            1, function(x) paste0(x[1], "[", x[2], "]"))
  }


  start <- object$final_iter[1]
  end <- object$n.iter
  thin <- object$thin

  out <- vector("list", n.chains)

  for(ii in 1:n.chains) {
    par.vals <- object$par.value[, , , ii, drop=FALSE]
    x <- t(matrix(c(par.vals), ncomp*npars, n.iter))
    colnames(x) <- par_names_long
    if (object$ncomp == 1) x <- x[, -1, drop = FALSE]
    out[[ii]] <- coda::mcmc(x, start = start, end = end, thin = thin)
  }

  coda::as.mcmc.list(out)
}




#' Fix label switching in angmcmc objects
#' @inheritParams pointest
#'
#' @param ... arguments other than \code{z, K, complete, mcmc, p}
#' and \code{data} passed to \link[label.switching]{label.switching}. See details.
#'
#' @details \code{fix_label} is a wrapper for  \link[label.switching]{label.switching} from
#' package \code{label.switching} for \code{angmcmc} objects. The arguments
#' \code{z, K, complete, mcmc, p} and \code{data} are appropriately filled in
#' from \code{object}. The \code{label.switching} argument \code{method} can
#' be a scalar or vector; for this wrapper it defaults to \code{"STEPHENS"} if the \code{angmcmc} was
#' created with permutation sampling (by setting perm_sampling = TRUE in
#' \link{fit_angmix}), and to \code{"DATA-BASED"} otherwise.
#'
#' @return Returns a single \code{angmcmc} object or a list of \code{angmcmc} objects (according as whether
#' the argument \code{method} is a scalar or vector) with label switchings corrected (after burn-in and thin)
#' according to the resulting permutation from \link[label.switching]{label.switching}.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # now apply fix_label
#' fit.vmsin.20.fix <- fix_label(fit.vmsin.20)
#'
#' @export

fix_label <- function(object, ...) {

  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  if (object$ncomp == 1)
    stop("Fitted model has only one component")


  # first collect all necessary arguments from all chains,
  # and collapse into a wide matrix matrix
  K <- ncomp <- object$ncomp
  n.data <- object$n.data
  n.iter.final <- object$n.iter.final
  n.chains <- object$n.chains

  final_iter <- object$final_iter

  npar <- length(object$par.name)

  ell <- list(...)

  if (is.null(ell$method)) {
    if (object$perm_sampling)
      ell$method <- "STEPHENS"
    else
      ell$method <- "DATA-BASED"
  }


  method <- ell$method

  iter_no <- matrix(1:(n.chains*n.iter.final), ncol = n.chains)

  z <- matrix(1, n.chains*n.iter.final, n.data)
  p <- array(1, c(n.chains*n.iter.final, n.data, ncomp))

  mcmc_par <- array(1, c(n.chains*n.iter.final, ncomp, npar))

  for(ii in 1:n.chains) {
    z[iter_no[, ii], ] <- t(object$clus.ind[, final_iter, ii])
    p[iter_no[, ii], , ] <- aperm(object$mem.prob[, , final_iter, ii], c(3, 1, 2))
    mcmc_par[iter_no[, ii], , ] <- aperm(object$par.value[, , final_iter, ii], c(3, 2, 1))
  }

  par_long <- aperm(mcmc_par, c(3, 2, 1))
  mem_prob_long <- aperm(p, c(2, 3, 1))
  clus_ind_long <- t(z)

  # browser()

  complete <- function(x, z, pars) {
    # x: data (size = n)
    # z: allocation vector (size = n)
    # pars: K x npars matrix of normal mixture parameters:
    # pars[k, 1] = k-th mixing proportion
    # pars[k, -1] = k-th component parameters
    pars.t <- t(pars)
    rownames(pars.t) <- object$par.name
    K <- ncol(pars.t)

    # pars.t[1, ] : mixing proportions
    # pars.t[-1, ] : model pars

    if (object$type == "bi")
      data_groups <- lapply(1:K, function(j) x[which(z == j), ])
    else
      data_groups <- lapply(1:K, function(j) x[which(z == j)])

    n.groups <- listLen(data_groups)
    res <- 0
    for(j in 1:K) {
      inargs <- list_by_row(pars.t[-1, j, drop = FALSE])
      inargs$x <- data_groups[[j]]
      if (object$model %in% c("wnorm", "wnorm2"))
        inargs$int.displ <- object$int.displ
      else if (object$model == "vmcos")
        inargs$qrnd_grid <- object$qrnd_grid
      res <- res + sum(
        signif(
          log(
            do.call(paste0("d", object$model), inargs)
          ),
          8
        )
      )
    }

    res + sum(n.groups * log(pars.t[1, ]))
  }

  labswitch_in <- addtolist(ell,
                            z = z, K = ncomp,
                            complete = complete,
                            mcmc = mcmc_par,
                            p = p, data = object$data)


  # lab_switch <-
  #   label.switching::label.switching(method = method,
  #                                    z = z, K = ncomp,
  #                                    complete = complete,
  #                                    mcmc = mcmc_par,
  #                                    p = p, data = object$data,
  #                                    ...)

  lab_switch <-
    do.call(label.switching::label.switching, labswitch_in)



  # now rearrage the component labels according to the rearrangement
  # obtained in lab_switch, and revert back to the original (1 higer dim
  # array) form of the outputs

  res <- vector("list", length(method))

  for(j in 1:length(method)) {
    tmp <- object
    method_curr <- method[j]
    perm_mat <- lab_switch[["permutations"]][[j]]

    par.value.mod <- tmp$par.value[, , final_iter, , drop=FALSE]
    mem.prob.mod <- tmp$mem.prob[, , final_iter, , drop=FALSE]
    clus.ind.mod <- tmp$clus.ind[, final_iter, , drop=FALSE]
    if (object$return_tune_param)
      tune_param_mod <- tmp$tune_param[, final_iter, , drop=FALSE]

    # relabel on these subarrays
    for(ii in 1:n.chains) {
      for(iter in 1:n.iter.final) {
        # browser()
        par.value.mod[ , , iter, ii] <- par.value.mod[, perm_mat[iter_no[iter, ii], ], iter, ii]
        mem.prob.mod[, , iter, ii] <- mem.prob.mod[, perm_mat[iter_no[iter, ii], ], iter, ii]
        clus.ind.mod[, iter, ii] <- perm_mat[iter_no[iter, ii], clus.ind.mod[, iter, ii]]
        if (object$return_tune_param) {
          tume_param_temp_mat <- matrix(tune_param_mod[, iter, ii], ncol=ncomp)
          tune_param_mod[, iter, ii] <- c(tume_param_temp_mat[, perm_mat[iter_no[iter, ii], ]])
        }
      }
    }

    tmp$par.value[, , final_iter, ] <- par.value.mod
    tmp$mem.prob[, , final_iter, ] <- mem.prob.mod
    tmp$clus.ind[, final_iter, ] <- clus.ind.mod
    if (object$return_tune_param)
      tmp$tune_param[, final_iter, ] <- tune_param_mod

    tmp$fixed.label <- TRUE
    res[[j]] <- tmp
  }

  if (length(method) == 1)
    res <- res[[1]]

  res
}

#' Point estimates for parameters from an angmcmc object
#' @param object angular MCMC object.
#' @param fn function, or a single character string specifying its name, to evaluate on MCMC samples to estimate
#' parameters.  Defaults to \code{mean}, which computes the estimated posterior mean.
#' Note that if \code{fn = "MODE"} (warning: not \code{"mode"}) or \code{fn = "MAP"}, then the maximum aposteriori estimate (MAP) is
#' calculated.
#' @param par.name vector of names of parameters for which point estimates are to be computed.  If \code{NULL}, results for all parameters are provided.
#' @param comp.label vector of component labels (positive integers, e.g., \code{1, 2, ...}) for which point estimates are to be computed.
#' If \code{NULL}, results for all components are provided.
#' @param chain.no vector of chain numbers whose samples are to be be used.
#' in the estimation. By default all chains are used.
#' @param ... additional arguments to be passed to the function.
#'
#' @return Returns a matrix of point estimates, or vector of point estimates if \code{length(par.name)==1} or \code{length(comp.label)==1}.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # estimate parameters by sample mean
#' (est_mean <- pointest(fit.vmsin.20))
#' # estimate parameters by sample median
#' (est_median <- pointest(fit.vmsin.20, fn = median))
#' # estimate parameters by MAP
#' (est_median <- pointest(fit.vmsin.20, fn = "MODE"))
#' @export

pointest <- function(object, fn = mean, par.name,
                     comp.label, chain.no, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")


  if (is.character(fn)) {
    if (fn == "MODE" | fn == "MAP") {
      do_MAP <- TRUE
    }
    else {
      do_MAP <- FALSE
      fn <- match.fun(fn)
    }
  } else {
    do_MAP <- FALSE
    fn <- match.fun(fn)
  }


  ell <- list(...)
  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if (missing(par.name)) {
    par.name <- object$par.name
  }  else if (any(!par.name %in% object$par.name)) {
    stop("invalid par.name")
  }

  if (missing(comp.label)) {
    comp.label <- 1:object$ncomp
  } else if (any(!comp.label %in% 1:object$ncomp)) {
    stop("invalid component label")
  }

  if (missing(chain.no)) {
    chain.no <- 1:object$n.chains
  } else if (any(!chain.no %in% 1:object$n.chains)) {
    stop("invalid chain number")
  }

  par.posit <- which(object$par.name %in% par.name)

  final_iter <- object$final_iter


  if (!do_MAP) {
    res <- apply(object$par.value[par.posit, comp.label, final_iter, chain.no, drop = FALSE],
                 c(1, 2), fn, drop = FALSE)
  }

  else {
    map_indic <- which(object$lpd[, chain.no, drop=FALSE]
                       == max(object$lpd[final_iter, chain.no]),
                       arr.ind = TRUE)[1, ]  # first entry of max
    res_mat <- as.matrix(object$par.value[ , , map_indic[1], map_indic[2]])

    res <- res_mat[par.posit, comp.label, drop=FALSE]
  }

  rownames(res) <- par.name
  colnames(res) <- comp.label
  res
}


#' Quantile estimates for parameters from an angmcmc object
#' @inheritParams pointest
#' @param x angmcmc object
#' @inheritParams stats::quantile
#' @return Returns a three dimensional array of quantiles, or a matrix (vector) of quantiles
#'  if one (or two) among \code{par.name},  \code{comp.label}, \code{probs} has length 1.
#' @param ... further arguments to pass to \code{quantile}.  In particular, \code{probs = seq(0, 1, 0.25)}
#' is the default vector of quantiles computed for each parameter.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # 0.025th quantiles
#' (quant_025 <- quantile(fit.vmsin.20, prob = 0.025))
#' # 0.975th quantiles
#' (quant_975 <- quantile(fit.vmsin.20, prob = 0.975))
#' # default quantiles
#' (quant_def <- quantile(fit.vmsin.20))
#'
#' @export


quantile.angmcmc <- function(x, par.name, comp.label, chain.no,
                             probs = seq(0, 1, 0.25), ...)
{

  object <- x

  ell <- list(...)

  if (!is.angmcmc(object))
    stop("\'x\' must be an angmcmc object")

  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if (missing(par.name)) {
    par.name <- object$par.name
  }  else if (any(!par.name %in% object$par.name)) {
    stop("invalid par.name")
  }

  if (missing(comp.label)) {
    comp.label <- 1:object$ncomp
  } else if (any(!comp.label %in% 1:object$ncomp)) {
    stop("invalid component label")
  }

  if (missing(chain.no)) {
    chain.no <- 1:object$n.chains
  } else if (any(!chain.no %in% 1:object$n.chains)) {
    stop("invalid chain number")
  }

  par.posit <- which(object$par.name %in% par.name)

  final_iter <- object$final_iter

  n_probs <- length(probs)

  res <- list()

  for(jj in 1:n_probs) {
    res[[jj]] <- apply(object$par.value[par.posit, comp.label, final_iter, chain.no, drop = FALSE],
                       c(1, 2),
                       function(y) quantile(y, probs = probs[jj]))
    rownames(res[[jj]]) <- par.name
    colnames(res[[jj]]) <- comp.label
  }

  names(res) <- paste0(round(probs*100), "%")

  res
}


#' Extract MCMC samples for parameters from an angmcmc object
#' @inheritParams pointest
#' @param object angular MCMC object
#' @param drop logical. Should the dimension of the output be dropped, if \code{par.name},
#' \code{comp.label} or \code{chain.no} has a single level?
#' @details The default for both \code{par.name} and \code{comp.label} are the all possible choices
#' available in \code{object}.
#' @return
#' Returns  a four dimensional array with
#'
#' dimension 1 - model parameters and mixing proportions
#' dimention 2 - components
#' dimension 3 - MCMC iterations
#' dimension 4 - chain number
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # extract Markov chain realizations for kappa1 from component 1
#' extr_kappa1_1 <- extractsamples(fit.vmsin.20, "kappa1", 1)
#' # for kappa1 from component from all components
#' extr_kappa1 <- extractsamples(fit.vmsin.20, "kappa1")
#' # for all parameters in component 1
#' extr_1 <- extractsamples(fit.vmsin.20, comp.label = 1)
#'
#' @export

extractsamples <- function(object, par.name, comp.label,
                           chain.no, drop = TRUE, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  ell <- list(...)

  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if (missing(par.name)) {
    par.name <- object$par.name
  }  else if (any(!par.name %in% object$par.name)) {
    stop("invalid par.name")
  }

  if (missing(comp.label)) {
    comp.label <- 1:object$ncomp
  } else if (any(!comp.label %in% 1:object$ncomp)) {
    stop("invalid component label")
  }

  if (missing(chain.no)) {
    chain.no <- 1:object$n.chains
  } else if (any(!chain.no %in% 1:object$n.chains)) {
    stop("invalid chain number")
  }

  par.posit <- which(object$par.name %in% par.name)
  final_iter <- object$final_iter
  out <- object$par.value[par.posit, comp.label, final_iter, chain.no, drop = drop]

  # browser()
  if (!drop)
    dimnames(out) <- list(par.name, comp.label, NULL, chain.no)

  out
}


#' Summary statistics for parameters from an angmcmc object
#' @inheritParams base::summary
#' @inheritParams pointest
#' @param object angular MCMC object.
#' @details Computes (after thinning and discarding burn-in) point estimates with 95\% posterior credible sets for all components and all parameters,
#' together with the sample averages of log likelihood and log posterior density.
#' @return Returns a list with elements \code{estimate, lower, upper, llik} and \code{lpd}. \code{estimate}
#' is itself a list with three elements: \code{mean}, \code{median} and \code{mode} providing the
#' sample mean, sample median and (sample) MAP estimates.
#'
#' Note that \code{summary.angmcmc} has its own print method, providing a table the estimated mean and 95\% credible intervals for each parameter
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' summary(fit.vmsin.20)
#'
#' @export

summary.angmcmc <- function(object, par.name, comp.label,
                            chain.no, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  ell <- list(...)

  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if (missing(par.name)) {
    par.name <- object$par.name
  }  else if (any(!par.name %in% object$par.name)) {
    stop("invalid par.name")
  }

  if (missing(comp.label)) {
    comp.label <- 1:object$ncomp
  } else if (any(!comp.label %in% 1:object$ncomp)) {
    stop("invalid component label")
  }

  if (missing(chain.no)) {
    chain.no <- 1:object$n.chains
  } else if (any(!chain.no %in% 1:object$n.chains)) {
    stop("invalid chain number")
  }


  est <- pointest(object, mean, par.name, comp.label, chain.no)
  map <- pointest(object, "MODE", par.name, comp.label, chain.no)
  quants <- quantile(object, par.name, comp.label, chain.no,
                     probs = c(0.5, 0.025, 0.975))

  med <- quants[[1]]
  low <- quants[[2]]
  up <- quants[[3]]

  llik <- as.numeric(logLik(object))
  lpd <- mean(object$lprior) + llik

  res <- list("estimate " = list(mean = est, median = med, mode = map),  "upper" = up, "lower" = low,
              "llik" = llik, "lpd" = lpd)

  class(res) <- "summary_angmcmc"
  res
}

#' @export
print.summary_angmcmc <- function(x, ...)
{
  res <- x
  print(est_ci(res$estimate$mean, res$lower, res$upper, digits=2))
}



# #' AIC and BIC for angmcmc objects -- not needed
# #' @inheritParams stats::AIC
# #' @inheritParams logLik.angmcmc
# #' @param ... additional argument to be passed to \link{logLik.angmcmc}
# #' @return AIC computes the AIC and BIC computes BIC for \code{angmcmc} objects.
# #'
# #' @details
# #' Note that \code{AIC.angmcmc} and \code{BIC.angmcmc} calls \link{logLik.angmcmc},
# #' which calculates Bayes estimate of the log-likelihood and *not* the maximum
# #' likelihood. As such, care needs to be taken while using theses quantities.
# #'
# #' \eqn{\hat{L}} is estimated by the sample maximum obtained from the MCMC realizations.
# #'
# #' @examples
# #' # illustration only - more iterations needed for convergence
# #' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
# #'                              n.chains = 1)
# #' AIC(fit.vmsin.20)
# #' BIC(fit.vmsin.20)
# #'
# #' @export
# #'
# #' AIC.angmcmc <- function(object, ..., k = 2)
# #' {
# #'   if (!is.angmcmc(object))
# #'     stop("\'object\' must be an angmcmc object")
# #'
# #'   ell <- list(...)
# #'
# #'   if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
# #'     warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")
# #'
# #'   llik <- logLik.angmcmc(object, ...)
# #'   k*attr(llik, "df") - 2*as.numeric(llik)
# #' }
# #'
# #'
# #'
# #' @rdname AIC.angmcmc
# #' @export
# #'
# #' BIC.angmcmc <- function(object, ...)
# #' {
# #'   AIC.angmcmc(object, ..., k=log(object$n.data))
# #' }


#' Deviance Information Criterion (DIC) for angmcmc objects
#' @inheritParams pointest
#' @param ... additional model specific arguments to be passed to \code{DIC}. For example, \code{int.displ}
#' specifies integer dispacement in wnorm and wnorm2 models. See \link{fit_wnormmix} and
#' \link{fit_wnorm2mix} for more details.
#' @param form form of DIC to use. Available choices are 1 and 2 (default). See details.
#' @return Computes the DIC for a given angmcmc object
#' @details Given a deviance function \eqn{D(\theta) = -2 log(p(y|\theta))}, and an estimate
#' \eqn{\theta* = (\sum \theta_i) / N} of the posterior mean
#' \eqn{E(\theta|y)}, where \eqn{y} denote the data, \eqn{\theta} are the unknown
#' parameters of the model, \eqn{\theta_1, ..., \theta_N} are MCMC samples from the posterior
#' distribution of \eqn{\theta} given \eqn{y} and \eqn{p(y|\theta)} is the likelihood function,
#' the (form 1 of) Deviance Infomation Criterion (DIC) is defined as
#' \deqn{DIC = 2 ( (\sum_{s=1}^N D(\theta_s)) / N - D(\theta*) )}
#' The second form for DIC is given by
#' \deqn{DIC = D(\theta*) - 4 \hat{var} \log p(y|\theta_s)}
#' where for \eqn{i = 1, ..., n}, \eqn{\hat{var} \log p(y|\theta)} denotes the estimated variance
#' of the log likelihood based on the realizations \eqn{\theta_1, ..., \theta_N}.
#'
#' Like AIC and BIC, DIC is an asymptotic approximation for large samples, and
#' is only valid when the posterior distribution is approximately normal.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' DIC(fit.vmsin.20)
#'
#' @export

DIC <- function(object, form = 2, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  ell <- list(...)

  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if(!form %in% 1:2)
    stop("form must be either 1 or 2")

  final_iter <- object$final_iter

  all_llik <- object$llik[final_iter, ]
  llik_estpar <-
    as.numeric(logLik.angmcmc(object, fn = mean, method = 1))
  D_all <- -2*all_llik
  Dbar <- mean(D_all)
  Dhat <- -2*llik_estpar

  if (form == 1) {
    pD <- Dbar - Dhat
  } else {
    pD <- 0.5*var(D_all)
  }

  c("pD" = pD, "DIC" = pD + Dbar)
}



#' Watanabe-Akaike Information Criterion (WAIC) for angmcmc objects
#'
#' @param x angmcmc object.
#' @param ... additional model specific arguments to be passed to \link[loo]{waic} from loo. For example, \code{int.displ}
#' specifies integer displacement in wnorm and wnorm2 models. See \link{fit_wnormmix} and
#' \link{fit_wnorm2mix} for more details.
#' @return Computes the WAIC for a given angmcmc object.
#' @details
#' Given a deviance function \eqn{D(\eta) = -2 \log(p(y|\eta))}, and an estimate
#' \eqn{\eta* = (\sum \eta_i) / n} of the posterior mean
#' \eqn{E(\eta|y)}, where \eqn{y = (y_1, ..., y_n)} denote the data, \eqn{\eta} is the unknown
#' parameter vector of the model, \eqn{\eta_1, ..., \eta_N} are MCMC samples from the posterior
#' distribution of \eqn{\eta} given \eqn{y} and \eqn{p(y|\eta)} is the likelihood function,
#' the Watanabe-Akaike Information Criterion (WAIC) is defined as
#' \deqn{WAIC = LPPD - p_W}
#' where
#' \deqn{LPPD  = \sum_{i=1}^n \log (N^{-1} \sum_{s=1}^N p(y_i|\eta_s) )}
#' and (form 1 of)
#' \deqn{p_W =  2 \sum_{i=1}^n [ \log (N^{-1} \sum_{s=1}^N p(y_i|\eta_s) ) - N^{-1} \sum_{s=1}^N \log \:p(y_i|\eta_s) ].}
#' An alternative form (form 2) for \eqn{p_W} is given by
#' \deqn{p_W = \sum_{i=1}^n \hat{var} \log p(y_i|\eta)}
#' where for \eqn{i = 1, ..., n}, \eqn{\hat{var} \log p(y_i|\eta)} denotes the estimated variance
#' of \eqn{\log p(y_i|\eta)} based on the realizations \eqn{\eta_1, ..., \eta_N}.
#'
#' Note that waic.angmcmc calls \link[loo]{waic} for computation. If the likelihood contribution of each data
#' point for each MCMC iteration is available in \code{object} (can be returned by setting \code{return_llik_contri = TRUE})
#' during \link{fit_angmix} call), \code{waic.array} is used; otherwise \code{waic.function} is
#' called. Computation is much faster if the likelihood contributions are available - however, they are very
#' memory intensive, and by default not returned in \link{fit_angmix}.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1, return_llik_contri = TRUE)
#' library(loo)
#' waic(fit.vmsin.20)
#'
#' @importFrom loo waic loo
#'
#' @export


waic.angmcmc <- function(x, ...)
{
  object <- x

  if (!is.angmcmc(object))
    stop("\'x\' must be an angmcmc object")

  ell <- list(...)

  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  final_iter <- object$final_iter


  if (object$return_llik_contri) {
    all_llik_contri <- object$llik.contri
    llik_for_loo <- aperm(all_llik_contri[, final_iter, , drop = FALSE], c(2, 3, 1))
    waic <- loo::waic.array(llik_for_loo,  ...)
  }

  else {
    samples <- object$par.value[, , final_iter, , drop=FALSE]
    chain_id <- unlist(lapply(1:object$n.chains, function(j) rep(j, length(final_iter))))

    dim_samples <- dim(samples)
    samples_long <- array(c(samples),
                          dim = c(dim_samples[1:2], dim_samples[3]*dim_samples[4]))
    dimnames(samples_long)[[1]] <- object$par.name

    calc_llik_contri <- function(data_i, draws) {
      N <- dim(draws)[3]
      out <- rep(NA,  N)
      for (j in 1:N) {
        input_pars <- list_by_row(draws[, , j])
        input_pars$log <- TRUE
        input_pars$x <- c(data_i)
        out[j] <- signif(
          do.call(paste0("d", object$model, "mix"), input_pars),
          8
        )
      }
      out
    }

    data <- as.matrix(object$data)

    waic <- loo::waic.function(calc_llik_contri,
                               data = data, draws = samples_long,
                               ...)


  }

  waic
}


#' Leave-one-out cross-validation (LOO) for angmcmc objects
#' @inheritParams waic.angmcmc
#'
#' @examples
#' \donttest{
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1, return_llik_contri = TRUE)
#' library(loo)
#' loo(fit.vmsin.20)
#' }
#' @export
#'
#' @details
#' Note that loo.angmcmc calls \link[loo]{loo} for computation. If the likelihood contribution of each data
#' point for each MCMC iteration is available in \code{object} (can be returned by setting \code{return_llik_contri = TRUE})
#' during \link{fit_angmix} call), \code{loo.array} is used; otherwise \code{loo.function} is
#' called. Computation is much faster if the likelihood contributions are available - however, they are very
#' memory intensive, and by default not returned in \link{fit_angmix}.
#'
#'

loo.angmcmc <- function(x, ...)
{
  object <- x


  if (!is.angmcmc(object))
    stop("\'x\' must be an angmcmc object")

  ell <- list(...)

  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  final_iter <- object$final_iter

  if (object$return_llik_contri) {
    all_llik_contri <- object$llik.contri
    llik_for_loo <- aperm(all_llik_contri[, final_iter, , drop = FALSE],
                          c(2, 3, 1))
    r_eff <- loo::relative_eff(exp(llik_for_loo), ...)
    looic <- loo::loo.array(llik_for_loo, r_eff = r_eff,
                            ...)
  }

  else {
    samples <- object$par.value[, , final_iter, , drop=FALSE]
    chain_id <- unlist(lapply(1:object$n.chains, function(j) rep(j, length(final_iter))))

    dim_samples <- dim(samples)
    samples_long <- array(c(samples),
                          dim = c(dim_samples[1:2], dim_samples[3]*dim_samples[4]))
    dimnames(samples_long)[[1]] <- object$par.name

    calc_llik_contri <- function(data_i, draws) {
      N <- dim(draws)[3]
      out <- rep(NA,  N)
      for (j in 1:N) {
        input_pars <- list_by_row(draws[, , j])
        input_pars$log <- TRUE
        input_pars$x <- c(data_i)
        out[j] <- signif(
          do.call(paste0("d", object$model, "mix"),
                  input_pars),
          8
        )
      }
      out
    }

    data <- as.matrix(object$data)


    r_eff <- loo::relative_eff(calc_llik_contri,
                               chain_id = chain_id,
                               data=data,
                               draws = samples_long)
    looic <- loo::loo.function(calc_llik_contri, r_eff = r_eff,
                               data = data, draws = samples_long,
                               ...)
  }

  looic
}



#' Density and random deviates from an angmcmc object
#' @inheritParams pointest
#' @inheritParams dvmsin
#' @param type Method of estimating density/generating random deviates. Possible choices are
#' \code{"post-pred"} and \code{"point-est"}. See details. Defaults to \code{"point-est"}.
#' @param object angular MCMC object. The dimension of the model must match with \code{x}.
#' @param x vector, if univariate or a two column matrix, if bivariate, with each row a 2-D vector, (can
#' also be a data frame of similar dimensions) of points where the
#' densities are to be computed.
#' @param n number of observations to be generated.
#' @return \code{d_fitted} gives a vector the densities computed at the given points  and \code{r_fitted}
#' creates a vector (if univariate) or a matrix (if bivariate) with each row being a 2-D point, of random deviates.
#'
#' @details
#'
#' If \code{type = 'point-est'}, density is evaluated/random samples are generated at a point estimate of
#' the parameter values.  To estimate the mixture density, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples (using the function \link{pointest}), yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}.
#' Then the mixture density \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by
#' \eqn{f(x|\hat{\eta})}. The random deviates are generated from the estimated mixture density \eqn{f(x|\hat{\eta})}.
#'
#' If \code{type == 'post-pred'}, posterior predictive samples and densities are returned. That
#' is, the average density \eqn{S^{-1} \sum_{s = 1}^S f(x | \eta_s)} is returned in \code{d_fitted},
#' where \eqn{\eta_1, \dots, \eta_S} is the set posterior MCMC samples obtained from \code{object}. In
#' \code{r_fitted}, first a random sub-sample \eqn{\eta_{(1)}, \dots, \eta_{(n)}} of size \code{n} from the
#' set of posterior samples \eqn{\eta_1, \dots, \eta_S} is drawn (with replacement if \code{n} > S). Then
#' the i-th posterior predictive data point is generated from the mixture density
#' \eqn{f(x|\eta_{(i)})} for i = 1,..., n.
#'
#' @examples
#' set.seed(1)
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' d_fitted(c(0,0), fit.vmsin.20, type = "post-pred")
#' d_fitted(c(0,0), fit.vmsin.20, type = "point-est")
#'
#' r_fitted(10, fit.vmsin.20, type = "post-pred")
#' r_fitted(10, fit.vmsin.20, type = "point-est")
#' @export


d_fitted <- function(x, object, type = "point-est", fn = mean, log=FALSE,
                     chain.no, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  ell <- list(...)
  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if (is.data.frame(x))
    x <- as.matrix(x)


  if(object$type == "bi") {
    if((length(dim(x)) < 2 && length(x) != 2) || (length(dim(x)) == 2 && tail(dim(x), 1) != 2)
       || (length(dim(x)) > 2)) stop("x must either be a bivariate vector or a two-column matrix")
  }

  if (!type %in% c("post-pred", "point-est")) {
    stop ("\'type\' must either be \'post-pred\' or \'point-est\'")
  }

  inargs <- list(x = x)
  if (object$model %in% c("wnorm", "wnorm2"))
    inargs$int.displ <- object$int.displ
  else if (object$model == "vmcos")
    inargs$qrnd_grid <- object$qrnd_grid

  if (type == "point-est") {
    est <- pointest(object, fn = fn)
    inargs <- c(list_by_row(est), inargs)
    inargs$log <- log
    out <- signif(suppressWarnings(do.call(paste0("d", object$model, "mix"), inargs)), 8)
  } else {
    samp <- extractsamples(object, chain.no = chain.no, drop = FALSE)
    dim_samp <- dim(samp)
    nsamp <- dim_samp[3]*dim_samp[4]
    samp_coll <- samp
    dim(samp_coll) <- c(dim_samp[1:2], nsamp)
    dimnames(samp_coll)[1:2] <- dimnames(samp)[1:2]
    inargs$log <- FALSE
    out_list <- lapply(seq_len(nsamp),
                       function(j) {
                         inargs1 <- c(list_by_row(as.matrix(samp_coll[, , j])), inargs)
                         signif(
                           c(suppressWarnings(do.call(paste0("d",
                                                             object$model,
                                                             "mix"), inargs1))),
                           8
                         )
                       })
    out <- Reduce('+', out_list)/nsamp

    if (log) out <- log(out)
  }

  out
}



#' @rdname d_fitted
#' @export
r_fitted <- function(n=1, object, type =  "point-est", fn = mean,
                     chain.no,  ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  ell <- list(...)
  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  if (!type %in% c("post-pred", "point-est")) {
    stop ("\'type\' must either be \'post-pred\' or \'point-est\'")
  }


  inargs <- list()

  if (object$model == "vmcos") {
    inargs$qrnd_grid <- object$qrnd_grid
  }


  if (object$model %in% c("vmcos", "vmsin")) {
    inargs$method <- "naive"
  }


  if (type == "point-est") {
    est <- pointest(object, fn = fn, chain.no = chain.no)
    inargs <- c(list_by_row(as.matrix(est)), inargs)
    inargs$n <- n
    out <- suppressWarnings(do.call(paste0("r", object$model, "mix"), inargs))
  } else {

    samp <- extractsamples(object, chain.no = chain.no, drop = FALSE)
    dim_samp <- dim(samp)
    nsamp <- dim_samp[3]*dim_samp[4]
    samp_coll <- samp
    dim(samp_coll) <- c(dim_samp[1:2], nsamp)
    dimnames(samp_coll)[1:2] <- dimnames(samp)[1:2]
    ids <- sample(seq_len(dim_samp[3]*dim_samp[4]), size = n,
                  replace = n > nsamp)

    if (object$type == "bi") {
      out_row <- numeric(2)
    } else {
      out_row <- numeric(1)
    }

    out <- vapply(ids,
                  function(j) {
                    inargs1 <- c(list_by_row(as.matrix(samp_coll[, , j])), inargs)
                    inargs1$n <- 1
                    c(suppressWarnings(do.call(paste0("r",
                                                      object$model,
                                                      "mix"), inargs1)))
                  },
                  out_row)

    if (object$type == "bi") out <- t(out)
  }

  out
}


#' Extract Log-Likelihood from angmcmc objects
#' @inheritParams pointest
#' @inheritParams stats::logLik
#' @param method interger specifying method of estimating the log likelihood. Must be 1 or 2. Defaults to 1. See details.
#' @param fn function to evaluate on the iteration-wise log-likelihood values obtained during MCMC run if \code{method = 1}; or,
#' if \code{method = 2}, function to evaluate on the MCMC samples for parameter estimation (passed to \link{pointest}).
#' Defaults to \code{max} if \code{method = 1} and \code{mean} if \code{method = 2}.
#'
#' @details
#'
#' There are two ways to estimate the log likelihood from the model. If \code{method = 1},
#' then log likelihood is estimated by applying \code{fn} (defaults to max, if method = 1)
#' direclty on the log likelihood values from observed during the MCMC run.
#' On the other hand, if \code{method == 2}, then  parameter estimates
#' are first computed using \code{pointest} with \code{fn}
#' (defaults to "MODE", if \code{method == 2}) applied on the MCMC samples,
#' and then then log likelihood is evaluated at the parameter estimates.
#'
#'
#' The degrees of the likelihood function is the total number of free parameters estimated in the mixture models,
#' which is equal to \eqn{6K - 1} for bivariate models (vmsin, vmcos and wnorm2), or \eqn{3K - 1} for univariate
#' models (vm and wnorm), where \eqn{K} denotes the number of components in the mixture model.
#'
#' @return Returns an object of class \link{logLik}. This is a number (the estimated log likelihood) with attributes "df"
#' (degrees of freedom) and "nobs" (number of observations).
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' logLik(fit.vmsin.20)
#' @export


logLik.angmcmc <- function(object, method = 1, fn, ...)
{
  if (!is.angmcmc(object))
    stop("\'object\' must be an angmcmc object")

  if (!method %in% c(1, 2))
    stop("method must be either 1 or 2")

  if (missing(fn)) {
    if (method == 1) fn <- "max"
    else fn <- "MODE"
  }

  do_MAP <- FALSE

  if (is.character(fn)) {
    if (fn %in% c("MODE", "MAP")) {
      if (method == 1)
        stop("fn can be \'MODE\' or \'MAP\' only when method = 2")
      else do_MAP <- TRUE
    }
  }

  if (!do_MAP) fn <- match.fun(fn)

  ell <- list(...)
  if (any(!is.null(ell$burnin),  !is.null(ell$thin)))
    warning("Use of burnin and thin are obsolete. Specify \'burnin.prop\' and \'thin\' during original MCMC run, or use \'add_burnin_thin\'.")

  final_iter <- object$final_iter

  if (method == 2) {
    all_den <- d_fitted(object$data, object, fn = fn, log=TRUE)
    llik <- sum(all_den)
  } else {
    llik <- fn(object$llik[final_iter, ])
  }


  if(object$type == "uni") {
    df <- 3*object$ncomp - 1
  } else if (object$cov.restrict == "ZERO") {
    df <- 5*object$ncomp - 1 #the cov terms aren't present
  } else {
    df <- 6*object$ncomp - 1
  }

  n_data <- object$n.data
  result <- llik
  attributes(result) <- list("df" = df, "nobs" = n_data)
  class(result) <- "logLik"
  result

}


#' Log Marginal Likelihood via Bridge Sampling for angmcmc objects
#' @param samples angmcmc object
#' @param ... additional argument passed to \link[bridgesampling]{bridge_sampler}. Note that default for
#' the argument \code{method} is \code{"warp3"}, (instead of \code{"normal"} as used in \code{bridgesampling} package)
#' to account for multi-modality of the posterior density.
#' @param ave_over_chains logical. Separately call \link[bridgesampling]{bridge_sampler} on
#' each chain in the angmcmc object and then take the average? Defaults to \code{TRUE}.
#' See details.
#'
#' @details
#' Marginal likelihood is calculated by first converting the \code{angmcmc} object \code{samples} to an
#' \code{mcmc.list} object, and then by passing the resulting \code{mcmc.list} object to \link[bridgesampling]{bridge_sampler}.
#' If variablity across multiple chains (if any) are very different,
#' then calling \link[bridgesampling]{bridge_sampler} separately for each chain
#' usually provides more stable results; the final log ML is computed by averaging over
#' chain specific MLs.
#'
#' @importFrom bridgesampling bridge_sampler
#'
#' @examples
#' \donttest{
#' library(future)
#' library(parallel)
#' # plan(multisession, gc = TRUE) # parallelize chains
#'
#' set.seed(100)
#' MC.fit <- fit_angmix("vmsin", tim8, ncomp=3, n.iter=5000,
#'                      n.chains = 3)
#'
#'
#' library(bridgesampling)
#' bridge_sampler(MC.fit)
#' }
#'
#' @export

bridge_sampler.angmcmc <- function(samples, ..., ave_over_chains = TRUE)
{

  object <- samples

  if (!is.angmcmc(object))
    stop("\'samples\' must be an angmcmc object")


  ell <- list(...)

  if (is.null(ell$method)) {
    ell$method <- "warp3"
  }

  n.chains <- object$n.chains

  if (object$model == "wnorm2") {
    kappa_tmp <- object$par.value
    sig_tmp <-
      apply(kappa_tmp, 2:4,
            function(x) c(x[1],
                          unname(kappas2sigmas_wnorm2(x[2], x[3], x[4])),
                          x[5:6]))
    names(sig_tmp) <- names(kappa_tmp)
    object$par.value <- sig_tmp
  }


  object_mcmc <- coda::as.mcmc.list(object)

  object_mcmc_matrix <- as.matrix(object_mcmc)

  gam.loc <- object$gam.loc
  gam.rate <- 1/object$gam.scale


  norm.var <- object$norm.var

  if (object$type == "bi" & length(norm.var) == 1) norm.var <- rep(norm.var, 3)


  if (object$ncomp == 1) {

    calc_lpd <- function(par_vec, data) {
      # reparametrize if wnorm2
      if (object$model == "wnorm2") {
        par_vec <- c(sigmas2kappas_wnorm2(par_vec[1], par_vec[2], par_vec[3]),
                     par_vec[4:5])
      }
      allpar_mat <- as.matrix(par_vec)

      rownames(allpar_mat) <- object$par.name[-1]
      inargs <- list_by_row(allpar_mat)
      inargs$x <- data
      if (grepl("wnorm", object$model))
        inargs$int.displ <- object$int.displ
      if (object$model == "vmcos")
        inargs$qrnd_grid <- object$qrnd_grid

      inargs$log <- TRUE
      llik <- sum(signif(do.call(paste0("d", object$model), inargs), 8))

      if (object$type == "bi") {
        lprior <- sum(-0.5*c((log(allpar_mat[1, ]))^2/norm.var[1],
                             (log(allpar_mat[2, ]))^2/norm.var[2],
                             ((allpar_mat[3, ]))^2/norm.var[3]))
      } else {
        lprior <- sum(-0.5*log(allpar_mat[1, ])^2/norm.var)
      }

      unname(llik+lprior)
    }


    # lower and upper limits for reparameterized wnorm2

    if (object$model == "wnorm2") {
      lower <- object$modelpar.lower
      upper <- object$modelpar.upper
      upper[1:2, ] <- 1e5
      lower[3, ] <- -(1-1e-7)
      upper[3, ] <- (1-1e-7)
      lower <- c(lower)
      upper <- c(upper)
    }
    else {
      lower <- c(object$modelpar.lower)
      upper <- c(object$modelpar.upper)
    }

    names(lower) <- names(upper) <- colnames(object_mcmc_matrix)

    if(n.chains == 1 || !ave_over_chains) {
      tmp_mat <- do.call(rbind, lapply(1:n.chains,
                                       function(j)
                                         as.matrix(object_mcmc[[j]])))

      do.call(bridgesampling::bridge_sampler,
              c(list(samples=tmp_mat,
                     log_posterior=calc_lpd,
                     data=object$data,
                     lb=lower, ub=upper),
                ell))
    }
    else {
      all_bridge_samp <- vector("list", n.chains)
      for (j in 1:n.chains) {
        # tmp <- coda::as.mcmc.list(select_chains(object, j))
        tmp_mat <- as.matrix(object_mcmc[[j]])
        if (!any(ell$silent)) {
          cat(paste0("\napplying \'bridge_sampler\' on chain ", j, "...\n\n"))
        }

        all_bridge_samp[[j]] <-
          do.call(bridgesampling::bridge_sampler,
                  c(list(samples=tmp_mat,
                         log_posterior=calc_lpd,
                         data=object$data,
                         lb=lower, ub=upper),
                    ell))
      }

      final_res <- all_bridge_samp[[1]]
      final_res$logml <- mean(vapply(all_bridge_samp,
                                     function(x) x$logml, 0))
      final_res$niter <- max(vapply(all_bridge_samp,
                                    function(x) x$niter, 0))
      final_res$method <- all_bridge_samp[[1]]$method
      final_res$q11 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q11))
      final_res$q12 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q12))
      final_res$q21 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q21))
      final_res$q22 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q22))

      final_res$all <- all_bridge_samp
      final_res
    }
  }


  else {

    calc_lpd <- function(par_vec, data) {
      allpar_mat <- matrix(par_vec,
                           nrow = length(object$par.name))
      # reparametrize if wnorm2
      if (object$model == "wnorm2") {
        allpar_mat <- apply(allpar_mat, 2,
                            function(x) c(x[1],
                                          unname(sigmas2kappas_wnorm2(x[2], x[3], x[4])),
                                          x[5:6]))
      }
      allpar_mat[1, ] <- allpar_mat[1, ]/sum(allpar_mat[1, ])
      rownames(allpar_mat) <- object$par.name
      inargs <- list_by_row(allpar_mat)
      inargs$x <- data
      if (grepl("wnorm", object$model))
        inargs$int.displ <- object$int.displ
      if (object$model == "vmcos")
        inargs$qrnd_grid <- object$qrnd_grid

      inargs$log <- TRUE
      llik <- suppressWarnings(
        sum(
          signif(do.call(paste0("d", object$model, "mix"), inargs)), 8
        )
      )


      par_mat <- allpar_mat[-1, ]
      pmix_vec <- allpar_mat[1, ]


      if (object$type == "bi") {
        lprior <- sum(-0.5*c((log(allpar_mat[1, ]))^2/norm.var[1],
                             (log(allpar_mat[2, ]))^2/norm.var[2],
                             ((allpar_mat[3, ]))^2/norm.var[3])) +
          sum(object$pmix.alpha*log(pmix_vec))
      } else {
        lprior <- sum(-0.5*log(allpar_mat[1, ])^2/norm.var) +
          sum(object$pmix.alpha*log(pmix_vec))
      }

      # if (object$type == "bi") {
      #   lprior <- sum((gam.loc - 1)*log(par_mat[1:2, ]) -
      #                   gam.rate*par_mat[1:2, ]) -
      #     0.5*sum(par_mat[3, ]^2)/object$norm.var +
      #     sum(object$pmix.alpha*log(pmix_vec))
      # } else {
      #   lprior <- sum((gam.loc - 1)*log(allpar_mat[1, ]) -
      #                   gam.rate*allpar_mat[1, ]) +
      #     sum(object$pmix.alpha*log(pmix_vec))
      # }

      unname(llik+lprior)
    }


    # lower and upper limits for reparameterized wnorm2

    if (object$model == "wnorm2") {
      lower <- rbind(0, object$modelpar.lower)
      upper <- rbind(Inf, object$modelpar.upper)
      upper[2:3, ] <- 1e5
      lower[4, ] <- -(1-1e-7)
      upper[4, ] <- (1-1e-7)
      lower <- c(lower)
      upper <- c(upper)
    }

    else {
      lower <- c(rbind(0, object$modelpar.lower))
      upper <- c(rbind(Inf, object$modelpar.upper))
    }

    names(lower) <- names(upper) <-
      colnames(object_mcmc_matrix)
    if (n.chains == 1 || !ave_over_chains) {
      tmp_mat <- do.call(rbind, lapply(1:n.chains,
                                       function(j)
                                         as.matrix(object_mcmc[[j]])))
      do.call(bridgesampling::bridge_sampler,
              c(list(samples=tmp_mat,
                     log_posterior=calc_lpd,
                     data=object$data,
                     lb=lower, ub=upper),
                ell))

    }
    else {
      all_bridge_samp <- vector("list", n.chains)
      for (j in 1:n.chains) {
        # tmp <- coda::as.mcmc.list(select_chains(object, j))
        tmp_mat <- as.matrix(object_mcmc[[j]])
        if (!any(ell$silent)) {
          cat(paste0("\napplying \'bridge_sampler\' on chain ", j, "...\n\n"))
        }

        all_bridge_samp[[j]] <-
          do.call(bridgesampling::bridge_sampler,
                  c(list(samples=tmp_mat,
                         log_posterior=calc_lpd,
                         data=object$data,
                         lb=lower, ub=upper),
                    ell))
      }

      final_res <- all_bridge_samp[[1]]
      final_res$logml <- mean(vapply(all_bridge_samp,
                                     function(x) x$logml, 0))
      final_res$niter <- max(vapply(all_bridge_samp,
                                    function(x) x$niter, 0))
      final_res$method <- all_bridge_samp[[1]]$method
      final_res$q11 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q11))
      final_res$q12 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q12))
      final_res$q21 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q21))
      final_res$q22 <- unlist(lapply(all_bridge_samp,
                                     function(x) x$q22))

      final_res$all <- all_bridge_samp
      final_res
    }
  }

}





#' Finding latent allocation (component indicators) from an angmcmc object
#' @inheritParams pointest
#' @param ... passed to \link{pointest} to estimate parameters.
#'
#' @details
#' In order to find the latent component indicators, estimates
#' of mixing proportions and model parameters are first computed via
#' pointest. Then, a data point is assigned label j, if the j-th
#' component gives highest density for that point.
#'
#' @return
#' Returns a vector of length n, where n is the length (if univariate) or
#' number of rows (if bivariate) of the data used in original fit.
#' i-th entry of the output vector provides component label for the i-th data point.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # now find latent allocation
#' latent_allocation(fit.vmsin.20)
#' @export
latent_allocation <- function(object, ...)
{

  if (!is.angmcmc(object))
    stop("object must be an \'angmcmc\' object")

  ell <- list(...)
  if (any(!is.null(c(ell$comp.label, ell$par.name)))) {
    warning("\'par.name\' and \'comp.label\' are not used")
    ell$comp.label <- NULL
    ell$par.name <- NULL
  }


  if (object$ncomp > 1) {
    parest <- pointest(object, ...)

    pi.mix <- parest[1, ]
    par.mat <- parest[-1, ]

    if (object$model == "vmsin") {
      l_c_vmsin <- log_const_vmsin_all(par.mat)
      mem_p_model <- mem_p_sin(object$data, par.mat, pi.mix, l_c_vmsin)
    }
    else if (object$model == "vmcos") {
      l_c_vmcos <- log_const_vmcos_all(par.mat, object$qrnd_grid)
      mem_p_model <- mem_p_cos(object$data, par.mat, pi.mix, l_c_vmcos)
    }
    else if (object$model == "wnorm2") {
      l_c_wnorm2 <- log_const_wnorm2_all(par.mat)
      mem_p_model <- mem_p_wnorm2(object$data, par.mat, pi.mix, l_c_wnorm2, object$omega.2pi)
    }
    else if (object$model == "vm") {
      l_c_vm <- log_const_univm_all(par.mat)
      mem_p_model <- mem_p_univm(object$data, par.mat, pi.mix, l_c_vm)
    }
    else {
      l_c_wnorm <- log_const_uniwnorm_all(par.mat)
      mem_p_model <- mem_p_uniwnorm(object$data, par.mat, pi.mix, l_c_wnorm, object$omega.2pi)
    }

    labs <- apply(mem_p_model, 1, which.max_entry1)
  }

  else {
    warning("\'object\' has only one component. A vector of all 1\'s will is returned")
    labs <- rep(1, object$n.data)

  }

  labs
}
