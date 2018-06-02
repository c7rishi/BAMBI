#' Fitting Bivariate and univariate angular mixture models
#'
#' @importFrom gtools rdirichlet
#'
#' @param model angular model whose mixtures are to be fitted. Available choices are \code{"vmsin", "vmcos"} and \code{"wnorm2"} for
#' bivariate data, and \code{"vm"} and \code{"wnorm"} for univariate data.
#' @param data data matrix (if bivarate, in which case it must have two columns) or vector. If outside, the values
#' are transformed into the scale \eqn{[0, 2\pi)}.
#' @param ncomp number of components in the mixture model. Must be a positive integer. vector values are not allowed.
#' If \code{comp == 1}, a single component model is fitted.
#' @param start_par list with elements \code{pmix} (ignored if \code{comp == 1}), together with \code{kappa1, kappa2, mu1} and \code{mu2},
#' for bivariate models, and \code{kappa} and \code{mu} for univariate models,
#' all being vectors of length same as \code{ncomp}.
#' These provides the starting values for the Markov chain; with \eqn{j}-th component of each vector corresponding to the \eqn{j}-th
#' component of the mixture distribution. If missing, the data is first clustered into \code{ncomp} groups either via k-means (after
#' projecting onto a unit sphere), or randomly, depending on \code{rand_start},  and then moment estimators for components are used as
#' the starting points. Note that a very wrong starting point can potentially lead the chain to get stuck at a wrong solution for thousands
#' of iterations. As such, we recommend using the default option, which is k-means followed by moment estimation.
#' @param cov.restrict Should there be any restriction on the covariance parameter for a bivariate model. Available choices are
#' \code{"POSITIVE", "NEGATIVE", "ZERO"} and "NONE". Note that \code{"ZERO"} fits a mixture with product components. Defaults to
#' \code{"NONE"}.
#' @param unimodal.component logical. Should each component in the mixture model be unimodal? Only used if \code{model} is either \code{"vmsin"}
#' or \code{"vmcos"}. Defaults to FALSE.
#' @param n.chains number of chains to run. Must be a positive integer.
#' @param chains_parallel logical. Should the chains be run in parallel? Defaluts to TRUE, and ignored if \code{n.chains} = 1.
#' Note that parallelization is implemented via \link{future_lapply} from package \code{future.apply} which
#' uses futures for this purpose, and thus provides a convenient way of parallelization across various OSs and computing environments.
#' However, a proper \link{plan} must be set for the parallization before running the chain. Otherwise the chains will run sequentially.
#' @param method MCMC strategy to be used for the model paramters:  \code{"hmc"} or \code{"rwmh"}.
#' @param perm_sampling logical. Should the permutation sampling algorithm of Fruhwirth-Schnatter (2001) be used?
#' If TRUE, at every iteration after burnin, once model parameters and mixing proportions are sampled,
#' a random permutation of 1, ..., ncomp is considered, and components are relabelled according
#' to this random permutation. This forced random label switchings may imporve the mixing rate of the chage. However, (automated) tuning
#' is very difficult with such a scheme, as there is no simple way of keeping track of the "original" component labels. This creates problem
#' with computing standard deviations of the generated model parameters, thus making the
#' scaling step used in tuning for \code{epsilon} or \code{paramscale} problematic as well. As such, \code{perm_sampling} is always turned
#' off during burn-in (even if \code{autotune = FALSE}), and turned on thereafter, if \code{TRUE}.
#' Defaults to and is set to \code{FALSE} if \code{pmix.alpha} is a vector with non-identical values, which indicates a non-exchangeable
#' Dirichlet prior for the mixing proportion.
#' @param int.displ absolute integer displacement for each coordinate for \code{wnorm} and \code{wnorm2} models (ignored otherwise). Default is 3.
#' Allowed minimum and maximum are 1 and 5 respectively.
#' @param epsilon,L  tuning parameters for HMC; ignored if \code{method = "rwmh"}. \code{epsilon} (step-size) is a single number,
#' or a vector of size \code{2*ncomp} for univariate models and \code{5*ncomp} for bivariate models. Note that the "mass matrix"
#' in HMC is assumed to be identity. As such, \code{epsilon}'s corresponding to different model parameters need to be in proper scale for
#' optimal acceptance rate. Can be autotuned during burnin. See \code{autotune}.
#' \code{L} (leapfrog steps) is a positive integer or a vector of positive integers of length \code{n.chains}.
#' If multiple chains are used, we suggest same \code{L} values acorss different chains to make the chains as homogenous as possible.
#'
#' @param epsilon.random logical. Should \code{epsilon*delta}, where \code{delta} is a random
#' number between \code{(1-epsilon.incr, 1+epsilon.incr)} be used instead of \code{epsilon} at each iteration?
#' Ignored if \code{method = "rwmh"}.
#' @param L.random logical. Should a random integer between \code{L.orig/exp(L.incr)} and \code{L.orig*exp(L.incr)}be used instead as \code{L}
#' at each iteration? Ignored if \code{method = "rwmh"}. Defaults to \code{TRUE}.
#' @param L.incr amount of randomness incorporated in L if \code{L.random = TRUE}.
#' @param epsilon.incr amount of randomness incorporated in \code{epsilon} if \code{epsilon.random = TRUE}.
#' @param propscale tuning parameters for RWMH; a vector of size 5 (for bivariate models) or 2 (for univariate models) representing
#' the variances for the proposal normal densities
#' for the model parameters. Ignored if \code{method = "hmc"}. Can be autotuned during burnin. See \code{autotune}.
#' @param n.iter number of iterations for the Markov Chain.
#' @param gam.loc,gam.scale location and scale (hyper-) parameters for the gamma prior for \code{kappa1} and \code{kappa2}. See
#' \link{dgamma}. Defaults are \code{gam.loc = 0, gam.scale = 1000} that makes the prior non-informative.
#' @param pmix.alpha concentration parameter(s) for the Dirichlet prior for \code{pmix}. Must either be a positive real number, or a vector
#' with positive entries and of length \code{ncomp}. The default is \eqn{(r+r(r+1)/2)/2+3}, where \eqn{r} is 1 or 2 according as whether
#' the model is univariate or bivariate. Note that it is recommended to use larger \code{alpha} values to ensure the a good posterior behavior,
#' especially when \link{fit_incremental_angmix} is used for model selection, which handles overfitting in "let two component-specific parameters be
#  identical", uses total number of components in the fitted model as the estimator for true component
#' size, and then penalizes for model complexity. See Fruhwirth-Schnatter (2011) for more details on this.
#' @param norm.var variance (hyper-) parameter in the normal prior for \code{kappa3}. (Prior mean is zero).
#' Default is 1000 that makes the prior non-informative. Ignored if the model is univariate or if cov.restrict = "ZERO".
#' @param burnin.prop proportion of iterations to used for burnin. Must be a be a number in [0, 1].
#' Default is 0.5.
#' @param thin thining size to be used. Must be a positive integer. If \code{thin = } n, then every nth iteration is reatained
#' in the final MCMC sample.
#' @param  autotune logical. Should the Markov chain auto-tune the parameter \code{epsilon} (in HMC) or
#' \code{propscale} (in RWMH) during burn-in?  Set to \code{TRUE} by default. An adaptive tuning strategy is implemented.
#' Here, at every 10th iteration during in burn-in, the acceptance ratio in the last \code{tune_ave_size}
#' iterations is calculated. Then the tuning parameter is decreased  (increased) by a factor of
#' \code{1-tune.incr} (\code{1+tune.incr}) if the calculated acceptance rate
#' falls below (above) \code{accpt.prob.lower} (\code{accpt.prob.upper}). In addditon, when \code{iter} is a multiple of
#' \code{tune_ave_size}, \code{epsilon} for each model parameter is rescaled via the standard deviation of
#' the corresponding parameter over the past \code{tune_ave_size} iterations.
#' @param tune.prop proportion of *\code{burnin}* used to tune the parameters (\code{epsilon} in HMC and
#' \code{propscale} in RWMH). Must be a number between 0 and 1; defaults to 1.  Ignored if \code{autotune == FALSE}.
#' @param show.progress logical. Should a progress bar be included?
#' @param accpt.prob.lower,accpt.prob.upper lower and upper limits of acceptance ratio to be maintained while tuning
#' during burn-in. Must be numbers between 0 and 1, which \code{accpt.prob.lower < accpt.prob.upper}. See \code{autotune}. Default to (0.6, 0,9) for HMC and  (0.3, 0.5) for RWMH.
#' Ignored if \code{autotune = FALSE}.
#' @param tune.incr how much should the tuning parameter be increased or decreased at each step while tuning during burn-in?
#' Must be a number between 0 and 1. See \code{autotune}. Defaults to 0.05. Ignored if \code{autotune = FALSE}.
#' @param rand_start logical. Should a random starting clustering be used? Must be either a scalar, or a vector of length \code{ncomp},
#' one for each chain. Ignored if \code{start_par} is supplied. See \code{start_par} for more details. Defaults to \code{FALSE}.
#' @param tune_ave_size number previous iterations used to compute the acceptance rate while tuning in burn-in. Must be a positive
#' integer. Defaults to 100.
#' @param qrnd_grid,n_qrnd Used only if \code{method="vmcos"}. See \link{dvmcos} for details.
#' @param kappa_upper,kappa_lower upper and lower bounds for the concentration and (absolute) association parameters. Must be a positive integers. Defaults to 150 and 1e-4,
#' and parameter with value above or below these limits rarely make sense in practice.
#' Warning: values much larger or smaller than the default are not recommended as they can cause numerical instability.
#' @param return_llik_contri logical. Should the log likelihood contribution of each data point for each MCMC iteration in each chain be returned? This makes
#' computation of \link{waic.angmcmc} and \link{loo.angmcmc} much faster. *Warning*: Depending on the length of data and \code{n.iter}, this can be
#' very memory intensive. We suggest setting \code{return_llik_contri = TRUE} only if \link{waic.angmcmc} and \link{loo.angmcmc} are aimed for. Defaults to
#' \code{FALSE}.
#' @param return_tune_param logical. Should the values of the tuning parameters used at each iteration in each chain be returned? Defaults to \code{FALSE}.
#' @param ... Unused.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_angmix("vmsin", tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' fit.vmsin.20
#'
#'
#' # Parallelization is implemented via \link{future_lapply} from the
#' # package future.apply. To parallelize, first provide a parallel
#' # plan() for futures. Otherwise the chains will run sequentially.
#' # Note that not all plan() might work on every OS, as they execute
#' # functions defined internally in fit_mixmodel. We suggest
#' # plan(multiprocess).
#'
#' \dontrun{
#' library(future)
#' library(parallel)
#' plan(multiprocess)
#'
#' set.seed(1)
#' MC.fit <- fit_angmix("vmsin", tim8, ncomp=3, n.iter=500,
#'                      n.chains = 3)
#'
#'
#' pointest(MC.fit)
#'
#' MC.fix <- fix_label(MC.fit)
#'
#' contour(MC.fit)
#' contour(MC.fix)
#' lpdtrace(MC.fit)
#' }
#'
#' @references
#' Fruhwirth-Schnatter, S. (2011). Label switching under model uncertainty. Mixtures: Estimation and Application, 213-239.
#'
#' Fruhwirth-Schnatter, S. (2001). Markov chain Monte Carlo estimation of classical and dynamic switching and mixture models. Journal of the American Statistical Association, 96(453), 194-209.
#'
#' @export

fit_angmix <- function(model = "vmsin",
                       data,
                       ncomp,
                       cov.restrict = "NONE",
                       unimodal.component = FALSE,
                       start_par = NULL,
                       rand_start = rep(FALSE, n.chains),
                       method="hmc",
                       perm_sampling = FALSE,
                       n.chains = 3,
                       chains_parallel = TRUE,
                       return_llik_contri = FALSE,
                       int.displ = 3,
                       epsilon = 0.1,
                       L = 10,
                       epsilon.random=TRUE,
                       L.random=FALSE,
                       burnin.prop = 0.5,
                       tune.prop = 1,
                       thin = 1,
                       propscale = 0.05,
                       n.iter = 500,
                       gam.loc = 0.001,
                       gam.scale = 1000,
                       pmix.alpha = NULL,
                       norm.var = 1000,
                       autotune = TRUE,
                       show.progress = TRUE,
                       accpt.prob.upper,
                       accpt.prob.lower,
                       epsilon.incr = 0.05,
                       L.incr = 0.075,
                       tune.incr = 0.05,
                       tune_ave_size = 100,
                       kappa_upper = 150,
                       kappa_lower = 1e-4,
                       return_tune_param = FALSE,
                       qrnd_grid = NULL,
                       n_qrnd = NULL, ...)
{

  # if(is.null(dim(data)) | !(mode(data) %in% c("list", "numeric", "data.frame") && ncol(data) == 2)) stop("non-compatible data")

  if (length(model) > 1) stop("\'model\' must be a scalar")

  if (model %in% c("vmsin", "vmcos", "wnorm2")) {
    type <- "bi"
  } else if (model %in% c("vm", "wnorm")) {
    type <- "uni"
  } else  {
    stop("non-compatible model")
  }

  if (missing(ncomp))
    stop("\'ncomp\' is missing, with no default.")

  if (length(ncomp) > 1)
    stop("\'ncomp\' must be a scalar")


  if (any(length(n.chains) > 1, length(n.chains) < 1, n.chains == 0))
    stop("Invalid n.chains")

  zero_cov <- FALSE

  if (type == "bi") {
    if (!(is.matrix(data) | is.data.frame(data)))
      stop("\'data\' must be a two column matrix for model = \'vmsin\', \'vmcos\' and \'wnorm2\'")

    if (ncol(data) != 2)
      stop("\'data\' must be a two column matrix for model = \'vmsin\', \'vmcos\' and \'wnorm2\'")

    if (length(cov.restrict) > 1 |
        !cov.restrict %in% c("NONE", "POSITIVE", "NEGATIVE", "ZERO"))
      stop("cov.restrict must be one of \"NONE\", \"POSITIVE\", \"NEGATIVE\" and  \"ZERO\"")


    data.rad <- rm_NA_rad(data)
    n.data <- nrow(data.rad)

    npar_1_comp <- 5
    par.names <- c("kappa1", "kappa2", "kappa3", "mu1", "mu2")

    par_lower <- replicate(ncomp, c(kappa_lower, kappa_lower, -kappa_upper, 0, 0))
    par_upper <- replicate(ncomp, c(kappa_upper, kappa_upper,
                                    kappa_upper, 2*pi, 2*pi))

    if (cov.restrict == "POSITIVE") {
      par_lower[3, ] <- 0

    } else if (cov.restrict == "NEGATIVE") {
      par_upper[3, ] <- 0

    } else if (cov.restrict == "ZERO") {
      par_lower[3, ] <- par_upper[3, ] <- 0
      zero_cov <- TRUE
    }

    if (missing(pmix.alpha))
      pmix.alpha <- 5.5

  }

  else {

    if ((is.matrix(data) | is.data.frame(data))) {
      data <- as.vector(as.matrix(data))
    }

    if (!is.numeric(data))
      stop("\'data\' must be a vector for \'model\' = \'vm\' and \'wnorm\'")

    data <- as.numeric(data)
    data.rad <- rm_NA_rad(data)
    n.data <- length(data.rad)

    npar_1_comp <- 2
    par.names <- c("kappa", "mu")

    par_lower <- replicate(ncomp, c(kappa_lower, 0 ))
    par_upper <- replicate(ncomp, c(kappa_upper, 2*pi))

    if (missing(pmix.alpha))
      pmix.alpha <- 4
  }

  if (any(length(perm_sampling) > 1, length(perm_sampling) < 1,
          !is.logical(perm_sampling) ))
    stop("Invalid perm_sampling")


  if (any(length(chains_parallel) > 1, length(chains_parallel) < 1,
          !is.logical(chains_parallel) ))
    stop("Invalid chains_parallel")

  if (any(length(chains_parallel) > 1, length(chains_parallel) < 1,
          !is.logical(chains_parallel) ))
    stop("Invalid chains_parallel")

  if (any(length(autotune) > 1, length(autotune) < 1,
          !is.logical(autotune) ))
    stop("Invalid autotune")

  if (any(length(unimodal.component) > 1, length(unimodal.component) < 1,
          !is.logical(unimodal.component) ))
    stop("Invalid unimodal.component")

  if (any(length(return_llik_contri) > 1, length(return_llik_contri) < 1,
          !is.logical(return_llik_contri) ))
    stop("Invalid return_llik_contri")

  if (any(length(return_tune_param) > 1, length(return_tune_param) < 1,
          !is.logical(return_tune_param) ))
    stop("Invalid return_tune_param")





  if (length(pmix.alpha) == 1)
    pmix.alpha <- rep(pmix.alpha, ncomp)
  else if (length(pmix.alpha) != ncomp)
    stop("length(pmix.alpha) and ncomp differ")

  if (any (diff(pmix.alpha) != 0))
    perm_sampling <- FALSE

  if (n.chains < 1)
    stop("\'n.chains\' must be a positive integer")

  if (n.chains == 1) {
    chains_parallel <- FALSE
  }




  if(length(burnin.prop) != 1 || burnin.prop < 0 || burnin.prop >= 1)
    stop("\"burnin.prop\" must be a number in [0, 1)")

  if(thin < 1)
    stop("\"thin\" must be a positive integer")

  if (tune.incr <= 0 | tune.incr >= 1)
    stop("\'tune.incr\' must be between 0 and 1")

  n.burnin <- ceiling(burnin.prop*n.iter)

  if(length(tune.prop) != 1 ||tune.prop < 0 || tune.prop > 1)
    stop("\"tune.prop\" must be in [0, 1]")

  iter.tune <- ceiling(burnin.prop*tune.prop*n.iter)

  if (iter.tune > n.burnin)
    iter.tune <- n.burnin

  n.iter.final <- n.iter - n.burnin
  burnin_iter <- seq_len(n.burnin)
  thin <- round(thin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter_set <- (seq_len(n.iter))[-burnin_iter][thin_filter]

  gam.rate <- 1/gam.scale

  curr.model <- model

  if (length(rand_start) == 1) {
    rand_start <- rep(rand_start, n.chains)
  }


  if (method == "hmc") {
    if (any(L < 1))
      stop("\'L\' must be a positive integer")

    L <- ceiling(L)

    if (length(L) == 1) {
      L <- rep(L, n.chains)
    } else if (length(L) != n.chains) {
      stop("\'L\' must be a vector of length \'n.chains\'")
    } else {
      rand_start <- rep(rand_start, n.chains)
    }

    if (missing(accpt.prob.upper)) {
      accpt.prob.upper <- 0.9
    }

    if (missing(accpt.prob.lower)) {
      accpt.prob.lower <- 0.6
    }

    if (length(epsilon) == ncomp*npar_1_comp) {
      tune_param <- matrix(c(epsilon), npar_1_comp, ncomp)
    } else if (length(epsilon) == 1) {
      tune_param <- matrix(epsilon, npar_1_comp, ncomp)
    }
    else {
      stop("epsilon must either be a scalar or a vector of length 2*ncomp (univariate) or 5*ncomp (bivariate)")
    }

  }

  #
  #   if (method == "hmc" & sum(rand_start) == 0 & n.chains > 1 & chains_parallel) {
  #     if (length(L) == 1 | max(L) == min(L)) {
  #       warning(paste("same L accross multiple chains while running them in",
  #                     "parllel with rand_start = FALSE will just make",
  #                     "identical copies of the same chain. L is changed"))
  #       L <- seq(ceiling(L/2), ceiling(2*L), length.out = n.chains)
  #     }
  #   }


  if (grepl(method, "rwmh")) # using rwmh
  {
    if (missing(accpt.prob.upper)) {
      accpt.prob.upper <- 0.5
    }

    if (missing(accpt.prob.lower)) {
      accpt.prob.lower <- 0.3
    }

    if (length(propscale) == ncomp*npar_1_comp) {
      tune_param <- matrix(c(propscale), npar_1_comp, ncomp)
    } else if (length(propscale) == 1) {
      tune_param <- matrix(propscale, npar_1_comp, ncomp)
    }
    else {
      stop("propscale must either be a scalar or a vector of length 2*ncomp (univariate) or 5*ncomp (bivariate)")
    }

  }

  nterms <- tune_ave_size

  if (!model %in% c("wnorm", "wnorm2"))
    int.displ <- omega.2pi <- NULL


  if (model != "vmcos") {
    qrnd_grid <- NULL
    n_qrnd <- NULL
  }
  # Now rename the model specific compiled llik and grad functions


  if (model == "vmsin") {

    if (unimodal.component) {
      dep_cov <- TRUE
      dep_cov_type <- "vmsin_unimodal"
    } else {
      dep_cov <- FALSE
      dep_cov_type <- NULL
    }

    # lpd and grad for all components - not required
    # lpd_grad_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   lpd_grad <- matrix(NA, 6, ncomp)
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       lpd_grad[, j] <- grad_llik_vmsin_C(data[obs_group[[j]], , drop=FALSE],
    #                                          par_mat[, j]) +
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1:2, j] - gam.rate,  -par_mat[3, j]/norm.var, 0, 0,
    #           # lprior
    #           sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #                 gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #         )
    #     } else {
    #       lpd_grad[, j] <-
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1:2, j] - gam.rate,  -par_mat[3, j]/norm.var, 0, 0,
    #           # lprior
    #           sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #                 gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #         )
    #     }
    #   }
    #   list(lpr = sum(lpd_grad[6, ]), grad = lpd_grad[1:5, ])
    # }
    #
    # lpd_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   res <- 0
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       # llik + prior
    #       res <- res +
    #         llik_vmsin_one_comp(data[obs_group[[j]], , drop=FALSE], par_mat[, j],
    #                             log(const_vmsin(par_mat[1, j],
    #                                             par_mat[2, j], par_mat[3, j]))) +
    #         sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #               gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #     } else{
    #       # only prior
    #       res <- res +
    #         sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #               gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #     }
    #   }
    #   res
    # }
    #



    lpd_grad_model_indep_1comp <-  function(data, par_vec,
                                            obs_group, n.clus) {
      lpd_grad <- matrix(NA, 6, 1)
      if (n.clus > 0) {
        lpd_grad <- grad_llik_vmsin_C(data[obs_group, , drop=FALSE],
                                      par_vec) +
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1:2] - gam.rate,  -par_vec[3]/norm.var, 0, 0,
            # lprior
            sum((gam.loc - 1)*log(par_vec[1:2])-
                  gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
          )
      } else {
        lpd_grad <-
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1:2] - gam.rate,  -par_vec[3]/norm.var, 0, 0,
            # lprior
            sum((gam.loc - 1)*log(par_vec[1:2])-
                  gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
          )
      }

      list(lpr = (lpd_grad[6]), grad = lpd_grad[1:5])
    }


    lpd_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      if (n.clus > 0) {
        # llik + prior
        res <-
          llik_vmsin_one_comp(data[obs_group, , drop=FALSE], par_vec,
                              log(const_vmsin(par_vec[1],
                                              par_vec[2], par_vec[3]))) +
          sum((gam.loc - 1)*log(par_vec[1:2])-
                gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
      } else{
        # only prior
        res <-
          sum((gam.loc - 1)*log(par_vec[1:2])-
                gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
      }
      res
    }




    llik_model_contri <- function(data, par_mat, pi_mix) {
      llik_vmsin_contri_C(data, par_mat, pi_mix,
                          log_const_vmsin_all(par_mat))
    }

    mem_p_model <- function(data, par_mat, pi_mix) {
      mem_p_sin(data, par_mat, pi_mix, log_const_vmsin_all(par_mat), 1)
    }
  }

  else if (model == "vmcos") {

    ell <- list(qrnd_grid = qrnd_grid, n_qrnd = n_qrnd)

    if (!is.null(ell$qrnd_grid)) {
      qrnd_grid <- ell$qrnd_grid
      dim_qrnd <- dim(qrnd_grid)
      if (!is.matrix(qrnd_grid) | is.null(dim_qrnd) |
          dim_qrnd[2] != 2)
        stop("qrnd_grid must be a two column matrix")
      n_qrnd <- dim_qrnd[1]
    } else if (!is.null(ell$n_qrnd)){
      n_qrnd <- round(ell$n_qrnd)
      if (n_qrnd < 1)
        stop("n_qrnd must be a positive integer")
      qrnd_grid <- sobol(n_qrnd, 2, FALSE)
    } else {
      n_qrnd <- 1e4
      qrnd_grid <- sobol(n_qrnd, 2, FALSE)
    }

    if (unimodal.component) {
      dep_cov <- TRUE
      dep_cov_type <- "vmcos_unimodal"
    } else {
      dep_cov <- FALSE
      dep_cov_type <- NULL
    }

    # full lpd and grad, not required
    # lpd_grad_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   lpd_grad <- matrix(NA, 6, ncomp)
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       lpd_grad[, j] <- grad_llik_vmcos_C(data[obs_group[[j]], , drop=FALSE],
    #                                          par_mat[, j], qrnd_grid) +
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1:2, j] - gam.rate,  -par_mat[3, j]/norm.var, 0, 0,
    #           # lprior
    #           sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #                 gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #         )
    #     } else {
    #       lpd_grad[, j] <-
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1:2, j] - gam.rate,  -par_mat[3, j]/norm.var, 0, 0,
    #           # lprior
    #           sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #                 gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #         )
    #     }
    #   }
    #   list(lpr = sum(lpd_grad[6, ]), grad = lpd_grad[1:5, ])
    # }
    #
    # lpd_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   res <- 0
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       # llik + prior
    #       res <- res +
    #         llik_vmcos_one_comp(data[obs_group[[j]], , drop=FALSE], par_mat[, j],
    #                             log(const_vmcos(par_mat[1, j],
    #                                             par_mat[2, j],
    #                                             par_mat[3, j],
    #                                             qrnd_grid))) +
    #         sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #               gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #     } else{
    #       # only prior
    #       res <- res +
    #         sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #               gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #     }
    #   }
    #   res
    # }


    lpd_grad_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      lpd_grad <- matrix(NA, 6, 1)
      if (n.clus > 0) {
        lpd_grad[] <- grad_llik_vmcos_C(data[obs_group, , drop=FALSE],
                                        par_vec[], qrnd_grid) +
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1:2] - gam.rate,  -par_vec[3]/norm.var, 0, 0,
            # lprior
            sum((gam.loc - 1)*log(par_vec[1:2])-
                  gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
          )
      } else {
        lpd_grad[] <-
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1:2] - gam.rate,  -par_vec[3]/norm.var, 0, 0,
            # lprior
            sum((gam.loc - 1)*log(par_vec[1:2])-
                  gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
          )
      }
      list(lpr = (lpd_grad[6]), grad = lpd_grad[1:5])
    }

    lpd_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {

      if (n.clus > 0) {
        # llik + prior
        res <-
          llik_vmcos_one_comp(data[obs_group, , drop=FALSE], par_vec[],
                              log(const_vmcos(par_vec[1],
                                              par_vec[2],
                                              par_vec[3],
                                              qrnd_grid))) +
          sum((gam.loc - 1)*log(par_vec[1:2])-
                gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
      } else{
        # only prior
        res <-
          sum((gam.loc - 1)*log(par_vec[1:2])-
                gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
      }
      res
    }



    llik_model_contri <- function(data, par_mat, pi_mix) {
      llik_vmcos_contri_C(data, par_mat, pi_mix,
                          log_const_vmcos_all(par_mat, qrnd_grid))
    }

    mem_p_model <- function(data, par_mat, pi_mix) {
      mem_p_cos(data, par_mat, pi_mix,
                log_const_vmcos_all(par_mat, qrnd_grid))
    }
  }

  else if (model == "wnorm2") {


    dep_cov <- TRUE
    dep_cov_type <- "wnorm2_bound"


    if(int.displ >= 5) int.displ <- 5
    else if(int.displ <= 1) int.displ <- 1

    int.displ <- floor(int.displ)
    omega.2pi.all <- expand.grid(-int.displ:int.displ,-int.displ:int.displ) * (2*pi) # 2pi * integer displacements
    omega.2pi <- as.matrix(omega.2pi.all)

    # lpd_grad_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   lpd_grad <- matrix(NA, 6, ncomp)
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       lpd_grad[, j] <- grad_llik_wnorm2_C(data[obs_group[[j]], , drop=FALSE],
    #                                           par_mat[, j], omega.2pi) +
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1:2, j] - gam.rate,  -par_mat[3, j]/norm.var, 0, 0,
    #           # lprior
    #           sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #                 gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #         )
    #     } else {
    #       lpd_grad[, j] <-
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1:2, j] - gam.rate,  -par_mat[3, j]/norm.var, 0, 0,
    #           # lprior
    #           sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #                 gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #         )
    #     }
    #   }
    #   list(lpr = sum(lpd_grad[6, ]), grad = lpd_grad[1:5, ])
    # }

    # lpd_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   res <- 0
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       # llik + prior
    #       res <- res +
    #         llik_wnorm2_one_comp(data[obs_group[[j]], , drop=FALSE], par_mat[, j],
    #                              l_const_wnorm2(par_mat[, j]),
    #                              omega.2pi) +
    #         sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #               gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #     } else{
    #       # only prior
    #       res <- res +
    #         sum((gam.loc - 1)*log(par_mat[1:2, j])-
    #               gam.rate*par_mat[1:2, j]) - 0.5*par_mat[3, j]^2/norm.var
    #     }
    #   }
    #   unname(res)
    # }


    lpd_grad_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      lpd_grad <- matrix(NA, 6, 1)
      if (n.clus > 0) {
        lpd_grad[] <- grad_llik_wnorm2_C(data[obs_group, , drop=FALSE],
                                         par_vec[], omega.2pi) +
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1:2] - gam.rate,  -par_vec[3]/norm.var, 0, 0,
            # lprior
            sum((gam.loc - 1)*log(par_vec[1:2])-
                  gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
          )
      } else {
        lpd_grad[] <-
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1:2] - gam.rate,  -par_vec[3]/norm.var, 0, 0,
            # lprior
            sum((gam.loc - 1)*log(par_vec[1:2])-
                  gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
          )
      }
      list(lpr = sum(lpd_grad[6]), grad = lpd_grad[1:5])
    }

    lpd_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      if (n.clus > 0) {
        # llik + prior
        res <-
          llik_wnorm2_one_comp(data[obs_group, , drop=FALSE], par_vec[],
                               l_const_wnorm2(par_vec[]),
                               omega.2pi) +
          sum((gam.loc - 1)*log(par_vec[1:2])-
                gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
      } else{
        # only prior
        res <-
          sum((gam.loc - 1)*log(par_vec[1:2])-
                gam.rate*par_vec[1:2]) - 0.5*par_vec[3]^2/norm.var
      }
      unname(res)
    }


    llik_model_contri <- function(data, par_mat, pi_mix) {
      llik_wnorm2_contri_C(data, par_mat, pi_mix,
                           log_const_wnorm2_all(par_mat), omega.2pi)
    }

    mem_p_model <- function(data, par_mat, pi_mix) {
      mem_p_wnorm2(data, par_mat, pi_mix, log_const_wnorm2_all(par_mat),
                   omega.2pi)
    }


  }

  else if (model == "vm") {

    dep_cov <- FALSE

    # lpd_grad_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   lpd_grad <- matrix(NA, 3, ncomp)
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       lpd_grad[, j] <- grad_llik_univm_C(data[obs_group[[j]]],
    #                                          par_mat[, j]) +
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1, j] - gam.rate,  0,
    #           # lprior
    #           (gam.loc - 1)*log(par_mat[1, j])-
    #             gam.rate*par_mat[1, j]
    #         )
    #
    #     } else {
    #       lpd_grad[, j] <-
    #         c( # grad for lprior
    #           (gam.loc - 1)/par_mat[1, j] - gam.rate,  0,
    #           # lprior
    #           (gam.loc - 1)*log(par_mat[1, j])-
    #             gam.rate*par_mat[1, j]
    #         )
    #     }
    #   }
    #   list(lpr = sum(lpd_grad[3, ]), grad = lpd_grad[1:2, ])
    # }
    #
    # lpd_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #   res <- 0
    #   for(j in 1:ncomp) {
    #     if (n.clus[j] > 0) {
    #       # llik + prior
    #       res <- res +
    #         llik_univm_one_comp(data[obs_group[[j]]], par_mat[, j],
    #                             log(const_univm(par_mat[1, j]))) +
    #         (gam.loc - 1)*log(par_mat[1, j]) - gam.rate*par_mat[1, j]
    #
    #     } else{
    #       # only prior
    #       res <- res +
    #         (gam.loc - 1)*log(par_mat[1, j]) - gam.rate*par_mat[1, j]
    #     }
    #   }
    #   unname(res)
    # }

    lpd_grad_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      lpd_grad <- matrix(NA, 3, 1)
      if (n.clus > 0) {
        lpd_grad[] <- grad_llik_univm_C(data[obs_group],
                                        par_vec[]) +
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1] - gam.rate,  0,
            # lprior
            (gam.loc - 1)*log(par_vec[1])-
              gam.rate*par_vec[1]
          )

      } else {
        lpd_grad[] <-
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1] - gam.rate,  0,
            # lprior
            (gam.loc - 1)*log(par_vec[1])-
              gam.rate*par_vec[1]
          )
      }
      list(lpr = (lpd_grad[3 ]), grad = lpd_grad[1:2 ])
    }

    lpd_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      if (n.clus > 0) {
        # llik + prior
        res <-
          llik_univm_one_comp(data[obs_group], par_vec[],
                              log(const_univm(par_vec[1]))) +
          (gam.loc - 1)*log(par_vec[1]) - gam.rate*par_vec[1]

      } else{
        # only prior
        res <-
          (gam.loc - 1)*log(par_vec[1]) - gam.rate*par_vec[1]
      }
      unname(res)
    }



    llik_model_contri <- function(data, par_mat, pi_mix) {
      llik_univm_contri_C(data, par_mat, pi_mix,
                          log_const_univm_all(par_mat))
    }

    mem_p_model <- function(data, par_mat, pi_mix) {
      mem_p_univm(data, par_mat, pi_mix, log_const_univm_all(par_mat))
    }



  }

  # else if (model == "wnorm")
  else    {

    dep_cov <- FALSE


    if(int.displ >= 5) int.displ <- 5
    else if(int.displ <= 1) int.displ <- 1

    int.displ <- floor(int.displ)
    omega.2pi <- (-int.displ):int.displ * (2*pi) # 2pi * 1d integer displacements

    #
    #     lpd_grad_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #       lpd_grad <- matrix(NA, 3, ncomp)
    #       for(j in 1:ncomp) {
    #         if (n.clus[j] > 0) {
    #           lpd_grad[, j] <- grad_llik_uniwnorm_C(data[obs_group[[j]]],
    #                                                 par_mat[, j], omega.2pi) +
    #             c( # grad for lprior
    #               (gam.loc - 1)/par_mat[1, j] - gam.rate,  0,
    #               # lprior
    #               (gam.loc - 1)*log(par_mat[1, j])-
    #                 gam.rate*par_mat[1, j]
    #             )
    #
    #         } else {
    #           lpd_grad[, j] <-
    #             c( # grad for lprior
    #               (gam.loc - 1)/par_mat[1, j] - gam.rate,  0,
    #               # lprior
    #               (gam.loc - 1)*log(par_mat[1, j])-
    #                 gam.rate*par_mat[1, j]
    #             )
    #         }
    #       }
    #       list(lpr = sum(lpd_grad[3, ]), grad = lpd_grad[1:2, ])
    #     }
    #
    #     lpd_model_indep <- function(data, par_mat, obs_group, n.clus) {
    #       res <- 0
    #       for(j in 1:ncomp) {
    #         if (n.clus[j] > 0) {
    #           # llik + prior
    #           res <- res +
    #             llik_uniwnorm_one_comp(data[obs_group[[j]]], par_mat[, j],
    #                                    l_const_uniwnorm(par_mat[1, j]),
    #                                    omega.2pi) +
    #             (gam.loc - 1)*log(par_mat[1, j]) - gam.rate*par_mat[1, j]
    #
    #         } else{
    #           # only prior
    #           res <- res +
    #             (gam.loc - 1)*log(par_mat[1, j]) - gam.rate*par_mat[1, j]
    #         }
    #       }
    #       unname(res)
    #     }

    lpd_grad_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      lpd_grad <- matrix(NA, 3, 1)
      if (n.clus > 0) {
        lpd_grad[] <- grad_llik_uniwnorm_C(data[obs_group],
                                           par_vec[], omega.2pi) +
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1] - gam.rate,  0,
            # lprior
            (gam.loc - 1)*log(par_vec[1])-
              gam.rate*par_vec[1]
          )

      } else {
        lpd_grad[] <-
          c( # grad for lprior
            (gam.loc - 1)/par_vec[1] - gam.rate,  0,
            # lprior
            (gam.loc - 1)*log(par_vec[1])-
              gam.rate*par_vec[1]
          )
      }

      list(lpr = sum(lpd_grad[3, ]), grad = lpd_grad[1:2, ])
    }

    lpd_model_indep_1comp <- function(data, par_vec, obs_group, n.clus) {
      if (n.clus > 0) {
        # llik + prior
        res <-
          llik_uniwnorm_one_comp(data[obs_group], par_vec[],
                                 l_const_uniwnorm(par_vec[1]),
                                 omega.2pi) +
          (gam.loc - 1)*log(par_vec[1]) - gam.rate*par_vec[1]

      } else{
        # only prior
        res <-
          (gam.loc - 1)*log(par_vec[1]) - gam.rate*par_vec[1]
      }

      unname(res)
    }



    llik_model_contri <- function(data, par_mat, pi_mix) {
      llik_uniwnorm_contri_C(data, par_mat, pi_mix,
                             log_const_uniwnorm_all(par_mat), omega.2pi)
    }

    mem_p_model <- function(data, par_mat, pi_mix) {
      mem_p_uniwnorm(data, par_mat, pi_mix, log_const_uniwnorm_all(par_mat), omega.2pi)
    }


    # grad_llik_uniwnorm_R <- function(data, par) {
    #   gr <- numDeriv::grad(function(par) llik_uniwnorm_one_comp(data, par,
    #                                                             l_const_uniwnorm(par[1]),
    #                                                             omega.2pi),
    #                        par)
    #   gr[3] <- llik_uniwnorm_one_comp(data, par,
    #                                   l_const_uniwnorm(par[1]), omega.2pi)
    #   gr
    # }

    # llik_uniwnorm_full(data.rad, par.mat, pi.mix, log_c = log_const_uniwnorm_all(par.mat), omega.2pi)
    #
    #
    # grad_llik_uniwnorm_C(data.rad, par.mat[, 1], omega.2pi)
    # grad_llik_uniwnorm_R(data.rad, par.mat[, 1])

  }




  if (!is.null(start_par)) {
    if (!is.list(start_par))
      stop("start_par must be a list")
    if (!is.list(start_par[[1]])) {
      if (length(setdiff(c(par.names, "pmix"), names(start_par))) > 0)
        stop("start_par does not have all parameters")
      start_par <- lapply(1:n.chains, function(ii) start_par)
    } else {
      if (length(start_par) != n.chains)
        stop("length(start_par) must be either 1 or equal to n.chains")
      for(jj in 1:n.chains) {
        if (length(setdiff(c(par.names, "pmix"), names(start_par[[jj]]))) > 0)
          stop(paste0("start_par[[", jj,   "]] does not have all parameters"))
      }
    }
  }

  if (any(is.null(start_par), !is.list(start_par[[1]]))
      & all(!rand_start) ) {
    starting_1 <- process_startpar(start_par,
                                   data.rad,
                                   ncomp,
                                   model,
                                   FALSE)
    starting <- lapply(1:n.chains, function(j) starting_1)
  } else {
    starting <- lapply(1:n.chains,
                       function(j)
                         process_startpar(start_par[[j]],
                                          data.rad,
                                          ncomp,
                                          model,
                                          rand_start[j]))
  }

  run_MC <- function(starting, L, chain_no) {

    # just to change the RNG state across chains, so that
    # no two c chains turn out to be identical
    change_rng_state <- runif(chain_no)

    if (method == "hmc") # using hmc
    {
      L_vec <- rep(L, n.iter)
      if (L.random)
        L.orig <- L
    }


    starting$par.mat[abs(starting$par.mat) >= kappa_upper/2] <- kappa_upper/2
    starting$par.mat[abs(starting$par.mat) <= 2*kappa_lower] <- 2*kappa_lower

    par.mat.all <- array(0, dim = c(npar_1_comp, ncomp, n.iter))
    pi.mix.all <- matrix(1, nrow = ncomp, ncol = n.iter)
    llik.all <- lprior.all <- lpd.all <- rep(-Inf, (n.iter))
    accpt.par.mat.all <- matrix(NA, ncomp, n.iter)
    modelpar.names <- par.names
    clus_ind.all <- matrix(1, nrow = n.data, ncol = n.iter)
    if (return_llik_contri) {
      llik_contri_all <- matrix(1, nrow = n.data, ncol = n.iter)
    } else {
      llik_contri_all <- NULL
    }
    mem_prob_all <- array(1, dim = c(n.data, ncomp, n.iter))



    epsilon_ave <- NULL
    L_ave <- NULL
    propscale_final <- NULL



    pi.mix <- starting$pi.mix
    par.mat <- as.matrix(starting$par.mat)
    if (zero_cov) {
      par.mat[3, ] <- 0
    } else if (cov.restrict == "POSITIVE") {
      par.mat[3, ] <- pmax(par.mat[3, ], 0)
    } else if (cov.restrict == "NEGATIVE") {
      par.mat[3, ] <- pmin(par.mat[3, ], 0)
    }


    # browser()

    if (ncomp == 1) {
      perm_sampling <- FALSE
      clus.ind <- clus.ind_curr <- rep(1, n.data)
      obs_group <- list(1:n.data)
      n.clus <- n.data
      post.wt <- matrix(1, n.data, 1)
      pi.mix <- 1
    }


    tcltk_fail <- FALSE

    if (show.progress & !exists("pb")) {
      pb <- tryCatch(tcltk::tkProgressBar(paste("Chain", chain_no),
                                          min = 1, max = n.iter),
                     error = function(e) "error")
      if (unlist(pb)[1] == "error") {
        show.progress <- FALSE
        tcltk_fail <- TRUE
      }
    }



    par.mat.all.order <- par.mat.all

    ntunes_up <- ntunes_down <- rep(0, ncomp)
    # tune_iter_no <- c()
    # ave_accpt_all <- c()

    tune_param_all <- matrix(NA, length(tune_param), n.iter)

    tune_status <- matrix(0, ncomp, n.iter)

    #  Run the Markov chain
    for(iter in  1:n.iter) {
      if (method == "hmc") {
        if (epsilon.random) {
          tune_param_final <- tune_param*runif(1, 1-epsilon.incr, 1+epsilon.incr)
        } else {
          tune_param_final <- tune_param
        }
        if (L.random)
          L_vec[iter] <- L <- sample(ceiling(L.orig/exp(L.incr)):ceiling(L.orig*exp(L.incr)), 1)
      }

      #----------------------------------------------------------------------------------
      # generating mixture proportions if ncomp > 1
      #----------------------------------------------------------------------------------


      if(ncomp > 1) {
        post.wt <- mem_p_model(data.rad, par.mat, pi.mix)
        clus.ind <-  cID(post.wt, ncomp, runif(n.data))
        # n.clus <- tabulate(clus.ind_curr, nbins = ncomp)
        obs_group <- lapply(1:ncomp, function(j) which(clus.ind == j))
        n.clus <- listLen(obs_group)
        pi.mix <-
          as.numeric(rdirichlet(1, (pmix.alpha + n.clus))) #new mixture proportions
      }




      #----------------------------------------------------------------------------------
      # generating par.mat
      #----------------------------------------------------------------------------------

      if (type == "bi") {



        if (method == "hmc") {

          # hmc_curr <-
          #   bounded_hmc_biv(lpr_grad =
          #                     function(par_mat)
          #                       lpd_grad_model_indep(data.rad, par_mat,
          #                                            obs_group, n.clus),
          #                   init =  par.mat,
          #                   lower = par_lower,
          #                   upper = par_upper,
          #                   dep_cov = dep_cov,
          #                   dep_cov_type = dep_cov_type,
          #                   zero_cov = zero_cov,
          #                   nsteps = L,
          #                   step = tune_param_final)
          #
          # par.mat <- hmc_curr$final
          # accpt.par.mat.all[iter] <- hmc_curr$acc

          # par_mat = par.mat
          # lpd_grad_model_indep(data.rad, par_mat,
          #                      obs_group, n.clus)
          #
          #
          # lpd_grad_model_indep_1comp(data.rad, par_mat[, 1],
          #                      obs_group[[1]], n.clus[1])
          # lpd_grad_model_indep_1comp(data.rad, par_mat[, 2],
          #                            obs_group[[2]], n.clus[2])
          #
          #
          # lpd_model_indep(data.rad, par_mat,
          #                 obs_group, n.clus)
          #
          # lpd_model_indep_1comp(data.rad, par_mat[, 1],
          #                 obs_group[[1]], n.clus[1])
          #

          # browser()

          for (j in 1:ncomp) {
            hmc_curr <-
              bounded_hmc_biv(lpr_grad =
                                function(par_vec)
                                  lpd_grad_model_indep_1comp(data.rad, par_vec,
                                                             obs_group[[j]], n.clus[j]),
                              init =  par.mat[, j, drop=FALSE],
                              lower = par_lower[, j, drop=FALSE],
                              upper = par_upper[, j, drop=FALSE],
                              dep_cov = dep_cov,
                              dep_cov_type = dep_cov_type,
                              zero_cov = zero_cov,
                              nsteps = L,
                              step = tune_param_final[, j, drop=FALSE])

            par.mat[, j] <- hmc_curr$final
            accpt.par.mat.all[j, iter] <- hmc_curr$acc
          }
        }

        else {

          # rwmh_curr <- bounded_rwmh_biv(lpr =
          #                                 function(par_mat)
          #                                   lpd_model_indep(data.rad, par_mat,
          #                                                   obs_group, n.clus),
          #                               init =  par.mat,
          #                               lower = par_lower,
          #                               upper = par_upper,
          #                               dep_cov = dep_cov,
          #                               dep_cov_type = dep_cov_type,
          #                               zero_cov = zero_cov,
          #                               step = tune_param)
          #
          # par.mat <- rwmh_curr$final
          # accpt.par.mat.all[iter] <- rwmh_curr$accpt

          for(j in 1:ncomp) {
            rwmh_curr <- bounded_rwmh_biv(lpr =
                                            function(par_vec)
                                              lpd_model_indep_1comp(data.rad, par_vec,
                                                                    obs_group[[j]], n.clus[j]),
                                          init =  par.mat[, j, drop=FALSE],
                                          lower = par_lower[, j, drop=FALSE],
                                          upper = par_upper[, j, drop=FALSE],
                                          dep_cov = dep_cov,
                                          dep_cov_type = dep_cov_type,
                                          zero_cov = zero_cov,
                                          step = tune_param[, j])

            par.mat[, j] <- rwmh_curr$final
            accpt.par.mat.all[j, iter] <- rwmh_curr$accpt
          }

        }


        lprior.all[iter] <-
          sum((gam.loc - 1)*log(par.mat[1:2, ])-
                gam.rate*par.mat[1:2, ]) -
          0.5*sum(par.mat[3, ]^2)/norm.var +
          sum(pmix.alpha*log(pi.mix))

      }

      # if type == "uni"
      else {

        if (method == "hmc") {
          # hmc_curr <-
          #   bounded_hmc_uni(lpr_grad =
          #                     function(par_mat)
          #                       lpd_grad_model_indep(data.rad, par_mat,
          #                                            obs_group, n.clus),
          #                   init =  par.mat,
          #                   lower = par_lower,
          #                   upper = par_upper,
          #                   nsteps = L,
          #                   step = tune_param_final)
          #
          # par.mat <- hmc_curr$final
          # accpt.par.mat.all[iter] <- hmc_curr$acc

          for (j in 1:ncomp) {
            hmc_curr <-
              bounded_hmc_uni(lpr_grad =
                                function(par_vec)
                                  lpd_grad_model_indep_1comp(data.rad, par_vec,
                                                             obs_group[[j]], n.clus[j]),
                              init =  par.mat[, j, drop=FALSE],
                              lower = par_lower[, j, drop=FALSE],
                              upper = par_upper[, j, drop=FALSE],
                              nsteps = L,
                              step = tune_param_final[, j, drop=FALSE])

            par.mat[, j] <- hmc_curr$final
            accpt.par.mat.all[j, iter] <- hmc_curr$acc
          }

        }

        else {

          # rwmh_curr <- bounded_rwmh_uni(lpr =
          #                                 function(par_mat)
          #                                   lpd_model_indep(data.rad, par_mat,
          #                                                   obs_group, n.clus),
          #                               init =  par.mat,
          #                               lower = par_lower,
          #                               upper = par_upper,
          #                               step = tune_param)
          #
          # par.mat <- rwmh_curr$final
          # accpt.par.mat.all[iter] <- rwmh_curr$accpt

          for(j in 1:ncomp) {
            rwmh_curr <- bounded_rwmh_uni(lpr =
                                            function(par_vec)
                                              lpd_model_indep_1comp(data.rad, par_vec,
                                                                    obs_group[[j]], n.clus[j]),
                                          init =  par.mat[, j, drop=FALSE],
                                          lower = par_lower[, j, drop=FALSE],
                                          upper = par_upper[, j, drop=FALSE],
                                          step = tune_param[, j])

            par.mat[, j] <- rwmh_curr$final
            accpt.par.mat.all[j, iter] <- rwmh_curr$accpt
          }
        }

        lprior.all[iter] <-
          sum((gam.loc - 1)*log(par.mat[1, ])-
                gam.rate*par.mat[1, ]) +
          sum(pmix.alpha*log(pi.mix))
      }




      # do permutation sampling only after tuning
      if (perm_sampling & iter > n.burnin) {
        rand_perm <- sample(1:ncomp)
        clus.ind <- rand_perm[clus.ind] # random label switch
        post.wt <- post.wt[, rand_perm, drop=FALSE]
        par.mat <- par.mat[, rand_perm, drop=FALSE]
        pi.mix <- pi.mix[rand_perm, drop=FALSE]

        par.mat.all.order[, , iter] <- par.mat[, order(rand_perm), drop=FALSE]
        # needed for tuning

        # if (method == "hmc")
        tune_param <- tune_param[, rand_perm, drop=FALSE]

      }

      # browser()

      llik_contri <- llik_model_contri(data.rad, par.mat, pi.mix)

      if (return_llik_contri)
        llik_contri_all[, iter] <- llik_contri

      llik.all[iter] <- sum(llik_contri)
      lpd.all[iter] <- llik.all[iter] + lprior.all[iter]

      par.mat.all[, , iter] <- par.mat
      pi.mix.all[, iter] <- pi.mix
      mem_prob_all[, , iter] <- post.wt
      clus_ind.all[, iter] <- clus.ind




      tune_param_all[, iter] <- c(tune_param)

      # tuning tune_param during burnin
      if (autotune & iter >= nterms & iter %% 10 == 0 &
          iter <= iter.tune) {
        ave_accpt_all <- rep(0, ncomp)
        for(j in 1:ncomp) {
          ave_accpt_all[j] <-
            ave_accpt <- sum(accpt.par.mat.all[j, (iter-nterms+1):iter])/nterms

          if (ave_accpt > accpt.prob.upper) {
            tune_param[, j] <- tune_param[, j] * (1 + tune.incr)

            ntunes_up[j] <- ntunes_up[j] + 1
            tune_status[j, iter] <- 1

          }
          else if (ave_accpt < accpt.prob.lower) {
            tune_param[, j] <- tune_param[, j] * (1 - tune.incr)

            ntunes_down[j] <- ntunes_down[j] + 1
            tune_status[j, iter] <- -1
          }

          if (iter %% nterms == 0) {

            par.sd <- apply(par.mat.all[, j, (iter-nterms+1):iter], 1, sd)
            mean_par.sd <- sum(par.sd)/(npar_1_comp)
            mean_tune_param_j <- sum(tune_param[, j])/(npar_1_comp)
            if(mean_par.sd > 0)
              tune_param[, j] <- 0.5*par.sd/mean_par.sd*mean_tune_param_j + 0.5*tune_param[, j]

            #   else {
            #   par.sd <- apply(par.mat.all.order[, , (iter-nterms+1):iter], c(1:2), sd)
            #   par.sd <- par.sd[, rand_perm]
            #   mean_par.sd <- sum(par.sd)/(npar_1_comp*ncomp)
            #   mean_tune_param <- sum(tune_param)/(npar_1_comp*ncomp)
            #   if(mean_par.sd > 0)
            #     tune_param <- 0.5*par.sd/mean_par.sd*mean_tune_param + 0.5*tune_param
            #
            #   cat(ave_accpt, tune_param[1,],lpd.all[iter], "\n")
            # }
          }
        }

        # if(iter %% nterms == 0)
        #   cat(paste("(", paste0(ave_accpt_all, collapse = ","), ")"),
        #       tune_param[1,], lpd.all[iter], "\n")
      }

      # if (method == "hmc" & autotune & iter >= nterms
      #     & iter %% 5 == 0 & iter <= n.burnin) {
      #   acr <- cor(lpd.all[(iter-nterms+2):iter],
      #              lpd.all[(iter-nterms+1):(iter-1)])
      #   if (acr < acr.lower) {
      #     L <- ceiling(L*(1 - tune.incr))
      #   } else if (acr > acr.upper) {
      #     L <- ceiling(L*(1 + tune.incr))
      #   }
      # }



      if (show.progress)
        if ((round(iter /  n.iter, 2) * 100) %% 5 == 0 || iter == n.iter + 1)
          tcltk::setTkProgressBar(pb, iter)

    }


    if(grepl(method, "hmc")) {
      epsilon_ave <- epsilon <- mean(tune_param_all)

      if(L.random) {
        L_ave <- sum(L_vec)/n.iter
      } else {
        L_ave <- L
      }
    }

    if(grepl(method, "rwmh")) {
      propscale_final <- propscale <- rowSums(tune_param_all)/n.iter
    }

    allpar_val <- array(1, dim = c(npar_1_comp+1, ncomp, n.iter))
    allpar_val[1, , ] <- pi.mix.all
    allpar_val[-1, , ] <- par.mat.all
    rm(pi.mix.all, par.mat.all)

    allpar_name <- c("pmix", modelpar.names)
    dimnames(allpar_val)[[1]] <- c("pmix", modelpar.names)

    if (!return_tune_param) {
      rm(tune_param_all)
      tune_param_all <- NULL
    }


    result <- list("par.value" = allpar_val, #[, , final_iter_set],
                   "par.name" = allpar_name,
                   "llik.contri" = llik_contri_all, #[, final_iter_set],
                   "llik" = llik.all, #[final_iter_set],
                   "lpd" = lpd.all,#[final_iter_set],
                   "lprior" = lprior.all,#[final_iter_set],
                   "accpt.modelpar" = accpt.par.mat.all,#[final_iter_set],
                   "clus.ind" = clus_ind.all,#[, final_iter_set],
                   "mem.prob" = mem_prob_all,#[, , final_iter_set],
                   "epsilon" = epsilon_ave,
                   "L" = L_ave,
                   "propscale" = propscale_final,
                   "tune_param" = tune_param_all,
                   "par.upper" = par_upper,
                   "par.lower" = par_lower,
                   "tcltk_fail" = tcltk_fail)

    if(show.progress) close(pb)

    result
  }


  # browser()
  # generate three chains in parallel, if possible
  if (chains_parallel) {

    res_list <- future.apply::future_lapply(1:n.chains,
                                            function(ii) run_MC(starting[[ii]],
                                                                L[ii], ii),
                                            future.seed = TRUE)
  } else {
    res_list <- lapply(1:n.chains,
                       function(ii) run_MC(starting[[ii]], L[ii], ii))
  }



  if (res_list[[1]]$tcltk_fail)
    warning("tcltk could not be loaded. \'show.progress\' was set to FALSE.")

  # combine the results from the lists
  allpar_val <- array(0, dim=c(npar_1_comp+1, ncomp, n.iter, n.chains))
  llik_all <- lprior_all <- lpd_all <- matrix(0, n.iter, n.chains)
  accpt.modelpar_all <- array(0, dim = c(ncomp, n.iter, n.chains))
  clus.ind_all <- array(0, dim=c(n.data, n.iter, n.chains))
  if (return_llik_contri) {
    llik_contri_all <- array(0, dim=c(n.data, n.iter, n.chains))
  } else {
    llik_contri_all <- NULL
  }

  mem_prob_all <- array(1, dim=c(n.data, ncomp, n.iter, n.chains))

  if (return_tune_param) {
    tune_param_all <- array(NA, dim=c(length(tune_param), n.iter, n.chains))
  } else {
    tune_param_all <- NULL
  }

  for(j in 1:n.chains) {
    allpar_val[, , , j] <- res_list[[j]]$par.value
    llik_all[, j] <- res_list[[j]]$llik
    lprior_all[, j] <- res_list[[j]]$lprior
    lpd_all[, j] <- res_list[[j]]$lpd
    accpt.modelpar_all[, , j] <- res_list[[j]]$accpt.modelpar
    clus.ind_all[, , j] <- res_list[[j]]$clus.ind
    if (return_llik_contri)
      llik_contri_all[, , j] <- res_list[[j]]$llik.contri
    mem_prob_all[, , , j] <- res_list[[j]]$mem.prob
    if (return_tune_param)
      tune_param_all[, , j] <- res_list[[j]]$tune_param
  }

  if(method == "hmc") {
    epsilon_final <- do.call(cbind, lapply(res_list, function(x) x$epsilon))
    L_final <- unlist(lapply(res_list, function(x) x$L))
    propscale_final <- NULL
  }
  else {
    propscale_final <- do.call(cbind, lapply(res_list, function(x) x$propscale))
    epsilon_final <- NULL
    L_final <- NULL
  }

  out <- list("par.value" = allpar_val,
              "clus.ind" = clus.ind_all,
              "par.name" = res_list[[1]]$par.name,
              "modelpar.lower" = res_list[[1]]$par.lower,
              "modelpar.upper" = res_list[[1]]$par.upper,
              "return_llik_contri" = return_llik_contri,
              "llik.contri" = llik_contri_all,
              "mem.prob" = mem_prob_all,
              "llik" = llik_all,
              "lpd" = lpd_all,
              "lprior" = lprior_all,
              "accpt.modelpar" = accpt.modelpar_all,
              "model" = curr.model,
              "method" = method,
              "perm_sampling" = perm_sampling,
              "epsilon.random" = epsilon.random,
              "epsilon" = epsilon_final,
              "L.random" = L.random,
              "L" = L_final,
              "iter.tune" = iter.tune,
              "propscale" = propscale_final,
              "tune_param" = tune_param_all,
              "return_tune_param" = return_tune_param,
              "type" = type,
              "data" = data.rad,
              "cov.restrict" = cov.restrict,
              "gam.loc" = gam.loc,
              "gam.scale" = gam.scale,
              "pmix.alpha" = pmix.alpha,
              "norm.var" = norm.var,
              "n.data" = n.data,
              "ncomp" = ncomp,
              "n.chains" = n.chains,
              "n.iter" = n.iter,
              "n.burnin" = n.burnin,
              "thin" = thin,
              "n.iter.final" = n.iter.final,
              "final_iter" = final_iter_set,
              "int.displ" = int.displ,
              "qrnd_grid" = qrnd_grid,
              "omega.2pi" = omega.2pi)

  class(out) <- "angmcmc"

  out
}
