#' Maximum likelihood estimation of bivariate von Mises parameters
#' @inheritParams fit_angmix
#' @param model Bivariate von Mises model. Must either be "vmsin" or "vmcos"
#' @param ... Additional arguments. See details.
#' @details The parameters \code{kappa1} and \code{kappa2} are optimized
#' in log scales. The method of optimization used (passed to \link{optim})
#' can be specified through \code{method} in \code{...}
#' (defaults to \code{"L-BFGS-B"}). Note, however, that
#' lower (0)  and upper (2*pi) bounds for \code{mu1} and \code{mu2}
#' are specified; so not all methods implemented in \link{optim} will work.
#' @return An object of class \link{mle-class}.
#' @examples
#' pars <- list(kappa1 = 3, kappa2 = 2, kappa3 = 1.5, mu1 = 0.5, mu2 = 1.5)
#' nsamp <- 2000
#' model <- "vmsin"
#' set.seed(100)
#' dat_gen <- do.call(paste0("r", model), c(list(n = nsamp), pars))
#'
#' est <- vm2_mle(dat_gen, model = model)
#' library(stats4)
#' coef(est)
#' vcov(est)
#' @export
vm2_mle <- function(data, model = c("vmsin", "vmcos"), ...) {

  model <- model[1]
  dots <- list(...)
  data <- data.matrix(data)
  call <- match.call()

  if (is.null(dots$method)) {
    dots$method <- "L-BFGS-B"
  }
  method <- dots$method

  if (model == "vmsin") {

    # if (unimodal.component) {
    #   dep_cov <- TRUE
    #   dep_cov_type <- "vmsin_unimodal"
    # } else {
    #   dep_cov <- FALSE
    #   dep_cov_type <- NULL
    # }


    # in log scale
    lpd_grad_model_indep_1comp <-  function(par_vec_lscale) {
      par_vec <- c(exp(par_vec_lscale[1:2]), par_vec_lscale[3:5])
      lpd_grad <- matrix(NA, 6, 1)

      lpd_grad <- signif(
        suppressWarnings(grad_llik_vmsin_C(data, par_vec))*
          c(par_vec[1:2], rep(1, 4)),
        8
      )
      list(lpr = (lpd_grad[6]), grad = lpd_grad[1:5])
    }

    start_par_gen <- start_par_vmsin

    hessian_fn <- function(par_vec) {
      numDeriv::hessian(
        func = function(par_vec) {
          -grad_llik_vmsin_C(data, par_vec)[6]
        },
        x = par_vec
      )
    }

  }

  else if (model == "vmcos") {

    ell <- dots[c("qrnd_grid", "n_qrnd")] #list(qrnd_grid = qrnd, n_qrnd = n_qrnd)

    if (!is.null(ell$qrnd)) {
      qrnd_grid <- ell$qrnd
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

    # in log scale
    lpd_grad_model_indep_1comp <- function(par_vec_lscale) {
      par_vec <- c(exp(par_vec_lscale[1:2]), par_vec_lscale[3:5])
      lpd_grad <- matrix(NA, 6, 1)
      lpd_grad[] <- signif(
        suppressWarnings(grad_llik_vmcos_C(data, par_vec, qrnd_grid)) *
          c(par_vec[1:2], rep(1, 4)),
        8
      )

      list(lpr = (lpd_grad[6]), grad = lpd_grad[1:5])
    }

    start_par_gen <- start_par_vmcos

    hessian_fn <- function(par_vec) {
      numDeriv::hessian(
        func = function(par_vec) {
          -grad_llik_vmcos_C(data, par_vec, qrnd_grid)[6]
        },
        x = par_vec
      )
    }

  }
  # browser()


  start <- start_par_gen(data)
  names(start) <- c("log_kappa1", "log_kappa2", "kappa3", "mu1", "mu2")

  start_lscale <- start
  start_lscale[c("log_kappa1", "log_kappa2")] <-
    log(start[c("log_kappa1", "log_kappa2")])

  opt <- optim(
    par = start_lscale,
    fn = function(par_lscale) {
      -lpd_grad_model_indep_1comp(par_lscale)$lpr
    },
    gr = function(par_lscale) {
      -lpd_grad_model_indep_1comp(par_lscale)$grad
    },
    lower = c(rep(-Inf, 3), 0, 0),
    upper = c(rep(Inf, 3), 2*pi, 2*pi),
    method = method
    # hessian = TRUE
  )

  est_par <- opt$par
  names(est_par)[1:2] <- c("kappa1", "kappa2")
  est_par[c("kappa1", "kappa2")] <- exp(est_par[c("kappa1", "kappa2")])

  hess <- hessian_fn(par_vec = est_par)
  dimnames(hess) <- list(names(est_par), names(est_par))

  res <- new(
    "mle",
    call = call,
    coef = est_par,
    fullcoef = unlist(est_par),
    vcov = solve(hess),
    min = opt$value,
    details = opt,
    minuslogl = function(kappa1, kappa2, kappa3, mu1, mu2) {
      par_lscale <- c(log(kappa1), log(kappa2), kappa3, mu1, mu2)
      -lpd_grad_model_indep_1comp(par_lscale)$lpr
    },
    nobs = nrow(data),
    method = method
  )

 res
}