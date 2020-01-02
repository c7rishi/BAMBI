
#' @export
vm2_mle <- function(data, model = c("vmsin", "vmcos", "indep"), ...) {
  model <- model[1]

  dots <- list(...)
  data <- data.matrix(data)
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
  names(start) <- c("kappa1", "kapp2", "kappa3", "mu1", "mu2")

  start_lscale <- start
  start_lscale[c("kappa1", "kapp2")] <- log(start[c("kappa1", "kapp2")])

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
    method = "L-BFGS-B"
    # hessian = TRUE
  )

  opt_adj <- opt
  opt_adj$par[c("kappa1", "kapp2")] <- exp(opt$par[c("kappa1", "kapp2")])

  hess <- hessian_fn(par_vec = opt_adj$par)
  dimnames(hess) <- list(names(opt_adj$par), names(opt_adj$par))
  opt_adj$hessian <- hess

  opt_adj
}
