#' Sample circular correlation coefficients
#' @param x two column matrix. NA values are not allowed.
#' @param type type of the circular correlation.
#' Must be one of "fl", "js", "tau1" and "tau2". See details.
#' @param alternative one of \code{"two.sided"}, \code{"less"} or \code{"greater"}
#' (defaults to \code{"two.sided"}).
#' Hypothesis test is performed only when \code{type} is either \code{"fl"} or \code{"js"},
#' in which case asymptotic standard error of the estimator is used to construct the test
#' statistic.
#' @param jackknife logical. Compute jackknifed estimate and standard error? Defaults to FALSE.
#' @param bootse logical. Compute bootstrap standard error? Defaults to FALSE.
#' @param n.boot number of bootstrapped samples to compute bootstrap standard error. Defaults to
#' 100. Ignored if \code{bootse} if FALSE.
#'
#' @details
#' \code{circ_cor} calculates the (sample) circular correlation between the columns of x.
#' Two parametric (the Jammalamadaka-Sarma (1988, equation 2.6) form \code{"js"}, and
#' the Fisher-Lee (1983, Section 3) form \code{"fl"})
#' and two non-parametric (two versions of Kendall's tau) correlation coefficients are considered.
#' The first version of Kendall's tau (\code{"tau1"}) is based on equation 2.1 in Fisher and Lee (1982),
#' whereas the second version (\code{"tau2"}) is computed using equations 6.7-6.8 in Zhan et al (2017).
#'
#' The cost-complexity for \code{"js"}, \code{"fl"}, \code{"tau2"} and \code{"tau1"} are \eqn{O(n), O(n^2), O(n^2)} and \eqn{O(n^3)}
#' respectively, where \eqn{n} denotes the number of rows in \code{x}. As such, for large \eqn{n} evaluation of
#' \code{"tau1"} will be slow.
#'
#'
#' @references
#'
#' Fisher, N. I. and Lee, A. J. (1982). Nonparametric measures of angular-angular association. Biometrika, 69(2), 315-321.
#'
#' Fisher, N. I. and Lee, A. J. (1983). A correlation coefficient for circular data. Biometrika, 70(2):327-332.
#'
#' Jammalamadaka, S. R. and Sarma, Y. (1988). A correlation coefficient for
#' angular variables. Statistical theory and data analysis II, pages 349-364.
#'
#' Zhan, X., Ma, T., Liu, S., & Shimizu, K. (2017). On circular correlation for data on the torus. Statistical Papers, 1-21.
#'
#'
#' @examples
#' # generate data from vmsin model
#' set.seed(1)
#' dat <- rvmsin(100, 2,3,-0.8,0,0)
#'
#' # now calculate circular correlation(s) between the 2 columns of dat
#' circ_cor(dat, type="js")
#' circ_cor(dat, type="fl")
#' circ_cor(dat, type="tau1")
#' circ_cor(dat, type="tau2")
#'
#'
#' @export


circ_cor <- function(x, type="js", alternative = "two.sided",
                     jackknife = FALSE, bootse = FALSE, n.boot = 100) {

  if (any(is.na(x)))
    stop("NA values in \'x\'")

  if (!is.matrix(x)) {
    stop("\'x\' must be a two-column matrix")
  }

  if (ncol(x)!=2) {
    stop("\'x\' must be a two-column matrix")
  }

  if (!type %in% c("fl", "js", "tau1", "tau2"))
    stop("\'type\' must be one of \'js\', \'fl\', \'tau1\' or \'tau2\'")

  x <- prncp_reg(x)
  n <- nrow(x)

  if (type == "fl") {
    calc_rho <- function(x) {

      rho <- calc_corr_fl(x)
      A_over_mu2 <- function(margin) {
        alpha <- sapply(1:2, function(p) sum(cos(p*margin))/n)
        beta <- sapply(1:2, function(p) sum(sin(p*margin))/n)
        A <- alpha[1]^2 + beta[1]^2 + alpha[2]*beta[1]^2 -
          alpha[1]^2*alpha[2] - 2*alpha[1]*beta[1]*beta[2]
        mu2 <- 0.5 * (1 - alpha[2]^2 - beta[2]^2)
        A/mu2
      }
      avar <- prod(apply(x, 2, A_over_mu2))

      se <- sqrt(avar)/sqrt(n)
      attr(rho, "se") <- se
      rho

    }
    # attr(rho_fl, "se") <- se
    # attr(rho_fl, "p.value") <- pval
    # rho
  } else if (type == "js") {
    calc_rho <- function(x) {
      sin_x_1_cent <- sin(x[, 1] - atan2(sum(sin(x[, 1])), sum(cos(x[, 1]))))
      sin_x_2_cent <- sin(x[, 2] - atan2(sum(sin(x[, 2])), sum(cos(x[, 2]))))
      num <- sum(sin_x_1_cent*sin_x_2_cent)
      den <- sqrt(sum(sin_x_1_cent^2)*sum(sin_x_2_cent^2))
      rho <- num/den

      # asymptotic variance
      # idx <- data.matrix(expand.grid(0:4, 0:4))
      idx <- rbind(
        c(2, 2),
        c(2, 0), c(0, 2),
        c(1, 3), c(3, 1),
        c(4, 0), c(0, 4)
      )
      rownames(idx) <- paste0(idx[, 1], idx[, 2])
      lambda <- apply(
        idx,
        1,
        function(ii) {
          sum(sin_x_1_cent^(ii[1]) * sin_x_2_cent^(ii[2]))/n
        }
      )
      avar <- unname(
        max(
          lambda["22"]/lambda["20"]*lambda["02"] -
            rho * (
              lambda["13"]/(lambda["20"]*sqrt(lambda["20"]*lambda["02"])) +
                lambda["31"]/(lambda["02"]*sqrt(lambda["20"]*lambda["02"]))
            ) +
            rho^2 / 4 * (
              1 +
                lambda["40"]/lambda["20"]^2 +
                lambda["04"]/lambda["02"]^2 +
                lambda["22"]/(lambda["20"]*lambda["02"])
            ),
          1e-10
        )
      )

      # browser()
      se <- sqrt(avar)/sqrt(n)
      # z <- rho_js/se
      attr(rho, "se") <- se
      # if (alternative == "two.sided") {
      #   pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
      # } else if (alternative == "less") {
      #   pval <- pnorm(z, lower.tail = TRUE)
      # } else if (alternative == "greater") {
      #   pval <- pnorm(z, lower.tail = FALSE)
      # }
      #
      # attr(rho_js, "se") <- se
      # attr(rho_js, "p.value") <- pval
      # rho_js
      rho
    }
  } else if (type == "tau1") {
    calc_rho <- calc_corr_tau_1
  } else {
    calc_rho <- calc_corr_tau_2
  }

  # browser()

  rho <- calc_rho(x)
  rho_attr <- attributes(rho)

  if (jackknife) {
    vals <- vapply(
      1:n,
      function(ii){
        c(calc_rho(x[-ii, , drop = FALSE]))
      },
      0
    )
    vals_adj <- n*rho - (n-1)*vals
    rho_attr$jackknife.est <- sum(vals_adj)/n
    rho_attr$jackknife.se <- sqrt(var(vals_adj)/n)
  }

  if (bootse) {
    boot_vals <- vapply(
      1:n.boot,
      function(ii) {
        idx <- sample(1:n, replace = TRUE)
        c(calc_rho(x[idx, , drop = FALSE]))
      },
      0
    )
    rho_attr$bootse <- sd(boot_vals)
  }

  if (type %in% c("js", "fl")) {
    z <- c(rho)/rho_attr$se
    if (alternative == "two.sided") {
      rho_attr$pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
    } else if (alternative == "less") {
      rho_attr$pval <- pnorm(z, lower.tail = TRUE)
    } else if (alternative == "greater") {
      rho_attr$pval <- pnorm(z, lower.tail = FALSE)
    }
  }

  attributes(rho) <- rho_attr
  rho

}


#' Analytic circular variances and correlations for bivariate angular models
#' @param model bivariate angular model. Must be one of \code{"vmsin"},
#' \code{"vmcos"}, or \code{"wnorm2"}.
#' @param kappa1,kappa2,kappa3 concentration and covariance parameters.
#' Recycled to the same size. kappa3^2 must be < kappa1*kappa2 in the wnorm2 model
#' (see \link{rwnorm2} for a detailed parameterization of \code{wnorm2}).
#' @param mu1,mu2 mean parameters. Ignored as they do not play any role in
#' the analytical formulas.
#' @param nsim Monte Carlo sample size. Ignored if all of \code{kappa1}, \code{kappa2}
#' and \code{abs(kappa3)} are < 150 or if model = \code{"wnorm2"}.
#' @inheritParams contour_model
#'
#' @return
#' Returns a list with elements \code{var1}, \code{var2} (circular variances for the
#' first and second coordinates), \code{rho_fl} and \code{rho_js} (circular correlations).
#' See details.
#'
#' @details
#' The function computes the analytic circular variances and correlations
#' (both Jammalamadaka-Sarma (JS) and Fisher-Lee (FL) forms) for von Mises sine,
#' von Mises cosine and bivariate wrapped normal distributions.
#'
#' For \code{wnorm2}, expressions for the circular variances,
#' JS and FL correlation coefficients can be found in Mardia and Jupp (2009),
#' Jammalamadaka and Sarma (1988) and Fisher and Lee (1983) respectively.
#' For \code{vmsin} and \code{vmcos} these expressions are provided in Chakraborty and Wong (2018).
#'
#' Because the analytic expressions in \code{vmsin} and \code{vmcos} models involve infinite sums
#' of product of Bessel functions,
#' if any of \code{kappa1}, \code{kappa2} and \code{abs(kappa3)} is larger
#' than or equal to 150, IID Monte Carlo with sample size \code{nsim} is used
#' to approximate \code{rho_js} for numerical stability.  From \code{rho_js},
#' \code{rho_fl} is computed using Corollary 2.2 in
#' Chakraborty and Wong (2018), which makes cost-complexity for
#' the \code{rho_fl} evaluation to be of order  O(\code{nsim}) for \code{vmsin}
#' and \code{vmcos} models. (In general,  \code{rho_fl} evaluation
#' is of order O(\code{nsim}^2)).
#'
#' In addition, for the \code{vmcos} model, when \code{-150 < kappa3 < -1}
#' or \code{50 < max(kappa1, kappa2, abs(kappa3)) <= 150}, the analytic formulas
#' in Chakraborty and Wong (2018) are used; however, the reciprocal of the normalizing
#' constant and its partial derivatives are all calculated numerically via (quasi) Monte carlo method for
#' numerical stability. These (quasi) random numbers can be provided through the
#' argument \code{qrnd}, which must be a two column matrix, with each element being
#' a  (quasi) random number between 0 and 1. Alternatively, if \code{n_qrnd} is
#' provided (and \code{qrnd} is missing), a two dimensional sobol sequence of size \code{n_qrnd} is
#' generated via the function \link{sobol} from the R package \code{qrng}. If none of \code{qrnd}
#' or \code{n_qrnd} is available, a two dimensional sobol sequence of size 1e4 is used.
#'
#'
#' @examples
#' circ_varcor_model("vmsin", kappa1= 1, kappa2 = 2, kappa3 = 3)
#'
#' # Monte Carlo approximation
#' set.seed(1)
#' dat <- rvmsin(1000, 1, 2, 3)
#' # sample circular variance
#' circ_var <- function(x)
#'   1 - mean(cos(x - atan2(mean(sin(x)), mean(cos(x))) ))
#' circ_var(dat[, 1])
#' circ_var(dat[, 2])
#' circ_cor(dat, "fl")
#' circ_cor(dat, "js")
#'
#' @references
#' Fisher, N. I. and Lee, A. (1983). A correlation coefficient for circular data. Biometrika, 70(2):327-332.
#'
#' Jammalamadaka, S. R. and Sarma, Y. (1988). A correlation coefficient for
#' angular variables. Statistical theory and data analysis II, pages 349-364.
#'
#' Mardia, K. and Jupp, P. (2009). Directional Statistics. Wiley Series in Probability and Statistics. Wiley.
#'
#' Chakraborty, S. and Wong, S, W.K. (2018). On the circular correlation coefficients
#' for bivariate von Mises distributions on a torus. arXiv e-print.
#'
#' @export

circ_varcor_model <- function(model = "vmsin", kappa1 = 1, kappa2 = 1, kappa3 = 0,
                              mu1 = 0, mu2 = 0, nsim = 1e4, ...)
{

  if(any(c(kappa1, kappa2) < 0))
    stop("kappa1 and kappa2 must be non-negative")

  if (length(model) != 1 | !model %in% c("vmsin", "vmcos", "wnorm2"))
    stop("\'model\' must be one of \"vmsin\", \"vmcos\" or \"wnorm2\"")

  if (nsim <= 0)
    stop("\'nsim\' must be a positive integer")

  ell <- list(...)
  nsim <- round(nsim)

  if (model == "vmcos") {


    if (!is.null(ell$qrnd)) {
      qrnd_grid <- ell$qrnd
      dim_qrnd <- dim(qrnd_grid)
      if (!is.matrix(qrnd_grid) | is.null(dim_qrnd) |
          dim_qrnd[2] != 2)
        stop("\'qrnd\' must be a two column matrix")
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
  }

  if(max(length(kappa1), length(kappa2), length(kappa3), length(mu1), length(mu2)) > 1) {
    expanded <- expand_args(kappa1, kappa2, kappa3, mu1, mu2)
    kappa1 <- expanded[[1]]
    kappa2 <- expanded[[2]]
    kappa3 <- expanded[[3]]
    if (model == "wnorm2" &
        any (kappa1*kappa2 - kappa3*kappa3 <= 1e-10))
      stop("abs(kappa3) must be less than sqrt(kappa1*kappa2) in wnorm2")
    lapply(1:length(kappa1),
           function(j) {
             inargs <- list(kappa1 = kappa1[j], kappa2 = kappa2[j],
                            kappa3 = kappa3[j])
             if (model == "vmcos") {
               inargs$qrnd_grid <- qrnd_grid
               # inargs$force_approx_const <- ell$force_approx_const
             }
             do.call(paste0(model, "_var_cor_singlepar"),
                     inargs)
           }
    )
  } else {
    if (model == "wnorm2" &
        (kappa1*kappa2 - kappa3*kappa3 <= 1e-10))
      stop("abs(kappa3) must be less than sqrt(kappa1*kappa2) in wnorm2")

    inargs <- list(kappa1 = kappa1, kappa2 = kappa2,
                   kappa3 = kappa3)
    if (model == "vmcos") inargs$qrnd_grid <- qrnd_grid
    do.call(paste0(model, "_var_cor_singlepar"),
            inargs)
  }

}

