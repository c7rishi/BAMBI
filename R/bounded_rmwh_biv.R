bounded_rwmh_biv <- function(lpr,
                             init,
                             lower,
                             upper,
                             dep_cov = FALSE,
                             dep_cov_type = NULL,
                             zero_cov = FALSE,
                             step)
{

  lpr.init <- lpr(init)
  ncomp <- ncol(init)
  prop <- init

  if (zero_cov) prop[3, ] <- 0

  prop[1:2, ] <-  pmax(init[1:2, ] + rnorm(2*ncomp, 0, step[1:2]), 1e-6)

  prop[4:5, ] <- prncp_reg(init[4:5, ] + rnorm(2*ncomp, 0, step[4:5]))

  if (!zero_cov) {
    prop[3, ] <-  init[3, ] + rnorm(ncomp, 0, step[3])

    if (dep_cov) {
      if (dep_cov_type %in% c("wnorm2_bound", "vmsin_unimodal")) {
        bd_k1k2 <- sqrt(prop[1, ]*prop[2, ])
        lower[3, ] <- pmax(-bd_k1k2, lower[3, ])
        upper[3, ] <- pmin(bd_k1k2, upper[3, ])
      } else {
        # dep_cov_type == "vmcos_unimodal"
        lower[3, ] <- pmax(-prop[1, ]*prop[2, ]/
                             (prop[1, ]+prop[2, ]),
                           lower[3, ])
      }
    }
  }

  bd_err <- any(c(prop-lower, upper-prop) < 0)

  if (bd_err) {
    lpr.prop <- -Inf
  } else {
    lpr.prop <- lpr(prop)
  }

  aprob <- min(1, exp(lpr.prop-lpr.init))

  accpt <- (runif(1) < aprob)

  # if (is.na(propcheck)) browser()

  if (accpt) {
    final <- prop
    lpr.final <- lpr.prop
  } else {
    final <- init
    lpr.final <- lpr.init
  }


  out <- list (final=final, lpr=lpr.final, step=step,
               aprob=aprob, accpt=accpt*1)

  out
}


bounded_rwmh_uni <- function(lpr,
                             init,
                             lower,
                             upper,
                             step)
{

  lpr.init <- lpr(init)
  ncomp <- ncol(init)
  prop <- init

  prop[1, ] <-  pmax(init[1, ] + rnorm(ncomp, 0, step[1]), 1e-6)
  prop[2, ] <- prncp_reg(init[2, ] + rnorm(ncomp, 0, step[2]))

  bd_err <- any(c(prop-lower, upper-prop) < 0)

  if (bd_err) {
    lpr.prop <- -Inf
  } else {
    lpr.prop <- lpr(prop)
  }

  aprob <- min(1, exp(lpr.prop-lpr.init))

  accpt <- (runif(1) < aprob)

  # if (is.na(propcheck)) browser()

  if (accpt) {
    final <- prop
    lpr.final <- lpr.prop
  } else {
    final <- init
    lpr.final <- lpr.init
  }


  out <- list (final=final, lpr=lpr.final, step=step,
               aprob=aprob, accpt=accpt*1)

  out
}
