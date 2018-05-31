bounded_hmc_biv <- function(lpr_grad,
                            init,
                            lower,
                            upper,
                            dep_cov = FALSE,
                            dep_cov_type = NULL,
                            zero_cov = FALSE,
                            nsteps = 1,
                            step)
{

  broken <- FALSE # will be true if breaks
  apr <- NULL #will  be replaced if doesn't break
  delta <- NULL #will be replaced if doesn't break

  lower1 <- lower
  upper1 <- upper

  ncomp <- 1
  npar <- 5
  kappa3_index <- 3
  # positions of kappa3 in vec(par_mat)

  lpr_grad.init <- lpr_grad(init)
  lpr.init <-   lpr_grad.init$lpr
  gr <- lpr_grad.init$grad

  init.p <- matrix(rnorm(npar), 5, ncomp)

  # Compute the kinetic energy at the start of the trajectory.
  kinetic.init <- sum(init.p^2) / 2

  # Compute the trajectory by the leapfrog method.
  q <- init
  p <- init.p

  if (zero_cov) {
    p[kappa3_index] <- 0
  }


  reflections <- 0

  # Make a half step for momentum at the beginning.

  p <- p + (step/2) * gr

  if (zero_cov) {
    p[kappa3_index] <- 0
  }

  for (i in 1:nsteps)
  {
    # Make a full step for the position.

    q <- q + step * p

    # check if broken
    if (any(is.nan(c(p, q)))) {
      broken <- TRUE
      #stop("Algorithm breaks. Try a smaller epsilon.")
      break
    }


    # Check for bound violations, and adjust position and momentum
    for(k in 1:npar) {
      # adjust the bound for the cov
      # param if dependent (for wnorm2)
      if(k == 3) {

        if (zero_cov) {
          next
        }
        else if (dep_cov) {

          if (any(q[1:2] <= 0)) {
            broken <- TRUE
            break
          }

          if (dep_cov_type %in% c("wnorm2_bound", "vmsin_unimodal")) {
            bd_k1k2 <- sqrt(q[1]*q[2])
            lower[k] <- max(lower1[k], -bd_k1k2)
            upper[k] <- min(upper1[k], bd_k1k2)
          }
          else {
            # dep_cov_type == "vmcos_unimodal"
            lower[k] <- max(lower1[k], -q[1]*q[2]/(q[1]+q[2]))
          }

        }

        # if (any(is.nan(c(upper[k], lower[k])))) {
        #   broken <- TRUE
        #   break
        # }
      }

      reflections_curr <- 0L

      while(all(q[k] < lower[k] || q[k] > upper[k],
                reflections_curr <= 100)) {
        if (q[k]<lower[k]) {
          q[k] <- lower[k] + (lower[k] - q[k])
          p[k] <- -p[k]
          reflections_curr <- reflections_curr + 1L
        }

        else if (q[k]>upper[k]) {
          q[k] <- upper[k] - (q[k] - upper[k])
          p[k] <- -p[k]
          reflections_curr <- reflections_curr + 1L
        }

      }

      reflections <- reflections + reflections_curr

      if (reflections_curr >= 100) {
        broken <- TRUE
        break
      }

    }

    # check for broken, bound violation, or nan
    # if so then break from i loop
    if (any(is.nan(c(p, q)),
            broken)) {
      broken <- TRUE
      #stop("Algorithm breaks. Try a smaller epsilon.")
      break
    }



    # Evaluate the gradient at the new position, provided not broken

    lpr_grad_current <- lpr_grad(q)
    lr <- lpr_grad_current$lpr
    gr <- lpr_grad_current$grad


    # Make a full step for the momentum, except when we're coming to the end of
    # the trajectory.

    if (i != nsteps) {
      p <- p + step * gr
      if (zero_cov) {
        p[kappa3_index] <- 0
      }
    }
  }


  # Make a half step for momentum at the end.

  p <- p + (step/2) * gr
  if (zero_cov) {
    p[kappa3_index] <- 0
  }

  # Negate momentum at end of trajectory to make the proposal symmetric.

  p <- -p


  # Look at log probability and kinetic energy at the end of the trajectory.
  lpr.prop <- lr
  kinetic.prop <- sum(p^2) / 2

  # Accept or reject the state at the end of the trajectory.
  H.init <- -lpr.init + kinetic.init
  H.prop <- -lpr.prop + kinetic.prop
  delta <- H.prop - H.init
  apr <- min(1,exp(-delta))


  if (any(is.nan(c(p, q, apr)))) {
    broken <- TRUE
    apr <- NULL
    delta <- NULL
    #stop("Algorithm breaks. Try a smaller epsilon.")
  }


  if(broken) # reject
  {
    final.q <- init
    final.p <- init.p
    lpr.final <- lpr.init
    acc <- 0
  }
  else if (runif(1) > apr) # reject
  {
    final.q <- init
    final.p <- init.p
    lpr.final <- lpr.init
    acc <- 0
  }
  else  # accept
  {
    final.q <- q
    final.p <- p
    lpr.final <- lr
    acc <- 1
  }

  # Return new state, its log probability and gradient, plus additional
  # information, including the trajectory, if requested.

  out <- list (final=final.q, final.p=final.p, lpr=lpr.final, step=step,
               apr=apr, accpt=acc, delta=delta, reflections=reflections, broken = broken)

  out
}


bounded_hmc_uni <- function(lpr_grad,
                            init,
                            lower,
                            upper,
                            nsteps = 1,
                            step)
{

  broken <- FALSE # will be true if breaks
  apr <- NULL #will  be replaced if doesn't break
  delta <- NULL #will be replaced if doesn't break

  ncomp <- ncol(init)
  npar <- 2*ncomp

  lpr_grad.init <- lpr_grad(init)
  lpr.init <-   lpr_grad.init$lpr
  gr <- lpr_grad.init$grad

  init.p <- matrix(rnorm(npar), 2, ncomp)

  # Compute the kinetic energy at the start of the trajectory.
  kinetic.init <- sum(init.p^2) / 2

  # Compute the trajectory by the leapfrog method.
  q <- init
  p <- init.p


  reflections <- 0

  # Make a half step for momentum at the beginning.

  p <- p + (step/2) * gr

  # browser()
  # Alternate full steps for position and momentum.

  for (i in 1:nsteps)
  {
    # Make a full step for the position.

    q <- q + step * p

    # check if broken
    if (any(is.nan(c(p, q)))) {
      broken <- TRUE
      #stop("Algorithm breaks. Try a smaller epsilon.")
      break
    }


    # Check for bound violations, and adjust position and momentum
    for(k in 1:npar) {
      reflections_curr <- 0L

      while(all(q[k] < lower[k] || q[k] > upper[k],
                reflections_curr <= 100)) {
        if (q[k]<lower[k]) {
          q[k] <- lower[k] + (lower[k] - q[k])
          p[k] <- -p[k]
          reflections_curr <- reflections_curr + 1L
        }

        else if (q[k]>upper[k]) {
          q[k] <- upper[k] - (q[k] - upper[k])
          p[k] <- -p[k]
          reflections_curr <- reflections_curr + 1L
        }

      }

      reflections <- reflections + reflections_curr

      if (reflections_curr >= 100) {
        broken <- TRUE
        break
      }
    }

    # check if broken
    if (any(is.nan(c(p, q)), broken)) {
      broken <- TRUE
      #stop("Algorithm breaks. Try a smaller epsilon.")
      break
    }


    # Evaluate the gradient at the new position.

    lpr_grad_current <- lpr_grad(q)
    lr <- lpr_grad_current$lpr
    gr <- lpr_grad_current$grad


    # Make a full step for the momentum, except when we're coming to the end of
    # the trajectory.

    if (i != nsteps) {
      p <- p + step * gr
    }
  }


  # Make a half step for momentum at the end.

  p <- p + (step/2) * gr

  # Negate momentum at end of trajectory to make the proposal symmetric.

  p <- -p


  # Look at log probability and kinetic energy at the end of the trajectory.
  lpr.prop <- lr
  kinetic.prop <- sum(p^2) / 2

  # Accept or reject the state at the end of the trajectory.
  H.init <- -lpr.init + kinetic.init
  H.prop <- -lpr.prop + kinetic.prop
  delta <- H.prop - H.init
  apr <- min(1,exp(-delta))


  if (any(is.nan(c(p, q, apr)))) {
    broken <- TRUE
    apr <- NULL
    delta <- NULL
    #stop("Algorithm breaks. Try a smaller epsilon.")
  }


  if(broken) # reject
  {
    final.q <- init
    final.p <- init.p
    lpr.final <- lpr.init
    acc <- 0
  }
  else if (runif(1) > apr) # reject
  {
    final.q <- init
    final.p <- init.p
    lpr.final <- lpr.init
    acc <- 0
  }
  else  # accept
  {
    final.q <- q
    final.p <- p
    lpr.final <- lr
    acc <- 1
  }

  # Return new state, its log probability and gradient, plus additional
  # information, including the trajectory, if requested.

  out <- list (final=final.q, final.p=final.p, lpr=lpr.final, step=step,
               apr=apr, accpt=acc, delta=delta, reflections=reflections, broken = broken)

  out
}
