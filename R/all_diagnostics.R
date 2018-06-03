#' Contour plot for angmcmc objects with bivariate data
#'
#' @inheritParams pointest
#' @param x angular MCMC object (with bivariate data).
#' @param show.data logical. Should the data points be added to the contour plot? Ignored if \code{object} is NOT supplied.
#' @param cex,col graphical parameters passed to \code{\link{points}} in graphics for plotting the data points.
#' Ignored if {show.data == FALSE}.
#' @inheritParams contour_model
#' @param ... additional arguments to be passed to the function \code{\link{contour}}.
#'
#' @details
#' \code{contour.angmcmc} is an S3 function for angmcmc objects that calls \code{\link{contour}} from graphics.
#'
#' To estimate the mixture density required to construct the contour plot, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples, yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}. Then the mixture density
#' \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by \eqn{f(x|\hat{\eta})}.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # now create a contour plot
#' contour(fit.vmsin.20)
#'
#' @export

contour.angmcmc <-  function(x, fn = mean, show.data = TRUE,
                             xpoints = seq(0, 2*pi, length.out = 100),
                             ypoints = seq(0, 2*pi, length.out = 100),
                             levels, nlevels = 20,
                             cex = 1, col = "red", ...)
{
  object <- x

  if (!is.angmcmc(object))
    stop("\'x\' must be an angmcmc object")

  if(x$type != "bi") stop("\"x\" is not a bivariate angmcmc object")

  if(missing(levels)) {
    levels <- exp(seq(-20,2, length.out = nlevels))
  }

  colnames_data <- colnames(x$data)

  if(is.null(colnames_data)) {
    xlab <- ylab <- ""
  } else {
    xlab <- colnames_data[1]
    ylab <- colnames_data[2]
  }

  if(x$ncomp > 1) {
    main <- paste("contour plot for fitted", x$ncomp, "component", x$model, "mixtures")
  } else {
    main <- paste("contour plot for fitted (single component)", x$model)
  }

  coords <- as.matrix(expand.grid(xpoints, ypoints))
  dens <- d_fitted(coords, x, fn = fn)
  contour(xpoints, ypoints, matrix(dens, nrow=length(xpoints)), levels=levels)

  if(show.data) points(x$data, col = col, cex = cex)

  title(main = main, xlab = xlab, ylab = ylab)
}


#
# panel_wireframe_cloud <- function(x, y, z, x2, y2, z2,...) {
#   panel.wireframe(x, y, z,...)
#   panel.cloud(x2, y2, z2,...)
# }


#' Density plots for angmcmc objects
#'
#' @description Plot fitted angular mixture model density surfaces or curves.
#' @inheritParams pointest
#' @param x angmcmc object.
#' @param plot logical. Should the density surface (if the fitted data is bivariate) or the density
#' curve (if univariate) be plotted?
#' @param log.density logical. Should log density be used for the plot?
#' @param ... additional arguments passed to \code{\link{persp}}, if
#' fitted data is bivariate, or to \link{hist} (if (\code{show.hist == TRUE})), if the fitted data is univariate
#' @param show.hist logical. Should a histogram for the data
#' points be added to the plot, if the fitted data is univariate? Ignored if data is
#' bivariate.
#' @param theta,phi,shade,expand,xlab,ylab,zlab,main grahpical parameters passed to \link{persp} (if
#' bivariate) or \link{plot} (if univariate). If the data is univariate,
#' \code{theta, phi, shade, expand} are ignored, and \code{zlab} plays the role of
#' \code{ylab}.
#' @param xpoints,ypoints Points on the  x and y coordinates (if bivariate) or only x coordinate
#' (if univariate) where the density is to be evaluated. Each defaults to seq(0, 2*pi, length.out=100).
#'
#' @details
#' When \code{plot==TRUE}, \code{densityplot.angmcmc} calls \link{persp} or
#' \link{plot} from graphics to draw the surface or curve.
#'
#' To estimate the mixture density, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples, yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}. Then the mixture density
#' \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by \eqn{f(x|\hat{\eta})}.
#'
#'
#' Note that \code{densityplot.angmcmc} \strong{does not} plot the kernel densitie estimates
#' of the MCMC parameters. (These plots can be obtained by first converting an \code{angmcmc}
#' object to an \code{mcmc} object via \link{as.mcmc.list}, and then
#' by using \code{densplot} from package coda on the resulting \code{mcmc.list} object. Instead,
#' \code{densityplot.angmcmc} returns the surface (if 2-D) or the curve (if 1-D)
#' of the fitted model density evaluated at the estimated parameter vector (obtain through \link{pointest}).
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # now create density surface with the default first 1/3 as burn-in and thin = 1
#' library(lattice)
#' densityplot(fit.vmsin.20)
#' # the viewing angles can be changed through the arguments theta and phi
#' # (passed to persp from graphics)
#' densityplot(fit.vmsin.20, theta = 45, phi = 45)
#'
#' # Now fit a vm mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vm.20 <- fit_vmmix(wind$angle, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' densityplot(fit.vm.20)
#'
#' @importFrom lattice densityplot
#'
#' @export

densityplot.angmcmc <- function(x, fn = mean, log.density = FALSE,
                                xpoints=seq(0, 2*pi, length.out=35),
                                ypoints=seq(0, 2*pi, length.out=35),
                                plot=TRUE,
                                show.hist=ifelse(log.density, FALSE, TRUE),
                                theta = 30, phi = 30, shade = 0.01,
                                expand = 0.5, xlab, ylab,
                                zlab = ifelse(log.density, "Log Density", "Density"),
                                main,
                                ...)
{
  if(!is.angmcmc(x))
    stop("\"x\" must be an angmcmc object")

  object <- x

  if (object$type == "bi") {

    colnames_data <- colnames(object$data)

    if (any(missing(xlab), missing(ylab))) {
      if (is.null(colnames_data)) {
        xlab <- ylab <- ""
      } else {
        xlab <- colnames_data[1]
        ylab <- colnames_data[2]
      }
    }

    if (missing(main)) {
      if(x$ncomp > 1) {
        main <- paste("contour plot for fitted", x$ncomp, "component", x$model, "mixtures")
      } else {
        main <- paste("contour plot for fitted (single component)", x$model)
      }
    }

    coords <- as.matrix(expand.grid(xpoints, ypoints))
    den <- d_fitted(coords, object, fn = fn)


    denmat <- matrix(den, nrow=length(xpoints))

    if(log.density) {
      denmat <- log(denmat)
    }


    out <- list(x=xpoints, y=ypoints, density=denmat)

    if (plot) {
      nrden <- nrow(denmat)
      ncden <- ncol(denmat)

      if(object$ncomp > 1) {
        main <- paste("Density surface for fitted", object$ncomp, "component",
                      object$model, "mixtures")
      } else {
        main <- paste("Density surface for fitted (single component)", object$model)
      }

      # Create a function interpolating colors in the range of specified colors
      jet.colors <- grDevices::colorRampPalette( c("blue", "green",
                                                   "yellow", "orange", "red") )
      # Generate the desired number of colors from this palette
      nbcol <- 500
      color <- jet.colors(nbcol)

      denfacet <- denmat[-1, -1] + denmat[-1, -ncden] +
        denmat[-nrden, -1] + denmat[-nrden, -ncden]
      # Recode facet z-values into color indices
      facetcol <- cut(denfacet, nbcol)

      persp(x=xpoints, y=ypoints, z=denmat, theta = theta, phi = phi, expand = expand, col = color[facetcol],
            ltheta = 120, shade = shade, ticktype = "detailed",
            xlab = xlab, ylab = ylab, zlab = zlab,
            main = main, ...) -> res

      # inargs <- list(...)
      # inargs$x <- denmat~x*y
      # inargs$data <- data.frame(x = xpoints,
      #                           y=rep(ypoints, each=length(xpoints)),
      #                           denmat=denmat)
      # inargs$outerbox <- FALSE
      # inargs$par.settings <- list(axis.line = list(col = 'transparent'))
      # if (is.null(inargs$xlab)) inargs$xlab <- xlab
      # if (is.null(inargs$ylab)) inargs$ylab <- ylab
      # if (is.null(inargs$colorkey)) inargs$colorkey <- FALSE
      # if (is.null(inargs$main)) inargs$main <- main
      # if (is.null(inargs$neval)) inargs$neval <- 100
      # if (is.null(inargs$aspect)) inargs$aspect <- c(61/87, 0.4)
      # if (is.null(inargs$zlab)) inargs$zlab <- list("Density", rot=90)
      # if (is.null(inargs$screen)) inargs$screen <- list(z=45, x=-45)
      # if (is.null(inargs$colorkey)) inargs$colorkey <- FALSE
      # if (is.null(inargs$scales))
      #   inargs$scales <- list(arrows=FALSE, col=1)
      # if (is.null(inargs$drape)) inargs$drape <- TRUE
      # if (is.null(inargs$light.source))
      #   inargs$light.source <- c(10,0,10)
      # if (is.null(inargs$col.regions))
      #   inargs$col.regions <- colorRampPalette(c("blue", "green",
      #                                            "yellow", "orange", "red"))(100)
      # if (is.null(inargs$par.settings))
      #   inargs$par.settings <- list(top.padding = 0,
      #                               bottom.padding = 0,
      #                               left.padding=0,
      #                               right.padding=0,
      #                               axis.line=list(col = 'transparent'))
      # do.call(wireframe, inargs)
    }

    invisible(out)
  }
  else {

    den <- d_fitted(xpoints, object, fn = fn)

    if (log.density) den <- log(den)

    out <- list(x=xpoints, density=den)

    if (plot) {
      if(show.hist){
        histplot <-  hist(object$data, plot = FALSE, ...)
      } else {
        histplot <- NULL
      }

      if (missing(main)) {
        if(object$ncomp > 1) {
          main <- paste("density plot for fitted", object$ncomp, "component", object$model, "mixtures")
        } else {
          main <- paste("density plot for fitted (single component)", object$model)
        }
      }

      if (missing(xlab))
        xlab <- "Angles in radian"

      y_max <- 1.1* max(den, histplot$density)
      plot(NULL, xlim=range(xpoints), ylim=c(0, y_max), xlab = xlab,
           ylab=zlab, main=main)
      points(xpoints, den, type = "l")

      if(show.hist) plot(histplot, freq = FALSE, add = TRUE, ...)


      title(main = main)
    }
  }
  invisible(out)
}



#' Trace plot for parameters from an angmcmc object
#' @inheritParams pointest
#' @param object angular MCMC object.
#' @param par parameter for which trace plot is to be created.
#' @param press.enter logical. Should the next plot in the series
#' be shown after you press "Enter"? Ignored if only a single plot
#' is to be created.
#' @param ... unused
#' @return
#' Returns a single plot if a single \code{par} and a single \code{comp.label} is supplied.
#' Otherwise, a series of plots is produced.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # trace plot for kappa1 in component 1
#' paramtrace(fit.vmsin.20, "kappa1", 1)
#' # for kappa1 in all components
#' paramtrace(fit.vmsin.20, "kappa1", press.enter = FALSE)
#' # for all parameters in component 1
#' paramtrace(fit.vmsin.20, comp.label = 1, press.enter = FALSE)
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @export

paramtrace <- function(object, par.name, comp.label, chain.no,
                       press.enter = TRUE, ...)
{
  if(!class(object) %in% "angmcmc") stop("\'object\' must be an angmcmc object")

  ell <- list(...)

  if (!is.null(ell$burnin))
    warning("Use of burnin is depriciated in postprocessing. Use \'burnin.prop\' during original MCMC run instead.")
  if (!is.null(ell$thin))
    warning("Use of thin is depriciated in postprocessing. Use \'thin\' during original MCMC run instead.")

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


  if(any(c(length(par.name), length(comp.label)) > 1)) {
    singleplot <- FALSE
  } else {
    singleplot <- TRUE
  }

  n.chains <- length(chain.no)
  col_brew <- brewer.pal(n = max(n.chains, 3), name = "Pastel2")

  if(singleplot) {
    val <- extractsamples(object, par.name, comp.label, chain.no, drop = FALSE)
    plot(NULL, xlim = c(1, object$n.iter.final), ylim = range(val),  ylab="", xlab = "Iteration")
    for(ch in 1:n.chains)
      points(val[, , , ch], type = "l", col = col_brew[ch])
    legend("bottomright", legend = paste("Chain", chain.no), col = col_brew,
           lty = 1)
    if(object$ncomp > 1) {
      ylab <- paste(par.name, "for component", comp.label)
      main <- paste("Traceplot for", ylab, "in", object$ncomp, "component", object$model, "mixtures")
    } else {
      ylab <- par.name
      main <- paste("Traceplot for", ylab, "in (single component)", object$model)
    }
    title(main = main, ylab = ylab)
  }

  # not singleplot
  else {
    nplots <- length(par.name) * length(comp.label)
    currplotno <- 1L
    for(par.curr in par.name) {
      for(comp.label.curr in comp.label) {
        val <- extractsamples(object, par.curr, comp.label.curr, chain.no, drop = FALSE)
        plot(NULL, xlim = c(1, object$n.iter.final), ylim = range(val),  ylab="", xlab = "Iteration")
        for(ch in 1:n.chains)
          points(val[, , , ch], type = "l", col = col_brew[ch])
        legend("bottomright", legend = paste("Chain", chain.no), col = col_brew,
               lty = 1)
        if(object$ncomp > 1) {
          ylab <- paste(par.curr, "for component", comp.label.curr)
          main <- paste("Traceplot for", ylab, "in ", object$ncomp, "component", object$model, "mixtures")
        } else {
          ylab <- par.curr
          main <- paste("Traceplot for", ylab, "in  (single component)", object$model)
        }

        title(main = main, ylab = ylab)
        if(currplotno < nplots) {
          if(press.enter) {
            press_enter()
            frame()
          }
        }
        currplotno <- currplotno + 1
      }
    }
  }
}


#' Trace and autocorrelation plots of log posterior density or log likelihood from an angmcmc object
#' @inheritParams paramtrace
#' @param object angular MCMC object.
#' @param use.llik logical. Should log likelihood be plotted instead of log posterior? Set
#' to \code{FALSE} by default.
#' @param plot.autocor logical. Should the autocorrelations be plotted as well?
#' @param lag.max maximum lag for autocorrelation.  Passed to \link{acf}. Ignored if
#' \code{plot.autocor = FALSE}.
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' # log posterior density trace
#' lpdtrace(fit.vmsin.20)
#' # log likelihood trace
#' lpdtrace(fit.vmsin.20, use.llik = TRUE)
#'
#' @export

lpdtrace <- function(object, chain.no, use.llik = FALSE,
                     plot.autocor = FALSE, press.enter = TRUE,
                     lag.max = NULL, ...)
{
  if(class(object) != "angmcmc") stop("lpdtrace can only be used for \'angmcmc\' objects")

  ell <- list(...)

  if (!is.null(ell$burnin))
    warning("Use of burnin is depriciated in postprocessing. Use \'burnin.prop\' during original MCMC run instead.")
  if (!is.null(ell$thin))
    warning("Use of thin is depriciated in postprocessing. Use \'thin\' during original MCMC run instead.")


  if (missing(chain.no)) {
    chain.no <- 1:object$n.chains
  } else if (any(!chain.no %in% 1:object$n.chains)) {
    stop("invalid chain number")
  }


  n.chains <- length(chain.no)
  final_iter_set <- object$final_iter
  col_brew <- brewer.pal(n = max(n.chains, 3), name = "Paired")

  # browser()

  if (use.llik) {
    val <- object$llik[final_iter_set, chain.no, drop = FALSE]

    plot(NULL, xlim = c(1, object$n.iter.final), ylim = range(val),
         ylab="Log Likelihood", xlab = "Iteration")
    for(ch in 1:n.chains)
      points(val[, ch], type = "l", col = col_brew[ch])
    legend("bottomright", legend = paste("Chain", chain.no), col = col_brew,
           lty = 1)

    if(object$ncomp > 1) {
      main <- paste("Log likelihood traceplot for ", object$ncomp, "component", object$model, "mixtures")
    } else {
      main <- paste("Log likelihood traceplot for (single component)", object$model)
    }

    title(main = main)



    if (press.enter & plot.autocor) {
      press_enter()
      frame()
    }

    if (plot.autocor) {
      all_autocors <- lapply(1:n.chains,
                             function(ch) acf(val[, ch], lag.max = lag.max, plot = FALSE))

      range_y <- range(unlist(lapply(all_autocors, function(j) j$acf)))
      for(ch in 1:n.chains) {
        acr <- all_autocors[[ch]]
        if (ch  == 1) {
          plot(acr$acf, type = "b", col = col_brew[ch],
               ylab="Log likelihood autocorrelation",
               xlab = "Lag", pch = 16, ylim = range_y)
        } else {
          points(acr$acf, type = "b", col = col_brew[ch], pch = 16)
        }
      }

      abline(h = 0, lty = 2)
      legend("topright", legend = paste("Chain", chain.no), col = col_brew,
             pch = 16, lty = 1)

      if(object$ncomp > 1) {
        main <- paste("Log likelihood autocorrelation plot for ", object$ncomp, "component", object$model, "mixtures")
      } else {
        main <- paste("Log likelihood autocorrelation plot for (single component)", object$model)
      }

      title(main = main)

    }



  } else {

    val <- object$lpd[final_iter_set, chain.no, drop = FALSE]

    plot(NULL, xlim = c(1, object$n.iter.final), ylim = range(val),
         ylab="Log Posterior Density", xlab = "Iteration")
    for(ch in 1:n.chains)
      points(val[, ch], type = "l", col = col_brew[ch])
    legend("bottomright", legend = paste("Chain", chain.no), col = col_brew,
           lty = 1)

    if(object$ncomp > 1) {
      main <- paste("Log Posterior Density traceplot for", object$ncomp, "component", object$model, "mixtures")
    } else {
      main <- paste("Log Posterior Density traceplot fitted (single component)", object$model)
    }

    title(main = main)


    if (press.enter & plot.autocor) {
      press_enter()
      frame()

    }

    if (plot.autocor) {

      all_autocors <- lapply(1:n.chains,
                             function(ch) acf(val[, ch], lag.max = lag.max, plot = FALSE))

      range_y <- range(unlist(lapply(all_autocors, function(j) j$acf)))
      for(ch in 1:n.chains) {
        acr <- all_autocors[[ch]]
        if (ch  == 1) {
          plot(acr$acf, type = "b", col = col_brew[ch],
               ylab="Log posterior autocorrelation",
               xlab = "Lag", pch = 16, ylim = range_y)
        } else {
          points(acr$acf, type = "b", col = col_brew[ch], pch = 16)
        }
      }

      abline(h = 0, lty = 2)
      legend("topright", legend = paste("Chain", chain.no), col = col_brew,
             pch = 16, lty = 1)

      if(object$ncomp > 1) {
        main <- paste("Log posterior autocorrelation plot for ", object$ncomp, "component", object$model, "mixtures")
      } else {
        main <- paste("Log posterior autocorrelation plot for (single component)", object$model)
      }

      title(main = main)

    }

  }



}


#' Summary plots for angmcmc objects
#' @inheritParams paramtrace
#' @inheritParams lpdtrace
#' @param do.paramtrace logical. Should the trace(s) for the
#' parameter(s) be plotted?
#' @param do.lpdtrace logical. Should the log posterior trace
#' be plotted?
#' @param use.llik logical. Should the log likelihood be plotted
#' instead? Ignored if \code{do.lpdtrace == FALSE}.
#' @param x angmcmc object
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              n.chains = 1)
#' plot(fit.vmsin.20, press.enter = FALSE)
#' @export

plot.angmcmc <- function(x, par.name, comp.label, chain.no,
                         do.paramtrace = TRUE,
                         do.lpdtrace = TRUE, use.llik = FALSE,
                         press.enter = TRUE, ...)
{
  if (!is.angmcmc(x))
    stop("\'x\' must be an angmcmc object")

  if (do.paramtrace)
    paramtrace(x, par.name, comp.label, chain.no, press.enter, ...)

  if (press.enter) {
    press_enter()
    frame()
  }

  if (do.lpdtrace)
    lpdtrace(x, chain.no, use.llik)
}




