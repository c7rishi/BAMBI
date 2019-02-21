#' Contourplot for bivariate angular mixture model densities
#' @param model bivariate angular model whose mixture is of interest. Must be one of
#' "vmsin", "vmcos" and "wnorm2".
#' @param kappa1,kappa2,kappa3,mu1,mu2,pmix model parameters and mixing
#' proportions. See the respective mixture model densities (\link{dvmsinmix}, \link{dvmcosmix},
#' \link{dwnorm2mix}) for more details.
#' @param levels numeric vector of levels at which to draw contour lines;
#' passed to the \link{contour} function in graphics.
#' @param nlevels	number of contour levels desired \strong{if} levels is not supplied;
#' passed to the \link{contour} function in graphics.
#' @param xpoints Points on the first (x-) coordinate where the density is to be evaluated.
#' Default to seq(0, 2*pi, length.out=100).
#' @param ypoints Points on the first (x-) coordinate where the density is to be evaluated.
#' Default to seq(0, 2*pi, length.out=100).
#' @param ... additional model specific argment
#' @param xlab,ylab,col,lty,main graphical parameters passed to \link{contour}.
#' @examples
#' contour_model('vmsin', 1, 1, 1.5, pi, pi)
#' contour_model('vmcos', 1, 1, 1.5, pi, pi)
#'
#'
#' @export

contour_model <- function(model = "vmsin", kappa1, kappa2, kappa3, mu1, mu2,
                          pmix = rep(1/length(kappa1), length(kappa1)),
                          xpoints = seq(0, 2*pi, length.out = 100),
                          ypoints = seq(0, 2*pi, length.out = 100),
                          levels, nlevels = 20,
                          xlab="x", ylab="y", col="black",
                          lty=1, main, ...)

{
  if (!model %in% c("vmsin", "vmcos", "wnorm2"))
    stop("model must be one of \"vmsin\", \"vmcos\" or \"wnorm2\"")

  if (missing(levels)) {
    levels <- exp(seq(-20,2, length.out = nlevels))
  }

  if (missing(main))
    main <- paste("Contour plot for", length(kappa1), "component",
                  model, "mixture density")

  coords <- as.matrix(expand.grid(xpoints, ypoints))
  inargs <- list(x = coords, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3,
                 mu1 = mu1, mu2 = mu2, pmix = pmix, ...)
  dens <- do.call(paste0("d", model, "mix"), inargs)

  contour(xpoints, ypoints, matrix(dens, nrow=length(xpoints)), levels=levels,
          xlab=xlab, ylab=ylab, col=col, lty=lty, main=main)

}



#' Surface for bivariate angular mixture model densities
#' @inheritParams contour_model
#' @param kappa1,kappa2,kappa3,mu1,mu2,pmix model parameters and mixing
#' proportions. See the respective mixture model densities (\link{dvmsinmix}, \link{dvmcosmix},
#' \link{dwnorm2mix}) for more details.
#' @param log.density logical. Should log density be used for the plot?
#' @param xlab,ylab,zlab,main graphical parameters passed to \code{lattice::wireframe}
#' @param ... additional arguments passed to \code{lattice::wireframe}
#' @examples
#' surface_model('vmsin', 1, 1, 1.5, pi, pi)
#' surface_model('vmcos', 1, 1, 1.5, pi, pi)
#'
#' @export

surface_model <- function(model = "vmsin", kappa1, kappa2, kappa3, mu1, mu2,
                          pmix = rep(1/length(kappa1), length(kappa1)),
                          xpoints = seq(0, 2*pi, length.out = 30),
                          ypoints = seq(0, 2*pi, length.out = 30),
                          log.density = FALSE, xlab="x", ylab="y",
                          zlab = ifelse(log.density, "Log Density", "Density"),
                          main, ...)
{
  if (!model %in% c("vmsin", "vmcos", "wnorm2"))
    stop("model must be one of \"vmsin\", \"vmcos\" or \"wnorm2\"")

  if (missing(main))
    main <- paste("Density surface for", length(kappa1), "component",
                  model, "mixture density")


  coords <- as.matrix(expand.grid(xpoints, ypoints))
  inargs <- list(x = coords, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3,
                 mu1 = mu1, mu2 = mu2, pmix = pmix, ...)
  dens <- do.call(paste0("d", model, "mix"), inargs)

  denmat <- matrix(dens, nrow = length(xpoints))

  if(log.density) {
    denmat <- log(denmat)
  }

  print(basic_surfaceplot(xpoints = xpoints, ypoints = ypoints,
                          denmat = denmat, xlab = xlab, ylab = ylab,
                          main = main, zlab = zlab, ...))

  # nrden <- nrow(denmat)
  # ncden <- ncol(denmat)
  #
  # denfacet <- denmat[-1, -1] + denmat[-1, -ncden] +
  #   denmat[-nrden, -1] + denmat[-nrden, -ncden]
  #
  # # Create a function interpolating colors in the range of specified colors
  # jet.colors <- grDevices::colorRampPalette( c("blue", "green",
  #                                              "yellow", "orange", "red") )
  # # Generate the desired number of colors from this palette
  # nbcol <- 500
  # color <- jet.colors(nbcol)
  #
  # denfacet <- denmat[-1, -1] + denmat[-1, -ncden] +
  #   denmat[-nrden, -1] + denmat[-nrden, -ncden]
  # # Recode facet z-values into color indices
  # facetcol <- cut(denfacet, nbcol)
  #
  # persp(x=xpoints, y=ypoints, z=denmat, theta = theta, phi = phi,
  #       expand = expand, col = color[facetcol],
  #       ltheta = 120, shade = shade, ticktype = "detailed",
  #       xlab = xlab, ylab = ylab, zlab = zlab,
  #       main = main, ...) -> res

  # wireframe(x=denmat~x*y,
  #           data=data.frame(x=xpoints,
  #                           y=rep(ypoints, each=length(xpoints)),
  #                           denmat=denmat),
  #           xlab=xlab, ylab=ylab, zlab=zlab,
  #           main=main, colorkey=colorkey, outerbox=FALSE,
  #           par.settings=list(axis.line=list(col = 'transparent')),
  #           neval=length(xpoints), aspect = c(61/87, 0.4),
  #           drape=TRUE, light.source=c(10,0,10),
  #           col.regions=colorRampPalette(c("blue", "green",
  #                                          "yellow", "orange",
  #                                          "red"))(100),
  #           scales=list(arrows=FALSE, col=1)
  # )
  # do.call(wireframe, inargs)
  # inargs$outerbox <- FALSE
  # inargs$par.settings <- list(axis.line = list(col = 'transparent'))
  # if (is.null(inargs$xlab)) inargs$xlab <- expression(Theta)
  # if (is.null(inargs$ylab)) inargs$ylab <- expression(Phi)
  # if (is.null(inargs$colorkey)) inargs$colorkey <- FALSE
  # if (is.null(inargs$main)) inargs$main <- list(maintext, font=2, cex=1.5)
  # if (is.null(inargs$neval)) inargs$neval <- 100
  # if (is.null(inargs$aspect)) inargs$aspect <- c(61/87, 0.4)
  # if (is.null(inargs$zlab)) inargs$zlab <- list("Density", rot=90)
  # if (is.null(inargs$screen)) inargs$screen <- list(z=45,x=-45)
  # if (is.null(inargs$colorkey)) inargs$colorkey <- TRUE
  # if (is.null(inargs$scales))
  #   inargs$scales <- list(
  #     arrows = FALSE,
  #     x=list(at = c(-pi, -pi/2, 0, pi, pi/2),
  #            labels = c(expression(-pi), expression(-pi/2), 0,
  #                       expression(pi), expression(pi/2))),
  #     y=list(at = c(-pi, -pi/2, 0, pi, pi/2),
  #            labels = c(expression(-pi), expression(-pi/2), 0,
  #                       expression(pi), expression(pi/2))),
  #     col = "black", font = 1, tck=1)
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
  #                               right.padding=0)
  #
  # inargs$par <- c(k1,k2,k3,0,0)
  #
  #
}
