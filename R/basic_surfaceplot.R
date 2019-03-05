basic_surfaceplot <- function(xpoints, ypoints, denmat,
                              xlab, ylab, zlab,
                              main = main, ...) {
  inargs <- list(...)
  inargs$x <- denmat~x*y
  inargs$data <- data.frame(x = xpoints,
                            y=rep(ypoints, each=length(xpoints)),
                            denmat=denmat)
  if (is.null(inargs$outerbox)) inargs$outerbox <- FALSE
  if (is.null(inargs$par.settings)) inargs$par.settings <- list(axis.line = list(col = 'transparent'))
  if (is.null(inargs$xlab)) inargs$xlab <- xlab
  if (is.null(inargs$ylab)) inargs$ylab <- ylab
  if (is.null(inargs$colorkey)) inargs$colorkey <- FALSE
  if (is.null(inargs$main)) inargs$main <- main
  if (is.null(inargs$neval)) inargs$neval <- 100
  if (is.null(inargs$aspect)) inargs$aspect <- c(61/87, 0.4)
  if (is.null(inargs$zlab)) inargs$zlab <- list("Density", rot=90)
  if (is.null(inargs$screen)) inargs$screen <- list(z=-30, x=-60)
  if (is.null(inargs$colorkey)) inargs$colorkey <- FALSE
  if (is.null(inargs$scales))
    inargs$scales <- list(arrows=FALSE, col=1)
  if (is.null(inargs$drape)) inargs$drape <- TRUE
  if (is.null(inargs$light.source))
    inargs$light.source <- c(10,0,10)
  if (is.null(inargs$col.regions))
    inargs$col.regions <- grDevices::colorRampPalette(c("steelblue", "green",
                                                        "yellow", "orange", "red"))(60)
  if (is.null(inargs$par.settings))
    inargs$par.settings <-
    list(axis.line = list(col = 'transparent'
    ),
    layout.heights = list(
      top.padding = 0,
      main.key.padding = 0,
      key.axis.padding = 0,
      axis.xlab.padding = 0,
      xlab.key.padding = 0,
      key.sub.padding = 0,
      bottom.padding = 0
    ),
    layout.widths = list(
      left.padding = 0,
      key.ylab.padding = 0,
      ylab.axis.padding = 0,
      axis.key.padding = 0,
      right.padding = 0
    ))

  if (is.null(inargs$zoom))
    inargs$zoom <- 0.85


  do.call(lattice::wireframe, inargs)
}
