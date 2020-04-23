# based on uniroot.all function from package rootSolve
# the package is now archived
uniroot.all <- function(
  f,
  interval,
  lower = min(interval),
  upper = max(interval),
  tol = .Machine$double.eps^0.2,
  maxiter = 1000,
  trace = 0,
  n = 100,
  ...
) {
  xseq <- seq(lower, upper, len = n + 1)
  f_xseq <- f(xseq, ...)
  out <- xseq[which(f_xseq == 0)]
  f_sign <- f_xseq[1:n] * f_xseq[2:(n + 1)]
  i_range <- which(f_sign < 0)
  for (i in i_range)  {
    out <- c(
      out,
      uniroot(f, lower = xseq[i],
              upper = xseq[i + 1],
              maxiter = maxiter,
              tol = tol,
              trace = trace,
              ...)$root
    )
  }
  out
}
