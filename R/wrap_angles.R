#' Wrap angles into \code{[-pi, pi]} or \code{[0, 2*pi]}
#' @param x numeric vector or matrix or data.frame.
#' @details
#' \code{minuspi_to_pi} wraps \code{x} into \code{[-pi, pi]},
#' while \code{zero_to_pi} wraps \code{x} into \code{[0, 2*pi]}.
#'
#' @examples
#' dat <- matrix(runif(100, -pi, pi), ncol=2)
#' dat1 <- zero_to_2pi(dat)
#' dat2 <- minuspi_to_pi(dat1)
#' all.equal(dat, dat2)
#' @export

zero_to_2pi <- function(x)
{
  x %% (2*pi)
}


#' @rdname zero_to_2pi
#' @export

minuspi_to_pi <- function(x)
{
  prncp_reg.minuspi.pi(x)
}
