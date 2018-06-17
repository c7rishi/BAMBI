expand_args <- function(...){
  ell <- list(...)
  max_length <- max(vapply(ell, length, 0))
  lapply(ell, rep, length.out = max_length)
}

ldgamanum <- function(x, loc, scale) {
  (loc-1)*log(x) - x/scale
}


# mean <- function(x) {
#   sum(x)/length(x)
# }


har_mean <- function(x) {
  if(sum(x == 0) > 0) stop("zero value in harmonic mean")
  1/mean(1/x)
}


sum_sq <- function(x) sum(x^2)


prncp_reg <- function(x) x %% (2*pi)
#makes and angle in [0, 2*pi]


prncp_reg.minuspi.pi <- function(x) {
  y <- (x + pi) %% (2*pi)
  y_neg <- which(y < 0)
  y[y_neg] <- 2*pi + y[y_neg]
  y - pi
}  #makes a single angle in [-pi, pi]


atan3 <- function(x) prncp_reg(atan(x[2]/x[1]))


sph2cart <- function(x)
  c(cos(x[1])*sin(x[2]), sin(x[1])*sin(x[2]), cos(x[2]))
#calculates unit vectors from pair of angles


listLen <- function(l)
  vapply(1:length(l), function(i) length(l[[i]]), 0)


rm_NA_rad <- function(data, rad = TRUE) {
  if(length(dim(data)) == 2) {
    phi.no.na <- data[,1]
    psi.no.na <- data[,2]
    na.phi.id <- NULL
    na.psi.id <- NULL
    is.na.phi <- is.na(data[,1])
    is.na.psi <- is.na(data[,2])
    if(sum(is.na.phi) > 0)
      na.phi.id <- which(is.na.phi)
    if(sum(is.na.psi) > 0)
      na.psi.id <- which(is.na.psi)
    na.id <- union(na.phi.id, na.psi.id)
    if(length(na.id) > 0){
      phi.no.na <- data[,1][-na.id]
      psi.no.na <- data[,2][-na.id]
    }
    if(rad)    res <- prncp_reg(cbind(phi.no.na, psi.no.na))
    else res <- prncp_reg(cbind(phi.no.na, psi.no.na) * pi / 180)
    colnames(res) <- colnames(data)
  } else {
    data.no.na <- data
    na.id <- NULL
    is.na.data <- is.na(data)
    if(sum(is.na.data) > 0){
      na.id <- which(is.na.data)
      data.no.na <- data[-na.id]
    }
    if(rad) res <- prncp_reg(data.no.na)
    else res <- prncp_reg(data.no.na * pi / 180)
  }
  res
} #removes NA and converts into radians


# rdirichlet <- function (n, alpha)  # random generation from dirichlet
# {
#   len <- length(alpha)
#   x <- matrix(rgamma(len * n, alpha), ncol = len, byrow = TRUE)
#   tot <- x %*% rep(1, len)
#   x/as.vector(tot)
# }


# rnorm2 <- function (n = 1, mu, Sigma)  # random generation from biv normal
# {
#   p <- 2L
#   eS <- eigen(Sigma, symmetric = TRUE)
#   ev <- eS$values
#   X <- matrix(rnorm(p * n), n)
#   X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
#   if (n == 1) drop(X)
#   else t(X)
# }


list_by_row <- function(mat, row_index)  # create a list with elements being rows of a matrix
{
  mat.list <- lapply(1:nrow(mat), function(j) mat[j, ])
  names(mat.list) <- rownames(mat)
  mat.list
}

addtolist <- function(list_in, ...)  # add element to a list
{
  ell <- list(...)
  c(list_in, ell)
}


press_enter <- function()  # waits for the user to press [enter]
{
  cat("Press [enter] to continue")
  line <- readline()
}


kappas2sigmas_wnorm2 <- function(kappa1, kappa2, kappa3) {
  den <- kappa1*kappa2 - kappa3^2
  sigma1 <- kappa2/den
  sigma2 <- kappa1/den
  rho <- -kappa3/sqrt(kappa1*kappa2)
  c(sigma11 = sigma1, sigma22 = sigma2, rho = rho)
}

sigmas2kappas_wnorm2 <- function(sigma11, sigma22, rho) {
  den <- sigma11*sigma22*(1-rho^2)
  kappa1 <- sigma22/den
  kappa2 <- sigma11/den
  kappa3 <- -rho/((1-rho^2)*sqrt(sigma11*sigma22))
  c(kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
}

which.max_entry1 <- function(x) {
  which.max(x)[1]
}


signif_or_round <- function(x, ...) {
  for(j in length(x)) {
    if (abs(x[j]) > 1) return(round(x[j], ...))
    else return(signif(x[j], ...))
  }
}


# print est (ci_lower, ci_upper) for each element
est_ci <- function(est, lower, upper, digits = 2)
{
  out_mat <- est
  for(j in 1:length(est)) {
    out_mat[j] <- paste0(format(signif_or_round(est[j], digits), nsmall = digits),
                         " (",
                         format(signif_or_round(lower[j], digits), nsmall = digits),
                         ", ",
                         format(signif_or_round(upper[j], digits), nsmall = digits),
                         ")")
  }

  as.data.frame(out_mat)
}
