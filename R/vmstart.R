#' @keywords internal
Ainv <- function(x) {
  if(x < 0.53) (2*x + x^3 + 5/6*x^5)
  else if(x >= 0.86 && x < 0.95) ((9 - 8*x + 3*x^2) / (8 * (1 - x)))
  else ((1.28 - 0.53*x^2) * tan(x*pi/2))
}


#' @keywords internal
A_bessel <- function(x){
  besselI(x, 1, TRUE)/besselI(x, 0, TRUE)
}


#' @keywords internal
start_par_vm <- function(data.sub) {
  x1 <- data.sub
  Sbar <- mean(sin(x1))
  Cbar <- mean(cos(x1))
  muhat <- atan(Sbar/Cbar) + pi * (Cbar < 0)
  Rbar <- sqrt(Sbar^2 + Cbar^2)
  k <- Ainv(Rbar)
  c(k, prncp_reg(muhat))
}  #starting parameters from  a dataset

#' @keywords internal
start_clus_kmeans_vm <- function(data.full, comp = 2, nstart = 10){
  data.full.cart <- t(sapply(data.full, function(x) c(cos(x), sin(x))))
  data.kmean <- kmeans(data.full.cart, centers = comp, nstart = nstart)
  ids <- data.kmean$cluster
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_vm(data.full[clust.ind[[m]]]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(order(apply(par, 2, function(x) sum_sq(x[1:3]))))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #kmeans, then start_par for each cluster

#' @keywords internal
start_clus_rand_vm <- function(data.full, comp = 2, nstart = 10) {
  rand.pi  <- runif(comp, 1/(2*comp), 2/comp)
  rand.pi <- rand.pi/sum(rand.pi)

  rand.multinom <- t(rmultinom(length(data.full), 1, rand.pi))
  ids <- apply(rand.multinom, 1, which.max)
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_vm(data.full[clust.ind[[m]]]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(order(apply(par, 2, function(x) sum_sq(x[1:3]))))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #random groups, then start_par for each cluster




#' @keywords internal
start_par_vmsin <- function(data.sub) {
  x1 <- data.sub[,1]; y1 <- data.sub[,2]
  Sbarphi <- mean(sin(x1))
  Cbarphi <- mean(cos(x1))
  phibar <- atan(Sbarphi/Cbarphi) + pi * (Cbarphi < 0)
  Rbarphi <- sqrt(Sbarphi^2 + Cbarphi^2)
  k1 <- Ainv(Rbarphi)

  Sbarpsi <- mean(sin(y1))
  Cbarpsi <- mean(cos(y1))
  psibar <- atan(Sbarpsi/Cbarpsi) + pi * (Cbarpsi < 0)
  Rbarpsi <- sqrt(Sbarpsi^2 + Cbarpsi^2)
  k2 <- Ainv(Rbarpsi)


  sindiffphi <- sin(outer(x1, x1, "-"))
  sindiffpsi <- sin(outer(y1, y1, "-"))
  rho <- sum(sindiffphi*sindiffpsi)/sum(sindiffphi^2)/sum(sindiffpsi^2)

  c(k1, k2, rho*sqrt(k1*k2), prncp_reg(phibar), prncp_reg(psibar))
}  #starting parameters from  a dataset

#' @keywords internal
start_clus_kmeans_vmsin <- function(data.full, comp = 2, nstart = 10){
  data.full.cart <- t(apply(data.full, 1, sph2cart))
  data.kmean <- kmeans(data.full.cart, centers = comp, nstart = nstart)
  ids <- data.kmean$cluster
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_vmsin(data.full[clust.ind[[m]],]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(order(apply(par, 2, function(x) sum_sq(x[1:3]))))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #kmeans, then start_par for each cluster

#' @keywords internal
start_clus_rand_vmsin <- function(data.full, comp = 2, nstart = 10) {
  rand.pi  <- runif(comp, 1/(2*comp), 2/comp)
  rand.pi <- rand.pi/sum(rand.pi)

  rand.multinom <- t(rmultinom(nrow(data.full), 1, rand.pi))
  ids <- apply(rand.multinom, 1, which.max)
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_vmsin(data.full[clust.ind[[m]],]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(order(apply(par, 2, function(x) sum_sq(x[1:3]))))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #random groups, then start_par for each cluster






#' @keywords internal
start_par_vmcos <- function(data.sub) {
  x1 <- data.sub[,1];
  y1 <- data.sub[,2]
  Sbarphi <- mean(sin(x1))
  Cbarphi <- mean(cos(x1))
  phibar <- atan(Sbarphi/Cbarphi) + pi * (Cbarphi < 0)
  Rbarphi <- sqrt(Sbarphi^2 + Cbarphi^2)
  k1 <- Ainv(Rbarphi)
  #k1 <- min(15, Ainv(Rbarphi))

  Sbarpsi <- mean(sin(y1))
  Cbarpsi <- mean(cos(y1))
  psibar <- atan(Sbarpsi/Cbarpsi) + pi * (Cbarpsi < 0)
  Rbarpsi <- sqrt(Sbarpsi^2 + Cbarpsi^2)
  k2 <- Ainv(Rbarpsi)
  #k2 <- min(15, Ainv(Rbarpsi))

  mu <- prncp_reg(phibar)
  nu <- prncp_reg(psibar)

  Sbarphi_psi <- mean(sin(x1-y1))
  Cbarphi_psi <- mean(cos(x1-y1))
  phi_psibar <- atan(Sbarphi_psi/Cbarphi_psi) + pi * (Cbarphi_psi < 0)
  Rbarphi_psi <- sqrt(Sbarphi_psi^2 + Cbarphi_psi^2)
  k3.unsgn <- Ainv(Rbarphi_psi)


  sindiffphi <- sin(outer(x1, x1, "-"))
  sindiffpsi <- sin(outer(y1, y1, "-"))
  rho <- sum(sindiffphi*sindiffpsi)/sum(sindiffphi^2)/sum(sindiffpsi^2)

  # Sbarphi_psi <- mean(sin(x1-y1+mu-nu))
  # Cbarphi_psi <- mean(cos(x1-y1+mu-nu))
  # Rbarphi_psi <- sqrt(Sbarphi_psi^2 + Cbarphi_psi^2)
  k3 <- sign(rho)*k3.unsgn



  c(k1, k2, k3, mu, nu)
}  #starting parameters from  a dataset

#' @keywords internal
start_clus_kmeans_vmcos <- function(data.full, comp = 2, nstart = 10){
  data.full.cart <- t(apply(data.full, 1, sph2cart))
  data.kmean <- kmeans(data.full.cart, centers = comp, nstart = nstart)
  ids <- data.kmean$cluster
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_vmcos(data.full[clust.ind[[m]],]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(order(apply(par, 2, function(x) sum_sq(x[1:3]))))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #kmeans, then start_par for each cluster

#' @keywords internal
start_clus_rand_vmcos <- function(data.full, comp = 2, nstart = 10){
  rand.pi  <- runif(comp, 1/(2*comp), 2/comp)
  rand.pi <- rand.pi/sum(rand.pi)

  rand.multinom <- t(rmultinom(nrow(data.full), 1, rand.pi))
  ids <- apply(rand.multinom, 1, which.max)
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_vmcos(data.full[clust.ind[[m]],]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(order(apply(par, 2,  sum_sq)))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #random groups, then start_par for each cluster

