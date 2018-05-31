start_par_wnorm <- function(data.sub) {
  x1 <- prncp_reg(data.sub)

  Sbar <- mean(sin(x1))
  Cbar <- mean(cos(x1))
  muhat <- prncp_reg(atan(Sbar/Cbar))
  Rbar <- sqrt(Sbar^2 + Cbar^2)

  c(1/(1-Rbar), muhat)
}  #starting parameters from  a dataset

start_clus_kmeans_wnorm <- function(data.full, comp = 2, nstart = 5){
  data.full.cart <- t(sapply(data.full, function(x) c(cos(x), sin(x))))
  data.kmean <- kmeans(data.full.cart, centers = comp, nstart = nstart)
  ids <- data.kmean$cluster
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_wnorm(data.full[clust.ind[[m]]]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(colSums(par^2))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #kmeans, then start_par for each cluster



start_clus_rand_wnorm <- function(data.full, comp = 2){
  rand.pi  <- runif(comp, 1/(2*comp), 2/comp)
  rand.pi <- rand.pi/sum(rand.pi)
  rand.multinom <- t(rmultinom(length(data.full), 1, rand.pi))
  ids <- apply(rand.multinom, 1, which.max)
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_wnorm(data.full[clust.ind[[m]]]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(colSums(par^2))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #random groups, then start_par for each cluster


start_par_wnorm2 <- function(data.sub) {
  data.sub.0.2pi <- prncp_reg(data.sub)
  mu <- prncp_reg(colMeans(data.sub.0.2pi))
  sigma.inv.vec <- solve(cov(data.sub.0.2pi))[c(1,4,2)]
  c(sigma.inv.vec, mu)
}  #starting parameters from  a dataset


start_clus_kmeans_wnorm2 <- function(data.full, comp = 2, nstart = 5){
  data.full.cart <- t(apply(data.full, 1, sph2cart))
  data.kmean <- kmeans(data.full.cart, centers = comp, nstart = nstart)
  ids <- data.kmean$cluster
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_wnorm2(data.full[clust.ind[[m]],]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(colSums(par^2))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #kmeans, then start_par for each cluster


start_clus_rand_wnorm2 <- function(data.full, comp = 2, nstart = 5){
  rand.pi  <- runif(comp, 1/(2*comp), 2/comp)
  rand.pi <- rand.pi/sum(rand.pi)
  rand.multinom <- t(rmultinom(nrow(data.full), 1, rand.pi))
  ids <- apply(rand.multinom, 1, which.max)
  clust.ind <- lapply(1:comp, function(i) which(ids == i))
  par <- sapply(1:length(clust.ind), function(m) start_par_wnorm2(data.full[clust.ind[[m]],]))
  pi.mix <- listLen(clust.ind)/length(unlist(clust.ind))
  order.conc <- order(colSums(par^2))
  list("par.mat" = par[,order.conc], "pi.mix" = pi.mix[order.conc], "clust.ind" = clust.ind[order.conc], "id" = ids)
}  #random groups, then start_par for each cluster
