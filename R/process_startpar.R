process_startpar <- function(start_par, data,
                             ncomp, model,
                             rand_start) {
  if (is.null(start_par) & !rand_start) {
    starting <-
      do.call(paste0("start_clus_kmeans_", model),
              list(data.full=data, comp=ncomp))
  }

  else if (is.null(start_par) & rand_start) {
    starting <-
      do.call(paste0("start_clus_rand_", model),
              list(data.full=data, comp=ncomp,nstart=5))
  }

  else {

    if (model %in% c("vmsin", "vmcos", "wnorm2")) {
      allpar <- start_par
      if (any(is.null(allpar$kappa1), is.null(allpar$kappa2), is.null(allpar$kappa3),
              is.null(allpar$mu1), is.null(allpar$mu2)) ) {
        stop("too few elements in start_par, with no default")
      }

      allpar1 <- list(allpar$kappa1, allpar$kappa2, allpar$kappa3, allpar$mu1, allpar$mu2)
      allpar_len <- listLen(allpar1)
      if (min(allpar_len) != max(allpar_len)){
        stop("component size mismatch: number of components of in the starting parameter vectors differ")
      }

      starting <- list(par.mat = rbind(start_par$kappa1, start_par$kappa2,
                                       start_par$kappa3,
                                       start_par$mu1,
                                       start_par$mu2))
    }
    else if (model %in% c("vm", "wnorm")) {
      allpar <- start_par
      if (any(is.null(allpar$kappa), is.null(allpar$mu)) ) {
        stop("too few elements in start_par, with no default")
      }

      allpar1 <- list(allpar$kappa, allpar$mu)
      allpar_len <- listLen(allpar1)
      if (min(allpar_len) != max(allpar_len)){
        stop("component size mismatch: number of components of in the starting parameter vectors differ")
      }

      starting <- list(par.mat = rbind(start_par$kappa,
                                       start_par$mu))
    }

    if (ncomp == 1) {
      starting$pi.mix <- 1
    }
    else {
      starting$pi.mix <- start_par$pmix
    }
  }

  starting
}
