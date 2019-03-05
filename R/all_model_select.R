#' Stepwise fitting of angular mixture models with incremental component sizes and optimum model selection
#' @inheritParams pointest
#' @inheritParams fit_angmix
#' @param start_ncomp starting component size. A single component model is fitted if \code{start_ncomp} is equal to one.
#' @param max_ncomp maximum number of components allowed in the mixture model.
#' @param crit model selection criteria, one of \code{"LOOIC", "WAIC", "AIC", "BIC", "DIC"} or \code{"LOGML"}. Default is
#' \code{"LOOIC"}.
#' @param L HMC tuning parameter (trajectory length) passed to \link{fit_angmix}. Can be a numeric vetor (or scalar), in which case
#' the same \code{L} is passed to all \link{fit_angmix} calls, or can be a list of length \code{max_ncomp-start_ncomp+1},
#' so that \code{L_list[[i]]} is passed as the argument \code{L} to \link{fit_angmix} call with \code{ncomp = max_ncomp+i-1}. See
#' \link{fit_angmix} for more details on \code{L} including its default values. Ignored if \code{method = "rwmh"}.
#' @param fn function to evaluate on MCMC samples to estimate parameters.
#' Defaults to \code{mean}, which computes the estimated posterior means. If \code{fn = max},
#' then MAP estimate is calculated from the MCMC run. Used only if \code{crit = "DIC"}, and ignored otherwise.
#' @param prev_par logical. Should the MAP estimated parameters from the model with \code{ncomp = K} be used in the model
#' with \code{ncomp = K+1} as the starting parameters, with the component with largest mixing proportion appearing twice in the
#' bigger model?
#' @param form form of crit to be used. Available choices are 1 and 2. Used only if \code{crit} is \code{"DIC"} and ignored otherwise.
#' @param ... additional arguments passed to \link{fit_angmix}.
#' @param logml_maxiter maximum number of iterations (\code{maxiter}) passed to \link{bridge_sampler} for calculating
#' \code{LOGML}. Ignored if \code{crit} is not \code{LOGML}.
#' @param fix_label logical. Should the label switchings on the current fit (only the corresponding "best chain" if \code{use_best_chain = TRUE})
#' be fixed before computing parameter estimates and model selection criterion? Defaults to \code{TRUE} if \code{perm_sampling} is true in
#' the \link{fit_angmix} call, or if \code{crit = "DIC"} and \code{form = 1}.
#' @param silent logical. Should the current status (such as what is the current component labels, which job is being done etc.)
#' be printed? Defaults to \code{TRUE}.
#' @param return_all logical. Should all angmcmc objects obtained during step-wise run be returned? *Warning*: depending on the
#' sizes of \code{n.iter}, \code{start_ncomp}, \code{max_ncomp} and \code{n.chains}, this can be very memory intesive. In such
#' cases, it is recommended that \code{return_all} be set to \code{FALSE}, and, if required, the intermediate fitted objects be
#' saved to file by setting \code{save_fits = TRUE}.
#' @param save_fits logical. Should the intermediate angmcmc objects obtained during step-wise run be saved
#'  to file using \link{save}? Defaults to TRUE. See \code{save_file} and \code{save_dir}.
#' @param save_file,save_dir \code{save_file} is a list of size \code{max_ncomp-start_ncomp+1},
#' with k-th entry providing the \code{file}
#' argument used to \link{save} the intermediate angmcmc object with \code{ncomp = k} (titled \code{"fit_angmcmc"}).
#' If not provided, then k-th element
#' of \code{save_file[[k]]} is taken to be \code{\link{paste}(save_dir, "comp_k", sep="/")}. Both are ignored if
#' \code{save_fits = FALSE}.
#' @param use_best_chain logical. Should only the "best" chain obtained duing each intermediate fit be used during
#' computation of model selection criterion? Here "best" means the chain
#' with largest (mean over iterations) log-posterior density. This can be helpful if one of the chains gets stuck at local optima. Defaults to TRUE.
#' @param return_llik_contri passed to \link{fit_angmix}. By default, set to \code{TRUE} if \code{crit} is either \code{"LOOIC"}
#' or \code{"WAIC"}, and to \code{FALSE} otherwise.
#'
#' @details
#' The goal is to fit an angular mixture model with an optimally chosen component size K.
#' To obtain an optimum K, mixture models with incremental component sizes
#' between \code{start_ncomp} and \code{max_ncomp} are fitted stepwise using \link{fit_angmix}.
#' The mixture model with the first minimum value of the model selection criterion \code{crit}
#' is taken as the best model.
#'
#' Note that in each intermediate fitted model, the total number of components (instead of the number of
#' "non-empty components") in the model is used to estimate of the true component
#' size, and then the fitted model is penalized for model complexity (via the model selection criterion).
#' This approach of selecting an optimal K follows the perspective "let two component specific parameters
#' be identical" for overfitting mixtures, and as such the  Dirichlet prior hyper-parameters \code{pmix.alpha}
#' (passed to \link{fit_angmix}) should be large. See  Fruhwirth-Schnatter (2011) and Mengerson for more deltails.
#'
#' Note that the stability of \link{bridge_sampler} used in marginal likelihood estimation heavily depends on stationarity of the
#' chains. As such, while using this criterion, we recommending running the chain long engouh, and setting \code{fix_label = TRUE}
#' for optimal performance.
#'
#' @references
#' Fruhwirth-Schnatter, S.: Label switching under model uncertainty. In: Mengerson, K., Robert, C., Titterington, D. (eds.) Mixtures:
#' Estimation and Application, pp. 213-239. Wiley, New York (2011).
#'
#'
#'
#'
#' @return Returns a named list (with class = \code{stepfit}) with the following seven elements:
#'
#' \code{fit.all} (if \code{return_all = TRUE}) - a list all angmcmc objects created at each component size;
#'
#' \code{fit.best} - angmcmc object corresponding to the optimum component size;
#'
#' \code{ncomp.best} - optimum component size (integer);
#'
#' \code{crit} - which model comparison criterion used (one of \code{"LOOIC", "WAIC", "AIC", "BIC", "DIC"} or \code{"LOGML"});
#'
#' \code{crit.all} - all \code{crit} values calculated (for all component sizes);
#'
#' \code{crit.best} - \code{crit} value for the optimum component size; and
#'
#' \code{maxllik.all} - maximum (obtained from MCMC iterations) log likelihood for all fitted models
#'
#' \code{maxllik.best} - maximum log likelihodd for the optimal model; and
#'
#' \code{check_min} - logical; is the optimum component size less than \code{max_ncomp}?
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' set.seed(1)
#' fit.vmsin.step.15 <- fit_incremental_angmix("vmsin", tim8, "BIC", start_ncomp = 1,
#'                                           max_ncomp = 3, n.iter = 15,
#'                                           n.chains = 1, save_fits=FALSE)
#' (fit.vmsin.best.15 <- bestmodel(fit.vmsin.step.15))
#' lattice::densityplot(fit.vmsin.best.15)
#'
#' @export


fit_incremental_angmix <- function(model, data,
                                   crit = "LOOIC",
                                   start_ncomp=1, max_ncomp=10,
                                   L = NULL,
                                   fn = mean,
                                   fix_label = NULL,
                                   form = 2,
                                   start_par = NULL,
                                   prev_par = TRUE,
                                   logml_maxiter = 1e4,
                                   return_all = FALSE,
                                   save_fits = FALSE,
                                   save_file = NULL,
                                   save_dir = "",
                                   silent = FALSE,
                                   return_llik_contri = (crit %in% c("LOOIC", "WAIC")),
                                   use_best_chain = TRUE,
                                   ...)
{
  if (length(model) > 1)
    stop("\'model\' must be a scalar")
  if(missing(model))
    stop("argument \"model\" is missing, with no default")
  if(!crit %in% c("AIC", "BIC", "DIC", "WAIC", "LOOIC", "LOGML"))
    stop("non-compatible criterion")
  if (model %in% c("vmsin", "vmcos", "wnorm2")) {
    type <- "bi"
  } else if (model %in% c("vm", "wnorm")) {
    type <- "uni"
  } else  {
    stop("non-compatible model")
  }

  if (!missing(start_par)) {
    if (min(listLen(start_par)) != max(listLen(start_par)))
      stop("Lengths of elements in start_par differ")
    else if((listLen(start_par))[1] != start_ncomp)
      stop("Number of components in start_par must be equal to \'start_ncomp\'")
  }


  if(start_ncomp > max_ncomp)
    stop("\'start_ncomp\' cannot be smaller than \'max_ncomp\'")

  all_ncomp <- start_ncomp:max_ncomp
  n_ncomp <- length(all_ncomp)


  if (is.null(save_file)) {
    save_file <- lapply(all_ncomp, function(j) paste0(save_dir, "/comp_", j, ".Rdata"))
  }
  else if (!is.list(save_file) | length(save_file) != n_ncomp)
    stop("\'save_file\' must be a list of length max_ncomp-start_ncomp+1")


  crit_print <- crit

  if (crit == "LOGML")
    crit_print <- "(negative) LOGML"


  if (type == "bi") {
    if (!(is.matrix(data) | is.data.frame(data)))
      stop("\'data\' must be a two column matrix for model = \'vmsin\', \'vmcos\' and \'wnorm2\'")

    if (ncol(data) != 2)
      stop("\'data\' must be a two column matrix for model = \'vmsin\', \'vmcos\' and \'wnorm2\'")

    data.rad <- rm_NA_rad(data)
    n.data <- nrow(data.rad)

  }
  else  {
    if (!is.numeric(data))
      stop("\'data\' must be a vector for \'model\' = \'vm\' and \'wnorm\'")
  }

  fit_all <- NULL

  if (return_all)
    fit_all <- vector("list", n_ncomp)


  ell <- list(...)

  if (is.null(ell$perm_sampling))
    ell$perm_sampling <- formals(fit_angmix)$perm_sampling

  if (is.null(fix_label)) {
    if(any(form == 1 & crit == "DIC", ell$perm_sampling & prev_par,
           ell$perm_sampling & crit == "LOGML")) {
      fix_label <- TRUE
    } else {
      fix_label <- FALSE
    }

  }


  # all_fit <- vector("list", n_ncomp)
  all_input <- list("data" = data, "model" = model,
                    return_llik_contri = return_llik_contri,
                    ...)

  formal_fit_angmix <- formals(fit_angmix)

  if (is.null(all_input$n.chains))
    all_input$n.chains <- formal_fit_angmix$n.chains

  n.chains <- all_input$n.chains


  # when missing L, get default value of L and make L_list
  if (is.null(L)) {
    L_list <- lapply(1:(max_ncomp-start_ncomp+1),
                     function(j) {
                       ncomp <- start_ncomp-j+1
                       eval(formal_fit_angmix$L)
                     })
  }
  # when L is not null, check if it's in correct format first
  else if (all(!is.list(L) & !is.numeric(L),
               is.list(L) & length(L) != (max_ncomp - start_ncomp + 1),
               !is.numeric(L)))
    stop("L must either be a list of length max_ncomp-start_ncomp+1 or a vector")

  else if (!is.list(L) & is.numeric(L)) {
    L_list <- lapply(1:(max_ncomp - start_ncomp + 1),
                     function(j) L)
    # all elements of L_list are the same
  }

  else if (is.list(L)) {
    L_list <- L
    # L_list is just L, when given properly
  }


  all_par_est <- vector("list", length = max_ncomp-start_ncomp+1)


  # all_crit <- rep(0, n_ncomp)
  all_crit <- vector("list", length = n_ncomp)
  all_maxllik <- rep(0, n_ncomp)

  if(!form %in% 1:2) form <- 1

  check_min <- FALSE


  for(j in seq_len(length(all_ncomp))) {
    all_input$ncomp <- all_ncomp[j]
    all_input$L <- L_list[[j]]

    # copy the previous fit as fit_prev if j > 1
    if (j > 1) {
      fit_prev <- fit_angmcmc
      rm(fit_angmcmc)
      gc()
    }

    # starting parameters for j > 1, ncomp >= 3
    if (j > 1 & prev_par & all_ncomp[j] > 2) {
      all_par <- all_par_est[[j-1]]
      # find the component with largest mix_prop
      copy_comp_id <- which.max(all_par[1, ])[1]
      new_comp <- all_par[, copy_comp_id]
      new_comp_id <- all_input$ncomp
      all_par <- cbind(all_par, new_comp)
      # distribute the weights between the "new" and "old" components
      all_par[1, c(copy_comp_id, new_comp_id)] <-
        all_par[1, copy_comp_id]/2

      colnames(all_par) <- 1:new_comp_id
      start_par <- list_by_row(all_par)
      all_input$start_par <- start_par
    }

    else if (j == 1 & !missing(start_par)) {
      all_input$start_par <- start_par
    }
    else {
      all_input$start_par <- NULL
    }







    if (!silent) {
      cat("**************\n")
      cat(paste("Fitting", model,
                "mixture model with ncomp = ",
                all_ncomp[j], "...\n") )

    }



    fit_angmcmc <- do.call(fit_angmix, all_input)


    all_maxllik[j] <- maxllik_curr <- max(fit_angmcmc$llik[fit_angmcmc$final_iter, ])

    if (!silent)
      cat(paste("\nMaximum log likelihood (from MCMC iterations) =",
                round(maxllik_curr, 3), "\n"))


    if (use_best_chain) {
      best.chain.id <- which.max(
        vapply(1:fit_angmcmc$n.chains,
               function (j) mean(fit_angmcmc$lpd[fit_angmcmc$final_iter, j]),
               0))
      fit_angmcmc_adj <- select_chains(fit_angmcmc, best.chain.id)
    } else {
      fit_angmcmc_adj <- fit_angmcmc
    }


    if (fix_label & all_ncomp[j] > 1) {
      if(!silent)
        cat("\nCalling fix_label with default settings ...\n")

      fit_angmcmc_adj <- fix_label(fit_angmcmc_adj)

      # replace fit_angmcmc by fit_angmcmc_adj if fix_label and !use_best_chain
      if (!use_best_chain) {
        fit_angmcmc <- fit_angmcmc_adj
      }

    }
    else if (!silent) {
      cat("\nSkipping fix_label call ...\n")
    }



    all_par_est[[j]] <- pointest(fit_angmcmc_adj, fn = "MODE")



    if(save_fits) {
      if (!silent)
        cat(paste0("\nSaving the output (titled \'fit_angmcmc\') with filename = \'",
                   save_file[[j]],  "\' ...\n"))
      save(fit_angmcmc, file=save_file[[j]])
    }


    if (return_all)
      fit_all[[j]] <- fit_angmcmc



    if(!silent)
      cat("\nComputing model selection criterion ...\n")

    if (crit == "WAIC") {
      curr_crit <- suppressWarnings(loo::waic(fit_angmcmc_adj))
      all_crit[[j]] <- curr_crit
    }

    else if (crit == "LOOIC") {
      curr_crit <- suppressWarnings(loo::loo(fit_angmcmc_adj))
      all_crit[[j]] <- curr_crit
    }

    else if (crit == "LOGML") {
      curr_crit <- tryCatch(bridgesampling::bridge_sampler(fit_angmcmc_adj, silent = TRUE, maxiter = logml_maxiter),
                            error = function(e) "error")

      if (unlist(curr_crit)[1] == "error")
        stop(paste0("log posterior too unstable with ncomp = ",
                    all_ncomp[j], " to calculate log ML. Try a different criterion."))
      all_crit[[j]] <- -curr_crit$logml
    }

    else if (crit == "DIC") {
      curr_crit <- DIC(fit_angmcmc_adj, form=form)
      all_crit[[j]] <- curr_crit["DIC"]
    }

    else if (crit == "AIC") {
      all_crit[[j]] <- AIC(fit_angmcmc_adj)
    }

    else {
      all_crit[[j]] <- BIC(fit_angmcmc_adj)
    }

    # if(j > start_ncomp) cat("\n")

    if(!silent) {
      crit_val_print <- ""
      if (!crit %in% c("LOOIC", "WAIC"))
        crit_val_print <- round(all_crit[[j]], 3)

      cat(paste("\t", "ncomp = ", all_ncomp[j], ",\t",
                crit_print, ":", crit_val_print))

      if (crit %in% c("LOOIC", "WAIC")) {
        cat("\n")
        cat(suppressWarnings(capture.output(all_crit[[j]])), sep = "\n")
      }
      cat("\n")
      cat("**************\n\n")
    }

    if (j > 1 ) {
      if (crit %in% c("LOOIC", "WAIC")) {
        compare_crit <- suppressWarnings(loo::compare(all_crit[[j-1]],
                                                      all_crit[[j]]))
        # test for signif improvement in elpd
        # H0: curr elpd - prev elpd <= 0 vs Ha: >
        zscore <- compare_crit["elpd_diff"]/compare_crit["se"]
        if (zscore <= 1.645) {
          # fail to reject null at alpha = 0.05 --
          # so no signific improvement in curr elpd compared to prev
          check_min <- TRUE
          j.best <- j-1
          cat("\nImprovement in predicitive accuracy not significant. Stopping...\n")
          fit_best <- fit_prev #previous fit is best
          break
        }
      } else if (all_crit[[j]] > all_crit[[j-1]]) {
        check_min <- TRUE
        j.best <- j-1
        cat("\nFirst minimum attained. Stopping...\n")
        fit_best <- fit_prev #previous fit is best
        break
      }
    }


    if (all_ncomp[j] == max_ncomp) {
      cat("\n\'max_ncomp\' reached. Stopping...\n")
      j.best <- max_ncomp
      fit_best <- fit_angmcmc
    }
  }

  result <- list("fit.all" = fit_all[1:j], "fit.best" = fit_best,
                 "ncomp.best" = all_ncomp[j.best], "crit" = crit,
                 "crit.all" = all_crit[1:j],
                 "crit.best" = all_crit[[j.best]],
                 "maxllik.all" = all_maxllik[1:j],
                 "maxllik.best" = all_maxllik[j.best],
                 "check_min" = check_min)
  class(result) <- "stepfit"

  result

}


#' Extracting angmcmc object corresponding to the best fitted model in stepwise fits
#'
#' @param step_object stepwise fitted object obtained from \link{fit_incremental_angmix}.
#'
#' @return Returns an angmcmc object corresponding to the the best fitted model in step_object.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' set.seed(1)
#' fit.vmsin.step.15 <- fit_incremental_angmix("vmsin", tim8, start_ncomp = 1,
#'                                             max_ncomp = 3, n.iter = 15,
#'                                             n.chains = 1)
#' fit.vmsin.best.15 <- bestmodel(fit.vmsin.step.15)
#' fit.vmsin.best.15
#'
#' @export

bestmodel <- function(step_object) {
  if(class(step_object) != "stepfit") stop("\'step_object\' is not a stepwise fitted object")
  step_object$fit.best
}
