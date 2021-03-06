# BAMBI 2.1.0
* fixed bugs that were causing some examples to fail 
(due to an update in package loo)
* updated documentation of dwnorm, dwnorm2, fit_angmix and
fit_incremental_angmix
* added an optional argument 'force_approx_const' to dvmcos
* added a function 'bestcriterion' to extract the value of best model selection criterion from stepwise fits.
* densityplot now uses lattice::wireframe instead of persp
* added alpha shading argument in contour.angmcmc.
* removed press.enter argument from plot.angmcmc, paramtrace and lpdtrace.
* changed BESSI, BESSI0, BESSI0_C and BESSI1. They're now all wrappers for R::bessel_i.
* changed "n"-cluster to "n"-component in print.angmcmc.
* Deleted AIC.angmcmc, BIC.angmcmc -- they can be directly computed using stats::AIC and stats::BIC.



# BAMBI 2.0.0
* Combined all fitting functions with known component sizes into one single funciton called 'fit_angmix'. The previous model specific fitting functions e.g. 'fit_vmsinmix' etc. are retained for backward compatibility, but they all now wrap to the same 'fit_angmix' function. The 'fit_angmix' function includes a number of new features, such as  multiple chains, an improved auto-tuning algorithm, permutation sampling, and option to impose restrictions on the component specific covariance parameters (such as fitting product mixtures).
* fit_stepwise_univariate and fit_stepwise_bivariate are now merged into fit_incremental_angmix.
* coda::as.mcmc.list method for mcmc objects is added.
* Added 'add_burnin_thin' as a post-processing function.
* Added 'bridge_sampler.angmcmc' method for bridgesampling::bridge_sampler
* WAIC computation is now done through the S3 method waic.angmcmc from package loo.
* Added loo.angmcmc function.
* fix_label now supports all methods for label.switching::label.switching
* A plot method for angmcmc objects is addedd.
* A naive rejection sampling method for vmsin and vmcos generation is added.


# BAMBI 1.2.1
* Added two functions, one for evaluation of sample circular correlations for paired circular data, and one for evaluating population circular correlations and variances for vmsin, vmcos and wnorm2 models.
* Corrected a few typos in documentation. 
* Corrected valgrind memory errors.


# BAMBI 1.1.1
* Added arXiv link to vignette in DESCRIPTION.
* Added native routine registration.
* Corrected a few typos in documentation. 


# BAMBI 1.1.0
* Added option to choose hyper parameters in the Dirichlet prior for mixing proportion.  
* Added a function to create logLik objects from angmcmc objects.
* Added a function to extract the best fitted model object from step-wise fits.


# BAMBI 1.0.1
* Corrected bugs that caused build to fail with clang, MacOS, and Solaris.
