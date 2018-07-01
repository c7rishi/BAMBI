# BAMBI v2.0.1
This is a resubmission. In this version we have fixed the following NOTE from R check:
"Namespace in Imports field not imported from: ‘future’"


# BAMBI v2.0.0

## Resubmission
This is a resubmission. A number of changes has been made, a summary of which are as follows:
* Changed the maintainer's mailing address
* Combined all fitting functions with known component sizes into one single funciton called 'fit_angmix'. The previous model specific fitting functions e.g. 'fit_vmsinmix' etc. are retained for backward compatibility, but they all now wrap to the same 'fit_angmix' function. The 'fit_angmix' function includes a number of new features, such as  multiple chains, an improved auto-tuning algorithm, permutation sampling, and option to impose restrictions on the component specific covariance parameters (such as fitting product mixtures).
* fit_stepwise_univariate and fit_stepwise_bivariate are now merged into fit_incremental_angmix.
* coda::as.mcmc.list method for mcmc objects is added.
* Added 'add_burnin_thin' as a post-processing function.
* Added 'bridge_sampler.angmcmc' method for bridgesampling::bridge_sampler
* WAIC computation is now done through the S3 method waic.angmcmc from package loo.
* Added loo.angmcmc function.
* fix_label now supports all methods for label.switching::label.switching
* added a plot method for angmcmc objects.
* added a naive rejection sampling method for vmsin and vmcos generation.
* Internal behaviors of a number of existing functions have been changed. 


## Test environments
* local Windows 10 install, R 3.5.0
* win-builder (devel)
* Red Hat Enterprise Linux Workstation release 6.8, R 3.4.2


## R CMD check results
There were no ERRORs or WARNINGs. 

There was one NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Saptarshi Chakraborty <chakra.saptarshi@gmail.com>'

New maintainer:
  Saptarshi Chakraborty <chakra.saptarshi@gmail.com>
Old maintainer(s):
  Saptarshi Chakraborty <c7rishi@ufl.edu>








# BAMBI v1.2.1

## Resubmission
This is a resubmission. In this version we have:
* Added two functions, one for evaluation of sample circular correlations for paired circular data, and one for evaluating population circular correlations and variances for vmsin, vmcos and wnorm2 models.
* Corrected a few typos in documentation. 
* Corrected Valgrind memory errors.

## Test environments
* local Windows 10 install, R 3.4.4
* win-builder (devel)
* Red Hat Enterprise Linux Workstation release 6.8, R 3.4.2

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE on Red Hat Enterprise Linux, but not in local Windows 10 nor in win-builder. The NOTE was on the size of installed package:

* checking installed package size ... NOTE
 installed size is  5.4Mb
  sub-directories of 1Mb or more:
    libs   5.0Mb


	
	
	
	
# BAMBI v1.1.1

## Resubmission
This is a resubmission. In this version we have:
* Added arXiv link to vignette in DESCRIPTION.
* Added native routine registration.
* Corrected a few typos in documentation. 


## Test environments
* local Windows 10 install, R 3.4.1
* Red Hat Enterprise Linux Workstation release 6.8, R 3.3.2
* win-builder (devel)


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE on local Windows, but not in Red Hat Enterprise Linux nor in win-builder. The NOTE was still there after adding native registration (using tools::package_native_routine_registration_skeleton):

* checking compiled code ... NOTE
File 'BAMBI/libs/x64/BAMBI.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

It is good practice to register native routines and to disable symbol
search.

See 'Writing portable packages' in the 'Writing R Extensions' manual.
