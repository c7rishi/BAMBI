# BAMBI v2.3.4

## Resubmission

This is a re-submission. In this version we have replaced all instances of '<<' redirect-based initialization of RcppArmadillo matrices or vectors with brace-initialization in the source C++ scripts. This follows the latest guidelines of RcppArmadillo, which notes that the '<<' redirect based initialization will be deprecated in a future release of the package.   


## Changelog:
* Replace '<<' redirect based initialization of armadillo vec/mat to brace-initializatiaon in c++ scripts


## Test environment:

* local Windows 10 install, R 4.2.1
* win-builder (devel, release, and old release)
* Ubuntu release 20.04, R 4.2.1
* experimental macbuilder at https://mac.r-project.org/

## R CMD check results

There were no ERRORs or WARNINGs in any platforms, and no NOTES on Windows - local or win-builder (devel, release, and old release).

There was a NOTE on Ubuntu 20.04 and on experimental macbuilder about the installed size of the package. Below is the NOTE from Ubuntu 20.04: 

installed size is  8.3Mb
  sub-directories of 1Mb or more:
    libs   7.8Mb



# BAMBI v2.3.3

Dear CRAN maintainers,

This re-submission aims to fully sync all references between the manuscript and the software. Our Journal of Statistical Software manuscript for BAMBI is now ready for publication, and we have recently received a DOI for the manuscript. The DOI displayed in the CITATION is for a new JSS publication that will be registered after publication on CRAN. We were unaware that we would be asked for the release of an updated package by JSS. We regret the need for making two concurrent submissions within the same week, and sincerely appreciate your understanding.     

## Resubmission
This is a resubmission. In this version we have made the following change.

## Changelog:
* Updated all citations to the manuscript for BAMBI with the DOI for the JSS article
* Created a new .Rd documentation for the package with a reference to the BAMBI manuscript.    

## Test environment:

* local Windows 10 install, R 4.1.1
* win-builder (devel, release, and old release)
* Ubuntu release 20.04, R 4.1.0

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE on the DOI (please see our comment above):
```
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.18637/jss.v099.i11
    From: DESCRIPTION
          man/BAMBI.Rd
    Status: 404
    Message: Not Found

Found the following (possibly) invalid DOIs:
  DOI: 10.18637/jss.v099.i11
    From: DESCRIPTION
          inst/CITATION
    Status: Not Found
    Message: 404
```

There was one NOTE on compiled code on R 4.1.1 on Windows 10 but not on win-builder:
```
> checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Note: information on .o files for x64 is not available
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  File 'BAMBI/libs/i386/BAMBI.dll':
    Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  File 'BAMBI/libs/x64/BAMBI.dll':
    Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  
  It is good practice to register native routines and to disable symbol
  search.
  
  See 'Writing portable packages' in the 'Writing R Extensions' manual.
```

There was one NOTE on R 4.1.0 on Ubuntu 20.04:
```
* checking installed package size ... NOTE
  installed size is  8.5Mb
  sub-directories of 1Mb or more:
    libs   7.9Mb

```



# BAMBI v2.3.2

## Resubmission
This is a resubmission. In this version we have made the following change.

## Changelog:
* Updated contour.angmcmc so that the ellipsis "..." arguments are properly passed to graphics::contour()

## Test environment:

* local Windows 10 install, R 4.1.1
* win-builder (devel, release, and old release)
* Ubuntu release 20.04, R 4.1.0

## R CMD check results
There were no ERRORs or WARNINGs.

There was a NOTE on compile code on R 4.1.1 on Windows 10 but not on win-builder:
```
> checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Note: information on .o files for x64 is not available
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  Warning in read_symbols_from_dll(so, rarch) :
    this requires 'objdump.exe' to be on the PATH
  File 'BAMBI/libs/i386/BAMBI.dll':
    Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  File 'BAMBI/libs/x64/BAMBI.dll':
    Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  
  It is good practice to register native routines and to disable symbol
  search.
  
  See 'Writing portable packages' in the 'Writing R Extensions' manual.
```

There was one NOTE on R 4.1.0 on Ubuntu 20.04:
```
* checking installed package size ... NOTE
  installed size is  8.5Mb
  sub-directories of 1Mb or more:
    libs   7.9Mb

```


# BAMBI v2.3.1

## Resubmission
This is a resubmission. In this version we have made the following changes.

## Changelog:
* Modified internal handling of floating points to ensure reproducibility across all platforms
* Updated loo::compare() (which is now deprecated) calls to loo::loo_compare()
* Introduced multiplicity adjustments to the sequential tests for incremental model fits.
* Minor bug fixes

## Test environments
* local Windows 10 install, R 4.1.0
* win-builder (devel and release)
* Ubuntu release 20.04, R 4.1.0

## R CMD check results
There were no ERRORs or WARNINGs.

There was a NOTE on compile code on R 4.1.0 on  Windows 10 but not on win-builder:
```
> checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Note: information on .o files for x64 is not available
  File 'D:/Documents/GitHub/BAMBI.Rcheck/BAMBI/libs/i386/BAMBI.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'D:/Documents/GitHub/BAMBI.Rcheck/BAMBI/libs/x64/BAMBI.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  
  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.
  
  See 'Writing portable packages' in the 'Writing R Extensions' manual. 
```

There was one NOTE on package size on R 4.1.0 on Ubuntu 20.04:

```
* checking installed package size ... NOTE 
  installed size is  8.5Mb
  sub-directories of 1Mb or more:
    libs   7.9Mb
```



# BAMBI v2.3.0
## Resubmission
This is a resubmission. In this version we have made the following changes.

## Changelog:
* added a function for maximum likelihood estimation of (single component) bivariate von Mises distribution
* fixed bug on computation of JS correlation coefficient 
* the default sampling algorithm from vmsin and vmcos are "naive" if n < 1e5
* removed dependency from rootSolve

## Test environments
* local Windows 10 install, R 4.0.0
* win-builder (devel and release)
* Ubuntu release 19.04, R 3.6.0

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 



# BAMBI v2.2.0

## Resubmission
This is a resubmission. In this version we have made the following changes.

### Changelog:
* added posterior predictive density evaluation in d_fitted and posterior predictive sampling in r_fitted
* fixed bugs that were creating some platform differences in fit_angmix
* added a check for NA values in data inside fit_angmix 

## Test environments
* local Windows 10 install, R 3.6.1
* win-builder (devel)
* Ubuntu release 19.04, R 3.6.0

## Test environments
* local Windows 10 install, R 3.5.3
* win-builder (devel)
* Ubuntu release 18.04, R 3.5.3

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 



# BAMBI v2.1.0

## Resubmission
This is a resubmission. In this version we have made the following changes.

### Changelog:
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



## Test environments
* local Windows 10 install, R 3.5.3
* win-builder (devel)
* Ubuntu release 18.04, R 3.5.3

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE on Ubuntu 18.04, but not in local Windows 10 nor in win-builder. The NOTE was on the size of installed package:

* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    libs   4.8Mb


# BAMBI v2.0.1

## Resubmission
This is a resubmission. In this version we have fixed the following NOTE:
"Namespace in Imports field not imported from: ‘future’"


## Test environments
* local Windows 10 install, R 3.5.0
* win-builder (devel)
* Ubuntu release 16.04, R 3.4.0


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
