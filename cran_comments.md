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