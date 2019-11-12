# Change Log
All notable changes to the **MIBRR** project (and some not so notable ones) will 
be documented in this file.

The format is based on [Keep a Changelog][kacl], and this project (attempts to) 
adhere to [Semantic Versioning][sv].

NOTE: 

- As of `mibrr` version 0.0.0.9000, I am returning to this project after
  approximately 2.5 years away because I was unable to devote any attention to
  this project during my postdoc. All initial development, therefore, is largely
  undocumented.
- On 2017-11-06 the package name was changed from `mibrr` to `MIBRR`, so the 
  version number was reset to 0.0.0.9000, as well

## 0.3.2.9000 - 2019-03-13 - ACTIVE

### Fixed
- Disabled multithreading in Eigen. We don't want Eigen to autonomously spawn
  its own threads. Doing so leads to unexpected overcommitting of CPU resources
  (especially when running in a cluster environment).
  
### Changed
- Improved the start-up message.

### Added
- A new testing script *visPriors.R* to generate plots visualizing the penalty 
  parameters' prior distributions. This functionality will be incorporated into 
  an exported function in the future.
  
- The option to restart failed optimizations of the penalty parameters in 
  `miben` with randomly perturbed starting values.
	  
	  - The number of restarts is set by the *optMaxRestarts* in the control 
	    list.

- A new control list option, *optStrict*, that dictates whether *optMaxRestarts* 
  failed optimizations of the penalty parameters results in a fatal error or in 
  a warning.
    
## 0.3.1.9000 - 2019-02-26

### Fixed
- A bug with the implementation of the `vanilla` imputation method.

### Changed
- Cleaned up the C++ source code a little bit.

### Added
- A control list option *useBetaMeans* that will tell the Gibbs sampler to use 
  the OLS estimates of *Beta* when updating *sigma* (as opposed to the default 
  option of using *Beta*s most recent sampled values).
- A function, `ppCheck`, to generate overlaid density plots of the observed and 
  posterior predicted samples in a *MibrrFit* object.
  
## 0.3.0.9000 - 2019-01-25

### Changed
- Removed all dependencies except for **mice**; the necessary functions are now
  imported from the previous dependencies.
- Updated the testing script: */testing/testPackage.R* to systematize the 
  testing implemented therein.

### Added
- A new R syntax file: *unitTests.R*.
- Unit tests for the random number samplers (`testMvn`, `testInvGamma`, 
  `testInvGauss`, `testGig`, `testInvChiSq`, `testIncGamma`, and `testSamplers`; 
  not currently exported).
- Unit test for the missing data indexing (`testMissIndex`; not currently 
  exported).
- Unit test for data manipulation/subsetting (`testDataProcessing`; not 
  currently exported).
- Unit test for data scaling (`testDataScaling`; not currently exported).
- Unit test for missing data filling (`testMissFill`; not currently exported).

## 0.2.0.9000 - 2019-01-21

### Fixed
- The scaling issues addressed in the changes listed below seem to have been a 
  large part of the poor performance exhibited in recent simulations. The full 
  extent of the improvement still needs to be evaluated via new simulations.
  
### Changed
- Beta samples are now returned in their raw metric
- Predictors are standardized on-the-fly with respect to the "training set" 
  used to estimate each elementary imputation model
- Outcomes are left in their raw metric (uncentered) during model estimation
- The "test sets" used to generate the imputations are standardized with 
  respect to the "training set" moments

## 0.1.0.9000 - 2019-01-15 

### Changed
- Updated version number

### NOTE
- This version shall be considered a "stable" reference version of the original 
  implementation of `MIBRR` that was used in Lang (2015). 
- This version is STILL BROKEN
- This version shall act as a known starting point from which I can fix the 
  broken implementation.
  
## 0.0.0.9007 - 2018-11-21

### Fixed
- Two bugs causing crashes when trying to impute univariate missing data

### Added
- A new exported function, `vanilla`, that implements a basic MI without 
  using any shrinkage priors.
- The option to specify known means and standard deviations for the data. This 
  option is mostly for debugging purposes.
  
## 0.0.0.9006 - 2018-06-07

### Fixed
- Bug causing crashes when initially imputing incomplete auxiliary variables 
  with mice.
- Removed redundant definition of `simpleImpute` function from 
  'helperFunctions.R'
- Removed unused/redundant definition of `smoothLambda` member function.

### Changed
- Improved the way random numbers are generated. Each subprocess now gets its 
  own R'Lecuyer RNG stream.
- Removed the `simpleIntercept` option from the control list.
  
### Added
- Option to specify an active `rlecuyer` RNG stream for the current session.
  This stream will be re-set as the active RNG stream when `MIBRR` returns.
  
## 0.0.0.9005 - 2018-05-04

### Fixed
- Bug keeping user from changing optimization algorithm used for MCEM
- Bug keeping 'lambdaHistory' from filling when using `mibl`.
- Added the correct 'LICENSE' file to the repository root directory

### Changed
- Removed restriction that forced 'L-BFGS-B' as the optimization method when 
  doing constrained optimization for MCEM
- Improved implementation of 'mibrrW' and 'mibrrL' functions.

## 0.0.0.9004 - 2018-02-15

### Fixed
- Removed some unused variables from the C++ source.
- Bug causing crashes when a Sigma's R-Hat fails convergence checks.

## 0.0.0.9003 - 2018-02-13

### Changed
- Updated all documentation files (they still need work, though).
- Updated the development branch README file.
- Improved the way random numbers are passed to C++ to seed the samplers class.
- Header file extension changed from '.hpp' to '.h' to appease CRAN checks.
- Cleaning up doc files
- `predictMibrr` renamed to `postPredict`
- Cleaned up the `initialize()` member function and implemented it in a less 
  stupid way

### Added 
- Fully Bayesian estimation of the Elastic Net and LASSO penalty parameters
- A `getImpData` function to replace missing values and generate the multiply
  imputed datasets after a run of `miben` or `mibl`.
- A `getParams` function to extract the posterior samples of the model parameters
  from a fitted `MibrrFit` object.
- An `getField` function to extract arbitrary fields from a `MibrrFit` object.
- Documentation for `getImpData`, `getField`, `getParams` functions and for the 
  `MibrrFit` class.
- Options to generate predictions for MAP scores and EAP scores in `postPredict`

## 0.0.0.9002 - 2018-02-08

### Changed
- Removed the `mibrr` function. The new `init`, `mcem`, and `postProcess` 
  subroutines are directly wrapped by `miben`, `mibl`, `ben`, `bl`.
- Pulled all main exported functions into new source file: 
  "exportedPrimaryFunctions.R"
- Broke core estimation subroutines into a new source files "subroutines.R"
- Pulled all samplers into a separate class `MibrrSamplers`

    - `MibrrGibbs` now extends `MibrrSamplers`

- Moved many R functions into the `MibrrFit` class as methods
- Broke optimization/Gibbs sampling routines out into a separate source file: 
  "02_EstimationMethods.R"
- Broke exported helper functions into a separate source file: 
  "exportedHelperFunctions.R"
	  
### Added
- A generalized inverse Gaussian sampler
- A metadata Reference Class: `MibrrFit`

    - Objects returned by `MIBRR` calls will now have class `MibrrFit`
	
## 0.0.0.9001 - 2018-02-06

### Added
- *Makevars.win* file allowing C++11 support on Windows
- Help file for `predictData` dataset
- Appropriate `importFrom` statements to pull functions from `mice`, `lavaan`, 
  `optimx`, `stats`, and `mvtnorm`

### Fixed
- Edited various message and comment text to remove references to LGPL-3 (now 
  GPL-3).

## 0.0.0.9000 - 2017-11-09

### Changed
- Changed package name from `mibrr` to `MIBRR`
- Reset version number
- Cleaned up the */testing* directory
- Made GitHub repository public

### Added
- A subdirectory: */builds/oldName* to hold the builds generated under the old 
  naming scheme

## RE-NAME and VERSION RESET

## 0.0.0.9007 - 2017-11-06

### Added
- New example data: `predictData`

### Changed
- Improved the `predictMibrr` function
- Updated documentation examples

## 0.0.0.9006 - 2017-11-03

### Added
- Help files for `predictMibrr`, `mibrrW`, and `mibrrL`

### Changed
- Updated docs to accurately reflect new default for the `smoothingWindow` 
  control parameter
  
## 0.0.0.9005 - 2017-11-02

### Changed
- Tweaked the iteration output printed to stdout
- Updated the README.md file
- Renamed 'testingSupportFunctions.R' to 'testingFunctions.R'
- Renamed 'mibrrHelperFunction.R' to 'helperFunctions.R'
- Renamed 'extraFunctions.cpp' to 'testingFunctions.cpp'
- Removed debugging/testing functions from the exports
- Updated default for 'smoothingWindow' control parameter

### Fixed
- A bug triggered by malformed lambda matrices in models with only 1 DV
- Added back support for lambda smoothing windows that was inadvertently lost 
  when moving optimization to the R layer

## 0.0.0.9004 - 2017-11-01

### Added
- Included a "documentation" directory
- Put the .pdf package doc in the "documentation" directory
- Help files for the `ben` and `bl` functions

### Changed
- Updated documentation
- Changed license from LGPL-3 to GPL-3 since `mibrr` no longer uses `nlopt` 

    - `nlopt` is distributed under the LGPL, which precluded releasing `mibrr` 
	  under the GPL.
	   
- Tweaked the implementation of `ben` and `bl` functions.

    - Both functions now work with incomplete data
	- Both functions only allow one DV per call
	
### Fixed
- Bug caused by trying to impute fully observed covariate matrices

## 0.0.0.9003 - 2017-10-27

### Changed
- Updated copyright holder's email address in the source code

### Added
- Exported functions to examine missing data indexing
- Confirmed that the missing data are being indexed correctly

## 0.0.0.9002 - 2017-10-25

### Changed
- Moved optimization of Bayesian LASSO's penalty parameter into the R layer

### Fixed
- Bug causing failures of penalty parameters' optimization in MIBEN when 
  lambda = 0 at some iterations
- Bug triggered when missing data was coded as `NA`.

## 0.0.0.9001 - 2017-09-30

### Added
- *archive* directory to hold frozen source files
- froze the current *src* directory and stored it as "archive/frozen_src-20170930.tar.gz"
- froze the current *R* directory and stored it as "archive/frozen_R-20170930.tar.gz"

### Changed
- Removed `NLopt`-based optimization for the MCEM steps
- Replaced the `NLopt`-based optimization with an `optimx`-based solution in the R layer
- Use a different way to index various Gibbs sample sizes

## 0.0.0.9000 - 2017-09-30

### Changed
- Adopted semantic versioning scheme
- Reset version number
- Cleaned up/restructured the project directory

### Added
- This changelog file (CHANGELOG.md).

[kacl]: http://keepachangelog.com/
[sv]:   http://semver.org/
[hw]:   http://r-pkgs.had.co.nz/
