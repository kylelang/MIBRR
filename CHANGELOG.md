# Change Log
All notable changes to the **MIBRR** project will be documented in this file.

The format is based on [Keep a Changelog][kacl], and this project (attempts to) 
adhere to [Semantic Versioning][sv].

NOTE: 

- As of `mibrr` version 0.0.0.9000, I am returning to this project after
  approximately 2.5 years away because I was unable to devote any attention to
  this project during my postdoc. All initial development, therefore, is largely
  undocumented.
- On 2017-11-06 the package name was changed from `mibrr` to `MIBRR`, so the 
  version number was reset to 0.0.0.9000, as well

## 0.0.0.9000 - XXXX-XX-XX

### Changed
- Changed package name from `mibrr` to `MIBRR`
- Reset version number

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
