# Change Log
All notable changes to the **mibrr** project will be documented in this file.

The format is based on [Keep a Changelog][kacl], and this project adheres to
[Semantic Versioning][sv].

NOTE: As of version 0.0.0.9000, I am returning to this project after
approximately 2.5 years away because I was unable to devote any attention to
this project during my postdoc. All initial development, therefore, is largely
undocumented.

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
- Moved optimization of Bayesian LASSO's penalty paramter into the R layer

### Fixed
- Bug causing failures of penalty parameters' optimization in MIBEN when 
  lambda = 0 at some iterations
- Bug triggered when missing data was coded as `NA`.

## 0.0.0.9001 - 2017-09-30

### Added
- *archive* directory to hold frozen source files
- froze the current *src* directory and stored it as "archive/frozen_src-20170930.tar.gz"
- froze the current *R* directory adn stored it as "archive/frozen_R-20170930.tar.gz"

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
