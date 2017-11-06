# MIBRR: Multiple Imputation with Bayesian Regularized Regression
This repository hosts development of the R package `MIBRR`.

- Licensing information is given in the [LICENSE][] file.
- Built tarballs of the `MIBRR` package are available in the [builds][] 
  directory.
- Stand-alone documentation is available in the [documentation][docs] directory.
- The source files for the most recent stable version of `MIBRR` are available 
  in the [source][src] directory.

`MIBRR` is alpha software, so please expect frequent---and dramatic---changes to 
the package's functionality and user interface. Please report any bugs that you 
encounter in the issues section of the project page. You may also leave requests 
for new features in the issues section.

Thank you for your interest in the MIBRR project! I hope you find my software
useful!

## Installation
The best way to install `MIBRR` is to use the `devtools::install_github` 
function.

1. First, make sure that you have `devtools` installed on your system
2. Next, execute the following lines:

		library(devtools)
		install_github("kylelang/MIBRR/source/MIBRR")
    
3. Finally, load `MIBRR` and enjoy:

		library(MIBRR)

If the `devtools`-based approach does not work, you can download one of the
built tar-balls from the [builds][] directory and manually install the package
from source by executing the following lines:

	install.packages(pkgs  = "/save_path/MIBRR_version.tar.gz",
	                 repos = NULL,
                     type  = "source")

Where *save_path* is replaced by the (relative or absolute) file path to the
location where you saved the tar-ball, and *version* is replaced with the correct
version number for the tar-ball that you downloaded.

## Examples

The `MIBRR` package contains four primary functions: `miben`, `mibl`, `ben`, and 
`bl`. 

- The `miben` and `mibl` functions do multiple imputation using the Bayesian 
  elastic net and Bayesian LASSO, respectively. 
- The `ben` and `bl` functions fit Bayesian elastic net and Bayesian LASSO
  models to incomplete data without returning and imputed datasets.
- Basic missing data treatments using `miben` or `mibl` might look like the 
  following:

		## Load some data:
		data(mibrrExampleData)

		## Create M = 100 multiply imputed datasets:
		mibenOut <- miben(data       = mibrrExampleData,
                          nImps      = 100,
                          targetVars = c("y", paste0("x", c(1 : 3))),
                          ignoreVars = "idNum")
			  
		miblOut <- mibl(data       = mibrrExampleData,
                        nImps      = 100,
                        targetVars = c("y", paste0("x", c(1 : 3))),
                        ignoreVars = "idNum")
				
		## Extract list of imputed datasets:
		mibenImps <- mibenOut$imps
		miblImps  <- miblOut$imps
		
- Use the following to fit models using `ben` or `bl`:

		## Load some data:
		data(predictData)

		trainData <- predictData$train
		testData  <- predictData$test
		
		## Estimate a Bayesian elastic net model:
		benOut <- ben(data = trainData,
                      y    = "agree",
                      X    = setdiff(colnames(trainData), "agree")
                      )
		   
		## Estimate a Bayesian LASSO model:
		blOut <- bl(data = trainData,
                    y    = "agree",
                    X    = setdiff(colnames(trainData), "agree")
                    )

		## Generate out-of-sample predictions:
	    benPred <- predictMibrr(object = benOut, newData = testData)
		blPred  <- predictMibrr(object = blOut, newData = testData)
		
- Posterior predictions can also be generated from `miben` and `mibl` models:

		## Load some data:
		data(predictData)

		missData <- predictData$incomplete
		testData <- predictData$test
		
		## Create M = 100 multiply imputed datasets:
		mibenOut <- miben(data = missData, nImps = 100)
		miblOut  <- mibl(data = missData, nImps  = 100)
		
		## Generate out-of-sample predictions:
	    mibenPred <- predictMibrr(object = mibenOut, newData = testData)
		miblPred  <- predictMibrr(object = miblOut, newData = testData)
		
		
[builds]:  https://github.com/kylelang/MIBRR/tree/master/builds/
[docs]:    https://github.com/kylelang/MIBRR/tree/master/documentation/
[src]:     https://github.com/kylelang/MIBRR/tree/master/source/MIBRR
[LICENSE]: https://github.com/kylelang/MIBRR/blob/master/LICENSE
