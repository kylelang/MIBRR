# mibrr: Multiple Imputation with Bayesian Regularized Regression
This repository hosts development of the R package `mibrr`.

- Licensing information is given in the [LICENSE][] file.
- Built tarballs of the `mibrr` package are available in the [builds][] 
  directory.
- Stand-alone documentation is available in the [documentation][docs] directory.
- The source files for the most recent stable version of `mibrr` are available 
  in the [source][src] directory.

`mibrr` is alpha software, so please expect frequent---and dramatic---changes to 
the package's functionality and user interface. Please report any bugs that you 
encounter in the issues section of the project page. You may also leave requests 
for new features in the issues section.

Thank you for your interest in the MIBRR project! I hope you find my software
useful!

## Installation
The best way to install `mibrr` is to use the `devtools::install_github` 
function.

1. First, make sure that you have `devtools` installed on your system
2. Next, execute the following lines:

		library(devtools)
		install_github("kylelang/mibrr/source/mibrr")
    
3. Finally, load `mibrr` and enjoy:

		library(mibrr)

If the `devtools`-based approach does not work, you can download one of the
built tar-balls from the [builds][] directory and manually install the package
from source by executing the following lines:

	install.packages(pkgs  = "/SAVE_PATH/mibrr_VERSION.tar.gz",
	                 repos = NULL,
					 type  = "source")

Where *SAVE_PATH* is replaced by the (relative or absolute) file path to the
location where you saved the tar-ball, and *VERSION* is replaced with the correct
version number for the tar-ball that you downloaded.

## Examples

The `mibrr` package contains four primary functions: `miben`, `mibl`, `ben`, and 
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
		
- Use the follwing to fit models using `ben` or `bl`:

		## Load some data:
		data(mibrrExampleData)

		trainData <- mibrrExampleData[1 : 175, ]
		testData  <- mibrrExampleData[176 : 200, ]
		
		## Estimate a Bayesian elastic net model:
		benOut <- ben(data = trainData,
		              y    = "y",
					  X    = paste0("x", c(1 : 3))
					  )
		   
		## Estimate a Bayesian LASSO model:
		blOut <- bl(data = trainData,
		            y    = "y",
					X    = paste0("x", c(1 : 3))
					)

		## Generate out-of-sample predictions:
	    benPred <- 
			predictMibrr(object  = benOut,
			             newData = as.matrix(testData[ , paste0("x", c(1 : 3))])
						 )

		blPred <- 
			predictMibrr(object  = blOut,
			             newData = as.matrix(testData[ , paste0("x", c(1 : 3))])
						 )

[builds]:  https://github.com/kylelang/mibrr/tree/develop/builds/
[docs]:    https://github.com/kylelang/mibrr/tree/develop/documentation/
[src]:     https://github.com/kylelang/mibrr/tree/develop/source/mibrr
[LICENSE]: https://github.com/kylelang/mibrr/blob/develop/LICENSE
