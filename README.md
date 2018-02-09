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
The best way to install (the development version of) `MIBRR` is to use the 
`devtools::install_github` function.

1. First, make sure that you have `devtools` installed on your system
2. Next, execute the following lines:

		library(devtools)
		install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")
    
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

- The underlying model in each of these four primary functions can be estimated
  using either Markov Chain Expectation Maximization (MCEM) or fully Bayesian
  modeling.
- The `miben` and `mibl` functions do multiple imputation using the Bayesian 
  elastic net and Bayesian LASSO, respectively.
  
    - A list of imputed datasets can be generated using the `complete` function.
	
- Basic missing data treatments using `miben` or `mibl` might look like the 
  following:

		## Load some data:
		data(mibrrExampleData)

		## Estimate the imputation models using MCEM:
		mibenOut <- miben(data       = mibrrExampleData,
                          targetVars = c("y", paste0("x", c(1 : 3))),
                          ignoreVars = "idNum")
			  
		miblOut <- mibl(data       = mibrrExampleData,
                        targetVars = c("y", paste0("x", c(1 : 3))),
                        ignoreVars = "idNum")
						
		## Estimate the imputation models using fully Bayesian modeling:
		mibenOut <- miben(data         = mibrrExampleData,
		                  targetVars   = c("y", paste0("x", c(1 : 3))),
                          ignoreVars   = "idNum",
                          doMcem       = FALSE,
                          sampleSizes  = c(500, 500),
                          lam1PriorPar = c(1.0, 0.1),
                          lam2PriorPar = c(1.0, 0.1)
                          )
			  
		miblOut <- mibl(data         = mibrrExampleData,
		                targetVars   = c("y", paste0("x", c(1 : 3))),
                        ignoreVars   = "idNum",
                        doMcem       = FALSE,
                        sampleSizes  = c(500, 500),
                        lam1PriorPar = c(1.0, 0.1)
                        )
				
		## Extract list of 100 imputed datasets:
		mibenImps <- complete(mibrrFit = mibenOut, nImps = 100)
		miblImps  <- complete(mibrrFit = miblOut, nImps = 100)
		
- The `ben` and `bl` functions fit Bayesian elastic net and Bayesian LASSO
  models to incomplete data without returning and imputed datasets.
- Use the following to fit models using `ben` or `bl`:

		## Load some data:
		data(predictData)

		trainData <- predictData$train
		testData  <- predictData$test
		
		## Estimate a Bayesian elastic net model using MCEM:
		benOut <- ben(data = trainData,
                      y    = "agree",
                      X    = setdiff(colnames(trainData), "agree")
                      )
		   
		## Estimate a Bayesian LASSO model using MCEM:
		blOut <- bl(data = trainData,
                    y    = "agree",
                    X    = setdiff(colnames(trainData), "agree")
                    )

		## Estimate a Bayesian elastic net model using full Bayes:
		benOut <- ben(data         = trainData,
                      y            = "agree",
                      X            = setdiff(colnames(trainData), "agree"),
                      doMcem       = FALSE,
                      sampleSizes  = c(500, 500),
                      lam1PriorPar = c(1.0, 0.1),
                      lam2PriorPar = c(1.0, 0.1)
                      )
		   
		## Estimate a Bayesian LASSO model using full Bayes:
		blOut <- bl(data         = trainData,
                    y            = "agree",
                    X            = setdiff(colnames(trainData), "agree"),
                    doMcem       = FALSE,
                    sampleSizes  = c(500, 500),
                    lam1PriorPar = c(1.0, 0.1)
                    )
					
		## Extract posterior parameter samples:
		benPars <- getParams(mibrrFit = benOut, target = "agree")
		blPars  <- getParams(mibrrFit = blOut, target = "agree")
		
		## Generate out-of-sample predictions:
	    benPred <- predictMibrr(mibrrFit = benOut, newData = testData)
		blPred  <- predictMibrr(mibrrFit = blOut, newData = testData)
		
- Posterior predictions can also be generated from `miben` and `mibl` models:

		## Load some data:
		data(predictData)

		missData <- predictData$incomplete
		testData <- predictData$test
		
		## Estimate the imputation models:
		mibenOut <- miben(data = missData)
		miblOut  <- mibl(data = missData)
		
		## Generate out-of-sample predictions:
	    mibenPred <- predictMibrr(mibrrFit = mibenOut, newData = testData)
		miblPred  <- predictMibrr(mibrrFit = miblOut, newData = testData)
		
		
[builds]:  https://github.com/kylelang/MIBRR/tree/develop/builds/
[docs]:    https://github.com/kylelang/MIBRR/tree/develop/documentation/
[src]:     https://github.com/kylelang/MIBRR/tree/develop/source/MIBRR
[LICENSE]: https://github.com/kylelang/MIBRR/blob/develop/LICENSE
