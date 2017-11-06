# MIBRR Development ToDo List
## Kyle M. Lang
## Last Update: 2017-NOV-06

## Add a check for completely empty rows
- Check for fully missing rows, remove them, and add back after model fitting

## Tweak the implementation of `predictMibrr` - DONE (2017-NOV-06)
- Need to make the *newData* specification more robust
    - Accept a data.frame
	- Can rely on user to keep track of the predictor ordering when we're 
	  re-ordering the data columns before fitting
	
## Update docs to reflect new default for smoothingWindow - DONE (2017-NOV-03)

## Add missing docs - DONE (2017-NOV-03)

## Make sure we're dealing with missing data on predictors in a sensible way

- Zero imputation will only work with MCAR missingness
- Can we implement a reasonably efficient single imputation-based treatment?

## Find a good way to compute data centers with missing values

- Need to center data for the BEN to work, but center computed from incomplete 
  data will be biased.
- How can we get unbiased centers to use for correctly re-centering the imputed 
  data?

## Implement MCEM for MIBL - DONE (2017-OCT-25)

- Still need to add the deterministic update step for BL lambda

## Assess the current state of the package

- What works, what does not?

## Remove Nlopt - DONE (2017-OCT-01)

- Installing `NLopt` is proving to be a huge barrier to portability

## Replace NLopt with the initial R-based `optimx` optimization - DONE (2017-OCT-01)

- We may see a performance hit, but we'll probably fold the MCEM tasks into the 
  Gibbs sampler in the future, anyway
- Using `optimx` won't induce any portability problems

## Check that we didn't break anything when moving to the new way of indexing missing data - DONE (2017-OCT-27)

- The results are looking dubious--we need to find out why
    - The initial missing data treatment is a prime suspect
    - As is the new indexing method and potential problems with the 
	  post-processing that it may have induced
    - If we can get a stable version running, we need to start some simulations 
	  before doing any more development work.

## Remove intercept regularization option - DONE
## Use FIML to get means for centering - DONE (and partially reverted)

- FIML will be cumbersome/infeasible when P > N
    - I tried implementing an imputation-based solution
    - Don't know if this worked, may be a cause of the poor performance that we 
	  need to track down.
