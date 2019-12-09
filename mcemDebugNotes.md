# MCEM Debugging Notes

## Last Updated: 2019-12-03

## Problem:
- The MCEM implementation seems to be underidentified
- Chains started at separate values will converge to stable, but unique, 
  equilibria

## ToDo:
- Isolate where this problem occurs:
  - YES:
	- MIBEN
	  - dissRerun simulation: pm = 10, exp1, n = 100, dense
	  - dissRerun simulation: pm = 10, exp1, n = 100, sparse
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0, nTargets = 4 (just a 
		little)
	  - dissRerun simulation: pm = 10, exp1, n = 100, sparse, zCor = 0.3
	  - dissRerun simulation: pm = 10, exp1, n = 100, sparse, zCor = 0.5
	- MIBL
	  - dissRerun simulation: pm = 10, exp1, n = 100, dense
	  - dissRerun simulation: pm = 10, exp1, n = 100, sparse
	  - simple mvn data: n = 100, p = 3, pm = 30, cor = 0.3, nTargets = 2
	  - simple mvn data: n = 100, p = 3, pm = 50, cor = 0.3, nTargets = 2
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0.3, nTargets
        = 4 (but not as bad as the p = 3 version)
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0.0, nTargets = 4
	  - dissRerun simulation: pm = 10, exp1, n = 100, sparse, zCor = 0.3
	  - dissRerun simulation: pm = 10, exp1, n = 100, sparse, zCor = 0.5
	- BEN
	- BL
  - NO:
	- MIBEN
	  - simple mvn data: n = 100, p = 3, pm = 10, cor = 0.3, nTargets = 2
	  - simple mvn data: n = 100, p = 3, pm = 30, cor = 0.3, nTargets = 2
	  - simple mvn data: n = 100, p = 3, pm = 50, cor = 0.3, nTargets = 2
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0.3, nTargets = 4
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0.5, nTargets = 4
	- MIBL
	  - simple mvn data: n = 100, p = 3, pm = 10, cor = 0.3, nTargets = 2
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0.5, nTargets = 4
	- BEN
	- BL
  
- See if using different optimization algorithms helps
  - Unlikely since both MIBEN and MIBL are affected

- See if chains started at the same value diverge.
  - Doesn't look like it. If we start the simple mvn case at the same values, 
	the chains will track identical trajectories.
	
- Check if the mismatch between the Park and Casella (2008) lambda estimate and
  mine is due to my intercept.
  
## Observations:
- When the MIBL chains don't mix, the sampled lambdas tend to be huge (about an
  order of magnitude larger than the well-behaved chains).
- Increasing the size of the burn-in Gibbs samples (i.e., n = 25, 50, 250),
  doesn't help the MIBL chains in the simple mvn data case.
- **Low correlations seem to cause problems**
  - Lambdas cannot be optimized in BEN
  - Lambdas explode and/or are not identified in BL
  - With BEN/BL, in the simple MVN case, r = 0.2 seems to be sufficient
- The following apply to a sparse version of the simple MVN case
  - Setup: 5 independent variables, 5 correlated at r = 0.5
  - With N = 25 for approximation iterations, 13 out of 15 reps fail
  - With N = 250 for approximation iterations, 9 out of 15 reps fail
  - With N = 500 for approximation iterations, 9 out of 15 reps fail
    (performance doesn't seem much better than N = 250).
  - With N = 1000 for approximation iterations, performance is pretty much the 
	same as N = 250 and N = 500.
  - Large differences between the trivial and non-trivial effect sizes seem to
    cause problems.
	- Setting trivial effects to r = 0 and non-trivial to r = 0.5 leads to poor
      convergence rates
	- Setting trivial effects to r = 0 and non-trivial to r = 0.25 produces
      perfect convergence and no identification issues
	- Setting trivial effects to r = 0 and non-trivial to r = 0.15 is too
      much. We get moderate convergence problems and and identification issues
  - When most effects are trivial (i.e., p < 3 non-trivial predictors out of
    10), things work well
  - When all all effects are non-trivial, things work well.
  - When 3 to 9 effects (out of 10) are non-trivial, problems occur
- I can mostly replicate the analysis of the *diabetes* data from Park & Casella 
  (2008).
  - The MCEM chains converge and mix nicely
  - I get the same posterior medians of beta
  - I get the same credible intervals for beta
  - I get the same L1 norm of beta relative to the least squares estimates
  - **BUT** my estimate of lambda is approximately 20 times larger than theirs
	- I get lambda ~ 5, they get lambda ~ 0.237
  - Using the Park & Casella (2008) prior parameterization, I get very different
    results
	- If I set the prior to exponential with mean of 10 times the ML estimate
      (as done in Park and Casella, 2008), I see the following relative to my ML
      results:
	  - About the same posterior medians of beta
	  - About the same credible intervals for beta
	  - About the same (relative) L1 norm for beta
	  - Different (much smaller) posterior median of lambda
- I can replicate the BEN and BL results for the *diabetes* data from Li and Lin
  (2010), if I scale the predictor data.
  - If the predictor data is not scaled, the optimization fails.
  - The predictors should be scaled internally, so this seems odd.

## Fixes
- The standardization of the beta samples was being reverted before they were
  stored.
  - This was done to return the posterior samples on the same scale as the raw
    data (since standardization is done on-the-fly).
  - However, this reversion was done before computing the MC expectations used
    to parameterize the loglikelihood function for lambda.
  - It seems that undoing the standardization before optimizing lambda was
    causing some serious issues, including:
	1. Convergence failures in the optimization

## Ideas
- The MCEM algorithm produces a Markov chain of lambdas. Maybe we should use 
  these chains in a more traditional way.
  - Right now, we're taking the final estimate of lambda and using that to 
	parameterize our posteriour samples
  - Why don't we average over a stationary chunk of the chain to get our final 
	estimate of lambda?
  - We'll need to check the chains for convergence, first.
	- Maybe just use two split chains to compute and R-Hat?
