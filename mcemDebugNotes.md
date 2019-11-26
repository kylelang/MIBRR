# MCEM Debugging Notes

## Last Updated: 2019-11-18

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
	  - simple mvn data: n = 100, p = 14, pm = 30, cor = 0, nTargets = 4 (just a little)
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
	
## Observations:
- When the MIBL chains don't mix, the sampled lambdas tend to be huge (about an
  order of magnitude larger than the well-behaved chains).
- Increasing the size of the burn-in Gibbs samples (i.e., n = 25, 50, 250),
  doesn't help the MIBL chains in the simple mvn data case.
- Low correlations seem to cause problems
  - Lambdas cannot be optimized in BEN
  - Lambdas explode and/or are not identified in BL
  - With BEN/BL, in the simple MVN case, r = 0.2 seems to be sufficient
