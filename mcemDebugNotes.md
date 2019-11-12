# MCEM Debugging Notes

## Last Updated: 2019-11-12

## Problem:
- The MCEM implementation seems to be underidentified
- Chains started at seperate values will converge to stable, but unique, equilibria

## ToDo:
- Isolate where this problem occurs:
  - MIBEN
	- dissRerun simulation: pm = 10, exp1, n = 100, dense
	- dissRerun simulation: pm = 10, exp1, n = 100, sparse
  - MIBL
	- dissRerun simulation: pm = 10, exp1, n = 100, dense
	- dissRerun simulation: pm = 10, exp1, n = 100, sparse
  - BEN
  - BL
  
- See if using different optimization algorithms helps
  - Unlikely since both MIBEN and MIBL are affected
