# MIBRR Development Notes
## Kyle M. Lang
## Last Updated: 2017-SEP-30

# 2016-NOV-05:
## The most recent branch is "moreControl"

## MibrrData2 returns predictor matrices that do not include a constant columns

## MibbrGibbs2 treats the intercept and slope components seperatly

## MibrrGibbs absorbs the intercept into the predictor matrix and avoid regularizing by setting _tau[0, 0] = 0

# 2016-NOV-10:
## Question: If the observed parts of the data are centered, is it ever possible to estimate a non-zero intercept in an imputation model?

## "classic mode" calls version 2 of all the class definitions

- This version appears to regularize the intercept

## Version 1 of the class definitions offers the option to regularize the intercepts

