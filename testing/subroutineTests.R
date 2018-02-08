### Title:    Subroutines for the MIBRR Package
### Author:   Kyle M. Lang
### Created:  2017-NOV-28
### Modified: 2018-FEB-07

##--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------##
##  Copyright (C) 2018 Kyle M. Lang <k.m.lang@uvt.nl>                          ##  
##                                                                             ##
##  This file is part of MIBRR.                                                ##
##                                                                             ##
##  This program is free software: you can redistribute it and/or modify it    ##
##  under the terms of the GNU General Public License as published by the      ##
##  Free Software Foundation, either version 3 of the License, or (at you      ##
##  option) any later version.                                                 ##
##                                                                             ##
##  This program is distributed in the hope that it will be useful, but        ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of                 ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General   ##
##  Public License for more details.                                           ##
##                                                                             ##
##  You should have received a copy of the GNU General Public License along    ##
##  with this program. If not, see <http://www.gnu.org/licenses/>.             ##
##-----------------------------------------------------------------------------##

#rm(list = ls(all = TRUE))

library(MIBRR)

load("../source/MIBRR/data/mibrrExampleData.RData")

#doBl           <- FALSE
#doImp          <- TRUE
#doMcem         <- TRUE
#data           <- mibrrExampleData
#nImps          <- 100
#targetVars     <- c("y", paste0("x", c(1 : 3)))
#ignoreVars     <- "idNum"
#iterations     <- c(30, 10)
#sampleSizes    <- list(rep(25, 2),
#                      rep(250, 2),
#                      rep(500, 2)
#                      )
#missCode       <- NA
#returnConvInfo <- TRUE
#returnParams   <- FALSE
#verbose        <- TRUE
#seed           <- NULL
#control        <- list()

dat1              <- mibrrExampleData
dat1[is.na(dat1)] <- -999

mibrrFit1 <- init(doBl           = FALSE,
                  doImp          = TRUE,
                  doMcem         = TRUE,
                  data           = dat1,
                  nImps          = 100,
                  targetVars     = paste0("x", c(1 : 3)), #c("y", paste0("x", c(1 : 3))),
                  ignoreVars     = "idNum",
                  iterations     = c(20, 5),
                  sampleSizes    = list(rep(25, 2),
                                        rep(250, 2),
                                        rep(500, 2)
                                        ),
                  missCode       = -999,
                  verbose        = TRUE,
                  seed           = NULL,
                  control        = list(center = FALSE)
                  )


tmp <- mibrrFit1$data
dat2 <- dat1[ , colnames(tmp)]

head(dat2)
head(tmp)

test <- dat2 == tmp

all.equal(!test, dat2 == -999)

head(mibrrFit1$data)
mibrrFit1$targetVars
mibrrFit1$missList
mibrrFit1$smoothingWindow

mibrrFit1 <- mcem(mibrrFit1)

mibrrFit1$gibbsOut$x1$imps
mibrrFit1$missList
mibrrFit1$data

    
mibrrFit1 <- postProcess(mibrrFit1)

mibrrFit1$gibbsOut$x1$imps

test <- mibrrFit1$getImpDataset()
test

impRowsPool <- mibrrFit1$impRowsPool
targetVars <- mibrrFit1$targetVars
gibbsOut <- mibrrFit1$gibbsOut

mibrrFit1$ignoredColumns
