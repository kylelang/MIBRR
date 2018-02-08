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

rm(list = ls(all = TRUE))

library(MIBRR)

load("../source/MIBRR/data/mibrrExampleData.RData")

dat1              <- mibrrExampleData
dat1[is.na(dat1)] <- -999

mibrrFit1 <- MIBRR::::init(doBl           = FALSE,
                           doImp          = TRUE,
                           doMcem         = TRUE,
                           data           = dat1,
                           nImps          = 100,
                           targetVars     = paste0("x", c(1 : 3)),
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
mibrrFit1 <- MIBRR:::mcem(mibrrFit1)
mibrrFit1 <- MIBRR:::postProcess(mibrrFit1)
