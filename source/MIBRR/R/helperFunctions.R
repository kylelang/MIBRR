### Title:    Helper Functions for MIBRR
### Author:   Kyle M. Lang
### Created:  2014-12-09
### Modified: 2019-12-11

##--------------------- COPYRIGHT & LICENSING INFORMATION --------------------##
##  Copyright (C) 2019 Kyle M. Lang <k.m.lang@uvt.nl>                         ##
##                                                                            ##
##  This file is part of MIBRR.                                               ##
##                                                                            ##
##  This program is free software: you can redistribute it and/or modify it   ##
##  under the terms of the GNU General Public License as published by the     ##
##  Free Software Foundation, either version 3 of the License, or (at you     ##
##  option) any later version.                                                ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful, but       ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of                ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General  ##
##  Public License for more details.                                          ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License along   ##
##  with this program. If not, see <http://www.gnu.org/licenses/>.            ##
##----------------------------------------------------------------------------##

## Print startup message:
.onAttach <- function(libname, pkgname) {
    tmp <- read.dcf(file = system.file("DESCRIPTION", package = pkgname))
    
    bDate <- gsub(".*(\\d{4}-\\d{2}-\\d{2}\\s\\d{1,2}:\\d{2}:\\d{2}\\sUTC).*",
                  "\\1",
                  tmp[ , "Built"])
    
    author <- gsub("\\s\\[.*\\]", "", tmp[ , "Author"])

    greet0 <- c(paste0("Package:   ", pkgname),
                paste0("Version:   ", tmp[ , "Version"]),
                paste0("Built:     ", bDate),
                paste0("Copyright: ",
                       format(Sys.time(), "%Y"),
                       " (",
                       author,
                       ")\n")
                )
    
    greet <-
        c(greet0,
          strwrap(
              paste0(pkgname,
                     " is distributed under Version 3 of the GNU General Public License (GPL-3); execute 'mibrrL()' for details. ",
                     pkgname,
                     " comes with ABSOLUTELY NO WARRANTY; execute 'mibrrW()' for details. ",
                     pkgname,
                     " is beta software. Please report any bugs to ",
                     tmp[ , "BugReports"],
                     "."),
              width = 81)
          )
    
    for(i in greet) packageStartupMessage(i)
}

###--------------------------------------------------------------------------###

## Print 'x' only if 'verbose = TRUE':
vcat <- function(x) if(parent.frame()$mibrrFit$verbose) cat(x)

###--------------------------------------------------------------------------###

## Estimate the mode of a continous vector:
numMode <- function(x) {
    dens <- density(x, na.rm = TRUE)
    dens$x[which.max(dens$y)]
}

###--------------------------------------------------------------------------###

## Cast object to given type:
cast <- function(obj, type)
    eval(call(paste0("as.", type), obj))

###--------------------------------------------------------------------------###

## Set the control list arguments for a particular class:
setControl <- function(x, where) {
    ## Get the fields for the current class:
    fields <- getRefClass(class(where))$fields()
    
    ## Assign the control list entries to the correct classes:
    for(n in names(x))
        if(n %in% names(fields))
            where$field(n, cast(x[n], fields[n]))
}
