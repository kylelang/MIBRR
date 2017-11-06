### Title:    Prepare Testing Data for MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2017-NOV-06

rm(list = ls(all = TRUE))
set.seed(235711)

library(psych)
library(devtools)

install_github("kylelang/SURF/source/SURF")
library(SURF)

data(bfi)
tmp <- na.omit(bfi)

keys <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
             consc = c("C1", "C2", "C3", "-C4", "-C5"),
             extra = c("-E1", "-E2", "E3", "E4", "E5"),
             neuro = c("N1", "N2", "N3", "N4", "N5"),
             open  = c("O1", "-O2", "O3", "O4", "-O5")
             )

scores <- scoreItems(keys = keys, items = tmp)$scores

ed.d           <- model.matrix(~factor(tmp$education))[ , -1]
colnames(ed.d) <-
    c("finish_hs", "some_college", "college_grad", "graduate_degree")

male            <- tmp$gender
male[male == 2] <- 0

bfi2           <- data.frame(scores, age = tmp$age, male, ed.d)
rownames(bfi2) <- NULL

marPreds <- c("age",
              "male",
              "finish_hs",
              "some_college",
              "college_grad",
              "graduate_degree")
targets  <- list(mar = setdiff(colnames(bfi2), marPreds),
                 mcar = NA,
                 mnar = NA)
pm       <- list(mar = 0.25)
snr      <- list(mar = 5)

dat1      <- bfi2[sample(c(1 : nrow(bfi2)), 550), ]
trainData <- dat1[1 : 500, ]
testData  <- dat1[501: 550, ]

missData <- SURF::imposeMissData(data    = trainData,
                                 targets = targets,
                                 preds   = marPreds,
                                 pm      = pm,
                                 snr     = snr)$data 

predictData <- list(train = trainData, test = testData, incomplete = missData)

save(predictData, file = "../source/mibrr/data/predictData.RData")
