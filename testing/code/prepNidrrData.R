### Title:    Prepare NIDRR Data
### Author:   Kyle M. Lang
### Created:  2016-MAY-05
### Modified: 2016-MAY-05

library(foreign)

dataDir2 <- "~/data/research/active/nomImps/appliedExample/data/"
fileName2 <- "NIDRR Follow Up Data.sav"

dat1 <- read.spss(paste0(dataDir2, fileName2),
                  to.data.frame = TRUE)

typeVec <- unlist(lapply(dat1, class))

dat2 <- dat1[ , typeVec == "numeric"]
colnames(dat2)

dropTerms <- c("total", "check", "date", "^id", "id$", "var")
dropVars <- colnames(dat2)[grep(paste0(dropTerms, collapse = "|"),
                                colnames(dat2),
                                ignore = TRUE)
                           ]

dat3 <- dat2[ , setdiff(colnames(dat2), dropVars)]

levelVec <- unlist(lapply(dat3, FUN = function(x) length(unique(na.omit(x)))))
constVars <- colnames(dat3)[levelVec == 0 | levelVec == 1]

dat4 <- dat3[ , setdiff(colnames(dat3), constVars)]

y1Names <- colnames(dat4)[-grep("Y2|Y3|Yr2|Yr3", colnames(dat4))]
y1Data <- dat4[ , y1Names]

sdsNames <- colnames(y1Data)[grep("sds", colnames(y1Data), ignore = TRUE)]
dat5 <- y1Data[ , setdiff(colnames(y1Data), sdsNames)]

airNames <- colnames(dat5)[grep("air", colnames(dat5), ignore = TRUE)]
airNames

dat6 <- dat5[ , setdiff(colnames(dat5), airNames)]

colnames(dat6)

tmp <- cor(dat6, use = "pairwise")

mean(tmp, na.rm = TRUE)
range(tmp, na.rm = TRUE)
quantile(tmp, c(0.05, 0.95), na.rm = TRUE)
