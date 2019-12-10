### Title:    Initialize the Default Control List
### Author:   Kyle M. Lang
### Created:  2019-12-10
### Modified: 2019-12-10

MIBRR_CONTROL <- list(
    checkConv         = TRUE,
    convThresh        = 1.1,
    usePcStarts       = FALSE,
    smoothingWindow   = 10,
    minPredCor        = 0.1,
    miceIters         = 10,
    miceRidge         = 0.001,
    miceMethod        = "pmm",
    preserveStructure = TRUE,
    optTraceLevel     = 0L,
    optCheckKkt       = TRUE,
    optMethod         = "L-BFGS-B",
    optBoundLambda    = TRUE,
    savePpSams        = FALSE,
    useBetaMeans      = FALSE,
    optMaxRestarts    = 5L,
    optRestartRatio   = 0.1,
    optStrict         = TRUE,
    centerType        = "median",
    dumpParamHistory  = FALSE,
    phHistoryLength   = 10L
)
