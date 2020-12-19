args <- base::commandArgs(trailingOnly=TRUE)
if ("mo" %in% args) {
  write("Doing multi-objective optimization.", stderr())
}
std
stop(42)

#library(renv)
#.libPaths(normalizePath("./renv/library/R-4.0/x86_64-w64-mingw32"))
#renv::load()

#f <- file("stdin")
#open(f)
#line <- readLines(con = f, n = 1)
#theta <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(trimws(line))))



# TO BE READ FROM THE JSON PASSED TO STDIN!
#theta <- c(0.3776, 0.6595, 0.8929, 0.5608, 0.6438, 0.2005, 0.5477, 0.2816, 0.396, 0.9446, 0.5169, 0.0926, 0.5961, 0.9379, 0.9794, 0.7635, 0.119, 0.7881, 0.0146, 0.9049, 0.8775, 0.4115, 0.3353)
theta <- c(0.08566332, 0.6240046, 0.8747402,
           0, 0.11, 0.11, 0.9, 0.005, # A
           0, 0.375, 0.375, 0.075, 0.01, # CP
           0, 0.525, 0.525, 1, 0.025, # FREQ
           0, 0.78, 0.725, 0.98, 0.775 # SCD
           )


##################################### Compute the MLMRC with FireDrill!
source(file = "./notebooks/common-funcs.R")
source(file = "./models/modelsR6.R")

referenceData <- readRDS("./optim/referenceData.rds")
referenceBoundaries <- readRDS("./optim/referenceBoundaries.rds")
queryData <- readRDS("./optim/queryData.rds")

################### Model:

mlm <- MultilevelModel$new(
  referenceData = referenceData,
  referenceBoundaries = referenceBoundaries)

mlm$setQueryData(series = "ALL", queryData = queryData)
mlm$sourceFiles <- c("./notebooks/common-funcs.R", "./models/modelsR6.R")

################### Stages:

stage1_special <- Stage1Rectifier$new(
  yLimitsRef = "self", yLimitsQuery = "self")

stage1_ordinary <- Stage1NoModel$new(
  approxRefFun = TRUE, approxQueryFun = TRUE)
stage1_ordinary$setYLimitsRef("self")$setYLimitsQuery("self")

stage2_special <- Stage2Rectifier$new(
  scoreWarpingFunction = c("mono_rel", "start_end", "f_warp_org-vs-f_warp_np_opt"),
  scoreSignals = c("f_ref-vs-f_query_np", "f_ref_warp-vs-f_query_warp"),
  numSamples = 2e3)

stage2_ordinary <- Stage2$new(numSamples = 2e3)
stage2_ordinary$addScoreMethod(
  name = "Score_LM", scoreMethod = function(f1, f2, numSamples) {
    stat_diff_2_functions_lm_score2(maxAngleBetween = 30, requireSign = FALSE, numSamples = numSamples)(f1, f2)^5 })

################## Sub-Models:

create_submodel_special <- function(varName, intervalName, weight) {
  sm <- SubModel$new(varName = varName, intervalName = intervalName, weight = weight)
  sm$setStage1(stage1_special$clone(deep = TRUE))
  sm$setStage2(stage2_special$clone(deep = TRUE))
}

create_submodel_ordinary <- function(varName, intervalName, weight) {
  sm <- SubModel$new(varName = varName, intervalName = intervalName, weight = weight)
  sm$setStage1(stage1_ordinary$clone(deep = TRUE))
  s2 <- stage2_ordinary$clone(deep = TRUE)
  # Doing it this way will come in handy should we ever use this
  # in another interval.
  s2$addScoreMethod(name = "Score_IntervalLength",
                    scoreMethod = stat_diff_custom_score(
                      callback = custom_interval_length_score, mlm = mlm, interval = intervalName, useA = FALSE))
  sm$setStage2(s2)
}

sm_A_Begin <- create_submodel_special("A", "Begin", .8)
sm_CP_Begin <- create_submodel_special("CP", "Begin", .8)
sm_FREQ_Begin <- create_submodel_special("FREQ", "Begin", .8)
sm_SCD_Begin <- create_submodel_special("SCD", "Begin", .8)

sm_A_LongStretch <- create_submodel_ordinary("A", "LongStretch", .6)
sm_CP_LongStretch <- create_submodel_ordinary("CP", "LongStretch", .6)
sm_FREQ_LongStretch <- create_submodel_ordinary("FREQ", "LongStretch", .6)
sm_SCD_LongStretch <- create_submodel_ordinary("SCD", "LongStretch", .6)

sm_A_FireDrill <- create_submodel_special("A", "FireDrill", 1)
sm_CP_FireDrill <- create_submodel_special("CP", "FireDrill", 1)
sm_FREQ_FireDrill <- create_submodel_special("FREQ", "FireDrill", 1)
sm_SCD_FireDrill <- create_submodel_special("SCD", "FireDrill", 1)

sm_A_Aftermath <- create_submodel_special("A", "Aftermath", 1)
sm_CP_Aftermath <- create_submodel_special("CP", "Aftermath", 1)
sm_FREQ_Aftermath <- create_submodel_special("FREQ", "Aftermath", 1)
sm_SCD_Aftermath <- create_submodel_special("SCD", "Aftermath", 1)

mlm$setAllSubModels(sm_A_Begin, sm_CP_Begin, sm_FREQ_Begin, sm_SCD_Begin)
mlm$setAllSubModels(sm_A_LongStretch, sm_CP_LongStretch, sm_FREQ_LongStretch, sm_SCD_LongStretch)
mlm$setAllSubModels(sm_A_FireDrill, sm_CP_FireDrill, sm_FREQ_FireDrill, sm_FREQ_FireDrill)
mlm$setAllSubModels(sm_A_Aftermath, sm_CP_Aftermath, sm_FREQ_Aftermath, sm_FREQ_Begin)

################# Constraints:

mlm$flushLinIneqConstraints()
mlm$addDefaultBoundaryConstraints()

# distance >= 0.025 and <= 0.85
mlm$addDefaultBoundaryDistanceConstraints(boundaryDistance = .025, op = "geq")
mlm$addDefaultBoundaryDistanceConstraints(boundaryDistance = .850, op = "leq")

mlm$constrainBoundaryInterval(boundaryIndexOrName = 1, value = .01, op = "geq")
mlm$constrainBoundaryInterval(boundaryIndexOrName = 3, value = .99, op = "leq")

validationResults <- readRDS("./results/validationResults.rds")
mlm$boundariesCalibrated[1, ] <- validationResults[[1]]$fitResult$optResult$par

################## MLMRC!

mlmrc <- MultilevelModelReferenceCalibrator$new(mlm = mlm)
#mlmrc$addDefaultVariableYConstraints() # THOSE ARE 40(!) EXTRA CONSTRAINTS, LET'S SKIP
mlmrc$copyBoundaryConstraints()

###############################################
###############################################
# NOTE THAT IN PAGMO, THE LINEAR INEQUALITIES ARE SATISFIED IF THEY
# ARE NEGATIVE, I.E., inequalities are in the form g(x)<=0.
###############################################
###############################################
res <- c()
mlmrc$setTheta(theta = theta)
mlmrc$theta[mlmrc$theta == 0] <- .Machine$double.eps
constr <- -1 * as.vector((mlmrc$getUi() %*% mlmrc$getTheta()) - mlmrc$getCi())
if (!mlmrc$validateLinIneqConstraints()) {
  # DO NOT PROCEED, but return the violated constraints!
  #res <- c(.Machine$double.xmax, constr)
  temp <- constr[constr >= 0]
  score <- 20 + (prod(1 + temp) - 1) * 1e2
  res <- c(score, constr)
} else {
  res <- c(
    -log(mlmrc$compute(forceSeq = TRUE)$aggregateUsing_Honel()),
    constr)
}


write(jsonlite::base64_enc(jsonlite::toJSON(res)), stdout())
