library(foreach)

#' TODO: Description
#' 
#' @param fireDrillPattern the readily processed Fire Drill
#' pattern, as a single data.frame. Needs to contain data for all
#' variables and all intervals. The expected columns are:
#' 'x', 'y', 't' (variable) and 'interval' (name of interval).
#' @param fireDrillProject the readily processed project that is
#' suspected to contain a Fire Drill. Should have the same format
#' as the data.frame 'fireDrillPattern'.
#' @param listOfSubModels a list where each entry is a sub-model
#' for a single variable in a single interval. Each sub-model is
#' expected to return a score within [0,1] (where 1 is best) and
#' is given the reference- and query-signals. The name of each
#' sub-model in this list must follow this pattern:
#' "VARIABLENAME_INTERVALNAME".
#' @param subModelWeights a named vector with linear weights for
#' each sub-model. Same naming convention as in the list of sub-
#' models!
create_fire_drill_model <- function(
  fireDrillPattern, fireDrillProject, listOfSubModels,
  subModelWeights = sapply(names(listOfSubModels), function(n) {
    vec <- c()
    vec[n] <- 1
    vec
  }, USE.NAMES = FALSE))
{ 
  # First, we divide the original pattern according to the
  # original boundaries, and store the reference signals.
  # We use the same names in the list as we expect for the
  # list of sub-models!
  referenceSignalData <- list()
  referenceSignalFuncs <- list()
  for (t in levels(fireDrillPattern$t)) {
    for (i in levels(fireDrillPattern$interval)) {
      varData <- fireDrillPattern[fireDrillPattern$t == t & fireDrillPattern$interval == i, ]
      
      referenceSignalData[[paste(t, i, sep = "_")]] <- varData
      
      referenceSignalFuncs[[paste(t, i, sep = "_")]] <- pattern_approxfun(
        yData = varData$y,
        xData = varData$x,
        yLimits = range(varData$y))
    }
  }
  
  
  #' @param x is a vector with the current boundaries.
  objectiveFunc <- function(x, returnAllScores = FALSE) {
    
    scores <- foreach::foreach(
      varAndInterval = names(listOfSubModels),
      .combine = c,
      .inorder = FALSE,
      .packages = c("dtw", "Metrics", "numDeriv",
                    "philentropy", "pracma", "rootSolve",
                    "SimilarityMeasures", "stats", "utils")
    ) %dopar% {
      boundaries <- sort(unique(c(0, 1, x)))
      sp <- strsplit(varAndInterval, "_")[[1]]
      vName <- sp[1]
      iName <- sp[2]
      
      boundaryStart <- if (iName == "Begin") {
        1 } else if (iName == "LongStretch") {
        2 } else if (iName == "FireDrill") {
        3 } else { 4 }
      boundaryEnd <- boundaryStart + 1
      
      dataRef <- referenceSignalData[[varAndInterval]]
      fnRef <- referenceSignalFuncs[[varAndInterval]]
      
      # The next step is to extract data from the project,
      # according to the current boundaries and interval.
      dataQuery <- fireDrillProject[
        fireDrillProject$t == vName &
        fireDrillProject$x >= boundaries[boundaryStart] &
        fireDrillProject$x < boundaries[boundaryEnd], ]
      
      # Let's make a function out of it:
      fnQuery <- pattern_approxfun(
        yData = dataQuery$y,
        xData = dataQuery$x#,
        #yLimits = range(dataRef$y) # make it is within reference
      )
      
      subModel <- listOfSubModels[[varAndInterval]]
      subModelWeight <- subModelWeights[[varAndInterval]]
      
      # Everything is prepared, let's call the model!
      vec <- c()
      vec[varAndInterval] <-
        subModelWeight * subModel(dataRef, fnRef, dataQuery, fnQuery)
      vec
    }
    
    scores <- 1 + penalizeScore(scores)
    
    if (returnAllScores) {
      # These are weighted and penalized!
      scores
    } else {
      prod(scores)
    }
  }
  
  return(objectiveFunc)
}



#' The input to this model is the reference- and query signal
#' and it will output the score.
sub_model_no_model <- function(dataRef, fnRef, dataQuery, fnQuery) {
  
  # JSD as upper bound log(2), so we normalize it.
  # 1 minus that gives us a score.
  # That to a large power increase sensitivity of the score,
  # since JSD tends to give us larger scores.
  return((1 - stat_diff_2_functions_symmetric_JSD(f1 = fnRef, f2 = fnQuery)$value / log(2))^4)
  
  # For ^5, the slope of the function becomes > 1 at ~0.6687,
  # for ^6, this happens at ~0.6988,
  # for ^4, at ~0.6299,
  # for ^3, at ~0.5774
  
}

