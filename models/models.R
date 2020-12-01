library(foreach)


#' @param patternData data.frame with columns 'x', 'y' and the
#' two factor-columns 't' (type of variable) and 'interval'.
#' @return list where each entry is a chunk of the patternData
#' that corresponds to one type and one interval. The names of
#' the entries in this list are in the format "t_interval".
create_reference_data <- function(patternData) {
  # First, we divide the original pattern according to the
  # original boundaries, and store the reference signals.
  # We use the same names in the list as we expect for the
  # list of sub-models!
  referenceSignalData <- list()
  for (t in levels(patternData$t)) {
    for (i in levels(patternData$interval)) {
      varData <- patternData[
        patternData$t == t & patternData$interval == i, ]
      
      referenceSignalData[[paste(t, i, sep = "_")]] <- varData
    }
  }
  
  referenceSignalData
}



#' TODO: Description
#' 
#' @param fireDrillProject the readily processed project that is
#' suspected to contain a Fire Drill. Should have the same format
#' as the data.frame 'fireDrillPattern'.
#' @param listOfSubModels a list where each entry is a sub-model
#' for a single variable in a single interval. Each sub-model is
#' expected to return a score within [0,1] (where 1 is best) and
#' is given the reference- and query-signals. The name of each
#' sub-model in this list must follow this pattern:
#' "VARIABLENAME_INTERVALNAME". Also, each submodel may carry a
#' list of properties to be found in its attributes. Currently,
#' only the attribute 'weight' is used.
create_fire_drill_model <- function(
  fireDrillProject, listOfSubModels)
{ 
  #' @param x is a vector with the current boundaries.
  objectiveFunc <- function(x, returnAllScores = FALSE) {
    
    print(42)
    
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
      
      # The next step is to extract data from the project,
      # according to the current boundaries and interval.
      dataQuery <- fireDrillProject[
        fireDrillProject$t == vName &
        fireDrillProject$x >= boundaries[boundaryStart] &
        fireDrillProject$x < boundaries[boundaryEnd], ]
      
      subModel <- listOfSubModels[[varAndInterval]]
      subModelAttr <- attributes(subModel)
      subModelWeight <- if ("weight" %in% names(subModelAttr)) subModelAttr$weight else 1
      
      # Everything is prepared, let's call the model!
      vec <- c()
      vec[varAndInterval] <-
        subModelWeight * penalizeScore(subModel(dataQuery))
      vec
    }
    
    scores <- 1 + scores
    
    if (returnAllScores) {
      # These are weighted and penalized!
      scores
    } else {
      prod(scores)
    }
  }
  
  return(objectiveFunc)
}



# sub_model <- function(dataRef, dataQuery, stage1, stage2) {
#   
#   pair_of_funcs_or_vectors <- stage1(dataRef, dataQuery)
#   
#   final_score <- stage2(pair_of_funcs_or_vectors)
#   
#   return(final_score)
# 
# }



#' Generic sub-model that takes the reference data and two
#' previously created stages, as well as a weight for this
#' sub-model. Returns the product of all scores as returned
#' by stage 2.
generic_sub_model <- function(weight = 1, dataRef, stage1, stage2) {
  temp <- function(dataQuery) {
    scores <- stage2(stage1(dataQuery = dataQuery))
    # for the single scores of this sub-model prod is OK!
    prod(scores)
  }
  
  attributes(temp) <- list(weight = weight)
  temp
}


SubModel <- setClass(
  Class = "SubModel",
  
  slots = c(
    weight = "numeric",
    dataRef = "data.frame",
    dataQuery = "data.frame",
    
    stage1 = "function",
    stage1Result = "list",
    stage2 = "function",
    stage2Result = "list"
  ),
  
  prototype = list(
    weight = 1
  )
)


setValidity("SubModel", function(object) {
  if (object@weight < 0 || object@weight > 1) {
    return("Weight must be 0 <= w <= 1.")
  }
  TRUE
})

setMethod("initialize", "SubModel", function(.Object, ...) {
  .Object <- callNextMethod() # call super.initialize()
  
  if (!all(c("stage1", "stage2") %in% names(list(...)))) {
    stop("stage1 and stage2 are required to be functions.")
  }
  
  .Object
})



setGeneric("updateStage1", def = function(.Object, df) {
  standardGeneric("updateStage1")
})

setMethod("updateStage1", "SubModel", function(.Object, df) {
 .Object@dataRef <- df
 .Object@stage1Result <- .Object@stage1(df)
 .Object
})




setGeneric("updateStage2", def = function(.Object) {
  standardGeneric("updateStage2")
})

setMethod("updateStage2", "SubModel", function(.Object) {
  .Object <- .Object@stage2(.Object@stage1Result)
  .Object
})



setGeneric("fit", def = function(.Object) {
  standardGeneric("fit")
})

setMethod("fit", "SubModel", function(.Object) {
  sm <- updateStage1(.Object, .Object@dataQuery)
  sm <- updateStage2(sm)
  
  scores <- sm@stage2Result
  
  scores <- l$stage2Result
  # for the single scores of this sub-model prod is OK!
  prod(scores)
})

setMethod("plot", "SubModel", function(.Object, x, y, foo) {
  print(43)
  print(foo)
})


create_stage1_No_Model <- function(
  dataRef,
  yLimitsRef = range(dataRef$y),
  approxRefFun = TRUE, approxQueryFun = TRUE
) {
  fnRef <- if (approxRefFun) pattern_approxfun(
    yData = dataRef$y,
    xData = dataRef$x,
    yLimits = yLimitsRef) else NULL
  
  return(function(dataQuery, yLimitsQuery = yLimitsRef) {
    
    fnQuery <- if (approxQueryFun) pattern_approxfun(
      yData = dataQuery$y,
      xData = dataQuery$x,
      yLimits = yLimitsQuery) else NULL
    
    list(
      dataRef = dataRef,
      dataQuery = dataQuery,
      fnRef = fnRef,
      fnQuery = fnQuery
    )
  })
}


#' Creates a second stage for a sub-model that operates on
#' two functions and use arbitrary many score functions.
#' 
#' @param erroneousScore Function to be called on erroneous
#' scores (less than 0, larger than 1, NaN, NA etc.)
#' @param aggregateSubScores Function to use to aggregate
#' multiple scores as returned by some score-methods.
#' @param ... Any number of score-methods that accept two
#' functions f1, f2 and return a plain score (or a plain
#' vector of sub-scores).
create_stage2_two_functions_using_scores <- function(
  ...,
  erroneousScore = function(type, score) stop(
    paste0("The score ", type, " produced the erroneous value ", score, ".")),
  aggregateSubScores = prod
) {
  scoreMethods <- unname(list(...))
  
  return(function(stage1Result) {
    
    # This stage requires a previous stage that produced also
    # 2 functions for the data.
    
    # These two always exist,
    d1 <- stage1Result$dataRef
    d2 <- stage1Result$dataQuery
    # .. and these two may:
    f1 <- stage1Result$fnRef
    f2 <- stage1Result$fnQuery
    
    if (!is.function(f1) || !is.function(f2)) {
      stop("Stage 1 did not produce all functions.")
    }
    
    checkScore <- function(type, score) {
      if (!all(is.numeric(score)) || any(is.na(score)) || any(score < 0) || any(score > 1)) {
        erroneousScore(type, score)
      }
      score
    }
    
    
    # Each produced score will be added to this vector. After
    # all metrics were computed, the vector is returned to the
    # sub-model, which then decides how to proceed (e.g., by
    # using some aggregation).
    singleScores <- c()
    
    
    for (i in 1:length(scoreMethods)) {
      scoreMethod <- scoreMethods[[i]]
      score <- checkScore(deparse(scoreMethod), scoreMethod(f1 = f1, f2 = f2))
      singleScores <- c(singleScores, aggregateSubScores(score))
    }
    
    singleScores
  })
}



#' create_sub_model_stage1 <- function(dataRef, dataQuery) {
#'   
#'   # dataRef is never changing so we can go ahead and
#'   # approximate its function directly.
#'   fnRef <- pattern_approxfun(
#'     yData = dataRef$y,
#'     xData = dataRef$x,
#'     yLimits = range(dataRef$y))
#'   
#'   
#'   
#' }
#' 
#' 
#' 
#' #' The input to this model is the reference- and query signal
#' #' and it will output the score.
#' sub_model_no_model <- function(dataRef, dataQuery) {
#'   
#'   fnRef <- pattern_approxfun(
#'     yData = dataRef$y,
#'     xData = dataRef$x,
#'     yLimits = range(dataRef$y))
#'   
#'   # Let's make a function out of it:
#'   fnQuery <- pattern_approxfun(
#'     yData = dataQuery$y,
#'     xData = dataQuery$x,
#'     yLimits = range(dataRef$y)) # should be within reference
#'   
#'   # JSD as upper bound log(2), so we normalize it.
#'   # 1 minus that gives us a score.
#'   # That to a large power increase sensitivity of the score,
#'   # since JSD tends to give us larger scores.
#'   #return((1 - stat_diff_2_functions_symmetric_JSD(f1 = fnRef, f2 = fnQuery)$value / log(2))^4)
#'   
#'   mi <- stat_diff_2_functions_mutual_information(f1 = fnRef, f2 = fnQuery)
#'   return(
#'     1 *
#'     #max(0, stat_diff_2_functions_cor(f1 = fnRef, f2 = fnQuery)$value) * # negative correlation is bad!
#'     (1 - area_diff_2_functions(f1 = fnRef, f2 = fnQuery)$value) *
#'     (1 - stat_diff_2_functions_frechet(f1 = fnRef, f2 = fnQuery, numSamples = 60)$value) *
#'     ## MI should not be used when one of the signals is entirely flat,
#'     ## as due to rounding errors the value may get out of bounds.
#'     ## If we later develop a score and encounter this case, we need to warn or stop.
#'     #((mi$entropy1 / mi$value) * (mi$entropy2 / mi$value)) * # 'symmetric MI'
#'     #(mi$entropy2 / mi$value) * # 'asymmetric MI' - how much entropy does the query explain in the MI
#'     1
#'   )
#'   
#'   # For ^5, the slope of the function becomes > 1 at ~0.6687,
#'   # for ^6, this happens at ~0.6988,
#'   # for ^4, at ~0.6299,
#'   # for ^3, at ~0.5774
#'   
#' }

