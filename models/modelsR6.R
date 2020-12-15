library(R6)
library(foreach)
library(ggplot2)


MultilevelModelReferenceCalibrator <- R6Class(
  "MultilevelModelReferenceCalibrator",
  
  lock_objects = FALSE,
  
  inherit = LinearInequalityConstraints,
  
  private = list(
    scoreAggCallback = NULL
  ),
  
  public = list(
    initialize = function(mlm, scoreAggCallback = function(sa) sa$aggregateUsing_Honel()) {
      stopifnot(inherits(mlm, "MultilevelModel") && R6::is.R6(mlm))
      stopifnot(is.matrix(mlm$boundariesCalibrated) && !any(is.na(mlm$boundariesCalibrated)))
      stopifnot(is.function(scoreAggCallback) && length(methods::formalArgs(scoreAggCallback)) == 1)
      
      self$mlm <- mlm
      # Let's copy it over. During optimization, we will do our affine
      # transformations and then set the data back to the MLM, before
      # we call its compute()-method.
      self$refData <- mlm$refData[, ]
      
      # Needed later in fit()
      private$scoreAggCallback <- scoreAggCallback
      
      theta <- mlm$getTheta()
      bounds <- sort(unique(c(0, 1, mlm$refBoundaries)))
      # For each variable, we add the current intercept with each boundary
      # from left to right, including the 0,1 boundaries.
      # Each variable has 'numIntervals' + 1 y-values we can adjust.
      for (t in levels(mlm$refData$t)) {
        tFunc <- (function() {
          d <- mlm$refData[mlm$refData$t == t, ]
          f <- stats::approxfun(x = d$x, y = d$y)
          return(function(x) {
            if (x < min(d$x)) d$y[1]
            if (x > max(d$x)) utils::tail(d$y, 1)
            f(x)
          })
        })()
        for (bIdx in 1:length(bounds)) {
          newC <- c()
          newC[paste(t, bIdx, sep = "_")] <- tFunc(bounds[bIdx])
          theta <- c(theta, newC)
        }
      }
      
      super$initialize(theta = theta)
      self$refTheta <- theta
      colnames(self$linIneqs) <- c(names(theta), "CI")
    },
    
    #' Transfers over the MLM's boundary constraints into our matrix.
    copyBoundaryConstraints = function() {
      for (rIdx in rownames(self$mlm$linIneqs)) {
        newC <- rep(0, ncol(self$linIneqs))
        newC[1:self$mlm$numBoundaries] <- self$mlm$linIneqs[rIdx, 1:self$mlm$numBoundaries]
        newC[length(newC)] <- self$mlm$linIneqs[rIdx, self$mlm$numBoundaries + 1]
        mlmrc$setLinIneqConstraint(name = rIdx, ineqs = newC)
      }
      invisible(self)
    },
    
    addDefaultVariableYConstraints = function() {
      idxOffset <- length(self$mlm$getTheta())
      bounds <- sort(unique(c(0, 1, mlm$refBoundaries)))
      lvls <- levels(self$mlm$refData$t)
      
      for (t in lvls) {
        lvlIdx <- which(t == lvls)
        for (bIdx in 1:length(bounds)) {
          newC <- rep(0, ncol(self$linIneqs))
          newC[idxOffset + (lvlIdx - 1) * length(bounds) + bIdx] <- 1
          newC[length(newC)] <- 0
          # >= 0:
          self$setLinIneqConstraint(name = paste0(t, "_", bIdx, "_geq_v"), ineqs = newC)
          # .. and <= 1:
          newC <- -1 * newC
          newC[length(newC)] <- -1
          self$setLinIneqConstraint(name = paste0("-", t, "_", bIdx, "_geq_v"), ineqs = newC)
        }
      }
      invisible(self)
    },
    
    
    computeRefData = function(limitYExtent = c(0,1), bbFlatIfBelow = 1e-3) { 
      # Here, we need the calibrated boundaries -- we use them as the reference!
      # I.e., these are superior to the user-set reference boundaries!
      stopifnot(is.matrix(self$mlm$boundariesCalibrated) && !any(is.na(self$mlm$boundariesCalibrated)))
      stopifnot(missing(limitYExtent) ||
                (is.numeric(limitYExtent) && length(limitYExtent) == 2 && !any(is.na(limitYExtent))))
      # Current values for all boundaries and y-intersections (in that order).
      theta <- self$getTheta()
      # Needed to re-slice the reference data.
      refBounds <- sort(unique(c(0, 1, self$mlm$boundariesCalibrated[1, ])))
      # The current bounds that are used to transform the reference data.
      newBounds <- sort(unique(c(0, 1, theta[1:self$mlm$numBoundaries])))
      refData <- NULL
      
      
      for (t in levels(self$mlm$refData$t)) {
        for (bIdx in 1:self$mlm$numIntervals) {
          xyData <- self$refData[
            self$refData$t == t & self$refData$interval == self$mlm$intervalNames[bIdx], ]
          xyData <- xyData[order(xyData$x), ]
          
          xRange_ref <- range(xyData$x)
          xExtent_ref <- xRange_ref[2] - xRange_ref[1]
          bStart <- newBounds[bIdx]
          bEnd <- newBounds[bIdx + 1]
          xRange <- c(bStart, bEnd)
          xExtent <- bEnd - bStart
          xData <- xyData$x - xRange_ref[1]
          if (xExtent_ref > 0) {
            xData <- xData / xExtent_ref
          }
          if (xExtent > 0) {
            xData <- xData * xExtent
          }
          
          # For re-fitting Y, we do not need X.
          yExtent_ref <- range(xyData$y)
          yStart_ref <- xyData$y[1]
          yEnd_ref <- xyData$y[length(xyData$y)]
          yDiff_ref <- yEnd_ref - yStart_ref
          
          yStart_new <- theta[paste0(t, "_", bIdx)]
          yEnd_new <- theta[paste0(t, "_", bIdx + 1)]
          yDiff_new <- yEnd_new - yStart_new
          yData <- xyData$y
          
          # If the reference bounding box is flat, the current
          # sub-pattern should be regarded as a flat line. Otherwise,
          # we may amplify small deviations of an otherwise flat line
          # way beyond its proportions. So if we encounter a flat
          # line, we change its slope rather doing a linear affine
          # transformation with the sub-pattern.
          if (abs(max(yData) - min(yData)) < bbFlatIfBelow) {
            # Note that if yDiff_new is flat, we only squeeze
            # the pattern, and the else-branch will work.
            yDiffBefore <- yData[length(yData)] - yData[1]
            yDiffAfter <- yEnd_new - yStart_new
            yDiff <- yDiffBefore - yDiffAfter
            yData <- yData - sapply(xData, function(x) {
              x / xExtent * yDiff
            })
          } else {
            yScale <- yDiff_new / yDiff_ref
            yData <- yData * yScale
          }
          yData <- yData - yData[1] + yStart_new
          
          # Do this ALWAYS last!
          xData <- xData + xRange[1]
          
          
          # bStart_ref <- refBounds[bIdx]
          # bEnd_ref <- refBounds[bIdx + 1]
          # xRange_ref <- range(xyData$x)
          # xExtent_ref <- xRange_ref[2] - xRange_ref[1]
          # yLeft_ref <- xyData$y[1]
          # yRight_ref <- tail(xyData$y, 1)
          # yRange_ref <- range(xyData$y)
          # yExtent_ref <- yRange_ref[2] - yRange_ref[1]
          # 
          # 
          # bStart <- newBounds[bIdx]
          # bEnd <- newBounds[bIdx + 1]
          # xRange <- c(bStart, bEnd)
          # xExtent <- bEnd - bStart
          # yLeft <- theta[paste0(t, "_", bIdx)]
          # yRight <- theta[paste0(t, "_", bIdx + 1)]
          # yRange <- range(yLeft, yRight)
          # yExtent <- yRange[2] - yRange[1]
          # 
          # xData <- xyData$x - xRange_ref[1]
          # if (xExtent_ref > 0) {
          #   xData <- xData / xExtent_ref
          # }
          # if (xExtent > 0) {
          #   xData <- xData * xExtent
          # }
          # 
          # 
          # yData <- xyData$y - xyData$y[1]
          # yDiffBefore <- yData[length(yData)] - yData[1]
          # yDiffAfter <- yRight - yLeft
          # yDiff <- yDiffBefore - yDiffAfter
          # 
          # 
          # yData <- yData + yLeft
          # yData <- yData - sapply(xData, function(x) {
          #   x / xExtent * yDiff
          # })
          # xData <- xData + xRange[1]
          # 
          # # Check if we need to limit Y:
          # if (is.numeric(limitYExtent)) {
          #   yData[yData < limitYExtent[1]] <- limitYExtent[1]
          #   yData[yData > limitYExtent[2]] <- limitYExtent[2]
          # }
          
          refData <- rbind(refData, data.frame(
            x = xData,
            y = yData,
            t = t,
            interval = self$mlm$intervalNames[bIdx],
            stringsAsFactors = FALSE
          ))
        }
      }
      
      refData$interval <- factor(
        x = refData$interval,
        levels = levels(self$refData$interval),
        ordered = is.ordered(self$refData$interval))
      
      refData$t <- factor(
        x = refData$t,
        levels = levels(self$refData$t),
        ordered = is.ordered(self$refData$t))
      
      refData
    },
    
    
    #' Using the current theta, we derive the reference data, which
    #' is then given to the MLM and its compute() is called.
    compute = function(forceSeq = NULL) {
      stopifnot(is.matrix(self$mlm$boundariesCalibrated) && !any(is.na(self$mlm$boundariesCalibrated)))
      self$mlm$refBoundaries <- self$mlm$boundariesCalibrated[1, ]
      
      refData <- self$computeRefData()
      self$mlm$setAllBoundaries(values = self$theta[1:self$mlm$numBoundaries])
      self$mlm$refData <- refData
      
      # .. then compute using those!
      cArgs <- list()
      if (!missing(forceSeq)) {
        cArgs[["forceSeq"]] <- forceSeq
      }
      do.call(what = self$mlm$compute, args = cArgs)
    },
    
    #' Alters the underlying MLMs reference data and boundaries to
    #' produce a best fit between the model and all of its data.
    #' This method is extremely expensive and should probably never
    #' be used directly, but rather called from an external optimization
    #' process that sets theta and then calls compute in a highly
    #' parallel fashion. We keep this method as a reference, but it is
    #' mostly here for academic purposes.
    fit = function(verbose = FALSE, reltol = sqrt(.Machine$double.eps), method = c("Nelder-Mead", "SANN")[1], callback = NULL) {
      stopifnot(self$validateLinIneqConstraints())
      # callback must take MLMRC, MLM, fitHist, score_raw
      stopifnot(missing(callback) ||
                (is.function(callback) && 4 == length(methods::formalArgs(callback))))
      
      histCols <- c("begin", "end", "duration", "score_raw", "score_log", names(self$theta))
      fitHist <- matrix(nrow = 0, ncol = length(histCols))
      colnames(fitHist) <- histCols
      
      beginOpt <- as.numeric(Sys.time())
      optR <- stats::constrOptim(
        control = list(
          reltol = reltol,
          maxit = 2e4
        ),
        method = method,
        ui = self$getUi(),
        theta = self$getTheta(),
        ci = self$getCi(),
        grad = NULL,
        f = function(x) {
          names(x) <- names(self$theta)
          if (verbose) {
            cat(paste0("Params (", length(x), "): ", paste0(sapply(x, function(p) {
              format(x, nsmall = 3, digits = 3)
            }), collapse = ", ")))
          }
          
          # Set all parameters!
          self$setTheta(x)
          
          begin <- as.numeric(Sys.time())
          
          scoreAgg <- self$compute()
          score_raw <- private$scoreAggCallback(scoreAgg)
          score_log <- -log(score_raw)
          
          finish <- as.numeric(Sys.time())
          fitHist <<- rbind(fitHist, c(
            begin, finish, finish - begin, score_raw, score_log, x
          ))
          
          if (verbose) {
            cat(paste0(" -- Value: ", format(score_log, digits = 10, nsmall = 5),
                       " -- Duration: ", format(finish - begin, digits = 2, nsmall = 2), "s\n"))
          }
          
          if (is.function(callback)) {
            callback(self, self$mlm, fitHist, score_raw)
          }
          
          score_log
        }
      )
      finishOpt <- as.numeric(Sys.time())
      
      list(
        begin = beginOpt,
        end = finishOpt,
        duration = finishOpt - beginOpt,
        fitHist = fitHist,
        optResult = optR
      )
    }
  )
)




LinearInequalityConstraints <- R6Class(
  "LinearInequalityConstraints",
  
  lock_objects = FALSE,
  
  public = list(
    initialize = function(theta) {
      stopifnot(is.numeric(theta) && length(theta) > 0)
      
      self$numParams <- length(theta)
      self$linIneqs <- matrix(nrow = 0, ncol = self$numParams + 1) # +1 for theta-column
      self$theta <- theta
    },
    
    hasLinIneqConstraint = function(name) {
      stopifnot(is.character(name) && nchar(name) > 0)
      name %in% rownames(self$linIneqs)
    },
    
    removeLinIneqConstraint = function(name) {
      stopifnot(self$hasLinIneqConstraint(name))
      rIdx <- which(rownames(self$linIneqs) == name)
      self$linIneqs <- self$linIneqs[-rIdx, ]
      invisible(self)
    },
    
    setLinIneqConstraint = function(name, ineqs) {
      stopifnot(is.vector(ineqs) && length(ineqs) == self$numParams + 1)
      stopifnot(all(is.numeric(ineqs)) && !any(is.na(ineqs)))
      
      # Set or replace semantics:
      if (!self$hasLinIneqConstraint(name)) {
        newRow <- matrix(ncol = self$numParams + 1, nrow = 1)
        rownames(newRow) <- name
        self$linIneqs <- rbind(self$linIneqs, newRow)
      }
      
      self$linIneqs[name, ] <- ineqs
      invisible(self)
    },
    
    flushLinIneqConstraints = function() {
      self$linIneqs <- self$linIneqs[-1:-nrow(self$linIneqs), ]
      invisible(self)
    },
    
    setTheta = function(theta) {
      stopifnot(is.numeric(theta) && length(theta) == length(self$theta))
      names(theta) <- names(self$theta)
      self$theta <- theta
      invisible(self)
    },
    
    getTheta = function() {
      self$theta
    },
    
    getUi = function() {
      self$linIneqs[, 1:self$numParams]
    },
    
    getCi = function() {
      self$linIneqs[, self$numParams + 1]
    },
    
    validateLinIneqConstraints = function() {
      theta <- self$getTheta()
      ui <- self$getUi()
      ci <- self$getCi()
      
      res <- ui %*% theta - ci
      !any(is.na(res)) && all(res >= 0)
    }
  )
)



#' The model we use to describe, fit and score arbitrary many time
#' series is a Multilevel model.
#' 
#' @source {https://en.wikipedia.org/wiki/Multilevel_model}
MultilevelModel <- R6Class(
  "MultilevelModel",
  
  inherit = LinearInequalityConstraints,
  
  lock_objects = FALSE,
  
  private = list(
    subModels = NULL,
    scoreAggCallback = NULL
  ),
  
  public = list(
    #' @param intervalNames ordered character of intervals, i.e.,
    #' the first interval's name is the first interval. If there
    #' are m intervals, there must be m-1 boundaries.
    initialize = function(
      referenceData,
      intervalNames = levels(referenceData$interval),
      referenceBoundaries = NULL,
      boundaryNames = if (missing(referenceBoundaries)) NULL else names(referenceBoundaries),
      scoreAggCallback = function(sa) sa$aggregateUsing_Honel()
    ) {
      
      stopifnot(is.data.frame(referenceData) && nrow(referenceData) > 0)
      stopifnot(is.numeric(referenceData$x) && is.numeric(referenceData$y))
      stopifnot(is.factor(referenceData$t))
      stopifnot(is.factor(referenceData$interval) && length(levels(referenceData$interval)) > 1)
      stopifnot(all(is.character(intervalNames) & nchar(intervalNames) > 0))
      stopifnot(length(intervalNames) == length(levels(referenceData$interval)))
      stopifnot(all(intervalNames %in% levels(referenceData$interval)))
      stopifnot(length(intervalNames) > 1)
      stopifnot(missing(referenceBoundaries) || (all(
        is.numeric(referenceBoundaries) &
        referenceBoundaries >= 0 &
        referenceBoundaries <= 1)))
      stopifnot(missing(boundaryNames) || (length(intervalNames) == 1 + length(boundaryNames)))
      stopifnot(is.function(scoreAggCallback) && length(methods::formalArgs(scoreAggCallback)) == 1)
      
      # Order by x ascending
      self$refData <- referenceData[order(referenceData$x), ]
      self$numVars <- length(levels(referenceData$t))
      
      # For each series, we may have different data.
      self$queryData <- list()
      
      self$refBoundaries <- referenceBoundaries
      
      self$intervalNames <- intervalNames
      self$numIntervals <- length(intervalNames)
      
      # Now for n intervals, there will be n-1 boundaries.
      # We do not initialize them, however.
      self$numBoundaries <- length(levels(referenceData$interval)) - 1
      self$boundaries <- matrix(nrow = 1, ncol = self$numBoundaries)
      colnames(self$boundaries) <- boundaryNames # NULL is OK
      if (!missing(referenceBoundaries)) {
        self$boundaries[1, ] <- referenceBoundaries
      }
      
      # Now that we know the amount of boundaries, we can
      # initialize a structure for the linear inequalities:
      super$initialize(theta = if (missing(referenceBoundaries)) rep(NA, self$numBoundaries) else referenceBoundaries)
      
      
      # Next step is to generate slots for the sub-models.
      # For each variable in each interval, there may or
      # may not be a sub-model. The entire Multilevel Model
      # can only be fit if at least one sub-model is present.
      # The sub-models must be set using the naming scheme
      # "VARIABLE_INTERVAL", with the same names as in the
      # reference data.
      private$subModels <- list()
      for (t in levels(referenceData$t)) {
        for (i in levels(referenceData$interval)) {
          private$subModels[[paste(t, i, sep = "_")]] <- NA
        }
      }
      
      
      # We also allow additional sub-models that do not have
      # to capture a tuple of variable and interval. These
      # are called meta-models and are computed exactly as
      # the regular models.
      private$metaSubModels <- list()
      
      
      # Those properties will be filled if the model is calibrated
      self$boundariesCalibrated <- matrix(nrow = 1, ncol = self$numBoundaries)
      colnames(self$boundariesCalibrated) <- boundaryNames
      # To be of type ScoreAggregator
      self$scoreCalibrated <- NA
      # Also, ScoreAggregator; only set if the previous is set.
      self$scoreRef <- NA
      
      private$scoreAggCallback <- scoreAggCallback
      
      
      # Can be written to from outside. These files will be
      # sourced in the parallel foreach loop.
      self$sourceFiles <- c()
    },
    
    calibrate = function(requireAbsoluteConverge = FALSE, verbose = FALSE) {
      stopifnot(self$validateLinIneqConstraints())
      
      # Note this is a list where the names are the series' names.
      queryDataOld <- self$queryData
      boundariesOld <- self$boundaries[, ]
      scoreBefore <- private$scoreAggCallback(self$compute())
      
      self$setQueryData(series = "ALL", queryData = self$refData)
      fitResult <- self$fit(verbose = verbose)
      
      # Revert the Query data:
      for (series in names(queryDataOld)) {
        self$setQueryData(series = series, queryData = queryDataOld[[series]])
      }
      
      # The optimization might not converge. However, it
      # may still find better boundaries than the reference
      # boundaries. Better means a higher score. if the
      # parameter 'requireAbsoluteConvergence' is FALSE, even an
      # optimization that failed may provide such better
      # boundaries, which we will then use.
      if (requireAbsoluteConverge && resDummy$optResult$convergence != 0) {
        stop("The model requires convergence but failed to do so.")
      }
      if (verbose) {
        cat(paste0("The model did ", (if (resDummy$optResult$convergence == 0) "" else "not"), " converge.\n"))
      }
      
      # OK, let's check if the found boundaries provide
      # better scores.
      
      best <- fitResult$fitHist[
        which.max(fitResult$fitHist[, "score_raw"]), ]
      bestBounds <- tail(best, self$numBoundaries)
      self$setAllBoundaries(bestBounds)
      scoreAfterAgg <- self$compute()
      scoreAfter <- private$scoreAggCallback(scoreAfterAgg)
      if (scoreAfter > scoreBefore) {
        if (verbose) {
          cat(paste0("More optimal than the reference boundaries were found: ",
                     paste0(self$boundaries[1, ], collapse = ", ")))
          cat(paste0(" -- Old Score: ", scoreBefore, " -- New Score: ", scoreAfter, "\n"))
        }
        self$boundariesCalibrated[1, ] <- self$boundaries[1, ]
        self$scoreCalibrated <- scoreAfterAgg
        
        boundsBefore <- self$boundaries[1, ]
        self$setAllBoundaries(self$refBoundaries)
        self$scoreRef <- self$compute()
        self$setAllBoundaries(boundsBefore)
      } else {
        if (verbose) {
          cat("The calibration failed to find more optimal boundaries\n.")
        }
      }
      
      invisible(fitResult)
    },
    
    ###### Some overrides:
    setTheta = function(theta) {
      mlm$setAllBoundaries(theta)
    },
    
    getTheta = function() {
      self$boundaries[1, ]
    },
    
    getUi = function() {
      self$linIneqs[, 1:self$numBoundaries]
    },
    
    getCi = function() {
      self$linIneqs[, self$numBoundaries + 1]
    },
    
    constrainBoundaryInterval = function(boundaryIndexOrName, value, op = c("leq", "geq")[1]) {
      stopifnot(self$hasBoundary(boundaryIndexOrName))
      stopifnot(is.numeric(value) && !is.nan(value))
      stopifnot(op %in% c("leq", "geq"))
      isNeg <- if (op == "leq") -1 else 1
      
      boundaryIdx <- if(is.numeric(boundaryIndexOrName)) boundaryIndexOrName else which(colnames(self$boundaries) == boundaryIndexOrName)
      newConstr <- matrix(nrow = 1, ncol = self$numBoundaries + 1, data = c(rep(0, self$numBoundaries), isNeg * value))
      newConstr[1, boundaryIdx] <- isNeg * 1
      
      b1 <- colnames(self$boundaries)[boundaryIdx]
      if (op == "leq") {
        b1 <- paste0("-", b1)
      }
      
      self$setLinIneqConstraint(
        name = paste0(b1, "_geq_v"),
        ineqs = c(newConstr[1, ]))
      invisible(self)
    },
    
    constrainBoundaryDistance = function(boundary1_IndexOrName, boundary2_IndexOrName, value, op = c("leq", "geq")[1]) {
      stopifnot(self$hasBoundary(boundary1_IndexOrName))
      stopifnot(self$hasBoundary(boundary2_IndexOrName))
      stopifnot(is.numeric(value) && !is.nan(value) && value >= 0)
      stopifnot(op %in% c("leq", "geq"))
      
      isNeg <- if (op == "leq") -1 else 1
      
      boundaryIdx_1 <- if(is.numeric(boundary1_IndexOrName)) boundary1_IndexOrName else which(colnames(self$boundaries) == boundary1_IndexOrName)
      boundaryIdx_2 <- if(is.numeric(boundary2_IndexOrName)) boundary2_IndexOrName else which(colnames(self$boundaries) == boundary2_IndexOrName)
      stopifnot(boundaryIdx_1 < boundaryIdx_2)
      
      newConstr <- matrix(nrow = 1, ncol = self$numBoundaries + 1, data = c(rep(0, self$numBoundaries), isNeg * value))
      newConstr[1, boundaryIdx_1] <- isNeg * -1
      newConstr[1, boundaryIdx_2] <- isNeg * 1
      
      b1 <- colnames(self$boundaries)[boundaryIdx_1]
      b2 <- colnames(self$boundaries)[boundaryIdx_2]
      cName <- if (op == "leq") paste0("+", b1, "-", b2) else paste0("-", b1, "+", b2)
      
      
      self$setLinIneqConstraint(
        name = paste0(cName, "_geq_v"),
        ineqs = c(newConstr[1, ]))
      invisible(self)
    },
    
    fixateBoundaryAt = function(indexOrName, fixateAt, threshold = sqrt(.Machine$double.eps)) {
      self$setBoundary(indexOrName = indexOrName, value = fixateAt)
      self$constrainBoundaryInterval(boundaryIndexOrName = indexOrName, value = fixateAt + threshold, op = "leq")
      self$constrainBoundaryInterval(boundaryIndexOrName = indexOrName, value = fixateAt - threshold, op = "geq")
      
      invisible(self)
    },
    
    #' Sets constraints that ascertain that:
    #' - each boundary is within [0,1]
    addDefaultBoundaryConstraints = function() {
      
      for (bIdx in 1:self$numBoundaries) {
        self$constrainBoundaryInterval(
          boundaryIndexOrName = bIdx, value = 0, op = "geq")
        self$constrainBoundaryInterval(
          boundaryIndexOrName = bIdx, value = 1, op = "leq")
      }
      invisible(self)
    },
    
    #' Sets constraints that ascertain that:
    #' - two neighbouring boundaries have a minimum or
    #'   maximum distance (this function can be called
    #'   once using 'geq' and once using 'leq')
    addDefaultBoundaryDistanceConstraints = function(
      boundaryDistance = .Machine$double.eps, op = c("geq", "leq")[1])
    {
      stopifnot(op %in% c("leq", "geq"))
      
      for (bIdx in 1:self$numBoundaries) {
        if (bIdx > 1) {
          self$constrainBoundaryDistance(
            boundary1_IndexOrName = bIdx - 1,
            boundary2_IndexOrName = bIdx,
            value = boundaryDistance, op = op)
        }
      }
      invisible(self)
    },
    
    getCurrentIntervalRange = function(intervalIndexOrName) {
      intIdx <- if (is.numeric(intervalIndexOrName)) intervalIndexOrName else which(self$intervalNames == intervalIndexOrName)
      stopifnot(intIdx >= 1 && intIdx <= length(self$intervalNames))
      
      b1 <- if (intIdx == 1) c("[start]" = 0) else self$boundaries[1, intIdx - 1]
      b2 <- if (intIdx == length(self$intervalNames)) c("[end]" = 1) else self$boundaries[1, intIdx]
      
      c(b1, b2)
    },
    
    getMinMaxBoundaryDistances = function(boundary1_IndexOrName, boundary2_IndexOrName) {
      
      boundaryIdx_1 <- if(is.numeric(boundary1_IndexOrName)) boundary1_IndexOrName else which(colnames(self$boundaries) == boundary1_IndexOrName)
      boundaryIdx_2 <- if(is.numeric(boundary2_IndexOrName)) boundary2_IndexOrName else which(colnames(self$boundaries) == boundary2_IndexOrName)
      stopifnot(boundaryIdx_1 + 1 == boundaryIdx_2)
      
      b1 <- colnames(self$boundaries)[boundaryIdx_1]
      b2 <- colnames(self$boundaries)[boundaryIdx_2]
      cNameLeq <- paste0("+", b1, "-", b2, "_geq_v")
      cNameGeq <- paste0("-", b1, "+", b2, "_geq_v")
      
      c(
        "min" = if (cNameGeq %in% rownames(self$linIneqs)) self$linIneqs[cNameGeq, self$numBoundaries + 1] else NA,
        # this is negative because if leq
        "max" = if (cNameLeq %in% rownames(self$linIneqs)) abs(self$linIneqs[cNameLeq, self$numBoundaries + 1]) else NA
      )
    },
    
    
    
    #' (Un-)sets Query data for a series. Omit 'queryData' to unset.
    setQueryData = function(series, queryData = NA) {
      stopifnot(is.character(series) && nchar(series) > 0)
      stopifnot(missing(queryData) ||
        (is.data.frame(queryData) && is.numeric(queryData$x) && is.numeric(queryData$y) && is.factor(queryData$t)))
      
      if (missing(queryData) || is.na(queryData)) {
        self$queryData[[series]] <- NULL
      } else {
        self$queryData[[series]] <- queryData[order(queryData$x), ]
      }
      invisible(self)
    },
    
    
    
    #' Sets a sub-model. Also sets self as MLM of the sub-model.
    setSubModel = function(model) {
      stopifnot(inherits(model, "SubModel") && R6::is.R6(model))
      stopifnot(model$isMetaModel() ||
                (model$name %in% names(private$subModels))) # remember that slots are pre-allocated!
      
      if (model$isMetaModel()) {
        private$metaSubModels[[model$name]] <- model
      } else {
        private$subModels[[model$name]] <- model
      }
      
      model$setMLM(self)
      invisible(self)
    },
    
    setAllSubModels = function(...) {
      for (sm in list(...)) {
        self$setSubModel(sm)
      }
      invisible(self)
    },
    
    removeSubModel = function(model) {
      stopifnot(inherits(model, "SubModel") && R6::is.R6(model))
      
      isMeta <- model$isMetaModel()
      inSm <- model$name %in% names(private$subModels)
      inMsm <- model$name %in% names(private$metaSubModels)
      
      stopifnot((isMeta && inMsm) || (!isMeta && inSm))
      stopifnot(inSm || inMsm)
      
      model$setMLM(NA)
      if (inSm) {
        private$subModels[[model$name]] <- NA
      } else {
        private$metaSubModels[[model$name]] <- NULL # Note that we're un-setting!
      }
      invisible(self)
    },
    
    getSubModel = function(name) {
      inSm <- name %in% names(private$subModels)
      inMsm <- name %in% names(private$metaSubModels)
      
      stopifnot(inSm || inMsm)
      
      if (inSm) private$subModels[[name]] else private$metaSubModels[[name]]
    },
    
    getSubModelsInUse = function(
      includeOrdinarySubModels = TRUE, includeMetaSubModels = FALSE
    ) {
      modelNames <- c()
      if (includeOrdinarySubModels) {
        modelNames <- names(private$subModels)
      }
      if (includeMetaSubModels) {
        modelNames <- c(modelNames, names(private$metaSubModels))
      }
      
      inUse <- sapply(modelNames, function(name) {
        sm <- self$getSubModel(name)
        inherits(sm, "SubModel") && R6::is.R6(sm)
      })
      
      modelNames[inUse]
    },
    
    
    
    hasBoundary = function(indexOrName) {
      (is.character(indexOrName) && indexOrName %in% colnames(self$boundaries)) ||
      (is.numeric(indexOrName) && (indexOrName >= 1 || indexOrName <= ncol(self$boundaries)))
    },
    
    getBoundary = function(indexOrName) {
      stopifnot(self$hasBoundary(indexOrName = indexOrName))
      
      self$boundaries[1, indexOrName]
    },
    
    setBoundary = function(indexOrName, value) {
      stopifnot(self$hasBoundary(indexOrName = indexOrName))
      stopifnot(value >= 0, value <= 1)
      
      self$boundaries[1, indexOrName] <- value
      invisible(self)
    },
    
    setAllBoundaries = function(values) {
      stopifnot(is.vector(values) && length(values) == self$numBoundaries)
      stopifnot(is.null(names(values)) || all(colnames(self$boundaries) %in% names(values)))
      
      for (i in 1:length(values)) {
        self$setBoundary(indexOrName = i, value = values[i])
      }
      invisible(self)
    },
    
    
    fit = function(verbose = FALSE, reltol = sqrt(.Machine$double.eps), method = c("Nelder-Mead", "SANN")[1]) {
      stopifnot(self$validateLinIneqConstraints())
      
      histCols <- c("begin", "end", "duration", "score_raw", "score_log", colnames(self$boundaries))
      fitHist <- matrix(nrow = 0, ncol = length(histCols))
      colnames(fitHist) <- histCols
      
      
      beginOpt <- as.numeric(Sys.time())
      optR <- stats::constrOptim(
        control = list(
          reltol = reltol
        ),
        method = method,
        ui = self$getUi(),
        theta = self$getTheta(),
        ci = self$getCi(),
        grad = NULL,
        f = function(x) {
          if (verbose) {
            cat(paste0("Boundaries: ", paste0(sapply(x, function(b) {
              format(b, digits = 10, nsmall = 10)
            }), collapse = ", ")))
          }
          
          x <- x + (rnorm(length(x))/10)^5
          for (idx in 1:length(x)) {
            self$setBoundary(indexOrName = idx, value = x[idx])
          }
          
          begin <- as.numeric(Sys.time())
          
          scoreAgg <- self$compute()
          score_raw <- private$scoreAggCallback(scoreAgg)
          score_log <- -log(score_raw)
          
          finish <- as.numeric(Sys.time())
          fitHist <<- rbind(fitHist, c(
            begin, finish, finish - begin, score_raw, score_log, x
          ))
          
          if (verbose) {
            cat(paste0(" -- Value: ", format(score_log, digits = 10, nsmall = 5),
                       " --  Duration: ", format(finish - begin, digits = 2, nsmall = 2), "s\n"))
          }
          
          score_log
        }
      )
      finishOpt <- as.numeric(Sys.time())
      
      list(
        begin = beginOpt,
        end = finishOpt,
        duration = finishOpt - beginOpt,
        fitHist = fitHist,
        optResult = optR
      )
    },
    
    
    #' Evaluates the entire MLM using the currently set boundaries.
    compute = function(forceSeq = NULL) {
      stopifnot(!any(is.na(self$boundaries)))
      
      foreachOp <- if (missing(forceSeq) || forceSeq != TRUE) foreach::`%dopar%` else foreach::`%do%`
      smAggs <- foreachOp(foreach::foreach(
        subModelName = self$getSubModelsInUse(),
        .inorder = FALSE,
        .export = c("self", "private"),
        .packages = c("dtw", "Metrics", "numDeriv",
                      "philentropy", "pracma", "rootSolve",
                      "SimilarityMeasures", "stats", "utils")
      ), {
        for (file in self$sourceFiles) {
          source(file = file)
        }
        
        sm <- self$getSubModel(name = subModelName)
        
        # We may need to calculate the score w.r.t. more than one
        # query data series, and that's what the nested aggregator
        # is for. This is the case when fitting one MLM to multiple
        # data series at once (average model).
        saSm <- ScoreAggregator$new(namePrefix = subModelName)
        
        if (sm$isMetaModel()) {
          score <- private$scoreAggCallback(sm$compute())
          saSm$setRawScore(
            rawScore = RawScore$new(name = sm$name, value = score))
        } else {
          intervalIdx <- which(self$intervalNames == sm$intervalName)
          
          refData <- self$refData[
            self$refData$t == sm$varName & self$refData$interval == sm$intervalName, ]
          
          sm$setReferenceData(referenceData = refData)
          
          delimStart <- if (intervalIdx == 1) 0 else self$boundaries[1, intervalIdx - 1]
          delimEnd <- if (intervalIdx == length(self$intervalNames)) 1 else self$boundaries[1, intervalIdx]
          
          
          for (series in names(self$queryData)) {
            queryData <- self$queryData[[series]]
            queryData <- queryData[
              queryData$t == sm$varName &
              queryData$x >= delimStart & queryData$x < delimEnd, ]
            
            sm$setQueryData(queryData = queryData)
            # The 2nd stage of any sub-model can return any number
            # of scores. However, usually these are scores based
            # on low-level metrics and all have weight=1, so doing
            # the product or mean ought to be fine.
            tempSa <- sm$compute()
            score <- private$scoreAggCallback(tempSa)
            
            saSm$setRawScore(
              rawScore = RawScore$new(
                name = paste(sm$name, series, sep = "_"),
                value = score))
          }
        }
        
        saSm
      })
      #   smAggs <- append(smAggs, saSm)
      # }
      
      
      sa <- ScoreAggregator$new()
      lapply(smAggs, function(smAgg) {
        sm <- self$getSubModel(smAgg$prefix)
        
        sa$setRawScore(rawScore = RawScore$new(
          name = smAgg$prefix, value = private$scoreAggCallback(smAgg), weight = sm$weight))
      })
      
      sa
    },
    
    plot = function() {
      plot_fullPattern <- ggplot(data = self$refData, aes(x = x, y = y, color = t)) +
        geom_line() +
        labs(color = "Variable") +
        scale_x_continuous(
          breaks = seq(0, 1, by = .05)
        ) +
        theme_light() +
        theme(axis.text.x = element_text(angle = -45, vjust = 0)) +
        scale_color_brewer(palette = "Set1")
      
      for (bName in colnames(self$boundaries)) {
        if (!is.na(self$boundaries[1, bName])) {
          plot_fullPattern <- plot_fullPattern +
            geom_vline(xintercept = self$boundaries[1, bName], color = "blue", size = .4)
        }
      }
      
      if (!any(is.na(self$boundariesCalibrated))) {
        for (bName in colnames(self$boundariesCalibrated)) {
          plot_fullPattern <- plot_fullPattern +
            geom_vline(xintercept = self$boundariesCalibrated[1, bName], color = "red", size = .6)
        }
      }
      
      
      
      plot_fullPattern
    }
  )
)



RawScore <- R6Class(
  "RawScore",
  
  lock_objects = FALSE,
  
  public = list(
    initialize = function(name, value, weight = 1, allowNA = FALSE) {
      stopifnot(is.character(name) && nchar(name) > 0)
      stopifnot(allowNA || (is.numeric(value) && value >= 0 && value <= 1))
      stopifnot(is.numeric(weight) && weight >= 0) # let's allow weights > 1..
      
      self$name <- name
      self$value <- value
      self$weight <- weight
      self$allowNA <- allowNA
      self$isNA <- is.na(value)
    }
  )
)


ScoreAggregator <- R6Class(
  "ScoreAggregator",
  
  lock_objects = FALSE,
  
  public = list(
    #' @param allowNAScores The number of raw scores that can have
    #' an NA value. These are ignored when aggregating.
    initialize = function(namePrefix = NULL, allowNAScores = 0) {
      if (is.character(namePrefix) && nchar(namePrefix) > 0) {
        self$prefix <- namePrefix
      } else {
        self$prefix <- NULL
      }
      
      stopifnot(is.numeric(allowNAScores) && allowNAScores >= 0)
      self$allowNAScores = allowNAScores
      
      self$rawScores <- list()
      # The default is actually to use linear transformation.
      self$setNonLinearTransform(base::identity)
    },
    
    flushRawScores = function() {
      self$rawScores <- list()
    },
    
    setNonLinearTransform = function(callback) {
      stopifnot(is.function(callback))
      self$nonLinearTransform <- callback
      invisible(self)
    },
    
    getRawScores = function(includeNAScores = FALSE) {
      if (includeNAScores) {
        return(self$rawScores)
      } else {
        Filter(f = function(rs) !rs$isNA, x = self$rawScores)
      }
    },
    
    getNumNAScores = function() {
      length(Filter(f = function(rs) rs$isNA, x = self$rawScores))
    },
    
    getWeights = function(includeNAScores = FALSE) {
      unlist(lapply(self$getRawScores(includeNAScores = includeNAScores),
        function(rs) self$nonLinearTransform(rs$weight)))
    },
    
    getScoreValues = function(skipNonLinearTransform = FALSE, includeNAScores = FALSE) {
      useTrans <- if (skipNonLinearTransform) base::identity else self$nonLinearTransform
      unlist(lapply(self$getRawScores(includeNAScores = includeNAScores),
        function(rs) useTrans(rs$value)))
    },
    
    
    
    aggregateUsing_custom = function(callback) {
      stopifnot(is.function(callback))
      callback(self$getWeights() * self$getScoreValues())
    },
    
    aggregateUsing_mean = function() {
      self$aggregateUsing_custom(callback = mean)
    },
    
    aggregateUsing_prod = function() {
      stopifnot(all(self$getWeights() == 1))
      stopifnot(length(self$rawScores) > 0)
      self$aggregateUsing_custom(callback = prod)
    },
    
    aggregateUsing_Honel = function() {
      weights <- self$getWeights()
      scores <- self$getScoreValues()
      upperBound <- prod(1 + weights) - 1 # lower bound is 1
      
      (prod(1 + weights * scores) - 1) / upperBound
    },
    
    
    setRawScore = function(rawScore) {
      stopifnot(inherits(rawScore, "RawScore") && R6::is.R6(rawScore))
      
      if (rawScore$isNA && self$getNumNAScores() == self$allowNAScores) {
        stop("No more NA raw scores allowed.")
      }
      
      self$rawScores[[paste0(self$prefix, rawScore$name)]] <- rawScore
      invisible(self)
    },
    
    removeRawScore = function(nameOrRawScore) {
      stopifnot(is.character(nameOrRawScore) || R6::is.R6(nameOrRawScore))
      
      name <- if (R6::is.R6(nameOrRawScore)) nameOrRawScore$name else nameOrRawScore
      nameAlt <- paste0(self$prefix, name)
      if (name %in% names(self$rawScores)) {
        self$rawScores[[name]] <- NULL
      } else if (nameAlt %in% names(self$rawScores)) {
        self$rawScores[[nameAlt]] <- NULL
      }
      invisible(self)
    }
  )
)




MetaSubModel <- R6Class(
  "MetaSubModel",
  
  inherit = SubModel,
  
  lock_objects = FALSE,
  
  public = list(
    initialize = function(name, weight = 1) {
      super$initialize(varName = "_META_", intervalName = name, weight = weight)
    },
    
    isMetaModel = function() {
      TRUE
    }
  )
)



SubModel <- R6Class(
  "SubModel",
  
  lock_objects = FALSE,
  
  private = list(
    mlm = NULL
  ),
  
  public = list(
    initialize = function(varName, intervalName, referenceData = NULL, weight = 1) {
      stopifnot(all(is.character(c(varName, intervalName))))
      stopifnot(is.numeric(weight) || weight >= 0 || weight <= 1)
      
      self$varName <- varName
      self$intervalName <- intervalName
      self$name <- paste(varName, intervalName, sep = "_")
      if (missing(referenceData)) {
        self$setReferenceData()
      } else {
        self$setReferenceData(referenceData)
      }
      
      self$weight <- weight
      
      # A reference to the MLM, to be set by it
      private$mlm <- NA
    },
    
    isMetaModel = function() {
      FALSE
    },
    
    setMLM = function(mlm = NA) {
      stopifnot((inherits(mlm, "MultilevelModel") && R6::is.R6(mlm)) ||
                (missing(mlm) || is.na(mlm)))
      private$mlm <- mlm
      invisible(self)
    },
    
    getMLM = function() {
      private$mlm
    },
    
    setStage1 = function(stage1) {
      stopifnot(inherits(stage1, "Stage1") && R6::is.R6(stage1))
      self$stage1 <- stage1
      stage1$setSubModel(self)
      invisible(self)
    },
    
    setStage2 = function(stage2) {
      stopifnot(inherits(stage2, "Stage2") && R6::is.R6(stage2))
      self$stage2 <- stage2
      stage2$setSubModel(self)
      invisible(self)
    },
    
    
    compute = function() {
      stopifnot(inherits(self$stage1, "Stage1") && R6::is.R6(self$stage1))
      stopifnot(inherits(self$stage2, "Stage2") && R6::is.R6(self$stage2))
      
      self$stage1$setRefData(dataRef = self$getReferenceData())
      self$stage1$setQueryData(dataQuery = self$getQueryData())
      
      stage1Result <- self$stage1$compute()
      self$stage2$computeScores(stage1Result = stage1Result)
    },
    
    
    setReferenceData = function(referenceData = NULL) {
      if (missing(referenceData)) {
        self$refData <- NA
        return(invisible(self))
      }
      
      stopifnot(is.data.frame(referenceData))
      stopifnot(is.numeric(referenceData$x) && is.numeric(referenceData$y))
      stopifnot(is.factor(referenceData$t) && is.factor(referenceData$interval))
      # One SubModel can only represent one variable in one interval.
      stopifnot(length(unique(referenceData$t)) == 1 && length(unique(referenceData$interval)) == 1)
      
      self$refData <- referenceData
      invisible(self)
    },
    
    getReferenceData = function() {
      self$refData
    },
    
    
    
    setQueryData = function(queryData = NULL) {
      if (missing(queryData)) {
        self$queryData <- NA
        return(invisible(self))
      }
      
      stopifnot(is.data.frame(queryData))
      stopifnot(is.numeric(queryData$x) && is.numeric(queryData$y))
      stopifnot(is.factor(queryData$t))
      stopifnot(nrow(queryData[queryData$t == self$varName, ]) > 0)
      
      self$queryData <- queryData
      invisible(self)
    },
    
    getQueryData = function() {
      self$queryData
    }
  )
)



Stage2 <- R6Class(
  "Stage2",
  
  lock_objects = FALSE,
  
  private = list(
    subModel = NULL
  ),
  
  public = list(
    initialize = function(namePrefix = NULL) {
      # Init sub-score aggregation to default:
      #self$setSubScoreAggregation()
      
      self$namePrefix <- namePrefix
      self$scoreMethods <- list()
      self$scoreAgg <- ScoreAggregator$new(namePrefix = namePrefix)
      
      private$subModel <- NA
    },
    
    setSubModel = function(subModel) {
      stopifnot(inherits(subModel, "SubModel") && R6::is.R6(subModel))
      private$subModel <- subModel
      invisible(self)
    },
    
    
    addScoreMethod = function(scoreMethod, name) {
      stopifnot(is.function(scoreMethod) && is.character(name) && nchar(name) > 0)
      stopifnot(!(name %in% names(self$scoreMethods)))
      
      self$scoreMethods[[name]] <- scoreMethod
      invisible(self)
    },
    
    removeScoreMethod = function(name) {
      stopifnot(self$hasScoreMethod(name))
      self$scoreMethods[[name]] <- NULL
      invisible(self)
    },
    
    hasScoreMethod = function(name) {
      name %in% names(self$scoreMethods)
    },
    
    setScoreMethod = function(scoreMethod, name) {
      if (self$hasScoreMethod(name)) {
        self$removeScoreMethod(name)
      }
      
      self$addScoreMethod(scoreMethod = scoreMethod, name = name)
    },
    
    computeScores = function(stage1Result) {
      stopifnot(length(self$scoreMethods) > 0)
      
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
      
      self$scoreAgg$flushRawScores()
      
      for (smName in names(self$scoreMethods)) {
        sm <- self$scoreMethods[[smName]]
        score <- RawScore$new(
          name = smName, value = sm(f1 = f1, f2 = f2))
        self$scoreAgg$setRawScore(score)
      }
      
      self$scoreAgg
    }
  )
)



Stage1 <- R6Class(
  "Stage1",
  
  lock_objects = FALSE,
  
  private = list(
    subModel = NULL
  ),
  
  public = list(
    initialize = function(yLimitsRef = c("01", "self")[1], yLimitsQuery = c("ref", "01", "self")[1]) {
      self$dataRef <- NA
      self$dataQuery <- NA
      self$setYLimitsRef(yLimits = yLimitsRef)
      self$setYLimitsQuery(yLimits = yLimitsQuery)
      
      private$subModel <- NA
    },
    
    setSubModel = function(subModel) {
      stopifnot(inherits(subModel, "SubModel") && R6::is.R6(subModel))
      private$subModel <- subModel
      invisible(self)
    },
    
    setYLimitsRef = function(yLimits = c("01", "self")[1]) {
      stopifnot(yLimits %in% c("01", "self"))
      self$yLimitsRef <- yLimits
      invisible(self)
    },
    
    setYLimitsQuery = function(yLimits = c("ref", "01", "self")[1]) {
      stopifnot(yLimits %in% c("ref", "01", "self"))
      self$yLimitsQuery <- yLimits
      invisible(self)
    },
    
    setRefData = function(dataRef) {
      stopifnot(is.data.frame(dataRef) && all(c("x", "y") %in% colnames(dataRef)))
      self$dataRef <- dataRef
      invisible(self)
    },
    
    setQueryData = function(dataQuery) {
      stopifnot(is.data.frame(dataQuery) && all(c("x", "y") %in% colnames(dataQuery)))
      self$dataQuery <- dataQuery
      invisible(self)
    },
    
    
    compute = function() {
      stopifnot(is.data.frame(self$dataRef) && is.data.frame(self$dataQuery))
      NA
    }
  )
)


Stage1NoModel <- R6Class(
  "Stage1NoModel",
  
  inherit = Stage1,
  
  lock_objects = FALSE,
  
  public = list(
    initialize = function(approxRefFun = TRUE, approxQueryFun = TRUE) {
      super$initialize()
      
      stopifnot(is.logical(approxRefFun) && is.logical(approxQueryFun))
      
      self$approxRefFun <- approxRefFun
      self$approxQueryFun <- approxQueryFun
      
      self$fnRef <- NA
      self$fnQuery <- NA
    },
    
    #' Does only return a list with all data and the approximated
    #' functions. Calls Stage1::compute() to ensure these data were
    #' set previously. The returned list always contains the keys
    #' for the approximated functions, even if they were not created.
    #' 
    #' @return list with keys 'dataRef', 'dataQuery', 'fnRef', 'fnQuery'
    compute = function() {
      super$compute() # Only performs checks - does not return anything
      
      
      if (self$approxRefFun) {
        self$fnRef <- pattern_approxfun(
          yData = self$dataRef$y,
          xData = self$dataRef$x,
          yLimits = if (self$yLimitsRef == "01") {
            c(0, 1)
          } else { range(self$dataRef$y) })
      }
      
      if (self$approxQueryFun) {
        if (self$yLimitsQuery == "ref" && !is.data.frame(self$dataRef)) {
          stop("Need reference data for y-limits of query.")
        }
        
        self$fnQuery <- pattern_approxfun(
          yData = self$dataQuery$y,
          xData = self$dataQuery$x,
          yLimits = if (self$yLimitsQuery == "ref") {
            range(self$dataRef$y)
          } else if (self$yLimitsQuery == "01") {
            c(0, 1)
          } else {
            range(self$dataQuery$y)
          })
      }

            
      list(
        dataRef = self$dataRef,
        dataQuery = self$dataQuery,
        fnRef = self$fnRef,
        fnQuery = self$fnQuery
      )
    }
  )
)



Stage1Rectifier <- R6Class(
  "Stage1Rectifier",
  
  inherit = Stage1NoModel,
  
  lock_objects = FALSE,
  
  public = list(
    initialize = function(yLimitsRef = "01", yLimitsQuery = "ref") {
      super$initialize(approxRefFun = TRUE, approxQueryFun = TRUE)
      
      # Note that these only apply to when a function is approximated!
      # That means that the data is NOT changed.
      self$setYLimitsRef(yLimitsRef)
      self$setYLimitsQuery(yLimitsQuery)
    },
    
    compute = function() {
      # This is the list with data and functions.
      s1Result <- super$compute()
      
      # Before we begin, we may have to transform the data. We sample
      # from the approximated functions exactly at the original x-
      # positions.
      normalizeX <- function(xData) {
        xData <- xData - min(xData)
        if (max(xData) > 0) {
          xData <- xData / max(xData)
        }
        xData
      }
      
      dataRef <- data.frame(
        x = normalizeX(s1Result$dataRef$x),
        y = sapply(normalizeX(s1Result$dataRef$x), s1Result$fnRef))
      
      dataQuery <- data.frame(
        x = normalizeX(s1Result$dataQuery$x),
        y = sapply(normalizeX(s1Result$dataQuery$x), s1Result$fnQuery))
      
      
      align <- dtw::dtw(
        # Note that x=QUERY and y=REF!
        x = dataQuery$y,
        y = dataRef$y,
        keep.internals = TRUE,
        open.begin = FALSE,
        open.end = FALSE)
      
      ex <- extract_signal_from_window(
        dtwAlign = align,
        window = dataQuery$y,
        throwIfFlat = FALSE,
        idxMethod = "smooth",
        smoothCnt = 3)
      
      paw <- pattern_approxfun_warp(
        dtwAlign = align,
        includeOptimizedAlign = TRUE,
        signalRef = dataRef$y,
        signalQuery = dataQuery$y)
      
      # Should not be used manually, use get_dtw_scorable_functions
      #ex_warp <- extract_warping_from_dtw(
      #  dtwAlign = align,
      #  signalRef = dataRef$y,
      #  signalQuery = dataQuery$y)
      
      dtwFuncs <- get_dtw_scorable_functions(
        dtwAlign = align,
        signalRef = dataRef$y,
        signalQuery = dataQuery$y)
      
      s1Result$align <- align
      s1Result$ex <- ex
      s1Result$paw <- paw
      s1Result$dtwFuncs <- dtwFuncs
      
      s1Result
    }
  )
)



Stage2Rectifier <- R6Class(
  "Stage2Rectifier",
  
  inherit = Stage2,
  
  lock_objects = FALSE,
  
  public = list(
    initialize = function(
      scoreWarpingFunction = c(
        "mono", "mono_rel", "start_end", "resid", "rel_resid",
        "f_warp_org-vs-f_warp_np", "f_warp_org-vs-f_warp_np_opt",
        "f_warp-vs-f_lm"),
      scoreSignals = c(
        "f_ref-vs-f_query", "f_ref-vs-f_query_np",
        "f_ref-vs-f_ref_warp", "f_ref-vs-f_query_warp",
        "f_query-vs-f_query_warp", "f_query-vs-f_ref_warp",
        "f_ref_warp-vs-f_query_warp")
    ) {
      super$initialize()
      
      self$scoreWarping <- scoreWarpingFunction
      self$scoreSignals <- scoreSignals
    },
    
    computeScores = function(stage1Result) {
      stopifnot("align" %in% names(stage1Result))
      
      
      # stage1Result comes from Stage1Rectifier, and contains A LOT of
      # functions and values that we can possibly score. The results is
      # a merger of some operations and contains these:
      
      # ex = extract_signal_from_window(..)
      # (ex) monotonicity, monotonicity_rel, warp_resid, warp_rel_resid,
      # start/end rel&ref and more..
      # (ex) data - that is the query without plateaus. In the beginning
      # we used this a lot to find good alignments, and it is quite good.
      # It was the first method developed to assess the goodness of a dtw
      # alignment according to only the warping function. The alternative
      # to this function is to use what is in the current window, which
      # is everything in the current interval when begin/end are open.
      
      # paw = pattern_approxfun_warp(..)
      # (paw) contains pairs/triplets of functions that can be scored
      # against each other: 
      # - f_warp vs. f_lm (slightly worse than the next triplet)
      # - f_warp_org vs. (f_warp_np OR (often slightly better) f_warp_np_opt)
      # (paw) also contains the two original warping functions as known
      # the 3-way plot: f_warp_ref, f_warp_query. Those haven't really
      # been used until now, but I guess it would be straightforward to
      # score them (not necessarily together, however).
      
      # dtwFuncs = get_dtw_scorable_functions(..)
      # (dtwFuncs) contains 4 functions that can be matched pairwise.
      # They were obtained from extract_warping_from_dtw(..), which should
      # NOT be used directly.
      # - f_ref vs. (f_ref_warp OR f_query_warp) [1 A/B]
      # - f_query vs. (f_query_warp OR f_ref_warp) [2 A/B]
      # (dtwFuncs) f_ref and f_ref_warp represent the original reference
      # signal and how it would look after applying the warping from the
      # query. f_query and f_query_warp are the same, just for the query.
      # As for which function against which other function: The examples
      # from above seem reasonable for assessing a score (esp. 1B, 2B),
      # but comparing any of these functions is not well tested at the
      # moment, esp. not 1A, 2A. Also note that none of these functions
      # remove plateaus.
      
      #######################################################################
      
      
      # This is a short list of what we can do (not complete). We can also
      # cherry-pick and score only some of these:
      #
      # 1. Score warping function: monotonicity, residuals, start&end etc.,
      #    this is the old-school way that works sufficiently well. We should
      #    also include the start/end (how much of the query matched).
      # 2. Score warping functions against itself: This also proved to work
      #    well, esp. the warping functions against its np-version or the
      #    np&optimized version.
      # 3. Use either of the two pairs from dtwFuncs. This is currently not
      #    well tested.
      
      #######################################################################
      
      
      
      ex <- stage1Result$ex
      paw <- stage1Result$paw
      dtwFuncs <- stage1Result$dtwFuncs
      
      ################################### Direct scores:
      
      if ("mono" %in% self$scoreWarping) {
        self$scoreAgg$setRawScore(
          RawScore$new(name = "mono", value = ex$monotonicity))
      }
      if ("mono_rel" %in% self$scoreWarping) {
        self$scoreAgg$setRawScore(
          RawScore$new(name = "mono_rel", value = ex$monotonicity_rel))
      }
      if ("start_end" %in% self$scoreWarping) {
        # We score how much of the query could be mapped
        self$scoreAgg$setRawScore(
          RawScore$new(name = "start_end", value = ex$end_rel - ex$start_rel))
      }
      if ("resid" %in% self$scoreWarping) {
        self$scoreAgg$setRawScore(
          RawScore$new(name = "resid", value = ex$warp_resid$score))
      }
      if ("rel_resid" %in% self$scoreWarping) {
        self$scoreAgg$setRawScore(
          RawScore$new(name = "rel_resid", value = ex$warp_rel_resid$score))
      }
      
      ################################### Warping Function scores:
      
      aggregateSubScores <- function(name, methods, f1, f2) {
        tempSa <- ScoreAggregator$new(allowNAScores = 1)
        k <- 0
        for (m in methods) {
          tempSa$setRawScore(RawScore$new(
            name = paste0(k), value = m(f1, f2), allowNA = TRUE
          ))
          k <- k + 1
        }
        
        self$scoreAgg$setRawScore(
          RawScore$new(name = name, value = tempSa$aggregateUsing_prod()))
      }
      
      warpScoreMethods <- list(
        # Use upper bound from data is safer here.
        area_diff_2_functions_score(useUpperBoundFromData = TRUE),
        stat_diff_2_functions_cor_score(allowReturnNA = TRUE),
        stat_diff_2_functions_arclen_score(),
        stat_diff_2_functions_sd_var_mae_rmse_score(
          use = "rmse", useUpperBoundFromData = TRUE)
      )
      
      if ("f_warp_org-vs-f_warp_np" %in% self$scoreWarping) {
        # These are in [0,1]
        f1 <- paw$f_warp_org
        f2 <- paw$f_warp_np
        
        aggregateSubScores("f_warp_org-vs-f_warp_np", warpScoreMethods, f1, f2)
      }
      if ("f_warp_org-vs-f_warp_np_opt" %in% self$scoreWarping) {
        f1 <- paw$f_warp_org
        f2 <- paw$f_warp_np_opt
        
        aggregateSubScores("f_warp_org-vs-f_warp_np_opt", warpScoreMethods, f1, f2)
      }
      if ("f_warp-vs-f_lm" %in% self$scoreWarping) {
        f1 <- paw$f_warp
        f2 <- paw$f_lm
        
        aggregateSubScores("f_warp-vs-f_lm", warpScoreMethods, f1, f2)
      }
      
      ################################### Ref vs. Query Function scores:
      
      refVsQueryMethods <- list(
        stat_diff_2_functions_symmetric_JSD_score(),
        area_diff_2_functions_score(useUpperBoundFromData = TRUE),
        stat_diff_2_functions_cor_score())
      
      if ("f_ref-vs-f_query" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_ref
        f2 <- dtwFuncs$f_query
        
        aggregateSubScores("f_ref-vs-f_query", refVsQueryMethods, f1, f2)
      }
      if ("f_ref-vs-f_query_np" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_ref
        f2 <- pattern_approxfun(
          yData = ex$data,
          yLimits = range(stage1Result$dataRef$y))
        
        aggregateSubScores("f_ref-vs-f_query_np", refVsQueryMethods, f1, f2)
      }
      if ("f_ref-vs-f_ref_warp" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_ref
        f2 <- dtwFuncs$f_ref_warp
        
        aggregateSubScores("f_ref-vs-f_ref_warp", refVsQueryMethods, f1, f2)
      }
      if ("f_ref-vs-f_query_warp" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_ref
        f2 <- dtwFuncs$f_query_warp
        
        aggregateSubScores("f_ref-vs-f_query_warp", refVsQueryMethods, f1, f2)
      }
      if ("f_query-vs-f_query_warp" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_query
        f2 <- dtwFuncs$f_query_warp
        
        aggregateSubScores("f_query-vs-f_query_warp", refVsQueryMethods, f1, f2)
      }
      if ("f_query-vs-f_ref_warp" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_query
        f2 <- dtwFuncs$f_ref_warp
        
        aggregateSubScores("f_query-vs-f_ref_warp", refVsQueryMethods, f1, f2)
      }
      if ("f_ref_warp-vs-f_query_warp" %in% self$scoreSignals) {
        f1 <- dtwFuncs$f_ref_warp
        f2 <- dtwFuncs$f_query_warp
        
        aggregateSubScores("f_ref_warp-vs-f_query_warp", refVsQueryMethods, f1, f2)
      }
      
      self$scoreAgg
    }
  )
)




