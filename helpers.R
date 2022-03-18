library(RMariaDB)


getExperimentConn = function(cnfFile = "../my.cnf") {
  return(dbConnect(
    RMariaDB::MariaDB(),
    default.file = cnfFile,
    group = "experiments")
  )
}

getDataset = function(dsName, removeUnwantedColums = TRUE) {
  conn <- getExperimentConn()
  result <- dbSendQuery(conn, paste("SELECT * FROM ", dsName))
  ds <- dbFetch(result)
  
  if (removeUnwantedColums) {
    removeNames <- c(
      #"SHA1",
      "RepoPathOrUrl",
      "AuthorName", "CommitterName", "AuthorTime",
      "CommitterTime", "MinutesSincePreviousCommit", "Message",
      "AuthorEmail", "CommitterEmail",
      "AuthorNominalLabel", "CommitterNominalLabel"
      #,"ParentCommitSHA1s"
      )
    
    ds <- ds[, !(names(ds) %in% removeNames)]
  }
  
  dbClearResult(result)
  dbDisconnect(conn)
  return(ds)
}

doWithParallelCluster <- function(expr, errorValue = NULL, numCores = parallel::detectCores()) {
  cl <- parallel::makePSOCKcluster(numCores)
  doSNOW::registerDoSNOW(cl)
  mev <- missing(errorValue)
  
  result <- tryCatch(expr, error=function(cond) {
    if (!mev) {
      return(errorValue)
    }
    return(cond)
  }, finally = {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    cl <- NULL
    gc()
  })
  return(result)
}

doWithParallelClusterExplicit <- function(cl, expr, errorValue = NULL, stopCl = TRUE) {
  doSNOW::registerDoSNOW(cl = cl)
  mev <- missing(errorValue)
  
  tryCatch(expr, error = function(cond) {
    if (!mev) {
      return(errorValue)
    }
    return(cond)
  }, finally = {
    if (stopCl) {
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
      gc()
    }
  })
}


#' Returns a list of seeds used in parallel training with caret. For
#' repeatability, we need deterministic seeds. The amount depends on
#' the amounts of hyperparamenters, and number of folds/repeats.
#' @param nh integer, the number of hyperparameters
#' @param amount integer, the number of seeds, usually this is number
#' of folds * number of repeats.
#' @param seed integer used in \code{set.seed()}. Given an identical
#' seed, this function produces the same seeds (idempotent).
#' @return list with seeds that can be used in caret's trainControl
get_seeds <- function(nh, amount, seed = 42) {
  set.seed(seed)
  
  seeds <- vector(mode = "list", length = amount + 1)
  for(i in 1:amount) seeds[[i]] <- sample.int(.Machine$integer.max, nh)
  # For the last model:
  seeds[[amount + 1]] <- sample.int(.Machine$integer.max, 1)
  return(seeds)
}

balanceDatasetSmote <- function(data, stateColumn) {
  lvls <- if (is.factor(data[[stateColumn]])) levels(data[[stateColumn]]) else NULL
  d <- table(data[[stateColumn]])
  m <- names(which.max(d))
  # We'll sample all other classes until we reach this for each:
  targetAmount <- d[[m]]
  
  # Get the other classes:
  otherClasses <- names(d)[!(names(d) %in% m)]
  
  # Add the over-represented class already to the final data:
  dataLargestClass <- data[data[[stateColumn]] == m, ]
  dataFinal <- dataLargestClass[, ]
  dataFinal[[stateColumn]] <- as.character(dataFinal[[stateColumn]])
  
  # Now, for each class, over-sample it and add to final frame:
  for (oc in otherClasses) {
    dataOtherClass <- data[data[[stateColumn]] == oc, ]
    temp <- rbind(dataLargestClass, dataOtherClass)
    
    # SMOTE requires factor-labels:
    temp[[stateColumn]] <- factor(temp[[stateColumn]])
    
    overSampled <- DMwR::SMOTE(
      form = formula(paste0(stateColumn, "~.")),
      data = temp,
      perc.over = 100 * ceiling(nrow(dataLargestClass) / nrow(dataOtherClass)),
      perc.under = 100
    )
    
    # Since we rounded up, let's only sample what we need:
    overSampled <- overSampled[overSampled[[stateColumn]] == oc, ]
    overSampled <- overSampled[sample(
      x = rownames(overSampled), size = min(nrow(overSampled), nrow(dataLargestClass))), ]
    
    # .. change to character again:
    overSampled[[stateColumn]] <- as.character(overSampled[[stateColumn]])
    dataFinal <- rbind(dataFinal, overSampled)
  }
  
  if (is.character(lvls)) {
    dataFinal[[stateColumn]] <- factor(dataFinal[[stateColumn]], levels = lvls)
  }
  
  return(dataFinal)
}


loadResultsOrCompute <- function(file, computeExpr) {
  file <- base::normalizePath(file, mustWork = FALSE)
  if (file.exists(file)) {
    return(base::readRDS(file))
  }
  
  res <- base::tryCatch(
    expr = computeExpr, error = function(cond) cond)
  
  # 'res' may have more than one class.
  if (any(class(res) %in% c("simpleError", "error", "condition"))) {
    print(traceback())
    stop(paste0("The computation failed: ", res))
  }
  
  base::saveRDS(res, file)
  return(res)
}

caretFitOneModeltoAllData <- function(method, tuneGrid, data) {
  set.seed(42)
  tr <- caret::trainControl(
    method = "none", p = 1, returnResamp = "all"
    , savePredictions = "all", classProbs = TRUE
    , number = 1)
  
  set.seed(43)
  caret::train(
    label ~ ., data = data, trControl = tr,
    tuneGrid = tuneGrid, preProcess = c("center", "scale"),
    method = method, verbose = FALSE)
}


saveAndPlotAsEPS <- function(ggplotInstance, fileName, width = 241.14749 / 72.27 * 2.54, height = 5 / 2.54) {
  ggplot2::ggsave(fileName, ggplotInstance,
         width = floor(width * 100) / 100,
         height = floor(height * 100) / 100,
         limitsize = F, device = cairo_pdf)
  ggplotInstance
}


curve2 <- function(func, from, to, col = "black", lty = 1, lwd = 1, add = FALSE, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL, ...) {
  f <- function(x) func(x)
  curve(expr = f, from = from, to = to, col = col, lty = lty, lwd = lwd, add = add, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main = main, ... = ...)
}



