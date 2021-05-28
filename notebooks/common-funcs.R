pattern_approxfun <- function(yData, smooth = FALSE, yLimits = c(0, 1), smoothSpan = .15, xData = NULL) {
  y1 <- yData[1]
  if (all(yData == y1)) {
    return(approxfun(
      x = c(0, 1),
      y = c(y1, y1)
    ))
  }
  
  if (!missing(xData)) {
    # There is xData, probably because yData was not uniformly sampled.
    # Or maybe yData is ordered according to xData..
    # It's fine, but we need to scale it to [0,1]
    xData <- xData - min(xData)
    if (max(xData) > 0) {
      xData <- xData / max(xData)
    }
  }
  
  if (smooth) {
    temp <- stats::loess.smooth(
      x = if (missing(xData)) 1:length(yData) else xData,
      y = yData,
      span = smoothSpan,
      evaluation = length(yData)
    )
    
    yData <- temp$y
    xData <- temp$x - min(temp$x)
    if (max(xData) > 0) {
      xData <- xData / max(xData)
    }
  } else {
    xData <- if (missing(xData)) seq(0, 1, by = 1 / (length(yData) - 1)) else xData
  }
  
  yData <- yData - min(yData)
  if (max(yData) > 0) {
    yData <- yData / max(yData)
  }
  
  if (!missing(yLimits)) {
    # Scale Y into the requested range:
    yData <- yData * (max(yLimits) - min(yLimits)) + min(yLimits)
  }
  
  f <- stats::approxfun(
    x = xData,
    y = yData
  )
  
  return(function(x) {
    if (x < 0 || x > 1) {
      warning("Function used out of bounds [0,1], returning 0.")
      return(0)
    }
    return(f(x))
  })
}


# Let's define a function to plot both patterns as functions:
plot_2_functions <- function(fReference, fQuery, interval = c(0, 1)) {
  ggplot2::ggplot(data.frame(x = interval), ggplot2::aes(x)) +
    ggplot2::stat_function(fun = fReference, ggplot2::aes(color="Reference")) +
    ggplot2::stat_function(fun = fQuery, ggplot2::aes(color="Query")) +
    ggplot2::scale_color_manual("Patterns", values = c("black", "red")) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank())
}


extract_signal_from_window <- function(
  dtwAlign, window, throwIfFlat = TRUE,
  idxMethod = c("discrete", "smooth"), smoothCnt = 3
) {
  idxMethod <- if (missing(idxMethod)) idxMethod[1] else idxMethod
  
  # First, we check whether the warping function is flat. This
  # test has two components: a linear regression, and a check
  # of the vertical co-domain of that function.
  tempLm <- stats::lm(
    formula = y ~ x, data = data.frame(
      x = 1:length(dtwAlign$index2),
      y = dtwAlign$index2
    )
  )
  # Note the slope of the warping function is never negative..
  if (stats::coef(tempLm)["x"] < .1) {
    # The slope is less than 0.1.
    # It is very flat, let's check the values in the co-domain:
    if (((max(dtwAlign$index2) - min(dtwAlign$index2)) / nrow(dtwAlign$reference)) < .1) {
      # Also, the values in the co-domain cover less than 10%
      # of the available range, indicating an actual flat line.
      # We will return two values, both 0, so that a function
      # can be approximated.
      wMsg <- "DTW warping function is flat."
      if (throwIfFlat) {
        stop(wMsg)
      }
      warning(wMsg)
      
      return(list(
        monotonicity = NA,
        monotonicity_rel = NA,
        warp_resid = list(
          lm = NULL,
          score = NA,
          sd = NA,
          var = NA,
          mae = NA,
          rmse = NA
        ),
        warp_rel_resid = list(
          lm = NULL,
          score = NA,
          sd = NA,
          var = NA,
          mae = NA,
          rmse = NA
        ),
        indices = c(0, 0),
        start = NA,
        start_ref = NA,
        start_rel = NA,
        start_rel_ref = NA,
        end = NA,
        end_ref = NA,
        end_rel = NA,
        end_rel_ref = NA,
        data = c(0, 0)
      ))
    }
  }
  
  # Also, we want to do a simple linear regression over the entire window.
  # This is often needed to, e.g., check the slope of the signal within it.
  windowLm <- stats::lm(
    formula = y ~ x, data = data.frame(
      x = 1:length(window),
      y = window
    )
  )
  windowLmCoef <- stats::coef(windowLm)
  
  if (idxMethod == "discrete") {
    # Here we store all indices that have a next index that has a value
    # greater than that of the current index. This means there is a slope
    # with a value > 0. In such a case, we store the index and its next
    # index in an array.
    indices <- c()
    for (idx in 1:(length(dtwAlign$index2) - 1)) {
      if (dtwAlign$index2[idx] < dtwAlign$index2[idx + 1]) {
        indices <- c(indices, c(idx, idx + 1))
      }
    }
    indices <- unique(indices)
  } else if (idxMethod == "smooth") {
    indices <- c()
    for (idx in 1:(length(dtwAlign$index2) - smoothCnt + 1)) {
      temp <- dtwAlign$index2[idx:(idx + smoothCnt - 1)]
      if (max(temp) > min(temp)) {
        indices <- c(indices, idx:(idx + smoothCnt - 1))
      }
    }
    indices <- unique(indices)
  }
  
  # Additional flatness-check:
  if (throwIfFlat && length(indices) == 0) {
    stop("Warping function is flat, no path between signals exist.")
  }
  
  # An additional metric is computed from the linear model of the warping
  # function, where we measure the residuals of it. We do this twice, once
  # for the entire function, once after cutting in and out.
  normalizeData <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }
  
  warpLm <- stats::lm(
    formula = y ~ x, data = data.frame(
      x = seq(0, 1, len = length(dtwAlign$index2)),
      y = normalizeData(dtwAlign$index2)
    ))
  warpLm_rel <- stats::lm(
    formula = y ~ x, data = data.frame(
      x = seq(0, 1, len = length(min(indices):max(indices))),
      y = normalizeData(
        dtwAlign$index2[min(indices):max(indices)])
    ))
  
  warp_res <- stats::residuals(warpLm)
  warp_rel_res <- stats::residuals(warpLm_rel)
  
  return(list(
    # We also are interested in the monotonicity or continuity of the
    # reference and the signal. Ideally, all indices are in the previous
    # list, indicating a strong monotonic behavior. For saddle points
    # however, indices will be missing. With monotonicity we are measuring
    # actually the continuity of the warping function.
    # Note also that the absolute monotonicity is equal to how much of the
    # query could be mapped to the reference.
    monotonicity = length(indices) / length(dtwAlign$index2),
    # Cut off flat start and end and measure again:
    monotonicity_rel = length(indices) / (max(indices) - min(indices) + 1),
    
    warp_resid = list(
      lm = warpLm,
      # times 2 because in the worst case, the difference between a
      # linear regression and its data in the unit square is 1/2.
      score = 1 - (mean(abs(warp_res)) * 2),
      sd = stats::sd(warp_res),
      var = stats::var(warp_res),
      mae = mean(abs(warp_res)),
      rmse = sqrt(mean(warp_res^2))
    ),
    warp_rel_resid = list(
      lm = warpLm_rel,
      score = 1 - (mean(abs(warp_rel_res)) * 2),
      sd = stats::sd(warp_rel_res),
      var = stats::var(warp_rel_res),
      mae = mean(abs(warp_rel_res)),
      rmse = sqrt(mean(warp_rel_res^2))
    ),
    
    window_info = list(
      lm = windowLm,
      mean = mean(window),
      slope = windowLmCoef["x"],
      intercept = windowLmCoef["(Intercept)"]
    ),
    
    indices = indices,
    start_ref = min(dtwAlign$index2[indices]),
    start_rel_ref = (min(dtwAlign$index2[indices]) - 1) / (dtwAlign$M - 1),
    end_ref = max(dtwAlign$index2[indices]),
    end_rel_ref = (max(dtwAlign$index2[indices]) - 1) / (dtwAlign$M - 1),
    
    start = min(dtwAlign$index1[indices]),
    start_rel = (min(dtwAlign$index1[indices]) - 1) / (dtwAlign$N - 1),
    end = max(dtwAlign$index1[indices]),
    end_rel = (max(dtwAlign$index1[indices]) - 1) / (dtwAlign$N - 1),
    
    # That data in the window that could be warped.
    data = window[dtwAlign$index1[indices]]
  ))
}


area_diff_2_functions <- function(f1, f2, useCubintegrate = TRUE) {
  # Find the intersections of both functions:
  intersections <- sort(rootSolve::uniroot.all(
    f = function(x) f1(x) - f2(x),
    interval = c(0, 1)))
  
  if (length(intersections) == 0) {
    # One function is complete below/above the other
    intersections <- c(0, 1)
  }
  
  # Check that lower/upper integration boundaries exist:
  if (intersections[1] > 0) {
    intersections <- c(0, intersections)
  }
  if (utils::tail(intersections, 1) < 1) {
    intersections <- c(intersections, 1)
  }
  
  integrateFn <- function(f, l, u) {
    if (useCubintegrate) {
      cubature::cubintegrate(
        f = f, lower = l, upper = u)$integral
    } else {
      stats::integrate(
        f = f, lower = l, upper = u, subdivisions = 1e5)$value
    }
  }
  
  
  # Now, for each pair of intersections, we integrate both
  # functions and sum up the areas.
  
  areas <- c()
  for (intsec in 1:(length(intersections) - 1)) {
    temp <- abs(
      integrateFn(
        f = f1, l = intersections[intsec], u = intersections[intsec + 1])
      -
      integrateFn(
        f = f2, l = intersections[intsec], u = intersections[intsec + 1])
    )
    
    areas <- c(areas, temp)
  }
  
  return(list(
    areas = areas,
    value = sum(areas),
    intersections = intersections
  ))
}


#' Computes the are between two function in the unit-square. That area
#' is directly anti-proportional to the resulting score, i.e., a low
#' area between corresponds to a high score. Returns a stat_diff-style
#' function with f1, f2.
#' 
#' Note: While this function will also work for when both function use
#' smaller or larger bounds than the unit-square by setting the para-
#' meter 'useUpperBoundFromData' to TRUE and using a sufficiently large
#' number of 'numSamples', it is generally recommeded to use the score
#' stat_diff_2_functions_sd_var_mae_rmse_score() with MAE instead in
#' these cases, also setting 'useUpperBoundFromData' to true. Otherwise,
#' this function should give more precise results as it uses true inte-
#' gration.
#' 
#' @param useUpperBoundFromData If FALSE (default), the two functions
#' are expected to make use of the bounds of the unit-square. Some-
#' times the two functions cover another rectangular area, and for
#' that case you can set this parameter to TRUE.
#' @param numSamples Only used when using the upper bound from data
#' to find the maximum/minimum extents for both functions.
area_diff_2_functions_score <- function(
  useUpperBoundFromData = FALSE,
  numSamples = 1e4
) {
  return(function(f1, f2) {
    yRange <- if (!useUpperBoundFromData) { c(0,1) } else {
      temp <- stat_diff_2_functions(
        f1 = f1, f2 = f2, numSamples = numSamples)
      range(temp$dataF1, temp$dataF2, na.rm = TRUE)
    }
    yExtent <- abs(yRange[2] - yRange[1])
    
    score <- 1 - (area_diff_2_functions(f1, f2)$value / yExtent)
    
    # We do this to make sure the score is within [0,1], because the
    # actual yExtent may be larger than what we found, due to a limited
    # number of samples we took.
    return(min(1, max(0, score)))
  })
}


#' Samples from both function within [0,1] and returns the
#' indices used and the values from either function. These
#' two vectors can then be compared statistically.
stat_diff_2_functions <- function(f1, f2, statFunc = function(a, b) {
  suppressWarnings({
    stats::cor(x = a, y = b, use = "pairwise.complete.obs")
  })
}, numSamples = 1e4) {
  #indices <- sort(stats::runif(n = numSamples, min = 0, max = 1))
  indices <- seq(0, 1, length.out = numSamples)
  
  d1 <- sapply(indices, f1)
  d2 <- sapply(indices, f2)
  
  return(list(
    dataF1 = d1,
    dataF2 = d2,
    indices = indices,
    value = statFunc(d1, d2)
  ))
}

stat_diff_2_functions_cov <- function(f1, f2, numSamples = 1e4) {
  return(stat_diff_2_functions(f1 = f1, f2 = f2, statFunc = function(x, y) {
    return(stats::cov(x = x, y = y, use = "pairwise.complete.obs"))
  }, numSamples = numSamples))
}

stat_diff_2_functions_cor <- function(f1, f2, numSamples = 1e4) {
  return(stat_diff_2_functions(f1 = f1, f2 = f2, statFunc = function(x, y) {
    return(stats::cor(x = x, y = y, use = "pairwise.complete.obs"))
  }, numSamples = numSamples))
}

#' Computes a score based on Pearson-, Kendall- or Spearman
#' correlation. Returns a stat_diff-style function with f1,
#' f2.
#' 
#' @param requiredSign must be one of 1, 0, -1. If 1, then
#' only positive correlations will lead to score > 0. Likewise,
#' for -1, only negative correlations will yield a score > 0.
#' If 0, then any correlation, positive or negative, produces
#' a score > 0.
#' @param corrType one of "pearson", "kendall", "spearman"
#' @param allowReturnNA cor() returns NA in some cases, for
#' example when at least one of the vectors has an SD of 0. An
#' NA score is useless, but it might be still worth detecting
#' and dealing with it. If NA is encountered and this is NOT
#' set to TRUE, this score will throw an error!
#' @param numSamples amount of samples to use
#' @return function with parameters f1, f2.
stat_diff_2_functions_cor_score <- function(
  requiredSign = c(1, 0, -1)[1],
  corrType = c("pearson", "kendall", "spearman")[1],
  allowReturnNA = FALSE,
  numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- suppressWarnings({
      switch (
        corrType,
        "pearson"  = stat_diff_2_functions_cor(f1 = f1, f2 = f2, numSamples = numSamples),
        "kendall"  = stat_diff_2_functions_cor_kendall(f1 = f1, f2 = f2, numSamples = numSamples),
        "spearman" = stat_diff_2_functions_cor_spearman(f1 = f1, f2 = f2, numSamples = numSamples),
        {
          stop(paste0("Don't know ", corrType, "."))
        }
      )$value
    })
    
    if (is.na(temp)) {
      if (allowReturnNA) {
        return(temp)
      } else {
        stop("Correlation resulted in NA.")
      }
    }
    
    switch (paste0(requiredSign),
      "1"  = max(0, temp),
      "-1" = abs(min(0, temp)),
      "0"  = abs(temp),
      {
        stop(paste0("Cannot handle ", requiredSign, "."))
      }
    )
  })
}

stat_diff_2_functions_cor_kendall <- function(f1, f2, numSamples = 1e4) {
  return(stat_diff_2_functions(f1 = f1, f2 = f2, statFunc = function(x, y) {
    return(stats::cor(x = x, y = y, use = "pairwise.complete.obs", method = "kendall"))
  }, numSamples = numSamples))
}

stat_diff_2_functions_cor_spearman <- function(f1, f2, numSamples = 1e4) {
  return(stat_diff_2_functions(f1 = f1, f2 = f2, statFunc = function(x, y) {
    return(stats::cor(x = x, y = y, use = "pairwise.complete.obs", method = "spearman"))
  }, numSamples = numSamples))
}

stat_diff_2_functions_var <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  temp$value <- stats::var(temp$dataF1 - temp$dataF2, na.rm = TRUE)
  return(temp)
}

#' Computes the sample standard deviation.
stat_diff_2_functions_sd <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  temp$value <- stats::sd(temp$dataF1 - temp$dataF2, na.rm = TRUE)
  return(temp)
}

stat_diff_2_functions_mae <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  temp$value <- mean(abs(temp$dataF1[idx] - temp$dataF2[idx]))
  return(temp)
}

stat_diff_2_functions_rmse <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  temp$value <- Metrics::rmse(temp$dataF1[idx], temp$dataF2[idx])
  return(temp)
}

#' Computes a score based on the sample standard deviation, or
#' the sample variance, or the MAE or RMSE. The returned score
#' is normalized w.r.t. the used method's upper bound in the unit-
#' square, so it is always in the range [0,1]. Returns a stat_diff-
#' style function with f1, f2.
#' 
#' @param use one of "sd", "var", "mae", "rmse"
#' @param useUpperBoundFromData If FALSE (default), the two functions
#' are expected to make use of the bounds of the unit-square. Some-
#' times the two functions cover another rectangular area, and for
#' that case you can set this parameter to TRUE.
#' @param numSamples number of samples to use
#' @return function with parameters f1, f2
stat_diff_2_functions_sd_var_mae_rmse_score <- function(
  use = c("sd", "var", "mae", "rmse")[1],
  useUpperBoundFromData = FALSE,
  numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- switch (use,
      "sd"   = stat_diff_2_functions_sd(f1 = f1, f2 = f2, numSamples = numSamples),
      "var"  = stat_diff_2_functions_var(f1 = f1, f2 = f2, numSamples = numSamples),
      "mae"  = stat_diff_2_functions_mae(f1 = f1, f2 = f2, numSamples = numSamples),
      "rmse" = stat_diff_2_functions_rmse(f1 = f1, f2 = f2, numSamples = numSamples),
      {
        stop(paste0("Don't know ", corrType, "."))
      }
    )
    
    yRange <- if (!useUpperBoundFromData) { c(0,1) } else {
      range(temp$dataF1, temp$dataF2, na.rm = TRUE)
    }
    yExtent <- abs(yRange[2] - yRange[1])
    
    upperBound <- switch (use,
      "sd"   = sqrt(1/2) * yExtent,
      "var"  = (sqrt(1/2) * yExtent)^2,
      "mae"  = yExtent,
      "rmse" = yExtent
    )
    
    1 - (temp$value / upperBound)
  })
}

#' Conveniently, this function's co-domain is [0,1] for f1,f2
#' defined in the unit-square. Computing a score can be done
#' by 1 - frechet-distance.
#' @note The Fréchet distance is very costly and its runtime
#' increases exponentially with 'numSamples'. Using a value
#' much larger than 100 is not recommended. 50 seems to be
#' an acceptable trade-off. However, for grid-searches a
#' value larger than 50 is recommended.
stat_diff_2_functions_frechet <- function(f1, f2, numSamples = 1e2) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  temp$value <- SimilarityMeasures::Frechet(
    traj1 = matrix(data = stats::na.exclude(temp$dataF1), ncol = 1),
    traj2 = matrix(data = stats::na.exclude(temp$dataF2), ncol = 1))
  return(temp)
}

#' Calculates the approximate arc-length of both functions using
#' @seealso {pracma::poly_length()}. Returns both lengths as well
#' as the ratio between f1's and f2's length. The ratio is < 0 if
#' f1's length is longer than f2's; > 0, otherwise. Ideally, the
#' value hence is 0. The arc-length itself has no upper bound.
#' 
#' Note that numSamples should ideally be >= 1e5!
stat_diff_2_functions_arclen <- function(f1, f2, numSamples = 1e5) {
  temp <- stat_diff_2_functions(f1, f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  x <- seq(0, 1, len = sum(idx))
  
  arcLen1 <- pracma::poly_length(x = x, y = temp$dataF1[idx])
  arcLen2 <- pracma::poly_length(x = x, y = temp$dataF2[idx])
  
  temp$arcLen1 <- arcLen1
  temp$arcLen2 <- arcLen2
  temp$value <- 1 - (arcLen1 / arcLen2)
  return(temp)
}

#' Computes a score based on comparing the arc-lengths of two
#' curves in the unit-square. The score is directly proportional
#' to the ratio of the functions' arc-lengths. A score of one
#' hence means that both have the same length. Returns a stat_diff-
#' style function with f1, f2.
#' 
#' @param requiredSign must be one of 1, 0, -1. This is the sign
#' as expected to be returned from stat_diff_2_functions_arclen.
#' If, f1's arc-length is longer than f2's, that function will
#' return a negative value. If f1's arc-length is shorter than
#' f2's, returns a positive value. If the expected signs differ,
#' the score is 0.
#' @param numSamples number of samples to use
#' @return function with parameters f1, f2
stat_diff_2_functions_arclen_score <- function(
  requiredSign = c(0, 1, -1)[1],
  numSamples = 1e5
) {
  return(function(f1, f2) {
    temp <- stat_diff_2_functions_arclen(
      f1 = f1, f2 = f2, numSamples = numSamples)
    
    ratio <- temp$arcLen1 / temp$arcLen2
    if (ratio > 1) {
      ratio <- 1 / ratio
    }
    
    switch(paste0(requiredSign),
      "0"  = ratio,
      "1"  = if (temp$value < 0) 0 else ratio,
      "-1" = if (temp$value > 0) 0 else ratio
    )
  })
}

#' Numerically approximates the gradient for both functions
#' over the sampled interval. The result is one vector for
#' each function with the slope at each evaluated x.
stat_diff_2_functions_grad <- function(f1, f2, numSamples = 1e4) {
  m <- matrix(ncol = 2, nrow = numSamples)
  s <- seq(0, 1, len = numSamples)
  
  for (i in 1:numSamples) {
    m[i, ] <- tryCatch({
      c(
        numDeriv::grad(f1, x = s[i]),
        numDeriv::grad(f2, x = s[i])
      )
    }, error = function(cond) c(NA, NA))
  }
  
  useRows <- stats::complete.cases(m)
  if ((1 - (sum(useRows) / numSamples)) > (1 / (numSamples / 10))) {
    stop("Cannot obtain sufficiently enough gradients.")
  }
  
  list(
    dataF1 = m[useRows, 1],
    dataF2 = m[useRows, 2]
  )
}


#' Performs linear regression on both functions. Returns
#' everything and already contains some additional properties:
#' - angle: the angle between the two LMs in radians. It
#'   is in the range [pi,-pi]. It is positive if f1 has a
#'   steeper slope than f2; negative, otherwise.
#' - intersect: NA if the two LMs do not intersect in the
#'   unit square. Otherwise, the x-coordinate of it.
#' - lm1, lm2: the LMs for f1 and f2. Those can be used
#'   for, e.g., extracting the residuals.
stat_diff_2_functions_lm <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  itv <- seq(0, 1, len = sum(idx))
  
  lm1 <- stats::lm(
    formula = y ~ x, data = list(x = itv, y = temp$dataF1[idx]))
  lm2 <- stats::lm(
    formula = y ~ x, data = list(x = itv, y = temp$dataF2[idx]))
  
  # The idea behind these two vectors is this:
  # - The two functions have a support of [0,1]
  # - if we let them go through 0,0 then we just need to
  #   calculate the value at 1 to get the slope
  vec1 <- c(1,
            stats::predict(lm1, newdata = list(x = c(1))) -
            stats::predict(lm1, newdata = list(x = c(0))))
  vec2 <- c(1,
            stats::predict(lm2, newdata = list(x = c(1))) -
            stats::predict(lm2, newdata = list(x = c(0))))
  
  # This would always be positive. However, we want to return
  # a negative value if slope(f1) < slope(f2)
  isNeg <- if (stats::coef(lm1)[["x"]] < stats::coef(lm2)[["x"]]) -1 else 1
  
  # Also include if the two models cross over (their lines intersect).
  # stats::uniroot() throws an error if 
  int <- sort(rootSolve::uniroot.all(function(z) {
    stats::predict(lm1, newdata = list(x = z)) - stats::predict(lm2, newdata = list(x = z))
  }, interval = c(0, 1)))
  
  arcCosAngle <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  # To increase numeric stability, as sometimes the calculated arcCos angle
  # maybe, e.g., 1+1e-16.
  arcCosAngle <- max(-1, min(1, arcCosAngle))
  
  temp$value <- list(
    lm1 = lm1,
    lm2 = lm2,
    angle = isNeg * acos(arcCosAngle),
    intersect = if (length(int) == 0) NA else int
  )
  
  return(temp)
}


#' @param expectAngleDeg NA if no expectation, otherwise an angle in degrees
#' that is in the range [-pi, pi]
#' @param expectIntersect NA if no expection/check, otherwise requires or not
#' requires an intersect of the two LM's linear functions within the unit square.
#' @param numSamples number of samples to use
#' @return function with parameters f1, f2.
stat_diff_2_functions_lm_score <- function(
  expectAngleDeg = NA,
  expectIntersect = c(NA, TRUE, FALSE)[1],
  numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- stat_diff_2_functions_lm(f1, f2, numSamples = numSamples)
    
    if (!is.na(expectIntersect) &&
        ((expectIntersect && is.na(temp$value$intersect)) ||
        (!expectIntersect && !is.na(temp$value$intersect)))) {
      return(0) # score is 0 if no intersect but expected (or vice versa)
    } else if (!is.na(expectAngleDeg)) {
      if (expectAngleDeg < -180 || expectAngleDeg > 180) {
        stop(paste0("The expectedAngle is not -180* <= a <= 180*: ", expectAngleDeg))
      }
      
      actualAngleDeg <- temp$value$angle / pi * 180
      diffAngleDeg <- expectAngleDeg - actualAngleDeg
      return(1 - abs(diffAngleDeg / 360))
    } else {
      stop("Need to either score intersect or angle between.")
    }
  })
}


#' Similar to \code{stat_diff_2_functions_lm_score}, this function
#' computes a score based on the two LMs of two functions. Instead
#' of expecting an angle, here we expect a maximum angle between
#' those two functions.
#' 
#' @param maxAngleBetween expected angle between the two LMs; must
#' be specified as degrees in the range [-180, 180]
#' @param requireSign If TRUE, then the deviation must have the same
#' sign as the maxAngleBetween. If TRUE and the sign is not met, the
#' resulting score is 0. If false, the sign is discarded, effectively
#' moving the scoring into the range [0, 180].
#' @param numSamples number of samples to use
#' @return function with parameters f1, f2.
stat_diff_2_functions_lm_score2 <- function(
  maxAngleBetween,
  requireSign = TRUE,
  numSamples = 1e4
) {
  return(function(f1, f2) {
    if (maxAngleBetween <= -180 || maxAngleBetween >= 180) {
      stop(paste0("The maxAngleBetween is not -180* < a < 180*: ", maxAngleBetween))
    }
    
    temp <- stat_diff_2_functions_lm(f1, f2, numSamples = numSamples)
    diffAngleDeg <- temp$value$angle / pi * 180
    
    if (!requireSign) {
      maxAngleBetween <- abs(maxAngleBetween)
      diffAngleDeg <- abs(diffAngleDeg)
    }
    
    if ((maxAngleBetween < 0 && (diffAngleDeg < maxAngleBetween || diffAngleDeg > 0)) ||
        (maxAngleBetween > 0 && (diffAngleDeg > maxAngleBetween || diffAngleDeg < 0))) {
      return(0)
    }
    
    return(1 - abs(diffAngleDeg / maxAngleBetween))
  })
}


#' This function should be preferred over the other one,
#' whenever there is doubt that the given functions f1,f2
#' do not integrate to exactly one over the support [0,1].
stat_diff_2_functions_symmetric_KL_sampled <- function(f1, f2, numSamples = 1e4) {
  temp1 <- stat_diff_2_functions_philentropy_sampled(
    f1 = f1, f2 = f2, numSamples = numSamples,
    method = "kullback-leibler", skipNormalize = FALSE)
  temp2 <- stat_diff_2_functions_philentropy_sampled(
    f1 = f2, f2 = f1, numSamples = numSamples,
    method = "kullback-leibler", skipNormalize = FALSE)
  
  temp2$value <- unname(temp2$value + temp1$value) # symmetric KL!
  temp2
  
  # temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  # idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  # 
  # vec1 <- temp$dataF1[idx]
  # vec1 <- vec1 / sum(vec1)
  # vec2 <- temp$dataF2[idx]
  # vec2 <- vec2 / sum(vec2)
  # 
  # temp$value <- suppressMessages({
  #   philentropy::KL(
  #     x = rbind(vec1, vec2), test.na = FALSE, unit = "log") +
  #   philentropy::KL(
  #     x = rbind(vec2, vec1), test.na = FALSE, unit = "log")
  # })
  # 
  # return(temp)
}


#' This is the symmetric KL-divergence
#' Note that both functions must return strictly
#' positive (including 0) values because of their
#' use in the logarithm.
#' 
#' @note WARNING: This function may not work with
#' z-normalized functions!
stat_diff_2_functions_symmetric_KL <- function(
  f1, f2, numSamples = 1e4, sampleOnError = TRUE,
  useCubintegrate = TRUE
) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2)
  tol <- if (useCubintegrate) 0 else 1e-8
  
  PQ <- function(x) {
    f1x <- f1(x)
    f2x <- f2(x)
    
    f1x[is.na(f1x) | f1x < 0] <- 0
    f2x[is.na(f2x) | f2x < 0] <- 0
    
    if (all(f1x == 0) || all(f2x == 0)) {
      return(rep(0, length(x)))
    }
    return(f1x * log(f1x / f2x))
  }
  QP <- function(x) {
    f1x <- f1(x)
    f2x <- f2(x)
    
    f1x[is.na(f1x) | f1x < 0] <- 0
    f2x[is.na(f2x) | f2x < 0] <- 0
    
    if (all(f1x == 0) || all(f2x == 0)) {
      return(rep(0, length(x)))
    }
    return(f2x * log(f2x / f1x))
  }
  
  
  temp$value <- tryCatch({
    if (useCubintegrate) {
      cubature::cubintegrate(
        f = PQ, lower = tol, upper = 1 - tol)$integral +
      cubature::cubintegrate(
        f = QP, lower = tol, upper = 1 - tol)$integral
    } else {
      stats::integrate(
        f = PQ, lower = tol, upper = 1 - tol,
        subdivisions = 10^log10(numSamples),
        stop.on.error = FALSE)$value +
      stats::integrate(
        f = QP, lower = tol, upper = 1 - tol,
        subdivisions = 10^log10(numSamples),
        stop.on.error = FALSE)$value
    }
  }, error = function(cond) {
    if (sampleOnError) {
      s <- seq(0, 1, len=numSamples)
      useVals <- sapply(s, PQ) + sapply(s, QP)
      
      if ((sum(is.na(useVals)) / numSamples) > (1 / (numSamples / 10))) {
        stop("Cannot sample sufficient amount of values for KL-divergence.")
      }
      return(mean(abs(stats::na.exclude(useVals))))
    }
    warning(cond)
    return(NA)
  })
  
  return(temp)
}


#' This function should be preferred over the other one,
#' whenever there is doubt that the given functions f1,f2
#' do not integrate to exactly one over the support [0,1].
stat_diff_2_functions_symmetric_JSD_sampled <- function(f1, f2, numSamples = 1e4) {
  return(stat_diff_2_functions_philentropy_sampled(
    f1 = f1, f2 = f2, numSamples = numSamples, method = "jensen-shannon"))
}


#' https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
#' 
#' @note WARNING: This function may not work with
#' z-normalized functions!
stat_diff_2_functions_symmetric_JSD <- function(
  f1, f2, numSamples = 1e4, sampleOnError = TRUE,
  useCubintegrate = TRUE
) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  tol <- if (useCubintegrate) 0 else 1e-8
  
  M <- function(x) 1/2 * (f1(x) + f2(x))
  PM <- function(x) {
    f1x <- f1(x)
    mx <- M(x)
    
    f1x[is.na(f1x) | f1x < 0] <- 0
    mx[is.na(mx) | mx < 0] <- 0
    
    # Note that log(0) is -Inf!
    # Also, we must not divide by 0..
    # We must not pass values < 0 to log either..
    if (all(f1x == 0) || all(mx == 0)) {
      return(rep(0, length(x)))
    }
    
    return(f1x * log(f1x / mx))
  }
  QM <- function(x) {
    f2x <- f2(x)
    mx <- M(x)
    
    f2x[is.na(f2x) | f2x < 0] <- 0
    mx[is.na(mx) | mx < 0] <- 0
    
    if (all(f2x == 0) || all(mx == 0)) {
      return(rep(0, length(x)))
    }
    
    return(f2x * log(f2x / mx))
  }
  
  JSD <- function(x) 1/2 * PM(x) + 1/2 * QM(x)
  
  temp$value <- tryCatch({
    if (useCubintegrate) {
      cubature::cubintegrate(
        f = JSD, lower = tol, upper = 1 - tol)$integral
    } else {
      stats::integrate(
        f = JSD, lower = tol, upper = 1 - tol, subdivisions = 10^log10(numSamples))$value
    }
  }, error = function(cond) {
    if (sampleOnError) {
      useVals <- sapply(seq(0, 1, len=numSamples), JSD)
      if ((sum(is.na(useVals)) / numSamples) > (1 / (numSamples / 10))) {
        stop("Cannot sample sufficient amount of values from JSD.")
      }
      return(mean(abs(stats::na.exclude(useVals))))
    }
    warning(cond)
    return(NA)
  })
  
  return(temp)
}


#' Computes a score based on the Jensen-Shannon divergence of two functions
#' in the unit-square (or at least with co-domain [0,1]), seen as proper
#' probability density functions. The divergence appears to capture more
#' properties than the other low-level scores and is often sufficient for
#' fitting models. Returns a stat_diff- style function with f1, f2, sampleOnError.
#' 
#' @param sensitivityExponent an exponent to exponentiate the computed score
#' by. The JS-divergence tends to produce relatively high scores, so that a
#' somewhat high exponent helps to spread them out more.
#' @param useSampledJSD if true, uses stat_diff_2_functions_symmetric_JSD_sampled
#' instead of the older and less robust stat_diff_2_functions_symmetric_JSD. Less
#' robust because it assumes that both functions integrate to 1, and it should
#' not be relied on the user having ascertained that. The newer JSD metric however
#' resamples both functions to ensure that property.
#' @param numSamples number of samples to use
#' @return function with parameters f1, f2
stat_diff_2_functions_symmetric_JSD_score <- function(
  sensitivityExponent = 5, useSampledJSD = TRUE, numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- if (useSampledJSD)
      stat_diff_2_functions_symmetric_JSD_sampled(f1 = f1, f2 = f2, numSamples = numSamples)
      else stat_diff_2_functions_symmetric_JSD(
        f1 = f1, f2 = f2,
        numSamples = numSamples, sampleOnError = if (useSampledJSD) NA else TRUE)
    
    # New: We pull the root of the JSD to make it a true
    # statistical distance metric. Note this also moves
    # its upper bound from log(2) to sqrt(log(2)) ~ 0.83.
    unname((1 - (sqrt(temp$value) / sqrt(log(2))))^sensitivityExponent)
  })
}


#' Samples the same values (range) from both f1,f2 and then normalizes
#' the obtained vectors to sum up to 1. If f1,f2 are proper PDFs, that
#' should be the case anyway (i.e., this function, like the others,
#' assumes that the given functions represent proper PDFs with the support
#' of [0,1]). It always uses the log-unit (ln).
#' 
#' @source {https://cran.r-project.org/web/packages/philentropy/vignettes/Information_Theory.html}
#' @param skipNormalize default FALSE, should be set to true when computing
#' distances that do not interpret the sampled data as probability vectors,
#' such as "euclidean" or "manhattan". Since we mostly use this method for
#' vectors of probabilities, FALSE is the default.
stat_diff_2_functions_philentropy_sampled <- function(
  f1, f2, numSamples = 1e4,
  method = philentropy::getDistMethods()[1],
  skipNormalize = FALSE
) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  
  vec1 <- temp$dataF1[idx]
  if (!skipNormalize) {
    vec1 <- vec1 - min(vec1) + .1 # if z-normalized, otherwise integral is 0!
    s1 <- sum(vec1)
    if (s1 > 0) {
      vec1 <- vec1 / s1
    } else {
      warning("Samples of f1 sum to 0.")
    }
  }
  
  vec2 <- temp$dataF2[idx]
  if (!skipNormalize) {
    vec2 <- vec2 - min(vec2) + .1
    s2 <- sum(vec2)
    if (s2 > 0) {
      vec2 <- vec2 / s2
    } else {
      warning("Samples of f2 sum to 0.")
    }
  }
  
  temp$value <- suppressMessages({
    philentropy::distance(
      x = matrix(data = c(vec1, vec2), nrow = 2, byrow = TRUE),
      method = method, unit = "log")
  })
  
  return(temp)
}


#' Calculates the discrete Cross-Entropy H(P,Q) by sampling from
#' f1,f2, then normalizing the sampled data to sum up to 1. The
#' value returned is in Shannon-bits (log2).
stat_diff_2_functions_cross_entropy <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2) & temp$dataF1 > 0 & temp$dataF2 > 0
  
  vec1 <- temp$dataF1[idx]
  vec2 <- temp$dataF2[idx]
  
  vec1 <- vec1 - min(vec1) # if z-normalized, otherwise integral is 0!
  s1 <- sum(vec1)
  if (s1 > 0) {
    vec1 <- vec1 / s1
  }
  
  vec2 <- vec2 - min(vec2)
  s2 <- sum(vec2)
  if (s2 > 0) {
    vec2 <- vec2 / s2
  }
  
  temp$value <- -sum(vec1 * log2(vec2))
  return(temp)
}


#' Returns the Mutual Information of both functions, in (Shannon) bits.
#' 
#' Note: Also returns the entropy for either function.
#' Note: Also returns the Joint Entropy.
#' 
#' If f1/X and f2/Y were to be mutually independent, the MI would be 0.
#' The closer either entropy is to the mutual information, the more like
#' the two functions are, as either could resemble the other. So it may
#' be a good idea to maximize the ratio between one of the entropies and
#' the mutual information, as it can maximally be 1.
#' 
#' Note that this function uses log2() internally, so that entropy and
#' Mutual Information are expressed as bits. Obtaining the amount of
#' information that can be encoded is done by 2^(entropy in bits).
stat_diff_2_functions_mutual_information <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2) & temp$dataF1 > 0 & temp$dataF2 > 0
  u <- "log2"
  
  vec1 <- temp$dataF1[idx]
  vec2 <- temp$dataF2[idx]
  
  vec1 <- vec1 - min(vec1) # if z-normalized, otherwise integral is 0!
  s1 <- sum(vec1)
  if (s1 > 0) {
    vec1 <- vec1 / s1
  }
  
  vec2 <- vec2 - min(vec2)
  s2 <- sum(vec2)
  if (s2 > 0) {
    vec2 <- vec2 / s2
  }
  
  vecJ <- vec1 * vec2
  vecJ <- vecJ / sum(vecJ)
  
  temp$entropy1 <- philentropy::H(vec1, unit = u)
  temp$entropy2 <- philentropy::H(vec2, unit = u)
  #temp$jointEntropy <- philentropy::JE(x = vecJ, unit = u)
  #temp$value <- philentropy::MI(x = vec1, y = vec2,  xy = vecJ, unit = u)
  temp$jointEntropy <- philentropy::MI(x = vec1, y = vec2,  xy = vecJ, unit = u)
  temp$value <- philentropy::JE(x = vecJ, unit = u)
  
  return(temp)
}

stat_diff_2_functions_mutual_information_manual <- function(f1, f2, numSamples = 1e3) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  
  vec1 <- temp$dataF1[idx]
  vec2 <- temp$dataF2[idx]
  
  vec1 <- vec1 - min(vec1) # if z-normalized, otherwise integral is 0!
  s1 <- sum(vec1)
  if (s1 > 0) {
    vec1 <- vec1 / s1
  } else stop()
  
  vec2 <- vec2 - min(vec2)
  s2 <- sum(vec2)
  if (s2 > 0) {
    vec2 <- vec2 / s2
  } else stop()
  
  # idx <- !((vec1 == 0) | (vec2 == 0))
  # vec1 <- vec1[idx]
  # vec2 <- vec2[idx]
  
  temp$vec1 <- vec1
  temp$vec2 <- vec2
  temp
}


#' Computes a score based on the mutual information, using the entropy
#' of either function. Returns a stat_diff-style function with f1, f2.
#' 
#' @param symmetry One of 0, -1 or 1. The entropy of the first function
#' divided by the mutual information is returned when using -1. The one
#' of the second divided by the mutual information if using 1. When
#' using the value 0, both of these are multiplied with each other, so
#' that the score becomes symmetric (commutative). The symmetric version
#' should be used whenever f1 and f2 are exchangeable. -1 should be used
#' when the goal is to find out how well f2 resembles f1, and 1 for fin-
#' ding out how well f1 resembles f2.
#' @param useBits The mutual information and the entropies are computed
#' using Shannon-bits, and the score then is computed by putting these
#' into relation. If you want to compare actual amounts of information,
#' set this to false. This then results in replacing the log2 bits-value
#' with 2 to the power of it. Comparisons of actual amounts of infor-
#' mation are more sensitive.
#' @param sensitivityExponent Like the JSD score, this score tends to
#' produce relatively high scores. To increase the sensitivity, these
#' scores can be spread out more by exponentiating them. If 'useBits'
#' was set to FALSE, this exponent should be lower, as the score is
#' already more sensitive. The default exponents are a somewhat good
#' match for when the entropy and mutual information have about ten
#' to fifteen bits.
#' @param numSamples number of samples to use. Careful! This score is
#' naturally very sensitive w.r.t. how many samples are drawn. While
#' the actual amount is of less importance, it is crucial to use the
#' same amount of samples across all comparisons made.
#' @return function with parameters f1, f2
stat_diff_2_functions_mutual_information_score <- function(
  symmetry = c(0, -1, 1)[1],
  useBits = FALSE,
  sensitivityExponent = if (useBits) 5 else 2,
  numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- stat_diff_2_functions_mutual_information(
      f1 = f1, f2 = f2, numSamples = numSamples)
    
    mi <- temp$value
    e1 <- temp$entropy1
    e2 <- temp$entropy2
    temp$mi <- mi
    
    if (!useBits) {
      mi <- 2^mi
      e1 <- 2^e1
      e2 <- 2^e2
    }
    
    r1 <- e1 / mi
    r2 <- e2 / mi
    
    temp$value <- switch (paste0(symmetry),
      "0"  = r1 * r2,
      "1"  = r1,
      "-1" = r2
    )^sensitivityExponent
    
    temp
  })
}


stat_diff_2_functions_mutual_information_score2 <- function(
  numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- stat_diff_2_functions_mutual_information(
      f1 = f1, f2 = f2, numSamples = numSamples)
    
    I_xy <- temp$value
    H_x <- temp$entropy1
    H_y <- temp$entropy2
    
    temp$C_xy <- I_xy / H_y
    temp$C_yx <- I_xy / H_x
    
    temp$R <- I_xy / (H_x + H_y)
    temp$R_max <- min(H_x, H_y) / (H_x + H_y)
    
    temp$U <- 2 * temp$R
    
    # total correlation
    temp$tc <- I_xy / min(H_x, H_y)
    
    temp
  })
}


#' Extract four paths from a DTW-alignment: The two warping-functions;
#' one for the reference, one for the query. Also extracts paths for
#' the warping function applied to the reference and to the query. These
#' last two can be compared to the original reference or query to assess
#' the goodness of match/fit.
#' 
#' Note that this function does not alter any data or ranges thereof, so
#' that if the data is transformed into functions or patterns, it is
#' recommended that it is scaled together, before it is put into a common
#' range.
extract_warping_from_dtw <- function(dtwAlign, signalRef, signalQuery) {
  return(list(
    # The warping function for the Reference.
    warpRef = suppressWarnings({
      dtw::warp(d = dtwAlign, index.reference = TRUE)
    }),
    # The warping function for the Query.
    warpQuery = suppressWarnings({
      dtw::warp(d = dtwAlign, index.reference = FALSE)
    }),
    # The Reference after warping the query to it.
    warpTransRef = signalRef[dtwAlign$index2],
    # The Query after warping the reference to it.
    warpTransQuery = signalQuery[dtwAlign$index1]
  ))
}


#' Estimates the warping function and overlays its linear regression.
#' Then scales both function together to be in the unit square. These
#' two functions can then be used with metrics such as stat_diff_2_functions_cor.
#'
#' The '_np' variants of the functions imply 'no plateaus', i.e.,
#' those are the same functions with their horizontal plateaus being
#' removed.
#'
#' You may request an optimized align between the original warping
#' function and its non-plateau version. The optimization moves the
#' non-plateau version on the y-axis to minimize the residual sum of squares.
#' 
#' If given the reference- and query-signal, also returns the two warping-
#' functions (one for reference, one for query) as known from the 3-way plot.
pattern_approxfun_warp <- function(dtwAlign, includeOptimizedAlign = TRUE, signalRef = NULL, signalQuery = NULL) {
  warpData <- dtwAlign$index2
  len <- length(warpData)
  
  warpLm <- stats::lm(
    formula = y ~ x, data = data.frame(
      x = 1:len, y = warpData))
  
  linData <- stats::predict(warpLm, newdata = data.frame(x = 1:len))
  
  # The two vectors 'warpData' and 'linData' need to be scaled together.
  ranges <- range(warpData, linData)
  extent <- ranges[2] - ranges[1]
  
  warpData <- warpData - ranges[1]
  warpData <- warpData / extent
  linData <- linData - ranges[1]
  linData <- linData / extent
  
  ex <- extract_signal_from_window(
    dtwAlign = dtwAlign, window = dtwAlign$index2, idxMethod = "smooth")
  
  f_warp_org <- pattern_approxfun(yData = dtwAlign$index2)
  f_warp_np <- pattern_approxfun(yData = ex$data)
  f_warp_np_opt <- NULL
  
  if (includeOptimizedAlign) {
    optFunc <- function(offsetLm) {
      tempFunc <- function(x) f_warp_np(x) + offsetLm
      area_diff_2_functions(f_warp_org, tempFunc)$value
    }
    
    f_warp_np_opt <- tryCatch({
      optR <- stats::optim(
        par = 0, fn = optFunc, method = "BFGS")
      
      if (optR$convergence != 0) {
        stop("The optimization did not converge.")
      }
      function(x) f_warp_np(x) + optR$par
    }, error = function(cond) {
      warning(cond)
      NULL
    })
  }
  
  res <- list(
    # These two functions have been scaled and translated together, as the
    # linear regression of the warping function may be either larger or
    # smaller than the warping function, and therefore the limits for the
    # warping function have to be adjusted. Only both functions together
    # have a domain of [0,1] and co-domain of [0,1].
    f_warp = pattern_approxfun(yData = warpData, yLimits = range(warpData)),
    f_lm = pattern_approxfun(yData = linData, yLimits = range(linData)),
    
    # The original warping function, squeezed into [0,1] (x&y). We need it
    # for comparisons against its non-plateau variant. Comparisons should
    # be either made against 'f_warp_np' or 'f_warp_np_opt', and nothing else.
    f_warp_org = f_warp_org,
    # The original warping function w/o plateaus, also in [0,1] for x&y
    # This is necessary because the length of this function is <= f_warp_org
    f_warp_np = f_warp_np,
    # The warping function w/o plateaus, and adjusted on the y-axis to
    # minimize the residual sum of squares between it and the original
    # warping function.
    f_warp_np_opt = f_warp_np_opt,
    
    # These two represent the warping functions as known from 3-way plot,
    # both for the reference and the query.
    f_warp_ref = NULL,
    f_warp_query = NULL
  )
  
  if (!missing(signalRef) && !missing(signalQuery)) {
    ex_dtw <- extract_warping_from_dtw(
      dtwAlign = dtwAlign, signalRef = signalRef, signalQuery = signalQuery)
    
    res$f_warp_ref <- pattern_approxfun(yData = ex_dtw$warpRef)
    res$f_warp_query <- pattern_approxfun(yData = ex_dtw$warpQuery)
  }
  
  return(res)
}


#' Uses 'warpTransRef' and 'warpTransQuery' from @seealso {extract_warping_from_dtw()}
#' to generate functions in the unit-square. Generates two pairs of functions, one
#' pair for the reference and one for the query. These pairs can be used to assess the
#' goodness of match.
get_dtw_scorable_functions <- function(dtwAlign, signalRef = NULL, signalQuery = NULL) {
  ex <- extract_warping_from_dtw(
    dtwAlign = dtwAlign, signalRef = signalRef, signalQuery = signalQuery)
  
  safe_divide <- function(numerator, divisor, threshold = 1e-20) {
    if (divisor < threshold) {
      return(numerator)
    }
    res <- numerator / divisor
    if (any(is.na(res))) {
      stop("Division produces NA.")
    }
    res
  }
  
  
  range_ref <- range(signalRef, ex$warpTransRef)
  extent_ref <- range_ref[2] - range_ref[1]
  data_ref <- safe_divide(signalRef - range_ref[1], extent_ref)
  data_ref_warp <- safe_divide(ex$warpTransRef - range_ref[1], extent_ref)
  
  f_ref <- pattern_approxfun(yData = data_ref, yLimits = c(0, 1))
  f_ref_warp <- pattern_approxfun(yData = data_ref_warp, yLimits = c(0, 1))
  
  range_q <- range(signalQuery, ex$warpTransQuery)
  extent_q <- range_q[2] - range_q[1]
  data_q <- safe_divide(signalQuery - range_q[1], extent_q)
  data_q_warp <- safe_divide(ex$warpTransQuery - range_q[1], extent_q)
  
  f_q <- pattern_approxfun(yData = data_q, yLimits = c(0, 1))
  f_q_warp <- pattern_approxfun(yData = data_q_warp, yLimits = c(0, 1))
  
  return(list(
    f_ref = f_ref,
    f_ref_warp = f_ref_warp,
    
    f_query = f_q,
    f_query_warp = f_q_warp
  ))
}


densitySafe <- function(data, ratio = 1, kernel = "gauss", bw = "SJ", from = 0, to = 1, safeVal = 0) {
  d <- stats::density(x = data, from = from, to = to, kernel = kernel, bw = bw, n = 2^max(10, ceiling(log2(length(data)))))
  f <- stats::approxfun(x = d$x, y = d$y)
  r <- range(d$x)
  f1 <- Vectorize(function(x) {
    if (x < r[1] || x > r[2]) safeVal else f(x) * ratio
  })
  attributes(f1) <- list(
    min = r[1], max = r[2], ratio = ratio, ymax = max(d$y) * ratio,
    x = d$x, y = d$y)
  f1
}


poly_autofit <- function(yData, xData = 1:length(yData), maxDegree = 5, method = c("AIC", "BIC", "Radj")[1], selectBy = c("findBest", "stopEarly")[1]) {
  
  # We want to minimize the cost, and while a higher
  # R.adjusted is better, lower AIC/BIC are better.
  minimizeMult <- if (method == "Radj") -1 else 1
  methodFunc <- function(model) {
    if (method == "Radj") {
      return(stats::summary.lm(model)$adj.r.squared)
    } else if (method == "AIC") {
      return(stats::AIC(model))
    } else {
      return(stats::BIC(model))
    }
  }
  
  startVal <- .Machine$double.xmax
  degree <- 1
  degreeBest <- 1
  tempModel <- NULL
  
  while (TRUE && degree <= maxDegree) {
    tempModel <- tryCatch(expr = {
      stats::lm(
        formula = yData ~ stats::poly(xData, degree = degree))
    }, error = function(cond) FALSE)
    
    if (is.logical(tempModel) && tempModel == FALSE) {
      break
    }
    
    newVal <- methodFunc(tempModel) * minimizeMult
    if (newVal < startVal) {
      startVal <- newVal
      degreeBest <- degree
    } else if (selectBy == "stopEarly") {
      return(tempModel)
    }
    
    degree <- degree + 1
  }
  
  # We tested all, so let's return the best:
  stats::lm(
    formula = yData ~ stats::poly(xData, degree = degreeBest))
}

poly_autofit_max <- function(x, y, startDeg = 10) {
  p <- NULL
  while (TRUE) {
    temp <- tryCatch({
      stats::lm(formula = y ~ stats::poly(x = x, degree = startDeg))
    }, error = function(cond) FALSE)
    
    if (is.logical(temp) && temp == FALSE) {
      break
    } else {
      p <- temp
      startDeg <- startDeg + 1
    }
  }
  p
}


# The following 4 functions were implemented from [@xi2000bearing].


RMS <- function(vec) {
  l <- length(vec)
  m <- mean(vec)
  sqrt(1 / l * sum((vec - m)^2))
}

Kurtosis <- function(vec) {
  l <- length(vec)
  m <- mean(vec)
  (1 / l * sum((vec - m)^4)) / RMS(vec)
}

Peak <- function(vec) {
  .5 * (max(vec) - min(vec))
}

ImpulseFactor <- function(vec) {
  l <- length(vec)
  Peak(vec) / (1 / l * sum(abs(vec)))
}


#' Treats both functions as signals and uniformly samples from
#' them. Then computes one or more signal measures for either
#' sampled signal and puts both into relation. Each relation is
#' moved into [0,1], where 1 means both signals have the same
#' relation for the computed property. Returns a stat_diff-
#' style function with f1, f2.
#' 
#' @param use One of "RMS", "Kurtosis", "Peak", "ImpulseFactor"
#' or "all". If "all", computes all of the previous and returns
#' a named vector instead of a single score.
#' @param requiredSign One of 0, 1, -1. Determines the required
#' ratio for each computed property. 1 means that the value for
#' f1 must be greater than the one for f2, and -1 vice versa. 0
#' means that it does not matter. If a required sign is not met,
#' the resulting score for the property is 0.
#' @param numSamples number of samples to use
#' @return function with parameters f1, f2
stat_diff_2_functions_signals_score <- function(
  use = c("all", "RMS", "Kurtosis", "Peak", "ImpulseFactor")[1],
  requiredSign = c(0, 1, -1)[1],
  numSamples = 1e4
) {
  return(function(f1, f2) {
    temp <- stat_diff_2_functions(
      f1 = f1, f2 = f2, numSamples = numSamples)
    idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
    d1 <- temp$dataF1[idx]
    d2 <- temp$dataF2[idx]
    
    fac1 <- c()
    fac2 <- c()
    useAll <- use == "all"
    if (useAll || "RMS" %in% use) {
      fac1 <- c(fac1, "RMS" = RMS(d1))
      fac2 <- c(fac2, "RMS" = RMS(d2))
    }
    if (useAll || "Kurtosis" %in% use) {
      fac1 <- c(fac1, "Kurtosis" = Kurtosis(d1))
      fac2 <- c(fac2, "Kurtosis" = Kurtosis(d2))
    }
    if (useAll || "Peak" %in% use) {
      fac1 <- c(fac1, "Peak" = Peak(d1))
      fac2 <- c(fac2, "Peak" = Peak(d2))
    }
    if (useAll || "ImpulseFactor" %in% use) {
      fac1 <- c(fac1, "ImpulseFactor" = ImpulseFactor(d1))
      fac2 <- c(fac2, "ImpulseFactor" = ImpulseFactor(d2))
    }
    
    frac <- fac1 / fac2
    temp <- 1 - frac
    if (requiredSign == 1) {
      temp[temp < 0] <- 0
    } else if (requiredSign == -1) {
      temp[temp > 0] <- 0
    }
    
    sapply(frac, function(r) if (r < 1) r else 1 / r)
  })
}


#' Creates a custom scoring method with a callback that is
#' passed the arguments specified by the three dots (...).
#' stat_diff-style scores accepts the two functions f1,f2,
#' but a custom score does not have to. If however it does
#' have these parameter names, then f1,f2 will be passed
#' as first two arguments, followed by the three dots.
#' Instead of specifying args using the 3 dots, you can
#' also make a closure over the values you need.
#' 
#' @param callback function that may accept any (even zero)
#' arguments. If it accpets f1,f2, then those are passed.
#' @param ... additional arguments that should be passed to
#' the callback. NOTE: these are evaluated immediately.
#' @return stat_diff-style function that accepts f1,f2.
stat_diff_custom_score <- function(callback, ...) {
  stopifnot(is.function(callback))
  argNames <- methods::formalArgs(callback)
  
  hasF1F2 <- all(c("f1", "f2") %in% argNames)
  # This is UTMOST important, as otherwise, the scope might not
  # be longer accessible at later run-time!
  useArgs <- list(...)
  stopifnot(length(useArgs) == 0 || all(argNames %in% names(useArgs)))
  
  
  return(function(f1, f2) {
    if (hasF1F2) {
      useArgs$f1 <- f1
      useArgs$f2 <- f2
    }
    return(do.call(what = callback, args = useArgs))
  })
}



#' Used to create a non-linear transformation R -> R for scores,
#' with both domain and co-domain [0,1].
#' 
#' @param t A threshold (0,1) until which (including) the function
#' will behave linearly. Values larger than it are exponentiated
#' using the parameter \code{k}. Values until \code{t} are calculated
#' as \code{x * t^(k - 1)}.
#' @param k An exponent for the non-linear transformation that is
#' applied for all x > \code{t}, i.e., \code{x^k}.
#' @return A vectorized function R -> R, that non-linearly trans-
#' formes the given score, into [0,1] again
create_penalizeScore <- function(t = .4, k = 2.2) {
  fac <- t^(k-1)
  return(Vectorize(function(x) {
    if (x <= t) x * fac else x^k
  }))
}


custom_interval_length_score <- function(mlm, interval, useA) {
  # Create a custom score based on the current interval length.
  # We have two options:
  # - A: measure the ratio of the lengths of the current interval to the pattern's interval:
  #      the closer it is to the reference, the higher the score.
  # - B: return a score that is higher if the chosen interval is longer:
  #      The longer, the better, that's it!
  
  curr <- mlm$getCurrentIntervalRange(intervalIndexOrName = interval)
  currExt <- curr[2] - curr[1]
  
  if (useA) {
    # We have to access the original boundaries first.
    intIdx <- which(mlm$intervalNames == interval)
    ref <- c(
      if (intIdx == 1) 0 else mlm$refBoundaries[intIdx - 1],
      if (intIdx == length(mlm$intervalNames)) 1 else mlm$refBoundaries[intIdx])
    
    refExt <- ref[2] - ref[1]
    ratio <- refExt / currExt
    ratio <- if (ratio > 1) 1 / ratio else ratio
    
    return(ratio)
  } else {
    # The best way would be to divide the current length
    # by the maximum theoretical length, based on the lin.
    # ineq. constraints. This is not necessarily trivial,
    # so we use the absolute theoretical maximum of 1.
    return(currExt)
  }
}



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



plot_project_data <- function(data, boundaries = c()) {
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_line(aes(color = t), size = .75) +
    theme_light() +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05), limits = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.10), limits = c(0, 1)) +
    theme(axis.text.x = element_text(angle = -45, vjust = 0)) +
    labs(color = "Activity") +
    scale_color_brewer(palette = "Set1")
  
  for (b in boundaries) {
    p <- p + geom_vline(xintercept = b, color = "blue", size = .5)
  }
  
  p
}



