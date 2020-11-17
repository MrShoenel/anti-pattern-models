pattern_approxfun <- function(yData, smooth = FALSE, yLimits = c(0, 1), smoothSpan = .15) {
  if (all(yData == 0)) {
    return(approxfun(
      x = c(0, 1),
      y = c(0, 0)
    ))
  }
  
  if (smooth) {
    temp <- stats::loess.smooth(
      x = 1:length(yData),
      y = yData,
      span = smoothSpan
    )
    
    yData <- temp$y
    xData <- temp$x - min(temp$x)
    xData <- xData / max(xData)
  } else {
    xData <- seq(0, 1, by = 1 / (length(yData) - 1))
  }
  
  yData <- yData - min(yData)
  yData <- yData / max(yData)
  
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


extract_signal_from_window <- function(dtwAlign, window, throwIfFlat = TRUE, idxMethod = c("discrete", "smooth"), smoothCnt = 3) {
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
    for (idx in 1:(min(length(window), length(dtwAlign$index2)) - 1)) {
      if (dtwAlign$index2[idx] < dtwAlign$index2[idx + 1]) {
        indices <- c(indices, c(idx, idx + 1))
      }
    }
    indices <- unique(indices)
  } else if (idxMethod == "smooth") {
    indices <- c()
    for (idx in 1:(min(length(window), length(dtwAlign$index2)) - smoothCnt + 1)) {
      temp <- dtwAlign$index2[idx:(idx + smoothCnt - 1)]
      if (max(temp) > min(temp)) {
        indices <- c(indices, idx:(idx + smoothCnt - 1))
      }
    }
    indices <- unique(indices)
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
      rmse = sqrt(sum(mean(warp_res^2)))
    ),
    warp_rel_resid = list(
      lm = warpLm_rel,
      score = 1 - (mean(abs(warp_rel_res)) * 2),
      sd = stats::sd(warp_rel_res),
      var = stats::var(warp_rel_res),
      mae = mean(abs(warp_rel_res)),
      rmse = sqrt(sum(mean(warp_rel_res^2)))
    ),
    
    window_info = list(
      lm = windowLm,
      mean = mean(window),
      slope = windowLmCoef["x"],
      intercept = windowLmCoef["(Intercept)"]
    ),
    
    indices = indices,
    start = min(indices),
    start_rel = min(indices) / dtwAlign$N,
    end = max(indices),
    end_rel = max(indices) / dtwAlign$N,
    start_ref = min(dtwAlign$index2),
    start_rel_ref = min(dtwAlign$index2) / dtwAlign$M,
    end_ref = max(dtwAlign$index2),
    end_rel_ref = max(dtwAlign$index2) / dtwAlign$M,
    data = window[indices]
  ))
}


area_diff_2_functions <- function(f1, f2) {
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
  
  
  
  # Now, for each pair of intersections, we integrate both
  # functions and sum up the areas.
  
  areas <- c()
  for (intsec in 1:(length(intersections) - 1)) {
    temp <- abs(
      stats::integrate(
        f = f1,
        lower = intersections[intsec],
        upper = intersections[intsec + 1],
        subdivisions = 1e5)$value
      -
        stats::integrate(
          f = f2,
          lower = intersections[intsec],
          upper = intersections[intsec + 1],
          subdivisions = 1e5)$value)
    
    areas <- c(areas, temp)
  }
  
  return(list(
    areas = areas,
    value = sum(areas),
    intersections = intersections
  ))
}


#' Samples from both function within [0,1] and returns the
#' indices used and the values from either function. These
#' two vectors can then be compared statistically.
stat_diff_2_functions <- function(f1, f2, statFunc = stats::cor, numSamples = 1e4) {
  #indices <- sort(stats::runif(n = numSamples, min = 0, max = 1))
  indices <- seq(0, 1, by = 1 / numSamples)
  
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

#' Conveniently, this function's co-domain is [0,1] for f1,f2
#' defined in the unit-square.
stat_diff_2_functions_frechet <- function(f1, f2, numSamples = 1e2) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples, statFunc = function(x, y) {
    
  })
  temp$value <- SimilarityMeasures::Frechet(
    traj1 = matrix(data = stats::na.exclude(temp$dataF1), ncol = 1),
    traj2 = matrix(data = stats::na.exclude(temp$dataF2), ncol = 1))
  return(temp)
}

#' Calculates the approximate arc-length of both functions using
#' @seealso {pracma::poly_length()}. Returns both lengths as well
#' as the ratio between f1's and f2's length. The ratio is < 0 if
#' f1's length is shorter than f2's; > 0, otherwise. Ideally, the
#' value hence is 0. The arc-length itself has no upper bound.
#' 
#' Note that numSamples should ideally be >= 1e5!
stat_diff_2_functions_arclen <- function(f1, f2, numSamples = 1e5) {
  temp <- stat_diff_2_functions(f1, f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  x <- seq(0, 1, len = length(idx))
  
  arcLen1 <- pracma::poly_length(x = x, y = temp$dataF1[idx])
  arcLen2 <- pracma::poly_length(x = x, y = temp$dataF2[idx])
  
  temp$arcLen1 <- arcLen1
  temp$arcLen2 <- arcLen2
  temp$value <- 1 - (arcLen1 / arcLen2)
  return(temp)
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
  int <- rootSolve::uniroot.all(function(z) {
    stats::predict(lm1, newdata = list(x = z)) - stats::predict(lm2, newdata = list(x = z))
  }, interval = c(0, 1))
  
  temp$value <- list(
    lm1 = lm1,
    lm2 = lm2,
    angle = isNeg * acos(sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))),
    intersect = if (length(int) == 0) NA else int
  )
  
  return(temp)
}


#' This function should be preferred over the other one,
#' whenever there is doubt that the given functions f1,f2
#' do not integrate to exactly one over the support [0,1].
stat_diff_2_functions_symmetric_KL_sampled <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  
  vec1 <- temp$dataF1[idx]
  vec1 <- vec1 / sum(vec1)
  vec2 <- temp$dataF2[idx]
  vec2 <- vec2 / sum(vec2)
  
  temp$value <- suppressMessages({
    philentropy::KL(
      x = rbind(vec1, vec2), test.na = FALSE, unit = "log") +
    philentropy::KL(
      x = rbind(vec2, vec1), test.na = FALSE, unit = "log")
  })
  
  return(temp)
}


#' This is the symmetric KL-divergence
#' Note that both functions must return strictly
#' positive (including 0) values because of their
#' use in the logarithm.
stat_diff_2_functions_symmetric_KL <- function(f1, f2, numSamples = 1e4, sampleOnError = TRUE) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2)
  tol <- 1e-8
  
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
    stats::integrate(f = PQ, lower = tol, upper = 1 - tol,
                     subdivisions = 10^log10(numSamples),
                     stop.on.error = FALSE)$value +
    stats::integrate(f = QP, lower = tol, upper = 1 - tol,
                     subdivisions = 10^log10(numSamples),
                     stop.on.error = FALSE)$value
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
stat_diff_2_functions_symmetric_JSD <- function(f1, f2, numSamples = 1e4, sampleOnError = TRUE) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  tol <- 1e-8
  
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
    stats::integrate(
      f = JSD, lower = tol, upper = 1 - tol, subdivisions = 10^log10(numSamples))$value
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


#' Samples the same values (range) from both f1,f2 and then normalizes
#' the obtained vectors to sum up to 1. If f1,f2 are proper PDFs, that
#' should be the case anyway (i.e., this function, like the others,
#' assumes that the given functions represent proper PDFs with the support
#' of [0,1]). It always uses the log-unit (ln).
stat_diff_2_functions_philentropy_sampled <- function(f1, f2, numSamples = 1e4, method = philentropy::getDistMethods()[1]) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2)
  
  vec1 <- temp$dataF1[idx]
  vec1 <- vec1 / sum(vec1)
  vec2 <- temp$dataF2[idx]
  vec2 <- vec2 / sum(vec2)
  
  temp$value <- suppressMessages({
    philentropy::distance(
      x = matrix(data = c(vec1, vec2), nrow = 2, byrow = TRUE),
      method = method, unit = "log")
  })
  
  return(temp)
}


#' Calculates the discrete Cross-Entropy H(P,Q) by sampling from
#' f1,f2, then normalizing the sampled data to sum up to 1.
stat_diff_2_functions_cross_entropy <- function(f1, f2, numSamples = 1e4) {
  temp <- stat_diff_2_functions(f1 = f1, f2 = f2, numSamples = numSamples)
  idx <- !is.na(temp$dataF1) & !is.na(temp$dataF2) & temp$dataF1 > 0 & temp$dataF2 > 0
  
  vec1 <- temp$dataF1[idx]
  vec1 <- vec1 / sum(vec1)
  vec2 <- temp$dataF2[idx]
  vec2 <- vec2 / sum(vec2)
  
  temp$value <- -sum(vec1 * log(vec2))
  return(temp)
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
#' function and its
pattern_approxfun_warp <- function(dtwAlign, includeOptimizedAlign = TRUE) {
  len <- length(dtwAlign$index2)
  
  warpLm <- stats::lm(
    formula = y ~ x, data = data.frame(
      x = 1:len, y = dtwAlign$index2))
  
  warpData <- dtwAlign$index2
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
  
  return(list(
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
    f_warp_np_opt = f_warp_np_opt
  ))
}

