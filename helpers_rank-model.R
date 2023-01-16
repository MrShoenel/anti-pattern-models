# File taken from https://raw.githubusercontent.com/MrShoenel/R-rank-model/master/R/helpers.R

curve2 <- function(func, from, to, col = "black", lty = 1, lwd = 1, add = FALSE, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL, ...) {
  f <- function(x) func(x)
  curve(expr = f, from = from, to = to, col = col, lty = lty, lwd = lwd, add = add, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main = main, ... = ...)
}


fit_cont_parametric <- function(data, test_distr = c("beta", "cauchy", "gamma", "logis", "lnorm", "norm", "unif", "weibull"), test_alpha = 0.05, jitter_data = TRUE, over_scale = 0) {
  # For the Beta distribution, data needs to be in range [0,1].
  # For the Gamma distribution, data needs to be in range (0, inf]
  
  data <- as.numeric(data)
  if (jitter_data) {
    expo <- round(log10(mean(data)) - 7)
    data <- data + runif(n = length(data), min = -10^expo, max = 10^expo)
  }
  best <- NULL
  
  e <- options("show.error.messages")
  options(show.error.messages = FALSE)
  for (distr in test_distr) {
    result <- tryCatch(expr = {
      suppressMessages({
        suppressWarnings(expr = {
          fitdistrplus::fitdist(data = data, distr = distr)
        })
      })
    }, error = function(cond) cond, finally = {
      options(show.error.messages = e)
    })
    
    if ("error" %in% class(result)) {
      # Cannot fit this distribution. Try next.
      next
    }
    
    dist_params <- setNames(as.list(unname(result$estimate)), names(result$estimate))
    if (over_scale > 0) {
      if ("scale" %in% names(dist_params)) {
        dist_params$scale <- dist_params$scale * (1 + over_scale)
      } else if ("min" %in% names(dist_params) && "max" %in% names(dist_params)) {
        ext <- dist_params$max - dist_params$min
        dist_params$min <- dist_params$min - over_scale / 2 * ext
        dist_params$max <- dist_params$max + over_scale / 2 * ext
      }
    }
    temp <- rlang::duplicate(dist_params)
    temp$x <- quote(data)
    temp$y <- paste0("p", distr) # E.g., "pnorm"
    test <- do.call(what = ks.test, args = temp)
    temp$p.value <- test$p.value
    temp$statistic <- test$statistic
    
    if (test$p.value >= test_alpha && (is.null(best) || temp$statistic < best$statistic)) {
      best <- temp
      best$distr <- distr
      best$dist_params <- dist_params
    }
  }
  
  if (is.null(best)) {
    stop("Cannot fit a parametric continuous distribution to the given data.")
  }
  
  pdf <- get(paste0("d", best$distr))
  cdf <- get(paste0("p", best$distr))
  ppf <- get(paste0("q", best$distr))
  
  to_dist_args <- function(x, use = c("x", "q", "p")) {
    templ <- rlang::duplicate(x = best$dist_params)
    templ[[match.arg(arg = use)]] <- x
    templ
  }
  
  list(
    "distr" = best$distr,
    "dist_params" = best$dist_params,
    "p_value" = best$p.value,
    "statistic" = best$statistic,
    "range" = range(data),
    
    "pdf" = function(x) {
      do.call(what = pdf, args = to_dist_args(x = x, use = "x"))
    },
    "cdf" = function(q) {
      do.call(what = cdf, args = to_dist_args(x = q, use = "q"))
    },
    "ppf" = function(p) {
      do.call(what = ppf, args = to_dist_args(x = p, use = "p"))
    }
  )
}


wasserstein_distance <- function(sample1, sample2, continuous = FALSE, discrete_samples = 1e4) {
  e1 <- stats::ecdf(sample1)
  e2 <- stats::ecdf(sample2)
  min_ <- min(c(sample1, sample2))
  max_ <- max(c(sample1, sample2))
  ext <- max_ - min_
  
  if (continuous) {
    return(cubature::cubintegrate(f = function(x) abs(e1(x) - e2(x)), lower = min_, upper = max_, maxEval = 2e4)$integral)
  }
  
  d <- 0
  for (x in seq(from = min_, to = max_, length.out = discrete_samples)) {
    d <- d + abs(e1(x) - e2(x))
  }
  d / discrete_samples * ext
}

kl_div_approx <- function(sample1, sample2) {
  d1 <- stats::density(sample1)
  d2 <- stats::density(sample2)
  y1 <- d1$y / sum(d1$y)
  y2 <- d2$y / sum(d2$y)
  
  sum(y1 * (log2(y1) - log2(y2)))
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


