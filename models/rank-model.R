# File taken from https://raw.githubusercontent.com/MrShoenel/R-rank-model/master/R/rank-model.R

# Based on https://arxiv.org/pdf/2111.04682.pdf
smooth_max <- function(x1, x2, alpha = 0.01) {
  (x1 + x2 + sqrt((x1 - x2)^2 + alpha)) / 2
}

smooth_min <- function(x1, x2, alpha = 0.01) {
  (x1 + x2 - sqrt((x1 - x2)^2 + alpha)) / 2
}

hard_sigmoid <- function(x) smooth_min(1, smooth_max(0, x))

sigmoid <- function(x) 1 / (1 + exp(-x))

swish <- function(x, beta=0.25) {
  x * sigmoid(beta * x)
}

make_safe_ppf <- function(ppf, x_min = 1e-220, x_max = 1-1e-16) {
  function(p) {
    p[p <= 0] <- x_min
    p[p >= 1] <- x_max
    ppf(p=p)
  }
}


make_smooth_ecdf <- function(values, slope = 0.025, inverse = FALSE) {
  r <- range(values)
  e <- stats::ecdf(values)
  x <- sort(unique(values))
  y <- e(x)
  if (slope > 0) {
    ext <- r[2] - r[1]
    # Add a slight slope before and after for numeric stability.
    x <- c(r[1] - ext, x, r[2] + ext)
    y <- c(0 - slope, y, 1 + slope)
  }
  
  # Note that the inversed ECDF (the EPPF,) will have an x-range of [0-slope, 1+slope].
  # We do it this way so that we allow the PPF to be called outside its range which may
  # be useful for new, unseen data that is outside of the known range.
  `attributes<-`(x = stats::approxfun(x = if (inverse) y else x, y = if (inverse) x else y, yleft = if (inverse) min(x) else y[1], yright = if (inverse) max(x) else y[length(y)]), value = list(
    "min" = min(values),
    "max" = max(values),
    "range" = range(values),
    
    "slope_min" = min(x),
    "slope_max" = max(x),
    "slope_range" = range(x)
  ))
}


create_model <- function(df_train, x_cols, y_col, cdf_type = c("auto", "gauss", "logis", "ecdf", "smooth", "fit")) {
  df_train <- as.data.frame(df_train)
  cdf_type <- match.arg(cdf_type)
  
  model_info <- list(
    "distr" = list()
  )
  
  make_gauss_cdf <- function(data, over_scale = 0) {
    mean_ <- mean(data)
    sd_ <- (1 + over_scale) * sd(data)
    function(q) pnorm(q = q, mean = mean_, sd = sd_)
  }
  
  make_logis_cdf <- function(data) {
    mean_ <- mean(data)
    sd_ <- mean(data)
    function(q) plogis(q = q, location = mean_, scale = sd_, log.p = FALSE)
  }
  
  ppf <- NULL
  cdfs <- list()
  if (cdf_type == "gauss") {
    for (x_col in x_cols) {
      cdfs[[x_col]] <- make_gauss_cdf(data = df_train[, x_col], over_scale = .1)
    }
    ppf_mean <- mean(df_train[, y_col])
    ppf_sd <- 1.1 * sd(df_train[, y_col]) # Note how we over-scale here as well
    ppf <- make_safe_ppf(ppf = function(p) qnorm(p = p, mean = ppf_mean, sd = ppf_sd))
  } else if (cdf_type == "logis") {
    for (x_col in x_cols) {
      cdfs[[x_col]] <- make_logis_cdf(data = df_train[, x_col])
    }
    ppf_mean <- mean(df_train[, y_col])
    ppf_sd <- sd(df_train[, y_col])
    ppf <- function(p) {
      p[p <= 0] <- 1e-220
      p[p >= 1] <- 1-1e-16
      qlogis(p = p, location = ppf_mean, scale = ppf_sd, log.p = FALSE)
    }
  } else if (cdf_type == "ecdf") {
    for (x_col in x_cols) {
      cdfs[[x_col]] <- stats::ecdf(x = df_train[, x_col])
    }
    ppf <- make_smooth_ecdf(values = df_train[, y_col], slope = 0, inverse = TRUE)
  } else if (cdf_type == "smooth") {
    for (x_col in x_cols) {
      # Create a smoothed ECDF that also has slopes to cope with previously unseen data.
      cdfs[[x_col]] <- make_smooth_ecdf(values = df_train[, x_col])
    }
    ppf <- make_smooth_ecdf(values = df_train[, y_col], inverse = TRUE)
  } else if (cdf_type == "fit") {
    for (x_col in x_cols) {
      the_fit <- fit_cont_parametric(data = df_train[, x_col], over_scale = .1)
      model_info$distr[[x_col]] <- list(
        "distr" = the_fit$distr,
        "dist_params" = the_fit$dist_params,
        "p_value" = the_fit$p_value,
        "statistic" = the_fit$statistic)
      cdfs[[x_col]] <- the_fit$cdf
    }
    the_fit <- fit_cont_parametric(data = df_train[, y_col], over_scale = .1)
    model_info$distr[[y_col]] <- list(
      "distr" = the_fit$distr,
      "dist_params" = the_fit$dist_params,
      "p_value" = the_fit$p_value,
      "statistic" = the_fit$statistic)
    ppf_ <- the_fit$ppf
    ppf <- make_safe_ppf(ppf = ppf_)
  } else if (cdf_type == "auto") {
    # Try fit, then fall back to smooth ECDF.
    for (x_col in x_cols) {
      cdfs[[x_col]] <- tryCatch({
        the_fit <- fit_cont_parametric(data = df_train[, x_col], over_scale = .1)
        model_info$distr[[x_col]] <- list(
          "distr" = the_fit$distr,
          "dist_params" = the_fit$dist_params,
          "p_value" = the_fit$p_value,
          "statistic" = the_fit$statistic)
        the_fit$cdf
      }, error = function(cond) {
        make_smooth_ecdf(values = df_train[, x_col])
      })
    }
    ppf <- tryCatch({
      the_fit <- fit_cont_parametric(data = df_train[, y_col], over_scale = .1)
      model_info$distr[[y_col]] <- list(
        "distr" = the_fit$distr,
        "dist_params" = the_fit$dist_params,
        "p_value" = the_fit$p_value,
        "statistic" = the_fit$statistic)
      the_fit$ppf
    }, error = function(cond) {
      make_smooth_ecdf(values = df_train[, y_col], inverse = TRUE)
    })
  }
  
  m <- function(x, df) {
    a_m <- x[1]
    b_m <- x[2] # Model output bias for Sigmoid
    # For each feature (x), we have a weight, and a scale and translate, like w*F(a+b*x)
    # First, there come the weights, then the a's, then the b's
    num_feats <- length(x_cols)
    weights <- x[2+(1:num_feats)]
    a_s <- x[2+num_feats+(1:num_feats)]
    b_s <- x[2+num_feats+num_feats+(1:num_feats)]
    
    res <- c()
    for (rn in rownames(df)) {
      data <- c()
      for (x_col in x_cols) {
        data[x_col] <- cdfs[[x_col]](df[rn, x_col])
      }
      res[rn] <- ppf(hard_sigmoid(a_m + b_m * weights %*% swish(a_s + b_s * as.numeric(data))))
    }
    res
  }
  
  `attributes<-`(m, list(
    "info" = model_info
  ))
}

model_loss <- function(model, x, df, y_col) {
  df <- as.data.frame(df)
  Metrics::rmse(actual = as.numeric(df[, y_col]), predicted = model(x = x, df = df))
}

model_loss_ws <- function(model, x, df, y_col) {
  sample1 <- as.numeric(df[, y_col])
  sample2 <- model(x = x, df = df)
  log(1 + wasserstein_distance(sample1 = sample1, sample2 = sample2, continuous = FALSE))
}

model_loss_kla <- function(model, x, df, y_col) {
  sample1 <- as.numeric(df[, y_col])
  sample2 <- model(x = x, df = df)
  kl_div_approx(sample1 = sample1, sample2 = sample2)
}
