M <- function(theta_b_org, theta_b, r, f, num_samples = 1e3) {
  stopifnot(length(theta_b_org) == length(theta_b))
  stopifnot(is.function(r) && is.function(f))
  stopifnot(num_samples >= length(theta_b))
  
  X <- seq(from = min(theta_b), to = max(theta_b), length.out = num_samples)
  y <- Vectorize(r)(X)
  
  # Contrary to the pseudo-code, we will create a list of
  # all intervals and scaled+translated functions to look
  # them up, that's much faster.
  
  f_primes <- list()
  for (iIdx in seq_len(length(theta_b) - 1)) {
    f_primes[[iIdx]] <- (function(int_idx) {
      start_org <- theta_b_org[int_idx]
      end_org <- theta_b_org[int_idx + 1]
      start <- theta_b[int_idx]
      end <- theta_b[int_idx + 1]
      extent <- end - start
      if (extent == 0) {
        extent <- .Machine$double.eps
      }
      
      function(x) {
        f( ((x - start) / extent) * (end_org - start_org) + start_org )
      }
    })(iIdx)
  }
  
  
  # Now we can sample from all intervals. Note that this
  # Model stays in an interval for as long as it is valid.
  # If a subsequent interval overlaps, it will not switch
  # into it.
  
  determine_interval_for_x <- function(x) {
    start <- min(theta_b)
    end <- max(theta_b)
    stopifnot(x >= start && x <= end)
    
    # Note that intervals are [start, end), except for the
    # last interval, which [start, end].
    start_last <- theta_b[length(theta_b) - 1]
    if (x >= start_last) {
      return(length(theta_b) - 1)
    }
    
    iIdx <- 1
    for (bIdx in seq_len(length(theta_b) - 1)) {
      if (x >= theta_b[bIdx] && x < theta_b[bIdx + 1]) {
        break
      } else {
        iIdx <- iIdx + 1
      }
    }
    
    if (iIdx > (length(theta_b) - 1)) {
      print(x)
      print(theta_b)
      stop("42!")
    }
    
    iIdx
  }
  
  y_hat <- rep(NA, num_samples)
  for (xIdx in seq_len(num_samples)) {
    iIdx <- determine_interval_for_x(X[xIdx])
    y_hat[xIdx] <- f_primes[[iIdx]](X[xIdx])
  }
  
  stopifnot(!any(is.na(c(y, y_hat))))
  list(
    X = X, # return the support used
    y = y,
    y_hat = y_hat
  )
}


M_grad <- function(theta_b_org, theta_b, r, f_grad, num_samples = 1e3) {
  M(theta_b_org = theta_b_org, theta_b = theta_b, r = r, f = f_grad, num_samples = num_samples)
}

M_grad2 <- function(theta_b_org, theta_b, r, f_grad2, num_samples = 1e3) {
  M(theta_b_org = theta_b_org, theta_b = theta_b, r = r, f = f_grad2, num_samples = num_samples)
}


M_new <- function(theta_b_org, theta_b, r, f, num_samples = 1e3) {
  stopifnot(length(theta_b_org) == length(theta_b))
  stopifnot(is.function(r) && is.function(f))
  stopifnot(num_samples >= length(theta_b))
  
  X <- seq(from = min(theta_b), to = max(theta_b), length.out = num_samples)
  y <- Vectorize(r)(X)
  
  # Contrary to the pseudo-code, we will create a list of
  # all intervals and scaled+translated functions to look
  # them up, that's much faster.
  
  f_primes <- list()
  for (iIdx in seq_len(length(theta_b) - 1)) {
    f_primes[[iIdx]] <- (function(q) {
      
      start_org <- theta_b_org[q]
      end_org <- theta_b_org[q + 1]
      start_q <- theta_b[q]
      end_q <- theta_b[q + 1]
      
      function(x) {
        d <- end_q - start_q
        if (d == 0) {
          return(0)
        }
        x_rel <- (x - start_q) / d
        x_o <- start_org + (x_rel * (end_org - start_org))
        f(x_o)
      }
    })(iIdx)
  }
  
  
  # Now we can sample from all intervals. Note that this
  # Model stays in an interval for as long as it is valid.
  # If a subsequent interval overlaps, it will not switch
  # into it.
  
  determine_interval_for_x <- function(x) {
    start <- min(theta_b)
    end <- max(theta_b)
    stopifnot(x >= start && x <= end)
    
    # Note that intervals are [start, end), except for the
    # last interval, which [start, end].
    start_last <- theta_b[length(theta_b) - 1]
    if (x >= start_last) {
      return(length(theta_b) - 1)
    }
    
    iIdx <- 1
    for (bIdx in seq_len(length(theta_b) - 1)) {
      if (x >= theta_b[bIdx] && x < theta_b[bIdx + 1]) {
        break
      } else {
        iIdx <- iIdx + 1
      }
    }
    
    if (iIdx > (length(theta_b) - 1)) {
      print(x)
      print(theta_b)
      stop("42!")
    }
    
    iIdx
  }
  
  y_hat <- rep(NA, num_samples)
  int_idx <- rep(NA, num_samples)
  for (xIdx in seq_len(num_samples)) {
    iIdx <- determine_interval_for_x(X[xIdx])
    y_hat[xIdx] <- f_primes[[iIdx]](X[xIdx])
    int_idx[xIdx] <- iIdx
  }
  
  stopifnot(!any(is.na(c(y, y_hat))))
  data.frame(
    X = X, # return the support used
    y = y,
    y_hat = y_hat,
    int_idx = int_idx
  )
}


