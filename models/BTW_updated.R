
M_updated <- function(theta_b_org, theta_b, r, f, num_samples = 1e3, zNormalize = FALSE) {
  stopifnot(length(theta_b_org) == length(theta_b))
  stopifnot(is.function(r) && is.function(f))
  stopifnot(num_samples >= length(theta_b))
  
  # The reference boundaries need to be sorted. Otherwise,
  # there would be intervals with negative length.
  stopifnot(all.equal(theta_b_org, sort(theta_b_org)))
  
  s <- theta_b[1]
  theta_l <- sapply(seq_len(length(theta_b) - 1), function(b) {
    theta_b[b + 1] - theta_b[b]
  })
  theta_l_org <- sapply(seq_len(length(theta_b_org) - 1), function(b) {
    theta_b_org[b + 1] - theta_b_org[b]
  })
  
  
  X <- seq(from = min(theta_b), to = max(theta_b), length.out = num_samples)
  y <- Vectorize(r)(X)
  
  # Contrary to the pseudo-code, we will create a list of
  # all intervals and scaled+translated functions to look
  # them up, that's much faster.
  
  f_primes <- list()
  for (iIdx in seq_len(length(theta_b) - 1)) {
    f_primes[[paste0(iIdx)]] <- (function(q) {
      
      lq <- if(theta_l[q] == 0) 1 else theta_l[q]
      phi_q <- if (q == 1) 0 else sum(theta_l[1:(q - 1)])
      b_os <- theta_b_org[q]
      de_o <- theta_l_org[q]
      
      Vectorize(function(x) {
        temp <- de_o * (x - s - phi_q) / lq + b_os
        f(temp)
      })
    })(iIdx)
  }
  
  
  # Now we can sample from all intervals. Note that this
  # Model stays in an interval for as long as it is valid.
  # If a subsequent interval overlaps, it will not switch
  # into it.
  start <- min(theta_b)
  end <- max(theta_b)
  start_last <- theta_b[length(theta_b) - 1]
  
  determine_interval_for_x <- function(x) {
    stopifnot(x >= start && x <= end)
    
    # Note that intervals are [start, end), except for the
    # last interval, which is [start, end].
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
      warning("44!")
      return("NA")
    }
    
    paste0(iIdx)
  }
  
  y_hat <- rep(NA, num_samples)
  int_idx <- rep(NA, num_samples)
  for (xIdx in seq_len(num_samples)) {
    iIdx <- determine_interval_for_x(X[xIdx])
    #y_hat[xIdx] <- f_primes[[iIdx]](X[xIdx])
    y_hat[xIdx] <- if (iIdx == "NA") NA_real_ else f_primes[[iIdx]](X[xIdx])
    int_idx[xIdx] <- iIdx
  }
  
  if (zNormalize) {
    y <- (y - mean(y)) / sd(y)
    y_hat <- (y_hat - mean(y_hat)) / sd(y_hat)
  }
  
  # TODO: We have currently allowed this - the error function needs to take care of this!
  #stopifnot(!any(is.na(c(y, y_hat))))
  data.frame(
    X = X, # return the support used
    y = y,
    y_hat = y_hat,
    int_idx = factor(x = int_idx, levels = paste0(1:(length(theta_b) - 1)), ordered = TRUE)
  )
}


L_updated_log <- function(
  theta_b_org,
  theta_b,
  r, f,
  weightErr = 1
) {
  loss <- 0
  res <- M_updated(theta_b_org = theta_b_org, theta_b = theta_b, r = r, f = f)
  
  ######## KL
  idx_not_NA <- !(is.na(res$y) | is.na(res$y_hat))
  if (sum(idx_not_NA) == 0) {
    stop("All data is NA.")
  }
  y_not_NA <- res$y[idx_not_NA]
  y_hat_not_NA <- res$y_hat[idx_not_NA]

  y_not_NA <- y_not_NA - min(y_not_NA)
  y_not_NA <- y_not_NA / sum(y_not_NA)
  y_hat_not_NA <- y_hat_not_NA - min(y_hat_not_NA)
  y_hat_not_NA <- y_hat_not_NA / sum(y_hat_not_NA)

  idx_not_0 <- !(y_not_NA == 0 | y_hat_not_NA == 0)
  if (sum(idx_not_0) == 0) {
    stop("All data is 0.")
  }

  y_not_NA <- y_not_NA[idx_not_0]
  y_hat_not_NA <- y_hat_not_NA[idx_not_0]

  dkl_pq <- sum(y_not_NA * log(y_not_NA / y_hat_not_NA))
  dkl_qp <- sum(y_hat_not_NA * log(y_hat_not_NA / y_not_NA))

  loss <- weightErr * (dkl_pq + dkl_qp)
  
  
  
  
  ######### 
  #numData <- sum(complete.cases(res$y_hat))
  #loss <- weightErr * log(1 + sum(na.omit(res$y - res$y_hat)^2) / numData)
  
  
  ######## RSS
  # loss <- loss + weightErr * log(1 + sum(na.omit(res$y - res$y_hat)^2))
  
  ######### RMSE
  # numData <- sum(complete.cases(res$y_hat))
  # loss <- loss + weightErr * log(1 + sum(na.omit(res$y - res$y_hat)^2))
  
  loss
}