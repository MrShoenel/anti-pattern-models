
M_updated <- function(theta_b_org, theta_b, r, f, num_samples = 1e3, zNormalize = FALSE, valueForZeroLenIntervals = NA_real_) {
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
  
  X_r <- range(theta_b_org)
  x_r <- seq(from = X_r[1], to = X_r[2], length.out = num_samples)
  Q <- seq_len(length.out = length(theta_l))
  
  y <- Vectorize(r)(x_r)
  y_hat <- rep(NA_real_, num_samples)
  q_idx <- rep(NA, num_samples)
  X_q <- c()

  i <- 1
  for (q in Q) {
    last_q <- q == rev(Q)[1]
    b_os <- theta_b_org[q]
    b_oe <- theta_b_org[q + 1]
    b_qs <- theta_b[q]
    b_qe <- theta_b[q + 1]
    
    x_r_supp <- if (last_q) {
      x_r[x_r >= b_os & x_r <= b_oe]
    } else {
      x_r[x_r >= b_os & x_r < b_oe]
    }
    supp_len <- length(x_r_supp)
    stopifnot(supp_len > 0)
    
    l_q <- theta_l[q]
    phi_q <- if (q == 1) 0 else sum(theta_l[1:(q - 1)])
    de_o <- theta_l_org[q]
    
    samp_func <- Vectorize(function(x) {
      temp <- de_o * (x - s - phi_q) / l_q + b_os
      f(temp)
    })
    
    # Let's check the new query interval.
    if (l_q == 0) {
      # It was requested to have a length of 0, but we need to
      # sample 'supp_len' number of elements, which we cannot.
      # We repeat 'valueForZeroLenIntervals' this element instead.
      y_hat[i:(i + supp_len - 1)] <- valueForZeroLenIntervals
      X_q <- c(X_q, rep(b_qs, supp_len))
    } else {
      # An interval of negative length means that it goes backward.
      # This case should be caught by the regularizers. However, we
      # can still return some meaningful values, by sampling from
      # the range of the interval, and then reversing the values.
      
      # This also works if b_qs > b_qe ..
      x_q <- seq(from = b_qs, to = b_qe, length.out = supp_len)
      X_q <- c(X_q, x_q) # append support used (before reversing)
      
      y_hat[i:(i + supp_len - 1)] <- samp_func(x_r_supp)
    }
    
    q_idx[i:(i + supp_len - 1)] <- q
    
    i <- i + supp_len
  }
  
  if (zNormalize) {
    y <- (y - mean(y)) / sd(y)
    y_hat <- (y_hat - mean(y_hat)) / sd(y_hat)
  }
  
  data.frame(
    x_r = x_r,
    X_q = X_q,
    y = y,
    y_hat = y_hat,
    q_idx = factor(x = q_idx, levels = paste0(1:length(theta_l)), ordered = TRUE)
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
  
  # ######## KL
  # idx_not_NA <- !(is.na(res$y) | is.na(res$y_hat))
  # if (sum(idx_not_NA) == 0) {
  #   stop("All data is NA.")
  # }
  # y_not_NA <- res$y[idx_not_NA]
  # y_hat_not_NA <- res$y_hat[idx_not_NA]
  # 
  # y_not_NA <- y_not_NA - min(y_not_NA)
  # y_not_NA <- y_not_NA / sum(y_not_NA)
  # y_hat_not_NA <- y_hat_not_NA - min(y_hat_not_NA)
  # y_hat_not_NA <- y_hat_not_NA / sum(y_hat_not_NA)
  # 
  # idx_not_0 <- !(y_not_NA == 0 | y_hat_not_NA == 0)
  # if (sum(idx_not_0) == 0) {
  #   stop("All data is 0.")
  # }
  # 
  # y_not_NA <- y_not_NA[idx_not_0]
  # y_hat_not_NA <- y_hat_not_NA[idx_not_0]
  # 
  # dkl_pq <- sum(y_not_NA * log(y_not_NA / y_hat_not_NA))
  # dkl_qp <- sum(y_hat_not_NA * log(y_hat_not_NA / y_not_NA))
  # 
  # loss <- weightErr * (dkl_pq + dkl_qp)
  
  
  
  
  ######### 
  #numData <- sum(complete.cases(res$y_hat))
  #loss <- weightErr * log(1 + sum(na.omit(res$y - res$y_hat)^2) / numData)
  
  
  ######## RSS
  numData <- sum(complete.cases(res$y_hat))
  loss <- loss + weightErr * log(1 + sum(na.omit(res$y - res$y_hat)^2))
  
  loss
}