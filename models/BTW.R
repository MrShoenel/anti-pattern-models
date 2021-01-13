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




M_final <- function(theta_b_org, theta_b, r, f, num_samples = 1e3) {
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
    f_primes[[paste0(iIdx)]] <- (function(q) {
      
      b_os <- theta_b_org[q]
      b_oe <- theta_b_org[q + 1]
      b_qs <- theta_b[q]
      b_qe <- theta_b[q + 1]
      de_o <- b_oe - b_os
      de_q <- b_qe - b_qs
      frac <- if (de_q == 0) de_o / .Machine$double.eps else de_o / de_q
      
      function(x) {
        temp <- (x - b_qs) * frac + b_os
        if (is.na(temp)) {
          stop(c(x, b_qs, frac, b_os))
        }
        f(temp)
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
    
    paste0(iIdx)
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
    int_idx = factor(x = int_idx, levels = paste0(1:(length(theta_b) - 1)), ordered = TRUE)
  )
}


L_final <- function(
  theta_b_org,
  theta_b,
  r, f,
  useR1 = TRUE,
  useR2 = TRUE
) {
  res <- M_final(theta_b_org = theta_b_org, theta_b = theta_b, r = r, f = f)
  
  loss <- sum((res$y - res$y_hat)^2)
  
  if (useR1 || useR2) {
    vts <- utils::head(theta_b, -1)
    vte <- utils::tail(theta_b, -1)
    X_r <- range(theta_b_org)
    X_q <- range(theta_b)
    
    
    if(useR1) {
      R <- Vectorize(function(x) if (x < 0) 0 else x)
      H <- Vectorize(function(x) if (x < 0) 0 else 1)
      p_phi <- function(v1, v2) H(R(v1 - v2)) * R(v1 - v2) + H(R(v2 - v1)) * R(v2 - v1)
      
      vts_o <- sort(vts)
      vts_or <- rev(vts_o)
      vte_o <- sort(vte)
      vte_or <- rev(vte_o)
      
      v_s <- p_phi(vts, vts_o)
      v_e <- p_phi(vte, vte_o)
      u_s <- p_phi(vts_o, vts_or)
      u_e <- p_phi(vte_o, vte_or)
      
      loss <- loss - log(1 - (sum(v_s + v_e) / sum(u_s + u_e)))
    }
    
    
    if (useR2) {
      loss <- loss - log((X_q[2] - X_q[1]) / (X_r[2] - X_q[1]))
    }
  }
  
  loss
}


L_final_grad <- function(
  theta_b_org,
  theta_b,
  r, f, f_p,
  useR1 = TRUE,
  useR2 = TRUE
) {
  res <- M_final(theta_b_org = theta_b_org, theta_b = theta_b, r = r, f = f)
  
  # \theta_s, \theta_e
  ts <- utils::head(theta_b_org, -1)
  te <- utils::tail(theta_b_org, -1)
  
  # \vartheta_s, \vartheta_e
  vts <- utils::head(theta_b, -1)
  vte <- utils::tail(theta_b, -1)
  Q <- as.numeric(levels(res$int_idx))
  
  theta_new <- c()
  
  for (q in Q) {
    grad_bqs <- 0
    grad_bqe <- 0
    
    bos <- ts[q]
    boe <- te[q]
    
    bqs <- vts[q]
    bqe <- vte[q]
    
    delta_o <- boe - bos
    x_q <- res$X[res$int_idx == q]
    y_q <- res$y[res$int_idx == q]
    
    for (i in seq_len(length(x_q))) {
      x_i <- x_q[i]
      temp_i <- if (bqs - bqe == 0) 1 else bqs - bqe # TODO: check what makes sense: 1 or .Machine$double.eps
      psi_i <- bos + ((delta_o * (bqs - x_i)) / temp_i)
      
      big_fac <- 2 * delta_o * (y_q[i] - f(psi_i)) * f_p(psi_i) / (temp_i)^2
      
      grad_bqs <- grad_bqs + big_fac * (bqe - x_i)
      grad_bqe <- grad_bqe + big_fac * (x_i - bqs)
    }
    
    if (q <= max(Q)) {
      theta_new <- c(theta_new, grad_bqs)
    }
    if (q == max(Q)) {
      theta_new <- c(theta_new, grad_bqe)
    }
  }
  
  theta_new
}











