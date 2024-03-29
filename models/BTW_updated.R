Stabilize <- function(f, lb, ub) {
  Vectorize(function(x) {
    if (x < lb) f(lb) else if (x > ub) f(ub) else f(x)
  })
}

R <- Vectorize(function(x) if (x < 0) 0 else x)

H <- Vectorize(function(x) if (x > 0) 1 else 0)


M_updated <- function(theta_b_org, theta_b, r, f, num_samples = 1e3, zNormalize = FALSE, valueForZeroLenIntervals = 0) {
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
  
  Q <- seq_len(length.out = length(theta_l))
  X_r <- range(theta_b_org)
  x_r <- seq(from = X_r[1], to = X_r[2], length.out = num_samples)
  X_q <- range(theta_b)
  x_q <- seq(from = X_q[1], to = X_q[2], length.out = num_samples)
  
  x_q_used <- c()
  
  y <- Vectorize(r)(x_q)
  y_hat <- rep(NA_real_, num_samples)
  q_idx <- rep(NA, num_samples)
  
  i <- 1
  for (q in Q) {
    last_q <- q == rev(Q)[1]
    b_qs <- theta_b[q]
    b_qe <- theta_b[q + 1]
    b_os <- theta_b_org[q]
    b_oe <- theta_b_org[q + 1]
    
    num_samples_in_q <- length(if (last_q) {
      x_r[x_r >= b_os & x_r <= b_oe]
    } else {
      x_r[x_r >= b_os & x_r < b_oe]
    })
    
    q_idx[i:(i + num_samples_in_q - 1)] <- q
    
    if ((b_qe - b_qs) == 0) {
      y_hat[i:(i + num_samples_in_q - 1)] <- valueForZeroLenIntervals
      x_q_used <- c(x_q_used, rep(b_qs, num_samples_in_q))
      i <- i + num_samples_in_q
      next
    }
    
    bq_range <- range(b_qs, b_qe)
    bq_rev <- b_qs > b_qe
    
    x_q_supp <- seq(
      from = bq_range[1], to = bq_range[2], length.out = num_samples_in_q)
    x_q_used <- c(x_q_used, x_q_supp)
    
    l_q <- theta_l[q]
    phi_q <- if (q == 1) 0 else sum(theta_l[1:(q - 1)])
    de_o <- theta_l_org[q]
    
    samp_func <- Vectorize(function(x) {
      use_l_q <- if (l_q == 0) 1 else l_q
      
      temp <- de_o * (x - s - phi_q) / use_l_q + b_os
      f(temp)
    })
    
    use_vals <- samp_func(x_q_supp)
    if (bq_rev) {
      use_vals <- rev(use_vals)
    }
    
    y_hat[i:(i + num_samples_in_q - 1)] <- use_vals
    
    i <- i + num_samples_in_q
  }
  
  f_hat <- stats::approxfun(x = x_q_used, y = y_hat, ties = mean, rule = 2)
  y_hat <- sapply(x_r, f_hat)
  
  if (zNormalize) {
    y <- (y - mean(y)) / sd(y)
    y_hat <- (y_hat - mean(y_hat)) / sd(y_hat)
  }
  
  data.frame(
    X = x_r,
    X_q = x_q_used,
    y = y,
    y_hat = y_hat,
    q_idx = factor(x = q_idx, levels = paste0(1:length(theta_l)), ordered = TRUE)
  )
}


M_final_no_NA <- function(
  theta_b_org, theta_b, r, f, num_samples = 1e3, zNormalize = FALSE
) {
  stopifnot(length(theta_b_org) == length(theta_b))
  stopifnot(is.function(r) && is.function(f))
  stopifnot(num_samples >= length(theta_b))
  
  X_r <- seq(from = min(theta_b_org), to = max(theta_b_org), length.out = num_samples)
  X <- seq(from = min(theta_b), to = max(theta_b), length.out = num_samples)
  y <- Vectorize(r)(X_r)
  
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
  
  temp <- sapply(seq_len(length.out = length(theta_b) - 1), function(t) {
    range(theta_b[t], theta_b[t + 1])
  })
  interval_supports <- matrix(data = temp, ncol = 2, byrow = TRUE)
  
  determine_interval_for_x <- function(x) {
    for (r in seq_len(nrow(interval_supports))) {
      is_last <- r == nrow(interval_supports)
      if (is_last) {
        return(paste0(r))
      } else {
        if (x >= interval_supports[r, 1] && x < interval_supports[r, 2]) {
          return(paste0(r))
        }
      }
    }
    stop("Should never get here!")
  }
  
  y_hat <- rep(NA, num_samples)
  int_idx <- rep(NA, num_samples)
  for (xIdx in seq_len(num_samples)) {
    iIdx <- determine_interval_for_x(X[xIdx])
    y_hat[xIdx] <- f_primes[[iIdx]](X[xIdx])
    int_idx[xIdx] <- iIdx
  }
  
  if (zNormalize) {
    y <- (y - mean(y)) / sd(y)
    y_hat <- (y_hat - mean(y_hat)) / sd(y_hat)
  }
  
  stopifnot(!any(is.na(y_hat)))
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
  model = M_updated,
  theta_w = NULL,
  weightErr = theta_w[1],
  weightR4 = theta_w[2],
  weightR5 = theta_w[3]
) {
  loss_raw <- 0
  loss <- 0
  res <- model(theta_b_org = theta_b_org, theta_b = theta_b, r = r, f = f)
  
  
  if (!missing(theta_w) && is.numeric(theta_w)) {
    # Regularizer for theta_w loss:
    # 1) Check that all weights are within [0,1]
    # 2) Check that all weights sum to 1
    
    # 1)
    lb <- (1 / length(theta_w))^2
    ub <- 1
    temp <- abs(R(lb - theta_w)) + abs(R(theta_w - ub))
    temp <- temp[temp > 0]
    # loss_raw <- loss_raw + sum(1 + temp)^(length(theta_w) * length(temp)) - 1
    # loss <- loss + log(sum(1 + temp)^(length(theta_w) * length(temp)))
    # loss_raw <- loss_raw + (length(theta_w) * sum(temp))^(length(theta_w) * length(temp))
    # loss <- loss - log(1 / (1 + (length(theta_w) * sum(temp))^(length(theta_w) * length(temp))))
    
    
    # ###### WORKS, but log(2) Rest..
    # loss_raw <- loss_raw + (length(theta_w) * sum(1 + temp))^(length(theta_w) * length(temp))
    # loss <- loss + log(1 + (length(theta_w) * sum(1 + temp))^(length(theta_w) * length(temp)))
    # ###### WORKS (alternative), no Rest..
    loss_raw <- loss_raw + (1 + length(theta_w))^(length(temp) * log(1 + sum(1 + temp)))
    loss <- loss + log((1 + length(theta_w))^(length(temp) * log(1 + sum(1 + temp))))
    
    # 2)
    temp <- abs(1 - sum(theta_w)) # ideally, this is 0
    # loss_raw <- loss_raw - log(1 / (1 + temp)^(1 + temp))
    # loss <- loss - log(1 / (1 + temp)^(1 + temp))
    loss_raw <- loss_raw + log(1 + length(theta_w) * temp)
    loss <- loss + log(1 + length(theta_w) * temp)
  } else {
    theta_w <- c(weightErr, weightR4, weightR5)
  }
  
  
  
  
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
  # loss <- loss + abs(theta_w[1]) * (dkl_pq + dkl_qp)
  
  
  
  ####### RSS
  loss_raw <- loss_raw + abs(theta_w[1]) * sum((res$y - res$y_hat)^2)
  loss <- loss + abs(theta_w[1]) * log(1 + sum((res$y - res$y_hat)^2))

  
  
  # ####### Stat-diff
  # f1 <- pattern_approxfun(
  #   yData = res$y,
  #   xData = seq(0, 1, length.out = length(res$y)))
  # f2 <- pattern_approxfun(
  #   yData = res$y_hat,
  #   xData = seq(0, 1, length.out = length(res$y_hat))
  # )
  # a <- area_diff_2_functions_score(useUpperBoundFromData = TRUE)(f1, f2)
  # loss <- loss - log(a)
  # 
  # c <- stat_diff_2_functions_cor_score(allowReturnNA = FALSE)(f1, f2)
  # loss <- loss - log(c)
  # 
  # # l <- stat_diff_2_functions_arclen_score()(f1, f2)
  # # loss <- loss - log(l)
  # # j <- stat_diff_2_functions_symmetric_JSD_score(useSampledJSD = FALSE)(f1, f2)
  # # loss <- loss - log(j)
  # # j <- stat_diff_2_functions_signals_score()(f1, f2)
  # # loss <- loss - log(prod(j))
  
  
  
  
  vts <- utils::head(theta_b, -1)
  vte <- utils::tail(theta_b, -1)
  
  s <- theta_b[1]
  theta_l <- sapply(seq_len(length(theta_b) - 1), function(b) {
    theta_b[b + 1] - theta_b[b]
  })
  theta_l_org <- sapply(seq_len(length(theta_b_org) - 1), function(b) {
    theta_b_org[b + 1] - theta_b_org[b]
  })
  
  # ####### R1:
  # p_phi <- function(v1, v2) H(R(v1 - v2)) * R(v1 - v2) + H(R(v2 - v1)) * R(v2 - v1)
  # 
  # vts_o <- sort(vts)
  # vts_or <- rev(vts_o)
  # vte_o <- sort(vte)
  # vte_or <- rev(vte_o)
  # 
  # v_s <- p_phi(vts, vts_o)
  # v_e <- p_phi(vte, vte_o)
  # u_s <- p_phi(vts_o, vts_or)
  # u_e <- p_phi(vte_o, vte_or)
  # 
  # eps <- .Machine$double.eps
  # temp <- sum(v_s + v_e) / sum(u_s + u_e)
  # temp <- temp * (1 - eps)
  # 
  # loss <- loss - log(1 - temp)
  
  
  # ####### R2:
  # X_r <- range(theta_b_org)
  # X_q <- range(theta_b)
  # eps <- .Machine$double.eps
  # temp <- (X_q[2] - X_q[1]) / (X_r[2] - X_r[1])
  # temp <- eps + temp * (1 - eps)
  # loss <- loss + abs(log(temp))
  
  
  # ###### R3:
  # X_r <- range(theta_b_org)
  # mu <- (X_r[2] - X_r[1]) / (length(theta_b_org) - 1)
  # loss <- loss + log(1 + sum((vte - vts - mu)^2))
  
  
  
  ###### R4 (Box bounds):
  lb <- min(theta_b_org)
  ub <- max(theta_b_org)
  temp <- abs(theta_b[theta_b < lb | theta_b > ub])
  loss_raw <- loss_raw + abs(theta_w[2]) * (sum(1 + temp)^length(temp) - 1)
  loss <- loss + abs(theta_w[2]) * log(sum(1 + temp)^length(temp))
  
  
  ###### R5 (neg Intervals):
  neg_l <- abs(theta_l[theta_l < 0])
  loss_raw <- loss_raw + abs(theta_w[3]) * (sum(1 + neg_l)^length(neg_l) - 1)
  loss <- loss + abs(theta_w[3]) * log(sum(1 + neg_l)^length(neg_l))
  
  
  
  
  
  
  loss
  # log(1 + loss_raw)
}




o_om_r1r2 <- function(x) {
  loss <- 0
  
  # # regularizer box-bounds:
  # temp <- q[, 1] + x
  # temp <- abs(temp[temp < 1/2*pi | temp > 3/2 * pi])
  # if (length(temp) > 0) {
  #   loss <- loss + log(1 + sum(temp)^length(temp))
  # }

  # deviation from expected distance:
  temp <- q[, 1] + x
  temp <- temp - temp[1]
  temp <- temp / (max(temp) - min(temp))
  loss <- loss + log(1 + sum((temp - qx_exp)^2))

  # Tikhonov:
  temp <- q[, 1] + x
  temp <- (temp - mean(temp)) / sd(temp)
  loss <- loss + log(1 + sum(temp^2))
  
  loss
}


o_om_neg_int <- function(x) {
  # Here, we penalize negative intervals.
  loss <- 0
  
  negs <- c()
  trans <- q[, 1] + x
  for (idx in seq_len(length.out = length(x) - 1)) {
    l <- trans[idx + 1] - trans[idx]
    if (l < 0) {
      negs <- c(negs, l)
    }
  }
  
  loss <- loss + log(1 + sum(negs^2)^length(negs))
  
  loss
}

o_om_illegal_trans <- function(x) {
  # Here, we check that each translation t_i > t_{i-1}
  loss <- 0
  
  vio <- c()
  trans <- q[, 1] + x
  for (idx in seq_len(length.out = length(x) - 1)) {
    if (trans[idx] > trans[idx + 1]) {
      vio <- c(vio, trans[idx] - trans[idx + 1])
    }
  }
  
  loss <- loss + log(1 + length(vio)^log(1 + sum((1 + vio)^2)))
  
  # loss <- loss + log(1 + sum((1 + vio)^2)^length(vio))
  
  loss
}


o_om <- function(x, isGrad = FALSE) {
  loss <- 0
  
  data_loss <- 0
  for (idx in seq_len(length.out = length(x))) {
    temp <- q[1:idx, 1] + x[1:idx]
    data_loss <- data_loss + (q[idx, 2] - r(max(temp)))^2
  }
  
  # data-loss:
  loss <- loss + log(1 + data_loss)
  
  # r1,r2:
  loss <- loss + o_om_r1r2(x)
  
  # # R: neg_int
  # loss <- loss + o_om_neg_int(x)
  
  # R: illegal_trans
  loss <- loss + o_om_illegal_trans(x)
  
  if (!isGrad) print(loss)
  loss
}


o_om_grad <- function(x) {
  grad <- rep(0, length(x))
  
  for (idx in seq_len(length.out = length(x))) {
    # First we check the conditions:
    tixi <- x[idx] + q[idx, 1]
    
    cond_pre <- if (idx == 1) TRUE else {
      tjxj <- x[1:(idx - 1)] + q[1:(idx - 1), 1]
      all(tjxj < tixi)
    }
    if (cond_pre) {
      # grad[idx] <- -2 * (q[idx, 2] - r(tixi)) * pracma::fderiv(f = r, x = tixi, method = "central")
      grad_denom <- -2 * q[idx, 2] * r(tixi) + r(tixi)^2 + q[idx, 2]^2 + 1
      if (grad_denom == 0) {
        grad_denom <- .Machine$double.eps
      }
      grad_num <- 2 * pracma::fderiv(f = r, x = tixi, method = "central") * (r(tixi) - q[idx, 2])
      grad[idx] <- grad_num / grad_denom
    }
  }
  
  grad_r1r2 <- pracma::grad(f = o_om_r1r2, x0 = x)
  grad <- grad + grad_r1r2
  
  # # neg_int:
  # grad_neg_int <- pracma::grad(f = o_om_neg_int, x0 = x)
  # grad <- grad + grad_neg_int
  
  # illegal_trans:
  grad_illegal_trans <- pracma::grad(f = o_om_illegal_trans, x0 = x)
  grad <- grad + grad_illegal_trans
  
  grad
}






