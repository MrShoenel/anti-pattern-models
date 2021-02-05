library(R6)



SRBTW <- R6Class(
  "SRBTW",
  
  lock_objects = FALSE,
  
  private = list(
    openBegin = NULL,
    openEnd = NULL,
    wp = NULL,
    wr = NULL,
    Q = NULL,
    subModels = NULL,
    gamma_bed = NULL,
    lambda = NULL,
    begin = NULL,
    end = NULL,
    thetaB = NULL,
    varthetaL = NULL,
    
    checkQ = function(q) {
      stopifnot(is.numeric(q) && !is.na(q) && q >= min(private$Q) && q <= max(private$Q))
    }
  ),
  
  public = list(
    initialize = function(wp, wc, theta_b, gamma_bed, lambda, begin, end, openBegin = FALSE, openEnd = FALSE) {
      stopifnot(is.function(wp) && is.function(wc))
      
      # Check theta_b:
      stopifnot(is.vector(theta_b) && is.numeric(theta_b) && !any(is.na(theta_b)) && length(theta_b) >= 2)
      stopifnot(all.equal(theta_b, sort(theta_b)) && length(unique(theta_b)) == length(theta_b))
      private$thetaB <- theta_b
      
      # Check gamma_bed (theta_b, theta_e, theta_d):
      stopifnot(is.vector(gamma_bed) && is.numeric(gamma_bed) && !any(is.na(gamma_bed)) && length(gamma_bed) == 3)
      stopifnot(gamma_bed[3] >= 0)
      stopifnot((gamma_bed[1] + gamma_bed[3]) <= gamma_bed[2])
      private$gamma_bed <- gamma_bed
      
      # Check lambda:
      stopifnot(is.vector(lambda) && is.numeric(lambda) && !any(is.na(lambda)) && all(lambda >= 0) && length(lambda) == (length(theta_b) - 1))
      private$lambda <- lambda
      
      # Check b,e: During construction, these need to be valid.
      # During optimization, they may be invalid.
      stopifnot(is.numeric(c(begin, end)) && !any(is.na(c(begin, end))))
      stopifnot(begin < end)
      private$begin <- begin
      private$end <- end
      
      private$wp <- wp
      private$wc <- wc
      
      private$Q <- seq_len(length.out = length(lambda))
      private$subModels <- list()
      
      self$setOpenBegin(ob = openBegin)
      self$setOpenEnd(oe = openEnd)
    },
    
    getQ = function() {
      private$Q
    },
    
    getLambda = function() {
      private$lambda
    },
    
    getLambda_q = function(q) {
      private$checkQ(q)
      private$lambda[q]
    },
    
    getWP = function() {
      private$wp
    },
    
    getWC = function() {
      private$wc
    },
    
    setParams = function(vartheta_l = NULL, begin = NULL, end = NULL) {
      if (!missing(vartheta_l)) {
        stopifnot(is.vector(vartheta_l) && is.numeric(vartheta_l) && !any(is.na(vartheta_l)) && length(vartheta_l) == (length(private$thetaB) - 1))
        private$varthetaL <- vartheta_l
      }
      
      if (!missing(begin)) {
        stopifnot(is.numeric(begin) && !is.na(begin))
        private$begin <- begin
      }
      
      if (!missing(end)) {
        stopifnot(is.numeric(end) && !is.na(end))
        private$end <- end
      }
      
      invisible(self)
    },
    
    getBegin = function() {
      stopifnot(!is.null(private$begin))
      private$begin
    },
    
    getEnd = function() {
      stopifnot(!is.null(private$end))
      private$end
    },
    
    getTb_q = function(q) {
      private$checkQ(q)
      private$thetaB[q]
    },
    
    getTe_q = function(q) {
      private$checkQ(q)
      private$thetaB[q + 1]
    },
    
    getLength_q = function(q) {
      private$checkQ(q)
      private$varthetaL[q]
    },
    
    getVarthetaL = function() {
      c(private$varthetaL)
    },
    
    setOpenBegin = function(ob = TRUE) {
      stopifnot(is.logical(ob))
      private$openBegin <- ob
      invisible(self)
    },
    
    isOpenBegin = function() {
      private$openBegin
    },
    
    setOpenEnd = function(oe = TRUE) {
      stopifnot(is.logical(oe))
      private$openEnd <- oe
      invisible(self)
    },
    
    isOpenEnd = function() {
      private$openEnd
    },
    
    getgamma_bed = function() {
      private$gamma_bed
    },
    
    getBeta_l = function() {
      gamma_bed <- self$getgamma_bed()
      min(gamma_bed[2] - gamma_bed[3],
          max(gamma_bed[1], min(private$begin, private$end)))
    },
    
    getBeta_u = function() {
      gamma_bed <- self$getgamma_bed()
      max(gamma_bed[1] + gamma_bed[3],
          min(gamma_bed[2], max(private$begin, private$end)))
    },
    
    getQForX = function(x) {
      b_l <- self$getBeta_l()
      b_u <- self$getBeta_u()
      stopifnot(x >= b_l && x <= b_u)
      
      for (q in self$getQ()) {
        if (x >= private$thetaB[q] && x < private$thetaB[q + 1]) {
          return(q)
        }
      }
      
      # x must be == b_u:
      return(max(self$getQ()))
    },
    
    getPhi = function() {
      fac1 <- self$getBeta_u() - self$getBeta_l()
      gamma_bed <- self$getgamma_bed()
      theta_d <- gamma_bed[3]
      
      lambda <- self$getLambda()
      fac2 <- max(theta_d, sum(sapply(X = self$getQ(), FUN = function(q) {
        max(lambda[q], self$getLength_q(q))
      })))
    },
    
    getPhi_q = function(q) {
      private$checkQ(q)
      
      if (q == 1) {
        return(0)
      }
      
      # Else: q > 1
      lambda <- self$getLambda()
      X <- seq_len(length.out = q - 1)
      sum(sapply(X = X, FUN = function(q) {
        max(lambda[q], self$getLength_q(q))
      })) / (self$getBeta_u() - self$getBeta_l())
    },
    
    
    getSubModel = function(q) {
      private$checkQ(q)
      
      qs <- paste0(q)
      if (!(qs %in% names(private$subModels))) {
        # Lazy init:
        private$subModels[[qs]] <- SRBTW_SubModel$new(srbtw = self, q = q)
      }
      
      private$subModels[[qs]]
    }
  )
)


SRBTW_SubModel <- R6Class(
  "SRBTW_SubModel",
  
  private = list(
    srbtw = NULL,
    q = NULL
  ),
  
  public = list(
    initialize = function(srbtw, q) {
      stopifnot(R6::is.R6(srbtw))
      stopifnot(is.numeric(q) && !is.na(q) && q >= 1)
      
      private$srbtw <- srbtw
      private$q <- q
    },
    
    getQ = function() {
      private$q
    },
    
    asTuple = function() {
      q <- private$q
      s <- private$srbtw
      vtl <- s$getVarthetaL()
      
      beta_l <- s$getBeta_l()
      beta_u <- s$getBeta_u()
      lambda <- s$getLambda()
      l_q <- s$getLength_q(q)
      l_prime_q <- max(lambda[q], max(-l_q, l_q))
      psi <- sum(sapply(s$getQ(), function(i) {
        max(lambda[i], max(-vtl[i], vtl[i]))
      }))
      l_q_c <- l_prime_q / psi * (beta_u - beta_l)
      
      phi_q <- if (q == 1) 0 else sum(sapply(seq_len(length.out = q - 1), function(i) {
        max(lambda[i], max(-vtl[i], vtl[i])) / psi * (beta_u - beta_l)
      }))
      sb_q <- beta_l + phi_q
      se_q <- sb_q + l_q_c
      tb_q <- s$getTb_q(q)
      te_q <- s$getTe_q(q)
      delta_t_q <- te_q - tb_q
      f <- s$getWC()
      
      list(
        q = q,
        beta_l = beta_l,
        beta_u = beta_u,
        lambda_q = lambda[q],
        l_q = l_q,
        l_prime_q = l_prime_q,
        l_q_c = l_q_c,
        psi = psi,
        phi_q = phi_q,
        delta_t_q = delta_t_q,
        sb_q = sb_q,
        se_q = se_q,
        tb_q = tb_q,
        te_q = te_q,
        mqc = function(x) {
          x <- (x - tb_q) * l_q_c / delta_t_q + sb_q
          x <- max(beta_l, min(beta_u, x))
          f(x)
        }
      )
    }
  )
)



# SRBTW_SubModel$debug("asTuple")
# SRBTW$debug("initialize")


