library(R6)


SRBTW <- R6Class(
  "SRBTW",
  
  lock_objects = TRUE,
  
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
    varthetaL = NULL
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
      stopifnot(is.vector(lambda) && is.numeric(lambda) && !any(is.na(lambda)) && length(lambda) == (length(theta_b) - 2))
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
      stopifnot(is.numeric(q) && !is.na(q))
      Q <- self$getQ()
      stopifnot(q >= 1 && q <= Q)
      private$lambda[q]
    },
    
    getWP = function() {
      private$wp
    },
    
    getWC = function() {
      private$wc
    },
    
    getDelta_q = function(q) {
      stopifnot(is.numeric(q) && !is.na(q))
      Q <- self$getQ()
      stopifnot(q >= 1 && q <= Q)
      private$thetaB[q + 1] - private$thetaB[q]
    },
    
    setParams = function(vartheta_l = NULL, begin = NULL, end = NULL) {
      if (!missing(vartheta_l)) {
        stopifnot(is.vector(vartheta_l) && is.numeric(vartheta_l) && !any(is.na(vartheta_l)) && length(vartheta_l) == self$getQ())
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
    
    getSb_q = function(q) {
      stopifnot(!is.null(private$varthetaL))
      stopifnot(is.numeric(q) && !is.na(q))
      Q <- self$getQ()
      stopifnot(q >= 1 && q <= Q)
      
      private$thetaB[q]
    },
    
    getLength_q = function(q) {
      stopifnot(!is.null(private$varthetaL))
      stopifnot(is.numeric(q) && !is.na(q))
      Q <- self$getQ()
      stopifnot(q >= 1 && q <= Q)
      
      private$varthetaL[q]
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
          max(gamma_bed[1], min(begin, end)))
    },
    
    getBeta_u = function() {
      gamma_bed <- self$getgamma_bed()
      max(gamma_bed[1] + gamma_bed[3],
          min(gamma_bed[2], max(begin, end)))
    },
    
    getPhi = function() {
      fac1 <- self$getBeta_u() - self$getBeta_l()
      gamma_bed <- self$getgamma_bed()
      theta_d <- gamma_bed[3]
      
      lambda <- self$getLambda()
      X <- seq_len(length.out = self$getQ())
      fac2 <- max(theta_d, sum(sapply(X = X, FUN = function(q) {
        max(lambda[q], self$getLength_q(q))
      })))
    },
    
    getPhi_q = function(q) {
      stopifnot(!is.null(private$varthetaL))
      stopifnot(is.numeric(q) && !is.na(q))
      Q <- self$getQ()
      stopifnot(q >= 1 && q <= Q)
      
      if (q == 1) {
        return(0)
      }
      
      # Else: q > 1
      lambda <- self$getLambda()
      X <- seq_len(length.out = q - 1)
      sum(sapply(X = X, FUN = function(q) {
        max(lambda[q], self$getLength_q(q))
      }))
    },
    
    
    getSubModel = function(q) {
      stopifnot(is.numeric(q) && !is.na(q))
      Q <- self$getQ()
      stopifnot(q >= 1 && q <= Q)
      
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
      
      beta_l <- s$getBeta_l()
      phi_q <- s$getPhi_q(q)
      phi <- s$getPhi()
      delta_q <- s$getDelta_q(q)
      l_q <- s$getLength_q(q)
      lambda_q <- s$getLambda_q(q)
      sb_q <- s$getSb_q(q)
      f <- s$getWC()
      
      list(
        q = q,
        beta_l = beta_l,
        phi_q = phi_q,
        phi = phi,
        delta_q = delta_q,
        l_q = l_q,
        sb_q = sb_q,
        se_q = sb_q + delta_q,
        mqc = function(x) {
          f((x - beta_l - phi_q) * phi * delta_q / max(lambda_q, l_q) + sb_q)
        }
      )
    }
  )
)



