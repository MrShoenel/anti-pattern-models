library(R6)



SRBTW <- R6Class(
  "SRBTW",
  
  lock_objects = FALSE,
  
  private = list(
    openBegin = NULL,
    openEnd = NULL,
    wp = NULL,
    wc = NULL,
    wc_d1 = NULL,
    wc_d2 = NULL,
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
    initialize = function(wp, wc, theta_b, gamma_bed, lambda, begin, end, openBegin = FALSE, openEnd = FALSE, wc_d1 = NULL, wc_d2 = NULL) {
      stopifnot(is.function(wp) && is.function(wc))
      stopifnot(missing(wc_d1) || is.null(wc_d1) || is.function(wc_d1))
      stopifnot(missing(wc_d2) || is.null(wc_d2) || is.function(wc_d2))
      
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
      private$wc_d1 <- wc_d1
      private$wc_d2 <- wc_d2
      
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
    
    getWC_d1 = function() {
      private$wc_d1
    },
    
    getWC_d2 = function() {
      private$wc_d2
    },
    
    setAllParams = function(params) {
      stopifnot(is.vector(params) && is.numeric(params) && !any(is.na(params)) && length(params) == (length(private$thetaB) - 1 + (if (self$isOpenBegin()) 1 else 0) + if (self$isOpenEnd()) 1 else 0))
      
      begin <- if (self$isOpenBegin()) params[length(private$thetaB)]
      end <- if (self$isOpenEnd()) params[length(private$thetaB) + 1]
      self$setParams(
        vartheta_l = params[seq_len(length.out = length(private$thetaB) - 1)],
        begin = begin, end = end)
    },
    
    setParams = function(vartheta_l = NULL, begin = NULL, end = NULL) {
      if (!missing(vartheta_l) && !is.null(vartheta_l)) {
        stopifnot(is.vector(vartheta_l) && is.numeric(vartheta_l) && !any(is.na(vartheta_l)) && length(vartheta_l) == (length(private$thetaB) - 1))
        private$varthetaL <- vartheta_l
      }
      
      if (!missing(begin) && !is.null(begin)) {
        stopifnot(is.numeric(begin) && !is.na(begin))
        private$begin <- begin
      }
      
      if (!missing(end) && !is.null(end)) {
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
      wp <- s$getWP()
      wc <- s$getWC()
      wc_d1 <- s$getWC_d1()
      wc_d2 <- s$getWC_d2()
      
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
        wp = wp,
        wc = wp,
        wc_d1 = wc_d1,
        wc_d2 = wc_d2,
        mqc = function(x) {
          x <- (x - tb_q) * l_q_c / delta_t_q + sb_q
          # For numeric stability, rounding sometimes makes
          # x slightly larger or smaller than the boundaries.
          x <- max(beta_l, min(beta_u, x))
          wc(x)
        },
        mqc_d1 = function(x) {
          # Gradient of mqc, needs wc, wc_d1
          stop("TODO")
        },
        mqc_d2 = function(x) {
          # 2nd-order Gradient of mqc, needs wc, wc_d1, wc_d2
          stop("TODO")
        }
      )
    }
  )
)


SRBTW_Loss <- R6Class(
  "SRBTW_Loss",
  
  lock_objects = FALSE,
  
  private = list(
    srbtw = NULL,
    params = NULL,
    intervals = NULL,
    weight = NULL,
    continuous = NULL,
    numSamples = NULL
  ),
  
  public = list(
    initialize = function(
      srbtw, intervals = c(), params = NULL, weight = 1, continuous = FALSE,
      numSamples = if (continuous) NA_integer_ else length(intervals) * 1e3)
    {
      stopifnot(R6::is.R6(srbtw))
      stopifnot(is.logical(continuous))
      stopifnot(continuous == TRUE || (is.numeric(numSamples) && !is.na(numSamples) && numSamples > 0))
      
      private$srbtw <- srbtw
      
      if (!missing(intervals)) {
        private$intervals <- intervals
      }
      if (!missing(params)) {
        self$setParams(params = params)
      }
      
      self$setWeight(weight = weight)
      private$continuous <- continuous
      private$numSamples <- numSamples
    },
    
    setParams = function(params) {
      private$params <- params
      private$srbtw$setAllParams(params = params)
      invisible(self)
    },
    
    getParams = function() {
      stopifnot(is.null(private$params))
      private$params
    },
    
    setWeight = function(weight) {
      stopifnot(is.numeric(weight) && !is.na(weight) && weight > 0)
      private$weight <- weight
      invisible(self)
    },
    
    compute = function() {
      stopifnot(length(private$intervals) > 0)
      
      lapply(private$intervals, function(q) {
        private$srbtw$getSubModel(q = q)$asTuple()
      })
    },
    
    computeGrad = function() {
      stop("Abstract method")
    },
    
    computeGradNumeric = function() {
      pracma::grad(f = function(x) {
        self$setParams(params = x)
        self$compute()
      }, x0 = self$getParams())
    },
    
    computeHess = function() {
      stop("Abstract method")
    },
    
    computeHessNumeric = function() {
      pracma::hessian(f = function(x) {
        self$setParams(params = x)
        self$compute()
      }, x0 = self$getParams())
    }
  )
)


SRBTW_DataLoss <- R6Class(
  "SRBTW_DataLoss",
  
  inherit = SRBTW_Loss,
  
  private = list(
    lossFunc = NULL,
    lossFuncGrad = NULL,
    lossFuncHess = NULL
  ),
  
  public = list(
    initialize = function(
      srbtw, intervals = c(), params = NULL, weight = 1, continuous = FALSE,
      numSamples = if (continuous) NA_integer_ else length(intervals) * 1e3,
      lossFunc = NULL, lossFuncGrad = NULL, lossFuncHess = NULL
    ) {
      super$initialize(
        srbtw = srbtw, intervals = intervals, params = params,
        weight = weight, continuous = continuous, numSamples = numSamples)
      
      self$setLossFunc(lossFunc = lossFunc)
      self$setLossFuncGrad(lossFuncGrad = lossFuncGrad)
      self$setLossFuncHess(lossFuncHess = lossFuncHess)
    },
    
    setLossFunc = function(lossFunc = NULL) {
      stopifnot(missing(lossFunc) || is.null(lossFunc) || is.function(lossFunc))
      private$lossFunc <- lossFunc
      invisible(self)
    },
    
    setLossFuncGrad = function(lossFuncGrad = NULL) {
      stopifnot(missing(lossFuncGrad) || is.null(lossFuncGrad) || is.function(lossFuncGrad))
      private$lossFuncGrad <- lossFuncGrad
      invisible(self)
    },
    
    setLossFuncHess = function(lossFuncHess = NULL) {
      stopifnot(missing(lossFuncHess) || is.null(lossFuncHess) || is.function(lossFuncHess))
      private$lossFuncHess <- lossFuncHess
      invisible(self)
    },
    
    getLossFunc = function() {
      stopifnot(is.function(private$lossFunc))
      private$lossFunc
    },
    
    compute = function() {
      lossFn <- self$getLossFunc()
      
      loss <- 0
      listOfSms <- super$compute()
      
      loss <- loss + do.call(
        what = lossFn, args = list(loss = self, listOfSms = listOfSms))
      
      loss
    },
    
    computeGrad = function() {
      res <- NULL
      if (is.function(private$lossFuncGrad)) {
        # Each sm has ::asTuple(), which includes m_q^c and its derivative m_q^c_d1
        listOfSms <- super$compute()
        res <- do.call(
          what = private$lossFuncGrad, args = list(loss = self, listOfSms = listOfSms))
      } else {
        res <- self$computeGradNumeric()
      }
      
      res
    },
    
    computeHess = function() {
      res <- NULL
      if (is.function(private$lossFuncHess)) {
        # Each sm has ::asTuple(), which includes m_q^c and its derivatives m_q^c_d1, m_q^c_d2
        listOfSms <- super$compute()
        res <- do.call(
          what = private$lossFuncHess, args = list(loss = self, listOfSms = listOfSms))
      } else {
        res <- self$computeHessNumeric()
      }
      
      res
    }
  )
)



# SRBTW_SubModel$debug("asTuple")
# SRBTW$debug("initialize")


