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
      stopifnot(is.vector(params) && is.numeric(params) && !any(is.na(params)) && length(params) == (length(private$thetaB) - 1 + (if (self$isOpenBegin()) 1 else 0) + (if (self$isOpenEnd()) 1 else 0)))
      
      ob <- self$isOpenBegin()
      oe <- self$isOpenEnd()
      begin <- if (ob) params[length(private$thetaB)]
      end <- if (oe) params[length(private$thetaB) + (if (ob) 1 else 0)]
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
      r <- range(private$thetaB)
      stopifnot(x >= r[1] && x <= r[2])
      
      for (q in self$getQ()) {
        if (x >= private$thetaB[q] && x < private$thetaB[q + 1]) {
          return(q)
        }
      }
      
      # x must be == b_u:
      return(max(self$getQ()))
    },
    
    getPsi = function() {
      lambda <- self$getLambda()
      vtl <- self$getVarthetaL()
      
      sum(sapply(self$getQ(), function(i) {
        max(lambda[i], max(-vtl[i], vtl[i]))
      }))
    },
    
    getPhi_q = function(q) {
      private$checkQ(q)
      lambda <- self$getLambda()
      vtl <- self$getVarthetaL()
      psi <- self$getPsi()
      beta_l <- self$getBeta_l()
      beta_u <- self$getBeta_u()
      
      if (q == 1) 0 else sum(sapply(seq_len(length.out = q - 1), function(i) {
        max(lambda[i], max(-vtl[i], vtl[i])) / psi * (beta_u - beta_l)
      }))
    },
    
    getSubModel = function(q) {
      private$checkQ(q)
      
      qs <- paste0(q)
      if (!(qs %in% names(private$subModels))) {
        # Lazy init:
        private$subModels[[qs]] <- SRBTW_SubModel$new(srbtw = self, q = q)
      }
      
      private$subModels[[qs]]
    },
    
    M = function(x) {
      q <- self$getQForX(x)
      sm <- self$getSubModel(q = q)
      t <- sm$asTuple()
      t$mqc(x)
    },
    
    plot_warp = function(WP = TRUE, WC = TRUE, M = TRUE) {
      funcsColors <- c("red", "black", "blue")
      funcs <- c()
      p <- ggplot2::ggplot()
      if (WP) {
        funcs <- c(funcs, "WP")
        p <- p + ggplot2::stat_function(
          fun = private$wp, mapping = ggplot2::aes(color = "WP"))
      }
      if (WC) {
        funcs <- c(funcs, "WC")
        p <- p + ggplot2::stat_function(
          fun = private$wc, mapping = ggplot2::aes(color = "WC"))
      }
      if (M) {
        funcs <- c(funcs, "M")
        p <- p + ggplot2::stat_function(
          fun = Vectorize(function(x) {
            self$M(x = x)
          }), mapping = ggplot2::aes(color = "M"))
      }
      
      p + ggplot2::scale_color_manual(
        paste(funcs, collapse = " / "),
        values = head(x = funcsColors, n = length(funcs))) +
        ggplot2::xlim(range(
          c(range(private$thetaB), self$getBeta_l(), self$getBeta_u())
        )) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank())
    },
    
    #' TODO: Make pretty following the example at
    #' https://stackoverflow.com/questions/17148679/construct-a-manual-legend-for-a-complicated-plot
    plot_dual_original = function(vOffset = 1, showGammaB = TRUE, showGammaE = TRUE, showBetaL = TRUE, showBetaU = TRUE) {
      plotFac <- factor(levels = c("WP", "WC", "M", "Beta_L", "Beta_U", "Gamma_B", "Gamma_E"), ordered = TRUE)
      getFac <- function(f) factor(x = f, levels = levels(plotFac), ordered = TRUE)
      
      p <- ggplot2::ggplot()
      
      wp <- Vectorize(function(x) {
        private$wp(x) + vOffset
      })
      
      wc <- Vectorize(function(x) {
        private$wc(x)
      })
      
      tb <- private$thetaB
      beta_l <- self$getBeta_l()
      beta_u <- self$getBeta_u()
      gamma_bed <- self$getgamma_bed()
      
      extWp <- range(tb)
      extWc <- range(
        c(beta_l, beta_u,
          (if (showGammaB) gamma_bed[1] else c()),
          (if (showGammaE) gamma_bed[2] else c())))
      
      # Add the two signals:
      p <- p + ggplot2::stat_function(
        fun = wp, mapping = ggplot2::aes(color = getFac("WP")), xlim = extWp)
      p <- p + ggplot2::stat_function(
        fun = wc, mapping = ggplot2::aes(color = getFac("WC")), xlim = extWc)
      
      
      tb_wc <- sapply(self$getQ(), function(q) beta_l + self$getPhi_q(q = q))
      tb_wc <- c(tb_wc, beta_u)
      
      y_wp <- wp(tb)
      y_wc <- wc(tb_wc)
      
      
      # Now we can add the lines:
      for (i in seq_len(length.out = length(tb))) {
        p <- p + geom_path(
          data = data.frame(x = c(tb[i], tb_wc[i]), y = c(y_wp[i], y_wc[i])),
          mapping = aes(x = x, y = y), color = "blue", linetype = "dashed", size = .25)
      }
      
      ob <- self$isOpenBegin()
      oe <- self$isOpenEnd()
      
      # Conditionally add lines for open begin and/or end:
      dfBeta <- NULL
      if (showBetaL && ob) {
        dfBeta <- rbind(dfBeta, data.frame(x = beta_l, type = getFac("Beta_L")))
      }
      if (showBetaU && oe) {
        dfBeta <- rbind(dfBeta, data.frame(x = beta_u, type = getFac("Beta_U")))
      }
      if (is.data.frame(dfBeta)) {
        p <- p + ggplot2::geom_vline(
          mapping = ggplot2::aes(xintercept = x, color = type),
          data = dfBeta, show.legend = TRUE)
      }
      
      # Conditionally add the polygons for when gamma_b < b or gamma_e > e:
      dfGamma <- NULL
      gammaNumX <- 200
      if (showGammaB && ob) {
        gammaB_x <- seq(gamma_bed[1], beta_l, length.out = gammaNumX)
        dfGamma <- rbind(dfGamma, data.frame(
          x = c(min(tb), gammaB_x),
          y = c(wp(min(tb)), sapply(X = gammaB_x, FUN = wc)),
          type = rep(getFac("Gamma_B"), 1 + gammaNumX)
        ))
      }
      if (showGammaE && oe) {
        gammaE_x <- seq(gamma_bed[2], beta_u, length.out = gammaNumX)
        dfGamma <- rbind(dfGamma, data.frame(
          x = c(max(tb), gammaE_x),
          y = c(wp(max(tb)), sapply(X = gammaE_x, FUN = wc)),
          type = rep(getFac("Gamma_E"), 1 + gammaNumX)
        ))
      }
      if (is.data.frame(dfGamma)) {
        p <- p + ggplot2::geom_polygon(
          mapping = ggplot2::aes(x = x, y = y, fill = type),
          data = dfGamma, alpha = .25, show.legend = TRUE)
      }
      
      p +
        scale_color_brewer(palette = "Set1")  +
        scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(color = "Variable", fill = "Unused") +
        ggplot2::xlim(range(c(extWp, extWc))) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank())
    }
  )
)


#' This model is SRBTW, wrapped in SRBAW, so we can do both optimizations.
SRBTWBAW <- R6Class(
  "SRBTWBAW",
  
  inherit = SRBTW,
  
  private = list(
    v = NULL,
    varthetaY = NULL,
    lambdaYmin = NULL,
    lambdaYmax = NULL
  ),
  
  public = list(
    initialize = function(
      wp, wc, theta_b, gamma_bed, lambda, begin, end, # SRBTW
      lambda_ymin, lambda_ymax, # SRBAW (!)
      openBegin = FALSE, openEnd = FALSE, wc_d1 = NULL, wc_d2 = NULL # SRBTW
    ) {
      super$initialize(
        wp = wp, wc = wc, theta_b = theta_b, gamma_bed = gamma_bed, lambda = lambda,
        begin = begin, end = end, openBegin = openBegin, openEnd = openEnd,
        wc_d1 = wc_d1, wc_d2 = wc_d2)
      
      stopifnot(is.numeric(c(lambda_ymin, lambda_ymax)) && !any(is.na(c(lambda_ymin, lambda_ymax))))
      stopifnot(length(lambda_ymin) == length(lambda_ymax) && length(lambda_ymin) == (length(theta_b) - 1))
      
      private$lambdaYmin <- lambda_ymin
      private$lambdaYmax <- lambda_ymax
    },
    
    #' For SRBTWBAW, params is vartheta_l, [,b [, e]], v, vartheta_y
    setAllParams = function(params) {
      stopifnot(is.vector(params) && is.numeric(params) && !any(is.na(params)))
      
      ob <- self$isOpenBegin()
      oe <- self$isOpenEnd()
      # The total length is vartheta_l + [+b] [+e] v + vartheta_y,
      # and the two varthetas have the same length.
      vtLen <- length(private$thetaB) - 1
      oboeLen <- (if (ob) 1 else 0) + (if (oe) 1 else 0)
      superLen <- vtLen + oboeLen
      requireLen <- vtLen + oboeLen + 1 + vtLen
      stopifnot(length(params) == requireLen)
      
      super$setAllParams(
        params = params[seq_len(length.out = superLen)])
      
      v <- params[superLen + 1]
      vartheta_y <- params[(superLen + 2):length(params)]
      
      self$setParams(v = v, vartheta_y = vartheta_y)
    },
    
    setParams = function(vartheta_l = NULL, begin = NULL, end = NULL, v = NULL, vartheta_y = NULL) {
      super$setParams(vartheta_l = vartheta_l, begin = begin, end = end)
      
      if (!missing(v) && !is.null(v)) {
        stopifnot(is.numeric(v) && !is.na(v))
        private$v <- v
      }
      if (!missing(vartheta_y) && !is.null(vartheta_y)) {
        stopifnot(is.numeric(vartheta_y) && !any(is.na(vartheta_y)))
        stopifnot(length(vartheta_y) == (length(private$thetaB) - 1))
        private$varthetaY <- vartheta_y
      }
      
      invisible(self)
    },
    
    getVarthetaY_q = function(q) {
      private$checkQ(q = q)
      private$varthetaY[q]
    },
    
    getV = function() {
      private$v
    },
    
    getPhi_y_q = function(q) {
      private$checkQ(q = q)
      
      if (q == 1) {
        return(0)
      }
      
      # Else: q > 1
      X <- seq_len(length.out = q - 1)
      sum(sapply(X = X, FUN = function(r) {
        self$getVarthetaY_q(q = r)
      }))
    },
    
    getLambdaYmin = function() {
      private$lambdaYmin
    },
    
    getLambdaYmax = function() {
      private$lambdaYmax
    },
    
    getSubModel = function(q) {
      private$checkQ(q)
      
      qs <- paste0(q)
      if (!(qs %in% names(private$subModels))) {
        # Lazy init:
        private$subModels[[qs]] <- SRBTWBAW_SubModel$new(srbtw = self, q = q)
      }
      
      private$subModels[[qs]]
    },
    
    #' Overridden because some funcs, such as plot_warp, call M()
    M = function(x) {
      q <- self$getQForX(x)
      sm <- self$getSubModel(q = q)
      t <- sm$asTuple()
      t$nqc(x)
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
      beta_u <- s$getBeta_u()
      lambda <- s$getLambda()
      l_q <- s$getLength_q(q)
      l_prime_q <- max(lambda[q], max(-l_q, l_q))
      psi <- s$getPsi()
      l_q_c <- l_prime_q / psi * (beta_u - beta_l)
      phi_q <- s$getPhi_q(q)
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
        wc = wc,
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




SRBTWBAW_SubModel <- R6Class(
  "SRBTWBAW_SubModel",
  
  inherit = SRBTW_SubModel,
  
  public = list(
    initialize = function(srbtw, q) {
      super$initialize(srbtw = srbtw, q = q)
    },
    
    asTuple = function() {
      q <- private$q
      s <- private$srbtw
      t <- super$asTuple()
      
      lambda_ymin <- s$getLambdaYmin()
      lambda_ymax <- s$getLambdaYmax()
      
      v <- s$getV()
      vty_q <- s$getVarthetaY_q(q = q)
      phi_y_q <- s$getPhi_y_q(q = q)
      a_q <- vty_q / t$l_q_c
      lambda_ymin_q <- lambda_ymin[q]
      lambda_ymax_q <- lambda_ymax[q]
      
      t$iota_q <- t$l_q_c
      t$v <- v
      t$vty_q <- vty_q
      t$a_q <- a_q
      t$phi_y_q <- phi_y_q
      t$lambda_ymin_q <- lambda_ymin_q
      t$lambda_ymax_q <- lambda_ymax_q
      
      tq <- function(x) {
        a_q * (x - t$phi_q) + v + phi_y_q
      }
      t$tq <- tq
      
      t$nqc <- function(x) {
        nq <- t$mqc(x = x) + tq(x = x)
        max(lambda_ymin_q, min(nq, lambda_ymax_q))
      }
      
      t$nqc_d1 <- function(x) {
        # Gradient of nqc, needs wc, wc_d1
        stop("TODO")
      }
      t$nqc_d2 <- function(x) {
        # 2nd-order Gradient of nqc, needs wc, wc_d1, wc_d2
        stop("TODO")
      }
      
      t
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
    
    isContinuous = function() {
      private$continuous
    },
    
    setParams = function(params) {
      private$params <- params
      private$srbtw$setAllParams(params = params)
      invisible(self)
    },
    
    getParams = function() {
      stopifnot(!is.null(private$params))
      private$params
    },
    
    setWeight = function(weight) {
      stopifnot(is.numeric(weight) && !is.na(weight) && weight > 0)
      private$weight <- weight
      invisible(self)
    },
    
    getWeight = function() {
      private$weight
    },
    
    compute = function() {
      stopifnot(length(private$intervals) > 0)
      
      listOfSms <- list()
      for (q in private$intervals) {
        listOfSms <- append(listOfSms, private$srbtw$getSubModel(q = q))
      }
      
      listOfSms
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


SRBTW_SingleObjectiveOptimization <- R6Class(
  "SRBTW_SingleObjectiveOptimization",
  
  lock_objects = FALSE,
  
  inherit = SRBTW_Loss,
  
  private = list(
    srbtw = NULL,
    objectives = NULL
  ),
  
  public = list(
    initialize = function(srbtw) {
      stopifnot(R6::is.R6(srbtw))
      private$srbtw <- srbtw
      
      private$objectives <- list()
    },
    
    setParams = function(params) {
      private$params <- params
      private$srbtw$setAllParams(params = params)
      
      for (obj in private$objectives) {
        obj$setParams(params = params)
      }
      
      invisible(self)
    },
    
    addObjective = function(obj) {
      stopifnot(R6::is.R6(obj) && inherits(obj, "SRBTW_Loss"))
      
      private$objectives <- append(private$objectives, obj)
      invisible(self)
    },
    
    compute = function() {
      loss <- 0
      
      for (obj in private$objectives) {
        loss <- loss + obj$getWeight() * obj$compute()
      }
      
      loss
    },
    
    computeGrad = function() {
      grad <- NULL
      
      for (obj in private$objectives) {
        temp <- obj$getWeight() * obj$computeGrad()
        grad <- if (is.null(grad)) temp else grad + temp
      }
      
      grad
    },
    
    computeHess = function() {
      # I need to check how to do this exactly before we can proceed.
      stop("TODO")
    }
  )
)


# SRBTW$debug("setAllParams")
# SRBTW_SubModel$debug("asTuple")
# SRBTW$debug("initialize")
# SRBTW_Loss$debug("compute")
# SRBTW_DataLoss$debug("compute")
# SRBTW_SingleObjectiveOptimization$debug("compute")


