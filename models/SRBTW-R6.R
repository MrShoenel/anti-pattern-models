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
    thetaB = NULL,
    begin = NULL,
    end = NULL,
    varthetaL = NULL,
    
    checkQ = function(q) {
      stopifnot(is.numeric(q) && !is.na(q) && q >= min(private$Q) && q <= max(private$Q))
    },
    
    requireParams = function() {
      stopifnot(!is.null(private$varthetaL))
      stopifnot(!is.null(private$begin))
      stopifnot(!is.null(private$end))
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
      private$requireParams()
      private$begin
    },
    
    getEnd = function() {
      private$requireParams()
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
      private$requireParams()
      private$checkQ(q)
      private$varthetaL[q]
    },
    
    getVarthetaL = function() {
      private$requireParams()
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
    lambdaYmax = NULL,
    
    requireParamsBaw = function() {
      private$requireParams()
      stopifnot(!is.null(private$v))
      stopifnot(!is.null(private$varthetaY))
    }
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
      private$requireParamsBaw()
      private$checkQ(q = q)
      private$varthetaY[q]
    },
    
    getV = function() {
      private$requireParamsBaw()
      private$v
    },
    
    getPhi_y_q = function(q) {
      private$requireParamsBaw()
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
    
    getLambdaYmin_q = function(q) {
      private$checkQ(q)
      private$lambdaYmin[q]
    },
    
    getLambdaYmax = function() {
      private$lambdaYmax
    },
    
    getLambdaYmax_q = function(q) {
      private$checkQ(q)
      private$lambdaYmax[q]
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
      
      v <- s$getV()
      vty_q <- s$getVarthetaY_q(q = q)
      phi_y_q <- s$getPhi_y_q(q = q)
      # a_q <- vty_q / t$l_q_c
      tb_q <- s$getTb_q(q = q)
      a_q <- vty_q / (s$getTe_q(q = q) - tb_q)
      lambda_ymin_q <- s$getLambdaYmin_q(q = q)
      lambda_ymax_q <- s$getLambdaYmax_q(q = q)
      
      t$iota_q <- t$l_q_c
      t$v <- v
      t$vty_q <- vty_q
      t$a_q <- a_q
      t$phi_y_q <- phi_y_q
      t$lambda_ymin_q <- lambda_ymin_q
      t$lambda_ymax_q <- lambda_ymax_q
      
      tq <- function(x) {
        a_q * (x - tb_q) + v + phi_y_q
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
      stop("Abstract method 1")
    },
    
    computeGradNumeric = function() {
      pracma::grad(f = function(x) {
        self$setParams(params = x)
        self$compute()
      }, x0 = self$getParams())
    },
    
    computeHess = function() {
      stop("Abstract method 2")
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
      stopifnot(R6::is.R6(obj) && inherits(obj, SRBTW_Loss$classname))
      
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



############################################# Below we start to define some
############################################# classes of the future architecture.


Differentiable <- R6Class(
  "Differentiable",
  
  private = list(
    paramNames = NULL,
    numOuts = NULL,
    outputNames = NULL
  ),
  
  public = list(
    initialize = function(
      paramNames = c(), numberOfOutputs = 1,
      outputNames = paste0("o_", seq_len(length.out = numberOfOutputs))
    ) {
      stopifnot((length(paramNames) == 0) || (!any(is.na(paramNames)) && is.character(paramNames)))
      stopifnot(is.numeric(numberOfOutputs) && !is.na(numberOfOutputs) && numberOfOutputs >= 1)
      stopifnot(numberOfOutputs == length(outputNames) && is.character(outputNames) && all(nchar(outputNames) > 0))
      
      if (length(paramNames) > 0 && is.unsorted(x = paramNames, strictly = TRUE)) {
        warning("The parameter names were not sorted, now they are.")
      }
      
      private$paramNames <- sort(paramNames)
      private$numOuts <- numberOfOutputs
      private$outputNames <- outputNames
    },
    
    #' Must return a named vector with all parameters. The names
    #' must be those the Differentiable was initialized with.
    getParams = function() {
      stop("Abstract method 3")
    },
    
    getParamNames = function() {
      private$paramNames
    },
    
    getNumParams = function() {
      length(self$getParamNames())
    },
    
    getOutputNames = function() {
      private$outputNames
    },
    
    getNumOutputs = function() {
      length(self$getOutputNames())
    },
    
    #' Returns a function that takes a single argument, x, that
    #' is the parameters of the Differentiable. That function,
    #' together with the parameters, computes the output.
    get0Function = function() {
      stop("Abstract method 4")
    },
    
    get1stOrderPd = function(name = c()) {
      stop("Abstract method 5")
    },
    
    get2ndOrderPd = function(name = c()) {
      stop("Abstract method 6")
    }
  )
)


Objective <- R6Class(
  "Objective",
  
  inherit = Differentiable,
  
  lock_objects = FALSE,
  
  private = list(
    zeroFunc = function(x) rep(0, length(x))
  ),
  
  public = list(
    initialize = function(
      paramNames = c(), numberOfOutputs = 1,
      outputNames = paste0("o", seq_len(length.out = numberOfOutputs))
    ) {
      super$initialize(
        paramNames = paramNames,
        numberOfOutputs = numberOfOutputs,
        outputNames = outputNames)
    },
    
    #' We need a setter for, e.g., computing the gradient/hessian.
    setParams = function(params) {
      stop("Abstract method 7")
    },
    
    #' Computes the underlying function with the given (or own) parameters.
    compute0 = function() {
      do.call(what = self$get0Function(), args = list())
    },
    
    #' Returns a named vector of computed 1st-order partial derivatives.
    #' The names of that vector correspond to the parameter derived for,
    #' and for the output they were calculated for. The naming scheme is
    #' "parameterName_outputName".
    compute1stOrderPd = function(name = c()) {
      stop("Abstract method 8")
    },
    
    #' Returns a named vector of computed 2nd-order partial derivatives.
    #' The names of that vector correspond to comma-separated lists of
    #' what parameters it was derived for, e.g., "a,b_o1" means that the
    #' partial derivative was first derived for "a", then for "b", and
    #' that it concerns output "o1".
    compute2ndOrderPd = function(name = c()) {
      stop("Abstract method 9")
    },
    
    #' Computes all 1st-order partial derivatives, which is equal to the
    #' gradient. The result is an ordered and named matrix with one row
    #' (according to the order of the parameters of the Differentiable)
    #' for scalar-valued functions, and a Jacobian-matrix for vector-
    #' valued functions (one row per parameter).
    computeGradient = function() {
      # temp <- do.call(what = self$compute1stOrderPd, args = list())
      # # TODO: compose list into matrix/tensor
      # # TODO: Check how I did it in computeGradient_numeric
      stop("not done implementing")
    },
    
    #' Computes all 2nd-order partial derivatives, which is equal to the
    #' Hessian. The result is a 3rd-order Tensor with the 3rd-dimension
    #' having a length of 1 for scalar-valued functions, and n for n-valued
    #' functions. Any of these are ordered according to the order of
    #' the parameters of the Differentiable.
    computeHessian = function() {
      # temp <- do.call(what = self$compute2ndOrderPd, args = list())
      # # TODO: compose list into matrix/tensor
      stop("not done implementing")
    },
    
    #' Computes all 1st-order partial derivatives, which is equal to the
    #' gradient. The result is an ordered and named matrix with one row
    #' (according to the order of the parameters of the Differentiable)
    #' for scalar-valued functions, and a Jacobian-matrix for vector-
    #' valued functions (one row per parameter). Computation is done
    #' numerically using \code{pracma::grad()}.
    computeGradient_numeric = function() {
      g <- matrix(nrow = self$getNumOutputs(), ncol = self$getNumParams())
      colnames(g) <- self$getParamNames()
      rownames(g) <- self$getOutputNames()
      
      for (oIdx in seq_len(length.out = self$getNumOutputs())) {
        g[oIdx, ] <- pracma::grad(f = function(x) {
          self$setParams(params = x)
          self$compute0()[oIdx]
        }, x0 = self$getParams())
      }
      g
    },
    
    #' Computes all 2nd-order partial derivatives, which is equal to the
    #' Hessian. The result is a 3rd-order Tensor with the 3rd-dimension
    #' having a length of 1 for scalar-valued functions, and n for n-valued
    #' functions. Any of these are ordered according to the order of
    #' the parameters of the Differentiable. Computation is done
    #' numerically using \code{pracma::hessian()}.
    computeHessian_numeric = function() {
      # TODO: like computeGradient_numeric, use pracma::hessian and always
      # TODO: return a 3rd-order Tensor that has dim(3) >= 1..
      stop("Not implemented")
    },
    
    #' Similar to \code{compute1stOrderPd()}, but computes each derivative
    #' numerically using \code{pracma::fderiv()}.
    compute1stOrderPd_numeric = function(name = c()) {
      # TODO: use pracma::fderiv() for every parameter, support vector-valued functions like:
      # pracma::fderiv(f = function(x) {
      #   c(x[1]^2, x[2]^3)
      # }, x = c(1,2))
      stop("Not implemented")
    },
    
    #' Similar to \code{compute2ndOrderPd()}, but computes each derivative
    #' numerically using \code{pracma::fderiv(n=2)}.
    compute2ndOrderPd_numeric = function(name = c()) {
      # if (length(names) == 0) {
      #   return(self$computeHessian_numeric())
      # }
      # # TODO: see above, use pracma::fderiv() with n=2
      stop("Not implemented")
    }
  )
)



MultiObjective <- R6Class(
  "MultiObjective",
  
  inherit = Objective,
  
  private = list(
    objectives = NULL,
    
    requireOneObjective = function() {
      stopifnot(length(private$objectives) > 0)
      private$objectives[[names(private$objectives)[1]]]
    }
  ),
  
  public = list(
    initialize = function() {
      # The parameters depend on the objectives added. All added
      # objectives must have the same number of parameters and
      # outputs, and their names must be equal.
      super$initialize(paramNames = c(), numberOfOutputs = 1)
      
      private$objectives <- list()
    },
    
    getNumObjectives = function() {
      length(private$objectives)
    },
    
    # region Objective
    getNumOutputs = function() {
      o <- private$requireOneObjective()
      o$getNumOutputs()
    },
    
    getNumParams = function() {
      o <- private$requireOneObjective()
      o$getNumParams()
    },
    
    getOutputNames = function() {
      o <- private$requireOneObjective()
      o$getOutputNames()
    },
    
    getParamNames = function() {
      o <- private$requireOneObjective()
      o$getParamNames()
    },
    # endregion Objective
    
    setObjective = function(name, obj) {
      stopifnot(is.character(name) && length(name) == 1 && nchar(name) > 0)
      stopifnot(R6::is.R6(obj) && inherits(obj, Objective$classname))
      
      if (self$getNumObjectives() > 0) {
        o <- private$requireOneObjective()
        # Check compatibility:
        stopifnot(o$getNumOutputs() == obj$getNumOutputs())
        stopifnot(o$getNumParams() == obj$getNumParams())
        stopifnot(all.equal(o$getOutputNames(), obj$getOutputNames()))
        stopifnot(all.equal(o$getParamNames(), obj$getParamNames()))
      } else {
        private$numOuts <- obj$getNumOutputs()
        private$outputNames <- obj$getOutputNames()
        private$paramNames <- obj$getParamNames()
      }
      private$objectives[[name]] <- obj
      invisible(self)
    },
    
    hasObjective = function(name) {
      stopifnot(is.character(name) && length(name) == 1 && nchar(name) > 0)
      name %in% names(private$objectives)
    },
    
    removeObjective = function(name) {
      if (!self$hasObjective(name = name)) {
        stop(paste0("The Objective ", name, " is not known."))
      }
      stopifnot(is.character(name) && length(name) == 1 && nchar(name) > 0)
      private$objectives[[name]] <- NULL
      invisible(self)
    }
  )
)


srBTAW_LossLinearScalarizer <- R6Class(
  "srBTAW_LossLinearScalarizer",
  
  inherit = MultiObjective,
  
  private = list(
    computeParallel = NULL,
    gradientParallel = NULL,
    computeWithinGradientParallel = NULL,
    
    isComputingGrad = NULL,
    
    computeFunc = NULL,
    
    progressCallback = NULL,
    
    reportProgress = function(n) {
      # 'n' is the number of done steps and the total amount of
      # steps depends on whether we called compute or the gradient.
      total <- if (private$isComputingGrad) {
        self$getNumParams()
      } else {
        length(private$objectives)
      }
      do.call(what = private$progressCallback, args = list(
        step = n, total = total, what = if (private$isComputingGrad) "Gradient" else "Compute"
      ))
    }
  ),
  
  public = list(
    initialize = function(
      computeParallel = TRUE,
      gradientParallel = TRUE,
      computeWithinGradientParallel = FALSE,
      progressCallback = function(what = c("Compute", "Gradient")[1], step, total) { }
    ) {
      stopifnot(is.logical(computeParallel))
      stopifnot(is.logical(gradientParallel))
      stopifnot(is.logical(computeWithinGradientParallel))
      stopifnot(is.function(progressCallback))
      stopifnot(TRUE == all.equal(c("what", "step", "total"), methods::formalArgs(progressCallback)))
      super$initialize()
      
      private$computeParallel <- computeParallel
      private$gradientParallel <- gradientParallel
      private$computeWithinGradientParallel <- computeWithinGradientParallel
      private$progressCallback <- progressCallback
      private$isComputingGrad <- FALSE
      
      private$computeFunc <- function() {
        private$requireOneObjective()
        
        `%parop%` <- if (private$computeParallel) foreach::`%dopar%` else foreach::`%do%`
        # w_1*L_1 + ... + w_n*L_n
        losses <- foreach::foreach(
          objName = names(private$objectives),
          .combine = c,
          .inorder = FALSE, # does not matter here
          .export = c("private"),
          .verbose = FALSE,
          .options.snow = list(progress = function(n) {
            private$reportProgress(n)
          })
        ) %parop% {
          source("./common-funcs.R")
          source("../models/SRBTW-R6.R")
          obj <- private$objectives[[objName]] # instance of 'srBTAW_Loss'
          obj$getWeight() * obj$compute0()
        }
        
        sum(losses)
      }
    },
    
    getParams = function() {
      private$requireOneObjective()$getParams()
    },
    
    setParams = function(params) {
      for (name in names(private$objectives)) {
        private$objectives[[name]]$getSrBtaw()$setParams(params = params)
      }
      invisible(self)
    },
    
    setComputeParallel = function(computeParallel) {
      stopifnot(is.logical(computeParallel))
      private$computeParallel <- computeParallel
      invisible(self)
    },
    
    setObjective = function(name, obj) {
      stopifnot(inherits(obj, srBTAW_Loss$classname))
      super$setObjective(name = name, obj = obj)
    },
    
    get0Function = function() {
      private$computeFunc
    },
    
    #' Overridden so that we can work in parallel!
    computeGradient_numeric = function() {
      private$isComputingGrad <- TRUE

      if (!private$gradientParallel) {
        return(super$computeGradient_numeric())
      }

      g <- matrix(nrow = self$getNumOutputs(), ncol = self$getNumParams())
      colnames(g) <- self$getParamNames()
      rownames(g) <- self$getOutputNames()

      `%parop%` <- if (private$computeParallel) foreach::`%dopar%` else foreach::`%do%`

      for (oIdx in seq_len(length.out = self$getNumOutputs())) {
        g[oIdx, ] <- foreach::foreach(
          pIdx = seq_len(length.out = self$getNumParams()),
          .combine = c,
          .inorder = TRUE,
          .packages = c("pracma"),
          .export = c("self"),
          .options.snow = list(progress = function(n) {
            private$reportProgress(n)
          })
        ) %parop% {
          source("./common-funcs.R")
          source("../models/SRBTW-R6.R")
          # The following is not very elegant but it works!
          # We need to make a clone of the model and this objective:
          copy <- self$clone()
          srbtaw <- private$requireOneObjective()$getSrBtaw()$clone()
          srbtaw$setObjective(obj = copy)
          
          # TODO: this requires nested parallelization..
          copy$setComputeParallel(computeParallel = private$computeWithinGradientParallel)
          params <- srbtaw$getParams()

          gr <- function(x) {
            params[pIdx] <- x # hold everything else constant
            srbtaw$setParams(params = params)
            copy$compute0()[oIdx]
          }

          pracma::fderiv(
            f = gr, x = params[pIdx], n = 1, method = "central")
        }
      }

      private$isComputingGrad <- FALSE
      g
    }
  )
)



#' Abstract super-class for all models. This class is meant
#' to represent models that were fit previously -- that means
#' that there are concrete values for all parameters, as well
#' as the likelihood.
Model <- R6Class(
  "Model",
  
  inherit = Differentiable,
  
  private = list(
    params = NULL
  ),
  
  public = list(
    initialize = function(
      paramNames = c(), numberOfOutputs = 1,
      outputNames = paste0("o_", seq_len(length.out = numberOfOutputs))
    ) {
      super$initialize(
        paramNames = paramNames, numberOfOutputs = numberOfOutputs, outputNames = outputNames)
    },
    
    #' Implements abstract method from \code{Differentiable::getParams()}.
    getParams = function() {
      private$params
    },
    
    setParams = function(params) {
      stopifnot(is.vector(params) && length(params) == self$getNumParams())
      stopifnot(is.numeric(params) && !any(is.na(params)))
      allEq <- all.equal(self$getParamNames(), names(params))
      stopifnot(is.logical(allEq) && allEq)
      private$params <- params
      invisible(self)
    },
    
    likelihood = function() {
      stop("Abstract method 10")
    },
    
    residuals = function(...) {
      stop("Abstract method 11")
    },
    
    coefficients = function() {
      self$getParams()
    },
    
    fit = function() {
      stop("Abstract method 12")
    },
    
    AIC = function() {
      # 2*k - 2*ln(L)
      2 * self$getNumParams() - 2 * log(self$likelihood())
    },
    
    BIC = function() {
      # k*ln(n) - 2*ln(L),
      k <- self$getNumParams()
      L <- self$likelihood()
      # n = the number of data points, the number of observations, or equivalently, the sample size ->
      # n probably only makes sense for when all associated objectives are finite and discrete..
      stop("Abstract method 13")
    }
  )
)


srBTAW_Loss <- R6Class(
  "srBTAW_Loss",
  
  inherit = Objective,
  
  private = list(
    srbtaw = NULL,
    wpName = NULL,
    wcName = NULL,
    weight = NULL,
    intervals = NULL,
    minimize = NULL
  ),
  
  #' Arguments for the names of the signals should be NULL if the loss
  #' is not specific to a pair of variables.
  public = list(
    initialize = function(wpName = NULL, wcName = NULL, weight = 1, intervals = c(), minimize = TRUE) {
      stopifnot(is.null(wpName) || (is.character(wpName) && nchar(wpName) > 0))
      stopifnot(is.null(wcName) || (is.character(wcName) && nchar(wcName) > 0))
      stopifnot(is.numeric(weight) && !is.na(weight) && weight > 0 && weight <= 1)
      stopifnot(is.numeric(intervals) && !any(is.na(intervals)) && all(intervals >= 1))
      stopifnot(!is.unsorted(intervals))
      stopifnot(is.logical(minimize))
      
      private$wpName <- wpName
      private$wcName <- wcName
      private$weight <- weight
      private$intervals <- intervals
      private$minimize <- minimize
    },
    
    setParams = function(params) {
      srbtaw <- private$srbtaw
      stopifnot(R6::is.R6(srbtaw) && inherits(srbtaw, srBTAW$classname))
      srbtaw$setParams(params = params)
      invisible(self)
    },
    
    setSrBtaw = function(srbtaw) {
      private$srbtaw <- srbtaw
      invisible(self)
    },
    
    getSrBtaw = function() {
      private$srbtaw
    },
    
    getWpName = function() {
      private$wpName
    },
    
    getWcName = function() {
      private$wcName
    },
    
    getWeight = function() {
      private$weight
    },
    
    getIntervals = function() {
      private$intervals
    },
    
    setMinimization = function(val) {
      stopifnot(is.logical(val))
      private$minimize <- val
      invisible(self)
    },
    
    isMinimization = function() {
      private$minimze
    },
    
    
    getParamNames = function() {
      self$getSrBtaw()$getParamNames()
    },
    
    getOutputNames = function() {
      self$getSrBtaw()$getOutputNames()
    },
    
    getNumParams = function() {
      self$getSrBtaw()$getNumParams()
    },
    
    getParams = function() {
      self$getSrBtaw()$getParams()
    }
  )
)


srBTAW_LossWrapper <- R6Class(
  "srBTAW_LossWrapper",
  
  inherit = srBTAW_Loss,
  
  private = list(
    objective = NULL
  ),
  
  public = list(
    initialize = function(objective, weight = 1, minimize = TRUE) {
      super$initialize(
        wpName = NULL, wcName = NULL, weight = weight, intervals = c(1))
      stopifnot(R6::is.R6(objective) && inherits(obj, Objective$classname))
      private$objective <- objective
      
      self$setMinimization(val = minimize)
    },
    
    get0Function = function() {
      private$objective$get0Function()
    },
    
    get1stOrderPd = function(name = c()) {
      private$objective$get1stOrderPd(name = name)
    },
    
    get2ndOrderPd = function(name = c()) {
      private$objective$get2ndOrderPd(name = name)
    },
    
    compute0 = function() {
      private$objective$compute0()
    },
    
    getNumOutputs = function() {
      private$objective$getNumOutputs()
    },
    
    getNumParams = function() {
      private$objective$getNumParams()
    },
    
    getOutputNames = function() {
      private$objective$getOutputNames()
    },
    
    getParamNames = function() {
      private$objective$getParamNames()
    },
    
    getParams = function() {
      private$objective$getParams()
    },
    
    setParams = function(params) {
      private$objective$setParams(params = params)
    }
  )
)



srBTAW_Loss2Curves <- R6Class(
  "srBTAW_Loss2Curves",
  
  inherit = srBTAW_Loss,
  
  private = list(
    continuous = NULL,
    numSamples = NULL,
    use = NULL
  ),
  
  
  public = list(
    initialize = function(
      srbtaw, wpName, wcName, weight = 1, intervals = c(),
      use = c("area", "rss", "corr", "arclen", "jsd"),
      continuous = FALSE, numSamples = rep(if (continuous) NA_real_ else 1e3, length(intervals))
    ) {
      super$initialize(
        wpName = wpName, wcName = wcName, weight = weight, intervals = intervals)
      self$setSrBtaw(srbtaw = srbtaw)
      
      stopifnot(is.logical(continuous) && is.numeric(numSamples))
      
      private$continuous <- continuous
      private$numSamples <- numSamples
      private$use <- match.arg(use)
    },
    
    getNumOutputs = function() {
      1
    },
    
    funcArea = function() {
      continuous <- private$continuous
      srbtaw <- private$srbtaw
      qs <- private$intervals
      
      err <- 0
      res <- srbtaw$residuals(loss = self)
      idx <- 1
      for (q in qs) {
        t <- res[[q]]
        wc <- if (srbtaw$isBawEnabled()) t$nqc else t$mqc
        
        if (continuous) {
          err <- err + cubature::cubintegrate(
            f = function(x) abs(t$wp(x) - wc(x)), lower = t$tb_q, upper = t$te_q)$integral
        } else {
          X <- seq(from = t$tb_q, to = t$te_q, length.out = private$numSamples[idx])
          y <- sapply(X = X, FUN = t$wp)
          y_hat <- sapply(X = X, FUN = wc)
          err <- err + (sum(abs(y - y_hat)) / private$numSamples[idx])
        }
        
        idx <- idx + 1
      }
      
      `names<-`(c(log(1 + err)), srbtaw$getOutputNames())
      
    },
    
    funcCorr = function() {
      continuous <- private$continuous
      stopifnot(!continuous) # not currently supported/used here
      
      srbtaw <- private$srbtaw
      qs <- private$intervals
      
      corrs <- c() # we will calculate the mean correlation
      res <- srbtaw$residuals(loss = self)
      idx <- 1
      for (q in qs) {
        t <- res[[q]]
        wc <- if (srbtaw$isBawEnabled()) t$nqc else t$mqc
        
        X <- seq(from = t$tb_q, to = t$te_q, length.out = private$numSamples[idx])
        y <- sapply(X = X, FUN = t$wp)
        y_hat <- sapply(X = X, FUN = wc)
        temp <- suppressWarnings(expr = {
          stats::cor(x = y, y = y_hat)
        })
        if (is.na(temp)) {
          # Happens when stddev is 0; this happens when one (or both)
          # signal(s) is constant (straight line). In this case, we
          # apply the rule of interpreting this as no correlation, as
          # it could be even worse (negative correlation).
          temp <- 0
        }
        corrs <- c(corrs, temp)
        
        idx <- idx + 1
      }
      
      err <- 0.5 * (1 - mean(corrs)) # \mapsto [0,1], where 1 is the highest possible loss
      # re-scale error so it is < 1 (otherwise we'll get infinity)
      err <- err * (1 - sqrt(.Machine$double.eps))
      
      `names<-`(c(-log(1 - err)), srbtaw$getOutputNames())
    },
    
    funcArclen = function() {
      arclen <- function(func, from, to) {
        func <- Vectorize(FUN = func)
        pracma::arclength(f = function(k) c(func(k), k), a = from, b = to)$length
      }
      
      srbtaw <- private$srbtaw
      qs <- private$intervals
      
      arclens <- c() # we will calculate the mean arclength ratio
      res <- srbtaw$residuals(loss = self)
      for (q in qs) {
        t <- res[[q]]
        wc <- if (srbtaw$isBawEnabled()) t$nqc else t$mqc
        
        arclen_wp <- arclen(func = t$wp, from = t$tb_q, to = t$te_q)
        arclen_wc <- arclen(func = wc, from = t$tb_q, to = t$te_q)
        temp <- 1 - (min(arclen_wp, arclen_wc) / max(arclen_wp, arclen_wc))
        # \mapsto [0,1], where 1 is the maximum possible loss (0 means a perfect ratio)
        
        arclens <- c(arclens, temp)
      }
      
      `names<-`(c(-log(1 - mean(arclens))), srbtaw$getOutputNames())
    },
    
    funcRss = function() {
      continuous <- private$continuous
      srbtaw <- private$srbtaw
      qs <- private$intervals
      
      err <- 0
      res <- srbtaw$residuals(loss = self)
      idx <- 1
      for (q in qs) {
        t <- res[[q]]
        wc <- if (srbtaw$isBawEnabled()) t$nqc else t$mqc
        
        if (continuous) {
          err <- err + cubature::cubintegrate(
            f = function(x) (t$wp(x) - wc(x))^2, lower = t$tb_q, upper = t$te_q)$integral
        } else {
          X <- seq(from = t$tb_q, to = t$te_q, length.out = private$numSamples[idx])
          y <- sapply(X = X, FUN = t$wp)
          y_hat <- sapply(X = X, FUN = wc)
          err <- err + (sum((y - y_hat)^2) / private$numSamples[idx])
        }
        
        idx <- idx + 1
      }
      
      `names<-`(c(log(1 + err)), srbtaw$getOutputNames())
    },
    
    get0Function = function() {
      if (private$use == "rss") {
        self$funcRss
      } else if (private$use == "area") {
        self$funcArea
      } else if (private$use == "corr") {
        self$funcCorr
      } else if (private$use == "arclen") {
        self$funcArclen
      } else if (private$use == "jsd") {
        self$funcJsd
      } else {
        stop(paste0("Should not get here.."))
      }
    }
  )
)


srBTAW_Loss_Rss <- R6Class(
  "srBTAW_Loss_Rss",
  
  inherit = srBTAW_Loss,
  
  private = list(
    continuous = NULL,
    numSamples = NULL
  ),
  
  public = list(
    initialize = function(
      wpName, wcName, weight = 1, intervals = c(),
      continuous = FALSE, numSamples = rep(if (continuous) NA_real_ else 1e3, length(intervals))
    ) {
      super$initialize(
        wpName = wpName, wcName = wcName, weight = weight, intervals = intervals)
      
      stopifnot(is.logical(continuous) && is.numeric(numSamples))
      
      private$continuous <- continuous
      private$numSamples <- numSamples
    },
    
    getNumOutputs = function() {
      1
    },
    
    get0Function = function() {
      function() {
        # Now for each interval, we compute the RSS over the
        # residuals from the MLM. The residuals currently
        # are just the sub-models as tuple so we can do
        # whatever we want.
        
        continuous <- private$continuous
        srbtaw <- private$srbtaw
        qs <- private$intervals
        
        err <- 0
        res <- srbtaw$residuals(loss = self)
        idx <- 1
        for (q in qs) {
          t <- res[[q]]
          wc <- if (srbtaw$isBawEnabled()) t$nqc else t$mqc
          
          if (continuous) {
            err <- err + cubature::cubintegrate(
              f = function(x) (t$wp(x) - wc(x))^2, lower = t$tb_q, upper = t$te_q)$integral
          } else {
            X <- seq(from = t$tb_q, to = t$te_q, length.out = private$numSamples[idx])
            y <- sapply(X = X, FUN = t$wp)
            y_hat <- sapply(X = X, FUN = wc)
            err <- err + (sum((y - y_hat)^2) / private$numSamples[idx])
          }
          
          idx <- idx + 1
        }
        
        `names<-`(c(log(1 + err)), srbtaw$getOutputNames())
      }
    }
  )
)



Signal <- R6Class(
  "Signal",
  
  inherit = Differentiable,
  
  private = list(
    name = NULL,
    func = NULL,
    support = NULL,
    isWp = NULL,
    
    func_d1 = NULL,
    func_d2 = NULL
  ),
  
  public = list(
    initialize = function(name, func, support, isWp) {
      stopifnot(is.character(name) && nchar(name) > 0)
      stopifnot(is.function(func) && length(methods::formalArgs(func)) == 1)
      stopifnot(is.numeric(support) && length(support) == 2 && !any(is.na(support)))
      stopifnot(is.logical(isWp))
      
      private$name <- name
      private$func <- Vectorize(func)
      private$support <- support
      private$isWp <- isWp
      
      func_deriv <- function(x, n) {
        m <- "central"
        eps <- sqrt(.Machine$double.eps)
        if (abs(x - support[1]) < eps) {
          m <- "forward"
        } else if (abs(support[2] - x) < eps) {
          m <- "backward"
        }
        pracma::fderiv(f = func, x = x, n = n, method = m)
      }
      
      private$func_d1 <- Vectorize(function(x) {
        func_deriv(x = x, n = 1)
      })
      private$func_d2 <- Vectorize(function(x) {
        func_deriv(x = x, n = 2)
      })
    },
    
    getName = function() {
      private$name
    },
    
    getSupport = function() {
      private$support
    },
    
    isWarpingPattern = function() {
      private$isWp
    },
    
    
    # region Differentiable
    get0Function = function() {
      private$func
    },
    
    get1stOrderPd = function(name = c()) {
      # 'name' must be an empty vector or 'x'
      stopifnot(length(name) == 0 || (length(name) == 1 && name == "x"))
      private$func_d1
    },
    
    get2ndOrderPd = function(name = c()) {
      # 'name' must be an empty vector or 'x'
      stopifnot(length(name) == 0 || (length(name) == 1 && name == "x"))
      private$func_d2
    },
    
    getParams = function() {
      # Signal does not store a value for x (i.e., there is no setter).
      c(x = NA_real_)
    },
    
    getParamNames = function() {
      "x"
    },
    
    getNumParams = function() {
      1
    },
    
    getOutputNames = function() {
      "y"
    },
    
    getNumOutputs = function() {
      1
    },
    # endregion Differentiable
    
    
    plot = function(showSignal = TRUE, show1stDeriv = FALSE, show2ndDeriv = FALSE) {
      funcsColors <- c(Signal = "black", Signal_d1 = "red", Signal_d2 = "blue")
      funcs <- c()
      g <- ggplot2::ggplot()
      
      if (showSignal) {
        funcs <- c(1)
        g <- g + ggplot2::stat_function(
          fun = self$get0Function(), xlim = private$support,
          mapping = ggplot2::aes(color = names(funcsColors)[1]))
      }
      if (show1stDeriv) {
        funcs <- c(funcs, 2)
        g <- g + ggplot2::stat_function(
          fun = self$get1stOrderPd(), xlim = private$support,
          mapping = ggplot2::aes(color = names(funcsColors)[2]))
      }
      if (show2ndDeriv) {
        funcs <- c(funcs, 3)
        g <- g + ggplot2::stat_function(
          fun = self$get2ndOrderPd(), xlim = private$support,
          mapping = ggplot2::aes(color = names(funcsColors)[3]))
      }
      
      g + ggplot2::scale_color_manual(
        private$name,
        values = funcsColors[funcs]) +
        ggplot2::theme(legend.position = "bottom")
    }
  )
)



srBTAW <- R6Class(
  "srBTAW",
  
  inherit = Model,
  
  private = list(
    objective = NULL,
    
    useAmplitudeWarping = NULL,
    isObjectiveLogarithmic = NULL,
    
    theta_b = NULL, gamma_bed = NULL, lambda = NULL,
    begin = NULL, end = NULL, openBegin = NULL, openEnd = NULL,
    lambda_ymin = NULL, lambda_ymax = NULL,
    
    lastFitResult = NULL,
    
    instances = NULL,
    
    signals_wp = NULL, signals_wc = NULL,
    
    data = NULL, dataWeights = NULL,
    
    requireOneInstance = function() {
      stopifnot(length(private$instances) > 0)
      useName <- strsplit(x = names(private$instances)[1], split = "|", fixed = TRUE)[[1]]
      self$getInstance(wpName = useName[1], wcName = useName[2])
    },
    
    requireObjective = function() {
      if (!self$hasObjective()) stop("No Objective present.")
    },
    
    addInstance = function(wpName, wcName, instance) {
      key <- paste(wpName, wcName, sep = "|")
      private$instances[[key]] <- instance
    },
    
    createInstance = function(wp, wc) {
      ctor <- if (private$useAmplitudeWarping) SRBTWBAW$new else SRBTW$new
      useArgs <- list(
        wp = wp, wc = wc, theta_b = private$theta_b, gamma_bed = private$gamma_bed,
        lambda = private$lambda, begin = private$begin, end = private$end,
        openBegin = private$openBegin, openEnd = private$openEnd)
      
      if (private$useAmplitudeWarping) {
        useArgs$lambda_ymin <- private$lambda_ymin
        useArgs$lambda_ymax <- private$lambda_ymax
      }
      
      do.call(what = ctor, args = useArgs)
    }
  ),
  
  public = list(
    initialize = function(
      theta_b, gamma_bed, lambda, begin, end, openBegin, openEnd, # required BTW/BAW params
      useAmplitudeWarping = FALSE,
      lambda_ymin = NULL, lambda_ymax = NULL # required if BAW
      
      ,paramNames = c(), numberOfOutputs = 1,
      outputNames = paste0("o_", seq_len(length.out = numberOfOutputs)),
      
      isObjectiveLogarithmic = TRUE
    ) {
      stopifnot(!useAmplitudeWarping || (is.vector(lambda_ymin) && is.vector(lambda_ymax)))
      
      private$instances <- list()
      private$signals_wp <- list()
      private$signals_wc <- list()
      
      private$useAmplitudeWarping <- useAmplitudeWarping
      private$theta_b <- theta_b
      private$gamma_bed <- gamma_bed
      private$lambda <- lambda
      private$begin <- begin
      private$end <- end
      private$openBegin <- openBegin
      private$openEnd <- openEnd
      private$lambda_ymin <- lambda_ymin
      private$lambda_ymax <- lambda_ymax
      
      private$isObjectiveLogarithmic <- isObjectiveLogarithmic
      
      paramNames <- paste0("vtl_", seq_len(length.out = length(private$theta_b) - 1)) # vartheta_l
      if (openBegin) {
        paramNames <- c(paramNames, "b")
      }
      if (openEnd) {
        paramNames <- c(paramNames, "e")
      }
      if (useAmplitudeWarping) {
        # v, vartheta_y
        paramNames <- c(paramNames, "v", paste0("vty_", seq_len(length.out = length(private$theta_b) - 1)))
      }
      super$initialize(paramNames = sort(paramNames))
    },
    
    isBawEnabled = function() {
      private$useAmplitudeWarping
    },
    
    setIsObjectiveLogarithmic = function(val) {
      stopifnot(is.logical(val))
      private$isObjectiveLogarithmic <- val
      invisible(self)
    },
    
    getIsObjectiveLogarithmic = function() {
      private$isObjectiveLogarithmic
    },
    
    setSignal = function(signal) {
      stopifnot(R6::is.R6(signal) && inherits(signal, Signal$classname))
      
      if (signal$isWarpingPattern()) {
        private$signals_wp[[signal$getName()]] <- signal
      } else {
        private$signals_wc[[signal$getName()]] <- signal
      }
      
      invisible(self)
    },
    
    getSignal = function(name) {
      stopifnot(is.character(name))
      
      if (name %in% names(private$signals_wp)) {
        return(private$signals_wp[[name]])
      } else if (name %in% names(private$signals_wc)) {
        return(private$signals_wc[[name]])
      }
      stop(paste0("The signal ", name, " was not previously added."))
    },
    
    removeSignal = function(nameOrSignal) {
      stopifnot(is.character(nameOrSignal) || (R6::is.R6(signal) && inherits(signal, Signal$classname)))
      
      name <- if (is.character(nameOrSignal)) nameOrSignal else nameOrSignal$getName()
      if (name %in% names(private$signals_wp)) {
        private$signals_wp[[name]] <- NULL
      } else if (name %in% names(private$signals_wc)) {
        private$signals_wc[[name]] <- NULL
      } else {
        stop(paste0("The signal ", name, " was not previously added."))
      }
      
      invisible(self)
    },
    
    setObjective = function(obj = NULL) {
      stopifnot(is.null(obj) || (R6::is.R6(obj) && inherits(obj, Objective$classname)))
      private$objective <- obj
      invisible(self)
    },
    
    getObjective = function() {
      private$objective
    },
    
    hasObjective = function() {
      R6::is.R6(private$objective) && inherits(private$objective, Objective$classname)
    },
    
    hasInstance = function(wpName, wcName) {
      key <- paste(wpName, wcName, sep = "|")
      key %in% names(private$instances)
    },
    
    getInstance = function(wpName, wcName) {
      stopifnot(self$hasInstance(wpName = wpName, wcName = wcName))
      key <- paste(wpName, wcName, sep = "|")
      private$instances[[key]]
    },
    
    addLoss = function(loss) {
      stopifnot(R6::is.R6(loss) && inherits(loss, srBTAW_Loss$classname))
      
      wpName <- loss$getWpName()
      wcName <- loss$getWcName()
      if (!is.null(wpName) && !is.null(wcName) && !self$hasInstance(wpName = wpName, wcName = wcName)) {
        wpSig <- self$getSignal(name = wpName)
        stopifnot(wpSig$isWarpingPattern())
        wcSig <- self$getSignal(name = wcName)
        stopifnot(!wcSig$isWarpingPattern())
        
        private$addInstance(
          wpName = wpName, wcName = wcName,
          instance = private$createInstance(
            wp = wpSig$get0Function(), wc = wcSig$get0Function()))
        
        # Also, make sure the new instance gets may existing parameters:
        self$setParams(params = self$getParams())
      }
      
      loss$setSrBtaw(srbtaw = self)
      invisible(self)
    },
    # region Model
    
    #' Overridden so that we can set parameters to all instances.
    setParams = function(params) {
      super$setParams(params = params) # called also for validation
      
      # 'params' is a named and sorted vector, and we need to deconstruct
      # it into vartheta_l, b, e, v, vartheta_y and set these explicitly,
      # if they are present.
      useArgs <- list()
      patterns <- c("vtl_" = "vartheta_l", "b" = "begin", "e" = "end",
                    "v$" = "v", "vty_" = "vartheta_y")
      for (p in names(patterns)) {
        idx <- grep(pattern = p, x = names(params))
        if (length(idx) > 0) {
          temp <- params[idx]
          useArgs[[patterns[p]]] <- unname(temp[sort(names(temp))])
        }
      }
      
      for (instName in names(private$instances)) {
        do.call(what = private$instances[[instName]]$setParams, args = useArgs)
      }
      
      invisible(self)
    },
    
    likelihood = function() {
      fr <- private$lastFitResult
      stopifnot(R6::is.R6(fr) && inherits(fr, FitResult$classname))
      
      # Likelihood and loss are anti-proportional, the lower
      # the loss, the higher the likelihood.
      bestLoss <- fr$getBest(paramName = "loss", lowest = TRUE)[, "loss"]
      # The following has limit \infty for loss -> 0 and limit 0
      # for loss -> \infty (what we want)
      # TODO: Check if Objective is logarithmic
      -log(1 - 1/(1 + bestLoss))
    },
    
    #' This method is the most important for computing Objectives and Losses,
    #' as these are based on the performed warping. This method expects one
    #' parameter, a loss which must be an instance of \code{srBTAW_Loss}.
    residuals = function(...) {
      private$requireObjective()
      
      params <- list(...)
      stopifnot(length(params) > 0)
      loss <- utils::head(params, 1)[[1]]
      stopifnot(R6::is.R6(loss) && inherits(loss, srBTAW_Loss$classname))
      
      qs <- loss$getIntervals()
      stopifnot(length(qs) > 0) # You can request multiple q!
      
      instance <- self$getInstance(
        wpName = loss$getWpName(), wcName = loss$getWcName())
      
      res <- list()
      for (q in qs) {
        res[[q]] <- instance$getSubModel(q = q)$asTuple()
      }
      
      res
    },
    
    fit = function(verbose = TRUE) {
      private$requireObjective()
      obj <- self$getObjective()
      
      vtl <- private$requireOneInstance()$getVarthetaL()
      if (sum(vtl) > 1) {
        stop("To begin fitting, vartheta_l should sum up to 1.")
      }
      
      fr <- FitResult$new(paramNames = c(self$getParamNames(), "loss"))
      fr$start()
      
      # Bounds: for SRBTWBAW, params is vartheta_l, [, b [, e [, v, vartheta_y]]]
      bb_lower <- c(private$lambda,
                    (if (private$openBegin) private$gamma_bed[1] else c()),
                    (if (private$openEnd) private$gamma_bed[1] + private$gamma_bed[3] else c()),
                    (if (private$useAmplitudeWarping)
                      c(rep(min(private$lambda_ymin) - 1e3 * (max(private$lambda_ymax) - min(private$lambda_ymin)),
                            length(private$theta_b))) else c()))
      bb_upper <- c(rep(1, length(private$lambda)), # remember vartheta_l are ratios ideally
                    (if (private$openBegin) private$gamma_bed[2] - private$gamma_bed[3] else c()),
                    (if (private$openEnd) private$gamma_bed[2] else c()),
                    (if (private$useAmplitudeWarping)
                      c(rep(max(private$lambda_ymax) + 1e3 * (max(private$lambda_ymax) - min(private$lambda_ymin)),
                            length(private$theta_b))) else c()))
      
      optR <- stats::optim(
        par = self$getParams(),
        method = "L-BFGS-B",
        lower = bb_lower,
        upper = bb_upper,
        fn = function(x) {
          self$setParams(params = x)
          if (inherits(obj, srBTAW_Loss$classname)) {
            obj$setParams(params = x)
          }
          # otherwise, the method is abstract and the Objective has to have
          # a way of obtaining parameters implemented.
          
          fr$startStep()
          loss <- obj$compute0()
          fr$stopStep(resultParams = c(x, loss), verbose = verbose)
          loss
        },
        gr = function(x) {
          self$setParams(params = x)
          if (inherits(obj, srBTAW_Loss$classname)) {
            obj$setParams(params = x)
          }
          
          # Note that computeGradient_numeric returns a matrix, and
          # we are dealing with scalar-valued losses only.
          obj$computeGradient_numeric()[1, ]
        }
      )
      
      fr$setOptResult(optResult = optR)
      fr$finish()
      private$lastFitResult <- fr
      fr
    },
    # endregion Model
    
    
    # region Differentiable
    getNumParams = function() {
      n <- length(private$theta_b) - 1 # vartheta_l
      if (private$openBegin) {
        n <- n + 1
      }
      if (private$openEnd) {
        n <- n + 1
      }
      if (private$useAmplitudeWarping) {
        n <- n + length(private$theta_b) # v, vartheta_y
      }
      n
    },
    
    get0Function = function() {
      private$requireObjective()
      do.call(what = self$getObjective()$get0Function, args = list())
    },
    
    get1stOrderPd = function(name = c()) {
      private$requireObjective()
      do.call(what = self$getObjective()$get1stOrderPd, args = list(name = name))
    },
    
    get2ndOrderPd = function(name = c()) {
      private$requireObjective()
      do.call(what = self$getObjective()$get2ndOrderPd, args = list(name = name))
    }
    # endregion Differentiable
  )
)



srBTAW_MultiVartype <- R6Class(
  "srBTAW_MultiVartype",
  
  inherit = Model,
  
  private = list(
    instances = NULL,
    
    requireOneInstance = function() {
      stopifnot(length(private$instances) > 0)
      private$instances[[names(private$instances)[1]]]
    }
  ),
  
  public = list(
    initialize = function() {
      super$initialize()
      
      private$instances <- list()
    },
    
    setSrbtaw = function(varName, srbtaw) {
      stopifnot(is.character(varName) && length(varName) == 1)
      stopifnot(grepl(pattern = "^[a-z][a-z0-9]*$", x = varName, ignore.case = TRUE))
      stopifnot(R6::is.R6(srbtaw) && inherits(srbtaw, srBTAW$classname))
      
      if (length(private$instances) > 0) {
        # All need to have the same amount of parameters, and
        # their names must be the same.
        temp <- private$requireOneInstance()
        allEq <- all.equal(temp$getParamNames(), srbtaw$getParamNames())
        stopifnot(is.logical(allEq) && allEq)
      }
      
      private$instances[[varName]] <- srbtaw
      invisible(self)
    },
    
    # region Differentiable
    getParamNames = function() {
      temp <- private$requireOneInstance()
      
      # We have equally many vtl and vty, and one additional v
      # per instance. All instances share vtl, but not v, vty_*
      pn <- temp$getParamNames()
      pn <- pn[2:(1 + ((length(pn) - 1) / 2))] # vtl's
      
      for (name in sort(names(private$instances))) {
        pn <- c(pn, paste(c(temp$getParamNames()[1],
          utils::tail(temp$getParamNames(), (temp$getNumParams() - 1) / 2)), name, sep = "_"))
      }
      
      sort(pn)
    },
    # endregion
    
    # region Objective
    setParams = function(params) {
      stopifnot(is.vector(params) && length(params) == self$getNumParams())
      stopifnot(is.numeric(params) && !any(is.na(params)))
      allEq <- all.equal(self$getParamNames(), names(params))
      stopifnot(is.logical(allEq) && allEq)
      
      # params must be an ordered and named vector, and we expect these params:
      # v_[varname], vtl_1 .. vtl_n, vty_1_[varname] .. vty_n_[varname]
      # (the vty-subvector is repeated for each instance)
      
      pAll <- params[grepl(pattern = "^vtl_\\d+", x = names(params))]
      pAll <- pAll[sort(names(pAll))]
      
      for (name in names(private$instances)) {
        p_v_vty <- params[grepl(pattern = paste0("_", name, "$"), x = names(params))]
        p_v_vty <- p_v_vty[sort(names(p_v_vty))]
        names(p_v_vty) <- gsub(pattern = paste0("_", name, "$"), replacement = "", x = names(p_v_vty))
        useParams <- c(pAll, p_v_vty)
        useParams <- useParams[order(names(useParams))]
        
        do.call(what = private$instances[[name]]$setParams,
                args = list(params = useParams))
      }
      
      private$params <- params
      
      invisible(self)
    }
    # endregion
  )
)



TimeWarpRegularization <- R6Class(
  "TimeWarpRegularization",
  
  inherit = srBTAW_Loss,
  
  private = list(
    use = NULL,
    
    exintKappa = NULL
  ),
  
  #' The names of the signals are required so that we know which
  #' model's parameters to retrieve.
  public = list(
    initialize = function(wpName, wcName, intervals = c(), weight = 1, use = c("exsupp", "exint"), exintKappa = NULL) {
      super$initialize(wpName = wpName, wcName = wcName, intervals = intervals, weight = weight)
      private$use <- match.arg(use)
      private$exintKappa <- exintKappa
      
      if (private$use == "exint") {
        stopifnot(is.character(wpName) && is.character(wcName))
        stopifnot(is.null(exintKappa) || (is.numeric(exintKappa) && !any(is.na(exintKappa))))
        stopifnot(length(intervals) > 0)
      } else {
        if (length(intervals) > 0) {
          warning("The regularizer for extreme supports does not use intervals.")
        }
      }
    },
    
    funcExsupp = function() {
      srbtaw <- private$srbtaw
      mod <- srbtaw$getInstance(wpName = self$getWpName(), wcName = self$getWcName())
      gamma_bed <- mod$getgamma_bed()
      beta_l <- mod$getBeta_l()
      beta_u <- mod$getBeta_u()
      
      `names<-`(
        c(-log((beta_u - beta_l) / (gamma_bed[2] - gamma_bed[1] - gamma_bed[3]))),
        srbtaw$getOutputNames())
    },
    
    funcExint = function() {
      srbtaw <- private$srbtaw
      mod <- srbtaw$getInstance(wpName = self$getWpName(), wcName = self$getWcName())
      vtl <- mod$getVarthetaL()
      
      exintKappa <- private$exintKappa
      if (is.null(exintKappa)) {
        gamma_bed <- mod$getgamma_bed()
        exintKappa <- rep((gamma_bed[2] - gamma_bed[1]) / length(vtl), length(vtl))
      }
      stopifnot(length(exintKappa) == length(vtl))
      
      qs <- private$intervals
      ints <- c()
      res <- srbtaw$residuals(loss = self)
      for (q in qs) {
        ints <- c(ints, (vtl[q] - exintKappa[q])^2)
      }
      
      `names<-`(c(log(1 + sum(ints))), srbtaw$getOutputNames())
    },
    
    getNumOutputs = function() {
      1
    },
    
    get0Function = function() {
      if (private$use == "exsupp") {
        self$funcExsupp
      } else {
        self$funcExint
      }
    }
  )
)



YTransRegularization <- R6Class(
  "YTransRegularization",
  
  inherit = srBTAW_Loss,
  
  private = list(
    use = NULL
  ),
  
  public = list(
    initialize = function(wpName, wcName, intervals = c(), weight = 1, use = c("trapezoid", "tikhonov", "avgout")) {
      super$initialize(wpName = wpName, wcName = wcName, intervals = intervals, weight = weight)
      
      private$use <- match.arg(use)
      
      # For Average-out regularization, wpName & wcName are required,
      # because it requires access to the underlying SRBTW/SRBTWBAW's M-function.
      # They are also required for trapezoidal regularization, as we need to get
      # the residuals, which always require a specific model (even though any
      # model would do here).
      if (private$use == "avgout" || private$use == "trapezoid") {
        stopifnot(is.character(wpName) && nchar(wpName) > 0)
        stopifnot(is.character(wcName) && nchar(wcName) > 0)
      } else {
        # For the others, they better be NULL, because the others are
        # NOT specific to any one model.
        stopifnot(is.null(wpName) && is.null(wcName))
      }
    },
    
    getNumOutputs = function() {
      1
    },
    
    funcTrapezoid = function() {
      srbtaw <- private$srbtaw
      qs <- private$intervals
      
      ratios <- 0
      res <- srbtaw$residuals(loss = self)
      params <- srbtaw$getParams()
      v <- params["v"]
      for (q in qs) {
        t <- res[[q]]
        wc <- if (srbtaw$isBawEnabled()) t$nqc else t$mqc
        o_q <- v + t$phi_y_q
        beta_l <- t$lambda_ymin_q
        beta_u <- t$lambda_ymax_q
        delta_max_q <- (beta_u - beta_l) * (t$te_q - t$tb_q)
        
        # We have 6 cases, let's check 1-by-1 and add to the ratios:
        delta_q <- 0
        if (o_q < beta_l) {
          if ((o_q + t$vty_q) > beta_l) {
            s <- t$tb_q + ((beta_l - o_q) / t$a_q)
            delta_q <- cubature::cubintegrate(
              f = function(x) min(beta_u, t$tq(x)) - beta_l, lower = s, upper = t$te_q)$integral
            if (delta_q < 0) warning(paste("1", delta_q))
          } else {
            delta_q <- min(delta_max_q, cubature::cubintegrate(
              f = function(x) beta_l - t$tq(x), lower = t$tb_q, upper = t$te_q)$integral)
            if (delta_q < 0) warning(paste("2", delta_q))
          }
        } else if (o_q > beta_u) {
          if ((o_q + t$vty_q) < beta_u) {
            t1 <- if ((o_q + t$vty_q) > beta_l) t$te_q else t$te_q - ((beta_l - (o_q + t$vty_q)) / t$a_q)
            delta_q <- cubature::cubintegrate(
              f = function(x) min(beta_u, t$tq(x)) - beta_l, lower = t$tb_q, upper = t1)$integral
            if (delta_q < 0) warning(paste("3", delta_q))
          } else {
            delta_q <- 0
          }
        } else {
          if ((o_q + t$vty_q) > beta_u) {
            t1 <- ((o_q + t$vty_q) - beta_u) / t$a_q
            delta_q <- cubature::cubintegrate(
              f = function(x) min(beta_u, t$tq(x)) - beta_l, lower = t$tb_q, upper = t1)$integral
            if (delta_q < 0) warning(paste("4", delta_q))
          } else {
            t1 <- if ((o_q + t$vty_q) >= beta_l) t$te_q else t$te_q - ((beta_l - (o_q + t$vty_q)) / t$a_q)
            delta_q <- cubature::cubintegrate(
              f = function(x) t$tq(x) - beta_l, lower = t$tb_q, upper = t1)$integral
            if (delta_q < 0) warning(paste("5", delta_q))
          }
        }
        
        ratios <- c(ratios, min(1, max(0, delta_q / delta_max_q)))
      }
      
      `names<-`(c(-log(1 - mean(ratios))), srbtaw$getOutputNames())
    },
    
    funcTikhonov = function() {
      srbtaw <- private$srbtaw
      v <- srbtaw$getParams()[grep(pattern = "v$", x = names(srbtaw$getParams()))]
      vty <- srbtaw$getParams()[grep(pattern = "vty_", x = names(srbtaw$getParams()))]
      `names<-`(c(log(1 + sqrt(sum(c(v, vty)^2)))), srbtaw$getOutputNames())
    },
    
    funcAvgout = function() {
      srbtaw <- private$srbtaw
      mod <- srbtaw$getInstance(wpName = self$getWpName(), wcName = self$getWcName())
      
      qs <- private$intervals
      avgs <- c()
      res <- srbtaw$residuals(loss = self)
      for (q in qs) {
        t <- res[[q]]
        lb <- t$lambda_ymin_q
        ub <- t$lambda_ymax_q
        mu <- cubature::cubintegrate(f = mod$M, lower = 0, upper = 1)$integral
        avgs <- c(avgs, 2 * sqrt((1/2 * (ub - lb) - mu)^2))
      }
      
      `names<-`(c(-log(1 - mean(avgs))), srbtaw$getOutputNames())
    },
    
    get0Function = function() {
      if (private$use == "trapezoid") {
        self$funcTrapezoid
      } else if (private$use == "tikhonov") {
        self$funcTikhonov
      } else {
        self$funcAvgout
      }
    }
  )
)


residuals.srBTAW <- function(model, loss) {
  model$residuals(loss)
}




# srBTAW_Loss2Curves$debug("funcCorr")
# srBTAW$debug("fit")
# srBTAW_Loss_Rss$debug("get0Function")
# srBTAW$debug("residuals")
# SRBTW$debug("setAllParams")
# SRBTW_SubModel$debug("asTuple")
# SRBTW$debug("initialize")
# SRBTW_Loss$debug("compute")
# SRBTW_DataLoss$debug("compute")
# SRBTW_SingleObjectiveOptimization$debug("compute")


