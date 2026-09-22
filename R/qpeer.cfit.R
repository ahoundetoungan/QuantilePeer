
#' @title Predicting Instruments Using Cross-Fitting
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as \code{y ~ x1 + x2},
#' where `y` is the outcome and `x1` and `x2` are control variables, which may include contextual variables such as averages or quantiles among peers.
#' @param excluded.instruments An object of class \link[stats]{formula} indicating the excluded instruments. It should be specified as \code{~ z1 + z2},
#' where `z1` and `z2` are excluded instruments for the quantile peer outcomes.
#' @param A The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `A` must be a list of subnets, with the `g`-th element being an
#' \eqn{n_g \times n_g} adjacency matrix, where \eqn{n_g} is the number of nodes in the `g`-th subnet.
#' @param tau A numeric vector specifying the quantile levels at which new instruments should be computed.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see \link[stats]{quantile}).
#' @param data An optional data frame, list, or environment (or an object that can be coerced via \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `qpeer` is called.
#' @param tol A tolerance value used in the QR factorization to identify columns of the explanatory-variable and instrument matrices that ensure full rank
#' (see \link[base]{qr}).
#' @param structural A logical value indicating whether the reduced-form or structural specification should be estimated (see Details).
#' @param fixed.effects A logical value or character string indicating whether the model includes subnet fixed effects. Accepted values are `"no"` or `"FALSE"`
#' (no fixed effects), `"join"` or `TRUE` (the same fixed effects for isolated and non-isolated nodes within each subnet), and `"separate"`
#' (different fixed effects for isolated and non-isolated nodes within each subnet). Note that `"join"` fixed effects are not applicable to structural models;
#' therefore, `"join"` and `TRUE` are automatically converted to `"separate"`.
#' @param estimator A character string specifying the estimator used for the training model. Available options are:
#' \itemize{
#'   \item `lm` or `"ols"` for a linear model estimated using OLS;
#'   \item `"rf"` or `"random forest"` for random forests;
#'   \item `"xgboost"` for gradient boosting;
#'   \item `"lasso"` for Lasso regression.
#' }
#' @param checkrank A logical value indicating whether the instrument matrix should be checked for full rank. If the matrix is not full rank, redundant columns
#' are removed to obtain a full-rank matrix.
#' @param nthreads A strictly positive integer indicating the number of threads used when computing peer quantiles or performing the bootstrap.
#' @param print A logical value indicating whether prediction progress should be displayed.
#' @param nfold A strictly positive integer specifying the number of folds used for cross-fitting in the estimation of the prediction model.
#' @param control A list of additional arguments passed to the estimation method. For `"ols"`, see \link[stats]{lm}; for `"lasso"`, see
#' \link[glmnet]{cv.glmnet}; for `"xgboost"`, see \link[xgboost]{xgb.params} and \link[xgboost]{xgb.train}; and for `"rf"`, see \link[ranger]{ranger}.
#' @description
#' `qpeer.cfit` predicts instruments for quantile peer-effect models. It uses the variables in `X` and the excluded instruments specified in
#' `excluded.instruments` to predict quantiles of `y` among friends at the quantile levels specified in `tau`.
#'
#' These predictions can then be used as instruments to estimate the quantile peer-effect model via \code{\link{qpeer}}.
#' The resulting predictions are obtained using cross-fitting and are therefore exogenous.
#' @return A matrix whose columns correspond to the predicted quantiles of `y` among friends.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' ngr  <- 50  # Number of subnets
#' nvec <- rep(30, ngr)  # Size of subnets
#' n    <- sum(nvec)
#' 
#' ### Simulating Data
#' ## Network matrix
#' A <- lapply(1:ngr, function(g) {
#'   Ag <- matrix(rbinom(nvec[g]^2, 1, 0.3), nvec[g], nvec[g])
#'   diag(Ag) <- 0
#'   # Adding isolated nodes (important for the structural model)
#'   niso <- sample(0:nvec[g], 1, prob = (nvec[g] + 1):1 / sum((nvec[g] + 1):1))
#'   if (niso > 0) {
#'     Ag[sample(1:nvec[g], niso), ] <- 0
#'   }
#'   Ag
#' })
#' 
#' tau <- seq(0, 1, 1/3)
#' X   <- cbind(rnorm(n), rpois(n, 2))
#' lam <- c(0.2, 0.15, 0.1, 0.2)
#' gam <- c(2, -0.5, 1)
#' eps <- rnorm(n, 0, 0.4)
#' 
#' ## Generating `y`
#' y <- qpeer.sim(formula = ~ X, A = A, tau = tau, lambda = lam, 
#'                gamma = gam, epsilon = eps)$y
#' 
#' ### Estimation
#' ## Computing instruments
#' Z1 <- qpeer.inst(formula = ~ X, A = A, tau = seq(0, 1, 0.1), 
#'                  max.distance = 2, checkrank = TRUE)
#' Z1 <- Z1$instruments
#' Z2 <- qpeer.cfit(formula = y ~ X, excluded.instruments = ~ Z1,
#'                  A = A, tau = tau, estimator = "ols")
#' 
#' ## Reduced-form model 
#' rest <- qpeer(formula = y ~ X, excluded.instruments = ~ Z2, A = A, tau = tau)
#' summary(rest)
#' summary(rest, diagnostic = TRUE)  # Summary with diagnostics
#' }
#' @importFrom stats lm predict
#' @importFrom ranger ranger
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom glmnet cv.glmnet
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
qpeer.cfit <- function(formula, excluded.instruments, A, tau, type = 7, data, 
                       structural = FALSE, fixed.effects = FALSE, estimator = "ols", 
                       nfold = 2, tol = 1e-10, checkrank = FALSE, nthreads = 1, 
                       control = NULL, print = TRUE) {
  
  nthreads   <- fnthreads(nthreads = nthreads)
  # Quantiles
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  ntau       <- length(tau)
  
  # Estimator for the instrument
  if (tolower(estimator) %in% c("lin", "linear", "ols", "lm")) {
    estimator <- "OLS"
  }  else if (tolower(estimator) %in% c("rf", "r-f", "random forest", "random-forest", "randomforest")) {
    estimator <- "Random Forest"
  } else if (tolower(estimator) %in% c("xgboost")) {
    estimator <- "XGBoost"
  } else if (tolower(estimator) %in% c("lasso")) {
    estimator <- "LASSO"
  } else {
    stop("This estimator is not available.")
  }
  
  # Fixed effects
  if (is.character(fixed.effects[1])) fixed.effects <- tolower(fixed.effects)
  stopifnot(fixed.effects %in% c(FALSE, "no", TRUE, "yes", "join", "separate"))
  if (fixed.effects == FALSE) fixed.effects <- "no"
  if (fixed.effects == TRUE | fixed.effects == "yes") fixed.effects <- "join"
  if (structural & fixed.effects != "no") fixed.effects <- "separate"
  FEnum <- (0:2)[fixed.effects == c("no", "join", "separate")]
  
  # Network
  if (!is.list(A)) {
    A  <- list(A)
  }
  d        <- fnetwork(A = A)
  G        <- d$G
  GIs      <- d$GIs
  GnIs     <- d$GnIs
  nvec     <- d$nvec
  n        <- d$n
  igr      <- d$igr
  group    <- d$group
  groupidx <- d$groupidx
  lIs      <- d$lIs
  Is       <- d$Is
  lnIs     <- d$lnIs
  nIs      <- d$nIs
  ld       <- d$ld
  d        <- d$d
  
  # Assign fold
  gseed   <- as.integer(runif(1, 0, 1e9))
  if (nfold > G) {
    warning("The number of folds cannot be larger than the number of subnetworks. 'nfold' is set to the number of subnetworks.")
    nfold <- G
  }
  fold    <- fassignfold(subnetwork = rep(0:(G -1), nvec), nfold = nfold, seed = gseed)
  ifold   <- split(seq_along(fold), fold)
  
  
  # y and X
  formula    <- as.formula(formula)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = (fixed.effects != "no"), 
                                simulations = FALSE) 
  y          <- f.t.data$y
  X          <- f.t.data$X
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  xint       <- f.t.data$intercept
  qy         <- fQtauy(y = y, A = A, d = d, igroup = igr, group = group, groupidx = groupidx,
                       nvec = nvec, stau = tau, ngroup = G, n = n, ntau = ntau, type = type, 
                       nthreads = nthreads)
  qy         <- as.matrix(qy)
  
  # Instruments
  inst       <- as.formula(excluded.instruments); excluded.instruments <- inst
  if(length(inst) != 2) stop("The `excluded.instruments` argument must be in the format `~ z1 + z2 + ....`.")
  f.t.data   <- formula.to.data(formula = inst, data = data, fixed.effects = (fixed.effects != "no"), 
                                simulations = TRUE)
  ins        <- f.t.data$X
  zename     <- f.t.data$xname
  if (xint) {
    ins      <- ins[, zename != "(Intercept)"]
    zename   <- zename[zename != "(Intercept)"]
  } else {
    ins      <- ins
  }
  X          <- cbind(X, ins)
  
  # Demean fixed effect models
  # save original data
  if (fixed.effects != "no") {
    if (fixed.effects == "join") {
      qy     <- c(Demean(qy, igroup = igr, ngroup = G))
      X      <- Demean(X, igroup = igr, ngroup = G)
    } else {
      qy     <- c(Demean_separate(qy, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n))
      X      <- Demean_separate(X, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n)
    }
  }
  
  # Remove useless columns
  if (checkrank) {
    X        <- X[, fcheckrank(X = X, tol = tol) + 1, drop = FALSE]
  }
  X           <- as.data.frame(X)
  colnames(X) <- paste0("X", 1:ncol(X))
  
  # Prediction
  qyhat      <- matrix(NA, n, ntau)
  pb         <- NULL
  if (print) {
    pb       <- txtProgressBar(min = 0, max = ntau * nfold, style = 3)
  }
  ii         <- 0
  for (s in 1:ntau) {
    for (k in 1:nfold) { # For each fold
      # ytrain
      qy_notk <- qy[-ifold[[k]], s]
      # xtrain
      X_notk  <- X[-ifold[[k]],,drop = FALSE]
      # xpred
      X_k     <- X[ifold[[k]],,drop = FALSE]
      
      # Estimation
      qyhatk  <- NULL
      if (estimator == "OLS") {
        
        defname     <- c("formula", "data")
        ARG         <- c(list(formula = qy_notk ~ ., data = X_notk),
                         control[!(names(control) %in% defname)])
        model_train <- do.call(lm, ARG) 
        qyhatk      <- unname(predict(model_train, newdata = X_k))
        
      } else if (estimator == "Random Forest") {
        
        defname     <- c("formula", "data")
        ARG         <- c(list(formula = qy_notk ~ ., data = X_notk),
                         control[!(names(control) %in% defname)])
        model_train <- do.call(ranger, ARG)  
        qyhatk      <- predict(model_train, data = X_k)$predictions
        
      } else if (estimator == "XGBoost") {
        
        # Training data
        dtrain      <- xgb.DMatrix(data = as.matrix(X_notk),
                                   label = qy_notk)
        # prediction
        dpred       <- xgb.DMatrix(data = as.matrix(X_k))
        # parameters
        dpar        <- c(list(objective = "reg:squarederror"), 
                         control$params[names(control$params) != "objective"])
        
        # Training arguments
        defname     <- c("formula", "objective", "params")
        ARG         <- c(list(params = dpar, data = dtrain),
                         control[!(names(control) %in% defname)])
        ARG$nrounds <- if(is.null(ARG$nrounds)) 100 else ARG$nrounds
        # Training
        model_train <- do.call(xgb.train, ARG)
        # Prediction 
        qyhatk      <- predict(model_train, newdata = dpred)
        
      } else if (estimator == "LASSO") {
        
        defname     <- c("y", "x", "alpha")
        ARG         <- c(list(y = qy_notk, x = as.matrix(X_notk), alpha = 1),
                         control[!(names(control) %in% defname)])
        fitlasso    <- do.call(cv.glmnet, ARG) 
        qyhatk      <- as.numeric(predict(fitlasso, newx = as.matrix(X_k), 
                                          s = "lambda.min"))
      }
      
      qyhat[ifold[[k]], s] <- qyhatk
      ii            <- ii + 1
      if (print) {
        setTxtProgressBar(pb, ii)
      }
    }
  }
  
  # Prediction acurassy
  if (print) {
    RMSE    <- sapply(1:ntau, \(s) sqrt(mean((qyhat[,s] - qy[,s])^2)))
    MAE     <- sapply(1:ntau, \(s) sqrt(mean(abs(qyhat[,s] - qy[,s]))))
    names(RMSE) <- names(MAE) <- paste0("tau", 1:ntau)
    
    cat("\nRMSE:\n")
    print(RMSE)
    cat("\nMAE:\n")
    print(MAE)
  }
  
  qyhat
}

