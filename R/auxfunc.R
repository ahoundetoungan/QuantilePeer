#' @importFrom formula.tools env
#' @importFrom stats model.frame
#' @importFrom stats terms
#' @importFrom stats model.response
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
formula.to.data <- function(formula,
                            data, 
                            simulations   = FALSE,
                            fixed.effects = FALSE) {
  
  ## Extract data from the formula
  if (missing(data)) {
    data      <- env(formula)
  }
  formula     <- as.formula(formula)
  yname       <- NULL
  
  if (simulations) {
    if(length(formula) != 2) stop("`formula` is not correct. To simulate data, the expected format is ~ X1 + X2 + ....")
  } else {
    if(length(formula) != 3) stop("`formula` is not correct. For the estimation, the expected format is y ~ X1 + X2 + ....")
    yname     <- all.vars(formula)[1]
  }
  
  ## call model.frame()
  mf          <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y           <- model.response(mf, "numeric")
  X           <- model.matrix(terms(formula, data = data, rhs = 1), mf)
  xname       <- colnames(X)
  intercept   <- "(Intercept)" %in% xname
  if(fixed.effects & intercept){
    X         <- X[, xname != "(Intercept)", drop = FALSE]
    xname    <- xname[xname != "(Intercept)"]
    intercept <- FALSE
  }
  
  list("formula"   = formula, 
       "X"         = X, 
       "y"         = y,
       "intercept" = intercept,
       "yname"     = yname,
       "xname"     = xname)
}


fnetwork <- function(Glist) {
  dg       <- unlist(lapply(Glist, rowSums))
  Is       <- which(dg == 0) - 1
  nIs      <- which(dg != 0) - 1
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  list(dg = dg, M = M, nvec = nvec, n = n, igr = igr, Is = Is, nIs = nIs)
}

fcoef           <- function(Estimate, cov) {
  coef           <- cbind(Estimate, sqrt(diag(cov)), 0, 0)
  coef[,3]       <- coef[,1]/coef[,2]
  coef[,4]       <- 2*(1 - pnorm(abs(coef[,3])))
  colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  coef
}

fprintcoeft <- function(coef) {
  pval      <- coef[,ncol(coef)]
  pval_pt   <- sapply(pval, function(s){ifelse(s < 2e-16, "<2e-16", format(s, digit = 3))})
  refprob   <- c(0.001, 0.01, 0.05, 0.1)
  refstr    <- c("***",  "**", "*", ".", "")
  str       <- sapply(pval, function(s) refstr[1 + sum(s > refprob)])
  out       <- data.frame(coef[,-ncol(coef)], "P" = pval_pt, "S" = str); colnames(out) <- c(colnames(coef), "")
  print(out)
}

#' @importFrom stats pf
fdiagnostic  <- function(estim, sar = FALSE) {
  ntau     <- estim$model.info$ntau
  FE       <- estim$model.info$fixed.effects
  struc    <- estim$model.info$structural
  y        <- as.matrix(estim$data$y)
  endo     <- as.matrix(estim$data[[ifelse(sar, "Gy", "qy")]])
  X        <- estim$data$X
  ins      <- estim$data$instruments
  
  idX1     <- estim$model.info$idX1
  idX2     <- estim$model.info$idX2
  
  if (struc) {
    xb     <- X[, idX1, drop = FALSE] %*% estim$gmm$Estimate[idX1 + ntau + 1]
    if (length(idX2) == 0) {
      X    <- as.matrix(xb)
    } else {
      X    <- cbind(xb, X[, idX2, drop = FALSE])
    }
    ins    <- cbind(xb, ins)
  }
  
  M        <- estim$model.info$ngroup
  nvec     <- estim$model.info$nvec
  n        <- estim$model.info$n
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  Is       <- estim$data$isolated - 1
  nIs      <- estim$data$non.isolated - 1
  
  if (FE == "join") {
    y      <- demean(y, igroup = igr, ngroup = M)
    endo   <- demean(endo, igroup = igr, ngroup = M)
    X      <- demean(X, igroup = igr, ngroup = M)
    ins    <- demean(ins, igroup = igr, ngroup = M)
  } else if(FE == "separate") {
    y      <- demean_separate(y, igroup = igr, Is = Is, ngroup = M, n = n)
    endo   <- demean_separate(endo, igroup = igr, Is = Is, ngroup = M, n = n)
    X      <- demean_separate(X, igroup = igr, Is = Is, ngroup = M, n = n)
    ins    <- demean_separate(ins, igroup = igr, Is = Is, ngroup = M, n = n)
  }
  
  if (struc) {
    y      <- y[nIs + 1,, drop = FALSE]
    endo   <- endo[nIs + 1,, drop = FALSE]
    X      <- X[nIs + 1,, drop = FALSE]
    ins    <- ins[nIs + 1,, drop = FALSE]
  }
  
  ## Weak instrument test
  tpF      <- fFstat(y = endo, Xc = X, Xu = ins)
  
  ## Endogeneity test
  tpend    <- fFstat(y = y, Xc = cbind(endo, X), Xu = cbind(tpF$ru, endo, X))
  
  out      <- cbind(df1        = c(rep(tpF$df1, length(tpF$F)) , tpend$df1, estim$gmm$Jtest["df"]),
                    df2        = c(rep(tpF$df2, length(tpF$F)), tpend$df2, NA),
                    statistic  = c(tpF$F, tpend$F, estim$gmm$Jtest["statistic"]),
                    "p-value"  = estim$gmm$Jtest["p-value"])
  out[-nrow(out), 4] <- pf(out[-nrow(out), 3], out[-nrow(out), 1], out[-nrow(out), 2], lower.tail = FALSE)
  rn            <- "Weak instruments"
  if (ntau > 1 & !sar) {
    rn          <- paste0(rn, " (", colnames(estim$data$qy), ")")
  }
  rn            <- c(rn, "Wu-Hausman", "Hansen's J-test")
  rownames(out) <- rn
  out
}