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
    if(length(formula) != 2) stop("The `formula` argument is invalid. For data simulation, the expected format is `~ X1 + X2 + ...`.")
  } else {
    if(length(formula) != 3) stop("The `formula` argument is invalid. For estimation, the expected format is `y ~ X1 + X2 + ...`.")
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

#' @importFrom utils head
#' @importFrom utils tail
fnetwork   <- function(Glist, isol = NULL) {
  # Isol is true isolated than can be removed. But this argument is no longer used. 
  # See the more general argument which is now drop
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  ncs      <- c(0, cumsum(nvec))
  group    <- rep(0:(M - 1), nvec)
  groupidx <- unlist(lapply(1:M, \(m) 0:(nvec[m] - 1)))
  
  ldg      <- lapply(Glist, rowSums)
  dg       <- unlist(ldg)
  MIs      <- NULL
  MnIs     <- NULL
  lIs      <- NULL
  Is       <- NULL
  lnIs     <- NULL
  nIs      <- NULL
  if (length(isol) == 0) {
    MIs    <- sum(sapply(ldg, function(s) any(s == 0)))
    MnIs   <- sum(sapply(ldg, function(s) any(s != 0)))
    lIs    <- lapply(1:M, function(m) which(ldg[[m]] == 0) - 1 + ncs[m])
    Is     <- unlist(lIs)
    lnIs   <- lapply(1:M, function(m) which(ldg[[m]] != 0) - 1 + ncs[m])
    nIs    <- unlist(lnIs)
  } else {
    if (any(!(isol %in% 0:1) | !is.finite(isol))) {
      stop("`isolated` must be a binary (0/1) variable.")
    }
    if (length(isol) != n) {
      stop("`isolated` must be a vector of length n.")
    }
    lisol  <- lapply(1:M, function(m) isol[(ncs[m] + 1):ncs[m + 1]])
    MIs    <- sum(sapply(lisol, function(s) any(s == 1)))
    MnIs   <- sum(sapply(lisol, function(s) any(s != 1)))
    lIs    <- lapply(1:M, function(m) which(lisol[[m]] == 1) - 1 + ncs[m])
    Is     <- unlist(lIs)
    lnIs   <- lapply(1:M, function(m) which(lisol[[m]] != 1) - 1 + ncs[m])
    nIs    <- unlist(lnIs)
  }
  
  list(dg = dg, ldg = ldg, M = M, nvec = nvec, n = n, igr = ncs, group = group,
       groupidx = groupidx, Is = Is, nIs = nIs, lIs = lIs, lnIs = lnIs, MIs = MIs, 
       MnIs = MnIs)
}

fdrop <- function(drop, ldg, nvec, M, lIs, lnIs, y, X, qy, ins) {
  n        <- sum(nvec)
  if (any(!(drop %in% 0:1) | !is.finite(drop))) {
    stop("`drop` must be a binary (0/1) variable.")
  }
  if (length(drop) != n) {
    stop("`drop` must be a vector of length n.")
  }
  ncs      <- c(0, cumsum(nvec))
  olIs     <- lapply(1:M, function(m) ldg[[m]] == 0)
  oIs      <- unlist(olIs)
  lkeep    <- lapply(1:M, function(m) drop[(ncs[m] + 1):ncs[m + 1]] != 1)
  keep     <- unlist(lkeep)
  gkeep    <- sapply(1:M, function(m) sum(lkeep[[m]]) >= 1) # Groups I keep
  ldg      <- lapply(1:M, function(m) ldg[[m]][lkeep[[m]]])[gkeep]
  dg       <- unlist(ldg)
  M        <- length(ldg)
  nvec     <- sapply(ldg, length)
  n        <- sum(nvec)
  ncs      <- c(0, cumsum(nvec))
  group    <- rep(0:(M - 1), nvec)
  groupidx <- unlist(lapply(1:M, \(m) 0:(nvec[m] - 1)))
  MIs      <- sum(sapply(ldg, function(s) any(s == 0)))
  MnIs     <- sum(sapply(ldg, function(s) any(s != 0)))
  lIs      <- lapply(1:M, function(m) which(ldg[[m]] == 0) - 1 + ncs[m])
  Is       <- unlist(lIs)
  lnIs     <- lapply(1:M, function(m) which(ldg[[m]] != 0) - 1 + ncs[m])
  nIs      <- unlist(lnIs)
  y        <- y[keep]
  X        <- X[keep, , drop = FALSE]
  qy       <- qy[keep, , drop = FALSE]
  ins      <- ins[keep, , drop = FALSE]
  list(dg = dg, ldg = ldg, M = M, nvec = nvec, n = n, igr = ncs, group = group,
       groupidx = groupidx, Is = Is, nIs = nIs, lIs = lIs, lnIs = lnIs, MIs = MIs, 
       MnIs = MnIs, y = y, X = X, qy = qy, ins = ins)
}

fcheckrank <- function(X, tol = 1e-10) {
  which(fcheckrankEigen(X, tol)) - 1
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
  pval_pt   <- sapply(pval, function(s){ifelse(is.na(s), "NA", ifelse(s < 2e-16, "<2e-16", format(s, digit = 4)))})
  refprob   <- c(0.001, 0.01, 0.05, 0.1)
  refstr    <- c("***",  "**", "*", ".", "")
  str       <- sapply(pval, function(s) ifelse(is.na(s), "", refstr[1 + sum(s > refprob)]))
  out       <- data.frame(coef[,-ncol(coef), drop = FALSE], "P" = pval_pt, "S" = str); 
  colnames(out) <- c(colnames(coef), "")
  print(out)
}

#' @importFrom stats pf
fdiagnostic  <- function(object, nendo, seed, boot, nthreads, print) {
  FE       <- object$model.info$fixed.effects
  struc    <- object$model.info$structural
  M        <- object$model.info$ngroup
  nvec     <- object$model.info$nvec
  n        <- object$model.info$n
  HAC      <- object$model.info$HAC
  HACnum   <- (0:3)[HAC == c("iid", "hetero", "cluster", "cluster-bootstrap")]
  ncs      <- c(0, cumsum(nvec))
  idX1     <- object$model.info$idXiso - 1
  idX2     <- object$model.info$idXniso - 1
  y        <- as.matrix(object$data$y)
  endo     <- as.matrix(object$data[[nendo]])
  cnendo   <- colnames(endo)
  ntau     <- ncol(endo)
  X        <- object$data$X
  ins      <- object$data$instruments
  dg       <- object$data$degree
  index    <- which(!(colnames(ins) %in% colnames(X))) - 1
  lIs      <- lapply(object$data$isolated, \(s) s - 1)
  lnIs     <- lapply(object$data$non.isolated, \(s) s - 1)
  Is       <- unlist(lIs)
  nIs      <- unlist(lnIs)
  nvc      <- sapply(lnIs, length)
  theta    <- object$gmm$Estimate
  
  ins      <- fdatadiagnostic(y = y, endo = endo, X = X, ins = ins, theta = theta, idX1 = idX1, idX2 = idX2, igroup = ncs, 
                              nIs = nIs, LIs = lIs, LnIs = lnIs, n = n, ngroup = M, ntau = ntau, struc = struc, FE = FE)
  y        <- ins$y
  endo     <- ins$endo
  X        <- ins$X
  ins      <- ins$ins
  
  if (struc) {
    nvc    <- nvc[nvc > 0]
    M      <- length(nvc)
    ncs    <- c(0, cumsum(nvc))
  }
  
  out      <- NULL
  cvKP     <- NULL
  if (object$model.info$estimator %in% c("JIVE", "JIVE2")) {
    ## Weak instrument test
    tpF    <- fFstat(y = endo, X = ins, index = index, igroup = ncs, ngroup = M, HAC = HACnum)
    tpKP   <- fKPstat(qy = endo, Z = ins, index = index, igroup = ncs, HAC = HACnum)
    ## Endogeneity test
    # Not implemented
    
    out    <- cbind(df1        = c(rep(tpF$df1, ntau), tpKP$df, object$gmm$Jtest["df"]),
                    df2        = c(rep(tpF$df2, ntau), NA, NA),
                    statistic  = c(tpF$F, tpKP$stat, object$gmm$Jtest["statistic"]),
                    "p-value"  = object$gmm$Jtest["p-value"])
    out[-ntau - 2, 4] <- pf(out[-ntau - 2, 3], out[-ntau - 2, 1], out[-ntau - 2, 2], lower.tail = FALSE)
    out[ntau + 1, 4]  <- pchisq(out[ntau + 1, 3], out[ntau + 1, 1], lower.tail = FALSE)
    rn            <- "Weak many instruments"
    if (ntau > 1) {
      rn          <- paste0(rn, " (", cnendo, ")")
    }
    rn            <- c(rn, "Kleibergen-Paap rk Wald", "Hansen's J-test")
    rownames(out) <- rn
  } else if (HACnum == 3) {
    ## Weak instrument test
    tpKP   <- fKPstat_boot(qy = endo, Z = ins, index = index, LnIs = lnIs, 
                           LIs = lIs, ngroup = M, boot = boot, nthreads = nthreads,
                           seed = seed, print = print)
    
    ## Endogeneity test
    # not implemented!
    
    out    <- cbind(df1        = c(tpKP$df, object$gmm$Jtest["df"]),
                    df2        = c(NA, NA),
                    statistic  = c(tpKP$stat, object$gmm$Jtest["statistic"]),
                    "p-value"  = object$gmm$Jtest["p-value"])
    out[1, 4]     <- pchisq(out[1, 3], out[1, 1], lower.tail = FALSE)
    rn            <- c("Kleibergen-Paap rk Wald", "Hansen J")
    rownames(out) <- rn
    
  } else {
    ## Weak instrument test
    tpF    <- fFstat(y = endo, X = ins, index = index, igroup = ncs, ngroup = M, HAC = HACnum)
    tpKP   <- fKPstat(qy = endo, Z = ins, index = index, igroup = ncs, HAC = HACnum)
    
    ## Endogeneity test
    tpend  <- fFstat(y = y, X = cbind(tpF$ru, endo, X), index = (0:(ntau - 1)), igroup = ncs, ngroup = M, HAC = HACnum)
    
    out    <- cbind(df1        = c(rep(tpF$df1, ntau), tpKP$df, tpend$df1, object$gmm$Jtest["df"]),
                    df2        = c(rep(tpF$df2, ntau), NA, tpend$df2, NA),
                    statistic  = c(tpF$F, tpKP$stat, tpend$F, object$gmm$Jtest["statistic"]),
                    "p-value"  = object$gmm$Jtest["p-value"])
    out[-ntau - 3, 4] <- pf(out[-ntau - 3, 3], out[-ntau - 3, 1], out[-ntau - 3, 2], lower.tail = FALSE)
    out[ntau + 1, 4]  <- pchisq(out[ntau + 1, 3], out[ntau + 1, 1], lower.tail = FALSE)
    rn            <- "Weak instruments"
    if (ntau > 1) {
      rn          <- paste0(rn, " (", cnendo, ")")
    }
    rn            <- c(rn, "Kleibergen-Paap rk Wald", "Wu-Hausman", "Hansen J")
    rownames(out) <- rn
  }
  list(diag = out, cvKP = cvKP)
}

## Create data to start optimization
fCESdatainit  <- function (y, z, G, nvec, M, ldg, lIs, lnIs, drop) {
  n           <- sum(nvec)
  if (length(drop) == 0) {
    drop      <- rep(0, n)
  }
  if (any(!(drop %in% 0:1) | !is.finite(drop))) {
    stop("`drop` must be a binary (0/1) variable.")
  }
  if (length(drop) != n) {
    stop("`drop` must be a vector of length n.")
  }
  ncs         <- c(0, cumsum(nvec))
  friendindex <- lapply(1:M, function(m) {
    lapply(1:nvec[m], function(s) {
      which(G[[m]][s,] > 0) - 1
    })})
  frzeroy     <- as.integer(unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      any(y[friendindex[[m]][[s]] + ncs[m] + 1] <= 0)
    })})))
  frzeroz     <- as.integer(unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      any(z[friendindex[[m]][[s]] + ncs[m] + 1] <= 0)
    })})))
  lsel        <- lapply(1:M, function(m) drop[(ncs[m] + 1):ncs[m + 1]] != 1)
  
  # Max and Min of friend y and z
  yFmax       <- unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, max(y[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  yFmin       <- unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, min(y[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  zFmax       <- unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, max(z[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  zFmin       <- unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, min(z[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  
  # In selection variables
  ldg         <- lapply(1:M, function(m) ldg[[m]][lsel[[m]]])
  lIs         <- lapply(1:M, function(m) lIs[[m]][lsel[[m]][lIs[[m]] - ncs[m] + 1]])
  lnIs        <- lapply(1:M, function(m) lnIs[[m]][lsel[[m]][lnIs[[m]] - ncs[m] + 1]])
  Is          <- unlist(lIs)
  nIs         <- unlist(lnIs)
  
  # In selection variables if empty groups are removed
  keepg       <- sapply(1:M, function(m) length(ldg[[m]]) > 0)
  ldg         <- ldg[keepg]
  M           <- length(ldg)
  MIs         <- sum(sapply(lIs, function(s) length(s) > 0))
  MnIs        <- sum(sapply(lnIs, function(s) length(s) > 0))
  
  list(friendindex = friendindex, frzeroy = frzeroy, frzeroz = frzeroz, M = M, MIs = MIs, MnIs = MnIs,
       ldg = ldg, dg = unlist(ldg), lIs = lIs, Is = Is, lnIs = lnIs, nIs = nIs, hasIso = (length(Is) > 0),
       yFmax = yFmax, yFmin = yFmin, zFmax = zFmax, zFmin = zFmin)
}
