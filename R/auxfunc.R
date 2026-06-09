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
fnetwork   <- function(A, binary = FALSE) {
  G        <- length(A)
  nvec     <- unlist(lapply(A, nrow))
  n        <- sum(nvec)
  ncs      <- c(0, cumsum(nvec))
  group    <- rep(0:(G - 1), nvec)
  groupidx <- unlist(lapply(1:G, \(g) 0:(nvec[g] - 1)))
  
  ld       <- lapply(A, rowSums)
  d        <- unlist(ld)
  GIs      <- sum(sapply(ld, function(g) any(g == 0)))
  GnIs     <- sum(sapply(ld, function(g) any(g != 0)))
  lIs      <- lapply(1:G, function(g) which(ld[[g]] == 0) - 1 + ncs[g])
  Is       <- unlist(lIs)
  lnIs     <- lapply(1:G, function(g) which(ld[[g]] != 0) - 1 + ncs[g])
  nIs      <- unlist(lnIs)
  
  list(d = d, ld = ld, G = G, nvec = nvec, n = n, igr = ncs, group = group,
       groupidx = groupidx, Is = Is, nIs = nIs, lIs = lIs, lnIs = lnIs, GIs = GIs, 
       GnIs = GnIs)
}

fdrop <- function(drop, ld, nvec, G, lIs, lnIs, y, X, qy, ins) {
  n        <- sum(nvec)
  if (any(!(drop %in% 0:1) | !is.finite(drop))) {
    stop("`drop` must be a binary (0/1) variable.")
  }
  if (length(drop) != n) {
    stop("`drop` must be a vector of length n.")
  }
  ncs      <- c(0, cumsum(nvec))
  olIs     <- lapply(1:G, function(g) ld[[g]] == 0)
  oIs      <- unlist(olIs)
  lkeep    <- lapply(1:G, function(g) drop[(ncs[g] + 1):ncs[g + 1]] != 1)
  keep     <- unlist(lkeep)
  gkeep    <- sapply(1:G, function(g) sum(lkeep[[g]]) >= 1) # Groups I keep
  ld       <- lapply(1:G, function(g) ld[[g]][lkeep[[g]]])[gkeep]
  d        <- unlist(ld)
  G        <- length(ld)
  nvec     <- sapply(ld, length)
  n        <- sum(nvec)
  ncs      <- c(0, cumsum(nvec))
  group    <- rep(0:(G - 1), nvec)
  groupidx <- unlist(lapply(1:G, \(g) 0:(nvec[g] - 1)))
  GIs      <- sum(sapply(ld, function(g) any(g == 0)))
  GnIs     <- sum(sapply(ld, function(g) any(g != 0)))
  lIs      <- lapply(1:G, function(g) which(ld[[g]] == 0) - 1 + ncs[g])
  Is       <- unlist(lIs)
  lnIs     <- lapply(1:G, function(g) which(ld[[g]] != 0) - 1 + ncs[g])
  nIs      <- unlist(lnIs)
  y        <- y[keep]
  X        <- X[keep, , drop = FALSE]
  qy       <- qy[keep, , drop = FALSE]
  ins      <- ins[keep, , drop = FALSE]
  list(d = d, ld = ld, G = G, nvec = nvec, n = n, igr = ncs, group = group,
       groupidx = groupidx, Is = Is, nIs = nIs, lIs = lIs, lnIs = lnIs, GIs = GIs, 
       GnIs = GnIs, y = y, X = X, qy = qy, ins = ins)
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
  G        <- object$model.info$ngroup
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
  d        <- object$data$degree
  index    <- which(!(colnames(ins) %in% colnames(X))) - 1
  lIs      <- lapply(object$data$isolated, \(g) g - 1)
  lnIs     <- lapply(object$data$non.isolated, \(g) g - 1)
  Is       <- unlist(lIs)
  nIs      <- unlist(lnIs)
  nvc      <- sapply(lnIs, length)
  theta    <- object$gmm$Estimate
  
  ins      <- fdatadiagnostic(y = y, endo = endo, X = X, ins = ins, theta = theta, idX1 = idX1, idX2 = idX2, igroup = ncs, 
                              nIs = nIs, LIs = lIs, LnIs = lnIs, n = n, ngroup = G, ntau = ntau, struc = struc, FE = FE)
  y        <- ins$y
  endo     <- ins$endo
  X        <- ins$X
  ins      <- ins$ins
  
  if (struc) {
    nvc    <- nvc[nvc > 0]
    G      <- length(nvc)
    ncs    <- c(0, cumsum(nvc))
  }
  
  out      <- NULL
  cvKP     <- NULL
  if (object$model.info$estimator %in% c("JIVE", "JIVE2")) {
    ## Weak instrument test
    tpF    <- fFstat(y = endo, X = ins, index = index, igroup = ncs, ngroup = G, HAC = HACnum)
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
    tpKP   <- fKPstat_boot(qy = endo, Z = ins, index = index, igroup = ncs,
                           LnIs = lnIs, LIs = lIs, ngroup = G, boot = boot, 
                           nthreads = nthreads, seed = seed, print = print)
    
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
    tpF    <- fFstat(y = endo, X = ins, index = index, igroup = ncs, ngroup = G, HAC = HACnum)
    tpKP   <- fKPstat(qy = endo, Z = ins, index = index, igroup = ncs, HAC = HACnum)
    
    ## Endogeneity test
    tpend  <- fFstat(y = y, X = cbind(tpF$ru, endo, X), index = (0:(ntau - 1)), igroup = ncs, ngroup = G, HAC = HACnum)
    
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
fCESdatainit  <- function (y, z, A, nvec, G, ld, lIs, lnIs, drop) {
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
  friendindex <- lapply(1:G, function(g) {
    lapply(1:nvec[g], function(i) {
      which(A[[g]][i,] > 0) - 1
    })})
  frzeroy     <- as.integer(unlist(lapply(1:G, function(g){
    lapply(1:nvec[g], function(i){
      any(y[friendindex[[g]][[i]] + ncs[g] + 1] <= 0)
    })})))
  frzeroz     <- as.integer(unlist(lapply(1:G, function(g){
    lapply(1:nvec[g], function(i){
      any(z[friendindex[[g]][[i]] + ncs[g] + 1] <= 0)
    })})))
  lsel        <- lapply(1:G, function(g) drop[(ncs[g] + 1):ncs[g + 1]] != 1)
  
  # Max and Min of friend y and z
  yFmax       <- unlist(lapply(1:G, function(g){
    lapply(1:nvec[g], function(i){
      ifelse(ld[[g]][i] > 0, max(y[friendindex[[g]][[i]] + ncs[g] + 1]), NA)
    })
  }))
  yFmin       <- unlist(lapply(1:G, function(g){
    lapply(1:nvec[g], function(i){
      ifelse(ld[[g]][i] > 0, min(y[friendindex[[g]][[i]] + ncs[g] + 1]), NA)
    })
  }))
  zFmax       <- unlist(lapply(1:G, function(g){
    lapply(1:nvec[g], function(i){
      ifelse(ld[[g]][i] > 0, max(z[friendindex[[g]][[i]] + ncs[g] + 1]), NA)
    })
  }))
  zFmin       <- unlist(lapply(1:G, function(g){
    lapply(1:nvec[g], function(i){
      ifelse(ld[[g]][i] > 0, min(z[friendindex[[g]][[i]] + ncs[g] + 1]), NA)
    })
  }))
  
  # In selection variables
  ld          <- lapply(1:G, function(g) ld[[g]][lsel[[g]]])
  lIs         <- lapply(1:G, function(g) lIs[[g]][lsel[[g]][lIs[[g]] - ncs[g] + 1]])
  lnIs        <- lapply(1:G, function(g) lnIs[[g]][lsel[[g]][lnIs[[g]] - ncs[g] + 1]])
  Is          <- unlist(lIs)
  nIs         <- unlist(lnIs)
  
  # In selection variables if empty groups are removed
  keepg       <- sapply(1:G, function(g) length(ld[[g]]) > 0)
  ld          <- ld[keepg]
  G           <- length(ld)
  GIs         <- sum(sapply(lIs, function(g) length(g) > 0))
  GnIs        <- sum(sapply(lnIs, function(g) length(g) > 0))
  
  list(friendindex = friendindex, frzeroy = frzeroy, frzeroz = frzeroz, G = G, GIs = GIs, GnIs = GnIs,
       ld = ld, d = unlist(ld), lIs = lIs, Is = Is, lnIs = lnIs, nIs = nIs, hasIso = (length(Is) > 0),
       yFmax = yFmax, yFmin = yFmin, zFmax = zFmax, zFmin = zFmin)
}
