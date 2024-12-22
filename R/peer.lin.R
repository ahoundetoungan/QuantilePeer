#' @rdname qpeer
#' @export
linpeer <- function(formula, excluded.instruments, Glist, data, gmm.weight = "IV", 
                    structural = FALSE, fixed.effects = FALSE, HAC = "iid", checkrank = FALSE, tol = 1e-10){
  # GMM weight
  gmm.weight <- tolower(gmm.weight)
  stopifnot(gmm.weight %in% c("iv", "optimal", "ident"))
  
  # Variance structure
  HAC        <- tolower(HAC[1])
  stopifnot(HAC %in% c("iid", "hetero", "cluster"))
  HACnum     <- (0:2)[HAC == c("iid", "hetero", "cluster")]
  
  # Fixed effects
  if (is.character(fixed.effects[1])) fixed.effects <- tolower(fixed.effects)
  stopifnot(fixed.effects %in% c(FALSE, "no", TRUE, "yes", "join", "separate"))
  if (fixed.effects == FALSE) fixed.effects <- "no"
  if (fixed.effects == TRUE | fixed.effects == "yes") fixed.effects <- "join"
  if (structural & fixed.effects != "no") fixed.effects <- "separate"
  FEnum = (0:2)[fixed.effects == c("no", "join", "separate")]
  
  
  # Network
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  dg       <- fnetwork(Glist = Glist)
  M        <- dg$M
  MIs      <- dg$MIs
  MnIs     <- dg$MnIs
  nvec     <- dg$nvec
  ncum     <- c(0, cumsum(nvec))
  n        <- dg$n
  igr      <- dg$igr
  Is       <- dg$Is
  nIs      <- dg$nIs
  dg       <- dg$dg
  
  # Data
  # y and X
  formula    <- as.formula(formula)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = (fixed.effects != "no"), 
                                simulations = FALSE) 
  y          <- f.t.data$y
  X          <- f.t.data$X
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  xint       <- f.t.data$intercept
  Gy         <- as.matrix(unlist(lapply(1:M, function(m) Glist[[m]] %*% y[(ncum[m] + 1):ncum[m + 1]])))
  
  # Instruments
  inst       <- as.formula(excluded.instruments)
  if(length(inst) != 2) stop("The `excluded.instruments` argument must be in the format ~ `z1 + z2 + ....`.")
  f.t.data   <- formula.to.data(formula = inst, data = data, fixed.effects = (fixed.effects != "no"), 
                                simulations = TRUE)
  ins        <- f.t.data$X
  zename     <- f.t.data$xname
  if (xint) {
    ins      <- ins[, zename != "(Intercept)"]
  } else {
    ins      <- ins
  }
  
  # Demean fixed effect models
  # save original data
  y0         <- y
  Gy0        <- Gy
  X0         <- X
  ins0       <- ins
  if (fixed.effects != "no") {
    if (fixed.effects == "join") {
      y      <- c(demean(as.matrix(y), igroup = igr, ngroup = M))
      Gy     <- demean(Gy, igroup = igr, ngroup = M)
      X      <- demean(X, igroup = igr, ngroup = M)
      ins    <- demean(ins, igroup = igr, ngroup = M)
    } else {
      y      <- c(demean_separate(as.matrix(y), igroup = igr, Is = Is, ngroup = M, n = n))
      Gy     <- demean_separate(Gy, igroup = igr, Is = Is, ngroup = M, n = n)
      X      <- demean_separate(X, igroup = igr, Is = Is, ngroup = M, n = n)
      ins    <- demean_separate(ins, igroup = igr, Is = Is, ngroup = M, n = n)
    }
    colnames(X)   <- xname
    colnames(ins) <- zename
  }
  
  # Remove useless columns
  idX1       <- 0:(ncol(X) - 1)
  tlm        <- idX1
  if (structural) {
    idX1     <- c(fcheckrank(X = X[Is + 1,], tol = tol))
    tlm      <- c(fcheckrank(X = X[nIs + 1,], tol = tol))
  } else {
    tlm      <- c(fcheckrank(X = X, tol = tol))
  }
  idX2       <- which(!(tlm %in% idX1)) - 1 
  idX1       <- which(tlm %in% idX1) - 1 
  X          <- X[, tlm + 1, drop = FALSE]
  X0         <- X0[, tlm + 1, drop = FALSE]
  xname      <- xname[tlm + 1]
  Kx         <- ncol(X)
  if (structural) {
    if (length(c(fcheckrank(X = cbind(Gy, X)[nIs + 1,], tol = tol))) != (1 + Kx)) stop("The design matrix is not full rank.")
  } else {
    if (length(c(fcheckrank(X = cbind(Gy, X), tol = tol))) != (1 + Kx)) stop("The design matrix is not full rank.")
  }
  

  if (structural) {
    ins      <- cbind(X[, idX2 + 1], ins)
    ins0     <- cbind(X0[, idX2 + 1], ins0)
    zename   <- c(xname[idX2 + 1], zename)
    if (checkrank) {
      tlm    <- c(fcheckrank(X = ins[nIs + 1,], tol = tol))
    }
  } else {
    ins      <- cbind(X, ins)
    ins0     <- cbind(X0, ins0)
    zename   <- c(xname, zename)
    if (checkrank) {
      tlm    <- c(fcheckrank(X = ins, tol = tol))
    }
  }
  if (checkrank) {
    ins      <- ins[, tlm + 1, drop = FALSE]
    ins0     <- ins0[, tlm + 1, drop = FALSE]
    zename   <- zename[tlm + 1]
  }
  Kins       <- ncol(ins)
  
  # GMM
  GMMe       <- list()
  iv         <- (gmm.weight %in% c("iv", "optimal"))
  estname    <- NULL
  Kest       <- NULL
  if (structural) {
    Kx1      <- length(idX1)
    Kx2      <- length(idX2)
    if (Kins < Kx2 + 1) stop("Insufficient number of instruments: the model is not identified.")
    Kest1    <- ifelse(FEnum == 0, Kx1, Kx1 + MIs)
    Kest2    <- ifelse(FEnum == 0, Kx2 + 2, Kx2 + 1 + MnIs)
    if (length(Is) <= Kest1) stop("Insufficient number of isolated nodes for estimating the structural model.")
    if (length(nIs) <= Kest2) stop("Insufficient number of nonisolated nodes for estimating the structural model.")
    Kest     <- Kest1 + Kest2
    if (HACnum == 2 && (Kx1 >= MIs || Kins + 1 >= MnIs)) {
      warning("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. HAC is set to 'hetero'.")
      HACnum <- 1
      HAC    <- "hetero"
    }
    estname  <- c(paste0(c("G(conformity):", "G(total):"), yname), xname)
    
    # Estimation
    GMMe     <- fgmm_struc(y = y, X = X, qy = Gy, ins = ins, W1 = diag(Kx1), W2 = diag(Kins + 1), idX1 = idX1, 
                           idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igroup = igr, nIs = nIs, Is = Is, ngroup = M, 
                           Kins = Kins, Kx = Kx, ntau = 1, Kest = Kest, n = n, HAC = HACnum, iv = iv)
    if (gmm.weight == "optimal" & HACnum != 0) {
      GMMe   <- fgmm_struc(y = y, X = X, qy = Gy, ins = ins, W1 = solve(GMMe$VF1), W2 = solve(GMMe$VF2), idX1 = idX1, 
                           idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igroup = igr, nIs = nIs, Is = Is, ngroup = M, Kins = Kins, 
                           Kx = Kx, ntau = 1, Kest = Kest, n = n, HAC = HACnum, iv = FALSE)
    }
    Vpa        <- fStructParam(param = c(GMMe$beta, GMMe$lambda), covp = GMMe$Vpa, idX1 = idX1, idX2 = idX2, 
                               ntau = 1, Kx = Kx, Kx1 = Kx1, Kx2 = Kx2)
    GMMe$parms <- c(Vpa$theta)
    GMMe$Vpa   <- Vpa$Vpa
  } else {
    if (Kins < Kx + 1) stop("Insufficient number of instruments: the model is not identified.")
    Kest     <- ifelse(FEnum == 0, Kx + 1, ifelse(FEnum == 1, Kx + 1 + M, Kx + 1 + MIs + MnIs))
    if (n <= Kest) stop("Insufficient number of observations.")
    if (HACnum == 2 && Kins >= M) {
      warning("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. HAC is set to 'hetero'.")
      HACnum <- 1
      HAC    <- "hetero"
    }
    estname  <- c(paste0("G:", yname), xname)
    V        <- cbind(Gy, X)
    
    # Estimation
    GMMe     <- fgmm_red(y = y, V = V, ins = ins, W = diag(Kins), igroup = igr, ngroup = M, 
                         Kx = Kx, Kins = Kins, ntau = 1, Kest = Kest, n = n, HAC = HACnum, iv = iv)
    if (gmm.weight == "optimal" & HACnum != 0) {
      GMMe   <- fgmm_red(y = y, V = V, ins = ins, W = solve(GMMe$VZe), igroup = igr, ngroup = M, Kx = Kx, 
                         Kins = Kins, ntau = 1, Kest = Kest, n = n, HAC = HACnum, iv = FALSE)
    }
  }
  fv         <- c(GMMe$yhat)
  res        <- y - fv
  rs         <- sum((fv - mean(fv))^2)/sum((y - mean(y))^2)
  ars        <- 1 - (1 -rs)*(n - 1)/(n - length(c(GMMe$parms)))
  sigma      <- sqrt(GMMe$sigma2);  if(is.na(sigma)) sigma = NULL
  GMMe       <- list(Estimate = c(GMMe$parms), cov = GMMe$Vpa, sigma = sigma, fitted.values = fv, 
                     residuals = res, rsquared = rs, adjusted.rsquared = ars, df.residual = n - Kest,
                     Jtest = c("statistic" = GMMe$Overident, "df" = as.integer(GMMe$df)))
  names(GMMe$Estimate)  <- estname
  colnames(GMMe$cov)    <- estname
  rownames(GMMe$cov)    <- estname
  GMMe$Jtest["p-value"] <- ifelse(GMMe$Jtest["df"] > 0, 1 - pchisq(GMMe$Jtest["statistic"], GMMe$Jtest["df"]), NA)
  
  out       <- list(model.info  = list(n = n, ngroup = M, nvec = nvec, structural = structural, formula = formula, 
                                       excluded.instruments = excluded.instruments, gmm.weight = gmm.weight, 
                                       fixed.effects = fixed.effects, idX1 = idX1 + 1, idX2 = idX2 + 1, HAC = HAC),
                    gmm         = GMMe,
                    data        = list(y = y0, Gy = c(Gy0), X = X0, instruments = ins0, isolated = Is + 1, 
                                       non.isolated = nIs + 1, degree = dg))
  class(out) <- "linpeer"
  out
}


#' @rdname summary.qpeer
#' @export
summary.linpeer <- function(object, diagnostic = FALSE, diagnostics = FALSE, ...) {
  stopifnot(inherits(object, "linpeer"))
  diagn          <- NULL
  if (diagnostic || diagnostics) {
    diagn        <- fdiagnostic(object, nendo = "Gy")
  }
  coef           <- fcoef(Estimate = object$gmm$Estimate, cov = object$gmm$cov)
  out            <- c(object["model.info"], 
                      list(coefficients = coef, diagnostics = diagn),
                      object["gmm"], list(...))
  class(out)     <- "summary.linpeer"
  out
}

#' @rdname summary.qpeer
#' @export
print.summary.linpeer <- function(x, ...) {
  gmmw <- x$model.info$gmm.weight
  gmmw <- ifelse(gmmw == "ident", "Identity Matrix", ifelse(gmmw == "optimal", "Optimal", "IV"))
  hete <- x$model.info$HAC
  hete <- ifelse(hete == "iid", "IID", ifelse(hete == "hetero", "Individual", "Cluster"))
  sig  <- x$gmm$sigma
  FE   <- x$model.info$fixed.effects
  cat("Formula: ", as.character(x$model.info$formula),
      "\nExcluded instruments: ", as.character(x$model.info$excluded.instruments), 
      "\n\nModel: ", ifelse(x$model.info$structural, "Structural", "Reduced Form"),
      "\nFixed effects: ", paste0(toupper(substr(FE, 1, 1)), tolower(substr(FE, 2, nchar(FE)))), "\n", sep = "")
  
  coef       <- x$coefficients
  coef[,1:2] <- round(coef[,1:2], 7)
  coef[,3]   <- round(coef[,3], 5)
  cat("\nCoefficients:\n")
  fprintcoeft(coef)
  
  if (!is.null(x$diagnostics)) {
    coef       <- x$diagnostics
    coef[,3]   <- round(coef[,3], 5)
    cat("\nDiagnostic tests:\n")
    fprintcoeft(coef) 
  }
  cat("---\nSignif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01 \u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1\n\n")
  cat("GMM weight: ", gmmw, ", HAC: ", hete, sep = "")
  if (!is.null(sig)) {
    cat(", sigma: ", format(sig, digits = 5), sep = "")
  }
  cat("\nR-Squared: ", format(x$gmm$rsquared, digits = 5), 
      ", Adjusted R-squared: ", format(x$gmm$adjusted.rsquared, digits = 5), 
      "\nDegree of freedoms of residuals: ", x$gmm$df.residual, "\n", sep = "")
  class(x) <- "print.summary.linpeer"
  invisible(x)
}

#' @rdname summary.qpeer
#' @export
print.linpeer <- function(x, ...) {
  print(summary(x))
}