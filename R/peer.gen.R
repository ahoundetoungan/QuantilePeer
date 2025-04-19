#' @rdname qpeer
#' @export
genpeer <- function(formula, excluded.instruments, endogenous.variables, Glist, data, estimator = "IV", 
                    structural = FALSE, drop = NULL, fixed.effects = FALSE, 
                    HAC = "iid", checkrank = FALSE, 
                    compute.cov = TRUE, tol = 1e-10){
  # Estimator
  estimator <- tolower(estimator)
  stopifnot(estimator %in% c("iv", "gmm.optimal", "gmm.identity", "jive", "jive2"))
  estimator <- c("IV", "GMM.optimal", "GMM.identity", "JIVE", "JIVE2")[estimator == c("iv", "gmm.optimal", "gmm.identity", "jive", "jive2")]
  
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
  if ((HACnum == 1) & (FEnum != 0)) {
    HACnum   <- 2
    HAC      <- "cluster"
  }
  
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
  lIs      <- dg$lIs
  Is       <- dg$Is
  lnIs     <- dg$lnIs
  nIs      <- dg$nIs
  ldg      <- dg$ldg
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
  
  # Endo
  EnVar      <- as.formula(endogenous.variables); endogenous.variables <- EnVar
  if (length(EnVar) != 2) stop("The `endogenous.variables` argument must be in the format `~ y1 + y2 + ...`.")
  f.t.data   <- formula.to.data(formula = EnVar, data = data, fixed.effects = TRUE, simulations = TRUE)
  endo       <- f.t.data$X
  enname     <- f.t.data$xname
  Kendo      <- ncol(endo)
  
  # Instruments
  inst       <- as.formula(excluded.instruments); excluded.instruments <- inst
  if(length(inst) != 2) stop("The `excluded.instruments` argument must be in the format `~ z1 + z2 + ....`.")
  f.t.data   <- formula.to.data(formula = inst, data = data, fixed.effects = (fixed.effects != "no"), 
                                simulations = TRUE)
  ins        <- f.t.data$X
  zename     <- f.t.data$xname
  if (xint) {
    ins      <- ins[, zename != "(Intercept)"]
  } else {
    ins      <- ins
  }
  
  # Drop false isolated
  if (!is.null(drop)) {
    dg       <- fdrop(drop = drop, ldg = ldg, nvec = nvec, M = M, lIs = lIs, 
                      lnIs = lnIs, y = y, X = X, qy = endo, ins = ins)
    M        <- dg$M
    MIs      <- dg$MIs
    MnIs     <- dg$MnIs
    nvec     <- dg$nvec
    n        <- dg$n
    igr      <- dg$igr
    lIs      <- dg$lIs
    Is       <- dg$Is
    lnIs     <- dg$lnIs
    nIs      <- dg$nIs
    ldg      <- dg$ldg
    y        <- dg$y
    X        <- dg$X
    endo     <- dg$qy
    ins      <- dg$ins
    dg       <- dg$dg
  }
  
  # Demean fixed effect models
  # save original data
  y0         <- y
  endo0      <- endo
  X0         <- X
  ins0       <- ins
  if (fixed.effects != "no") {
    if (fixed.effects == "join") {
      y      <- c(demean(as.matrix(y), igroup = igr, ngroup = M))
      endo   <- demean(endo, igroup = igr, ngroup = M)
      X      <- demean(X, igroup = igr, ngroup = M)
      ins    <- demean(ins, igroup = igr, ngroup = M)
    } else {
      y      <- c(demean_separate(as.matrix(y), igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n))
      endo   <- demean_separate(endo, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n)
      X      <- demean_separate(X, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n)
      ins    <- demean_separate(ins, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n)
    }
    colnames(X)   <- xname
    colnames(ins) <- zename
  }
  
  # Remove useless columns
  idX1       <- 0:(ncol(X) - 1)
  tlm        <- idX1
  if (structural) {
    idX1     <- fcheckrank(X = X[Is + 1,], tol = tol)
    tlm      <- fcheckrank(X = X[nIs + 1,], tol = tol)
  } else {
    tlm      <- fcheckrank(X = X, tol = tol)
  }
  idX2       <- which(!(tlm %in% idX1)) - 1 
  idX1       <- which(tlm %in% idX1) - 1 
  X          <- X[, tlm + 1, drop = FALSE]
  X0         <- X0[, tlm + 1, drop = FALSE]
  xname      <- xname[tlm + 1]
  Kx         <- ncol(X)
  if (structural) {
    if (length(fcheckrank(X = cbind(endo, X)[nIs + 1,], tol = tol)) != (Kendo + Kx)) stop("The design matrix is not full rank.")
  } else {
    if (length(fcheckrank(X = cbind(endo, X), tol = tol)) != (Kendo + Kx)) stop("The design matrix is not full rank.")
  }
  
  
  if (structural) {
    ins      <- cbind(X[, idX2 + 1], ins)
    ins0     <- cbind(X0[, idX2 + 1], ins0)
    zename   <- c(xname[idX2 + 1], zename)
    if (checkrank) {
      tlm    <- fcheckrank(X = ins[nIs + 1,], tol = tol)
    }
  } else {
    ins      <- cbind(X, ins)
    ins0     <- cbind(X0, ins0)
    zename   <- c(xname, zename)
    if (checkrank) {
      tlm    <- fcheckrank(X = ins, tol = tol)
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
  iv         <- (estimator %in% c("IV", "GMM.optimal"))
  estname    <- NULL
  Kest       <- NULL
  if (structural) {
    Kx1      <- length(idX1)
    Kx2      <- length(idX2)
    if (Kins < Kx2 + Kendo) stop("Insufficient number of instruments: the model is not identified.")
    Kest1    <- ifelse(FEnum == 0, Kx1, Kx1 + MIs)
    Kest2    <- ifelse(FEnum == 0, Kx2 + Kendo + 1, Kx2 + Kendo + MnIs)
    if (length(Is) <= Kest1) stop("Insufficient number of isolated nodes for estimating the structural model.")
    if (length(nIs) <= Kest2) stop("Insufficient number of nonisolated nodes for estimating the structural model.")
    Kest     <- Kest1 + Kest2
    if (HACnum == 2 && (Kx1 >= MIs || Kins + 1 >= MnIs) && estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
      stop("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. HAC is set to 'iid' or 'hetero'.")
    }
    estname  <- c("Peers(conformity)", enname, xname)
    
    # Estimation
    GMMe     <- fstruct(y = y, X = X, qy = endo, ins = ins, idX1 = idX1, idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igr = igr, 
                        nIs = nIs, Is = Is, lnIs = lnIs, lIs = lIs, M = M, MnIs = MnIs, Kins = Kins, Kx = Kx, ntau = Kendo, 
                        Kest1 = Kest1, Kest2 = Kest2, n = n, HACnum = HACnum, iv = iv, estimator = estimator, 
                        compute.cov = compute.cov, estname = estname)
  } else {
    if (Kins < Kx + Kendo) stop("Insufficient number of instruments: the model is not identified.")
    Kest     <- ifelse(FEnum == 0, Kx + Kendo, ifelse(FEnum == 1, Kx + Kendo + M, Kx + Kendo + MIs + MnIs))
    if (n <= Kest) stop("Insufficient number of observations.")
    if (HACnum == 2 && Kins >= M && estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
      stop("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. HAC is set to 'iid' or 'hetero'.")
    }
    estname  <- c(enname, xname)
    V        <- cbind(endo, X)
    
    # Estimation
    GMMe     <- freduce(y = y, V = V, ins = ins, igr = igr, nvec = nvec, M = M, Kins = Kins, Kx = Kx, ntau = Kendo, 
                        Kest = Kest, n = n, HACnum = HACnum, iv = iv, estimator = estimator, compute.cov = compute.cov, 
                        estname = estname)
  }
  
  out       <- list(model.info  = list(n = n, ngroup = M, nvec = nvec, structural = structural, formula = formula, 
                                       endogenous.variables = endogenous.variables, excluded.instruments = excluded.instruments, 
                                       estimator = estimator, fixed.effects = fixed.effects, idXiso = idX1 + 1, idXniso = idX2 + 1, HAC = HAC,
                                       yname = yname, xnames = xname, znames = zename, endonames = enname),
                    gmm         = GMMe,
                    data        = list(y = y0, endogenous.variables = endo0, X = X0, instruments = ins0, isolated = Is + 1, 
                                       non.isolated = nIs + 1, degree = dg))
  class(out) <- "genpeer"
  out
}

#' @rdname summary.qpeer
#' @export
summary.genpeer <- function(object, fullparameters = TRUE, diagnostic = FALSE, diagnostics = FALSE, ...) {
  stopifnot(inherits(object, "genpeer"))
  if (is.null(object$gmm$cov)) {
    stop("The covariance matrix is not estimated.")
  }
  diagn          <- NULL
  cvKP           <- NULL
  if (diagnostic || diagnostics) {
    diagn        <- fdiagnostic(object, nendo = "endogenous.variables")
    cvKP         <- diagn$cvKP
    diagn        <- diagn$diag
  }
  
  coef           <- fcoef(Estimate = object$gmm$Estimate, cov = object$gmm$cov)
  out            <- c(object["model.info"], 
                      list(coefficients = coef, diagnostics = diagn, KP.cv = cvKP),
                      object["gmm"], list(...))
  class(out)     <- "summary.genpeer"
  out
}

#' @rdname summary.qpeer
#' @export
print.summary.genpeer <- function(x, ...) {
  esti <- x$model.info$estimator
  esti <- ifelse(esti == "GMM.identity", "GMM (Weight: Identity Matrix)", 
                 ifelse(esti == "GMM.optimal", "GMM (Weight: Optimal)",
                        ifelse(esti == "IV", "IV", 
                               ifelse(esti == "JIVE", "JIVE", "JIVE2"))))
  hete <- x$model.info$HAC
  hete <- ifelse(hete == "iid", "IID", ifelse(hete == "hetero", "Individual", "Cluster"))
  sig  <- x$gmm$sigma
  sig1 <- x$gmm$sigma1
  sig2 <- x$gmm$sigma2
  FE   <- x$model.info$fixed.effects
  cat("Formula: ", deparse(x$model.info$formula),
      "\nEndogenous variables: ", deparse(x$model.info$endogenous.variables), 
      "\nExcluded instruments: ", deparse(x$model.info$excluded.instruments), 
      "\n\nModel: ", ifelse(x$model.info$structural, "Structural", "Reduced Form"),
      "\nEstimator: ", esti,
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
  cat("---\nSignif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01 \u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1\n")
  
  cat("\nHAC: ", hete, sep = "")
  if (x$model.info$structural) {
    if (!is.null(sig1)) {
      if (!is.null(sig2)) {
        cat(", sigma (isolated): ", format(sig1, digits = 5), ", (non-isolated): ", format(sig2, digits = 5), sep = "")
      } else {
        cat(", sigma (isolated): ", format(sig1, digits = 5), sep = "")
      }
    }
  } else {
    if (!is.null(sig)) {
      cat(", sigma: ", format(sig, digits = 5), sep = "")
    }
  }
  cat("\nR-Squared: ", format(x$gmm$rsquared, digits = 5), 
      ", Adjusted R-squared: ", format(x$gmm$adjusted.rsquared, digits = 5), 
      "\nDegree of freedoms of residuals: ", x$gmm$df.residual, "\n", sep = "")
  class(x) <- "print.summary.genpeer"
  invisible(x)
}

#' @rdname summary.qpeer
#' @export
print.genpeer <- function(x, ...) {
  print(summary(x))
}