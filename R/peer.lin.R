#' @rdname qpeer
#' @export
linpeer <- function(formula, excluded.instruments, A, data, estimator = "IV", 
                    structural = FALSE, drop = NULL, fixed.effects = FALSE, HAC = "iid", 
                    checkrank = FALSE, compute.cov = TRUE, boot = 5e2, nthreads = 1, tol = 1e-10,
                    print = TRUE){
  nthreads   <- fnthreads(nthreads = nthreads)
  
  # Estimator
  estimator <- tolower(estimator)
  stopifnot(estimator %in% c("iv", "gmm.optimal", "gmm.identity", "jive", "jive2"))
  estimator <- c("IV", "GMM.optimal", "GMM.identity", "JIVE", "JIVE2")[estimator == c("iv", "gmm.optimal", "gmm.identity", "jive", "jive2")]
  
  # Variance structure
  HAC        <- tolower(HAC[1])
  stopifnot(HAC %in% c("iid", "hetero", "cluster", "cluster-bootstrap", "cluster-boot", "cboot", "c-boot"))
  HACnum     <- c(0:2, rep(3, 4))[HAC == c("iid", "hetero", "cluster", "cluster-bootstrap", "cluster-boot", "cboot", "c-boot")]
  HAC        <- c("iid", "hetero", "cluster", "cluster-bootstrap")[HACnum + 1]
  
  if ((HACnum == 3) & !(estimator %in% c("IV", "GMM.identity"))) {
    stop("Clustered bootstrap is implemented only for IV and for the GMM estimator with the identity weighting matrix.")
  }
  
  if ((HACnum == 3) & structural) {
    stop("Clustered bootstrap is not implemented for the structural model.")
  }
  
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
  if (!is.list(A)) {
    A  <- list(A)
  }
  d        <- fnetwork(A = A)
  G        <- d$G
  GIs      <- d$GIs
  GnIs     <- d$GnIs
  nvec     <- d$nvec
  ncum     <- c(0, cumsum(nvec))
  n        <- d$n
  igr      <- d$igr
  lIs      <- d$lIs
  Is       <- d$Is
  lnIs     <- d$lnIs
  nIs      <- d$nIs
  ld       <- d$ld
  d        <- d$d
  
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
  Ay         <- as.matrix(unlist(lapply(1:G, function(m) A[[m]] %*% y[(ncum[m] + 1):ncum[m + 1]])))
  
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
  
  # drop
  if (!is.null(drop)) {
    d        <- fdrop(drop = drop, ld = ld, nvec = nvec, G = G, lIs = lIs, 
                      lnIs = lnIs, y = y, X = X, qy = Ay, ins = ins)
    G        <- d$G
    GIs      <- d$GIs
    GnIs     <- d$GnIs
    nvec     <- d$nvec
    n        <- d$n
    igr      <- d$igr
    lIs      <- d$lIs
    Is       <- d$Is
    lnIs     <- d$lnIs
    nIs      <- d$nIs
    ld       <- d$ld
    y        <- d$y
    X        <- d$X
    Ay       <- d$qy
    ins      <- d$ins
    d        <- d$d
  }
  
  # Demean fixed effect models
  # save original data
  y0         <- y
  Ay0        <- Ay
  X0         <- X
  ins0       <- ins
  if (fixed.effects != "no") {
    if (fixed.effects == "join") {
      y      <- c(Demean(as.matrix(y), igroup = igr, ngroup = G))
      Ay     <- Demean(Ay, igroup = igr, ngroup = G)
      X      <- Demean(X, igroup = igr, ngroup = G)
      ins    <- Demean(ins, igroup = igr, ngroup = G)
    } else {
      y      <- c(Demean_separate(as.matrix(y), igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n))
      Ay     <- Demean_separate(Ay, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n)
      X      <- Demean_separate(X, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n)
      ins    <- Demean_separate(ins, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n)
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
    if (length(fcheckrank(X = cbind(Ay, X)[nIs + 1,], tol = tol)) != (1 + Kx)) stop("The design matrix is not full rank.")
  } else {
    if (length(fcheckrank(X = cbind(Ay, X), tol = tol)) != (1 + Kx)) stop("The design matrix is not full rank.")
  }
  
  
  ins        <- cbind(X, ins)
  ins0       <- cbind(X0, ins0)
  zename     <- c(xname, zename)
  if (checkrank) {
    tlm      <- fcheckrank(X = ins, tol = tol)
  }
  if (checkrank) {
    ins      <- ins[, tlm + 1, drop = FALSE]
    ins0     <- ins0[, tlm + 1, drop = FALSE]
    zename   <- zename[tlm + 1]
  }
  Kins       <- ncol(ins)
  
  # GMM
  seed       <- round(runif(1, 0, 1e9))
  GMMe       <- list()
  iv         <- (estimator %in% c("IV", "GMM.optimal"))
  estname    <- NULL
  Kest       <- NULL
  if (structural) {
    Kx1      <- length(idX1)
    Kx2      <- length(idX2)
    if (Kins < Kx2 + 1) stop("Insufficient number of instruments: the model is not identified.")
    Kest1    <- ifelse(FEnum == 0, Kx1, Kx1 + GIs)
    Kest2    <- ifelse(FEnum == 0, Kx2 + 2, Kx2 + 1 + GnIs)
    if (length(Is) <= Kest1) stop("Insufficient number of isolated nodes for estimating the structural model.")
    if (length(nIs) <= Kest2) stop("Insufficient number of nonisolated nodes for estimating the structural model.")
    Kest     <- Kest1 + Kest2
    if (HACnum == 2 && (Kx1 >= GIs || Kins + 1 >= GnIs) && estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
      stop("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. Set HAC to 'iid' or 'hetero'.")
    }
    estname  <- c(paste0(c("A(conformity):", "A(total):"), yname), xname)
    
    # Estimation
    GMMe     <- fstruct(y = y, X = X, qy = Ay, ins = ins, idX1 = idX1, idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igr = igr, 
                        nIs = nIs, Is = Is, lnIs = lnIs, lIs = lIs, G = G, GnIs = GnIs, Kins = Kins, Kx = Kx, ntau = 1, 
                        Kest1 = Kest1, Kest2 = Kest2, n = n, HACnum = HACnum, iv = iv, estimator = estimator, 
                        compute.cov = compute.cov, estname = estname)
  } else {
    if (Kins < Kx + 1) stop("Insufficient number of instruments: the model is not identified.")
    Kest     <- ifelse(FEnum == 0, Kx + 1, ifelse(FEnum == 1, Kx + 1 + G, Kx + 1 + GIs + GnIs))
    if (n <= Kest) stop("Insufficient number of observations.")
    if (HACnum == 2 && Kins >= G && estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
      stop("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. Set HAC to 'iid' or 'hetero', or use Bootstrap.")
    }
    estname  <- c(paste0("A:", yname), xname)
    V        <- cbind(Ay, X)
    
    # Estimation
    GMMe     <- freduce(y = y, V = V, ins = ins, igr = igr, nvec = nvec, G = G, Kins = Kins, Kx = Kx, ntau = 1, 
                        Kest = Kest, n = n, HACnum = HACnum, iv = iv, estimator = estimator, compute.cov = compute.cov, 
                        estname = estname, LnIs = lnIs, LIs = lIs, boot = boot, nthreads = nthreads, seed = seed,
                        print = print)
  }
  
  out       <- list(model.info  = list(n = n, ngroup = G, nvec = nvec, structural = structural, formula = formula, 
                                       excluded.instruments = excluded.instruments, estimator = estimator, 
                                       fixed.effects = fixed.effects, idXiso = idX1 + 1, idXniso = idX2 + 1, HAC = HAC, 
                                       yname = yname, xnames = xname, znames = zename, seed = seed, boot = boot),
                    gmm         = GMMe,
                    data        = list(y = y0, Ay = c(Ay0), X = X0, instruments = ins0, isolated = lapply(lIs, \(s) s + 1), 
                                       non.isolated = lapply(lnIs, \(s) s + 1), degree = ld))
  class(out) <- "linpeer"
  out
}


#' @rdname summary.qpeer
#' @export
summary.linpeer <- function(object, fullparameters = TRUE, diagnostic = FALSE, diagnostics = FALSE, 
                            boot = NULL, nthreads = 1L, print = TRUE, ...) {
  stopifnot(inherits(object, "linpeer"))
  if (is.null(object$gmm$cov)) {
    stop("The covariance matrix is not estimated.")
  }
  diagn          <- NULL
  cvKP           <- NULL
  seed           <- round(runif(1, 0, 1e9))
  boot           <- ifelse(is.null(boot), object$model.info$boot, boot)
  nthreads       <- fnthreads(nthreads = nthreads)
  if (diagnostic || diagnostics) {
    diagn        <- fdiagnostic(object, nendo = "Ay", seed = seed, boot = boot, 
                                nthreads = nthreads, print = print)
    cvKP         <- diagn$cvKP
    diagn        <- diagn$diag
  }
  
  if (fullparameters & object$model.info$structural) {
    yname        <- object$model.info$yname
    xnames       <- object$model.info$xnames
    est          <- object$gmm$Estimate
    covt         <- object$gmm$cov
    Kx1          <- length(object$model.info$idXiso)
    Kx2          <- length(object$model.info$idXniso)
    tp                  <- fStructParamFull(param = est, covp = covt, ntau = 1, Kx1 = Kx1, Kx2 = Kx2, quantile = 0) 
    tp$theta            <- c(tp$theta)
    names(tp$theta)     <- colnames(tp$Vpa) <- rownames(tp$Vpa) <- c(paste0(c("A(spillover):", "A(conformity):", "A(total):"), yname), xnames)
    object$gmm$Estimate <- tp$theta
    object$gmm$cov      <- tp$Vpa
  }
  coef           <- fcoef(Estimate = object$gmm$Estimate, cov = object$gmm$cov)
  out            <- c(object["model.info"], 
                      list(coefficients = coef, diagnostics = diagn, KP.cv = cvKP),
                      object["gmm"], list(...))
  class(out)     <- "summary.linpeer"
  out
}

#' @rdname summary.qpeer
#' @export
print.summary.linpeer <- function(x, ...) {
  esti <- x$model.info$estimator
  esti <- ifelse(esti == "GMM.identity", "GMM (Weight: Identity Matrix)", 
                 ifelse(esti == "GMM.optimal", "GMM (Weight: Optimal)",
                        ifelse(esti == "IV", "IV", 
                               ifelse(esti == "JIVE", "JIVE", "JIVE2"))))
  hete <- x$model.info$HAC
  hete <- ifelse(hete == "iid", "IID", ifelse(hete == "hetero", "Individual", ifelse(hete == "cluster", "Cluster", "Cluster (Bootstrap)")))
  sig  <- x$gmm$sigma
  sig1 <- x$gmm$sigma1
  sig2 <- x$gmm$sigma2
  FE   <- x$model.info$fixed.effects
  cat("Formula: ", deparse(x$model.info$formula),
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
  class(x) <- "print.summary.linpeer"
  invisible(x)
}

#' @rdname summary.qpeer
#' @export
print.linpeer <- function(x, ...) {
  print(summary(x, ...))
}


#' @title Simulating Linear Peer Effect Models
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param A The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `A` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\lambda', \gamma')'}, where  
#' \eqn{\lambda} is either the peer effect parameter for the reduced-form specification or a 2-vector with the first component being conformity peer effects and the second component representing total peer effects. 
#' The parameters \eqn{\lambda} and \eqn{\gamma} can also be specified separately using the arguments `lambda`, and `gamma`.
#' @param lambda The true value of the vector \eqn{\lambda}.
#' @param gamma The true value of the vector \eqn{\gamma}.
#' @param epsilon A vector of idiosyncratic error terms. If not specified, it will be simulated from a standard normal distribution. 
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `linpeer.sim` is called.
#' @param nthreads A strictly positive integer indicating the number of threads to use.
#' @param structural A logical value indicating whether simulations should be performed using the structural model. The default is the reduced-form model (see the Details section of \code{\link{qpeer}}).
#' @description
#' `linpeer.sim` simulates linear peer effect models.
#' @seealso \code{\link{qpeer}}, \code{\link{linpeer}}
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{Ay}{the average of y among friends.}
#' @examples 
#' set.seed(123)
#' ngr  <- 50
#' nvec <- rep(30, ngr)
#' n    <- sum(nvec)
#' A    <- lapply(1:ngr, function(z){
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z])
#'   diag(Gz) <- 0
#'   Gz/rowSums(Gz) # Row-normalized network
#' })
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' lam  <- 0.5
#' gam  <- c(2, -0.5, 1)
#' 
#' out  <- linpeer.sim(formula = ~ X, A = A, lambda = lam, gamma = gam)
#' summary(out$y)
#' @export
linpeer.sim   <- function(formula, A, parms, lambda, gamma, epsilon, 
                          structural = FALSE, nthreads = 1, data){
  if (!is.list(A)) {
    A  <- list(A)
  }
  d        <- fnetwork(A = A)
  G        <- d$G
  nvec     <- d$nvec
  n        <- d$n
  igr      <- d$igr
  Is       <- d$Is
  nIs      <- d$nIs
  d        <- d$d
  if (length(Is) <= 1 & structural) warning("The structural model requires isolated nodes.")
  
  # Data
  f.t.data <- formula.to.data(formula = formula, data = data, simulations = TRUE, fixed.effects = FALSE)
  formula  <- f.t.data$formula
  X        <- f.t.data$X
  if (nrow(X) != n) stop("The number of observations does not match the number of nodes in the network.")
  Kx       <- ncol(X)
  eps      <- NULL
  if(missing(epsilon)){
    eps    <- rnorm(n)
  } else{
    eps    <- c(epsilon)
    if (!(length(eps) %in% c(1, n))) stop("`epsilon` must be either a scalar or an n-dimensional vector.")
    if (length(eps) == 1) eps <- rep(eps, n)
  }
  
  # parameters
  lamst    <- NULL
  lam      <- NULL
  gam      <- NULL
  if (missing(parms)) {
    if (missing(lambda) | missing(gamma)) {
      stop("Define either `parms` or `lambda` and `gamma`.")
    }
    if (structural) {
      if (length(lambda) != 2){
        stop("length(lambda) is different from 2. See details on the structural model.")
      }
      lamst <- lambda[1]
      lam   <- lambda[2]
    } else {
      if (length(lambda) != 1){
        stop("lambda must be a scalar for the reduced-form model.")
      }
      lam   <- lambda
    }
    if (length(gamma) != Kx) stop("length(gamma) is different from ncol(X).")
    gam     <- gamma
  } else{
    if (!missing(lambda) | !missing(gamma)) {
      stop("Define either `parms` or `lambda` and `gamma`.")
    }
    if (structural) {
      if (length(parms) != (2 + Kx)) stop("length(parms) is different from 2 + ncol(X).")
      lamst <- parms[1]
      lam   <- parms[2]
    } else {
      if (length(parms) != (1 + Kx)) stop("length(parms) is different from 1 + ncol(X).")
      lam   <- parms[1]
    }
    gam     <- tail(parms, Kx)
  }
  if (sum(abs(lam)) >= 1) {
    warning("The absolute value of the total peer effects is greater than or equal to one, which may lead to multiple or no equilibria.")
  }
  if (structural && abs(lamst) >= 1) {
    stop("The absolute value of conformity peer effects must be strictly less than 1.")
  }
  
  # Solving the game
  ## talpha
  talpha <- c(X %*% gam + eps)
  if (structural) talpha[nIs + 1] <- talpha[nIs + 1]*(1 - lamst)
  
  y      <- fylim(A = A, talpha = talpha, igroup = igr, nvec = nvec, 
                  ngroup = G, lambda = lam, n = n, nthreads = nthreads)
  Ay     <- y$Ay
  y      <- y$y
  
  out    <- list("y"  = y,
                 "Ay" = Ay)
  out
}