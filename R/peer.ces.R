#' @title Estimation of CES-Based Peer Effects Models
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as \code{y ~ x1 + x2}, 
#' where `y` is the outcome and `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param instrument An object of class \link[stats]{formula} indicating the excluded instrument. It should be specified as \code{~ z},  
#' where `z` is the excluded instrument for the outcome. Following Boucher et al. (2024) \doi{10.3982/ECTA21048}, it can be an OLS exogenous prediction of `y`.  
#' This prediction is used to compute instruments for the CES function of peer outcomes.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `cespeer` is called.
#' @param tol A tolerance value used in the QR factorization to identify columns of explanatory variable and instrument matrices that ensure a full-rank matrix (see the \link[base]{qr} function).
#' The same tolerance is also used in the to minimize the concentrated GMM objective function (see \link[stats]{optimise}).
#' @param drop A dummy vector of the same length as the sample, indicating whether an observation should be dropped. 
#' This can be used, for example, to remove false isolates or to estimate the model only on non-isolated agents.
#' These observations cannot be directly removed from the network by the user because they may still be friends with other agents.
#' @param structural A logical value indicating whether the reduced-form or structural specification should be estimated (see details).
#' @param fixed.effects A logical value or string specifying whether the model includes subnet fixed effects. The fixed effects may differ between isolated and non-isolated nodes. Accepted values are `"no"` or `"FALSE"` (indicating no fixed effects), 
#' `"join"` or `TRUE` (indicating the same fixed effects for isolated and non-isolated nodes within each subnet), and `"separate"` (indicating different fixed effects for isolated and non-isolated nodes within each subnet). Note that `"join"` fixed effects are not applicable for structural models; 
#' `"join"` and `TRUE` are automatically converted to `"separate"`.
#' @param compute.cov A logical value indicating whether the covariance matrix of the estimator should be computed.
#' @param HAC A character string specifying the correlation structure among the idiosyncratic error terms for covariance computation. Options are `"iid"` for independent errors, `"hetero"` for heteroskedastic non-autocorrelated errors, and `"cluster"` for heteroskedastic errors with potential within-subnet correlation.
#' @param set.rho A fixed value for the CES substitution parameter to estimate a constrained model. Given this value, the other parameters can be estimated.  
#' @param grid.rho A finite grid of values for the CES substitution parameter \eqn{\rho} (see Details).  
#' This grid is used to obtain the starting value and define the GMM weight.  
#' It is recommended to use a finely subdivided grid.  
#' @param radius The radius of the subset in which the estimate for \eqn{\rho} is determined.  
#' The subset is a segment centered at the optimal \eqn{\rho} found using `grid.rho`.  
#' For better numerical optimization performance, use a finely subdivided `grid.rho` and a small `radius`.  
#' @description
#' `cespeer` estimates the CES-based peer effects model introduced by Boucher et al. (2024) \doi{10.3982/ECTA21048}. See Details.
#' @details 
#' Let \eqn{\mathcal{N}} be a set of \eqn{n} agents indexed by the integer \eqn{i \in [1, n]}.  
#' Agents are connected through a network characterized by an adjacency matrix \eqn{\mathbf{G} = [g_{ij}]} of dimension \eqn{n \times n}, where \eqn{g_{ij} = 1} if agent \eqn{j} is a friend of agent \eqn{i}, and \eqn{g_{ij} = 0} otherwise.  
#' In weighted networks, \eqn{g_{ij}} can be a nonnegative variable (not necessarily binary) that measures the intensity of the outgoing link from \eqn{i} to \eqn{j}. The model can also accommodate such networks.  
#' Note that the network generally consists of multiple independent subnets (e.g., schools).  
#' The `Glist` argument is the list of subnets. In the case of a single subnet, `Glist` should be a list containing one matrix.\cr
#' 
#' The reduced-form specification of the CES-based peer effects model is given by:
#' \deqn{y_i = \lambda\left(\sum_{j = 1}^n g_{ij}y_j^{\rho}\right)^{1/\rho} + \mathbf{x}_i^{\prime}\beta + \varepsilon_i,}
#' where \eqn{\varepsilon_i} is an idiosyncratic error term, \eqn{\lambda} captures the effect of the social norm \eqn{\left(\sum_{j = 1}^n g_{ij}y_j^{\rho}\right)^{1/\rho}},  
#' and \eqn{\beta} captures the effect of \eqn{\mathbf{x}_i} on \eqn{y_i}. The parameter \eqn{\rho} determines the form of the social norm in the model.  
#' - When \eqn{\rho > 1}, individuals are more sensitive to peers with high outcomes.  
#' - When \eqn{\rho < 1}, individuals are more sensitive to peers with low outcomes.  
#' - When \eqn{\rho = 1}, peer effects are uniform across peer outcome values.\cr
#' 
#' The structural specification of the model differs for isolated and non-isolated individuals.  
#' For an **isolated** \eqn{i}, the specification is similar to a standard linear-in-means model without social interactions, given by:
#' \deqn{y_i = \mathbf{x}_i^{\prime}\beta + \varepsilon_i.}
#' If node \eqn{i} is **non-isolated**, the specification is:
#' \deqn{y_i = \lambda\left(\sum_{j = 1}^n g_{ij}y_j^{\rho}\right)^{1/\rho} + (1 - \lambda^*)\mathbf{x}_i^{\prime}\beta + \varepsilon_i,}
#' where \eqn{\lambda^*} determines whether preferences exhibit conformity or complementarity/substitution.  
#' Identification of \eqn{\beta} and \eqn{\lambda^*} requires the network to include a sufficient number of isolated individuals.
#' @seealso \code{\link{qpeer}}, \code{\link{linpeer}}
#' @references Boucher, V., Rendall, M., Ushchev, P., & Zenou, Y. (2024). Toward a general theory of peer effects. Econometrica, 92(2), 543-565, \doi{10.3982/ECTA21048}.
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, including the number of subnets, the number of observations, and other key details.}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#'     \item{first.search}{A list containing initial estimations on the grid of values for \eqn{\rho}.}
#' @examples 
#' \donttest{
#' set.seed(123)
#' ngr  <- 50  # Number of subnets
#' nvec <- rep(30, ngr)  # Size of subnets
#' n    <- sum(nvec)
#' 
#' ### Simulating Data
#' ## Network matrix
#' G <- lapply(1:ngr, function(z) {
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z], nvec[z])
#'   diag(Gz) <- 0
#'   # Adding isolated nodes (important for the structural model)
#'   niso <- sample(0:nvec[z], 1, prob = (nvec[z] + 1):1 / sum((nvec[z] + 1):1))
#'   if (niso > 0) {
#'     Gz[sample(1:nvec[z], niso), ] <- 0
#'   }
#'   # Row-normalization
#'   rs   <- rowSums(Gz); rs[rs == 0] <- 1
#'   Gz/rs
#' })
#' 
#' X   <- cbind(rnorm(n), rpois(n, 2))
#' l   <- 0.55
#' b   <- c(2, -0.5, 1)
#' rho <- -2
#' eps <- rnorm(n, 0, 0.4)
#' 
#' ## Generating `y`
#' y <- cespeer.sim(formula = ~ X, Glist = G, rho = rho, lambda = l,
#'                  beta = b, epsilon = eps)$y
#' 
#' ### Estimation
#' ## Computing instruments
#' z <- fitted.values(lm(y ~ X))
#' 
#' ## Reduced-form model
#' rest <- cespeer(formula = y ~ X, instrument = ~ z, Glist = G, fixed.effects = "yes",
#'                 radius = 5, grid.rho = seq(-10, 10, 1))
#' summary(rest)
#' 
#' ## Structural model
#' sest <- cespeer(formula = y ~ X, instrument = ~ z, Glist = G, fixed.effects = "yes",
#'                 radius = 5, structural = TRUE, grid.rho = seq(-10, 10, 1))
#' summary(sest)
#' 
#' ## Quantile model
#' z    <- qpeer.inst(formula = ~ X, Glist = G, tau = seq(0, 1, 0.1), max.distance = 2, 
#'                    checkrank = TRUE)$instruments
#' qest <- qpeer(formula = y ~ X, excluded.instruments  = ~ z, Glist = G, 
#'               fixed.effects = "yes", tau = seq(0, 1, 1/3), structural = TRUE)
#' summary(qest)
#' }
#' @export
#' @importFrom stats optimise
cespeer <- function(formula, instrument, Glist, structural = FALSE, fixed.effects = FALSE,
                    set.rho = NULL, grid.rho = seq(-400, 400, radius), radius = 5, tol = 1e-8, 
                    drop = NULL, compute.cov = TRUE, HAC = "iid", data) {
  if (!is.null(set.rho)) {
    stopifnot(set.rho != 0)
  }
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
  # MIs      <- dg$MIs
  # MnIs     <- dg$MnIs
  nvec     <- dg$nvec
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
  formula  <- as.formula(formula)
  f.t.data <- formula.to.data(formula = formula, data = data, fixed.effects = (fixed.effects != "no"),
                              simulations = FALSE)
  y        <- f.t.data$y
  X        <- f.t.data$X
  xname    <- f.t.data$xname
  yname    <- f.t.data$yname
  Kx       <- ncol(X)
  # xname      <- colnames(X)
  # yname      <- "y"
  
  
  # Instrument
  inst     <- as.formula(instrument); instrument <- inst
  if(length(inst) != 2) stop("The `excluded.instruments` argument must be in the format `~ z` (a single instrument).")
  f.t.data <- formula.to.data(formula = inst, data = data, fixed.effects = (fixed.effects != "no"), 
                              simulations = TRUE)
  z       <- f.t.data$X
  zename  <- f.t.data$xname
  z       <- z[, zename != "(Intercept)"]
  zename  <- zename[zename != "(Intercept)"]
  if (!is.null(dim(z))) {
    stop("Only one excluded instrument can be used.")
  }
  
  # Create additional variables)
  hasIso   <- fCESdatainit(y = y, z = z, G = Glist, nvec = nvec, M = M, ldg = ldg, lIs = lIs, lnIs = lnIs, drop = drop)
  frindex  <- hasIso$friendindex
  frzeroy  <- hasIso$frzeroy
  frzeroz  <- hasIso$frzeroz
  ldg_st   <- hasIso$ldg
  dg_st    <- hasIso$dg
  M_st     <- hasIso$M
  MIs_st   <- hasIso$MIs
  MnIs_st  <- hasIso$MnIs
  yFMiMa   <- cbind(hasIso$yFmin, hasIso$yFmax)
  zFMiMa   <- cbind(hasIso$zFmin, hasIso$zFmax)
  lIs      <- hasIso$lIs # After selection
  Is       <- hasIso$Is # After selection
  lnIs     <- hasIso$lnIs # After selection
  nIs      <- hasIso$nIs # After selection
  hasIso   <- hasIso$hasIso
  niso     <- length(Is)
  nniso    <- length(nIs)
  n_st     <- niso + nniso
  sel      <- sort(c(Is, nIs))
  if ((1 %in% frzeroy) | (1 %in% frzeroz)) {
    stop("The outcome `y` and the instrument `z` must be strictly positive.")
  }
  
  # Index for Isolated to form a full rank matrix
  tp      <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, igroup = igr, frzeroy = frzeroy, 
                      frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                      n = n, Kx = Kx, ngroup = M, rho = 1, FEnum = FEnum, deriv = FALSE)
  idXiso   <- fcheckrank(X = tp[Is + 1, 1:Kx, drop = FALSE], tol = tol)
  Kxiso    <- length(idXiso)
  if (length(idXiso) == 0 & structural) {
    stop("The regressor matrix X for isolated nodes does not have sufficiently independent nonlinear columns for estimating the structural model.")
  }
  idXniso  <- fcheckrank(X = tp[nIs + 1, 1:Kx, drop = FALSE], tol = tol)
  idXniso  <- setdiff(idXniso, idXiso)
  Kxniso   <- length(idXniso)
  if ((Kxiso + Kxniso) != Kx) {
    stop("The regressor matrix X is not full rank (perhaps it includes a constant variable when fixed effects are required?).")
  }
  
  # Number of parameters estimated (for degree of freedom)
  Kest       <- NA
  Kendo      <- NA #Parameters associated with y
  nKendo     <- NA # names
  Kest1      <- NA
  Kest2      <- NA
  if (structural) {
    Kest1    <- ifelse(FEnum == 0, Kx, Kx + MIs_st)
    Kest2    <- ifelse(FEnum == 0, Kx + 3, Kx + 3 + MnIs_st)
    Kest     <- Kest
    Kendo    <- 3
    nKendo   <- c("rho", paste0(c("G(conformity):", "G(total):"), yname))
    if (niso <= Kest1) stop("Insufficient number of isolated nodes for estimating the structural model.")
    if (nniso <= Kest2) stop("Insufficient number of nonisolated nodes for estimating the structural model.")
  } else {
    Kest     <- ifelse(FEnum == 0, Kx + 2, ifelse(FEnum == 1, Kx + 2 + M, Kx + 2 + MIs_st + MnIs_st))
    if (nniso <= Kest) stop("Insufficient number of isolated nodes.")
    Kendo    <- 2
    nKendo   <- c("rho", paste0("G:", yname))
  }
  
  # beta1 and W
  b1    <- numeric()
  W     <- NULL
  if (structural) {
    b1  <- fOLS(data = tp, idX1 = idXiso, Is = Is, Kx = Kx)
    W   <- diag(Kxniso + 3)
  } else {
    W   <- diag(Kx + 2)
  }
  
  # Optimization
  opt   <- list()
  fes  <- NULL
  if (is.null(set.rho)) {
    grid.rho <- grid.rho[grid.rho != 0]
    # First estimation when W = I
    Fsear    <- sapply(grid.rho, function(s){
      tp1    <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, igroup = igr, frzeroy = frzeroy, 
                         frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, nvec = nvec, n = n, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                         Kx = Kx, ngroup = M, rho = s, FEnum = FEnum, deriv = compute.cov)
      tp2    <- fCESgmmrhoparms(rho = s, beta1 = b1, data = tp1, sel = sel, nIs = nIs, Is = Is, idX1 = idXiso, 
                                idX2 = idXniso, igroup = igr, ngroup = M, Kx = Kx, Kx2 = Kxniso, nniso = nniso, niso = niso, n = n, 
                                nst = n_st, Kest1 = Kest1, Kest2 = Kest2, Kest = Kest, structural = structural, COV = compute.cov)
      fCESgmmrhoobj(theta = tp2$theta, data = tp1, sel = sel, nIs = nIs, idX1 = idXiso, 
                    idX2 = idXniso, Kx = Kx, Kx2 = Kxniso, nniso = nniso, nst = n_st, structural = structural)
    })
    rhoi <- grid.rho[which.min(Fsear)]
    fes  <- fCESgmmparms(rho = rhoi, beta1 = b1, idX1 = idXiso, idX2 = idXniso, X = X, y = y, z = z, G = Glist, 
                         friendindex = frindex, igroup = igr, frzeroy = frzeroy, frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, 
                         nIs = nIs, Is = Is, sel = sel, W = W, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, n = n, nniso = nniso, 
                         niso = niso, nst = n_st, Kx = Kx, Kx2 = Kxniso, ngroup = M, FEnum = FEnum, Kest1 = Kest1, Kest2 = Kest2, 
                         Kest = Kest, structural = structural, HAC = HACnum, COV = FALSE)$theta
    
    # Second estimation
    tp  <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, igroup = igr, frzeroy = frzeroy, 
                    frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, nvec = nvec, n = n, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                    Kx = Kx, ngroup = M, rho = fes[1], FEnum = FEnum, deriv = FALSE)
    W   <- fCESWeight(theta = fes, data = tp, sel = sel,
                      nIs = nIs, idX1 = idXiso, idX2 = idXniso, structural = structural, Kx = Kx, Kx2 = Kxniso,
                      nniso = nniso, nst = n_st)
    tp  <- optimise(fCESgmmobj, tol = tol, interval = c(fes[1] - radius, fes[1] + radius),
                    beta1 = b1, idX1 = idXiso, idX2 = idXniso, X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                    igroup = igr, frzeroy = frzeroy, frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, nIs = nIs, sel = sel,
                    W = W, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, n = n, nniso = nniso, nst = n_st, Kx = Kx, 
                    Kx2 = Kxniso, ngroup = M, FEnum = FEnum, structural = structural)
    opt$objective <- tp$objective
    Est <- fCESgmmparms(rho = tp$minimum, beta1 = b1, idX1 = idXiso, idX2 = idXniso, X = X, y = y, z = z, G = Glist, 
                        friendindex = frindex, igroup = igr, frzeroy = frzeroy, frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, 
                        nIs = nIs, Is = Is, sel = sel, W = W, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, n = n, nniso = nniso, 
                        niso = niso, nst = n_st, Kx = Kx, Kx2 = Kxniso, ngroup = M, FEnum = FEnum, Kest1 = Kest1, Kest2 = Kest2, 
                        Kest = Kest, structural = structural, HAC = HACnum, COV = compute.cov)
    
    # print(Est$theta)
    opt <- c(opt, fCESParam(param = Est$theta, covp = Est$Vpa, idX1 = idXiso, idX2 = idXniso, Kx = Kx, Kx2 = Kxniso, estimrho = 1, 
                            structural = structural, COV = compute.cov), Est[c("sigma21", "sigma22", "sigma2")])
    opt$theta <- c(opt$theta)
    fes <- c(fCESParam(param = fes, covp = matrix(0), idX1 = idXiso, idX2 = idXniso, Kx = Kx, Kx2 = Kxniso, estimrho = 1, 
                       structural = structural, COV = FALSE)$theta)
    names(fes) <- c(nKendo, xname)
    fes        <- list(grid.rho = grid.rho, objective = Fsear, theta = fes)
  } else {
    rhoinf <- ifelse(is.finite(set.rho), 0, 1)
    tp  <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, igroup = igr, frzeroy = frzeroy, 
                    frzeroz = frzeroz, lIs = lIs, lnIs = lnIs, nvec = nvec, n = n, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                    Kx = Kx, ngroup = M, rho = set.rho, FEnum = FEnum, deriv = compute.cov)
    Est  <- fCESgmmrhoparms(rho = set.rho, beta1 = b1, data = tp, sel = sel, nIs = nIs, Is = Is, idX1 = idXiso, 
                            idX2 = idXniso, igroup = igr, ngroup = M, Kx = Kx, Kx2 = Kxniso, nniso = nniso, niso = niso, n = n, 
                            nst = n_st, Kest1 = Kest1, Kest2 = Kest2, Kest = Kest, rhoinf = rhoinf, structural = structural, 
                            COV = compute.cov)
    opt$objective <- fCESgmmrhoobj(theta = Est$theta, data = tp, sel = sel, nIs = nIs, idX1 = idXiso, 
                                   idX2 = idXniso, Kx = Kx, Kx2 = Kxniso, nniso = nniso, nst = n_st, structural = structural, rhoinf)
    opt  <- c(opt, fCESParam(param = Est$theta, covp = Est$Vpa, idX1 = idXiso, idX2 = idXniso, Kx = Kx, Kx2 = Kxniso, estimrho = 0, 
                             structural = structural, COV = compute.cov), Est[c("sigma21", "sigma22", "sigma2")])
  }
  opt$theta           <- c(opt$theta)
  names(opt$theta)    <- c(nKendo, xname)
  if (compute.cov) {
    colnames(opt$Vpa) <- rownames(opt$Vpa) <- c(nKendo, xname)
  } else{
    opt$Vpa           <- NULL
  }
  sigmaiso    <- sqrt(opt$sigma21);  if(is.na(sigmaiso)) sigmaiso = NULL
  sigmaniso   <- sqrt(opt$sigma22);  if(is.na(sigmaniso)) sigmaniso = NULL
  sigma       <- sqrt(opt$sigma2);  if(is.na(sigma)) sigma = NULL
  SIGMA       <- NULL
  if (structural) {
    SIGMA     <- list(sigma1 = sigmaiso, sigma2 = sigmaniso)
  } else {
    SIGMA     <- list(sigma = sigma)
  }
  Est         <- c(list(Estimate = opt$theta, cov = opt$Vpa, objective = opt$objective),  SIGMA,
                   list(counts = opt$counts, convergence = opt$convergence, message = opt$message))
  out         <- list(model.info   = list(n = n_st, ngroup = M, nvec = nvec, structural = structural, formula = formula, 
                                          instrument = instrument, fixed.effects = fixed.effects, idXiso = idXiso + 1, idXniso = idXniso + 1, 
                                          HAC = HAC, set.rho = set.rho, yname = yname, xnames = xname, zname = zename),
                      gmm          = Est,
                      first.search = fes)
  class(out)  <- "cespeer"
  out
}


#' @title Summary for the Estimation of CES-based Peer Effects Models
#' @param object An object of class \code{\link{cespeer}}.
#' @param ... Further arguments passed to or from other methods.
#' @param x An object of class \code{\link{summary.cespeer}} or \code{\link{cespeer}}.
#' @param fullparameters A logical value indicating whether all parameters should be summarized (may be useful for the structural model).
#' @description Summary and print methods for the class \code{\link{cespeer}}.
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, such as the number of subnets, number of observations, and other key details.}
#'     \item{coefficients}{A summary of the estimates, standard errors, and p-values.}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#' @export
summary.cespeer <- function(object, fullparameters = TRUE, ...) {
  stopifnot(inherits(object, "cespeer"))
  if (is.null(object$gmm$cov)) {
    stop("The covariance matrix is not estimated.")
  }
  if (fullparameters & object$model.info$structural) {
    yname        <- object$model.info$yname
    xnames       <- object$model.info$xnames
    est          <- object$gmm$Estimate
    covt         <- object$gmm$cov
    Kx1          <- length(object$model.info$idXiso)
    Kx2          <- length(object$model.info$idXniso)
    tp                  <- fStructParamFull(param = est, covp = covt, ntau = 1, Kx1 = Kx1, Kx2 = Kx2, quantile = 0, ces = TRUE) 
    tp$theta            <- c(tp$theta)
    names(tp$theta)     <- colnames(tp$Vpa) <- rownames(tp$Vpa) <- c(c("rho", paste0(c("G(spillover):", "G(conformity):", "G(total):"), yname)), xnames)
    object$gmm$Estimate <- tp$theta
    object$gmm$cov      <- tp$Vpa
  }
  coef           <- fcoef(Estimate = object$gmm$Estimate, cov = object$gmm$cov)
  if (!is.null(object$model.info$set.rho)) {
    coef[1, -1]  <- NA
  }
  out            <- c(object["model.info"], 
                      list(coefficients = coef),
                      object["gmm"], list(...))
  class(out)     <- "summary.cespeer"
  out
}

#' @rdname summary.cespeer
#' @export
print.summary.cespeer <- function(x, ...) {
  hete <- x$model.info$HAC
  hete <- ifelse(hete == "iid", "IID", ifelse(hete == "hetero", "Individual", "Cluster"))
  sig  <- x$gmm$sigma
  sig1 <- x$gmm$sigma1
  sig2 <- x$gmm$sigma2
  FE   <- x$model.info$fixed.effects
  cat("Formula: ", deparse(x$model.info$formula),
      "\nExcluded instrument: ", deparse(x$model.info$instrument), 
      "\n\nModel: ", ifelse(x$model.info$structural, "Structural", "Reduced Form"),
      "\nFixed effects: ", paste0(toupper(substr(FE, 1, 1)), tolower(substr(FE, 2, nchar(FE)))), "\n", sep = "")
  
  coef       <- x$coefficients
  coef[,1:2] <- round(coef[,1:2], 7)
  coef[,3]   <- round(coef[,3], 5)
  
  cat("\nCoefficients:\n")
  fprintcoeft(coef)
  
  cat("---\nSignif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01 \u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1\n\n")
  cat("HAC: ", hete, sep = "")
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
  cat("\nCES parameter -- testing whether rho = 1: prob = ", 
      ifelse(is.null(x$model.info$set.rho), round(2*(1 - pnorm(abs((coef[1, 1] - 1)/coef[1, 2]))), 5), NA), "\n")
  class(x) <- "print.summary.cespeer"
  invisible(x)
}

#' @rdname summary.cespeer
#' @export
print.cespeer <- function(x, ...) {
  print(summary(x))
}


#' @title Simulating Peer Effect Models with a CES Social Norm
#' @param formula A formula object (\link[stats]{formula}): a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\rho, \lambda', \beta')'}, where \eqn{\rho} is the substitution parameter of the CES function and 
#' \eqn{\lambda} is either the peer effect parameter for the reduced-form specification or a 2-vector with the first component being conformity peer effects and the second component representing total peer effects. 
#' The parameters \eqn{\rho}, \eqn{\lambda}, and \eqn{\beta} can also be specified separately using the arguments `rho`, `lambda`, and `beta` (see the Details section of \code{\link{cespeer}}).
#' @param rho The true value of the substitution parameter of the CES function.
#' @param lambda The true value of the peer effect parameter \eqn{\lambda}. It must include conformity and total peer effects for the structural model.
#' @param beta The true value of the vector \eqn{\beta}.
#' @param epsilon A vector of idiosyncratic error terms. If not specified, it will be simulated from a standard normal distribution (see the model specification in the Details section of \code{\link{cespeer}}). 
#' @param maxit The maximum number of iterations for the Fixed Point Iteration Method.
#' @param data An optional data frame, list, or environment containing the model variables. If a variable is not found in `data`, it is retrieved from \code{environment(formula)}, typically the environment from which `sim.cespeer` is called.
#' @param tol The tolerance value used in the Fixed Point Iteration Method to compute the outcome `y`. The process stops if the \eqn{\ell_1}-distance 
#' between two consecutive values of `y` is less than `tol`.
#' @param init An optional initial guess for the equilibrium.
#' @param structural A logical value indicating whether simulations should be performed using the structural model. The default is the reduced-form model (see the Details section of \code{\link{cespeer}}).
#' @description
#' `cespeer.sim` simulates peer effect models with a Constant Elasticity of Substitution (CES) based social norm (Boucher et al., 2024).
#' @seealso \code{\link{cespeer}}, \code{\link{qpeer.sim}}
#' @references Boucher, V., Rendall, M., Ushchev, P., & Zenou, Y. (2024). Toward a general theory of peer effects. *Econometrica, 92*(2), 543-565, \doi{10.3982/ECTA21048}.
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{epsilon}{The idiosyncratic error.}
#'     \item{init}{The initial guess.}
#'     \item{iteration}{The number of iterations before convergence.}
#' @examples 
#' set.seed(123)
#' ngr  <- 50
#' nvec <- rep(30, ngr)
#' n    <- sum(nvec)
#' G    <- lapply(1:ngr, function(z){
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z])
#'   diag(Gz) <- 0
#'   Gz/rowSums(Gz) # Row-normalized network
#' })
#' tau  <- seq(0, 1, 0.25)
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' l    <- 0.55
#' rho  <- 3
#' b    <- c(4, -0.5, 1)
#' 
#' out  <- cespeer.sim(formula = ~ X, Glist = G, rho = rho, lambda = l, beta = b)
#' summary(out$y)
#' out$iteration
#' @export
cespeer.sim <- function(formula, Glist, parms, rho, lambda, beta, epsilon, structural = FALSE, 
                        init, tol = 1e-10, maxit = 500, data){
  # Network
  
  dg       <- fnetwork(Glist = Glist)
  M        <- dg$M
  nvec     <- dg$nvec
  n        <- dg$n
  igr      <- dg$igr
  Is       <- dg$Is
  nIs      <- dg$nIs
  ldg      <- dg$ldg
  dg       <- dg$dg
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
  b        <- NULL
  if (missing(parms)) {
    if (missing(lambda) | missing(beta) | missing(rho)) {
      stop("Define either `parms` or `rho`, `lambda`, and `beta`.")
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
    if (length(beta) != Kx) stop("length(beta) is different from ncol(X).")
    b      <- beta
  } else{
    if (!missing(rho) | !missing(lambda) | !missing(beta)) {
      stop("Define either `parms` or `rho`, `lambda`, and `beta`.")
    }
    if (structural) {
      if (length(parms) != (3 + Kx)) stop("length(parms) is different from 3 + ncol(X). See details on the structural model.")
      rho   <- parms[1]
      lamst <- parms[2]
      lam   <- parms[3]
    } else {
      if (length(parms) != (2 + Kx)) stop("length(parms) is different from 2 + ncol(X).")
      rho   <- parms[1]
      lam   <- parms[2]
    }
    b      <- tail(parms, Kx)
  }
  if (sum(abs(lam)) >= 1) {
    warning("The absolute value of the total peer effects is greater than or equal to one, which may lead to multiple or no equilibria.")
  }
  if (structural && abs(lamst) >= 1) {
    stop("The absolute value of conformity peer effects must be strictly less than 1.")
  }
  if (rho == 0) {
    stop("`rho` cannot be zero.")
  }
  
  # Solving the game
  ## talpha
  talpha   <- X %*% b + eps
  if (structural) talpha[nIs + 1] <- talpha[nIs + 1]*(1 - lamst)
  if (any(talpha <= 0)) {
    stop("`x*beta + epsilon` is not strictly positive.")
  }
  
  ## init
  if (missing(init)) {
    init   <- rep(max(talpha)/(1 - lam))
  }
  if (length(init) == 1){
    init   <- rep(init, n)
  } else if (length(init) != n) {
    stop("`init` is not an n-vector.")
  }
  y        <- unlist(init)

  ## other variables
  ncs      <- c(0, cumsum(nvec))
  friendindex <- lapply(1:M, function(m) {
    lapply(1:nvec[m], function(s) {
      which(Glist[[m]][s,] > 0) - 1
    })})
  frzeroy  <- as.integer(unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      any(y[friendindex[[m]][[s]] + ncs[m] + 1] <= 0)
    })})))
  yFmax    <- unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, max(y[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  yFmin    <- unlist(lapply(1:M, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, min(y[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  yFMiMa   <- cbind(yFmin, yFmax)
  if (any(frzeroy == 1)) {
    stop("`init` is not strictly positive.")
  }
  
  # Compute equilibrium
  t        <- fNashECES(y = y, G = Glist, talpha = talpha, lambda = lam, rho = rho, friendindex = friendindex, 
                        igroup = igr, frzeroy = frzeroy, nvec = nvec, yFMiMa = yFMiMa, ngroup = M, n = n, tol = tol, 
                        maxit = maxit)
  # Output
  list("y"         = c(y),
       "epsilon"   = eps,
       "init"      = init,
       "iteration" = t)
} 

