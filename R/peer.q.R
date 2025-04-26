#' @title Estimating Peer Effects Models
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as \code{y ~ x1 + x2}, 
#' where `y` is the outcome and `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param excluded.instruments An object of class \link[stats]{formula} to indicate excluded instruments. It should be specified as \code{~ z1 + z2}, 
#' where `z1` and `z2` are excluded instruments for the quantile peer outcomes.
#' @param endogenous.variables An object of class \link[stats]{formula} that allows specifying endogenous variables. It is used to indicate the peer variables whose effects will be estimated. These can include average peer variables, quantile peer variables, 
#' or a combination of multiple variables. It should be specified as \code{~ y1 + y2}, where `y1` and `y2` are the endogenous peer variables.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param tau A numeric vector specifying the quantile levels.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see the \link[stats]{quantile} function).
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `qpeer` is called.
#' @param tol A tolerance value used in the QR factorization to identify columns of explanatory variable and instrument matrices that ensure a full-rank matrix (see the \link[base]{qr} function).
#' @param drop A dummy vector of the same length as the sample, indicating whether an observation should be dropped. 
#' This can be used, for example, to remove false isolates or to estimate the model only on non-isolated agents.
#' These observations cannot be directly removed from the network by the user because they may still be friends with other agents.
#' @param structural A logical value indicating whether the reduced-form or structural specification should be estimated (see Details).
#' @param fixed.effects A logical value or string specifying whether the model includes subnet fixed effects. The fixed effects may differ between isolated and non-isolated nodes. Accepted values are `"no"` or `"FALSE"` (indicating no fixed effects), 
#' `"join"` or `TRUE` (indicating the same fixed effects for isolated and non-isolated nodes within each subnet), and `"separate"` (indicating different fixed effects for isolated and non-isolated nodes within each subnet). Note that `"join"` fixed effects are not applicable for structural models; 
#' `"join"` and `TRUE` are automatically converted to `"separate"`.
#' @param estimator A character string specifying the estimator to be used. The available options are:
#' `"IV"` for the standard instrumental variable estimator,
#' `"gmm.identity"` for the GMM estimator with the identity matrix as the weight,
#' `"gmm.optimal"` for the GMM estimator with the optimal weight matrix,
#' `"JIVE"` for the Jackknife instrumental variable estimator, and
#' `"JIVE2"` for the Type 2 Jackknife instrumental variable estimator.
#' @param compute.cov A logical value indicating whether the covariance matrix of the estimator should be computed.
#' @param HAC A character string specifying the correlation structure among the idiosyncratic error terms for covariance computation. Options are `"iid"` for independent errors, `"hetero"` for heteroskedastic non-autocorrelated errors, and `"cluster"` for heteroskedastic errors with potential within-subnet correlation.
#' @param checkrank A logical value indicating whether the instrument matrix should be checked for full rank. If the matrix is not of full rank, unimportant columns will be removed to obtain a full-rank matrix.
#' @description
#' `qpeer` estimates the quantile peer effect models introduced by Houndetoungan (2025). In the \code{\link{linpeer}} function, quantile peer variables are replaced with the average peer variable, and they can be replaced with other peer variables in the \code{\link{genpeer}} function.
#' @details 
#' Let \eqn{\mathcal{N}} be a set of \eqn{n} agents indexed by the integer \eqn{i \in [1, n]}. 
#' Agents are connected through a network that is characterized by an adjacency matrix \eqn{\mathbf{G} = [g_{ij}]} of dimension \eqn{n \times n}, where \eqn{g_{ij} = 1} if agent \eqn{j} is a friend of agent \eqn{i}, and \eqn{g_{ij} = 0} otherwise. 
#' In weighted networks, \eqn{g_{ij}} can be a nonnegative variable (not necessarily binary) that measures the intensity of the outgoing link from \eqn{i} to \eqn{j}. The model can also accommodate such networks. Note that the network is generally constituted in many independent subnets (eg: schools). 
#' The `Glist` argument is the list of subnets. In the case of a single subnet, `Glist` will be a list containing one matrix.\cr
#' 
#' Let \eqn{\mathcal{T}} be a set of quantile levels. The reduced-form specification of quantile peer effect models is given by:
#' \deqn{y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + \mathbf{x}_i^{\prime}\beta + \varepsilon_i,}
#' where \eqn{\mathbf{y}_{-i} = (y_1, \ldots, y_{i-1}, y_{i+1}, \ldots, y_n)^{\prime}} is the vector of outcomes for other units, and \eqn{q_{\tau,i}(\mathbf{y}_{-i})} is the 
#' sample \eqn{\tau}-quantile of peer outcomes. The term \eqn{\varepsilon_i} is an idiosyncratic error term, \eqn{\lambda_{\tau}} captures the effect of the \eqn{\tau}-quantile of peer outcomes on \eqn{y_i}, 
#' and \eqn{\beta} captures the effect of \eqn{\mathbf{x}_i} on \eqn{y_i}. For the definition of the sample \eqn{\tau}-quantile, see Hyndman and Fan (1996). 
#' If the network matrix is weighted, the sample weighted quantile can be used, where the outcome for friend \eqn{j} of \eqn{i} is weighted by \eqn{g_{ij}}. It can be shown that
#' the sample \eqn{\tau}-quantile is a weighted average of two peer outcomes. For more details, see the \link[stats]{quantile} and \code{\link{qpeer.instruments}} functions. \cr
#' 
#' The quantile \eqn{q_{\tau,i}(\mathbf{y}_{-i})} can be replaced with the average peer variable in \code{\link{linpeer}} or with other measures in \code{\link{genpeer}} through the `endogenous.variables` argument. 
#' In \code{\link{genpeer}}, it is possible to specify multiple peer variables, such as male peer averages and female peer averages. Additionally, both quantiles and averages can be included (\code{\link{genpeer}} 
#' is general and encompasses \code{\link{qpeer}} and \code{\link{linpeer}}). See examples. \cr
#' 
#' One issue in linear peer effect models is that individual preferences with conformity and complementarity/substitution lead to the same reduced form. 
#' However, it is possible to disentangle both types of preferences using isolated individuals (individuals without friends). 
#' The structural specification of the model differs between isolated and nonisolated individuals.
#' For isolated \eqn{i}, the specification is similar to a standard linear-in-means model without social interactions, given by:
#' \deqn{y_i = \mathbf{x}_i^{\prime}\beta + \varepsilon_i.}
#' If node \eqn{i} is non-isolated, the specification is given by:
#' \deqn{y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + (1 - \lambda_2)(\mathbf{x}_i^{\prime}\beta  + \varepsilon_i),}
#' where \eqn{\lambda_2} determines whether preferences exhibit conformity or complementarity/substitution. In general, \eqn{\lambda_2 > 0} and this means that that preferences are conformist (anti-conformity may be possible in some models when \eqn{\lambda_2 < 0}). 
#' In contrast, when \eqn{\lambda_2 = 0}, there is complementarity/substitution between individuals depending on the signs of the \eqn{\lambda_{\tau}} parameters.
#' It is obvious that \eqn{\beta} and \eqn{\lambda_2} can be identified only if the network includes enough isolated individuals.
#' @seealso \code{\link{qpeer.sim}}, \code{\link{qpeer.instruments}}
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, such as the number of subnets, number of observations, and other key details.}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#'     \item{data}{A list containing the outcome, outcome quantiles among peers, control variables, and excluded instruments used in the model.}
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
#'   Gz
#' })
#' 
#' tau <- seq(0, 1, 1/3)
#' X   <- cbind(rnorm(n), rpois(n, 2))
#' l   <- c(0.2, 0.15, 0.1, 0.2)
#' b   <- c(2, -0.5, 1)
#' eps <- rnorm(n, 0, 0.4)
#' 
#' ## Generating `y`
#' y <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = l, 
#'                beta = b, epsilon = eps)$y
#' 
#' ### Estimation
#' ## Computing instruments
#' Z <- qpeer.inst(formula = ~ X, Glist = G, tau = seq(0, 1, 0.1), 
#'                 max.distance = 2, checkrank = TRUE)
#' Z <- Z$instruments
#' 
#' ## Reduced-form model 
#' rest <- qpeer(formula = y ~ X, excluded.instruments = ~ Z, Glist = G, tau = tau)
#' summary(rest)
#' summary(rest, diagnostic = TRUE)  # Summary with diagnostics
#' 
#' ## Structural model
#' sest <- qpeer(formula = y ~ X, excluded.instruments = ~ Z, Glist = G, tau = tau,
#'               structural = TRUE)
#' summary(sest, diagnostic = TRUE)
#' # The lambda^* parameter is y_q (conformity) in the outputs.
#' # There is no conformity in the data, so the estimate will be approximately 0.
#' 
#' ## Structural model with double fixed effects per subnet using optimal GMM 
#' ## and controlling for heteroskedasticity
#' sesto <- qpeer(formula = y ~ X, excluded.instruments = ~ Z, Glist = G, tau = tau,
#'                structural = TRUE, fixed.effects = "separate", HAC = "hetero", 
#'                estimator = "gmm.optimal")
#' summary(sesto, diagnostic = TRUE)
#' 
#' ## Average peer effect model
#' # Row-normalized network to compute instruments
#' Gnorm <- lapply(G, function(g) {
#'   d <- rowSums(g)
#'   d[d == 0] <- 1
#'   g / d
#' })
#' 
#' # GX and GGX
#' Gall <- Matrix::bdiag(Gnorm)
#' GX   <- as.matrix(Gall %*% X)
#' GGX  <- as.matrix(Gall %*% GX)
#' 
#' # Standard linear model
#' lpeer <- linpeer(formula = y ~ X + GX, excluded.instruments = ~ GGX, Glist = Gnorm)
#' summary(lpeer, diagnostic = TRUE)
#' # Note: The normalized network is used here by definition of the model.
#' # Contextual effects are also included (this is also possible for the quantile model).
#' 
#' # The standard model can also be structural
#' lpeers <- linpeer(formula = y ~ X + GX, excluded.instruments = ~ GGX, Glist = Gnorm,
#'                   structural = TRUE, fixed.effects = "separate")
#' summary(lpeers, diagnostic = TRUE)
#' 
#' ## Estimation using `genpeer`
#' # Average peer variable computed manually and included as an endogenous variable
#' Gy     <- as.vector(Gall %*% y)
#' gpeer1 <- genpeer(formula = y ~ X + GX, excluded.instruments = ~ GGX, 
#'                   endogenous.variables = ~ Gy, Glist = Gnorm, structural = TRUE, 
#'                   fixed.effects = "separate")
#' summary(gpeer1, diagnostic = TRUE)
#' 
#' # Using both average peer variables and quantile peer variables as endogenous,
#' # or only the quantile peer variable
#' # Quantile peer `y`
#' qy <- qpeer.inst(formula = y ~ 1, Glist = G, tau = tau)
#' qy <- qy$qy
#' 
#' # Model estimation
#' gpeer2 <- genpeer(formula = y ~ X + GX, excluded.instruments = ~ GGX + Z, 
#'                   endogenous.variables = ~ Gy + qy, Glist = Gnorm, structural = TRUE, 
#'                   fixed.effects = "separate")
#' summary(gpeer2, diagnostic = TRUE)}
#' @importFrom stats pchisq
#' @export
qpeer <- function(formula, excluded.instruments, Glist, tau, type = 7, data, 
                  estimator = "IV", structural = FALSE, fixed.effects = FALSE, 
                  HAC = "iid", checkrank = FALSE, drop = NULL,
                  compute.cov = TRUE, tol = 1e-10){
  # Quantiles
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  ntau       <- length(tau)
  
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
  FEnum <- (0:2)[fixed.effects == c("no", "join", "separate")]
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
  qy         <- fQtauy(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                       ngroup = M, n = n, ntau = ntau, type = type)
  
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
  
  # drop
  if (!is.null(drop)) {
    dg       <- fdrop(drop = drop, ldg = ldg, nvec = nvec, M = M, lIs = lIs, 
                      lnIs = lnIs, y = y, X = X, qy = qy, ins = ins)
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
    qy       <- dg$qy
    ins      <- dg$ins
    dg       <- dg$dg
  }
  
  # Demean fixed effect models
  # save original data
  y0         <- y
  qy0        <- qy
  X0         <- X
  ins0       <- ins
  if (fixed.effects != "no") {
    if (fixed.effects == "join") {
      y      <- c(demean(as.matrix(y), igroup = igr, ngroup = M))
      qy     <- demean(qy, igroup = igr, ngroup = M)
      X      <- demean(X, igroup = igr, ngroup = M)
      ins    <- demean(ins, igroup = igr, ngroup = M)
    } else {
      y      <- c(demean_separate(as.matrix(y), igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n))
      qy     <- demean_separate(qy, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n)
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
    if (length(fcheckrank(X = cbind(qy, X)[nIs + 1,], tol = tol)) != (ntau + Kx)) stop("The design matrix is not full rank.")
  } else {
    if (length(fcheckrank(X = cbind(qy, X), tol = tol)) != (ntau + Kx)) stop("The design matrix is not full rank.")
  }
  
  
  if (structural) {
    ins      <- cbind(X[, idX2 + 1], ins)
    ins0     <- cbind(X0[, idX2 + 1], ins0)
    zename   <- c(xname[idX2 + 1], zename)
    if (checkrank) {
      tlm      <- fcheckrank(X = ins[nIs + 1,], tol = tol)
    }
  } else {
    ins      <- cbind(X, ins)
    ins0     <- cbind(X0, ins0)
    zename   <- c(xname, zename)
    if (checkrank) {
      tlm      <- fcheckrank(X = ins, tol = tol)
    }
  }
  if (checkrank) {
    ins        <- ins[, tlm + 1, drop = FALSE]
    ins0       <- ins0[, tlm + 1, drop = FALSE]
    zename     <- zename[tlm + 1]
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
    if (Kins < Kx2 + ntau) stop("Insufficient number of instruments: the model is not identified.")
    Kest1    <- ifelse(FEnum == 0, Kx1, Kx1 + MIs)
    Kest2    <- ifelse(FEnum == 0, Kx2 + ntau + 1, Kx2 + ntau + MnIs)
    if (length(Is) <= Kest1) stop("Insufficient number of isolated nodes for estimating the structural model.")
    if (length(nIs) <= Kest2) stop("Insufficient number of nonisolated nodes for estimating the structural model.")
    Kest     <- Kest1 + Kest2
    if (HACnum == 2 && (Kx1 >= MIs || Kins + 1 >= MnIs) && estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
      stop("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. HAC is set to 'iid' or 'hetero'.")
    }
    estname  <- c(paste0(yname, paste0("_q", c("(conformity)", 1:ntau))), xname)
    
    # Estimation
    GMMe     <- fstruct(y = y, X = X, qy = qy, ins = ins, idX1 = idX1, idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igr = igr, 
                        nIs = nIs, Is = Is, lnIs = lnIs, lIs = lIs, M = M, MnIs = MnIs, Kins = Kins, Kx = Kx, ntau = ntau, 
                        Kest1 = Kest1, Kest2 = Kest2, n = n, HACnum = HACnum, iv = iv, estimator = estimator, 
                        compute.cov = compute.cov, estname = estname)
  } else {
    if (Kins < Kx + ntau) stop("Insufficient number of instruments: the model is not identified.")
    Kest     <- ifelse(FEnum == 0, Kx + ntau, ifelse(FEnum == 1, Kx + ntau + M, Kx + ntau + MIs + MnIs))
    if (n <= Kest) stop("Insufficient number of observations.")
    if (HACnum == 2 && Kins >= M && estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
      stop("Heteroskedasticity at the group (cluster) level is not possible because the number of groups is small. HAC is set to 'iid' or 'hetero'.")
    }
    estname  <- c(paste0(yname, paste0("_q", 1:ntau)), xname)
    V        <- cbind(qy, X)
    
    # Estimation
    GMMe     <- freduce(y = y, V = V, ins = ins, igr = igr, nvec = nvec, M = M, Kins = Kins, Kx = Kx, ntau = ntau, 
                        Kest = Kest, n = n, HACnum = HACnum, iv = iv, estimator = estimator, compute.cov = compute.cov, 
                        estname = estname)
  }
  if (ntau == 1) {
    colnames(qy0) <- paste0(yname, "_q")
  } else {
    colnames(qy0) <- paste0(yname, "_q", 1:ntau)
  }
  
  out       <- list(model.info  = list(n = n, ngroup = M, nvec = nvec, structural = structural, formula = formula, 
                                       excluded.instruments = excluded.instruments, tau = tau, ntau = ntau, type = type, 
                                       estimator = estimator, fixed.effects = fixed.effects, idXiso = idX1 + 1,  
                                       idXniso = idX2 + 1, HAC = HAC, yname = yname, xnames = xname, znames = zename),
                    gmm         = GMMe,
                    data        = list(y = y0, qy = qy0, X = X0, instruments = ins0, isolated = Is + 1, 
                                       non.isolated = nIs + 1, degree = dg))
  class(out) <- "qpeer"
  out
}

#' @title Summary for the Estimation of Quantile Peer Effects Models
#' @param object An object of class \code{\link{qpeer}}.
#' @param diagnostics,diagnostic Logical. Should diagnostic tests for the instrumental-variable regression be performed? 
#' These include an F-test of the first-stage regression for weak instruments, a Wu-Hausman test for endogeneity, 
#' and a Hansen's J-test for overidentifying restrictions (only if there are more instruments than regressors).
#' @param ... Further arguments passed to or from other methods.
#' @param x An object of class \code{\link{summary.qpeer}} or \code{\link{qpeer}}.
#' @param fullparameters A logical value indicating whether all parameters should be summarized (may be useful for the structural model).
#' @description Summary and print methods for the class \code{\link{qpeer}}.
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, such as the number of subnets, number of observations, and other key details.}
#'     \item{coefficients}{A summary of the estimates, standard errors, and p-values.}
#'     \item{diagnostics}{A summary of the diagnostic tests for the instrumental-variable regression if requested.}
#'     \item{KP.cv}{Critical values for the Kleibergen–Paap Wald test (5% level).}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @export
summary.qpeer <- function(object, fullparameters = TRUE, diagnostic = FALSE, diagnostics = FALSE, ...) {
  stopifnot(inherits(object, "qpeer"))
  if (is.null(object$gmm$cov)) {
    stop("The covariance matrix is not estimated.")
  }
  diagn          <- NULL
  cvKP           <- NULL
  if (diagnostic || diagnostics) {
    diagn        <- fdiagnostic(object, nendo = "qy")
    cvKP         <- diagn$cvKP
    diagn        <- diagn$diag
  }
  if (fullparameters) {
    yname        <- object$model.info$yname
    xnames       <- object$model.info$xnames
    ntau         <- object$model.info$ntau
    est          <- object$gmm$Estimate
    covt         <- object$gmm$cov
    Kx1          <- length(object$model.info$idXiso)
    Kx2          <- length(object$model.info$idXniso)
    if (object$model.info$structural) {
      tp                  <- fStructParamFull(param = est, covp = covt, ntau = ntau, Kx1 = Kx1, Kx2 = Kx2, quantile = 1) 
      tp$theta            <- c(tp$theta)
      names(tp$theta)     <- colnames(tp$Vpa) <- rownames(tp$Vpa) <- c(paste0(yname, paste0("_q", c("(spillover)", "(conformity)", "(total)", 1:ntau))), xnames)
      object$gmm$Estimate <- tp$theta
      object$gmm$cov      <- tp$Vpa
    } else {
      tp                  <- fParamFull(param = est, covp = covt, ntau = ntau, Kx1 = Kx1, Kx2 = Kx2)
      tp$theta            <- c(tp$theta)
      names(tp$theta)     <- colnames(tp$Vpa) <- rownames(tp$Vpa) <- c(paste0(yname, paste0("_q", c("(total)", 1:ntau))), xnames)
      object$gmm$Estimate <- tp$theta
      object$gmm$cov      <- tp$Vpa
    }
  }
  
  coef           <- fcoef(Estimate = object$gmm$Estimate, cov = object$gmm$cov)
  out            <- c(object["model.info"], 
                      list(coefficients = coef, diagnostics = diagn, KP.cv = cvKP),
                      object["gmm"], list(...))
  class(out)     <- "summary.qpeer"
  out
}

#' @rdname summary.qpeer
#' @export
print.summary.qpeer <- function(x, ...) {
  ntau <- x$model.info$ntau
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
      "\nExcluded instruments: ", deparse(x$model.info$excluded.instruments), 
      "\n\nModel: ", ifelse(x$model.info$structural, "Structural", "Reduced Form"),
      "\nEstimator: ", esti,
      "\nFixed effects: ", paste0(toupper(substr(FE, 1, 1)), tolower(substr(FE, 2, nchar(FE)))), "\n", sep = "")
  
  cat("Quantile type: ", x$model.info$type,
      "\nQuantile levels (", ntau, "): ", sep = "")
  cat(x$model.info$tau, "\n")
  
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
  class(x) <- "print.summary.qpeer"
  invisible(x)
}

#' @rdname summary.qpeer
#' @export
print.qpeer <- function(x, ...) {
  print(summary(x))
}

#' @title Simulating Linear Models with Quantile Peer Effects
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\lambda', \beta')'}, where \eqn{\lambda} is a vector of \eqn{\lambda_{\tau}} for each quantile level \eqn{\tau}. 
#' The parameters \eqn{\lambda} and \eqn{\beta} can also be specified separately using the arguments `lambda` and `beta`. For the structural model, 
#' \eqn{\lambda = (\lambda_2, \lambda_{\tau_1}, \lambda_{\tau_2}, \dots)^{\prime}} (see the Details section of \code{\link{qpeer}}).
#' @param lambda The true value of the vector \eqn{\lambda}.
#' @param beta The true value of the vector \eqn{\beta}.
#' @param tau The vector of quantile levels.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see the \link[stats]{quantile} function).
#' @param epsilon A vector of idiosyncratic error terms. If not specified, it will be simulated from a standard normal distribution (see the model specification in the Details section of \code{\link{qpeer}}). 
#' @param maxit The maximum number of iterations for the Fixed Point Iteration Method.
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `sim.qpeer` is called.
#' @param tol The tolerance value used in the Fixed Point Iteration Method to compute the outcome `y`. The process stops if the \eqn{\ell_1}-distance 
#' between two consecutive values of `y` is less than `tol`.
#' @param details A logical value indicating whether to save the indices and weights of the two peers whose weighted average determines the quantile.
#' @param structural A logical value indicating whether simulations should be performed using the structural model. The default is the reduced-form model (see the Details section of \code{\link{qpeer}}).
#' @description
#' `qpeer.sim` simulates the quantile peer effect models developed by Houndetoungan (2025).
#' @seealso \code{\link{qpeer}}, \code{\link{qpeer.instruments}}
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{qy}{Quantiles of the simulated variable among peers.}
#'     \item{epsilon}{The idiosyncratic error.}
#'     \item{index}{The indices of the two peers whose weighted average gives the quantile.}
#'     \item{weight}{The weights of the two peers whose weighted average gives the quantile.}
#'     \item{iteration}{The number of iterations before convergence.}
#' @examples 
#' set.seed(123)
#' ngr  <- 50
#' nvec <- rep(30, ngr)
#' n    <- sum(nvec)
#' G    <- lapply(1:ngr, function(z){
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z])
#'   diag(Gz) <- 0
#'   Gz
#' }) 
#' tau  <- seq(0, 1, 0.25)
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' l    <- c(0.2, 0.1, 0.05, 0.1, 0.2)
#' b    <- c(2, -0.5, 1)
#' 
#' out  <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = l, beta = b)
#' summary(out$y)
#' out$iteration
#' @importFrom stats rnorm
#' @importFrom utils head
#' @importFrom utils tail
#' @export
qpeer.sim <- function(formula, Glist, tau, parms, lambda, beta, epsilon, structural = FALSE, 
                      type = 7, tol = 1e-10, maxit = 500, details = TRUE, data){
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  # Network
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  dg       <- fnetwork(Glist = Glist)
  M        <- dg$M
  nvec     <- dg$nvec
  n        <- dg$n
  igr      <- dg$igr
  Is       <- dg$Is
  nIs      <- dg$nIs
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
  ntau     <- length(tau)
  ltst     <- NULL
  lt       <- NULL
  b        <- NULL
  if (missing(parms)) {
    if (missing(lambda) | missing(beta)) {
      stop("Define either `parms` or both `lambda` and `beta`.")
    }
    if (structural) {
      if (length(lambda) != (ntau + 1)){
        stop("length(lambda) is different from length(tau) + 1. See details on the structural model.")
      }
      ltst <- lambda[1]
      lt   <- tail(lambda, ntau)
    } else {
      if (length(lambda) != ntau){
        stop("length(lambda) is different from length(tau).")
      }
      lt   <- lambda
    }
    if (length(beta) != Kx) stop("length(beta) is different from ncol(X).")
    b      <- beta
  } else{
    if (!missing(lambda) | !missing(beta)) {
      stop("Define either `parms` or both `lambda` and `beta`.")
    }
    if (structural) {
      if (length(parms) != (ntau + Kx + 1)) stop("length(parms) is different from length(tau) + ncol(X) + 1. See details on the structural model.")
      ltst <- parms[1]
      lt   <- parms[2:(ntau + 1)]
    } else {
      if (length(parms) != (ntau + Kx)) stop("length(parms) is different from length(tau) + ncol(X).")
      lt   <- head(parms, ntau)
    }
    b      <- tail(parms, Kx)
  }
  if (sum(abs(lt)) >= 1) warning("The sum of the absolute values of lambda_tau is greater than or equal to one, the Nash Equilibrium may not be stable.")
  if (structural) {
    if (abs(ltst) >= 1) {
      stop("The absolute value of lambda[1] (the parameter that captures whether preferences indicate complementarity/substitution or conformism) must be less than 1.")
    }
  }
  
  # Solving the game
  talpha   <- X %*% b + eps
  if (structural) talpha[nIs + 1] <- talpha[nIs + 1]*(1 - ltst)
  y        <- rep(0, n)
  t        <- fNashE(y = y, G = Glist, d = dg, talpha = talpha, lambdatau = lt, igroup = igr, 
                     nvec = nvec, stau = tau, ngroup = M, n = n, ntau = ntau, type = type, 
                     tol = tol, maxit = maxit)
  # Quantile
  qy       <- NULL
  if (details) {
    qy     <- fQtauyWithIndex(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                              ngroup = M, n = n, ntau = ntau, type = type)
  } else {
    qy     <- fQtauy(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                     ngroup = M, n = n, ntau = ntau, type = type)
  }
  
  W        <- list()
  if (details) {
    W      <- lapply(1:ntau, function(s) fIndexMat(pi1 = qy$pi1[,s], pi2 = qy$pi2[,s], w1 = qy$w1[,s], w2 = qy$w2[,s], n = n))
    qy     <- qy$qy
  }
  
  if (ntau == 1){
    qy            <- c(qy)
    if (details) {
      names(W)    <- "tau"
    }
  } else {
    colnames(qy)  <- paste0("y_q", 1:ntau)
    if (details) {
      names(W)    <- paste0("tau_", 1:ntau)
    }
  }
  
  # Output
  list("y"         = y,
       "qy"        = qy, 
       "epsilon"   = eps,
       "weight"    = W,
       "iteration" = t)
} 


#' @title Specification Tests for Peer Effects Models
#'
#' @param model1,model2 Objects of class \code{\link{qpeer}}, \code{\link{linpeer}}, or \code{\link{genpeer}}.
#' @param which A character string indicating the type of test to be implemented. 
#' The value must be one of `"uniform"`, `"increasing"`, `"decreasing"`, `"wald"`, `"sargan"`, and `"encompassing"` (see Details).
#' @param boot An integer indicating the number of bootstrap replications to use for computing `p-values` in the `"increasing"` and `"decreasing"` tests.
#' @param maxit,eps_f,eps_g Control parameters for the `optim_lbfgs` solver used to optimize the objective function in the `"increasing"` and `"decreasing"` tests (see Kodde and Palm, 1986). 
#' The `optim_lbfgs` function is provided by the \pkg{RcppNumerical} package and is based on the `L-BFGS` method.
#' @param full A Boolean indicating whether the parameters associated with the exogenous variables should be compared. This is used
#' for tests that compare two competing parameter sets.
#' 
#' @references Hayashi, F. (2000). *Econometrics*. Princeton University Press.
#' @references Kodde, D. A., & Palm, F. C. (1986). Wald criteria for jointly testing equality and inequality restrictions. *Econometrica*, 54(5), 1243–1248.  
#'
#' @description
#' `qpeer.test` performs specification tests on peer effects models. These include monotonicity tests on quantile peer effects, as well as tests for instrument validity when an alternative set of instruments is available.
#'
#' @details
#' The monotonicity tests evaluate whether the quantile peer effects \eqn{\lambda_{\tau}} are constant, increasing, or decreasing. In this case, `model1` must be an object of class \code{\link{qpeer}}, and the `which` argument specifies the null hypothesis: `"uniform"`, `"increasing"`, or `"decreasing"`.  
#' For the `"uniform"` test, a standard Wald test is performed. For the `"increasing"` and `"decreasing"` tests, the procedure follows Kodde and Palm (1986). \cr
#'
#' The instrument validity tests assess whether a second set of instruments \eqn{Z_2} is valid, assuming that a baseline set \eqn{Z_1} is valid. In this case, both `model1` and `model2` must be objects of class \code{\link{qpeer}}, \code{\link{linpeer}}, or \code{\link{genpeer}}. 
#' The test compares the estimates obtained using each instrument set. If \eqn{Z_2} nests \eqn{Z_1}, it is recommended to compare the overidentification statistics from both estimations (see Hayashi, 2000, Proposition 3.7).
#' If \eqn{Z_2} does not nest \eqn{Z_1}, the estimates themselves are compared. To compare the overidentification statistics, set the `which` argument to `"sargan"`. To compare the estimates directly, set the `which` argument to `"wald"`.\cr
#' 
#' Given two competing models, it is possible to test whether one is significantly worse using an encompassing test by setting `which` to `"encompassing"`. The null hypothesis is that `model1` is not worse.
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
#'   Gz
#' })
#' 
#' tau <- seq(0, 1, 1/3)
#' X   <- cbind(rnorm(n), rpois(n, 2))
#' l   <- c(0.2, 0.15, 0.1, 0.2)
#' b   <- c(2, -0.5, 1)
#' eps <- rnorm(n, 0, 0.4)
#' 
#' ## Generating `y`
#' y <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = l, 
#'                beta = b, epsilon = eps)$y
#' 
#' ### Estimation
#' ## Computing instruments
#' Z1 <- qpeer.inst(formula = ~ X, Glist = G, tau = seq(0, 1, 0.1), 
#'                  max.distance = 2, checkrank = TRUE)$instruments
#' Z2 <- qpeer.inst(formula = y ~ X, Glist = G, tau = seq(0, 1, 0.1), 
#'                  max.distance = 2, checkrank = TRUE)$instruments
#' 
#' ## Reduced-form model 
#' rest1 <- qpeer(formula = y ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau)
#' summary(rest1, diagnostic = TRUE)  
#' rest2 <- qpeer(formula = y ~ X, excluded.instruments = ~ Z1 + Z2, Glist = G, tau = tau)
#' summary(rest2, diagnostic = TRUE)  
#' 
#' qpeer.test(model1 = rest1, which = "increasing")
#' qpeer.test(model1 = rest1, which = "decreasing")
#' qpeer.test(model1 = rest1, model2 = rest2, which = "sargan")
#' 
#' #' A model with a mispecified tau
#' rest3 <- qpeer(formula = y ~ X, excluded.instruments = ~ Z1 + Z2, Glist = G, tau = c(0, 1))
#' summary(rest3)
#' #' Test is rest3 is worse than rest1
#' qpeer.test(model1 = rest3, model2 = rest1, which = "encompassing")
#' }
#' @importFrom MASS ginv
#' @importFrom stats qchisq
#' @importFrom stats uniroot
#' @export
qpeer.test <- function(model1, model2 = NULL, which, full = FALSE,
                       boot = 1e4, maxit = 1e6, eps_f = 1e-9, eps_g = 1e-9) {
  which  <- tolower(which[1])
  stopifnot(which %in% c("increasing", "decreasing", "uniform", "wald", "sargan", "encompassing"))
  
  struc1 <- model1$model.info$structural
  struc2 <- model2$model.info$structural
  qy1    <- model1$data$qy
  qy2    <- model2$data$qy
  X1     <- model1$data$X
  X2     <- model2$data$X
  Z1     <- model1$data$instruments
  Z2     <- model2$data$instruments
  W11    <- model1$gmm$Wiso
  W21    <- model1$gmm$Wniso
  W1     <- model1$gmm$W
  W12    <- model2$gmm$Wiso
  W22    <- model2$gmm$Wniso
  W2     <- model2$gmm$W
  e1     <- model1$gmm$residuals
  e2     <- model2$gmm$residuals
  theta1 <- model1$gmm$Estimate
  theta2 <- model2$gmm$Estimate
  HAC1   <- (0:2)[model1$model.info$HAC == c("iid", "hetero", "cluster")]
  HAC2   <- (0:2)[model2$model.info$HAC == c("iid", "hetero", "cluster")]  
  idX11  <- model1$model.info$idXiso - 1
  idX21  <- model2$model.info$idXniso - 1
  idX12  <- model1$model.info$idXiso - 1
  idX22  <- model2$model.info$idXniso - 1
  nIs1   <- model1$data$non.isolated - 1
  Is1    <- model1$data$isolated - 1
  nIs2   <- model2$data$non.isolated - 1
  Is2    <- model2$data$isolated - 1
  ngr1   <- model1$model.info$ngroup
  ngr2   <- model2$model.info$ngroup
  FE1    <- model1$model.info$fixed.effects
  FE2    <- model2$model.info$fixed.effects
  nvec1  <- model1$model.info$nvec
  nvec2  <- model2$model.info$nvec
  n1     <- model1$model.info$n
  n2     <- model2$model.info$n
  tau1   <- model1$model.info$tau
  tau2   <- model2$model.info$tau
  
  Tval   <- NULL
  op     <- NULL
  lt1    <- NULL
  lt2    <- NULL
  dtheta <- NULL
  if (which %in% c("increasing", "decreasing", "uniform")) {
    stopifnot(class(model1) == "qpeer")
    
    if (!is.null(model2)) {
      stop('Testing whether peer effects are "increasing", "decreasing", or "uniform" is only possible on `model1`.')
    }
    
    ntau <- model1$model.info$ntau
    if (ntau == 1) {
      stop("There is only one peer effect parameter.")
    }
    
    struc <- model1$model.info$structural
    lt1  <- model1$gmm$Estimate[(struc + 1):(struc + ntau)]
    
    if (is.null(model1$gmm$cov)) {
      stop("The covariance matrix has not been estimated.")
    }
    
    CThe   <- model1$gmm$cov[(struc + 1):(struc + ntau), (struc + 1):(struc + ntau)]
    if (which == "uniform") {
      R    <- t(diag(1, ntau, ntau - 1))
      R    <- R - cbind(0, R[, 1:(ntau - 1)])
      stat <- c(t(R %*%  lt1) %*% solve(R %*% CThe %*% t(R), R %*%  lt1))
      df   <- ntau - 1
      pval <- pchisq(stat, df, lower.tail = FALSE)
      Tval <- c("statistic" = stat, "df" = df, "p-value" = pval)
    } else {
      tp   <- NULL
      ts   <- t(chol(CThe)) %*% matrix(rnorm(ntau * boot), ntau, boot) +  lt1
      suppressWarnings(
        tp <- fTestMonotone(thetahat =  lt1, Sigma = CThe, a = rep(ifelse(which == "increasing", 1, -1), ntau - 1), 
                            thetasimu = ts, Boot = boot, maxit = maxit,
                            eps_f = eps_f, eps_g = eps_g)
      )
      stat <- tp$optim$minimum
      op   <- tp$optim
      Pi   <- sapply(0:(ntau - 1), function(s) mean(tp$count == s))
      pval <- uniroot(function(p){
        qq   <- qchisq(p, df = (ntau - 1):0, lower.tail = FALSE)
        qq[qq == Inf] <- 1e10
        mean(qq*Pi) - stat
      }, interval = c(0, 1), tol = eps_f)$root
      Tval  <- c("statistic" = stat, "p-value" = pval)
    }
  } else if (which %in% c("wald", "cox", "encompassing")) {
    stopifnot(class(model1) %in% c("qpeer", "linpeer", "genpeer"))
    stopifnot(class(model2) %in% c("qpeer", "linpeer", "genpeer"))
    if (is.null(model1$gmm$cov) | is.null(model2$gmm$cov)) {
      stop("The covariance matrix has not been estimated.")
    }
    
    if (!(tolower(model1$model.info$estimator) %in% c("iv", "gmm.optimal", "gmm.identity") &&
          tolower(model2$model.info$estimator) %in% c("iv", "gmm.optimal", "gmm.identity"))) {
      stop("This test is only valid for IV and GMM estimators.")
    }
    
    tp     <- (struc1 == struc2) &&
      (length(e1) == length(e2)) &&
      (HAC1 == HAC2) &&
      all(nIs1 == nIs2) &&
      all(Is1 == Is2) &&
      (ngr1 == ngr2) &&
      (FE1 == FE2) &&
      all(nvec1 == nvec2) &&
      (n1 == n2)
    
    if (!tp) {
      stop("`model1` and `model2` differ in data or model specification.")
    }
    igr     <- matrix(c(cumsum(c(0, nvec1[-ngr1])), cumsum(nvec1) - 1), ncol = 2)
    ncs     <- c(0, cumsum(nvec1))
    LIs     <- lapply(1:ngr1, function(m) {
      Is1[Is1 %in% ncs[m]:(ncs[m + 1] - 1)]
    })
    LnIs    <- lapply(1:ngr1, function(m) {
      nIs1[nIs1 %in% ncs[m]:(ncs[m + 1] - 1)]
    })
    FEnum   <- (0:2)[FE1 == c("no", "join", "separate")]
    
    if (FE1 == "join") {
      qy1  <- demean(qy1, igr, ngr1)
      qy2  <- demean(qy2, igr, ngr1)
      X1   <- demean(X1, igr, ngr1)
      X2   <- demean(X2, igr, ngr1)
      Z1   <- demean(Z1, igr, ngr1)
      Z2   <- demean(Z2, igr, ngr1)
    } else if(FE1 == "separate") {
      qy1  <- demean_separate(qy1, igr, LIs, LnIs, ngr1, n1)
      qy2  <- demean_separate(qy2, igr, LIs, LnIs, ngr1, n1)
      X1   <- demean_separate(X1, igr, LIs, LnIs, ngr1, n1)
      X2   <- demean_separate(X2, igr, LIs, LnIs, ngr1, n1)
      Z1   <- demean_separate(Z1, igr, LIs, LnIs, ngr1, n1)
      Z2   <- demean_separate(Z2, igr, LIs, LnIs, ngr1, n1)
    }
    
    MIs    <- sum(sapply(LIs, function(s) length(s) > 0))
    MnIs   <- sum(sapply(LnIs, function(s) length(s) > 0))
    
    if (which == "wald") {
      tp     <- (length(qy1) == length(qy2)) &&
        (length(X1) == length(X2)) &&
        (length(W11) == length(W12)) &&
        (length(theta1) == length(theta2)) &&
        (length(idX11) == length(idX12)) &&
        (length(idX21) == length(idX22)) &&
        all(colnames(qy1) == colnames(qy2)) &&
        all(colnames(X1) == colnames(X2)) &&
        all(tau1 == tau2)
      
      if (!tp) {
        stop("`model1` and `model2` differ in data or model specification.")
      }
      
      CTT    <- NULL
      ntau   <- length(tau1)
      if (struc1) {
        K1   <- length(idX11)
        K2   <- length(idX21)
        Kest1 <- ifelse(FEnum == 0, K1, K1 + MIs)
        Kest2 <- ifelse(FEnum == 0, K2 + ntau + 1, K2 + ntau + MnIs)
        CTT   <- Cov2ThetaStruc(X1 = X1, qy1 = qy1, Z1 = Z1, W11 = W11, W21 = W21, e1 = e1, theta1 = theta1, 
                                Kest11 = Kest1, Kest21 = Kest2, X2 = X2, qy2 = qy2, Z2 = Z2, W12 = W12, W22 = W22,
                                e2 = e2, theta2 = theta2, Kest12 = Kest1, Kest22 = Kest2, idX1 = idX11, idX2 = idX21, 
                                nIs = nIs1, Is = Is1, ngroup = ngr1, cumsn = ncs, HAC = HAC1, full = full)
      } else {
        K     <- length(theta1)
        Kest  <- ifelse(FEnum == 0, K, ifelse(FEnum == 1, K + ngr1, K + MIs + MnIs))
        CTT   <- Cov2ThetaRed(X1 = X1, qy1 = qy1, Z1 = Z1, W1 = W1, e1 = e1, theta1 = theta1, Kest1 = Kest,
                              X2 = X2, qy2 = qy2, Z2 = Z2, W2 = W2, e2 = e2, theta2 = theta2, Kest2 = Kest, 
                              ngroup = ngr1, cumsn = ncs, HAC = HAC1, full = full)
      }
      stat    <- c(CTT$stat)
      df      <- c(CTT$df)
      pval    <- pchisq(stat, df, lower.tail = FALSE)
      Tval    <- c("statistic" = stat, "df" = df, "p-value" = pval)
      lt1     <- c(CTT$theta1)
      lt2     <- c(CTT$theta2)
    } else if (which == "encompassing") {
      ntau1   <- length(tau1)
      ntau2   <- length(tau2)
      CTT     <- NULL
      if (struc1) {
        K11   <- length(idX11)
        K21   <- length(idX21)
        K12   <- length(idX12)
        K22   <- length(idX22)
        Kest11<- ifelse(FEnum == 0, K11, K11 + MIs)
        Kest21<- ifelse(FEnum == 0, K21 + ntau1 + 1, K21 + ntau1 + MnIs)
        Kest12<- ifelse(FEnum == 0, K12, K12 + MIs)
        Kest22<- ifelse(FEnum == 0, K22 + ntau2 + 1, K22 + ntau2 + MnIs)
        CTT   <- list(F = fEncompassingStruc(X1 = X1, qy1 = qy1, Z1 = Z1, W11 = W11, W21 = W21, e1 = e1, theta1 = theta1,
                                             idX11 = idX11, idX21 = idX21, Kest11 = Kest11, Kest21 = Kest21,
                                             X2 = X2, qy2 = qy2, Z2 = Z2, W12 = W12, W22 = W22, e2 = e2, theta2 = theta2,
                                             idX12 = idX12, idX22 = idX22, Kest12 = Kest12, Kest22 = Kest22,
                                             nIs = nIs1, Is = Is1, ngroup = ngr1, cumsn = ncs, HAC = HAC1, full = full), 
                      
                      KP = fEncompassingStrucKP(X1 = X1, qy1 = qy1, Z1 = Z1, W11 = W11, W21 = W21, e1 = e1, theta1 = theta1,
                                                idX11 = idX11, idX21 = idX21, Kest11 = Kest11, Kest21 = Kest21,
                                                X2 = X2, qy2 = qy2, Z2 = Z2, W12 = W12, W22 = W22, e2 = e2, theta2 = theta2,
                                                idX12 = idX12, idX22 = idX22, Kest12 = Kest12, Kest22 = Kest22,
                                                nIs = nIs1, Is = Is1, ngroup = ngr1, cumsn = ncs, HAC = HAC1, full = full))
      } else {
        K1    <- length(theta1)
        K2    <- length(theta2)
        Kest1 <- ifelse(FEnum == 0, K1, ifelse(FEnum == 1, K1 + ngr1, K1 + MIs + MnIs))
        Kest2 <- ifelse(FEnum == 0, K2, ifelse(FEnum == 1, K2 + ngr1, K2 + MIs + MnIs))
        CTT   <- list(F = fEncompassingRed(X1 = X1, qy1 = qy1, Z1 = Z1, W1 = W1, e1 = e1, theta1 = theta1, Kest1 = Kest1,
                                  X2 = X2, qy2 = qy2, Z2 = Z2, W2 = W2, e2 = e2, theta2 = theta2, Kest2 = Kest2,
                                  ngroup = ngr1, cumsn = ncs, HAC = HAC1, full = full),
                      KP = fEncompassingRedKP(X1 = X1, qy1 = qy1, Z1 = Z1, W1 = W1, e1 = e1, theta1 = theta1, Kest1 = Kest1,
                                            X2 = X2, qy2 = qy2, Z2 = Z2, W2 = W2, e2 = e2, theta2 = theta2, Kest2 = Kest2,
                                            ngroup = ngr1, cumsn = ncs, HAC = HAC1, full = full))
      }
      Fstat   <- c(CTT$F$stat)
      Fdf1    <- c(CTT$F$df1)
      Fdf2    <- c(CTT$F$df2)
      Fpval   <- pf(Fstat, Fdf1, Fdf2, lower.tail = FALSE)
      KPstat  <- c(CTT$KP$stat)
      KPdf    <- c(CTT$KP$df)
      KPpval  <- pchisq(KPstat, KPdf, lower.tail = FALSE)
      Tval    <- cbind("statistic" = c(Fstat, KPstat), 
                       "df1"       = c(Fdf1, KPdf), 
                       "df2"       = c(Fdf2, NA), 
                       "p-value"   = c(Fpval, KPpval))
      rownames(Tval) <- c("robust F", "KP Wald rank")
      dtheta         <- c(CTT$diff)
    } else if (which == "cox") { #TO not available yet
    }
    
  } else if (which == "sargan") {
    stopifnot(class(model1) %in% c("qpeer", "linpeer", "genpeer"))
    stopifnot(class(model2) %in% c("qpeer", "linpeer", "genpeer"))
    if (!(tolower(model1$model.info$estimator) %in% c("iv", "gmm.optimal", "gmm.identity") &&
          tolower(model2$model.info$estimator) %in% c("iv", "gmm.optimal", "gmm.identity"))) {
      stop("This test is only valid for IV and GMM estimators.")
    }
    df     <- model2$gmm$Jtest[2] - model1$gmm$Jtest[2]
    if (df <= 0) {
      stop("sargan test required instruments in model2 to neast instruments in model1")
    }
    stat   <- model2$gmm$Jtest[1] - model1$gmm$Jtest[1]
    pval   <- pchisq(stat, df, lower.tail = FALSE)
    Tval  <- c("statistic" = stat, "df" = df, "p-value" = pval)
    struc  <- model1$model.info$structural
    ntau   <- model1$model.info$ntau
    lt1    <- model1$gmm$Estimate[1:(ntau + struc)]
    lt2    <- model2$gmm$Estimate[1:(ntau + struc)]
  }
out         <- list("test" = Tval, 
                    "lambda1" = lt1, "lambda2" = lt2,
                    which = which, boot = boot)
class(out)  <- "qpeer.test"
out
}

#' @title Printing Specification Tests for Peer Effects Models
#' @param x an object of class \code{\link{qpeer.test}}
#' @param ... Further arguments passed to or from other methods.
#' @description A print method for the class \code{\link{qpeer.test}}.
#' @export
print.qpeer.test <- function(x, ...) {
  stopifnot(class(x) == "qpeer.test")
  if (x$which %in% c("increasing", "decreasing", "uniform")) {
    cat("Testing quantile peer effect monotonicity\n\n")
    cat("Quantile peer effects:\n")
    cat("  lambda_tau:", x$lambda1, "\n")
    cat("  Null hypothesis: lambda_tau is", x$which, "\n")
    cat("  Statistic:", x$test["statistic"])
    cat(" -- p-value:", ifelse(x$test["p-value"] < 2e-16, "< 2e-16", format(x$test["p-value"], digits = 4)), "\n")
  } else if (x$which %in% c("wald", "sargan")) {
    cat("Testing instrument validity (", ifelse(x$which == "wald", "Wald", "J-Sargan"), " Test)", "\n\n", sep = "")
    cat("Quantile peer effects:\n")
    cat("  Model 1 (Z1):", x$lambda1, "\n")
    cat("  Model 2 (Z2):", x$lambda2, "\n")
    cat("  Null hypothesis: Z2 is exogenous\n")
    cat("  Statistic:", x$test["statistic"])
    cat(" -- p-value:", ifelse(x$test["p-value"] < 2e-16, "< 2e-16", format(x$test["p-value"], digits = 4)), "\n")
  } else if (x$which == "encompassing") {
    cat("Encompassing Test\n")
    cat("  Null hypothesis: Model 1 is not worse\n")
    cat("  Robust F statistic:", x$test["robust F", "statistic"])
    cat(" -- p-value:", ifelse(x$test["robust F", "p-value"] < 2e-16, "< 2e-16", format(x$test["robust F", "p-value"], digits = 4)), "\n")
    cat("  KP Wald rank statistic:", x$test["KP Wald rank", "statistic"])
    cat(" -- p-value:", ifelse(x$test["KP Wald rank", "p-value"] < 2e-16, "< 2e-16", format(x$test["KP Wald rank", "p-value"], digits = 4)), "\n")
  }
  invisible(x)
}
