#' @title Simulating Linear Models with Quantile Peer Effects
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\lambda', \beta')'}, where \eqn{\lambda} is a vector of \eqn{\lambda_{\tau}} for each quantile level \eqn{\tau} (see model specification in the details). 
#' The parameters \eqn{\lambda} and \eqn{\beta} can also be specified separately using the arguments `lambda` and `beta`. For the structural model, 
#' \eqn{\lambda = (\lambda^{*}, \lambda_{\tau_1}, \lambda_{\tau_2}, \dots)^{\prime}} (see details).
#' @param lambda The true value of the vector \eqn{\lambda}.
#' @param beta The true value of the vector \eqn{\beta}.
#' @param tau The vector of quantile levels.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see the \link[stats]{quantile} function).
#' @param epsilon A vector of idiosyncratic error terms. If not specified, it will be simulated from a standard normal distribution (see the model specification in details). 
#' @param maxit The maximum number of iterations for the Fixed Point Iteration Method.
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `qpeer.sim` is called.
#' @param tol The tolerance value used in the Fixed Point Iteration Method to compute the outcome `y`. The process stops if the \eqn{\ell_1}-distance 
#' between two consecutive values of `y` is less than `tol`.
#' @param details A logical value indicating whether to save the indices and weights of the two peers whose weighted average determines the quantile.
#' @param structural A logical value indicating whether simulations should be performed using the structural model. The default is the reduced-form model (see details).
#' @description
#' `qpeer.sim` simulates the quantile peer effect models developed by Houndetoungan (2025).
#' @details 
#' Let \eqn{\mathcal{T}} be a set of quantile levels. The reduced-form specification of quantile peer effect models is given by:
#' \deqn{y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + \mathbf{x}_i^{\prime}\beta + \varepsilon_i,}
#' where \eqn{\mathbf{y}_{-i} = (y_1, \ldots, y_{i-1}, y_{i+1}, \ldots, y_n)^{\prime}} is the vector of outcomes for other units, and \eqn{q_{\tau,i}(\mathbf{y}_{-i})} is the 
#' sample \eqn{\tau}-quantile of peer outcomes. The term \eqn{\varepsilon_i} is an idiosyncratic error term, \eqn{\lambda_{\tau}} captures the effect of the \eqn{\tau}-quantile of peer outcomes on \eqn{y_i}, 
#' and \eqn{\beta} captures the effect of \eqn{\mathbf{x}_i} on \eqn{y_i}. For the definition of the sample \eqn{\tau}-quantile, see Hyndman and Fan (1996). The network matrices in `Glist` can be weighted or unweighted. 
#' If weighted, the sample weighted quantile is computed, where the outcome for friend \eqn{j} of \eqn{i} is weighted by \eqn{g_{ij}}, the \eqn{(i, j)} entry of the network matrix. It can be shown that
#' the sample \eqn{\tau}-quantile is a weighted average of two peer outcomes. For more details, see the \link[stats]{quantile} and \code{\link{qpeer.instruments}} functions. \cr
#' 
#' The specification of the structural model depends on whether node \eqn{i} is isolated or not. For isolated \eqn{i}, the specification is similar to a standard linear-in-means model without social interactions, given by:
#' \deqn{y_i = \mathbf{x}_i^{\prime}\beta + \varepsilon_i.}
#' If node \eqn{i} is non-isolated, the specification is:
#' \deqn{y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + (1 - \lambda^*)\mathbf{x}_i^{\prime}\beta  + \varepsilon_i,}
#' where \eqn{\lambda^*} captures whether preferences describe complementarity/substitution or conformism.
#' @seealso \code{\link{qpeer.estim}}, \code{\link{qpeer.instruments}}
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{qy}{Quantiles of the simulated variable among peers.}
#'     \item{epsilon}{The idiosyncratic error.}
#'     \item{index}{The indices of the two peers whose weighted average gives the quantile.}
#'     \item{weight}{The weights of the two peers whose weighted average gives the quantile.}
#'     \item{iteration}{Number of iterations performed by the sub-network in the Fixed Point Iteration Method.}
#' @examples 
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
      if (length(parms) != (ntau + Kx)) stop("length(parms) is different from length(tau) + ncol(X) + 1. See details on the structural model.")
      ltst <- parms[1]
      lt   <- parms[2:(ntau + 1)]
    } else {
      if (length(parms) != (ntau + Kx)) stop("length(parms) is different from length(tau) + ncol(X).")
      lt   <- head(parms, ntau)
    }
    b      <- tail(parms, Kx)
  }
  if (sum(abs(lt)) >= 1) warning("The sum of the absolute values of lambda_tau is greater than or equal to one, the Nash Equilibrium may not be stable")
  if (structural) {
    if (abs(ltst) >= 1) {
      stop("The absolute value of lambda[1] (the parameter that captures whether preferences indicate complementarity/substitution or conformism) must be less than 1.")
    }
  }
  
  # Solving the game
  talpha   <- X %*% b
  if (structural) talpha[nIs + 1] <- talpha[nIs + 1]*(1 - ltst)
  talpha   <- talpha + eps
  y        <- rep(0, n)
  t        <- fNashE(y = y, G = Glist, d = dg, talpha, lambdatau = lt, igroup = igr, 
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
  
  i        <- NULL
  w        <- NULL
  if (details) {
    i      <- list(pi1 = qy$pi1, pi2 = qy$pi2)
    w      <- list(w1 = qy$w1, w2 = qy$w2)
    qy     <- qy$qy
  }
  
  if (ntau == 1){
    qy       <- c(qy)
    if (details) {
      i$pi1  <- c(i$pi1)
      i$pi2  <- c(i$pi2)
      w$w1   <- c(w$w1)
      w$w2   <- c(w$w2)
    }
  } else {
    colnames(qy)      <- paste0("y_q", 1:ntau)
    if (details) {
      colnames(i$pi1) <- paste0("pi1_q", 1:ntau)
      colnames(i$pi2) <- paste0("pi2_q", 1:ntau)
      colnames(w$w1)  <- paste0("w1_q", 1:ntau)
      colnames(w$w2)  <- paste0("w2_q", 1:ntau)
    }
  }
  
  # Output
  list("y"         = y,
       "qy"        = qy, 
       "epsilon"   = eps,
       "index"     = i,
       "weight"    = w,
       "iteration" = t)
} 

#' @importFrom gmm gmm
#' @export
qpeer.estim <- function(formula, instruments, Glist, tau, type = 7, data, optimal.instruments = FALSE, 
                     gmm.weight = "IV", structural = FALSE, fixed.effects = FALSE, tol = 1e-10, maxit = 500, ...){
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  ntau       <- length(tau)
  gmm.weight <- tolower(gmm.weight)
  stopifnot(gmm.weight %in% c("iv", "optimal", "ident"))
  
  if (is.character(fixed.effects[1])) fixed.effects <- tolower(fixed.effects)
  stopifnot(fixed.effects %in% c(FALSE, "no", TRUE, "yes", "join", "separate"))
  if (fixed.effects == TRUE | fixed.effects == "yes") fixed.effects <- ifelse(structural, "separate", "join")
  if (fixed.effects == FALSE) fixed.effects <- "no"
  
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
  if (length(Is) <= 1 & structural) stop("The structural model requires isolated nodes.")
  
  # linf and lsup
  linf     <- NULL
  lsup     <- NULL
  if (structural) {
    linf   <- rep(0, ntau + 1)
    lsup   <- rep(1, ntau + 1)
  }
  
  # Data
  # y and X
  formula    <- as.formula(formula)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = (fixed.effects != "no"), 
                                simulations = FALSE) 
  y          <- f.t.data$y
  X          <- f.t.data$X
  Kx         <- ncol(X)
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  xint       <- f.t.data$intercept
  qy         <- fQtauy(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                       ngroup = M, n = n, ntau = ntau, type = type)
  
  # Instruments
  inst       <- as.formula(instruments)
  if(length(inst) != 2) stop("Expected format for instruments is ~ X1 + X2 + ....")
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
  X0         <- X #Save initial X
  if (fixed.effects != "no") {
    if (fixed.effects == "join") {
      y      <- c(demean(as.matrix(y), igroup = igr, ngroup = M))
      qy     <- demean(qy, igroup = igr, ngroup = M)
      X      <- demean(X, igroup = igr, ngroup = M)
      ins    <- demean(ins, igroup = igr, ngroup = M)
    } else {
      y      <- c(demean_separate(as.matrix(y), igroup = igr, Is = Is, ngroup = M, n = n))
      qy     <- demean_separate(qy, igroup = igr, Is = Is, ngroup = M, n = n)
      X      <- demean_separate(X, igroup = igr, Is = Is, ngroup = M, n = n)
      ins    <- demean_separate(ins, igroup = igr, Is = Is, ngroup = M, n = n)
    }
  }
  
  # Construct the final instrument matrix
  XStruc     <- NULL
  if (structural) {
    Xiso     <- X; Xiso[nIs + 1, ] <- 0
    Xniso    <- X; Xniso[Is + 1, ] <- 0
    XStruc   <- cbind(Xiso, Xniso)
    XStruc   <- XStruc[, fcheckrank(XStruc) == 1, drop = FALSE]
    ins      <- cbind(XStruc, ins)
  } else {
    ins      <- cbind(X, ins)
  }
  
  
  
  # GMM weight
  weight     <- gmm.weight
  if (gmm.weight == "iv"){
    weight   <- solve(crossprod(ins)/n)
  }
  
  # GMM
  GMM1      <- NULL
  if (structural) {
    x       <- list(y = y, qy = qy, X = X, ins = ins, Is = Is, nIs = nIs, 
                    linf = linf, lsup = lsup, n = n, Kx = Kx, ntau = ntau, Kins = ncol(ins))
    if (gmm.weight == "iv") {
      GMM1  <- gmm(g, x = x, grad = dg, weightsMatrix = weight, 
                   t0 = rep(0, ntau + Kx + 1), ...)
    } else {
      GMM1  <- gmm(g, x = x, grad = dg, wmatrix = weight, 
                   t0 = rep(0, ntau + Kx + 1), ...)
    }
  } else {
    if (gmm.weight == "iv") {
      GMM1  <- gmm(y ~ -1 + qy + X, x = ins, weightsMatrix = weight)
    } else {
      GMM1  <- gmm(y ~ -1 + qy + X, x = ins, wmatrix = weight, ...)
    }
    if (ntau == 1) {
      names(GMM1$coefficients) <- c(paste0(yname, "_q"), xname)
    } else {
      names(GMM1$coefficients) <- c(paste0(yname, "_q", 1:ntau), xname)
    }
  }
  
  # Using the optimal instrument
  GMM2      <- list()
  if (optimal.instruments) {
    if (structural) {
      lt       <- flambda(lambdatilde = head(GMM1$coefficients, ntau + 1), linf = linf, lsup = lsup, ntau = ntau)
      ltst     <- lt[1]
      lt       <- lt[-1]
      b        <- tail(GMM1$coefficients, Kx)
      talpha   <- X0 %*% b
      talpha[nIs + 1] <- talpha[nIs + 1]*(1 - ltst)
      Ey       <- rep(0, n)
      t        <- fNashE(y = Ey, G = Glist, d = dg, talpha = talpha, lambdatau = lt, igroup = igr,
                        nvec = nvec, stau = tau, ngroup = M, n = n, ntau = ntau, type = type,
                        tol = tol, maxit = maxit)
      ins      <- fQtauy(y = Ey, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                        ngroup = M, n = n, ntau = ntau, type = type)
      if (fixed.effects != "no") {
        if (fixed.effects == "join") {
          ins   <- demean(ins, igroup = igr, ngroup = M)
        } else {
          ins   <- demean_separate(ins, igroup = igr, Is = Is, ngroup = M, n = n)
        }
      }
      ins      <- cbind(XStruc, ins)
      # GMM weight
      weight   <- gmm.weight
      if (gmm.weight == "iv"){
        weight <- solve(crossprod(ins))
      }
      
      # Estimation
      x       <- list(y = y, qy = qy, X = X, ins = ins, Is = Is, nIs = nIs, 
                      linf = linf, lsup = lsup, n = n, Kx = Kx, ntau = ntau, Kins = ncol(ins))
      if (gmm.weight == "iv") {
        GMM2  <- gmm(g, x = x, grad = dg, weightsMatrix = weight, 
                     t0 = GMM1$coefficients, ...)
      } else {
        GMM2  <- gmm(g, x = x, grad = dg, wmatrix = weight, 
                     t0 = GMM1$coefficients, ...)
      }
    } else {
      lt      <- head(GMM1$coefficients, ntau)
      b       <- tail(GMM1$coefficients, Kx)
      talpha  <- X0 %*% b
      Ey      <- rep(0, n)
      t       <- fNashE(y = Ey, G = Glist, d = dg, talpha, lambdatau = lt, igroup = igr,
                        nvec = nvec, stau = tau, ngroup = M, n = n, ntau = ntau, type = type,
                        tol = tol, maxit = maxit)
      ins     <- fQtauy(y = Ey, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                        ngroup = M, n = n, ntau = ntau, type = type)
      if (fixed.effects != "no") {
        if (fixed.effects == "join") {
          ins   <- demean(ins, igroup = igr, ngroup = M)
        } else {
          ins   <- demean_separate(ins, igroup = igr, Is = Is, ngroup = M, n = n)
        }
      }
      ins     <- cbind(X, ins)
      # GMM weight
      weight     <- gmm.weight
      if (gmm.weight == "iv"){
        weight   <- solve(crossprod(ins))
      }
      if (gmm.weight == "iv") {
        GMM2  <- gmm(y ~ -1 + qy + X, x = ins, weightsMatrix = weight, ...)
      } else {
        GMM2  <- gmm(y ~ -1 + qy + X, x = ins, wmatrix = weight, ...)
      }
      if (ntau == 1) {
        names(GMM2$coefficients) <- c(paste0(yname, "_q"), xname)
      } else {
        names(GMM2$coefficients) <- c(paste0(yname, "_q", 1:ntau), xname)
      }
    }
  }
  
  out       <- list(model.info  = list(n = n, ngroup = M, tau = tau, formula = formula, instruments = instruments,
                                      type = type, gmm.weight = gmm.weight, optimal.instruments = optimal.instruments,
                                      fixed.effects = fixed.effects),
                    gmm         = GMM1,
                    gmm.opt.ins = GMM2)
  class(out) <- c("qpeer.estim")
  out
}

