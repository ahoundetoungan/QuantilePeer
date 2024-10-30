#' @title Simulating Linear Models with Quantile Peer Effects
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\lambda', \beta')'}, where \eqn{\lambda} is a vector of \eqn{\lambda_{\tau}} for each quantile level \eqn{\tau} (see model specification in the details). 
#' The parameters \eqn{\lambda} and \eqn{\beta} can also be specified separately using the arguments `lambda` and `beta`.
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
#' @description
#' `qpeer.sim` simulates the quantile peer effect models developed by Houndetoungan (2025).
#' @details 
#' Let \eqn{\mathcal{T}} be a set of quantile levels. The quantile peer effect model is given by:
#' \deqn{y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + \mathbf{x}_i^{\prime}\beta + \varepsilon_i,}
#' where \eqn{\mathbf{y}_{-i} = (y_1, \ldots, y_{i-1}, y_{i+1}, \ldots, y_n)^{\prime}} is the vector of outcomes for other units and \eqn{q_{\tau,i}(\mathbf{y}_{-i})} is the 
#' sample \eqn{\tau}-quantile of peer outcomes. The term \eqn{\varepsilon_i} is an idiosyncratic error term, \eqn{\lambda_{\tau}} captures the effect of the \eqn{\tau}-quantile of peer outcomes on \eqn{y_i}, 
#' and \eqn{\beta} captures the effect of \eqn{\mathbf{x}_i} on \eqn{y_i}. For the definition of the sample \eqn{\tau}-quantile, see Hyndman and Fan (1996). The network matrices in `Glist` can be weighted or unweighted. 
#' If weighted, the sample weighted quantile is computed, where the outcome for friend \eqn{j} of \eqn{i} is weighted by \eqn{g_{ij}}, the \eqn{(i, j)} entry of the network matrix. It can be shown that
#' the sample \eqn{\tau}-quantile is a weighted average of two peer outcomes. For more details, see the \link[stats]{quantile} and \code{\link{quantvars}} functions.
#' @seealso \code{\link{qpeer.lim}}, \code{\link{qpeer.instruments}}
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{Qy}{Quantiles of the simulated variable among peers.}
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
qpeer.sim <- function(formula, Glist, tau, parms, lambda, beta, epsilon, type = 7, tol = 1e-10, maxit = 500, details = TRUE, data){
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  # Network
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  dg       <- unlist(lapply(Glist, rowSums))
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
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
  }
  
  # parameters
  ntau     <- length(tau)
  lt       <- NULL
  b        <- NULL
  if (missing(parms)) {
    if (missing(lambda) | missing(beta)) {
      stop("To compute optimal instruments, please define either `parms` or both `lambda` and `beta`.")
    }
    lt     <- lambda
    b      <- beta
    if (length(lt) != ntau) stop("length(lambda) is different from length(tau)")
    if (length(b) != Kx) stop("length(beta) is different from ncol(X)")
  } else{
    if (!missing(lambda) | !missing(beta)) {
      stop("Define either `parms` or both `lambda` and `beta`.")
    }
    if (length(parms) != (ntau + Kx)) stop("length(parms) is different from length(tau) + ncol(X)")
    lt     <- head(parms, ntau)
    b      <- tail(parms, Kx)
  }
  if (sum(abs(lt)) >= 1) warning("The sum of the absolute values of lambda_tau is greater than or equal to one, the Nash Equilibrium may not be stable")
  
  # Solving the game
  talpha   <- X %*% b + eps
  y        <- rep(0, n)
  t        <- fNashE(y = y, G = Glist, d = dg, talpha, lambdatau = lt, igroup = igr, 
                     nvec = nvec, stau = tau, ngroup = M, n = n, ntau = ntau, type = type, 
                     tol = tol, maxit = maxit)
  
  # Quantile
  Qy       <- NULL
  if (details) {
    Qy     <- fQtauyWithIndex(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                              ngroup = M, n = n, ntau = ntau, type = type)
  } else {
    Qy     <- fQtauy(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                     ngroup = M, n = n, ntau = ntau, type = type)
  }
  
  i        <- NULL
  w        <- NULL
  if (details) {
    i      <- list(pi1 = Qy$pi1, pi2 = Qy$pi2)
    w      <- list(w1 = 1 - Qy$w2, w2 = Qy$w2)
    Qy     <- Qy$Qy
  }
  
  if (ntau == 1){
    Qy       <- c(Qy)
    if (details) {
      i$pi1  <- c(i$pi1)
      i$pi2  <- c(i$pi2)
      w$w1   <- c(w$w1)
      w$w2   <- c(w$w2)
    }
  } else {
    colnames(Qy)      <- paste0("y_q", 1:ntau)
    if (details) {
      colnames(i$pi1) <- paste0("pi1_q", 1:ntau)
      colnames(i$pi2) <- paste0("pi2_q", 1:ntau)
      colnames(w$w1)  <- paste0("w1_q", 1:ntau)
      colnames(w$w2)  <- paste0("w2_q", 1:ntau)
    }
  }
  
  # Output
  list("y"         = y,
       "Qy"        = Qy, 
       "epsilon"   = eps,
       "index"     = i,
       "weight"    = w,
       "iteration" = t)
} 

#' @title Computing Instruments for Linear Models with Quantile Peer Effects
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be specified as, for example, \code{y ~ x1 + x2}, 
#' where `x1` and `x2` are variables for which the quantiles will be computed and `y` is the dependent variable. The quantiles of `x1` and `x2` are computed by ranking observations according to the values of `y` (see details). 
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param tau The vector of quantile levels.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see the \link[stats]{quantile} function).
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `quantvars` is called.
#' @param power The maximum network distance of friends to consider in computing instruments.
#' @param details A logical indicating whether to output the indices and weights of the two peers used to calculate the quantile as a weighted average.
#' @description
#' `qpeer.instruments` computes quantile peer variables. 
#' @details
#' The sample quantile is computed as a weighted average of two peer outcomes (see Hyndman and Fan, 1996). Specifically:
#'  \deqn{q_{\tau,i}(x_{-i}) = (1 - \omega_i)x_{i,(\pi_i)} + \omega_ix_{i,(\pi_i+1)},}
#' where \eqn{x_{i,(1)}, x_{i,(2)}, x_{i,(3)}, \ldots} are the order statistics of the outcome within \eqn{i}'s peers, and \eqn{q_{\tau,i}(x_{-i})} represents the sample \eqn{\tau}-quantile 
#' of the outcome within \eqn{i}'s peer group. To compute the instruments, the ranks \eqn{\pi_i} and the weights \eqn{\omega_i} for the variables in `X` are determined based on `y`.
#' The network matrices in `Glist` can be weighted or unweighted. If weighted, the sample weighted quantile is computed, where the outcome for friend \eqn{j} of \eqn{i} is weighted by \eqn{g_{ij}}, the \eqn{(i, j)} entry of the network matrix.
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @seealso \code{\link{qpeer.lim}}, \code{\link{qpeer.sim}}
#' @return A matrix including quantile peer variables
#' @return A list containing:
#'     \item{qy}{Quantiles of peer variable y.}
#'     \item{instruments}{Matrix of instruments.}
#'     \item{index}{The indices of the two peers whose weighted average gives the quantile.}
#'     \item{weight}{The weights of the two peers whose weighted average gives the quantile.}
#' @examples 
#' ngr  <- 50
#' nvec <- rep(30, ngr)
#' n    <- sum(nvec)
#' G    <- lapply(1:ngr, function(z){
#'   Gz <- matrix(rbinom(sum(nvec[z]*(nvec[z] - 1)), 1, 0.3), nvec[z])
#'   diag(Gz) <- 0
#'   Gz
#' }) 
#' tau  <- seq(0, 1, 0.25)
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' l    <- c(0.2, 0.1, 0.05, 0.1, 0.2)
#' b    <- c(2, -0.5, 1)
#' y    <- qpeer.sim(formula = ~X, Glist = G, tau = tau, lambda = l, beta = b)$y
#' Inst <- qpeer.instruments(formula = y ~ X, Glist = G, tau = tau, power = 2)
#' qy   <- Inst$qy
#' summary(qy)
#' Inst <- Inst$instruments
#' summary(Inst)
#' @export
qpeer.instruments <- function(formula, Glist, tau, type = 7, data, power = 1, details = FALSE){
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  stopifnot(power >= 1)
  # Network
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  dg       <- unlist(lapply(Glist, rowSums))
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  # Data
  formula    <- as.formula(formula)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = TRUE, simulations = FALSE) #Intercept is not necessary
  y          <- f.t.data$y
  X          <- f.t.data$X
  Kx         <- ncol(X)
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  
  # quantiles
  ntau       <- length(tau)
  qy       <- fQtauyIndex(y = y, G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                      ngroup = M, n = n, ntau = ntau, type = type)
  i        <- list(pi1 = qy$pi1, pi2 = qy$pi2)
  w        <- list(w1 = 1 - qy$w2, w2 = qy$w2)
  W        <- lapply(1:ntau, function(s) fIndexMat(pi1 = i$pi1[,s], pi2 = i$pi2[,s], w2 = w$w2[,s], n = n))
  qy       <- do.call(cbind, lapply(W, function(s) fProdWVI(W = s, V = as.matrix(y))))
  ins      <- do.call(cbind, lapply(W, function(s) fProdWVI(W = s, V = as.matrix(X), power = power)))
  
  suffix1  <- NULL
  suffix2  <- NULL
  suffiy   <- NULL
  if (details) {
    i$pi1  <- i$pi1 + 1
    i$pi2  <- i$pi2 + 1
    if (ntau == 1) {
      i$pi1      <- c(i$pi1)
      i$pi2      <- c(i$pi2)
      w$w1       <- c(w$w1)
      w$w2       <- c(w$w2)
      suffix1    <- "_q"
      suffix2    <- rep(paste0("_p", 1:power), each = Kx)
      suffiy     <- "_q"
    } else {
      colnames(i$pi1) <- paste0("pi1_q", 1:ntau)
      colnames(i$pi2) <- paste0("pi2_q", 1:ntau)
      colnames(w$w1)  <- paste0("w1_q", 1:ntau)
      colnames(w$w2)  <- paste0("w2_q", 1:ntau)
      suffix1         <- rep(paste0("_q", 1:ntau), each = Kx*power)
      suffix2         <- rep(rep(paste0("_p", 1:power), each = Kx), ntau)
      suffiy          <- paste0("_q", 1:ntau)
    }
  } else {
    i         <- NULL
    w         <- NULL
    if (ntau == 1) {
      suffix1 <- "_q"
      suffix2 <- rep(paste0("_p", 1:power), each = Kx)
      suffiy  <- "_q"
    } else {
      suffix1 <- rep(paste0("_q", 1:ntau), each = Kx*power)
      suffix2 <- rep(rep(paste0("_p", 1:power), each = Kx), ntau)
      suffiy  <- paste0("_q", 1:ntau)
    }
  }
  colnames(qy)  <- paste0(yname, suffiy)
  colnames(ins) <- paste0(xname, suffix1, suffix2)
  if (ncol(ins) == 1) {
    ins         <- c(ins)
  } 
  if (ncol(qy) == 1) {
    qy          <- c(qy)
  }

  
  list("qy"          = qy,
       "instruments" = ins, 
       "index"       = i,
       "weight"      = w)
  
}

# qpeer.lm <- function(formula, Glist, tau, type = 7, data, power = 1, stuctural = FALSE){
#   
# }