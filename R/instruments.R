#' @title Computing Instruments for Linear Models with Quantile Peer Effects
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be specified as, for example, \code{~ x1 + x2} or \code{y ~ x1 + x2}, 
#' where `x1` and `x2` are variables for which the quantiles will be computed and `y` is the dependent variable. If `y` is specified, then the quantiles of `x1` and `x2` are computed by ranking observations according to the values of `y` (see details). 
#' @param A The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `A` must be a list of subnets, with the `g`-th element being an \eqn{n_g \times n_g} adjacency matrix, where \eqn{n_g} is the number of nodes in the `g`-th subnet.
#' @param tau The vector of quantile levels.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see the \link[stats]{quantile} function).
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `qpeer.instruments` is called.
#' @param max.distance The maximum network distance of friends to consider in computing instruments.
#' @param checkrank A logical value indicating whether the instrument matrix should be checked for full rank. If the matrix is not of full rank, unimportant columns will be removed to obtain a full-rank matrix.
#' @param tol A tolerance value used in the QR factorization to identify columns that ensure a full-rank matrix (see the \link[base]{qr} function).
#' @param nthreads A strictly positive integer indicating the number of threads to use when computing the quantiles of peer variables.
#' @description
#' `qpeer.instruments` computes quantile peer variables. 
#' @details
#' The sample quantile is computed as a weighted average of two peer outcomes (see Hyndman and Fan, 1996). Specifically:
#'  \deqn{q_{\tau,i}(x_{-i}) = (1 - \omega_i)x_{i,(\pi_i)} + \omega_ix_{i,(\pi_i+1)},}
#' where \eqn{x_{i,(1)}, x_{i,(2)}, x_{i,(3)}, \ldots} are the order statistics of the outcome within \eqn{i}'s peers, and \eqn{q_{\tau,i}(x_{-i})} represents the sample \eqn{\tau}-quantile 
#' of the outcome within \eqn{i}'s peer group. If `y` is specified, then the ranks \eqn{\pi_i} and the weights \eqn{\omega_i} for the variables in `X` are determined based on `y`.
#' The network matrices in `A` can be weighted or unweighted. If weighted, the sample weighted quantile is computed, where the outcome for friend \eqn{j} of \eqn{i} is weighted by \eqn{g_{ij}}, the \eqn{(i, j)} entry of the network matrix.
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @seealso \code{\link{qpeer}}, \code{\link{qpeer.sim}}, \code{\link{linpeer}}
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
#' A    <- lapply(1:ngr, function(g){
#'   Ag <- matrix(rbinom(sum(nvec[g]*(nvec[g] - 1)), 1, 0.3), nvec[g])
#'   diag(Ag) <- 0
#'   Ag
#' }) 
#' tau  <- seq(0, 1, 0.25)
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' lam  <- c(0.2, 0.1, 0.05, 0.1, 0.2)
#' gam  <- c(2, -0.5, 1)
#' y    <- qpeer.sim(formula = ~X, A = A, tau = tau, lambda = lam, gamma = gam)$y
#' Inst <- qpeer.instruments(formula = ~ X, A = A, tau = tau, max.distance = 2)$instruments
#' summary(Inst)
#' @export
qpeer.instruments <- function(formula, A, tau, type = 7, data, max.distance = 1, 
                              checkrank = FALSE, nthreads = 1, tol = 1e-10){
  stopifnot(all((tau >= 0) & (tau <= 1)))
  stopifnot(type %in% 1:9)
  stopifnot(max.distance >= 1)
  nthreads <- fnthreads(nthreads = nthreads)
  
  # Network
  if (!is.list(A)) {
    A  <- list(A)
  }
  d        <- unlist(lapply(A, rowSums))
  G        <- length(A)
  nvec     <- unlist(lapply(A, nrow))
  n        <- sum(nvec)
  igr      <- c(0, cumsum(nvec))
  group    <- rep(0:(G - 1), nvec)
  groupidx <- unlist(lapply(1:G, \(g) 0:(nvec[g] - 1)))
  
  # Data
  formula    <- as.formula(formula)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = TRUE,
                                simulations = (length(formula) == 2)) #Intercept is not necessary
  y          <- f.t.data$y
  X          <- f.t.data$X
  Kx         <- ncol(X)
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  
  # quantiles
  ntau       <- length(tau)
  qy         <- NULL
  ins        <- NULL
  if (!is.null(yname)) {
    tp       <- fQtauyWithIndex(y = y, A = A, d = d, igroup = igr, group = group,
                                groupidx = groupidx, nvec = nvec, stau = tau, ngroup = G, 
                                n = n, ntau = ntau, type = type, nthreads = nthreads)
    W        <- lapply(1:ntau, function(s) fIndexMat(pi1 = tp$pi1[,s], pi2 = tp$pi2[,s], w1 = tp$w1[,s], w2 = tp$w2[,s], n = n))
    qy       <- tp$qy
    ins      <- do.call(cbind, lapply(W, function(s) fProdWVI(W = s, V = as.matrix(X), distance = max.distance)))
  } else {
    W        <- lapply(1:Kx, function(k){
      tp     <- fQtauyIndex(y = X[,k], A = A, d = d, igroup = igr, group = group,
                            groupidx = groupidx, nvec = nvec, stau = tau, 
                            ngroup = G, n = n, ntau = ntau, type = type, nthreads = nthreads)
      lapply(1:ntau, function(s) fIndexMat(pi1 = tp$pi1[,s], pi2 = tp$pi2[,s], w1 = tp$w1[,s], w2 = tp$w2[,s], n = n))})
    ins      <- do.call(cbind, lapply(1:Kx, function(s1) do.call(cbind, lapply(1:ntau, function(s2) fProdWVI(W = W[[s1]][[s2]], V = as.matrix(X[,s1]), distance = max.distance)))))
  }
 
  # Column name
  suffiins  <- NULL
  suffiy    <- NULL
  if (ntau == 1) {
    suffiins <- rep(paste0("_q_p", 1:max.distance), Kx)
    suffiy   <- "_q"
  } else {
    suffiins <- rep(paste0(rep(paste0("_q", 1:ntau), each = max.distance), rep(paste0("_p", 1:max.distance), ntau)), Kx)
    suffiy   <- paste0("_q", 1:ntau)
  }
  if (!is.null(yname)){
    colnames(qy)  <- paste0(yname, suffiy)
    if (ncol(qy) == 1) {
      qy  <- c(qy)
    } 
  }
  colnames(ins)   <- paste0(rep(xname, each = ntau*max.distance), suffiins)
  
  # Checking rank
  if (checkrank) {
    ins <- ins[, fcheckrank(X = ins, tol = tol) + 1, drop = FALSE]
  }
  if (ncol(ins) == 1) {
    ins <- c(ins)
  } 
  
  list("qy"          = qy,
       "instruments" = ins)
}


#' @rdname qpeer.instruments
#' @export
qpeer.instrument <- function(formula, A, tau, type = 7, data, max.distance = 1, 
                             checkrank = FALSE, nthreads = 1, tol = 1e-10){
  qpeer.instruments(formula, A, tau, type, data, max.distance, checkrank, 
                    nthreads, tol)
}

#' @rdname qpeer.instruments
#' @export
qpeer.inst <- function(formula, A, tau, type = 7, data, max.distance = 1, 
                       checkrank = FALSE, nthreads = 1, tol = 1e-10){
  qpeer.instruments(formula, A, tau, type, data, max.distance, checkrank, 
                    nthreads, tol)
}

#' @rdname qpeer.instruments
#' @export
qpeer.insts <- function(formula, A, tau, type = 7, data, max.distance = 1, 
                        checkrank = FALSE, nthreads = 1, tol = 1e-10){
  qpeer.instruments(formula, A, tau, type, data, max.distance, checkrank, 
                    nthreads, tol)
}


#' @title Simulating Optimal (or "Good") Instruments for Quantile Peer Effect Models
#'
#' @description
#' `qpeer.optimal.inst` computes optimal (or "good") instruments for quantile peer effect models using a bootstrap approach.
#'
#' @param boot A strictly positive integer indicating the number of bootstrap replications.
#' @param seed The random number generator (RNG) state used for random number generation. 
#' This can also be set using the \code{\link{set.seed}} function.
#' @param A The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), 
#' `A` must be a list of subnets, with the `g`-th element being an \eqn{n_g \times n_g} adjacency matrix, 
#' where \eqn{n_g} is the number of nodes in the `g`-th subnet.
#' @param maxit The maximum number of iterations for the fixed-point iteration method.
#' @param tol The tolerance value used in the fixed-point iteration method to compute the outcome `y`. 
#' The process stops if the \eqn{\ell_1}-distance between two consecutive values of `y` is less than `tol`.
#' @param tau A numeric vector specifying the quantile levels for the instrument matrix. 
#' The default value is the vector of quantile levels used to estimate the model.
#' @param model An object of class \code{\link{qpeer}} that contains an initial estimation of the model 
#' for which the "optimal" instruments will be computed.
#' @param nthreads A strictly positive integer indicating the number of threads to use when bootstrapping.
#'
#' @return A matrix of instruments, where the k-th column corresponds to the instrument for the k-th endogenous variable.
#' @examples
#' \donttest{
#' set.seed(123)
#' ngr  <- 30  # Number of subnets
#' nvec <- rep(30, ngr)  # Size of subnets
#' n    <- sum(nvec)
#' 
#' ### Simulating Data
#' ## Network matrix
#' A <- lapply(1:ngr, function(g) {
#'   Ag <- matrix(rbinom(nvec[g]^2, 1, 0.3), nvec[g], nvec[g])
#'   diag(Ag) <- 0
#'   # Adding isolated nodes (important for the structural model)
#'   niso <- sample(0:nvec[g], 1, prob = (nvec[g] + 1):1 / sum((nvec[g] + 1):1))
#'   if (niso > 0) {
#'     Ag[sample(1:nvec[g], niso), ] <- 0
#'   }
#'   Ag
#' })
#' 
#' tau <- seq(0, 1, 1/3)
#' X   <- cbind(rnorm(n), rpois(n, 2))
#' lam <- c(0.2, 0.15, 0.1, 0.2)
#' gam <- c(2, -0.5, 1)
#' eps <- rnorm(n, 0, 0.4)
#' 
#' ## Generating `y`
#' y <- qpeer.sim(formula = ~ X, A = A, tau = tau, lambda = lam,
#'                gamma = gam, epsilon = eps)$y
#' 
#' ### Estimation
#' ## Computing instruments
#' Z <- qpeer.inst(formula = ~ X, A = A, tau = seq(0, 1, 0.1),
#'                 max.distance = 2, checkrank = TRUE)
#' Z <- Z$instruments
#' 
#' ## Reduced-form model
#' rest <- qpeer(formula = y ~ X, excluded.instruments = ~ Z, A = A, tau = tau)
#' summary(rest, diagnostic = TRUE)  
#' 
#' # Estimation using the optimal instrument
#' Ired <- qpeer.optimal.insts(rest, A = A, boot = 100, nthreads = 2)
#' summary(qpeer(formula = y ~ X, excluded.instruments = ~ Ired, A = A, tau = tau), 
#'         diagnostic = TRUE)  
#' 
#' ## Structural model
#' sest <- qpeer(formula = y ~ X, excluded.instruments = ~ Z, A = A, tau = tau,
#'               structural = TRUE)
#' summary(sest, diagnostic = TRUE)
#' 
#' # Estimation using the optimal instrument
#' Istr <- qpeer.optimal.insts(sest, A = A, boot = 100, nthreads = 2)
#' summary(qpeer(formula = y ~ X, excluded.instruments = ~ Istr, A = A, tau = tau), 
#'         diagnostic = TRUE)
#' }
#' @importFrom stats runif
#' @export
qpeer.optimal.instruments <- function(model, 
                                      A,
                                      tau,
                                      boot     = 100L, 
                                      nthreads = 1L, 
                                      seed, 
                                      tol      = 1e-10,
                                      maxit    = 500) {
  stopifnot(class(model) == "qpeer")
  nthreads <- fnthreads(nthreads = nthreads)
  if (missing(tau)) {
    tau    <- model$model.info$tau
  }
  stopifnot(all((tau >= 0) & (tau <= 1)))
  if (missing(seed)) {
    seed   <- runif(1, 0, 1e6)
  }
  FE       <- ifelse(model$model.info$fixed.effects == "no", 0, 1)
  if (!is.list(A)) {
    A  <- list(A)
  }
  d        <- unlist(lapply(A, rowSums))
  if (length(d) != model$model.info$n) {
    stop("`A` seems to be different from the one used for the estimation.")
  }
  nvec     <- model$model.info$nvec
  ngroup   <- model$model.info$ngroup
  igroup   <- c(0, cumsum(nvec))
  group    <- rep(0:(ngroup - 1), nvec)
  groupidx <- unlist(lapply(1:ngroup, \(g) 0:(nvec[g] - 1)))
  lIs      <- lapply(model$data$isolated, \(s) s - 1)
  lnIs     <- lapply(model$data$non.isolated, \(s) s - 1)
  Is       <- unlist(lIs)
  nIs      <- unlist(lnIs)
  
  out      <- simInstrqpeer(y = model$data$y, qy = as.matrix(model$data$qy), X = model$data$X, A = A, 
                            d = d, igroup = igroup, group = group, groupidx = groupidx, estimate = model$gmm$Estimate, 
                            nIs = nIs, nvec = nvec, stau = model$model.info$tau, stauInst = tau,
                            boot = boot, fixedeffects = FE, structural =  model$model.info$structural, nthreads = nthreads, 
                            seed = seed, type = model$model.info$type, tol = tol, maxit = maxit)
  colnames(out) <- paste0("Instrument", 1:length(tau))
  out
}


#' @rdname qpeer.optimal.instruments
#' @export
qpeer.optimal.instrument <- function(model, 
                                     A,
                                     tau,
                                     boot     = 100L, 
                                     nthreads = 1L, 
                                     seed, 
                                     tol      = 1e-10,
                                     maxit    = 500) {
  qpeer.optimal.instruments(model, A, tau, boot, nthreads, seed, tol, maxit)
}

#' @rdname qpeer.optimal.instruments
#' @export
qpeer.optimal.insts <- function(model, 
                                A,
                                tau, 
                                boot     = 100L, 
                                nthreads = 1L, 
                                seed, 
                                tol      = 1e-10,
                                maxit    = 500) {
  qpeer.optimal.instruments(model, A, tau, boot, nthreads, seed, tol, maxit)
}

#' @rdname qpeer.optimal.instruments
#' @export
qpeer.optimal.inst <- function(model, 
                               A,
                               tau, 
                               boot     = 100L, 
                               nthreads = 1L, 
                               seed, 
                               tol      = 1e-10,
                               maxit    = 500) {
  qpeer.optimal.instruments(model, A, tau, boot, nthreads, seed, tol, maxit)
}

