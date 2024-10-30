#' @title Computing Quantile Peer Variables
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are variables for which the quantiles will be computed (see details). 
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' @param tau The vector of quantile levels.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms used to compute peer quantiles (see the \link[stats]{quantile} function).
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `quantvars` is called.
#' @param details A logical indicating whether to output the indices and weights of the two peers used to calculate the quantile as a weighted average.
#' @description
#' `quantvars` computes quantile peer variables. 
#' @details
#' The sample quantile is computed as a weighted average of two peer outcomes (see Hyndman and Fan, 1996). Specifically:
#'  \deqn{q_{\tau,i}(x_{-i}) = (1 - \omega_i)x_{i,(\pi_i)} + \omega_ix_{i,(\pi_i+1)},}
#' where \eqn{x_{i,(1)}, x_{i,(2)}, x_{i,(3)}, \ldots} are the order statistics of the outcome within \eqn{i}'s peers, and \eqn{q_{\tau,i}(x_{-i})} represents the sample \eqn{\tau}-quantile 
#' of the outcome within \eqn{i}'s peer group. The network matrices in `Glist` can be weighted or unweighted. If weighted, the sample weighted quantile is computed, where the outcome for friend \eqn{j} of \eqn{i} is weighted by \eqn{g_{ij}}, the \eqn{(i, j)} entry of the network matrix.
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @seealso \code{\link{qpeer.instruments}}, \code{\link{qpeer.lim}}.
#' @return A matrix including quantile peer variables
#' @return A list containing:
#'     \item{qx}{Quantiles of peer variable X.}
#'     \item{index}{The indices of the two peers whose weighted average gives the quantile.}
#'     \item{weight}{The weights of the two peers whose weighted average gives the quantile.}
#' @examples 
#' ngr  <- 50
#' nvec <- rep(30, ngr)
#' n    <- sum(nvec)
#' G    <- lapply(1:ngr, function(z){
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z])
#'   diag(Gz) <- 0
#'   Gz
#' })
#' X    <- cbind(rnorm(n), rpois(n, 2))
#' out1 <- quantvars(formula = ~ X, Glist = G, tau = seq(0, 1, 0.1), type = 2)
#' out2 <- quantvars(formula = ~ X, Glist = G, tau = seq(0, 1, 0.1), type = 7)
#' @export
quantvars  <- function(formula, Glist, tau, type = 7, data, details = FALSE){
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
  xname      <- NULL
  formula    <- as.formula(formula)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = TRUE,
                                simulations = (length(formula) == 2)) #Intercept is not necessary
  y          <- f.t.data$y
  X          <- f.t.data$X
  Kx         <- ncol(X)
  xname      <- f.t.data$xname
  
  # quantiles
  ntau       <- length(tau)
  qy         <- NULL
  qx         <- list()
  i          <- list()
  w          <- list()
  if (details) {
    qx       <- lapply(1:Kx, function(k){
      tp     <- fQtauyWithIndex(y = X[,k], G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                                ngroup = M, n = n, ntau = ntau, type = type)
      
      tpi    <- list(pi1 = tp$pi1 + 1, pi2 = tp$pi2 + 1)
      tpw    <- list(w1 = 1 - tp$w2, w2 = tp$w2)
      tp     <- tp$Qy
      if(ntau == 1){
        colnames(tp) <- paste(xname[k], "q", sep = "_")
        tpi$pi1      <- c(tpi$pi1)
        tpi$pi2      <- c(tpi$pi2)
        tpw$w1       <- c(tpw$w1)
        tpw$w2       <- c(tpw$w2)
      } else {
        colnames(tp)      <- paste(xname[k], paste0("q", 1:ntau), sep = "_")
        colnames(tpi$pi1) <- paste0("pi1_q", 1:ntau)
        colnames(tpi$pi2) <- paste0("pi2_q", 1:ntau)
        colnames(tpw$w1)  <- paste0("w1_q", 1:ntau)
        colnames(tpw$w2)  <- paste0("w2_q", 1:ntau)
      }
      list(q = tp, i = tpi, w = tpw)
    })
    i      <- lapply(1:Kx, function(k) qx[[k]]$i); names(i) <- xname
    w      <- lapply(1:Kx, function(k) qx[[k]]$w); names(w) <- xname
    qx     <- do.call(cbind, lapply(1:Kx, function(k) qx[[k]]$q))
  } else {
    qx  <- do.call(cbind, lapply(1:Kx, function(k){
      tp   <- fQtauy(y = X[,k], G = Glist, d = dg, igroup = igr, nvec = nvec, stau = tau, 
                     ngroup = M, n = n, ntau = ntau, type = type)
      if(ntau == 1){
        colnames(tp) <- paste(xname[k], "q", sep = "_")
      } else {
        colnames(tp) <- paste(xname[k], paste0("q", 1:ntau), sep = "_")
      }
      tp
    }))
  }
  
  if (ncol(qx) == 1) {
    qx  <- c(qx)
  } 
  
  list("qx"     = qx, 
       "index"  = i,
       "weight" = w)
}

