#' @title Demeaning Variables
#' 
#' @param X A matrix or vector to demean.
#' 
#' @param A The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `A` must be a list of subnets, with the `g`-th element being an \eqn{n_g \times n_g} adjacency matrix, where \eqn{n_g} is the number of nodes in the `g`-th subnet.
#' 
#' @param drop A logical vector of the same length as the sample, indicating whether an observation should be dropped. 
#' This can be used, for example, to remove false isolates or to estimate the model only on non-isolated agents.
#' These observations cannot be directly removed from the network by the user, as they may still be connected to other agents.
#' 
#' @param separate A logical value specifying whether variables should be demeaned separately for isolated and non-isolated individuals.
#' This is similar to setting `fixed.effects = "separate"` in \code{\link{qpeer}}.
#' 
#' @description
#' `demean` demeans variables by subtracting the within-subnetwork average. In each subnetwork, this transformation can be performed separately
#' for isolated and non-isolated nodes.
#' 
#' @return A matrix or vector with the same dimensions as \code{X}, containing the demeaned values.
#' @export
demean <- function(X, A, separate = FALSE, drop = NULL) {
  X    <- as.matrix(X)
  cn   <- colnames(X)
  # Network
  if (!is.list(A)) {
    A  <- list(A)
  }
  d        <- fnetwork(A = A)
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
  d        <- d$d
  
  # drop
  if (!is.null(drop)) {
    d        <- fdrop(drop = drop, ld = ld, nvec = nvec, G = G, lIs = lIs, lnIs = lnIs, 
                      y = rep(0, n), X = X, qy = matrix(0, n, 1), ins = matrix(0, n, 1))
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
    qy       <- d$qy
    ins      <- d$ins
    d        <- d$d
  }
  if (separate) {
    X        <- Demean_separate(X, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = G, n = n)
  } else {
    X        <- Demean(X, igroup = igr, ngroup = G)
  }
  
  if (ncol(X) == 1) {
    X        <- c(X)
  }
  colnames(X)<- cn
  X
}