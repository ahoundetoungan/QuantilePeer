#' @title Demeaning Variables
#' 
#' @param X A matrix or vector to demean.
#' 
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
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
#' @export
demean <- function(X, Glist, separate = FALSE, drop = NULL) {
  X    <- as.matrix(X)
  cn   <- colnames(X)
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
  
  # drop
  if (!is.null(drop)) {
    dg       <- fdrop(drop = drop, ldg = ldg, nvec = nvec, M = M, lIs = lIs, lnIs = lnIs, 
                      y = rep(0, n), X = X, qy = matrix(0, n, 1), ins = matrix(0, n, 1))
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
  if (separate) {
    X        <- Demean_separate(X, igroup = igr, LIs = lIs, LnIs = lnIs, ngroup = M, n = n)
  } else {
    X        <- Demean(X, igroup = igr, ngroup = M)
  }
  
  if (ncol(X) == 1) {
    X        <- c(X)
  }
  colnames(X)<- cn
  X
}