# Structural
fstruct <- function(y, X, qy, ins, idX1, idX2, Kx1, Kx2, igr, nIs, Is, lnIs, lIs, M, Kins, Kx, 
                    ntau, Kest1, Kest2,  n, HACnum, iv, estimator, compute.cov, estname) {
  Kest        <- Kest1 + Kest2
  GMMe        <- NULL
  if (estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
    GMMe      <- fgmm_struc(y = y, X = X, qy = qy, ins = ins, W1 = diag(Kx1), W2 = diag(Kins + 1), idX1 = idX1, 
                            idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igroup = igr, nIs = nIs, Is = Is, ngroup = M, 
                            Kins = Kins, Kx = Kx, ntau = ntau, Kest1 = Kest1, Kest2 = Kest2, n = n, HAC = HACnum, 
                            iv = iv)
    
    if (estimator == "GMM.optimal" & HACnum != 0) {
      GMMe    <- fgmm_struc(y = y, X = X, qy = qy, ins = ins, W1 = solve(GMMe$VF1), W2 = solve(GMMe$VF2), 
                            idX1 = idX1, idX2 = idX2, Kx1 = Kx1, Kx2 = Kx2, igroup = igr, nIs = nIs, Is = Is, 
                            ngroup = M, Kins = Kins, Kx = Kx, ntau = ntau, Kest1 = Kest1, Kest2 = Kest2, n = n, 
                            HAC = HACnum, iv = FALSE)
    }
  } else if (estimator == "JIVE") {
    if (HACnum == 2) {
      nvecnIs <- sapply(lnIs, length)
      nvecIs  <- sapply(lIs, length)
      hasnIs  <- as.integer(nvecnIs > 0)
      hasIs   <- as.integer(nvecIs > 0)
      GMMe    <- fJIVE_strucClu(y = y, X = X, qy = qy, ins = ins, idX1 = idX1, idX2 = idX2, nIs = nIs, Is = Is, 
                                LnIs = lnIs, LIs = lIs, igroup = igr, nvecnIs = nvecnIs, hasnIs = hasnIs, 
                                hasIs = hasIs, ngroup = M, Kx1 = Kx1, Kx2 = Kx2, Kins = Kins, ntau = ntau, n = n, 
                                COV = compute.cov)
    } else {
      GMMe    <- fJIVE_strucInd(y = y, X = X, qy = qy, ins = ins, idX1 = idX1, idX2 = idX2, nIs = nIs, 
                                Is = Is, Kx1 = Kx1, Kx2 = Kx2, Kins = Kins, ntau = ntau, n = n, Kest1 = Kest1, 
                                HAC = HACnum, COV = compute.cov)
    }
  } else {
    if (HACnum == 2) {
      nvecnIs <- sapply(lnIs, length)
      nvecIs  <- sapply(lIs, length)
      hasnIs  <- as.integer(nvecnIs > 0)
      hasIs   <- as.integer(nvecIs > 0)
      GMMe    <- fJIVE2_strucClu(y = y, X = X, qy = qy, ins = ins, idX1 = idX1, idX2 = idX2, nIs = nIs, Is = Is, 
                                 LnIs = lnIs, LIs = lIs, igroup = igr, hasnIs = hasnIs, 
                                 hasIs = hasIs, ngroup = M, Kx1 = Kx1, Kx2 = Kx2, Kins = Kins, ntau = ntau, n = n, 
                                 COV = compute.cov)
    } else {
      GMMe    <- fJIVE2_strucInd(y = y, X = X, qy = qy, ins = ins, idX1 = idX1, idX2 = idX2, nIs = nIs, 
                                 Is = Is, Kx1 = Kx1, Kx2 = Kx2, Kins = Kins, ntau = ntau, n = n, Kest1 = Kest1, 
                                 HAC = HACnum, COV = compute.cov)
    }
  }
  if (!compute.cov) {
    GMMe$Vpa  <- matrix(NA, 1 + ntau + Kx1 + Kx2, 1 + ntau + Kx1 + Kx2)
  }
  Vpa         <- fStructParam(param = c(GMMe$beta, GMMe$lambda), covp = GMMe$Vpa, idX1 = idX1, idX2 = idX2, 
                              ntau = ntau, Kx = Kx, Kx1 = Kx1, Kx2 = Kx2, COV = compute.cov)
  GMMe$parms  <- c(Vpa$theta)
  GMMe$Vpa    <- Vpa$Vpa
  # print(GMMe$Vpa)
  
  ###
  fv          <- c(GMMe$yhat)
  res         <- y - fv
  res[nIs + 1]<- res[nIs + 1]/(1 - GMMe$parms[1])
  rs          <- sum((fv - mean(fv))^2)/sum((y - mean(y))^2)
  ars         <- 1 - (1 - rs)*(n - 1)/(n - length(c(GMMe$parms)))
  sigmaiso    <- sqrt(GMMe$sigma21);  if(is.na(sigmaiso)) sigmaiso = NULL
  sigmaniso   <- sqrt(GMMe$sigma22);  if(is.na(sigmaniso)) sigmaniso = NULL
  GMMe        <- list(Estimate = c(GMMe$parms), cov = GMMe$Vpa, sigma1 = sigmaiso, sigma2 = sigmaniso, fitted.values = fv, 
                      residuals = res, rsquared = rs, adjusted.rsquared = ars, df.residual = n - Kest,
                      Jtest = c("statistic" = ifelse(is.null(GMMe$Overident), NA, GMMe$Overident), 
                                "df" = ifelse(is.null(GMMe$df), NA, as.integer(GMMe$df))), Wiso = GMMe$W1, Wniso = GMMe$W2)
  names(GMMe$Estimate)  <- colnames(GMMe$cov) <- rownames(GMMe$cov)    <- estname
  GMMe$Jtest["p-value"] <- ifelse(GMMe$Jtest["df"] > 0, 1 - pchisq(GMMe$Jtest["statistic"], GMMe$Jtest["df"]), NA)
  GMMe
}


# Reduced form
freduce <- function(y, V, ins, igr, nvec, M, Kins, Kx, ntau, Kest, n, HACnum, iv, 
                    estimator, compute.cov, estname) {
  GMMe        <- NULL
  if (estimator %in% c("IV", "GMM.optimal", "GMM.identity")) {
    GMMe      <- fgmm_red(y = y, V = V, ins = ins, W = diag(Kins), igroup = igr, ngroup = M, 
                          Kx = Kx, Kins = Kins, ntau = ntau, Kest = Kest, n = n, HAC = HACnum, iv = iv)
    if (estimator == "GMM.optimal" & HACnum != 0) {
      GMMe    <- fgmm_red(y = y, V = V, ins = ins, W = solve(GMMe$VZe), igroup = igr, ngroup = M, Kx = Kx, 
                          Kins = Kins, ntau = ntau, Kest = Kest, n = n, HAC = HACnum, iv = FALSE)
    }
  } else if (estimator == "JIVE") {
    if (HACnum == 2) {
      GMMe    <- fJIVE_redClu(y = y, V = V, ins = ins, igroup = igr, nvec = nvec, ngroup = M, Kx = Kx, 
                              Kins = Kins, ntau = ntau, n = n, Kest = Kest, COV = compute.cov)
    } else {
      GMMe    <- fJIVE_redInd(y = y, V = V, ins = ins, Kx = Kx, Kins = Kins, ntau = ntau, n = n, Kest = Kest, 
                              HAC = HACnum, COV = compute.cov)
    }
  } else {
    if (HACnum == 2) {
      GMMe    <- fJIVE2_redClu(y = y, V = V, ins = ins, igroup = igr, nvec = nvec, ngroup = M, Kx = Kx, 
                               Kins = Kins, ntau = ntau, n = n, Kest = Kest, COV = compute.cov)
    } else {
      GMMe    <- fJIVE2_redInd(y = y, V = V, ins = ins, Kx = Kx, Kins = Kins, ntau = ntau, n = n, Kest = Kest, 
                               HAC = HACnum, COV = compute.cov)
    }
  }
  
  ###
  if (!compute.cov) {
    GMMe$Vpa  <- matrix(NA, ntau + Kx, ntau + Kx)
  }
  fv          <- c(GMMe$yhat)
  res         <- y - fv
  rs          <- sum((fv - mean(fv))^2)/sum((y - mean(y))^2)
  ars         <- 1 - (1 -rs)*(n - 1)/(n - length(c(GMMe$parms)))
  sigma       <- sqrt(GMMe$sigma2);  if(is.na(sigma)) sigma = NULL
  GMMe        <- list(Estimate = c(GMMe$parms), cov = GMMe$Vpa, sigma = sigma, fitted.values = fv, 
                      residuals = res, rsquared = rs, adjusted.rsquared = ars, df.residual = n - Kest,
                      Jtest = c("statistic" = ifelse(is.null(GMMe$Overident), NA, GMMe$Overident), 
                                "df" = ifelse(is.null(GMMe$df), NA, as.integer(GMMe$df))), W = GMMe$W)
  names(GMMe$Estimate)  <- colnames(GMMe$cov) <- rownames(GMMe$cov)    <- estname
  GMMe$Jtest["p-value"] <- ifelse(GMMe$Jtest["df"] > 0, 1 - pchisq(GMMe$Jtest["statistic"], GMMe$Jtest["df"]), NA)
  GMMe
}