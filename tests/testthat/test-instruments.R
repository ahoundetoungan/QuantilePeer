test_that("qpeer.instruments returns expected output", {
  ngr  <- 5
  nvec <- rep(30, ngr)
  n    <- sum(nvec)
  ## Network matrix
  G <- lapply(1:ngr, function(z) {
    Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z], nvec[z])
    diag(Gz) <- 0
    # Adding isolated nodes (important for the structural model)
    niso <- sample(0:nvec[z], 1, prob = (nvec[z] + 1):1 / sum((nvec[z] + 1):1))
    if (niso > 0) {
      Gz[sample(1:nvec[z], niso), ] <- 0
    }
    Gz
  })
  
  iso    <- unlist(lapply(G, function(x) rowSums(x) == 0))
  tau    <- seq(0, 1, 0.25)
  X      <- cbind(rnorm(n), rpois(n, 2))
  l      <- c(0.2, 0.1, 0.05, 0.1, 0.2)
  b      <- c(2, -0.5, 1)
  
  ftest  <- function(type) {
  y      <- qpeer.sim(formula = ~X, Glist = G, tau = tau, lambda = l, beta = b, type = type)$y
  qy     <- qpeer.instruments(formula = y ~ X, Glist = G, tau = tau, max.distance = 1, type = type)$qy
  Inst1  <- qpeer.instruments(formula = ~ X, Glist = G, tau = tau, max.distance = 1, type = type)$instruments
  Inst2  <- qpeer.instruments(formula = y ~ X, Glist = G, tau = tau, max.distance = 1, type = type)
  
  expect_equal(rowSums(abs(qy[iso,])), rep(0, sum(iso))) # Testing quantiles of isolated peers
  
  # Testing for consistency with qpeer.inst qpeer.insts qpeer.instruments
  expect_equal(c(round(Inst1, 10)),
               c(round(qpeer.instruments(formula = ~ X, Glist = G, tau = tau, 
                                         max.distance = 1, type = type)$instruments, 10)))
  expect_equal(c(round(Inst1, 10)),
               c(round(qpeer.inst(formula = ~ X, Glist = G, tau = tau, 
                                  max.distance = 1, type = type)$instruments, 10)))
  expect_equal(c(round(Inst1, 10)),
               c(round(qpeer.insts(formula = ~ X, Glist = G, tau = tau, 
                                   max.distance = 1, type = type)$instruments, 10)))
  
  #check quantile for non isolates
  ncum   <- c(0, cumsum(nvec))
  n      <- sum(nvec) 
  
  # friend indices
  friend <- do.call(c, lapply(1:ngr, function(m) lapply(1:nvec[m], function(i) which(G[[m]][i,] > 0) + ncum[m])))
  
  # Testing quantiles of y
  qytest <- t(sapply(1:n, function(i) quantile(y[friend[[i]]], probs = tau, type = type)))
  expect_equal(c(round(qy[!iso,], 10)), c(round(qytest[!iso,], 10)))
  
  # Testing quantiles of X
  qXtest  <- cbind(t(sapply(1:n, function(i) quantile(X[friend[[i]], 1], probs = tau, type = type))), 
                   t(sapply(1:n, function(i) quantile(X[friend[[i]], 2], probs = tau, type = type))))
  expect_equal(c(round(qXtest[!iso,], 10)), c(round(Inst1[!iso,], 10)))
  
  # Testing quantile of X based on ranks of y
  expect_true(all(Inst2$qy == qy))
  NULL
  }
  lapply(1:9, ftest) # Test for all types of quantiles
})

# Need to test for max.distance