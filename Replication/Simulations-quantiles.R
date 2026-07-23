################################################################################
################################################################################
################ Identification of Heterogeneous Peer Effects ##################
########### Eyo I. Herstad, Aristide Houndetoungan, Myungkou Shin ##############
################################################################################
################################################################################

# Last updated: 2026-07-08

# This script reproduces the Monte Carlo simulations (Tables 4.2 and 4.3).
# Estimation results will be exported to a single Excel file (`Simulation.xlsx`)
# located in `OutResPath`. The Excel file contains the following sheets:
#  - Table 4.2 (reproduces Table 4.2),
#  - Table 4.3 (reproduces Table 4.3),
#  - Full results (not included in the paper).

rm(list = ls())
library(QuantilePeer)
library(CDatanet) 
library(parallel)
library(dplyr)
library(openxlsx)

# Where results should be saved
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED"
OutResPath  <- "~/personal/Quantile Peer Effects/Package/QuantilePeer/Results"
OutResPath  <- "~/Dropbox/Academy/1.Papers/Quantile Peer Effects/Package/QuantilePeer/Replication"

# Sample size:
ngroup     <- 50  # Number of subnetworks
nvec       <- rep(50, ngroup)  # Size of each subnetwork
n          <- sum(nvec)  # Total sample size

# tau
tau        <- seq(0, 1, 1/3)
ntau       <- length(tau)

# Contextual 
contextual <- TRUE

# Function: Takes lambda, simulates data, estimates the model, and returns estimates.
# lambda is the vector of peer effect parameters
# fixed.effects indicates whether the model includes fixed effects
# linear indicates whether data should be simulated using tje standar linear model
festim     <- function(lambda, fixed.effects = TRUE, linear = FALSE, 
                       nthreads = 7) {
  
  # We consider 2 explanatory variables and control for contextual variables
  # parameter gamma
  gamma    <- c(1.5, -0.8, -0.5, 1.2)
  if (!contextual) {
    gamma  <- c(1.5, -0.8)
  }
  
  # parameter sigma
  sigma    <- 1
  
  # Network
  N   <- sum(nvec*(nvec - 1))
  P   <- vec.to.mat(plogis(-3 + rnorm(N)), nvec)
  A   <- simnetwork(P)
  summary(unlist(lapply(A, rowSums)))
  
  # Fixed effects
  eff      <- rep(runif(ngroup, 8, 10), nvec)
  # X
  X        <- cbind(rnorm(n, eff, 1), rpois(n, 2)) # Make the effects correlated with X
  # AX AAX
  Anorm    <- norm.network(A) # Row-normalized network
  AX       <- peer.avg(Anorm, X)
  
  # Outcome y
  resim    <- TRUE
  Zces     <- NULL
  y        <- NULL
  while (resim) {
    form   <- NULL
    if (contextual) {
      form <- as.formula(~ -1 + eff + X + AX)
    } else {
      form <- as.formula(~ -1 + eff + X)
    }
    
    if (linear) {
      # If linear lambda is set to the sum of lambda
      y    <- linpeer.sim(formula = form, A = Anorm, lambda = sum(lambda), 
                          gamma = c(ifelse(fixed.effects, 1, 0), gamma), structural = FALSE, 
                          epsilon = rnorm(n, 0, sigma))$y 
    } else {
      y    <- qpeer.sim(formula = form, A = A, tau = tau, lambda = lambda, 
                        gamma = c(ifelse(fixed.effects, 1, 0), gamma), structural = FALSE, 
                        epsilon = rnorm(n, 0, sigma))$y 
    }
    
    # Instrument for CES model: exogenous prediction of y
    if (fixed.effects) {
      
      form <- as.formula(y ~ X + AX + grp)
      grp  <- as.factor(rep(1:ngroup, nvec))
      Zces <- fitted.values(lm(form))
      
    } else {
      
      form <- as.formula(y ~ X + AX)
      Zces <- fitted.values(lm(form))
      
    }
    
    # Resimulate
    ## This is because the CES does not allow negative outcome
    ## But with the large values of eff, this should not happen often
    ## it happen 2 times in 60K simulations.
    resim  <- (any(y < 0) | any(Zces < 0)) 
  }
  
  # Instruments for LiM and Quantile
  # linear model: A^2X
  AAX      <- peer.avg(Anorm, AX) #AX where A is row-normalized
  
  # Quantile model
  form     <- as.formula(~ X + AX)
  Z        <- qpeer.inst(formula = form, A = A, tau = seq(0, 1, 1/9), 
                         max.distance = 2, checkrank = TRUE)$instruments 
  
  # Estimation
  # Linear model
  
  form     <- NULL
  if (contextual) {
    form   <- as.formula(y ~ X + AX)
  } else {
    form   <- as.formula(y ~ X)
  }
  
  inst     <- NULL
  if (contextual) {
    inst   <- as.formula(~ AAX)
  } else {
    inst   <- as.formula(~ AX + AAX)
  }
  
  LIM     <- linpeer(formula = form, 
                     excluded.instruments = inst, 
                     A = Anorm,
                     structural = FALSE,
                     fixed.effects = fixed.effects)
  
  # Quantile-based model using Z1 as instruments
  Quant   <- qpeer(formula = form, 
                   excluded.instruments = ~ Z, 
                   A = Anorm,
                   structural = FALSE,
                   tau = tau,
                   fixed.effects = fixed.effects)
  
  # CES-based model
  Ces     <- cespeer(formula = form, instrument = ~ Zces, A = Anorm, 
                     structural = FALSE, fixed.effects = fixed.effects, radius = 5,
                     grid.rho = seq(-100, 100, 1))
  
  # Quantile model with 2 quantiles
  Quant2  <- qpeer(formula = form, 
                   excluded.instruments = ~ Z, 
                   A = Anorm,
                   structural = FALSE,
                   tau = seq(0, 1, 1),
                   fixed.effects = fixed.effects)
  
  # Quantile model with 3 quantiles
  Quant3  <- qpeer(formula = form, 
                   excluded.instruments = ~ Z, 
                   A = Anorm,
                   structural = FALSE,
                   tau = seq(0, 1, 1/2),
                   fixed.effects = fixed.effects)
  
  # Quantile model with 5 quantiles
  Quant5  <- qpeer(formula = form, 
                   excluded.instruments = ~ Z, 
                   A = Anorm,
                   structural = FALSE,
                   tau = seq(0, 1, 1/4),
                   fixed.effects = fixed.effects)
  
  # Encompassing test
  Etest    <- c(qpeer.test(Quant2, Quant3, which = "encompassing", boot = 1e3, 
                           print = FALSE, nthreads = nthreads)$test["KP Wald rank", "p-value"], # testing where 4 quantiles is better than 3 quantiles
                qpeer.test(Quant3, Quant, which = "encompassing", boot = 1e3, 
                           print = FALSE, nthreads = nthreads)$test["KP Wald rank", "p-value"], # testing where 3 quantiles is better than 4 quantiles
                qpeer.test(Quant, Quant5, which = "encompassing", boot = 1e3, 
                           print = FALSE, nthreads = nthreads)$test["KP Wald rank", "p-value"]) # testing where 5 quantiles is better than 4 quantiles
  
  LIM           <- summary(LIM, diagnostic = TRUE)
  LIM           <- c(LIM$gmm$Estimate, KPstat = LIM$diagnostics[2, 3], Jpvalue = LIM$diagnostics[4, 4],
                     Jpvalue05 = LIM$diagnostics[4, 4] < 0.05, Jpvalue10 = LIM$diagnostics[4, 4] < 0.1)
  names(LIM)    <- paste0("LIM.", names(LIM))
  
  Quant         <- summary(Quant, diagnostic = TRUE)
  Quant         <- c(Quant$gmm$Estimate, KPstat = Quant$diagnostics[ntau + 1, 3], 
                     Jpvalue = Quant$diagnostics[ntau + 3, 4], Jpvalue05 = Quant$diagnostics[ntau + 3, 4] < 0.05,
                     Jpvalue10 = Quant$diagnostics[ntau + 3, 4] < 0.1)
  names(Quant)  <- paste0("Q.", names(Quant))
  
  Ces           <- summary(Ces)$gmm$Estimate
  names(Ces)    <- paste0("CES.", names(Ces))
  
  Etest         <- c(Etest, Etest < 0.05, Etest < 0.1)
  names(Etest)  <- paste0(c("ET.ntau=2.3", "ET.ntau=3.4", "ET.ntau=4.5",
                            "ET05.ntau=2.3", "ET05.ntau=3.4", "ET05.ntau=4.5", 
                            "ET10.ntau=2.3", "ET10.ntau=3.4", "ET10.ntau=4.5"))
  
  Quant2        <- summary(Quant2)
  Quant2        <- c(Quant2$gmm$Estimate, KPstat = Quant2$diagnostics[ntau, 3], 
                     Jpvalue = Quant2$diagnostics[ntau + 2, 4], Jpvalue05 = Quant2$diagnostics[ntau + 2, 4] < 0.05, 
                     Jpvalue10 = Quant2$diagnostics[ntau + 2, 4] < 0.1)
  names(Quant2) <- paste0("Q2.", names(Quant2))
  
  Quant3        <- summary(Quant3)
  Quant3        <- c(Quant3$gmm$Estimate, KPstat = Quant3$diagnostics[ntau + 2, 3], 
                     Jpvalue = Quant3$diagnostics[ntau + 4, 4], Jpvalue05 = Quant3$diagnostics[ntau + 4, 4] < 0.05,
                     Jpvalue10 = Quant3$diagnostics[ntau + 4, 4] < 0.1)
  names(Quant3) <- paste0("Q3.", names(Quant3))
  
  Quant5        <- summary(Quant5)
  Quant5        <- c(Quant5$gmm$Estimate, KPstat = Quant5$diagnostics[ntau + 2, 3], 
                     Jpvalue = Quant5$diagnostics[ntau + 4, 4], Jpvalue05 = Quant5$diagnostics[ntau + 4, 4] < 0.05,
                     Jpvalue10 = Quant5$diagnostics[ntau + 4, 4] < 0.1)
  names(Quant5) <- paste0("Qb.", names(Quant5))
  
  
  c(LIM, Quant, Ces, Etest, Quant2, Quant3, Quant5)
}

# Number of simulations
nsim      <- 3

# Several values will be tested for peer effects with conformity parameter set to 0.2
# Estimation is done in parallel 
nthreads  <- 30
set.seed(2026)  # Set global seed for reproducibility

# Increasing lambda
lambda1    <- c(0, 0.05, 0.2, 0.3)
cat("========== DGP 1\n")
Est1       <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda1, TRUE, FALSE, nthreads)
}))

# Decreasing lambda
lambda2    <- c(0.3, 0.2, 0.05, 0)
cat("========== DGP 2\n")
Est2       <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda2, TRUE, FALSE, nthreads)
}))

# Concave lambda
lambda3    <- c(0, 0.275, 0.275, 0)
cat("========== DGP 3\n")
Est3      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda3, TRUE, FALSE, nthreads)
}))

# Convex lambda
lambda4    <- c(0.275, 0, 0, 0.275)
cat("========== DGP 4\n")
Est4       <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda4, TRUE, FALSE, nthreads)
}))

# Concave lambda
lambda5    <- c(-0.05, 0.35, 0.15, 0.1)
cat("========== DGP 5\n")
Est5       <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda5, TRUE, FALSE, nthreads)
}))

# Constant lambda with data simulated using the standard LIM model
lambda6    <- 0.55
cat("========== DGP 6\n")
Est6       <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda6, TRUE, TRUE, nthreads)
}))

save(Est1, Est2, Est3, Est4, Est5, Est6,
     file = paste0(OutResPath, "/Simulations", ifelse(contextual, "", "NoContext"),  ".Rda"))

# Summary
# Summary function
Sumfunc    <- function(x) {
  c(mean = mean(x), sd = sd(x))
}

load(paste0(OutResPath, "/Simulations", ifelse(contextual, "", "NoContext"),  ".Rda"))
Est        <- as.data.frame(rbind(cbind(apply(Est1, 1, Sumfunc)),
                                  cbind(apply(Est2, 1, Sumfunc)),
                                  cbind(apply(Est3, 1, Sumfunc)),
                                  cbind(apply(Est4, 1, Sumfunc)),
                                  cbind(apply(Est5, 1, Sumfunc)),
                                  cbind(apply(Est6, 1, Sumfunc))))

# Formatting to Excel
# Put standard error into parenthesis
Est <- Est %>% mutate(across(everything(), ~ sprintf("%.3f", .))) %>%  # Ensures exactly 3 decimals
  mutate(across(everything(), ~ ifelse(row_number() %% 2 == 0, paste0("(", trimws(.), ")"), .)))  # Parentheses without spaces 

# Insert blank rows between DGP
tp           <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(Est)))
colnames(tp) <- colnames(Est)
Est          <- Est %>% add_row(tp, .before = 11) %>%
  add_row(tp, .before = 9) %>%
  add_row(tp, .before = 7) %>%
  add_row(tp, .before = 5) %>%
  add_row(tp, .before = 3) %>%
  add_row(tp, .before = 1)

# Export to Excel
# Create the workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "Table 4.3") 
addWorksheet(wb, "Table 4.2") 
addWorksheet(wb, "Full results")

# Table 4.3
T4.3 <- Est %>% select(all_of(c(paste0("Q.y_q", 1:4), "LIM.A:y", "CES.rho", "CES.A:y")))
T4.3[seq(1, 16, 3), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1, collapse = ", "), ")$"),
                            paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2, collapse = ", "), ")$"),
                            paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3, collapse = ", "), ")$"),
                            paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4, collapse = ", "), ")$"),
                            paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5, collapse = ", "), ")$"),
                            paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$"))
T4.3 <- cbind(T4.3[,c("Q.y_q1", "Q.y_q2", "Q.y_q3", "Q.y_q4")], "V1" = NA,
              T4.3[,c("LIM.A:y")], "V2" = NA,
              T4.3[,c("CES.rho", "CES.A:y")])
# write
writeData(wb, "Table 4.3", T4.3, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)
# first row
fr    <- c(paste("lmda_q", 1:4), "", "lmda", "", "rho", "lmda")
for (i in 1:length(fr)) {
  writeData(wb, "Table 4.3", fr[i], startRow = 1, startCol = i)
}
# Merge cells
for (i in seq(1, nrow(T4.3), 3) + 1) {
  mergeCells(wb, "Table 4.3", cols = 1:ncol(T4.3), rows = i)
}
# Formatting
addStyle(wb, "Table 4.3", style = createStyle(halign = "center"), 
         rows = 1:(nrow(T4.3) + 1), cols = 1:ncol(T4.3), gridExpand = TRUE)

# Table 4.2
T4.2 <- Est %>% select(all_of(c("ET05.ntau=2.3", "ET10.ntau=2.3",
                                "ET05.ntau=3.4", "ET10.ntau=3.4", 
                                "ET05.ntau=4.5", "ET10.ntau=4.5" # Encompassing tests
))) %>% 
  filter((row_number() - 1) %% 3 != 2)
T4.2[seq(1, 11, 2), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1, collapse = ", "), ")$"),
                            paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2, collapse = ", "), ")$"),
                            paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3, collapse = ", "), ")$"),
                            paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4, collapse = ", "), ")$"),
                            paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5, collapse = ", "), ")$"),
                            paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$"))
T4.2 <- cbind(T4.2[,c("ET05.ntau=2.3", "ET10.ntau=2.3")], 
              V1 = NA,
              T4.2[,c("ET05.ntau=3.4", "ET10.ntau=3.4")], 
              V2 = NA,
              T4.2[,c("ET05.ntau=4.5", "ET10.ntau=4.5")])
# write
writeData(wb, "Table 4.2", T4.2, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)
# first row
fr    <- c("2vs3", "2vs3", "", "3vs4", "3vs4", "", "4vs5", "4vs5")
for (i in 1:length(fr)) {
  writeData(wb, "Table 4.2", fr[i], startRow = 1, startCol = i)
}
# Merge cells
for (i in seq(1, nrow(T4.2), 2) + 1) {
  mergeCells(wb, "Table 4.2", cols = 1:ncol(T4.2), rows = i)
}
# Formatting
addStyle(wb, "Table 4.2", style = createStyle(halign = "center"), 
         rows = 1:(nrow(T4.2) + 1), cols = 1:ncol(T4.2), gridExpand = TRUE)

# Full results
tp <- Est
tp[seq(1, 16, 3), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1, collapse = ", "), ")$"),
                          paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2, collapse = ", "), ")$"),
                          paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3, collapse = ", "), ")$"),
                          paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4, collapse = ", "), ")$"),
                          paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5, collapse = ", "), ")$"),
                          paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$"))
writeData(wb, "Full results", tp, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, paste0(OutResPath, "/Simulations", ifelse(contextual, "", "NoContext"),  ".xlsx"), 
             overwrite = TRUE)
