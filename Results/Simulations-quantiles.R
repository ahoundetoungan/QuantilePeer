################################################################################
################################################################################
################ Identification of Heterogeneous Peer Effects ##################
########### Eyo I. Herstad, Aristide Houndetoungan, Myungkou Shin ##############
################################################################################
################################################################################

# Last updated: 2026-3-25

# This script reproduces the Monte Carlo simulations (Tables 5.3 and 5.4).
# Estimation results will be exported to a single Excel file (`Simulation.xlsx`)
# located in `OutResPath`. The Excel file contains the following sheets:
#  - Table 5.4 (reproduces Table 5.4),
#  - Table 5.3 (reproduces Table 5.3),
#  - Full results (not included in the paper).

rm(list = ls())
library(QuantilePeer)
library(CDatanet) 
library(parallel)
library(dplyr)
library(openxlsx)

# Where results should be saved
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED"

# I consider 2 explanatory variables and control for contextual variables
# parameter beta
beta       <- c(1.5, -0.8, -0.5, 1.2)
# parameter sigma
sigma      <- 1

# Sample size:
ngroup     <- 50  # Number of subnetworks
nvec       <- rep(50, ngroup)  # Size of each subnetwork
n          <- sum(nvec)  # Total sample size

# Network degree distribution: I use the empirical distribution of AddHealth
# Degree possible values (support):
dvalues    <- 0:10  # Agents can have up to 10 friends
# Degree distribution:
ddvalues   <- c(0.22175143, 0.09047220, 0.10325461, 0.11262459, 0.11128805, 0.10039670,
                0.08578010, 0.07272753, 0.05633362, 0.03411014, 0.01126104)  # Degree distribution

# tau
tau        <- seq(0, 1, 1/3)
ntau       <- length(tau)

# Function: Takes lambda, simulates data, estimates the model, and returns estimates.
# lambda is the vector of peer effect parameters
# fixed.effects indicates whether the model includes fixed effects
# linear indicates whether data should be simulated using tje standar linear model
festim     <- function(lambda, fixed.effects = TRUE, linear = FALSE, nthreads = 7) {
  # Simulated degree for each agent in each subnetwork
  degree   <- lapply(1:ngroup, function(s) sample(dvalues, nvec[s], replace = TRUE, prob = ddvalues))
  
  # Peer group for each agent in each subnetwork
  peerg    <- lapply(1:ngroup, function(s){
    lapply(1:nvec[s], function(i) sample((1:nvec[s])[-i], degree[[s]][i]))
  })    
  
  # Network matrix
  G        <- list()
  for (s in 1:ngroup) {
    Gs     <- matrix(0, nvec[s], nvec[s])
    for (i in 1:nvec[s]) {
      Gs[i, peerg[[s]][[i]]] <- 1
    }
    G[[s]] <- Gs
  }
  
  # Fixed effects
  eff      <- rep(runif(ngroup, 8, 10), nvec)
  # X
  X        <- cbind(rnorm(n, eff, 1), rpois(n, 2)) # Make the effects correlated with X
  # GX GGX
  Gnorm    <- norm.network(G) # Row-normalized network
  GX       <- peer.avg(Gnorm, X)
  
  # Outcome y
  y        <- NULL
  if (linear) {
    # If linear lambda is set to the sum of lambda
    y      <- linpeer.sim(formula = ~ -1 + eff + X + GX, Glist = Gnorm, lambda = sum(lambda), 
                          beta = c(ifelse(fixed.effects, 1, 0), beta), structural = FALSE, 
                          epsilon = rnorm(n, 0, sigma))$y 
  } else {
    y      <- qpeer.sim(formula = ~ -1 + eff + X + GX, Glist = G, tau = tau, lambda = lambda, 
                        beta = c(ifelse(fixed.effects, 1, 0), beta), structural = FALSE, 
                        epsilon = rnorm(n, 0, sigma))$y 
  }
  
  # Instruments
  # linear model: G^2X
  GGX      <- peer.avg(Gnorm, GX) #GX where G is row-normalized
  # Quantile model
  Z        <- qpeer.inst(formula = ~ X + GX, Glist = G, tau = seq(0, 1, 1/9), 
                         max.distance = 2, checkrank = TRUE)$instruments 
  # CES model: exogenous prediction of y
  Zces     <- NULL
  if (fixed.effects) {
    grp    <- as.factor(rep(1:ngroup, nvec))
    Zces   <- fitted.values(lm(y ~ X + GX + grp))
  } else {
    Zces   <- fitted.values(lm(y ~ X + GX))
  }
  
  # Linear model
  LIM     <- linpeer(formula = y ~ X + GX, 
                     excluded.instruments = ~ GGX, 
                     Glist = Gnorm,
                     structural = FALSE,
                     fixed.effects = fixed.effects)
  
  # Quantile-based model using Z1 as instruments
  Quant   <- qpeer(formula = y ~ X + GX, 
                   excluded.instruments = ~ Z, 
                   Glist = Gnorm,
                   structural = FALSE,
                   tau = tau,
                   fixed.effects = fixed.effects)
  
  # CES-based model
  Ces     <- cespeer(formula = y ~ X + GX, instrument = ~ Zces, Glist = Gnorm, 
                     structural = FALSE, fixed.effects = fixed.effects, radius = 5,
                     grid.rho = seq(-100, 100, 1))
  
  # Quantile model with 2 quantiles
  Quant2  <- qpeer(formula = y ~ X + GX, 
                   excluded.instruments = ~ Z, 
                   Glist = Gnorm,
                   structural = FALSE,
                   tau = seq(0, 1, 1),
                   fixed.effects = fixed.effects)
  
  # Quantile model with 3 quantiles
  Quant3  <- qpeer(formula = y ~ X + GX, 
                   excluded.instruments = ~ Z, 
                   Glist = Gnorm,
                   structural = FALSE,
                   tau = seq(0, 1, 1/2),
                   fixed.effects = fixed.effects)
  
  # Quantile model with 5 quantiles
  Quant5  <- qpeer(formula = y ~ X + GX, 
                   excluded.instruments = ~ Z, 
                   Glist = Gnorm,
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
  names(LIM)    <- paste0(ifelse(fixed.effects, "FE.", ""), "LIM.", names(LIM))
  
  Quant         <- summary(Quant, diagnostic = TRUE)
  Quant         <- c(Quant$gmm$Estimate, KPstat = Quant$diagnostics[ntau + 1, 3], 
                     Jpvalue = Quant$diagnostics[ntau + 3, 4], Jpvalue05 = Quant$diagnostics[ntau + 3, 4] < 0.05,
                     Jpvalue10 = Quant$diagnostics[ntau + 3, 4] < 0.1)
  names(Quant)  <- paste0(ifelse(fixed.effects, "FE.", ""), "Q.", names(Quant))
  
  Ces           <- summary(Ces)$gmm$Estimate
  names(Ces)    <- paste0(ifelse(fixed.effects, "FE.", ""), "CES.", names(Ces))
  
  Etest         <- c(Etest, Etest < 0.05, Etest < 0.1)
  names(Etest)  <- paste0(ifelse(fixed.effects, "FE.", ""), c("ET.ntau=2.3", "ET.ntau=3.4", "ET.ntau=4.5",
                                                              "ET05.ntau=2.3", "ET05.ntau=3.4", "ET05.ntau=4.5", 
                                                              "ET10.ntau=2.3", "ET10.ntau=3.4", "ET10.ntau=4.5"))
  
  Quant2        <- summary(Quant2)
  Quant2        <- c(Quant2$gmm$Estimate, KPstat = Quant2$diagnostics[ntau, 3], 
                     Jpvalue = Quant2$diagnostics[ntau + 2, 4], Jpvalue05 = Quant2$diagnostics[ntau + 2, 4] < 0.05, 
                     Jpvalue10 = Quant2$diagnostics[ntau + 2, 4] < 0.1)
  names(Quant2) <- paste0(ifelse(fixed.effects, "FE.", ""), "Q2.", names(Quant2))
  
  Quant3        <- summary(Quant3)
  Quant3        <- c(Quant3$gmm$Estimate, KPstat = Quant3$diagnostics[ntau + 2, 3], 
                     Jpvalue = Quant3$diagnostics[ntau + 4, 4], Jpvalue05 = Quant3$diagnostics[ntau + 4, 4] < 0.05,
                     Jpvalue10 = Quant3$diagnostics[ntau + 4, 4] < 0.1)
  names(Quant3) <- paste0(ifelse(fixed.effects, "FE.", ""), "Q3.", names(Quant3))
  
  Quant5        <- summary(Quant5)
  Quant5        <- c(Quant5$gmm$Estimate, KPstat = Quant5$diagnostics[ntau + 2, 3], 
                     Jpvalue = Quant5$diagnostics[ntau + 4, 4], Jpvalue05 = Quant5$diagnostics[ntau + 4, 4] < 0.05,
                     Jpvalue10 = Quant5$diagnostics[ntau + 4, 4] < 0.1)
  names(Quant5) <- paste0(ifelse(fixed.effects, "FE.", ""), "Qb.", names(Quant5))
  
  
  c(LIM, Quant, Ces, Etest, Quant2, Quant3, Quant5)
}

# Number of simulations
nsim      <- 1e3

# Several values will be tested for peer effects with conformity parameter set to 0.2
# Estimation is done in parallel 
nthreads  <- 30
set.seed(2025)  # Set global seed for reproducibility

# Increasing lambda
lambda1    <- c(0, 0.05, 0.2, 0.3)
cat("DGP A without fixed effects\n")
Est11      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda1, FALSE, FALSE, nthreads)
}))
cat("DGP A with fixed effects\n")
Est12      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda1, TRUE, FALSE, nthreads)
}))

# Decreasing lambda
lambda2    <- c(0.3, 0.2, 0.05, 0)
cat("DGP B without fixed effects\n")
Est21      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda2, FALSE, FALSE, nthreads)
}))
cat("DGP B with fixed effects\n")
Est22      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda2, TRUE, FALSE, nthreads)
}))

# Concave lambda
lambda3    <- c(0, 0.275, 0.275, 0)
cat("DGP C without fixed effects\n")
Est31      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda3, FALSE, FALSE, nthreads)
}))
cat("DGP C with fixed effects\n")
Est32      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda3, TRUE, FALSE, nthreads)
}))

# Convex lambda
lambda4    <- c(0.275, 0, 0, 0.275)
cat("DGP D without fixed effects\n")
Est41      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda4, FALSE, FALSE, nthreads)
}))
cat("DGP D with fixed effects\n")
Est42      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda4, TRUE, FALSE, nthreads)
}))

# Concave lambda
lambda5    <- c(-0.05, 0.35, 0.15, 0.1)
cat("DGP E without fixed effects\n")
Est51      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda5, FALSE, FALSE, nthreads)
}))
cat("DGP E with fixed effects\n")
Est52      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda5, TRUE, FALSE, nthreads)
}))

# Constant lambda with data simulated using the standard LIM model
lambda6    <- 0.55
cat("DGP F without fixed effects\n")
Est61      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda6, FALSE, TRUE, nthreads)
}))
cat("DGP F with fixed effects\n")
Est62      <- do.call(cbind, lapply(1:nsim, \(i) {
  cat("Iteration:", i, "\n")
  festim(lambda6, TRUE, TRUE, nthreads)
}))

save(Est11, Est12, Est21, Est22, Est31, Est32, Est41, Est42, Est51, Est52, Est61, Est62,
     file = paste0(OutResPath, "/Simulations.Rda"))

# Summary
# Summary function
Sumfunc    <- function(x) {
  c(mean = mean(x), sd = sd(x))
}

load(paste0(OutResPath, "/Simulations.Rda"))
Est        <- as.data.frame(rbind(cbind(apply(Est11, 1, Sumfunc), apply(Est12, 1, Sumfunc)),
                                  cbind(apply(Est21, 1, Sumfunc), apply(Est22, 1, Sumfunc)),
                                  cbind(apply(Est31, 1, Sumfunc), apply(Est32, 1, Sumfunc)),
                                  cbind(apply(Est41, 1, Sumfunc), apply(Est42, 1, Sumfunc)),
                                  cbind(apply(Est51, 1, Sumfunc), apply(Est52, 1, Sumfunc)),
                                  cbind(apply(Est61, 1, Sumfunc), apply(Est62, 1, Sumfunc))))

# Formatting to Excel
# Put standard error into parenthesis
Est <- Est %>% mutate(across(everything(), ~ sprintf("%.3f", .))) %>%  # Ensures exactly 3 decimals
  mutate(across(everything(), ~ ifelse(row_number() %% 2 == 0, paste0("(", trimws(.), ")"), .)))  # Parentheses without spaces 

# Insert blank rows between GDP
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
addWorksheet(wb, "Table 5.4") 
addWorksheet(wb, "Table 5.3") 
addWorksheet(wb, "Full results")

# Table 5.4
T5.4 <- Est %>% select(all_of(c(paste0("FE.Q.y_q", 1:4), "FE.LIM.G:y", "FE.CES.rho", "FE.CES.G:y")))
T5.4[seq(1, 16, 3), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1, collapse = ", "), ")$"),
                            paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2, collapse = ", "), ")$"),
                            paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3, collapse = ", "), ")$"),
                            paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4, collapse = ", "), ")$"),
                            paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5, collapse = ", "), ")$"),
                            paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$"))
T5.4 <- cbind(T5.4[,c("FE.Q.y_q1", "FE.Q.y_q2", "FE.Q.y_q3", "FE.Q.y_q4")], "V1" = NA,
              T5.4[,c("FE.LIM.G:y")], "V2" = NA,
              T5.4[,c("FE.CES.rho", "FE.CES.G:y")])
# write
writeData(wb, "Table 5.4", T5.4, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)
# first row
fr    <- c(paste("lmda_q", 1:4), "", "lmda", "", "rho", "lmda")
for (i in 1:length(fr)) {
  writeData(wb, "Table 5.4", fr[i], startRow = 1, startCol = i)
}
# Merge cells
for (i in seq(1, nrow(T5.4), 3) + 1) {
  mergeCells(wb, "Table 5.4", cols = 1:ncol(T5.4), rows = i)
}
# Formatting
addStyle(wb, "Table 5.4", style = createStyle(halign = "center"), 
         rows = 1:(nrow(T5.4) + 1), cols = 1:ncol(T5.4), gridExpand = TRUE)

# Table 5.3
T5.3 <- Est %>% select(all_of(c("FE.ET05.ntau=2.3", "FE.ET10.ntau=2.3",
                                "FE.ET05.ntau=3.4", "FE.ET10.ntau=3.4", 
                                "FE.ET05.ntau=4.5", "FE.ET10.ntau=4.5" # Encompassing tests
))) %>% 
  filter((row_number() - 1) %% 3 != 2)
T5.3[seq(1, 11, 2), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1, collapse = ", "), ")$"),
                            paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2, collapse = ", "), ")$"),
                            paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3, collapse = ", "), ")$"),
                            paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4, collapse = ", "), ")$"),
                            paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5, collapse = ", "), ")$"),
                            paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$"))
T5.3 <- cbind(T5.3[,c("FE.ET05.ntau=2.3", "FE.ET10.ntau=2.3")], 
              V1 = NA,
              T5.3[,c("FE.ET05.ntau=3.4", "FE.ET10.ntau=3.4")], 
              V2 = NA,
              T5.3[,c("FE.ET05.ntau=4.5", "FE.ET10.ntau=4.5")])
# write
writeData(wb, "Table 5.3", T5.3, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)
# first row
fr    <- c("2vs3", "2vs3", "", "3vs4", "3vs4", "", "4vs5", "4vs5")
for (i in 1:length(fr)) {
  writeData(wb, "Table 5.3", fr[i], startRow = 1, startCol = i)
}
# Merge cells
for (i in seq(1, nrow(T5.3), 2) + 1) {
  mergeCells(wb, "Table 5.3", cols = 1:ncol(T5.3), rows = i)
}
# Formatting
addStyle(wb, "Table 5.3", style = createStyle(halign = "center"), 
         rows = 1:(nrow(T5.3) + 1), cols = 1:ncol(T5.3), gridExpand = TRUE)

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
saveWorkbook(wb, paste0(OutResPath, "/simulations.xlsx"), overwrite = TRUE)
