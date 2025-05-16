##############################################################################################################
##############################################################################################################
########################### Quantile Peer Effect Models by Aristide Houndetoungan ############################
##############################################################################################################

# Last updated: 2025-05-12

# This script reproduces the Monte Carlo simulations (Tables 5.1 and 5.2).
# Estimation results will be exported to a single Excel file (`Simulation.xlsx`)
# located in `OutResPath`. The Excel file contains the following sheets:
#  - Table 5.1 (reproduces Table 5.1),
#  - Table 5.2 (reproduces Table 5.2),
#  - Full results (not included in the paper).

rm(list = ls())
library(QuantilePeer)
library(CDatanet) #Another of my package to simulate and estimate peer effects models (count data, linear model, Tobit, ...)
library(parallel)
library(dplyr)
library(openxlsx)

# Where results should be saved
OutResPath  <- "~/Dropbox/Academy/1.Papers/Quantile Peer Effects/Simulations"

# I consider 2 explanatory variables and control for contextual variables
# parameter beta
beta       <- c(4, -0.5, 1, -0.2, 0.6)
# parameter sigma
sigma      <- 0.7

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
festim     <- function(lambda, fixed.effects = TRUE, linear = FALSE) {
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
  
  # X
  X        <- cbind(rnorm(n), rpois(n, 2))
  # GX GGX
  Gnorm    <- norm.network(G) # Row-normalized network
  GX       <- peer.avg(Gnorm, X)
  
  # Outcome y
  y        <- NULL
  if (linear) {
    # If linear lambda is set to the sum of lambda
    y      <- linpeer.sim(formula = ~ X + GX, Glist = Gnorm, lambda = sum(lambda), 
                          beta = beta, structural = FALSE, epsilon = rnorm(n, 0, sigma))$y 
  } else {
    y      <- qpeer.sim(formula = ~ X + GX, Glist = G, tau = tau, lambda = lambda, 
                          beta = beta, structural = TRUE, epsilon = rnorm(n, 0, sigma))$y 
  }
  
  # Instruments
  # linear model: G^2X
  GGX      <- peer.avg(Gnorm, GX) #GX where G is row-normalized
  # Quantile model
  Z1       <- qpeer.inst(formula = ~ X + GX, Glist = G, tau = seq(0, 1, 1/9), 
                         max.distance = 2, checkrank = TRUE)$instruments 
  Z2       <- qpeer.inst(formula = y ~ X + GX, Glist = G, tau = tau, 
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
                             structural = TRUE,
                             fixed.effects = fixed.effects)
  
  # Quantile-based model using Z1 as instruments
  Quant1  <- qpeer(formula = y ~ X + GX, 
                           excluded.instruments = ~ Z1, 
                           Glist = Gnorm,
                           structural = TRUE,
                           tau = tau,
                           fixed.effects = fixed.effects)
  
  # Quantile-based model using Z2 as instruments
  Quant2  <- qpeer(formula = y ~ X + GX, 
                           excluded.instruments = ~ Z2, 
                           Glist = Gnorm,
                           structural = TRUE,
                           tau = tau,
                           fixed.effects = fixed.effects)
  
  # Quantile-based model using Z1 nd Z2 as instruments
  Quant3  <- qpeer(formula = y ~ X + GX, 
                           excluded.instruments = ~ Z1 + Z2, 
                           Glist = Gnorm,
                           structural = TRUE,
                           tau = tau,
                           fixed.effects = fixed.effects)
  
  # CES-based model
  Ces     <- cespeer(formula = y ~ X + GX, instrument = ~ Zces, Glist = Gnorm, 
                             structural = TRUE, fixed.effects = fixed.effects, radius = 5,
                             grid.rho = seq(-100, 100, 1))
  
  # Quantile model with 3 quantiles
  Quant1a <- qpeer(formula = y ~ X + GX, 
              excluded.instruments = ~ Z1, 
              Glist = Gnorm,
              structural = TRUE,
              tau = seq(0, 1, 1/2),
              fixed.effects = fixed.effects)
  
  # Quantile model with 5 quantiles
  Quant1b <- qpeer(formula = y ~ X + GX, 
              excluded.instruments = ~ Z1, 
              Glist = Gnorm,
              structural = TRUE,
              tau = seq(0, 1, 1/4),
              fixed.effects = fixed.effects)
  
  # Validity of Z2 instruments
  valins   <-  c(Q2 = qpeer.test(Quant1, Quant2, which = "sargan")$test["p-value"],
                 Q3 = qpeer.test(Quant1, Quant3, which = "sargan")$test["p-value"])
  
  # Encompassing test
  Etest    <- c(qpeer.test(Quant1a, Quant1, which = "encompassing")$test["KP Wald rank", "p-value"], # testing where 4 quantiles is better than 3 quantiles
                qpeer.test(Quant1, Quant1a, which = "encompassing")$test["KP Wald rank", "p-value"], # testing where 3 quantiles is better than 4 quantiles
                qpeer.test(Quant1, Quant1b, which = "encompassing")$test["KP Wald rank", "p-value"]) # testing where 5 quantiles is better than 4 quantiles
  
  LIM           <- summary(LIM, diagnostic = TRUE)
  LIM           <- c(LIM$gmm$Estimate, KPstat = LIM$diagnostics[2, 3], Jpvalue = LIM$diagnostics[4, 4],
                     Jpvalue90 = LIM$diagnostics[4, 4] > 0.1, Jpvalue95 = LIM$diagnostics[4, 4] > 0.05)
  names(LIM)    <- paste0(ifelse(fixed.effects, "FE.", ""), "LIM.", names(LIM))
  
  Quant1        <- summary(Quant1, diagnostic = TRUE)
  Quant1        <- c(Quant1$gmm$Estimate, KPstat = Quant1$diagnostics[ntau + 1, 3], 
                     Jpvalue = Quant1$diagnostics[ntau + 3, 4], Jpvalue90 = Quant1$diagnostics[ntau + 3, 4] > 0.1,
                     Jpvalue95 = Quant1$diagnostics[ntau + 3, 4] > 0.05)
  names(Quant1) <- paste0(ifelse(fixed.effects, "FE.", ""), "Q1.", names(Quant1))
  
  Quant2        <- summary(Quant2, diagnostic = TRUE)
  Quant2        <- c(Quant2$gmm$Estimate, KPstat = Quant2$diagnostics[ntau + 1, 3], 
                     Jpvalue = Quant2$diagnostics[ntau + 3, 4], Jpvalue90 = Quant2$diagnostics[ntau + 3, 4] > 0.1,
                     Jpvalue95 = Quant2$diagnostics[ntau + 3, 4] > 0.05)
  names(Quant2) <- paste0(ifelse(fixed.effects, "FE.", ""), "Q2.", names(Quant2))
  
  Quant3        <- summary(Quant3, diagnostic = TRUE)
  Quant3        <- c(Quant3$gmm$Estimate, KPstat = Quant3$diagnostics[ntau + 1, 3], 
                     Jpvalue = Quant3$diagnostics[ntau + 3, 4], Jpvalue90 = Quant3$diagnostics[ntau + 3, 4] > 0.1, 
                     Jpvalue95 = Quant3$diagnostics[ntau + 3, 4] > 0.05)
  names(Quant3) <- paste0(ifelse(fixed.effects, "FE.", ""), "Q3.", names(Quant3))
  
  Ces           <- summary(Ces)$gmm$Estimate
  names(Ces)    <- paste0(ifelse(fixed.effects, "FE.", ""), "CES.", names(Ces))
  
  valins        <- c(valins, valins > 0.1, valins > 0.05)
  names(valins) <- paste0(ifelse(fixed.effects, "FE.", ""), c("Val.Q2", "Val.Q3",
                                                              "Val90.Q2", "Val90.Q3", 
                                                              "Val95.Q2", "Val95.Q3"))
    
  Etest         <- c(Etest, Etest > 0.1, Etest > 0.05)
  names(Etest)  <- paste0(ifelse(fixed.effects, "FE.", ""), c("ET.ntau=3.4", "ET.ntau=4.3", "ET.ntau=4.5",
                                                              "ET90.ntau=3.4", "ET90.ntau=4.3", "ET90.ntau=4.5", 
                                                              "ET95.ntau=3.4", "ET95.ntau=4.3", "ET95.ntau=4.5"))
  
  Quant1a       <- summary(Quant1a)
  Quant1a       <- c(Quant1a$gmm$Estimate, KPstat = Quant1a$diagnostics[ntau, 3], 
                     Jpvalue = Quant1a$diagnostics[ntau + 2, 4], Jpvalue90 = Quant1a$diagnostics[ntau + 2, 4] > 0.1, 
                     Jpvalue95 = Quant1a$diagnostics[ntau + 2, 4] > 0.05)
  names(Quant1a)<- paste0(ifelse(fixed.effects, "FE.", ""), "Q1a.", names(Quant1a))
  
  Quant1b       <- summary(Quant1b)
  Quant1b       <- c(Quant1b$gmm$Estimate, KPstat = Quant1b$diagnostics[ntau + 2, 3], 
                     Jpvalue = Quant1b$diagnostics[ntau + 4, 4], Jpvalue90 = Quant1b$diagnostics[ntau + 4, 4] > 0.1,
                     Jpvalue95 = Quant1b$diagnostics[ntau + 4, 4] > 0.05)
  names(Quant1b)<- paste0(ifelse(fixed.effects, "FE.", ""), "Q1b.", names(Quant1b))
  
  c(LIM, Quant1, Quant2, Quant3, Ces, valins, Etest, Quant1a, Quant1b)
}

# Number of simulations
nsim      <- 1e3

# Several values will be tested for peer effects with conformity parameter set to 0.2
# Estimation is done in parallel 
# mc.cores is the number of cored used (should be 1 for Windows users)
# set see
RNGkind("L'Ecuyer-CMRG")
set.seed(2025)  # Set global seed for reproducibility

# Increasing lambda
lambda1    <- c(0.2, 0, 0.05, 0.2, 0.3)
Est11      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda1, FALSE), mc.cores = 10))
Est12      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda1, TRUE), mc.cores = 10))

# Decreasing lambda
lambda2    <- c(0.2, 0.3, 0.2, 0.05, 0)
Est21      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda2, FALSE), mc.cores = 10))
Est22      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda2, TRUE), mc.cores = 10))

# Concave lambda
lambda3    <- c(0.2, 0, 0.275, 0.275, 0)
Est31      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda3, FALSE), mc.cores = 10))
Est32      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda3, TRUE), mc.cores = 10))

# Convex lambda
lambda4    <- c(0.2, 0.275, 0, 0, 0.275)
Est41      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda4, FALSE), mc.cores = 10))
Est42      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda4, TRUE), mc.cores = 10))

# Concave lambda
lambda5    <- c(0.2, -0.05, 0.35, 0.15, 0.1)
Est51      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda5, FALSE), mc.cores = 10))
Est52      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda5, TRUE), mc.cores = 10))

# Constant lambda with data simulated using the standard LIM model
lambda6    <- 0.55
Est61      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda6, FALSE, TRUE), mc.cores = 10))
Est62      <- do.call(cbind, mclapply(1:nsim, function(i) festim(lambda6, TRUE, TRUE), mc.cores = 10))

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
addWorksheet(wb, "Table 5.1") 
addWorksheet(wb, "Table 5.2") 
addWorksheet(wb, "Full results")

# Table 5.1
T5.1 <- Est %>% select(all_of(c(paste0("FE.Q1.y_q", 1:4), "FE.Q1.y_q(conformity)",
                              paste0("FE.Q3.y_q", 1:4), "FE.Q3.y_q(conformity)",
                              "FE.LIM.G(total):y", "FE.LIM.G(conformity):y", 
                              "FE.CES.rho", "FE.CES.G(total):y", "FE.CES.G(conformity):y")))
T5.1[seq(1, 16, 3), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda1[-1]) - lambda1[1], "$, ", 
                                 "$\\lambda_2 = ", lambda1[1], "$"),
                          paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda2[-1]) - lambda2[1], "$, ", 
                                 "$\\lambda_2 = ", lambda2[1], "$"),
                          paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda3[-1]) - lambda3[1], "$, ", 
                                 "$\\lambda_2 = ", lambda3[1], "$"),
                          paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda4[-1]) - lambda4[1], "$, ", 
                                 "$\\lambda_2 = ", lambda4[1], "$"),
                          paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda5[-1]) - lambda5[1], "$, ", 
                                 "$\\lambda_2 = ", lambda5[1], "$"),
                          paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$, $\\lambda_2 = 0$"))
T5.1 <- cbind(T5.1[,c("FE.Q1.y_q1", "FE.Q1.y_q2", "FE.Q1.y_q3", "FE.Q1.y_q4", "FE.Q1.y_q(conformity)")], "V1" = NA,
              T5.1[,c("FE.Q3.y_q1", "FE.Q3.y_q2", "FE.Q3.y_q3", "FE.Q3.y_q4", "FE.Q3.y_q(conformity)")], "V2" = NA,
              T5.1[,c("FE.LIM.G(total):y", "FE.LIM.G(conformity):y")], "V2" = NA,
              T5.1[,c("FE.CES.rho", "FE.CES.G(total):y", "FE.CES.G(conformity):y")])
# write
writeData(wb, "Table 5.1", T5.1, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)
# first row
fr    <- c(rep(c(paste("ld_q", 1:4), "ld2", ""), 2), "ld", "ld2", "", "rho", "ld", "ld2")
for (i in 1:length(fr)) {
  writeData(wb, "Table 5.1", fr[i], startRow = 1, startCol = i)
}
# Merge cells
for (i in seq(1, nrow(T5.1), 3) + 1) {
  mergeCells(wb, "Table 5.1", cols = 1:ncol(T5.1), rows = i)
}
# Formatting
addStyle(wb, "Table 5.1", style = createStyle(halign = "center"), 
         rows = 1:(nrow(T5.1) + 1), cols = 1:ncol(T5.1), gridExpand = TRUE)

# Table 5.2
T5.2 <- Est %>% select(all_of(c("FE.ET90.ntau=3.4", "FE.ET95.ntau=3.4", "FE.ET90.ntau=4.5", "FE.ET95.ntau=4.5", # Encompassing tests
                                "FE.Q1.KPstat", "FE.Q3.KPstat", # KP stat
                                "FE.Q3.Jpvalue90", "FE.Q3.Jpvalue95"))) %>% # Validity of Z2
  mutate(across(.cols = -c(FE.Q1.KPstat, FE.Q3.KPstat),
                .fns = ~ ifelse((row_number() - 1) %% 3 == 2, NA, .)))
T5.2[seq(1, 16, 3), 1] <- c("DGP A", "DGP B", "DGP C", "DGP D", "DGP E", "DGP F (LIM model)")
T5.2 <- cbind(T5.2[,c("FE.ET90.ntau=3.4", "FE.ET95.ntau=3.4", 
                      "FE.ET90.ntau=4.5", "FE.ET95.ntau=4.5")], "V1" = NA,
              T5.2[,c("FE.Q1.KPstat", "FE.Q3.KPstat", "FE.Q3.Jpvalue90", "FE.Q3.Jpvalue95")])
# write
writeData(wb, "Table 5.2", T5.2, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)
# first row
fr    <- c("3vs4", "3vs4", "4vs5", "4vs5", "", "KP Z1", "KP Z1, Z2", "Val Z2 5%", "Val Z2 10%")
for (i in 1:length(fr)) {
  writeData(wb, "Table 5.2", fr[i], startRow = 1, startCol = i)
}
# Merge cells
for (i in seq(1, nrow(T5.2), 3) + 1) {
  mergeCells(wb, "Table 5.2", cols = 1:ncol(T5.2), rows = i)
}
# Formatting
addStyle(wb, "Table 5.2", style = createStyle(halign = "center"), 
         rows = 1:(nrow(T5.2) + 1), cols = 1:ncol(T5.2), gridExpand = TRUE)

# Full results
tp <- Est
tp[seq(1, 16, 3), 1] <- c(paste0("DGP A: $\\boldsymbol\\lambda = (", paste0(lambda1[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda1[-1]) - lambda1[1], "$, ", 
                                 "$\\lambda_2 = ", lambda1[1], "$"),
                          paste0("DGP B: $\\boldsymbol\\lambda = (", paste0(lambda2[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda2[-1]) - lambda2[1], "$, ", 
                                 "$\\lambda_2 = ", lambda2[1], "$"),
                          paste0("DGP C: $\\boldsymbol\\lambda = (", paste0(lambda3[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda3[-1]) - lambda3[1], "$, ", 
                                 "$\\lambda_2 = ", lambda3[1], "$"),
                          paste0("DGP D: $\\boldsymbol\\lambda = (", paste0(lambda4[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda4[-1]) - lambda4[1], "$, ", 
                                 "$\\lambda_2 = ", lambda4[1], "$"),
                          paste0("DGP E: $\\boldsymbol\\lambda = (", paste0(lambda5[-1], collapse = ", "), ")$, ", 
                                 # "$\\lambda_1 = ", sum(lambda5[-1]) - lambda5[1], "$, ", 
                                 "$\\lambda_2 = ", lambda5[1], "$"),
                          paste0("DGP F (LIM model): $\\lambda = ", lambda6, "$, $\\lambda_2 = 0$"))
writeData(wb, "Full results", tp, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, paste0(OutResPath, "/simulations.xlsx"), overwrite = TRUE)
