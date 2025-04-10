##############################################################################################################
##############################################################################################################
########################### Quantile Peer Effect Models by Aristide Houndetoungan ############################
##############################################################################################################

# Last updated: 2025-04-03

# This script replicates the counterfactual analysis (Figure 6.1)
rm(list = ls()) 
library(dplyr)
library(QuantilePeer)
library(PartialNetwork)
library(ggplot2)

OutDataPath <- "~/Dropbox/Data/AHdata/CleanData/" # Where prepared data for each outcome are saved (/ at the end is important)
OutResPath  <- "~/Dropbox/Academy/1.Papers/Quantile Peer Effects/Application/Outputs/" # Where results should be saved


# List of outcome variables
depvar  <- c("gpa", "academiceffort", "nclubs", "futureperception", "trouble", "smoke", "drink", "risky",
             "selfesteem", "physicalexercise", "fight")

# The following function takes an outcome and compute the impact of an increase in the intercept.
fsim    <- function(outcome) {
  cat("Outcome: ", outcome, "\n", sep = "")
  ########################################################
  ################### Data Preparation ###################
  ########################################################
  load(file = paste0(OutDataPath, outcome, ".Rda")) # Load saved data using with the script 0-Data.R
  # Exogenous characteristics (excluding reference variables for identification)
  exovar  <- c("age", "age2", "grade", "grade2", "female", "hispanic", "racewhite", "raceblack", 
               "raceasian", "melhigh", "memhigh", "memiss", "mjprof", "mjother",  "mjmiss")
  y       <- data$y + 1 # + 1 because the CES-based model does not allow y <= 0
  X       <- as.matrix(data[,exovar])
  nmatch  <- data$nmatch # Number of unmatched friends
  match   <- data$match  # Number of matched friends
  iso     <- match == 0 # Isolated in of observed network
  trueiso <- (match == 0 & nmatch == 0) # True isolated
  Gnorm   <- norm.network(G) # Row-normalized network
  GX      <- peer.avg(Gnorm, X) 
  X       <- cbind(X, nmatch = nmatch, match = match) # Adding additional variables
  
  # I test two vectors for tau
  tau1    <- seq(0, 1, 1/3) # Four quantiles
  tau2    <- seq(0, 1, 1/4) # Five quantiles
  
  # Load results
  load(paste0(OutResPath, outcome, ".Rda"))
  
  # Intercepts 1 and 2
  I1      <- rep(0, nrow(X))
  I2      <- I1 + 1
  
  #############################################################
  ######################## Simulations ########################
  #############################################################
  
  #### Linear model
  LIMsim  <- linpeer.sim(formula = ~ -1 + I2 + X + GX, 
                         Glist = Gnorm,
                         lambda = SLIM$gmm$Estimate[c("G(conformity):y", "G(total):y")],
                         beta = c(1, SLIM$gmm$Estimate[c(paste0("X", colnames(X)),
                                                         paste0("GX", colnames(GX)))]),
                         structural = TRUE,
                         epsilon = 0)$y - 
    linpeer.sim(formula = ~ -1 + I1 + X + GX, 
                Glist = Gnorm,
                lambda = SLIM$gmm$Estimate[c("G(conformity):y", "G(total):y")],
                beta = c(1, SLIM$gmm$Estimate[c(paste0("X", colnames(X)),
                                                paste0("GX", colnames(GX)))]),
                structural = TRUE,
                epsilon = 0)$y
  
  #### Quantile model
  # because this model is not linear in y, it is important to estimate the residual
  # tau1
  Xmat     <- cbind(X, GX)
  # multiply isolated individual X by 1 - lambda2
  Xmat[match > 0,] <- Xmat[match > 0,]*(1 - SQ11$gmm$Estimate["y_q(conformity)"])
  # qy
  qy       <- qpeer.instrument(y ~ 1, tau = tau1, Glist = Gnorm)$qy
  # residuals
  res      <- y - cbind(qy, Xmat) %*% SQ11$gmm$Estimate[c(paste0("y_q", 1:length(tau1)),
                                                          paste0("X", colnames(X)),
                                                          paste0("GX", colnames(GX)))]
  res[match > 0] <- res[match > 0]/(1 - SQ11$gmm$Estimate["y_q(conformity)"])
  
  Qsim1    <- qpeer.sim(formula = ~ -1 + I2 + X + GX, 
                        Glist = Gnorm,
                        tau = tau1,
                        lambda = SQ11$gmm$Estimate[c("y_q(conformity)", paste0("y_q", 1:length(tau1)))],
                        beta = c(1, SQ11$gmm$Estimate[c(paste0("X", colnames(X)),
                                                        paste0("GX", colnames(GX)))]),
                        structural = TRUE,
                        epsilon = res)$y - y # because y is the outcome before the increase in the intercept
  
  # tau2
  Xmat     <- cbind(X, GX)
  # multiply isolated individual X by 1 - lambda2
  Xmat[match > 0,] <- Xmat[match > 0,]*(1 - SQ21$gmm$Estimate["y_q(conformity)"])
  # qy
  qy       <- qpeer.instrument(y ~ 1, tau = tau2, Glist = Gnorm)$qy
  # residuals
  res      <- y - cbind(qy, Xmat) %*% SQ21$gmm$Estimate[c(paste0("y_q", 1:length(tau2)),
                                                          paste0("X", colnames(X)),
                                                          paste0("GX", colnames(GX)))]
  res[match > 0] <- res[match > 0]/(1 - SQ21$gmm$Estimate["y_q(conformity)"])
  
  Qsim2    <- qpeer.sim(formula = ~ -1 + I2 + X + GX, 
                        Glist = Gnorm,
                        tau = tau2,
                        lambda = SQ21$gmm$Estimate[c("y_q(conformity)", paste0("y_q", 1:length(tau2)))],
                        beta = c(1, SQ21$gmm$Estimate[c(paste0("X", colnames(X)),
                                                        paste0("GX", colnames(GX)))]),
                        structural = TRUE,
                        epsilon = res)$y - y # because y is the outcome before the increase in the intercept
  
  #### CES model
  # because this model is not linear in y, it is important to estimate the residual
  Xmat     <- cbind(X, GX)
  # multiply isolated individual X by 1 - lambda2
  Xmat[match > 0,] <- Xmat[match > 0,]*(1 - SCES$gmm$Estimate["G(conformity):y"])
  # social norm
  Gy       <- cespeer.data(y, Glist = Gnorm, rho = SCES$gmm$Estimate["rho"])[,"ces(y, rho)"]
  # residuals
  res      <- y - cbind(Gy, Xmat) %*% SCES$gmm$Estimate[c("G(total):y",
                                                          paste0("X", colnames(X)),
                                                          paste0("GX", colnames(GX)))]
  res[match > 0] <- res[match > 0]/(1 - SCES$gmm$Estimate["G(conformity):y"])
  
  CESsim   <- cespeer.sim(formula = ~ -1 + I2 + X + GX, 
                          Glist = Gnorm,
                          rho = SCES$gmm$Estimate["rho"],
                          lambda = SCES$gmm$Estimate[c("G(conformity):y", "G(total):y")],
                          beta = c(1, SCES$gmm$Estimate[c(paste0("X", colnames(X)),
                                                          paste0("GX", colnames(GX)))]),
                          structural = TRUE,
                          epsilon = res)$y - y # because y is the outcome before the increase in the intercept
  saveRDS(data.frame(LIM = LIMsim, Q1 = Qsim1, Q2 = Qsim2, CES = CESsim),
          file = paste0(OutResPath, outcome, ".RDS"))
}

# Estimation
set.seed(2025)
sapply(depvar, fsim) 

# Label of Outcome
OUTCOME <- c("Academic achievements", "Academic effort", "Extracurricular activities", "Future perception", 
             "Trouble at school", "Smoking", "Drinking", "Risky behaviors", "Self-esteem", "Physical exercise", "Fighting")

# Scatterplot
# data to plot
dataplot <- do.call(rbind, lapply(1:length(depvar), function(k) {
  cat("Outcome: ", OUTCOME[k], "\n", sep = "")
  load(file = paste0(OutDataPath, depvar[k], ".Rda")) # Load saved data using with the script 0-Data.R
  SIM     <-  readRDS(paste0(OutResPath, depvar[k], ".RDS")) # Load results
  match   <- data$match
  SIM     <- SIM[match > 0,] #only for non-isolated because the social multiplier is one for isolated
  data.frame(outcome  = k,
             model    = rep(1:2, each = nrow(SIM)),
             Quantile = rep(SIM$Q1, 2), 
             Other    = c(SIM$LIM, CES = SIM$CES))
})) %>% mutate(outcome = factor(outcome, levels = 1:length(depvar), labels = OUTCOME),
               model = factor(model, levels = 1:2, labels = c("Quantile vs LIM", "Quantile vs CES")))

# Select only 10%
dataplot_sampled <- dataplot %>%
  group_by(outcome, model) %>%
  sample_frac(0.05)

(graph <- ggplot(dataplot_sampled, aes(x = Quantile, y = Other, colour = model, shape = model)) +
    geom_point(size = 1.5, alpha = 1) +
    facet_wrap(~ outcome, ncol = 2, scales = "free") +
    scale_colour_manual(values = c("#c66", "#66c")) +
    scale_shape_manual(values = c(2, 4)) +
    theme_minimal(base_size = 12, base_family = "Palatino") +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(hjust = 1),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    labs(
      x = "Quantile Model",
      y = "LIM and CES Models",
      colour = NULL,
      shape = NULL
    ))
ggsave("CFScatter.pdf", path = OutResPath, plot = graph, device = "pdf", width = 7.5, height = 10)

# Boxplot
# data to plot
dataplot <- do.call(rbind, lapply(1:length(depvar), function(k) {
  cat("Outcome: ", OUTCOME[k], "\n", sep = "")
  load(file = paste0(OutDataPath, depvar[k], ".Rda")) # Load saved data using with the script 0-Data.R
  SIM     <-  readRDS(paste0(OutResPath, depvar[k], ".RDS")) # Load results
  match   <- data$match
  SIM     <- SIM[match > 0,] #only for non-isolated because the social multiplier is one for isolated
  sim     <- c(SIM$Q1, SIM$LIM, SIM$CES)
  model   <- rep(1:3, each = nrow(SIM))
  data.frame(outcome = k, model, sim)
})) %>% mutate(outcome = factor(outcome, levels = 1:length(depvar), labels = OUTCOME),
               model = factor(model, levels = 1:3, labels = c("Quantile", "LIM", "CES")))

(graph <- ggplot(dataplot, aes(x = model, y = sim)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +
    facet_wrap(~ outcome, ncol = 3, scales = "free_y") +
    theme_minimal(base_size = 12, base_family = "Palatino") +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(hjust = 1),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      x = "Model",
      y = "Social Multiplier",
    ))

ggsave("CFBox.pdf", path = OutResPath, plot = graph, device = "pdf", width = 7, height = 7)

