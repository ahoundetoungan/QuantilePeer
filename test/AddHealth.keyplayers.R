##############################################################################################################
##############################################################################################################
########################### Quantile Peer Effect Models by Aristide Houndetoungan ############################
##############################################################################################################

# Last updated: 2025-05-16

# This script reproduces the counterfactual analysis (Figure 6.1)

rm(list = ls()) 
library(dplyr)
library(QuantilePeer)
library(PartialNetwork)
library(foreach)
library(doParallel)
library(doRNG)
library(ggplot2)

OutDataPath <- "~/Dropbox/Data/AHdata/CleanData/" # Where prepared data for each outcome are saved (/ at the end is important)
OutResPath  <- "~/Dropbox/Academy/1.Papers/Quantile Peer Effects/Application/Outputs/" # Where results should be saved


# List of outcome variables
depvar  <- c("gpa", "academiceffort", "nclubs", "futureperception", "trouble", "smoke", "drink", "risky",
             "selfesteem", "physicalexercise", "fight")

# The following function takes an outcome and compute the impact of an increase in the intercept.
fsim    <- function(outcome) {
  if (file.exists(paste0(OutResPath, outcome, ".txt"))) {
    file.remove(paste0(OutResPath, outcome, ".txt"))
  }
  sink(file = paste0(OutResPath, outcome, ".txt"))
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
  
  # Load results
  load(paste0(OutResPath, outcome, ".Rda"))
  
  #############################################################
  ######################## Simulations ########################
  #############################################################
  nvec     <- sapply(G, nrow)
  nsch     <- length(nvec)
  ncum     <- c(0, cumsum(nvec))
  sch      <- (1:nsch)[nvec < 100]
  j        <- 0
  KPlayer  <- as.data.frame(matrix(NA, sum(nvec[sch]), 7))
  colnames(KPlayer) <- c("School", "indegree", "outdegree", "LIM", "Q3", "Q4", "CES")
  
  for (k in 1:length(sch)) {
    s          <- sch[k]
    for (i in 1:nvec[s]) {
      cat("School: ", k, "/", length(sch), " ** Student: ", i, "/", nvec[s], "\n", sep = "")
      Gs       <- G[[s]]
      idi      <- sum(Gs[,i])
      odi      <- sum(Gs[i,])
      Gs[i,]   <- 0
      Gs[,i]   <- 0
      Gs       <- list(norm.network(Gs))
      
      #### Linear model
      Xmat     <- cbind(X, GX)
      # multiply isolated individual X by 1 - lambda2
      Xmat[match > 0,] <- Xmat[match > 0,]*(1 - SLIM$gmm$Estimate["G(conformity):y"])
      # Gy
      Gy       <- peer.avg(Gnorm, y)
      # residuals
      res      <- y - cbind(Gy, Xmat) %*% SLIM$gmm$Estimate[c("G(total):y",
                                                              paste0("X", colnames(X)),
                                                              paste0("GX", colnames(GX)))]
      res[match > 0] <- res[match > 0]/(1 - SLIM$gmm$Estimate["G(conformity):y"])
      # key player status
      KPLIM    <- mean(y[(ncum[s] + 1):ncum[s + 1]] - 
                         linpeer.sim(formula = ~ -1 + X[(ncum[s] + 1):ncum[s + 1],] + GX[(ncum[s] + 1):ncum[s + 1],], 
                                     Glist = Gs,
                                     lambda = SLIM$gmm$Estimate[c("G(conformity):y", "G(total):y")],
                                     beta = c(SLIM$gmm$Estimate[c(paste0("X", colnames(X)),
                                                                  paste0("GX", colnames(GX)))]),
                                     structural = TRUE,
                                     epsilon = res[(ncum[s] + 1):ncum[s + 1]])$y)
      
      #### Quantile model
      # because this model is not linear in y, it is important to estimate the residual
      # tau3
      tau3     <- SQ13$model.info$tau
      Xmat     <- cbind(X, GX)
      # multiply isolated individual X by 1 - lambda2
      Xmat[match > 0,] <- Xmat[match > 0,]*(1 - SQ13$gmm$Estimate["y_q(conformity)"])
      # qy
      qy       <- qpeer.instrument(y ~ 1, tau = tau3, Glist = Gnorm)$qy
      # residuals
      res      <- y - cbind(qy, Xmat) %*% SQ13$gmm$Estimate[c(paste0("y_q", 1:length(tau3)),
                                                              paste0("X", colnames(X)),
                                                              paste0("GX", colnames(GX)))]
      res[match > 0] <- res[match > 0]/(1 - SQ13$gmm$Estimate["y_q(conformity)"])
      # key player status
      KPQ3     <- mean(y[(ncum[s] + 1):ncum[s + 1]] - 
                         qpeer.sim(formula = ~ -1 + X[(ncum[s] + 1):ncum[s + 1],] + GX[(ncum[s] + 1):ncum[s + 1],], 
                                   Glist = Gs,
                                   tau = tau3,
                                   lambda = SQ13$gmm$Estimate[c("y_q(conformity)", paste0("y_q", 1:length(tau3)))],
                                   beta = c(SQ13$gmm$Estimate[c(paste0("X", colnames(X)),
                                                                paste0("GX", colnames(GX)))]),
                                   structural = TRUE,
                                   init = y[(ncum[s] + 1):ncum[s + 1]],
                                   epsilon = res[(ncum[s] + 1):ncum[s + 1]])$y)
      
      # tau4
      tau4     <- SQ14$model.info$tau
      Xmat     <- cbind(X, GX)
      # multiply isolated individual X by 1 - lambda2
      Xmat[match > 0,] <- Xmat[match > 0,]*(1 - SQ14$gmm$Estimate["y_q(conformity)"])
      # qy
      qy       <- qpeer.instrument(y ~ 1, tau = tau4, Glist = Gnorm)$qy
      # residuals
      res      <- y - cbind(qy, Xmat) %*% SQ14$gmm$Estimate[c(paste0("y_q", 1:length(tau4)),
                                                              paste0("X", colnames(X)),
                                                              paste0("GX", colnames(GX)))]
      res[match > 0] <- res[match > 0]/(1 - SQ14$gmm$Estimate["y_q(conformity)"])
      # key player status
      KPQ4     <- mean(y[(ncum[s] + 1):ncum[s + 1]] - 
                         qpeer.sim(formula = ~ -1 + X[(ncum[s] + 1):ncum[s + 1],] + GX[(ncum[s] + 1):ncum[s + 1],], 
                                   Glist = Gs,
                                   tau = tau4,
                                   lambda = SQ14$gmm$Estimate[c("y_q(conformity)", paste0("y_q", 1:length(tau4)))],
                                   beta = c(SQ14$gmm$Estimate[c(paste0("X", colnames(X)),
                                                                paste0("GX", colnames(GX)))]),
                                   structural = TRUE,
                                   init = y[(ncum[s] + 1):ncum[s + 1]],
                                   epsilon = res[(ncum[s] + 1):ncum[s + 1]])$y)
      
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
      # key player status
      KPCES    <- mean(y[(ncum[s] + 1):ncum[s + 1]] - 
                         cespeer.sim(formula = ~ -1 + X[(ncum[s] + 1):ncum[s + 1],] + GX[(ncum[s] + 1):ncum[s + 1],], 
                                     Glist = Gs,
                                     rho = SCES$gmm$Estimate["rho"],
                                     lambda = SCES$gmm$Estimate[c("G(conformity):y", "G(total):y")],
                                     beta = c(SCES$gmm$Estimate[c(paste0("X", colnames(X)),
                                                                  paste0("GX", colnames(GX)))]),
                                     structural = TRUE,
                                     init = y[(ncum[s] + 1):ncum[s + 1]],
                                     epsilon = res[(ncum[s] + 1):ncum[s + 1]])$y, na.rm = TRUE)
      
      KPlayer[j,] <- c(s, idi, odi, KPLIM, KPQ3, KPQ4, KPCES)
      j           <- j + 1
    }
  }
  sink()
  
  saveRDS(KPlayer, file = paste0(OutResPath, outcome, ".RDS"))
}

# Set up parallel backend
ncores  <- 11
cl      <- makeCluster(ncores)
registerDoParallel(cl)

# Estimation
foreach(outcome = depvar, .packages = c("dplyr", "QuantilePeer", "PartialNetwork"), .options.RNG = 2025) %dorng% {
  fsim(outcome)
}

# Stop cluster
stopCluster(cl)

# # Label of Outcome
# OUTCOME <- c("Academic achievements", "Academic effort", "Extracurricular activities", "Future perception", 
#              "Trouble at school", "Smoking", "Drinking", "Risky behaviors", "Self-esteem", "Physical exercise", "Fighting")
# 
# # Scatterplot
# # data to plot
# dataplot <- do.call(rbind, lapply(1:length(depvar), function(k) {
#   cat("Outcome: ", OUTCOME[k], "\n", sep = "")
#   load(file = paste0(OutDataPath, depvar[k], ".Rda")) # Load saved data using with the script 0-Data.R
#   SIM     <-  readRDS(paste0(OutResPath, depvar[k], ".RDS")) # Load results
#   match   <- data$match
#   SIM     <- SIM[match > 0,] #only for non-isolated because the social multiplier is one for isolated
#   data.frame(outcome  = k,
#              model    = rep(1:2, each = nrow(SIM)),
#              Quant3   = rep(SIM$Q3, 2), 
#              Quant4   = rep(SIM$Q4, 2),
#              Other    = c(SIM$LIM, CES = SIM$CES))
# })) %>% mutate(outcome = factor(outcome, levels = 1:length(depvar), labels = OUTCOME),
#                model = factor(model, levels = 1:2, labels = c("Quantile vs LIM", "Quantile vs CES")))
# 
# # Select only 10%
# dataplot_sampled <- dataplot %>%
#   group_by(outcome, model) %>%
#   sample_frac(0.05)
# 
# # tau3
# (graph <- ggplot(dataplot_sampled, aes(x = Quant3, y = Other, colour = model, shape = model)) +
#     geom_point(size = 1.5, alpha = 1) +
#     facet_wrap(~ outcome, ncol = 2, scales = "free") +
#     scale_colour_manual(values = c("#c66", "#66c")) +
#     scale_shape_manual(values = c(2, 4)) +
#     theme_minimal(base_size = 12, base_family = "Palatino") +
#     theme(
#       strip.text = element_text(face = "bold"),
#       axis.text.x = element_text(hjust = 1),
#       panel.spacing = unit(1, "lines"),
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "bottom"
#     ) +
#     labs(
#       x = "Quantile Model",
#       y = "LIM and CES Models",
#       colour = NULL,
#       shape = NULL
#     ))
# ggsave("CFScatter3.pdf", path = OutResPath, plot = graph, device = "pdf", width = 7.5, height = 10)
# 
# # tau3
# (graph <- ggplot(dataplot_sampled, aes(x = Quant4, y = Other, colour = model, shape = model)) +
#     geom_point(size = 1.5, alpha = 1) +
#     facet_wrap(~ outcome, ncol = 2, scales = "free") +
#     scale_colour_manual(values = c("#c66", "#66c")) +
#     scale_shape_manual(values = c(2, 4)) +
#     theme_minimal(base_size = 12, base_family = "Palatino") +
#     theme(
#       strip.text = element_text(face = "bold"),
#       axis.text.x = element_text(hjust = 1),
#       panel.spacing = unit(1, "lines"),
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "bottom"
#     ) +
#     labs(
#       x = "Quantile Model",
#       y = "LIM and CES Models",
#       colour = NULL,
#       shape = NULL
#     ))
# ggsave("CFScatter4.pdf", path = OutResPath, plot = graph, device = "pdf", width = 7.5, height = 10)
# 
# # Boxplot
# # data to plot
# dataplot <- do.call(rbind, lapply(1:length(depvar), function(k) {
#   cat("Outcome: ", OUTCOME[k], "\n", sep = "")
#   load(file = paste0(OutDataPath, depvar[k], ".Rda")) # Load saved data using with the script 0-Data.R
#   SIM     <-  readRDS(paste0(OutResPath, depvar[k], ".RDS")) # Load results
#   match   <- data$match
#   SIM     <- SIM[match > 0,] #only for non-isolated because the social multiplier is one for isolated
#   sim     <- c(SIM$Q3, SIM$Q4, SIM$LIM, SIM$CES)
#   model   <- rep(1:4, each = nrow(SIM))
#   data.frame(outcome = k, model, sim)
# })) %>% mutate(outcome = factor(outcome, levels = 1:length(depvar), labels = OUTCOME),
#                modnum  = model,
#                model   = factor(model, levels = 1:4, 
#                                 labels = c(expression(paste("Q(", d[tau], " = 3)")),
#                                            expression(paste("Q(", d[tau], " = 4)")),
#                                            "LIM", "CES")))
# 
# (graph <- ggplot(dataplot, aes(x = model, y = sim)) +
#     geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +
#     facet_wrap(~ outcome, ncol = 3, scales = "free_y") +
#     theme_minimal(base_size = 12, base_family = "Palatino") +
#     theme(
#       strip.text = element_text(face = "bold"),
#       axis.text.x = element_text(hjust = 0.5),
#       panel.spacing = unit(1, "lines"),
#       plot.title = element_text(hjust = 0.5)
#     ) +
#     labs(
#       x = "Model",
#       y = "Social Multiplier",
#     ) +
#     scale_x_discrete(labels = scales::label_parse())) 
# 
# ggsave("CFBox.pdf", path = OutResPath, plot = graph, device = "pdf", width = 10, height = 6)
# 
