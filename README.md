# QtNet: An R Package for estimating Models with Quantile Peer Effects

Aristide Houndetoungan

## Introduction

The **QtNet** package includes all functions for the replication of the results in Houndetoungan (2025). The exact replication codes are located in the folder [**test**](https://github.com/ahoundetoungan/QtNet/tree/master/test). Below, we also provide detailed examples of how to use the estimators described in the paper.

## Installation

Installation is possible from this GitHub repository using the following code:

``` r
library(remotes)
install_github("ahoundetoungan/QtNet")
```

**Important:** Windows users must first install [**Rtools**](https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html) to enable installation from a GitHub repository.

## A brief description of the model

Let $\mathcal{N}$ be a set of $n$ agents indexed by the integer $i \in [1, ~n]$. Agents are connected through a network that is characterized by an adjacency matrix $\mathbf{G} = [g_{ij}]$ of dimension $n \times n$, where $g_{ij} = 1$ if agent $j$ is a friend of agent $i$, and $g_{ij} = 0$ otherwise. The network may be directed, meaning that $g_{ij} = 1$ does not imply that $g_{ji} = 1$. Self-connection is impossible; that is, $g_{ii} = 0$ for all $i$. In weighted networks, $g_{ij}$ can be a nonnegative variable (not necessarily binary) that measures the intensity of the outgoing link from $i$ to $j$. The model can also accommodate such networks. Note that the network is generally constituted in many independent subnets (eg: schools).

Let $\mathcal{T}$ be a set of quantile levels. The reduced-form specification of quantile peer effect models is given by:

``` math
y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + \boldsymbol{x}_i^{\prime}\beta + \varepsilon_i,
```

where $\mathbf{y} = (y_1, \ldots, y_{i-1}, y_{i+1}, \ldots, y_n)^{\prime}$ is the vector of outcomes for other units, and $q_{\tau,i}(\mathbf{y}_{-i})$ is the sample $\tau$-quantile of peer outcomes. The term $\varepsilon_i$ is an idiosyncratic error term, $\lambda_{\tau}$ captures the effect of the $\tau$-quantile of peer outcomes on $y_i$, and $\beta$ captures the effect of the exogenous variable $\boldsymbol{x}_i$ on $y_i$. For the definition of the sample $\tau$-quantile, see [Hyndman and Fan (1996)](https://doi.org/10.2307/2684934). If the network matrix is weighted, the sample weighted quantile can be used, where the outcome for friend $j$ of $i$ is weighted by $g_{ij}$.

One issue in linear peer effect models is that individual preferences with conformity and complementarity/substitution lead to the same reduced form. However, it is possible to disentangle both types of preferences using isolated individuals (see [Boucher and Fortin, 2016](https://doi.org/10.1093/oxfordhb/9780199948277.013.22), and [Boucher et al., 2024](https://doi.org/10.3982/ECTA21048)). Isolated individuals are those who have no friends. The structural specification of the model differs between isolated and nonisolated individuals.

For isolated $i$, the specification is similar to a standard linear-in-means model without social interactions, given by: $$y_i = \mathbf{x}_i^{\prime}\beta + \varepsilon_i.$$ If $i$ is non-isolated, the specification is given by: $$y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + (1 - \lambda^*)\mathbf{x}_i^{\prime}\beta  + \varepsilon_i,$$ where $\lambda^*$ determines whether preferences exhibit conformity or complementarity/substitution. In general, $\lambda^* > 0$ means that preferences are conformist (anti-conformity may be possible in some models when $\lambda^* < 0$). In contrast, when $\lambda^* = 0$, there is complementarity/substitution between individuals depending on the signs of the $\lambda_{\tau}$ parameters.

## How to use the **QtNet** package

In this section, I present the main functions of the package and how they can be used through examples.\
The main functions of the package include: - `qpeer.sim`: simulating data from models with quantile peer effects; - `qpeer.inst`: computing instruments for models with quantile peer effects; - `qpeer`, `linpeer`, and `genpeer`: estimating models with quantile peer effects.

Most of these functions are also classes that have `summary` and `print` methods.

Throughout this section, I use simulated data. To begin, I first create a network matrix `G` and two exogenous variables `X1` and `X2`.

library(QtNet)

ngr \<- 50 \# Number of subnets nvec \<- rep(30, ngr) \# Size of subnets n \<- sum(nvec)

### Simulating Data

## Network matrix

G \<- lapply(1:ngr, function(z) { Gz \<- matrix(rbinom(nvec[z]\^2, 1, 0.3), nvec[z], nvec[z]) diag(Gz) \<- 0 \# Adding isolated nodes (important for the structural model) niso \<- sample(0:nvec[z], 1, prob = (nvec[z] + 1):1 / sum((nvec[z] + 1):1)) if (niso \> 0) { Gz[sample(1:nvec[z], niso), ] \<- 0 } Gz })

tau \<- seq(0, 1, 1/3) X \<- cbind(rnorm(n), rpois(n, 2)) l \<- c(0.2, 0.15, 0.1, 0.2) b \<- c(2, -0.5, 1) eps \<- rnorm(n, 0, 0.4)

## Generating `y`

y \<- qpeer.sim(formula = \~ X, Glist = G, tau = tau, lambda = l, beta = b, epsilon = eps)\$y

### Estimation

## Computing instruments

Z \<- qpeer.inst(formula = \~ X, Glist = G, tau = seq(0, 1, 0.1), max.distance = 2, checkrank = TRUE) Z \<- Z\$instruments

## Reduced-form model

rest \<- qpeer(formula = y \~ X, excluded.instruments = \~ Z, Glist = G, tau = tau) summary(rest) summary(rest, diagnostic = TRUE) \# Summary with diagnostics

## Structural model

sest \<- qpeer(formula = y \~ X, excluded.instruments = \~ Z, Glist = G, tau = tau, structural = TRUE) summary(sest, diagnostic = TRUE) \# The lambda\^\* parameter is y_q (conformity) in the outputs. \# There is no conformity in the data, so the estimate will be approximately 0.

## Structural model with double fixed effects per subnet using optimal GMM

## and controlling for heteroskedasticity

sesto \<- qpeer(formula = y \~ X, excluded.instruments = \~ Z, Glist = G, tau = tau, structural = TRUE, fixed.effects = "separate", HAC = "hetero", gmm.weight = "optimal") summary(sesto, diagnostic = TRUE)

## Average peer effect model

# Row-normalized network to compute instruments

Gnorm \<- lapply(G, function(g) { d \<- rowSums(g) d[d == 0] \<- 1 g / d })

# GX and GGX

Gall \<- Matrix::bdiag(Gnorm) GX \<- as.matrix(Gall %*% X) GGX \<- as.matrix(Gall %*% GX)

# Standard linear model

lpeer \<- linpeer(formula = y \~ X + GX, excluded.instruments = \~ GGX, Glist = Gnorm) summary(lpeer, diagnostic = TRUE) \# Note: The normalized network is used here by definition of the model. \# Contextual effects are also included (this is also possible for the quantile model).

# The standard model can also be structural

lpeers \<- linpeer(formula = y \~ X + GX, excluded.instruments = \~ GGX, Glist = Gnorm, structural = TRUE, fixed.effects = "separate") summary(lpeers, diagnostic = TRUE)

## Estimation using `genpeer`

# Average peer variable computed manually and included as an endogenous variable

Gy \<- as.vector(Gall %\*% y) gpeer1 \<- genpeer(formula = y \~ X + GX, excluded.instruments = \~ GGX, endogenous.variables = \~ Gy, Glist = Gnorm, structural = TRUE, fixed.effects = "separate") summary(gpeer1, diagnostic = TRUE)

# Using both average peer variables and quantile peer variables as endogenous,

# or only the quantile peer variable

# Quantile peer `y`

qy \<- qpeer.inst(formula = y \~ 1, Glist = G, tau = tau) qy \<- qy\$qy

# Model estimation

gpeer2 \<- genpeer(formula = y \~ X + GX, excluded.instruments = \~ GGX + Z, endogenous.variables = \~ Gy + qy, Glist = Gnorm, structural = TRUE, fixed.effects = "separate") summary(gpeer2, diagnostic = TRUE)
