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

Let $\mathcal{N}$ be a set of $n$ agents indexed by the integer $i \in [1, ~n]$. Agents are connected through a network that is characterized by an adjacency matrix $\mathbf{G} = [g_{ij}]$ of dimension $n \times n$, where $g_{ij} = 1$ if agent $j$ is a friend of agent $i$, and $g_{ij} = 0$ otherwise. In weighted networks, $g_{ij}$ can be a nonnegative variable (not necessarily binary) that measures the intensity of the outgoing link from $i$ to $j$. The model can also accommodate such networks. Note that the network is generally constituted in many independent subnets (eg: schools).

Let $\mathcal{T}$ be a set of quantile levels. The reduced-form specification of quantile peer effect models is given by:

``` math
y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + \boldsymbol{x}_i^{\prime}\beta + \varepsilon_i,
```

where $`\mathbf{y}_{-i} = (y_1, \ldots, y_{i-1}, y_{i+1}, \ldots, y_n)^{\prime}`$ is the vector of outcomes for other units, and $`q_{\tau,i}(\mathbf{y}_{-i})`$ is the sample $\tau$-quantile of peer outcomes. The term $`\varepsilon_i`$ is an idiosyncratic error term, $`\lambda_{\tau}`$ captures the effect of the $`\tau`$-quantile of peer outcomes on $`y_i`$, and $`\beta`$ captures the effect of the exogenous variable $`\boldsymbol{x}_i`$ on $`y_i`$. For the definition of the sample $`\tau`$-quantile, see [Hyndman and Fan (1996)](https://doi.org/10.2307/2684934). If the network matrix is weighted, the sample weighted quantile can be used, where the outcome for friend $`j`$ of $`i`$ is weighted by $`g_{ij}`$.

One issue in linear peer effect models is that individual preferences with conformity and complementarity/substitution lead to the same reduced form. However, it is possible to disentangle both types of preferences using isolated individuals (see [Boucher and Fortin, 2016](https://doi.org/10.1093/oxfordhb/9780199948277.013.22), and [Boucher et al., 2024](https://doi.org/10.3982/ECTA21048)). Isolated individuals are those who have no friends. The structural specification of the model differs between isolated and nonisolated individuals.

For isolated $`i`$, the specification is similar to a standard linear-in-means model without social interactions, given by: 
``` math
y_i = \mathbf{x}_i^{\prime}\beta + \varepsilon_i.
```
If $`i`$ is non-isolated, the specification is given by: 
``` math
y_i = \sum_{\tau \in \mathcal{T}} \lambda_{\tau} q_{\tau,i}(\mathbf{y}_{-i}) + (1 - \lambda^*)\mathbf{x}_i^{\prime}\beta  + \varepsilon_i,
```
where $`\lambda^*`$ determines whether preferences exhibit conformity or complementarity/substitution. In general, $`\lambda^* > 0`$ and this means that that preferences are conformist (anti-conformity may be possible in some models when $`\lambda^* < 0`$). In contrast, when $`\lambda^* = 0`$, there is complementarity/substitution between individuals depending on the signs of the $`\lambda_{\tau}`$ parameters. It is obvious that $`\beta`$ and $`\lambda^*`$ can be identified only if the network includes enough isolated individuals.

## How to use the **QtNet** package

In this section, I present the main functions of the package and how they can be used through examples. The main functions of the package include:
- `qpeer.sim`: simulating data from models with quantile peer effects;
- `qpeer.inst`: computing instruments for models with quantile peer effects;
- `qpeer`, `linpeer`, and `genpeer`: estimating models with quantile peer effects.

Most of these functions are also classes that have `summary` and `print` methods.

### Data simulation
Throughout this section, I use simulated data. To begin, I first create a network matrix `G` and two exogenous variables, `X1` and `X2`. Importantly, I impose some isolated nodes for the identification of the structural model. Otherwise, only the reduced-form parameters can be identified.
```R
library(QtNet)
ngr  <- 50  # Number of subnets
nvec <- rep(30, ngr)  # Size of subnets
n    <- sum(nvec)

# Network matrix
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

X   <- cbind(rnorm(n), rpois(n, 2)); colnames(X) <- c("X1", "X2")
```

Using the network matrix and the exogenous variables, I can now generate the dependent variables. I consider the quantile levels `seq(0, 1, 1/3)` (which defines four quantiles) and two dependent variables. The first one is based on the reduced-form model (where there is no conformity parameter), while the second is based on the structural model. The following code assigns values to the parameters of the model, followed by simulating the dependent variables. 
```R
tau       <- seq(0, 1, 1/3) #quantile level
lambdatau <- c(0.1, 0.25, 0.2, 0.15) #lambda_tau
lambdast  <- 0.2 # lambda_start
beta      <- c(2, -0.5, 1)

# First dependent variable (reduced form without conformity)
y1        <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = lambdatau, 
                       structural = FALSE, beta = beta, epsilon = rnorm(n, 0, 0.4)) 
y1        <- y1$y #qpeer.sim returns a list of several object including y

# Second dependent variable (structural form with conformity)
y2        <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = c(lambdast, lambdatau), 
                       structural = TRUE, beta = beta, epsilon = rnorm(n, 0, 0.4)) 
y2        <- y2$y 
```
Note that we can also include contextual variables, such as averages of `X` among peers, as additional exogenous variables. In real-life situations, the practitioner would also have their networks, exogenous variables `X`, and dependent variables. In the next section, I illustrate how quantile peer effects can be estimated. The same process can also be applied to real-life data instead of the simulated data I am using.


### Estimation of quantile peer effects
As discussed in the paper, there are two possible instrument sets for quantile peer outcome. The first instrument set includes quantile of `X` among peers. However, one does not need to keep the same quantile level as in the model. To enhance the instrument strength, I recommend using finer quantile levels than those used in the model. The second instrument set is still quantile of `X` among peers, but with the key difference that the values of `X` of peers are ordered using the values of their dependent variable. 

Both instrument sets can be computed using `qpeer.inst`. In the `formula` argument, if a dependent variable is specified, then the second instrument set is computed. 
```R
# First instrument set
Z1  <- qpeer.inst(formula = ~ X, Glist = G, tau = seq(0, 1, 0.1),  
                  max.distance = 2, checkrank = TRUE) # finer subdivision of quantile levels
Z1  <- Z1$instruments #qpeer.inst returns a list of several object including instruments 

# Second instrument set: y1 is used to order X
Z21 <- qpeer.inst(formula = y1 ~ X, Glist = G, tau = seq(0, 1, 0.1),  
                  max.distance = 2, checkrank = TRUE) 
qy1 <- Z21$qy #quantile of y among peers                 
Z21 <- Z21$instruments

# Second instrument set: y2 is used to order X
Z22 <- qpeer.inst(formula = y2 ~ X, Glist = G, tau = seq(0, 1, 0.1),  
                  max.distance = 2, checkrank = TRUE) 
qy2 <- Z22$qy #quantile of y among peers                 
Z22 <- Z22$instruments
```
As in the standard linear model, it is possible to use the quantile of direct friends' and long-distance `X` (such as friends' friends) to strengthen the instruments. I set `max.distance = 2`, which means that I use the quantiles of `X` among direct friends and friends' friends. The `checkrank` argument ensures that the resulting instrument set is a full-rank matrix by removing columns that are linear combinations of others.
