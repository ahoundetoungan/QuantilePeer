---
output:
  html_document: default
  pdf_document: default
---
# QtNet: An R Package for estimating Models with Quantile Peer Effects

Aristide Houndetoungan

## Introduction

The **QtNet** package simulates and estimates quantile peer effect models that were introduced by Houndetoungan (2025). The exact replication codes of the results in the paper are located in the folder [**test**](https://github.com/ahoundetoungan/QtNet/tree/master/test). Below, I also provide detailed examples of how to use the package.

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
y1        <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = lambdatau, beta = beta, 
                       structural = FALSE, epsilon = rnorm(n, 0, 0.4)) 
y1        <- y1$y #qpeer.sim returns a list of several object including y

# Second dependent variable (structural form with conformity)
y2        <- qpeer.sim(formula = ~ X, Glist = G, tau = tau, lambda = c(lambdast, lambdatau), 
                       beta = beta, structural = TRUE, epsilon = rnorm(n, 0, 0.4)) 
y2        <- y2$y 
```
Note that we can also include contextual variables, such as averages of `X` among peers, as additional exogenous variables. 

### Instruments

As discussed in the paper, there are two possible instrument sets for quantile peer outcomes. The first type of instruments is the quantiles of `X` among peers. The second type is also the quantiles of `X` among peers, but with a key distinction that the values of `X` for peers are ordered using the values of their dependent variable. The second type of instruments can lead to the most efficient estimators.

Instruments can be computed using `qpeer.inst`. The choice of the instrument type can be specified through the `formula` argument. If `formula` is defined without a dependent variable (i.e., an expression of the type `~ X1 + X2 + ...`), then the first type is computed. In contrast, if `formula` is defined with a dependent variable (i.e., an expression of the type `y ~ X1 + X2 + ...`), the second type is computed. Importantly, it is not necessary to use the same quantile level as in the model, especially for the first instrument type. To enhance instrument strength, I recommend using finer quantile levels than those employed in the model.

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

### Estimation of quantile peer effects
Quantile peer effects are estimated using the general method of moments (GMM). The estimates can be obtained by employing the `qpeer` function. I begin with the reduced-form specification using both types of instruments for each dependent variable.
```R
M1 <- qpeer(formula = y1 ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau)
summary(M1, diagnostic = TRUE)

M2 <- qpeer(formula = y1 ~ X, excluded.instruments = ~ Z21, Glist = G, tau = tau)
summary(M2, diagnostic = TRUE)

M3 <- qpeer(formula = y2 ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau)
summary(M3, diagnostic = TRUE)

M4 <- qpeer(formula = y2 ~ X, excluded.instruments = ~ Z21, Glist = G, tau = tau)
summary(M4, diagnostic = TRUE)

```
Note that the `formula` argument does not include the quantile of peers, but only exogenous variables. The quantiles will be directly computed by the function given the quantile levels in the `tau` argument. Additionally, `excluded.instruments` should not include all instruments, but only those not used as exogenous variables.

The output of the `qpeer` function is a `class` to which one can apply a `summary` method. An important argument of the `summary` function is `diagnostics`, a logical value to specify whether diagnostic tests for the instrumental-variable regression should be performed. These tests include an F-test of the first-stage regression for weak instruments, a Wu-Hausman test for endogeneity, and a Hansen's J-test for overidentifying restrictions (only if there are more instruments than regressors).

The estimation results for `y2` are likely biased (even for the exogenous variable) because, for this dependent variable, preferences exhibit both conformity and complementarity, which make the coefficients of the exogenous variable different for isolated and non-isolated individuals. However, in the above code, the reduced-form specification is considered without distinction between the coefficient of `X` for isolated and non-isolated individuals (the underlying assumption is that preferences exhibit either conformity or complementarity). The following code replicates the same estimation using the structural specification.

```R
M5 <- qpeer(formula = y1 ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau,
            structural = TRUE)
summary(M5, diagnostic = TRUE)

M6 <- qpeer(formula = y1 ~ X, excluded.instruments = ~ Z21, Glist = G, tau = tau,
            structural = TRUE)
summary(M6, diagnostic = TRUE)

M7 <- qpeer(formula = y2 ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau,
            structural = TRUE)
summary(M7, diagnostic = TRUE)

M8 <- qpeer(formula = y2 ~ X, excluded.instruments = ~ Z22, Glist = G, tau = tau,
            structural = TRUE)
summary(M8, diagnostic = TRUE)
```
In the new results, the estimates seem reliable. For the dependent variable `y1`, the conformity parameter is not significant because the data are simulated by assuming complementarity.

The `qpeer` function offers several useful options, including changing the weight of the GMM estimator, controlling for subnet fixed effects, and accounting for heteroskedasticity. The GMM weight can be controlled through the `gmm.weight` argument. The default value `"IV"` corresponds to the standard instrumental variable (IV) weight. It is also possible to use the identity matrix (`gmm.weight = "ident"`) and the optimal GMM weight (`gmm.weight = "optimal"`).

The `fixed.effects` argument can be used to specify how to control for subnet fixed effects. The default value is `FALSE` or `"no"` to indicate that there are no fixed effects. Two levels of subnet fixed effects are possible: a single fixed effect per subnet (`fixed.effects = "join"`) and double fixed effects per subnet (`fixed.effects = "separate"`), each for isolated and non-isolated individuals (see [Houndetoungan et al., 2024](https://doi.org/10.48550/arXiv.2405.06850)). For the structural specification, the fixed effects are necessarily double per subnet.

The `HAC` argument can indicate the covariance structure of errors. The default value assumes homoscedasticity (`HAC = "iid"`). To control for heteroskedasticity at the individual level, `HAC` can be set to `"hetero"`. To control for heteroskedasticity at the subnet level, where the errors associated with individuals in the same subnet can be correlated, `HAC` can be set to `"cluster"`.

```R
M9  <- qpeer(formula = y1 ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau,
             structural = FALSE, gmm.weight = "optimal", fixed.effects = "separate", 
             HAC = "hetero")
summary(M9, diagnostic = TRUE)

M10 <- qpeer(formula = y2 ~ X, excluded.instruments = ~ Z1, Glist = G, tau = tau,
             structural = TRUE, gmm.weight = "optimal", fixed.effects = "separate", 
             HAC = "hetero")
summary(M10, diagnostic = TRUE)
```


### Other specifications
The package also provides a method for estimating the standard linear peer effect model, where the average of the outcome among peers serves as an explanatory variable. The `linpeer` function performs this estimation. The function arguments are similar to the `qpper` function (except for the quantile-specific arguments such as `tau` and `type`).  

For the standard linear model, instruments can be `GX` when the model does not include contextual effects or `G^2X` otherwise, where `G` is row-normalized so that `GX` can be interpreted as averages of `X` among peers (see [Bramoullé et al., 2009](http://dx.doi.org/10.1016/j.jeconom.2008.12.021)). It is also possible to include higher powers of `G` to enhance instrument strength (see [Houndetoungan and Maoude, 2024](https://doi.org/10.48550/arXiv.2402.05030)). To normalize the network and compute `GX`, I use the functions `norm.network` and `peer.avg` from the [**PartialNetwork**](https://cran.r-project.org/package=PartialNetwork) package. In the following code, I use both `GX` and `G^2X` as instruments.  
```R
library(PartialNetwork)
Gn  <- norm.network(G)
GX  <- peer.avg(Gn, X)
GGX <- peer.avg(Gn, GX)
M11 <- linpeer(formula = y2 ~ X, excluded.instruments = ~ GX + GGX, Glist = Gn, 
               structural = TRUE, gmm.weight = "optimal", fixed.effects = "separate", 
               HAC = "hetero")
summary(M11, diagnostic = TRUE)
```
The conformity parameter appears consistent, as the estimate is similar to what was obtained from the previous estimations using the quantile model. However, the total effect seems to correspond to the sum of the estimates across different quantiles from the previous estimations.  

In the `qpeer` and `linpeer` functions, the endogenous variables are not defined by the user. For example, in the quantile model, one only needs to set the quantile levels. This feature can limit the estimation of certain specifications. To address this issue, the `QtNet` package provides the generic function `genpeer`, which allows users to define their endogenous variables. This option is useful for estimating the effects of quantiles among girl and boy peers, combining average peer effects with quantile peer effects in the same model, or specifying other social norms.  

The `endogenous.variables` argument of the function must be defined as a formula where the endogenous variables are specified. Below is an example where average peer effects and quantile peer effects are estimated.  
```R
Gy2 <- peer.avg(Gn, y2)
qy2 <- qpeer.inst(formula = y2 ~ 1, Glist = G, tau = tau[1])$qy 
M12 <- genpeer(formula = y2 ~ X, excluded.instruments = ~ Z1 + GX + GGX, 
               endogenous.variables = ~ Gy2 + qy2, Glist = G, structural = TRUE, 
               gmm.weight = "optimal", fixed.effects = "separate", HAC = "hetero")
summary(M12, diagnostic = TRUE)
```
## References
- Boucher, V., & Fortin, B. (2016). Some challenges in the empirics of the effects of networks. *Handbook on the Economics of Networks*, 45, 48, <[doi:10.1093/oxfordhb/9780199948277.013.22](https://doi.org/10.1093/oxfordhb/9780199948277.013.22)>
- Boucher, V., Rendall, M., Ushchev, P., & Zenou, Y. (2024). Toward a general theory of peer effects. *Econometrica*, 92(2), 543-565. <[doi:10.3982/ECTA21048](https://doi.org/10.3982/ECTA21048)>
- Bramoullé, Y., Djebbari, H., & Fortin, B. (2009). Identification of peer effects through social networks. Journal of econometrics, 150(1), 41-55. <[doi:10.1016/j.jeconom.2008.12.021](https://doi.org/10.1016/j.jeconom.2008.12.021)>
- Houndetoungan, A. (2025). Quantile Peer Effect Models, *arXiv preprint arXiv:*. <[doi:]()>
- Houndetoungan, A., Kouame, C., & Vlassopoulos, M. (2024). Identifying Peer Effects in Networks with Unobserved Effort and Isolated Students, *arXiv preprint arXiv:2405.06850*. <[doi:10.48550/arXiv.2405.06850](https://doi.org/10.48550/arXiv.2405.06850)>
- Houndetoungan, A., & Maoude, A. H. (2024). Inference for Two-Stage Extremum Estimators, *arXiv preprint arXiv:2402.05030*. <[doi:10.48550/arXiv.2402.05030](https://doi.org/10.48550/arXiv.2402.05030)>
- Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. *The American Statistician*, 50(4), 361-365. <[doi:10.1080/00031305.1996.10473566](https://doi.org/10.1080/00031305.1996.10473566)>


