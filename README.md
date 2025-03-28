# QuantilePeer: An R Package for estimating Models with Quantile Peer Effects

Aristide Houndetoungan

## Introduction

The **QuantilePeer** package simulates and estimates quantile peer effect models introduced by [Houndetoungan (2025)](https://). Exact replication code for the results in the paper is available in the folder [**test**](https://github.com/ahoundetoungan/QuantilePeer/tree/main/test).  

See the [vignettes](https://github.com/ahoundetoungan/QuantilePeer/tree/main/vignettes) for detailed examples showing how to use the package.

The package also includes functions to simulate and estimate the CES-based peer effect model developed by [Boucher et al. (2024)](https://doi.org/10.3982/ECTA21048).

## Installation

Installation is possible from this GitHub repository using the following code:

``` r
library(remotes)
install_github("ahoundetoungan/QuantilePeer")
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

## How to use the **QuantilePeer** package

See [vignette](https://nbviewer.jupyter.org/github/ahoundetoungan/QuantilePeer/blob/master/doc/introduction_to_QuantilePeer.pdf).

## References
- Boucher, V., & Fortin, B. (2016). Some challenges in the empirics of the effects of networks. *Handbook on the Economics of Networks*, 45, 48, <[doi:10.1093/oxfordhb/9780199948277.013.22](https://doi.org/10.1093/oxfordhb/9780199948277.013.22)>
- Boucher, V., Rendall, M., Ushchev, P., & Zenou, Y. (2024). Toward a general theory of peer effects. *Econometrica*, 92(2), 543-565. <[doi:10.3982/ECTA21048](https://doi.org/10.3982/ECTA21048)>
- Houndetoungan, A. (2025). Quantile Peer Effect Models, *arXiv preprint arXiv:*. <[doi:]()>
- Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. *The American Statistician*, 50(4), 361-365. <[doi:10.1080/00031305.1996.10473566](https://doi.org/10.1080/00031305.1996.10473566)>


