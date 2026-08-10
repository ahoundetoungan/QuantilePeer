## QuantilePeer: An R package for Simulating and Estimation Quantile Peer Effect Models

**Aristide Houndetoungan**


 <!-- badges: start -->
  [![Lifecycle: stable](https://img.shields.io/badge/Lifecycle-Experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![R-CMD-check](https://github.com/ahoundetoungan/QuantilePeer/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ahoundetoungan/QuantilePeer/actions/workflows/R-CMD-check.yml)
  

  [![R-universe](https://ahoundetoungan.r-universe.dev/badges/QuantilePeer)](https://ahoundetoungan.r-universe.dev/QuantilePeer)
  [![CRAN](https://www.r-pkg.org/badges/version/QuantilePeer)](https://CRAN.R-project.org/package=QuantilePeer)
  [![DOI](https://img.shields.io/badge/DOI-10.32614%2FCRAN.package.QuantilePeer-blue)](https://doi.org/10.32614/CRAN.package.QuantilePeer)
  [![CRAN Downloads](https://img.shields.io/endpoint?url=https://ahoundetoungan.github.io/cranlogs/badges/QuantilePeer.json)](https://cran.r-project.org/package=QuantilePeer)


  [![Vignette](https://img.shields.io/badge/Vignette-blue.svg)](https://docs.google.com/viewer?url=https://github.com/ahoundetoungan/QuantilePeer/raw/main/doc/introduction_to_QuantilePeer.pdf)
<!-- badges: end -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/QuantilePeer)](https://cran.r-project.org/package=QuantilePeer) -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](DOI: 10.32614/CRAN.package) -->

### Overview

The **QuantilePeer** package implements estimation and simulation routines for quantile peer effect models introduced in [Houndetoungan (2025)](https://doi.org/10.48550/arXiv.2506.12920).
Replication code for all results in the paper is available in the [Replication](https://github.com/ahoundetoungan/QuantilePeer/tree/main/Replication) folder.

The package also includes functions to simulate and estimate the CES-based peer effect model developed by [Boucher et al. (2024)](https://doi.org/10.3982/ECTA21048).

## Installation
### CRAN version
**QuantilePeer** can be directly installed from CRAN.
```R
install.packages("QuantilePeer")
```

### GitHub version
It may be possible that we updated the package without submitting the new version to CRAN. The latest version (*but not necessary stable*) of **QuantilePeer** can be installed from this GitHub repos.
```R
remotes::install_github("ahoundetoungan/QuantilePeer", build_vignettes = TRUE)
```

### Getting Started
See the [vignettes](https://docs.google.com/viewer?url=https://github.com/ahoundetoungan/QuantilePeer/raw/main/vignettes/introduction_to_QuantilePeer.pdf) for detailed examples demonstrating how to use the package.
