# QuantilePeer: An R Package for Estimating Models with Quantile Peer Effects

**Aristide Houndetoungan**

 <!-- badges: start -->
  [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![.github/workflows/R-CMD-check.yaml](https://github.com/ahoundetoungan/QuantilePeer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ahoundetoungan/QuantilePeer/actions/workflows/R-CMD-check.yaml)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  ![CRAN](https://img.shields.io/badge/CRAN-not%20yet-lightgrey)
  ![R-universe](https://img.shields.io/badge/R--universe-not%20yet-lightgrey)

  [![vignette](https://img.shields.io/badge/vignette-introduction-blue.svg)](https://nbviewer.org/github/ahoundetoungan/QuantilePeer/blob/main/doc/introduction_to_QuantilePeer.pdf)

<!-- badges: end -->



<!-- [![CRAN status](https://www.r-pkg.org/badges/version/QuantilePeer)](https://cran.r-project.org/package=QuantilePeer) -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](DOI: 10.32614/CRAN.package) -->

## Overview

The **QuantilePeer** package implements estimation and simulation routines for quantile peer effect models introduced in [Houndetoungan (2025)](https://). Replication code for all results in the paper is available in the [`test`](https://github.com/ahoundetoungan/QuantilePeer/tree/main/test) folder.

The package also includes functions to simulate and estimate the CES-based peer effect model developed by [Boucher et al. (2024)](https://doi.org/10.3982/ECTA21048).

## Installation

You can install the package directly from this GitHub repository using:

```r
library(remotes)
install_github("ahoundetoungan/QuantilePeer")
```


**Note for Windows users:** Make sure [**Rtools**](https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html) is installed before installing from GitHub.


## Getting Started
See the [vignettes](https://nbviewer.org/github/ahoundetoungan/QuantilePeer/blob/main/doc/introduction_to_QuantilePeer.pdf) for detailed examples demonstrating how to use the package.
