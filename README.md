# Divisive cover algorithm

## Description

An R package for analyzing data with the divisive cover algorithm described in:

> N. Blaser, M. Brun (2017). Filtered covers.

## Installation

To install the stable version of this R package from CRAN:

    install.packages("dca", dependencies=TRUE)

To install the latest version of this R package directly from github:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("nello.blaser/dca")
    require(dca)