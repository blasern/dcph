# Divisive Cover Persistent Homology

## Description

An R package for calculating persistent homology with the divisive cover algorithm described in:

> N. Blaser, M. Brun. Mathematics in Computer Science (2018). [Divisive Cover](https://doi.org/10.1007/s11786-018-0352-6).

## Installation

To install the latest version of this R package directly from github:

    install.packages("devtools")
    devtools::install_github("blasern/dcph")

## Usage

```{R}
# load package
require(dcph)

# generate sample data
rcircle <- function(N, r, sd){
  radius <- rnorm(N, r, sd)
  angle <- runif(N, 0, 2 * pi)
  data.frame(x = radius * cos(angle), 
             y = radius * sin(angle))
}
data_matrix <- rcircle(200, 1, .1)

# calculate divisive cover
dc <- divisive_cover(data = data_matrix,
                     division_fct = relative_gap_division(0.05), 
                     stop_fct = stop_relative_filter(0.01))

# calculate persistent homology
pers <- persistent_homology(dc)

# plot persistence diagram and barcode
plot_persistence(pers, mode = "diag")
plot_persistence(pers, mode = "bars")
```
