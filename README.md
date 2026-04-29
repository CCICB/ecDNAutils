
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ecDNAutils

<!-- badges: start -->

<!-- badges: end -->

The goal of ecDNAutils is to …

## Installation

You can install the development version of ecDNAutils from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("CCICB/ecDNAutils")
```

## Quick Start

``` r
library(ecDNAutils)

# Specify directory containing linx files
linx_dir <- system.file(package="ecDNAutils", "CHP212")

# Parse Key linx files
links <- read_links(linx_dir, sample = "CHP212")
#> Searching for file: /private/var/folders/d9/x2yygv_13_15dw5f8fspdn880000gp/T/Rtmp45lEEy/temp_libpath94562165630a/ecDNAutils/CHP212/CHP212.linx.links.tsv
#>     > File found
ecdna <- read_ecdna(linx_dir, sample = "CHP212")
#> Searching for file: /private/var/folders/d9/x2yygv_13_15dw5f8fspdn880000gp/T/Rtmp45lEEy/temp_libpath94562165630a/ecDNAutils/CHP212/CHP212.linx.ecdna.csv
#>     > File found

# Process to derive key-insights, like what genes are amplified on each ecDNA
ampgenes <- ecdna_dataframe_to_ampgenes(ecdna)
```
