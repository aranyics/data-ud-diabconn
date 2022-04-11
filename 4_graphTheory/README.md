
# Steps for 4_graphTheory

# 1. Run swp_batch_compute.m Matlab script

This computes small-world propensity based on (Muldoon 2016).
Muldoon SF, Bridgeford EW, Bassett DS. Small-World Propensity and Weighted Brain Networks. Sci Rep. 2016 Feb 25;6:22057.

Codes for SWP calculations can be downloaded from (https://complexsystemsupenn.com/codedata).

The script produces SWP_diab.csv and SWP_obes.csv files.

# 2. Run conngraphstat.R script in R

1. Specify the path where diabconntable.R is found in line 8 of conngraphstat.R
YOUR.SCRIPT.PATH = ''

2. Run script conngraphstat.R

Notes:
- required R packages
library(stringr)
library(igraph)
library(assortnet)
library(scales)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

- The network strength, assortativity and balance measures are implemented in calc_stat.R

- Plotting functions are found in connectivity_plots.R
