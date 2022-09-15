## Packages to load
library(here)
library(tidyverse)
library(magrittr)
library(KFAS)
library(splines2)
library(invgamma)
library(mvtnorm)
# For parallelization
library(doParallel)
# For visualization
library(ggh4x)
library(ggrepel)
library(patchwork)
library(cowplot)
library(latex2exp)
library(ggflags)
library(plotly)

# The path of this file should be the working directory

# This code reports the order to execute the scripts of the project. 
# Some of the scripts may take sometime to run.
# Single execution of the script is possible as well.

### Create Rdata datasets
source(here('source','data_processing.R'))

### Fit BSP model and smoothing
source(here('application-countries','ITA_BSP.R'))
source(here('application-countries','SWE_BSP.R'))
source(here('application-countries','UK_BSP.R'))
source(here('application-countries','US_BSP.R'))
source(here('application-countries','smoothing.R'))
# Generate plots
source(here('application-countries','plot.R'))

### Rolling window and model comparison
# BSP model
source(here('application-countries','model-comparison','ITA_BSP_for.R'))
source(here('application-countries','model-comparison','SWE_BSP_for.R'))
source(here('application-countries','model-comparison','UK_BSP_for.R'))
source(here('application-countries','model-comparison','US_BSP_for.R'))
# Alternative models for comparison
# (Computationally intensive part due to rolling window)
# (these requires further libraries loaded in the respective scripts)
source(here('application-countries','model-comparison','LC_for.R'))
source(here('application-countries','model-comparison','APC_for.R'))
source(here('application-countries','model-comparison','CBD_for.R'))
source(here('application-countries','model-comparison','PLAT_for.R'))
source(here('application-countries','model-comparison','RH_for.R'))
source(here('application-countries','model-comparison','HU_for.R'))
source(here('application-countries','model-comparison','CP_for.R'))
# Table of results
source(here('application-countries','model-comparison','results.R'))

