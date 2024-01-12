# This code reports the order to execute the scripts
# to generate the results for the extra folder
# It requires to have installed all the packages in the related scripts

### Create Rdata datasets
source(here('source','data_processing_extra.R'))

### Rolling window and model comparison
# BSP model
source(here('application-countries','model-comparison','extra','FRA_BSP_for.R'))
source(here('application-countries','model-comparison','extra','DNK_BSP_for.R'))
source(here('application-countries','model-comparison','extra','CZE_BSP_for.R'))
# Alternative models for comparison
# (Computationally intensive part due to rolling window)
# (these requires further libraries loaded in the respective scripts)
source(here('application-countries','model-comparison','extra','LC_for.R'))
source(here('application-countries','model-comparison','extra','APC_for.R'))
source(here('application-countries','model-comparison','extra','CBD_for.R'))
source(here('application-countries','model-comparison','extra','PLAT_for.R'))
source(here('application-countries','model-comparison','extra','RH_for.R'))
source(here('application-countries','model-comparison','extra','HU_for.R'))
source(here('application-countries','model-comparison','extra','CP_for.R'))
# Post-processing 
source(here('application-countries','model-comparison','extra','results_extra.R'))

# Tables aggregating the results for the extra countries are generated in
# BSP-mortality/application-countries/model-comparison/results.R
# by setting flag_extra_countries <- TRUE