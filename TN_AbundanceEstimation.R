###############################################################################L

#   TN Abundance Estimation Analysis
#       
#       Running this script will use harvest data to estimate abundance of deer 
#       in Tennessee. The steps and dependencies involved include: 
#           1. Sourcing an R script (Rscripts/TN_harvest_data_prep.R) to format 
#                and prep the harvest data (Data/dkmaster.RData; Data/dmus.csv)
#           2. Sourcing an R script (Rscripts/Functions_AbundanceEstimation.R) 
#                 to load functions used in the analysis 
#           3. Create priors and specify other Bayesian requirements
#           4. Run the IPM model file (TN_IPM_RecRate.txt)
#   
#   Allison C Keever
#   TNTech University
#   Last Updated: 10/17/24

###############################################################################L

# Prep work space---------------------------------------------------------------

# Load packages needed in analyses
library(tidyverse)
library(jagsUI)
library(abind)

# Run script to load data and create data frames for analyses
source("Rscripts/TN_harvest_data_prep.R")

# Source functions needed for analyses
source("Rscripts/Functions_AbundanceEstimation.R")


# Bayesian requirements---------------------------------------------------------

# Estimate initial values for population size based on reported harvest, 
# reporting rates, and harvest rates. 
init_N <- est_init_pop(y.H = y_H_DMU, y.Af = y_Af_DMU, 
                       y.Am = y_Am_DMU, harvest.rate = c(0.15, 0.35), 
                       report.rate = c(0.85, 0.8))

state_init_N <- apply(init_N, c(1,2), FUN = function(x) sum(x))


# Generate priors for natural survival (NS), harvest survival (HS), fawn survival
# (S0), productivity (M), and reporting rate (report). Values reported as fawn to
# adult females then fawn to adult males on logit scale. 
ns_shape <- apply(array(t(data.frame(means = c(2.44, 2.44, 2.09, 2.59, 2.59, 1.5), 
                                     sds = c(0.25, 0.15, 0.15, 0.25, 0.2, 0.15))), 
                        dim = c(2, 3, 2)), c(1, 3), t)

hs_shape <- apply(array(t(data.frame(means = c(2.94, 1.9, 1.8, 2.94, 1, 0.75), 
                                     sds = c(0.35, 0.35, 0.35, 0.35, 0.35, 0.3))), 
                        dim = c(2, 3, 2)), c(1, 3), t)

rec_shape <- data.frame(means = c(-2.4, -0.35, -0.30), 
                        sds = 0.15)

report_shape <- data.frame(means = c(1.39, 1.15), sds = c(0.25, 0.25)) 

sd_years_HS <- c(shape = 3.92, rate = 28)

sd_years_NS <- c(shape = 0.98, rate = 14)

sd_years_report <- c(shape = 2, rate = 20)

sd_dmus_HS <- c(shape = 3.92, rate = 28)

sd_dmus_NS <- c(shape = 0.98, rate = 14)

sd_dmus_report <- c(shape = 2, rate = 20)

sd_years_S0 <- sd_dmus_S0 <- matrix(c(5.12, 32), ncol = 2)

# Set up some data for JAGS
O_dmu <- as.data.frame(apply(y_Af_DMU, c(2, 3), sum)) %>% 
  rownames_to_column(var = "year") %>% 
  pivot_longer(cols = contains("Unit"), names_to = "DMU", values_to = "fN") %>% 
  full_join(as.data.frame(apply(y_Am_DMU, c(2, 3), sum)) %>% 
              rownames_to_column(var = "year") %>% 
              pivot_longer(cols = contains("Unit"), names_to = "DMU", 
                           values_to = "mN")) %>% 
  split(.$DMU) %>% 
  map(~ .x %>% select(- year, - DMU)) %>% 
  abind::abind(., along = 3)

# Periods with different reporting requirements from 2005 to 2023
report.periods <- c(rep(1, 4), rep(2, 7), rep(3, 4), rep(4,4))


# Jags data 
jags_data_statewide <- list("nyears" = dim(y_H_DMU)[1], "nDMUs" = dim(y_H_DMU)[3],
                            "nages" = nrow(y_Am_statewide) + 1, 
                            "Period" = report.periods,
                            "initN" = init_N, "ns.shape" = ns_shape, 
                            "nperiods" = length(unique(report.periods)),
                            "hs.shape" = hs_shape,  
                            "rec.shape" = rec_shape,
                            "report.shape" = report_shape, 
                            "sd.yr.HS" = sd_years_HS, "sd.dmu.HS" = sd_dmus_HS, 
                            "sd.yr.NS" = sd_years_NS, "sd.dmu.NS" = sd_dmus_NS, 
                            "sd.yr.report" = sd_years_report, 
                            "sd.dmu.report" = sd_dmus_report, 
                            "sd.yr.rec" = sd_years_S0, "sd.dmu.rec" = sd_dmus_S0, 
                            "y.H" = y_H_statewide, "y.Af" = y_Af_statewide, 
                            "y.Am" = y_Am_statewide, "y.H.dmu" = y_H_DMU, 
                            "y.Af.dmu" = y_Af_DMU, "y.Am.dmu" = y_Am_DMU, 
                            "y.Hest" = h.est, "y.Hest.dmu" = h.est.dmu,
                            "O.aged" = cbind(apply(y_Af_statewide, 2, sum), 
                                             apply(y_Am_statewide, 2, sum)), 
                            "O.aged.dmu" = O_dmu, "npredyrs" = 0)


# Initial values to facilitate convergence statewide and by DMU
inits_statewide <- function() {list(N = init_pop_values(y.H = y_H_DMU, 
                                                        y.Af = y_Af_DMU, 
                                                        y.Am = y_Am_DMU))}


# Parameters to track
params <- c("N", "NS", "HS", "Report", "Total", "mu.NS", "mu.rec",
            "mu.HS", "mu.Report", "rec.rate", "lambda", "mlam",
            "Total.dmu", "rec.rate.dmu", "lambda.dmu", "mlam.dmu")



# Jags requirements
nt <- 5
nc <- 3
ni <- 400000
nb <- 200000
max.iter <- 1200000


# Run model---------------------------------------------------------------------

# Run IPM model
TN_Rec_statewide <- autojags(data = jags_data_statewide, inits = inits_statewide, 
                             parameters.to.save = params,
                             model.file = "TN_IPM_RecRate.txt", 
                             n.chains = nc, n.thin = nt, n.burnin = nb, 
                             iter.increment = ni, parallel = TRUE, 
                             max.iter = max.iter, 
                             codaOnly = c("lambda.dmu", "mlam.dmu"))





