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
#           4. Run the IPM model file (TN_IPM_Base.txt)
#   
#   Allison C Keever
#   TNTech University
#   Last Updated: 5/26/23

###############################################################################L

# Prep work space---------------------------------------------------------------

# Load packages needed in analyses
library(tidyverse)
library(R2jags)
library(tidybayes)
library(abind)

# Run script to load data and create data frames for analyses
source("Rscripts/TN_harvest_data_prep.R")

# Source functions needed for analyses
source("Rscripts/Functions_AbundanceEstimation.R")


# Bayesian requirements---------------------------------------------------------

# Estimate initial values for population size based on reported harvest, 
# reporting rates, and harvest rates. 
init_N <- est_init_pop(y.H = y_H_DMU, y.Af = y_Af_DMU, 
                       y.Am = y_Am_DMU, harvest.rate = c(0.12, 0.25), 
                       report.rate = c(0.8, 0.76))

state_init_N <- apply(init_N, c(1,2), FUN = function(x) sum(x))


# Generate priors for natural survival (NS), harvest survival (HS), fawn survival
# (S0), productivity (M), and reporting rate (report). Values reported as fawn to
# adult females then fawn to adult males on logit scale. 
ns_shape <- apply(array(t(data.frame(means = c(2.44, 2.44, 2.09, 2.59, 2.59, 1.32), 
                                   sds = c(0.3, 0.2, 0.2, 0.3, 0.2, 0.2))), 
                        dim = c(2, 3, 2)), c(1, 3), t)

hs_shape <- apply(array(t(data.frame(means = c(2.94, 1.9, 1.8, 2.94, 1, 0.75), 
                                     sds = c(0.35, 0.4, 0.4, 0.35, 0.35, 0.3))), 
                        dim = c(2, 3, 2)), c(1, 3), t)

s0_shape <- data.frame(means = c(-0.66, -0.66), sds = rep(0.35, 2))

m_shape <- generate_priors(values = data.frame(means = c(0.24, 1.42, 1.77),
                                               vars = rep(0.005, 3)), 
                           dist.type = "gamma")

report_shape <- apply(array(t(data.frame(means = rep(c(rep(1.39, 3), rep(1.15, 3)), 3), 
                                         sds = c(rep(0.3, 18)))), 
                            dim = c(2, 3, 2, 3)), c(1, 3, 4), t)

sd_years_HS <- apply(array(t(data.frame(shape = c(0.5, 1.28, 1.62, 0.5, 8.82, 3.92), 
                                     rate = c(10, 16, 18, 10, 42, 28))), 
                        dim = c(2, 3, 2)), c(1, 3), t)

sd_years_NS <- apply(array(t(data.frame(shape = c(2.88, 0.02, 0.08, 2.88, 0.18, 0.98), 
                                        rate = c(24, 2, 4, 24, 6, 14))), 
                           dim = c(2, 3, 2)), c(1, 3), t)

sd_years_report <- apply(array(t(data.frame(shape = c(rep(2, 6)), 
                                        rate = c(rep(20, 6)))), 
                           dim = c(2, 3, 2)), c(1, 3), t)

sd_dmus_HS <- apply(array(t(data.frame(shape = c(0.5, 1.28, 1.62, 0.5, 8.82, 3.92), 
                                       rate = c(10, 16, 18, 10, 42, 28))), 
                          dim = c(2, 3, 2)), c(1, 3), t)

sd_dmus_NS <- apply(array(t(data.frame(shape = c(2.88, 0.02, 0.08, 2.88, 0.18, 0.98), 
                                       rate = c(24, 2, 4, 24, 6, 14))), 
                          dim = c(2, 3, 2)), c(1, 3), t)

sd_dmus_report <- apply(array(t(data.frame(shape = c(rep(2, 6)), 
                                           rate = c(rep(20, 6)))), 
                              dim = c(2, 3, 2)), c(1, 3), t)

sd_years_S0 <- sd_dmus_S0 <- matrix(c(5.12, 32), ncol = 2)


# Set up some data for JAGS
O_dmu <- as.data.frame(apply(y_Af_DMU, c(2, 3), sum)) %>% 
  rownames_to_column(var = "year") %>% 
  pivot_longer(cols = 2:10, names_to = "DMU", values_to = "fN") %>% 
  full_join(as.data.frame(apply(y_Am_DMU, c(2, 3), sum)) %>% 
              rownames_to_column(var = "year") %>% 
              pivot_longer(cols = 2:10, names_to = "DMU", values_to = "mN")) %>% 
  split(.$DMU) %>% 
  map(~ .x %>% select(- year, - DMU)) %>% 
  abind::abind(., along = 3)

# Periods with different reporting requirements from 2005 to 2019
report.periods <- c(rep(1, 4), rep(2, 7), rep(3, 4))


# Jags data 
jags_data_statewide <- list("nyears" = dim(y_H_DMU)[1], "nDMUs" = dim(y_H_DMU)[3],
                            "nages" = nrow(y_Am_statewide) + 1, 
                            "Period" = report.periods,
                            "initN" = init_N, "ns.shape" = ns_shape, 
                            "nperiods" = length(unique(report.periods)),
                            "hs.shape" = hs_shape, "s0.shape" = s0_shape, 
                            "m.shape" = m_shape, "report.shape" = report_shape, 
                            "sd.yr.HS" = sd_years_HS, "sd.dmu.HS" = sd_dmus_HS, 
                            "sd.yr.NS" = sd_years_NS, "sd.dmu.NS" = sd_dmus_NS, 
                            "sd.yr.report" = sd_years_report, "sd.dmu.report" = sd_dmus_report, 
                            "sd.yr.S0" = sd_years_S0, "sd.dmu.S0" = sd_dmus_S0, 
                            "y.H" = y_H_statewide, "y.Af" = y_Af_statewide, 
                            "y.Am" = y_Am_statewide, "y.H.dmu" = y_H_DMU, 
                            "y.Af.dmu" = y_Af_DMU, "y.Am.dmu" = y_Am_DMU, 
                            "y.Hest" = h.est, "y.Hest.dmu" = h.est.dmu,
                            "O.aged" = cbind(apply(y_Af_statewide, 2, sum), 
                                             apply(y_Am_statewide, 2, sum)), 
                            "O.aged.dmu" = O_dmu, "npredyrs" = 0)


# Initial values to facilitate convergence statewide and by DMU
inits_statewide <- function() {list(N = init_pop_values(y.H = y_H_DMU, 
                                                        y.Af = y_Af_DMU[-4,,], 
                                                        y.Am = abind(y_Af_DMU[4,,], 
                                                                     y_Am_DMU, 
                                                                     along = 1)))}


# Parameters to track
params <- c("N", "NS", "HS", "S0", "Report", "Hest", "Total", "mu.NS",
            "mu.HS", "mu.Report", "mu.S0", "m", "sd.NS.yr", "sd.HS.yr",
            "sd.Report.yr", "sd.NS.dmu", "sd.HS.dmu", "sd.Report.dmu", "sd.S0.yr",
            "sd.S0.dmu",  "fawns", "does", "bucks", "H.total",
            "H.antlerless", "H.does", "H.yr.bucks", "H.ad.bucks", "rec.rate",
            "lambda", "mlam", "Total.dmu", "fawns.dmu", "does.dmu", "bucks.dmu",
            "rec.rate.dmu", "H.total.dmu", "H.antlerless.dmu", "H.yr.bucks.dmu",
            "H.ad.bucks.dmu", "lambda.dmu", "mlam.dmu")



# Jags requirements
nt <- 5
nc <- 3
ni <- 1000000
nb <- 100000


# Run model---------------------------------------------------------------------

# Run IPM model
TN_statewide <- run_model("TN_IPM_Base.txt", ni = ni, nb = nb, nt = nt, 
                          nc = nc, win.data = jags_data_statewide, 
                          params = params, inits = inits_statewide)





