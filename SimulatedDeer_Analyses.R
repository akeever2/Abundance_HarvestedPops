###############################################################################L

#   Simulated deer analyses
#       
#       Running this script will simulate deer populations and monitoring data
#       from those populations. Then those data will be used in an IPM to 
#       estimate deer abundance. 
#   
#   Allison C Keever
#   TNTech University
#   Last Updated: 6/5/23

###############################################################################L

# Prep work space----------------------------------------------------------------

# Load packages required for analyses
library(abind)
library(R2jags)
library(tidyverse)
library(tidybayes)


# Source code or data needed for analyses 
source("Rscripts/SimPopFunctions.R")


# Run function to perform analyses----------------------------------------------

# Run the wrapper function to simulate deer populations and harvest data, then 
# run IPM to estimate abundance of simulated deer. This is for a single scenario.

# props <- c(0.02, 0.06, 0.10)
props <- 0.02

for(i in seq_along(props)){
  
  sim_res <- sim_deer_analysis(init.N = rand_deer_init_pop(ndmus = 9, mu.init = 70000), 
                              mu.NS = matrix(data = c(0.92, 0.92, 0.88, 0.93, 0.93, 0.79), nrow = 3), 
                              sd.NS.dmu = matrix(rgamma(6, shape = 5.78, rate = 34), nrow = 3), 
                              sd.NS.yr = matrix(rgamma(6, shape = 5.78, rate = 34), nrow = 3), 
                              mu.HS = matrix(data = c(0.95, 0.87, 0.85, 0.95, 0.73, 0.67), nrow = 3), 
                              sd.HS.dmu = matrix(rgamma(6, shape = 5.78, rate = 34),nrow = 3), 
                              sd.HS.yr = matrix(rgamma(6, shape = 5.78, rate = 34), nrow = 3),
                              mu.report = matrix(data = c(0.80, 0.80, 0.80, 0.76, 0.76, 0.76), nrow = 3), 
                              sd.report.dmu = matrix(rgamma(6, shape = 5.78, rate = 34), nrow = 3), 
                              sd.report.yr = matrix(rgamma(6, shape = 5.78, rate = 34), nrow = 3), 
                              correct.sex = matrix(data = c(1, 1, 1, 1, 1, 1),  nrow = 3), 
                              mu.S0 = 0.34, sd.S0.dmu = 0.2, sd.S0.yr = 0.34, 
                              m = c(0.24, 1.42, 1.77), nyears = 10, nsims = 100, 
                              ndmus = 9, prop.f = props[i], prop.m = props[i], ni = 100000, 
                              nb = 25000, model.scripts = "IPM_HarvestRecon_Base.txt")
  
  save(sim_res, file = file.path(paste0("Results/sim02_", j, ".RData")))
  rm(sim_res)
  
}



