##############################################################################L

#   Functions for simulated deer monitoring analyses
#   Allison C Keever
#   TNTech University
#   Last updated: 6/5/23

##############################################################################L

# Function to simulate deer populations with harvest
deerpop_sim <- function(init.N, mu.NS, sd.NS.dmu, sd.NS.yr, mu.HS, sd.HS.dmu, 
                         sd.HS.yr, mu.S0, sd.S0.dmu, sd.S0.yr, m, nyears, 
                         nsims, ndmus){
  
  # Initial set up ----------------------------------
  
  # cat("Simulating deer populations: \n")
  
  # Progress bar set up
  pb <- txtProgressBar(min = 0, max = nsims, style = 3, width = nsims, char = "-")
  init_pb <- numeric(nsims)
  end_pb <- numeric(nsims)
  
  # Set up arrays to hold all the information
  N <- NS <- HS <- Hest <- array(NA, dim = c(3, nyears, 2, ndmus, nsims))
  eps.NS.dmu <- eps.HS.dmu <- array(NA, dim = c(3, 2, ndmus, nsims))
  eps.NS.yr <- eps.HS.yr <- array(NA, dim = c(3, nyears, 2, nsims))
  eps.S0.dmu <- array(NA, dim = c(ndmus, nsims))
  eps.S0.yr <- array(NA, dim = c(nyears, nsims))
  Fer <- array(NA, dim = c(3,nyears, ndmus, nsims))
  S0 <- array(NA, dim = c(nyears, 2, ndmus, nsims))
  lambda.dmu <- logla.dmu <- array(NA, dim = c(nyears, ndmus, nsims))
  lambda <- logla <- array(NA, dim = c(nyears, nsims))
  mlam.dmu <- array(NA, dim = c(ndmus, nsims))
  mlam <- array(NA, dim = c(nsims))
  
  # Throw warnings if arguments aren't set up properly 
  if(!all(dim(mu.NS) == c(3, 2))){ 
    stop("mu.NS does not have proper dimensions [age, sex]")
  } 
  
  if(!all(dim(mu.HS) == c(3, 2))){ 
    stop("mu.HS does not have proper dimensions [age, sex]")
  } 
  
  if(!all(dim(init.N) == c(3, 2, ndmus))){ 
    stop("init.N does not have proper dimensions [age, sex, DMU]")
  } 

  if(!all(dim(sd.NS.dmu) == c(3, 2)) | !all(dim(sd.NS.yr) == c(3, 2))){ 
    stop("sd.NS.dmu OR sd.NS.yr does not have proper dimensions [age, sex]")
  } 

  if(!all(dim(sd.HS.dmu) == c(3, 2)) | !all(dim(sd.HS.yr) == c(3, 2))){ 
    stop("sd.HS.dmu OR sd.HS.yr does not have proper dimensions [age, sex]")
  } 
  
  if(length(m) != 3){ 
    stop("m does not have proper dimensions [age]")
  }


  # Code to simulate populations  ----------------------------------
  
  
  # Loop through each simulated population. Each nsim pop will run for nyears. 
  for(sim in 1:nsims){
    
    # Time at start of each simulation for progress bar
    init_pb[sim] <- Sys.time()
    
    # Draw random effects for year and site
    for(a in 1:3){
      for(s in 1:2){
        eps.NS.dmu[a,s,,sim] <- rnorm(n = ndmus, 0, sd.NS.dmu[a,s])
        eps.HS.dmu[a,s,,sim] <- rnorm(n = ndmus, 0, sd.HS.dmu[a,s])
        eps.NS.yr[a,,s,sim] <- rnorm(n = nyears, 0, sd.NS.yr[a,s])
        eps.HS.yr[a,,s,sim] <- rnorm(n = nyears, 0, sd.HS.yr[a,s]) 
      }
    }
    eps.S0.dmu[,sim] <- rnorm(n = ndmus, 0, sd.S0.dmu)
    eps.S0.yr[,sim] <- rnorm(n = nyears, 0, sd.S0.yr)
    
    # Calculate demographic rates
    for(yr in 1:nyears){
      for(i in 1:ndmus){
        for(a in 1:3){
          for(s in 1:2){
            NS[a,yr,s,i,sim] <- plogis(qlogis(mu.NS[a,s]) + eps.NS.dmu[a,s,i,sim] +
                                         eps.NS.yr[a,yr,s,sim])
            HS[a,yr,s,i,sim] <- plogis(qlogis(mu.HS[a,s]) + eps.HS.dmu[a,s,i,sim] +
                                         eps.HS.yr[a,yr,s,sim])
          } # s
          
          Fer[a,yr,i,sim] <- (NS[a,yr,1,i,sim] ^ (8 / 12)) * HS[a,yr,1,i,sim] * m[a]
        } # a
        
        for(s in 1:2){
          S0[yr,s,i,sim] <- plogis(qlogis(mu.S0) + eps.S0.dmu[i,sim] + 
                                     eps.S0.yr[yr,sim])
        }
      }
    }
    
    # Set pop size in year 1
    for(a in 1:3){
      for(s in 1:2){
        for(i in 1:ndmus){
          N[a,1,s,i,sim] <- round(init.N[a,s,i])
          Hest[a,1,s,i,sim] <- rbinom(1, N[a,1,s,i,sim], 1 - HS[a,1,s,i,sim])
        } # i
      } # s
    } # a
      
      
    # Loop through each year for each simulated population
    for(yr in 2:nyears){
      # Calculate number harvested from each age and sex class
      for(i in 1:ndmus){
        for(s in 1:2){
          # Project the population for each age and sex class
          N[1,yr,s,i,sim] <- rpois(1, N[1,yr-1,1,i,sim] * Fer[1,yr-1,i,sim] * 0.5 * S0[yr-1,s,i,sim] +
                                     N[2,yr-1,1,i,sim] * Fer[2,yr-1,i,sim] * 0.5 * S0[yr-1,s,i,sim] +
                                     N[3,yr-1,1,i,sim] * Fer[3,yr-1,i,sim] * 0.5 * S0[yr-1,s,i,sim])
          
          N[2,yr,s,i,sim] <- rbinom(1, N[1,yr-1,s,i,sim] - Hest[1,yr-1,s,i,sim], NS[1,yr-1,s,i,sim])
          
          N[3,yr,s,i,sim] <- rbinom(1, N[2,yr-1,s,i,sim] - Hest[2,yr-1,s,i,sim], NS[2,yr-1,s,i,sim]) + 
            rbinom(1, N[3,yr-1,s,i,sim] - Hest[3,yr-1,s,i,sim], NS[3,yr-1,s,i,sim])
          
          for(a in 1:3){
            Hest[a,yr,s,i,sim] <- rbinom(1, N[a,yr,s,i,sim], 1 - HS[a,yr,s,i,sim])
          }
        }
        lambda.dmu[yr - 1,i,sim] <- sum(N[,yr,,i,sim]) / sum(N[,yr - 1,,i,sim])
        logla.dmu[yr - 1,i,sim] <- log(lambda.dmu[yr-1,i,sim])
      }
      
      lambda[yr - 1,sim] <- sum(N[,yr,,,sim]) / sum(N[,yr - 1,,,sim])
      logla[yr - 1,sim] <- log(lambda[yr-1,sim])
      
    } # end year loop
    
    for(i in 1:ndmus){
      mlam.dmu[i,sim] <- exp((1 / (nyears - 1)) * sum(logla.dmu[1:(nyears - 1),i,sim]))
    }
    
    mlam[sim] <- exp((1 / (nyears - 1)) * sum(logla[1:(nyears - 1),sim]))
    
  # Wrap up the function ----------------------------------------
    
    
    # End time for each simulation for progress bar
    end_pb[sim] <- Sys.time()
    
    # Progress bar update
    setTxtProgressBar(pb, sim)
    
    # Estimate remaining time
    time_pb <- round(lubridate::seconds_to_period(sum(end_pb - init_pb)), 0)
    est <- nsims * (mean(end_pb[end_pb != 0] - init_pb[init_pb != 0])) - time_pb
    rem_pb <- round(lubridate::seconds_to_period(est), 0)
    
    # Paste remaining time
    # cat(paste("  Estimated time remaining:", rem_pb, sep = ""), "")
    
    cat("\r  Estimated time remaining:", rem_pb, "secs")
    
    } # end simulation loop
  
  close(pb)
  
  return(list(N = N, Hest = Hest, S0 = S0, NS = NS, HS = HS, mu.NS = mu.NS, 
              sd.NS.dmu = sd.NS.dmu, sd.NS.yr = sd.NS.yr, mu.HS = mu.HS, 
              sd.HS.dmu = sd.HS.dmu, sd.HS.yr = sd.HS.yr, mu.S0 = mu.S0, 
              sd.S0.dmu = sd.S0.dmu, sd.S0.yr = sd.S0.yr, m = m, nyears = nyears, 
              nsims = nsims, ndmus = ndmus, lambda = lambda, logla = logla, 
              mlam = mlam, lambda.dmu = lambda.dmu, logla.dmu = logla.dmu, 
              mlam.dmu = mlam.dmu))
  
  
}


# Function to introduce error in the harvest data
harvest_data_scramble <- function(sim.results, mu.report, sd.report.dmu, 
                                  sd.report.yr, correct.sex, prop.f, prop.m){
  
  # Initial set up ----------------------------------
  
  
  # Set values from sim.results to make things easier
  nyears <- sim.results$nyears
  ndmus <- sim.results$ndmus
  nsims <- sim.results$nsims
  
  # Progress bar set up
  pb <- txtProgressBar(min = 0, max = nsims, style = 3, width = nsims, char = "-")
  init_pb <- numeric(nsims)
  end_pb <- numeric(nsims)
  
  # Set up arrays to hold all the information
  Report <- Hrep <- Hcorr <- array(NA, dim = c(3, nyears, 2, ndmus, nsims))
  eps.report.dmu <- array(NA, dim = c(3, 2, ndmus, nsims))
  eps.report.yr <- array(NA, dim = c(3, nyears, 2, nsims))
  y.H <- array(NA, dim = c(nyears, 2, nsims))
  y.H.dmu <- array(NA, dim = c(nyears, 2, ndmus, nsims))
  y.Af.dmu <- y.Am.dmu <- array(NA, dim = c(3, nyears, ndmus, nsims))
  y.Af <- y.Am <- array(NA, dim = c(3, nyears, nsims))
  
  # Throw warnings if arguments aren't set up properly 
  if(!all(dim(mu.report) == c(3, 2))){ 
    stop("mu.report does not have proper dimensions [age, sex]")
  } 
  
  if(!all(dim(sd.report.dmu) == c(3, 2)) | !all(dim(sd.report.yr) == c(3, 2))){ 
    stop("sd.report.dmu OR sd.report.yr does not have proper dimensions [age, sex]")
  } 
  
  
  # Code to simulate harvest data  ----------------------------------
  
  
  # Loop through each simulated population. Each nsim pop will run for nyears. 
  for(sim in 1:nsims){
    
    # Time at start of each simulation for progress bar
    init_pb[sim] <- Sys.time()
    
    # Draw random effects for year and site
    for(a in 1:3){
      for(s in 1:2){
        eps.report.dmu[a,s,,sim] <- rnorm(n = ndmus, 0, sd.report.dmu[a,s])
        eps.report.yr[a,,s,sim] <- rnorm(n = nyears, 0, sd.report.yr[a,s]) 
      }
    }
    
    # Calculate demographic rates
    for(yr in 1:nyears){
      for(i in 1:ndmus){
        for(a in 1:3){
          for(s in 1:2){
            Report[a,yr,s,i,sim] <- plogis(qlogis(mu.report[a,s]) + eps.report.dmu[a,s,i,sim] +
                                         eps.report.yr[a,yr,s,sim])
          } # s
        } # a
      } # i
    } # yr
    
    # Get reported harvest, including reporting and sexing errors
    for(yr in 1:nyears){
      for(i in 1:ndmus){
        for(a in 1:3){
          for(s in 1:2){
            # Reporting rate 
            Hrep[a,yr,s,i,sim] <- rbinom(1, sim.results$Hest[a,yr,s,i,sim], 
                                         Report[a,yr,s,i,sim])
          } # s
          
          # Sexing errors
          correct_f <- rbinom(1, Hrep[a,yr,1,i,sim], correct.sex[a,1])
          correct_m <- rbinom(1, Hrep[a,yr,2,i,sim], correct.sex[a,2])
          wrong_f <- Hrep[a,yr,1,i,sim] - correct_f
          wrong_m <- Hrep[a,yr,2,i,sim] - correct_m
          Hcorr[a,yr,1,i,sim] <- correct_f + wrong_m
          Hcorr[a,yr,2,i,sim] <- correct_m + wrong_f
        }
      }
    }
    
    # Format data for IPM models
    for(yr in 1:nyears){
      # Statewide bag harvest data
      y.H[yr,1,sim] <- sum(Hcorr[,yr,1,,sim]) + sum(Hcorr[1,yr,2,,sim])
      y.H[yr,2,sim] <- sum(Hcorr[2:3,yr,2,,sim])
      
      # DMU bag harvest data
      for(i in 1:ndmus){
        y.H.dmu[yr,1,i,sim] <- sum(Hcorr[,yr,1,i,sim]) + Hcorr[1,yr,2,i,sim]
        y.H.dmu[yr,2,i,sim] <- sum(Hcorr[2:3,yr,2,i,sim])
      }
      
      # DMU aged/sexed harvest data
      for(a in 1:3){
        for(i in 1:ndmus){
          y.Af.dmu[a,yr,i,sim] <- round(Hcorr[a,yr,1,i,sim] * prop.f)
          y.Am.dmu[a,yr,i,sim] <- round(Hcorr[a,yr,2,i,sim] * prop.m)
        }
        y.Af[a,yr,sim] <- sum(y.Af.dmu[a,yr,,sim])
        y.Am[a,yr,sim] <- sum(y.Am.dmu[a,yr,,sim])
      }
      
      
    }
    
    
    
    # Wrap up the function ----------------------------------------
    
    
    # End time for each simulation for progress bar
    end_pb[sim] <- Sys.time()
    
    # Progress bar update
    setTxtProgressBar(pb, sim)
    
    # Estimate remaining time
    time_pb <- round(lubridate::seconds_to_period(sum(end_pb - init_pb)), 0)
    est <- nsims * (mean(end_pb[end_pb != 0] - init_pb[init_pb != 0])) - time_pb
    rem_pb <- round(lubridate::seconds_to_period(est), 0)
    
    # Paste remaining time
    cat("\r  Estimated time remaining:", rem_pb, "secs")
    
  } # end simulation loop
  
  close(pb)
  
  sim.results$Report <- Report
  sim.results$Hrep <- Hcorr
  sim.results$mu.report <- mu.report
  sim.results$sd.report.dmu <- sd.report.dmu
  sim.results$sd.report.yr <- sd.report.yr
  sim.results$y.H <- y.H 
  sim.results$y.H.dmu <- y.H.dmu
  sim.results$y.Af.dmu <- y.Af.dmu
  sim.results$y.Am.dmu <- y.Am.dmu
  sim.results$y.Af <- y.Af
  sim.results$y.Am <- y.Am
  
  return(sim.results)
  
}


# Function to run the IPM models
run_IPM <- function(sim.results, mu.NS, mu.HS, mu.report, mu.S0, m, ni, nb, 
                    model.scripts){

  
  # Initial set up ----------------------------------

  # Set values from sim.results to make things easier
  nyears <- sim.results$nyears
  ndmus <- sim.results$ndmus
  nsims <- sim.results$nsims
  
  # Progress bar set up
  pb <- txtProgressBar(min = 0, max = nsims, style = 3, width = nsims, char = "-")
  init_pb <- numeric(nsims)
  end_pb <- numeric(nsims)
  
  # Generate priors for natural survival (NS), harvest survival (HS), fawn survival
  # (S0), productivity (M), and reporting rate (report). Values reported as fawn to
  # adult females then fawn to adult males on logit scale. 
  ns_shape <- apply(array(t(data.frame(means = qlogis(as.vector(mu.NS)), 
                                       sds = c(0.3, 0.25, 0.25, 0.3, 0.25, 0.25))), 
                          dim = c(2, 3, 2)), c(1, 3), t)
  
  hs_shape <- apply(array(t(data.frame(means = qlogis(as.vector(mu.HS)), 
                                       sds = c(0.35, 0.4, 0.4, 0.35, 0.35, 0.35))), 
                          dim = c(2, 3, 2)), c(1, 3), t)
  
  s0_shape <- data.frame(means = rep(qlogis(mu.S0), 2), sds = rep(0.35, 2))
  
  capture.output(m_shape <- generate_priors(values = data.frame(means = m,
                                                 vars = rep(0.005, 3)), 
                             dist.type = "gamma"))
  
  report_shape <- apply(array(t(data.frame(means = qlogis(as.vector(mu.report)), 
                                           sds = c(rep(0.5, 6)))), 
                              dim = c(2, 3, 2)), c(1, 3), t)
  
  capture.output(sd_years <- generate_priors(values = data.frame(means = rep(0.17, 4),
                                                  vars = rep(0.005, 4)), 
                              dist.type = "gamma"))
  
  capture.output(sd_dmus <- generate_priors(values = data.frame(means = rep(0.17, 4),
                                                 vars = rep(0.005, 4)), 
                             dist.type = "gamma"))
  
  # Set up lists to hold the information and set constant JAGS requirements
  out <- list()
  params <- c("N", "NS", "HS", "S0", "Report", "Hest", "Hrep", "Total", "mu.NS",
              "mu.HS", "mu.Report", "mu.S0", "m", "sd.NS.yr", "sd.HS.yr", 
              "sd.Report.yr", "sd.NS.dmu", "sd.HS.dmu", "sd.Report.dmu", "sd.S0.yr", 
              "sd.S0.dmu", "fawns", "does", "bucks", "H.total", 
              "H.antlerless", "H.yr.bucks", "H.ad.bucks", "rec.rate", 
              "lambda", "mlam", "Total.dmu", "fawns.dmu", "does.dmu", "bucks.dmu",
              "rec.rate.dmu", "H.total.dmu", "H.antlerless.dmu", "H.yr.bucks.dmu", 
              "H.ad.bucks.dmu", "lambda.dmu", "mlam.dmu")
  nt <- 5
  nc <- 3  
  
  
  
  # Code to run the models  ----------------------------------
  
  for(i in 1:nsims){

    # Time at start of each simulation for progress bar
    init_pb[i] <- Sys.time()
    
    # Estimate initial values for population size based on reported harvest, 
    # reporting rates, and harvest rates. 
    suppressMessages(init_N <- est_init_pop(y.H = sim.results$y.H.dmu[,,,i], 
                           y.Af = sim.results$y.Af.dmu[,,,i], 
                           y.Am = sim.results$y.Am.dmu[,,,i], 
                           harvest.rate = c(0.09, 0.20), 
                           report.rate = c(0.8, 0.76)))
    
    init_vals <- sim.results$N[,,,,i]
    init_vals[3,,,] <- NA
    init_vals[,1,,] <- NA
    
    # Set up some data for JAGS
    suppressMessages(O_dmu <- as.data.frame(apply(sim.results$y.Af.dmu[,,,i], c(2, 3), sum)) %>% 
      rownames_to_column(var = "year") %>% 
      pivot_longer(cols = 2:10, names_to = "DMU", values_to = "fN") %>% 
      full_join(as.data.frame(apply(sim.results$y.Am.dmu[,,,i], c(2, 3), sum)) %>% 
                  rownames_to_column(var = "year") %>% 
                  pivot_longer(cols = 2:10, names_to = "DMU", values_to = "mN")) %>% 
      split(.$DMU) %>% 
      map(~ .x %>% select(- year, - DMU)) %>% 
      abind::abind(., along = 3))

    
    jags.data <- list("nyears" = nyears, "nages" = 3, "nDMUs" = ndmus, 
                      "initN" = init_N, "ns.shape" = ns_shape, 
                      "hs.shape" = hs_shape, "s0.shape" = s0_shape, 
                      "m.shape" = m_shape, "report.shape" = report_shape, 
                      "sd.yr.gam" = sd_years, "sd.dmu.gam" = sd_dmus, 
                      "y.H" = sim.results$y.H[,,i], "y.Af" = sim.results$y.Af[,,i], 
                      "y.Am" = sim.results$y.Am[,,i], 
                      "y.H.dmu" = sim.results$y.H.dmu[,,,i], 
                      "y.Af.dmu" = sim.results$y.Af.dmu[,,,i], 
                      "y.Am.dmu" = sim.results$y.Am.dmu[,,,i],
                      "O.aged" = cbind(apply(sim.results$y.Af[,,i], 2, sum),
                                       apply(sim.results$y.Am[,,i], 2, sum)), 
                      "O.aged.dmu" = O_dmu, "npredyrs" = 0)
    
  
  inits <- function() {list(N = init_vals)}
  

  
  # Run the IPM model 
 
    base <- run_model(model.scripts, ni = ni, nb = nb, nt = nt, nc = nc, 
                      win.data = jags.data, params = params, inits = inits)
    
    out[[i]] <- base
    
    # Wrap up the function ----------------------------------------
    
    
    # End time for each simulation for progress bar
    end_pb[i] <- Sys.time()
    
    # Progress bar update
    setTxtProgressBar(pb, i)
    
    # Estimate remaining time
    time_pb <- round(lubridate::seconds_to_period(sum(end_pb - init_pb)), 0)
    est <- nsims * (mean(end_pb[end_pb != 0] - init_pb[init_pb != 0])) - time_pb
    rem_pb <- round(lubridate::seconds_to_period(est), 0)
    
    # Paste remaining time
    cat("\r  Sim #", i, "; Estimated time remaining:", rem_pb, "secs;")
    
  } # End simulations
  
  close(pb)
  
  cat("\n", "", "\n Estimated time per model run:", (mean(end_pb[end_pb != 0] - init_pb[init_pb != 0])))
  
  return(out)
  
}

# Function to calculate performance metrics
calc_pm <- function(sim.results, abund.est){
  
  # Initial set up ----------------------------------
  
  
  # Set values from sim.results to make things easier
  nyears <- sim.results$nyears
  ndmus <- sim.results$ndmus
  nsims <- sim.results$nsims
  
  
  # Code to calculate performance metrics  ----------------------------------
  
  # First set up a results data frame to hold everything, then add the estimates
  results <- list()
  for(i in 1:nsims){
    out.mcmc <- as.mcmc(abund.est[[i]])
    results[[i]] <- out.mcmc %>% 
      gather_draws(N[age,yr,sex,dmu], NS[age,yr,sex,dmu], HS[age,yr,sex,dmu], 
                   Report[age,yr,sex,dmu], Total[yr], H.total[yr], rec.rate[yr],
                   lambda[yr], mlam, Total.dmu[dmu,yr], H.total.dmu[dmu,yr], 
                   rec.rate.dmu[dmu,yr], lambda.dmu[dmu,yr], mlam.dmu[dmu]) %>%  
      mean_qi(.width = c(0.90)) %>% select(-.width, -.point, -.interval) %>%
      mutate(sim = i)
  }
  
  results <- bind_rows(results)
  
  # To simplify matters, get truth from simulations organized
  truth <- list()
  for(i in 1:nsims){
    N <- array_tree(sim.results$N[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr",
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "N")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2), sim = i)
    
    NS <- array_tree(sim.results$NS[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "NS")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2), sim = i)
    
    HS <- array_tree(sim.results$HS[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr",
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "HS")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2), sim = i)
    
    Report <- array_tree(sim.results$Report[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr",
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "Report")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2), sim = i)
    
    H.total <- array_tree(sim.results$Hest[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "Hest")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2)) %>% 
      group_by(yr) %>% 
      summarise(truth = sum(truth)) %>% 
      mutate(.variable = "H.total", sim = i)
    
    H.total.dmu <- array_tree(sim.results$Hest[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "Hest")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2)) %>% 
      group_by(yr, dmu) %>% 
      summarise(truth = sum(truth)) %>% 
      mutate(.variable = "H.total.dmu", sim = i)
    
    rec.rate <- array_tree(sim.results$N[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "N")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2)) %>% 
      group_by(age, sex, yr) %>% 
      summarise(tot = sum(truth)) %>% ungroup() %>% 
      group_by(yr) %>% 
      summarise(truth = sum(tot[age == 1]) / sum(tot[age == 2 | age == 3 & sex == 1])) %>%
      mutate(.variable = "rec.rate", sim = i)
    
    rec.rate.dmu <- array_tree(sim.results$N[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "N")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2)) %>% 
      group_by(age, sex, yr, dmu) %>% 
      summarise(tot = sum(truth)) %>% ungroup() %>% 
      group_by(yr, dmu) %>% 
      summarise(truth = sum(tot[age == 1]) / sum(tot[age == 2 | age == 3 & sex == 1])) %>%
      mutate(.variable = "rec.rate.dmu", sim = i)
    
    Total <- array_tree(sim.results$N[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "N")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2)) %>% 
      group_by(yr) %>% 
      summarise(truth = sum(truth)) %>% 
      mutate(.variable = "Total", sim = i)
    
    Total.dmu <- array_tree(sim.results$N[,,,,i], 4) %>%
      map(~ as.data.frame(do.call(rbind, array_tree(.x, 3))) %>%
            mutate(age = rep(1:3, 2), sex = rep(1:2, each = 3)) %>% 
            pivot_longer(cols = contains("V"), names_to = "yr", 
                         names_transform = list(yr = as.integer),
                         names_prefix = "V", values_to = "truth") %>%
            mutate(.variable = "N")) %>%
      bind_rows() %>% 
      mutate(dmu = rep(1:ndmus, each = nyears * 3 * 2)) %>% 
      group_by(yr, dmu) %>% 
      summarise(truth = sum(truth)) %>% 
      mutate(.variable = "Total.dmu", sim = i)
    
    lambda.dmu <- as.data.frame(sim.results$lambda.dmu[,,i]) %>% 
      slice_head(n = nyears - 1) %>%
      mutate(yr = 1:(nyears - 1), sim = i) %>%
      pivot_longer(cols = contains("V"), names_to = "dmu", 
                   names_transform = list(dmu = as.integer),
                   names_prefix = "V", values_to = "truth") %>%
      mutate(.variable = "lambda.dmu")
    
    lambda <- data.frame(truth = sim.results$lambda[,i]) %>% 
      slice_head(n = nyears - 1) %>%
      mutate(yr = 1:(nyears - 1), sim = i, .variable = "lambda")
    
    mlam <- data.frame(truth = sim.results$mlam[i]) %>% 
      mutate(sim = i, .variable = "mlam")
    
    mlam.dmu <- data.frame(truth = sim.results$mlam.dmu[,i]) %>% 
      mutate(dmu = 1:ndmus, sim = i, .variable = "mlam.dmu")
    
    
    truth[[i]] <- bind_rows(N, HS, NS, Report, H.total.dmu, H.total, rec.rate.dmu, 
                            rec.rate, Total.dmu, Total, lambda, lambda.dmu, 
                            mlam, mlam.dmu)
    
    rm(N, HS, NS, Report, H.total.dmu, H.total, rec.rate.dmu, rec.rate, Total.dmu, 
       Total, lambda, lambda.dmu, mlam, mlam.dmu)
    
  }
  
  truth <- bind_rows(truth)
  
  # Add in truth from the simulation results
  results <- results %>% left_join(truth)
  
  
  # Performance metrics: coefficient of error, bias, and coverage
  if(nsims > 1){
    pm_datum <- bind_rows(
    results %>%
      mutate(MSE = (.value - truth)^2, 
             bias = (.value - truth) / truth, 
             cover = ifelse(truth >= .lower & truth <= .upper, 1, 0)) %>% 
      filter(!(.variable %in% c("mlam", "mlam.dmu"))) %>%
      group_by(sim, age, sex, dmu, .variable) %>%
      summarise(MSEyr = sum(MSE) / (n() - 1), 
                biasyr = sum(bias), 
                coveryr = sum(cover), 
                totalN = sum(truth)) %>%
      ungroup() %>%
      group_by(age, sex, dmu, .variable) %>%
      summarise(coverage = sum(coveryr) / (nsims * nyears), 
                bias = sum(biasyr) / (nsims * nyears), 
                MSE = sum(MSEyr) / nsims, 
                totalN = sum(totalN)) %>%
      mutate(CE = sqrt(MSE) / (totalN / (nsims * nyears))) %>%
      select(age, sex, dmu, .variable, coverage, bias, CE),
    
    results %>%
      mutate(MSE = (.value - truth)^2, 
             bias = (.value - truth) / truth, 
             cover = ifelse(truth >= .lower & truth <= .upper, 1, 0)) %>% 
      filter((.variable %in% c("mlam", "mlam.dmu"))) %>%
      group_by(age, sex, dmu, .variable) %>%
      summarise(coverage = sum(cover) / (nsims), 
                bias = sum(bias) / (nsims), 
                MSE = sum(MSE) / nsims, 
                totalN = sum(truth)) %>%
      mutate(CE = sqrt(MSE) / (totalN / (nsims))) %>%
      select(age, sex, dmu, .variable, coverage, bias, CE)) 
  
  
  return(list(results = results, pm_datum = pm_datum))
    
  } else {
    
    return(list(results = results, pm_datum = NULL))
  }
  
  
}

# Wrapper function 
sim_deer_analysis <- function(init.N, mu.NS, sd.NS.dmu, sd.NS.yr, mu.HS, sd.HS.dmu, 
                    sd.HS.yr, mu.S0, sd.S0.dmu, sd.S0.yr, m, nyears, nsims, 
                    ndmus, mu.report, sd.report.dmu, sd.report.yr, correct.sex, 
                    prop.f, prop.m, ni, nb, model.scripts){
  
  # Simulate the deer population
  cat("Simulating deer populations... \n", "", "\n")
  
  sim_res <- deerpop_sim(init.N = init.N, mu.NS = mu.NS, sd.NS.dmu = sd.NS.dmu,
                         sd.NS.yr = sd.NS.yr, mu.HS = mu.HS, sd.HS.dmu = sd.HS.dmu, 
                         sd.HS.yr = sd.HS.yr, mu.S0 = mu.S0, sd.S0.dmu = sd.S0.dmu,
                         sd.S0.yr = sd.S0.yr, m = m, nyears = nyears, nsims = nsims, 
                         ndmus = ndmus)
  
  # Generate harvest data with errors
  cat("\n", "", "\n", "Generating and formatting harvest data... \n", "", "\n")
  
  harv_dat <- harvest_data_scramble(sim.results = sim_res, mu.report = mu.report,
                                    sd.report.dmu = sd.report.dmu, 
                                    sd.report.yr = sd.report.yr, 
                                    correct.sex = correct.sex, prop.f = prop.f,
                                    prop.m = prop.m)
  
  # Run the models for each simulated deer population
  cat("\n", "", "\n", "Running models to estimate abundance... \n", "", "\n")
  
  abund_est <- run_IPM(sim.results = harv_dat, mu.NS = mu.NS, mu.HS = mu.HS, 
                       mu.report = mu.report, mu.S0 = mu.S0, m = m, ni = ni, 
                       nb = nb, model.scripts = model.scripts)
  
  # Calculate performane metrics
  cat("\n", "", "\n", "Calculating performance metrics... \n", "", "\n")
  
  suppressMessages(pm_dat <- calc_pm(sim.results = harv_dat, abund.est = abund_est))
  
  
  
  return(list(sim_datum = harv_dat, ests = abund_est, pm_datum = pm_dat))
  
}

# Function to generate random initial values for deer population
rand_deer_init_pop <- function(ndmus, mu.init){
  init_N <- array(NA, dim = c(3, 2, ndmus))
  total_N <- NULL
  prop_f <- c(0.14, 0.13, 0.34)#c(0.15, 0.12, 0.27)
  prop_m <- c(0.14, 0.12, 0.13)#c(0.15, 0.10, 0.21)
  
  for(i in 1:ndmus){
    total_N[i] <- rnorm(1, mu.init, sd = 0.35 * mu.init)
    for(a in 1:3){
      init_N[,1,i] <- ceiling(total_N[i] * prop_f)
      init_N[,2,i] <- ceiling(total_N[i] * prop_m)
    }
  }
  
  return(init_N)
}

# Function to run model in JAGS. Useful if running many different models with 
# different data
run_model <- function(script, ni, nb, nt, nc, win.data, params, inits) {
  # Run the model
  suppressMessages(out <- jags(win.data, inits, params, script, 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
              jags.module = c("glm", "dic")))
  
  if(any(out$BUGSoutput$summary[,"Rhat"] > 1.05)){
    suppressMessages(out <- autojags(out, n.iter = 200000, n.thin = nt, n.update = 5,
                                   progress.bar = "none"))
  }
  
  return(out)
  
}

# Function to estimate/calculate initial population size for deer based on 
# harvest data (total and sexed/aged data), reporting rate, and harvest rate. 
est_init_pop <- function(y.H, y.Af, y.Am, harvest.rate = c(0.09, 0.20), 
                         report.rate = c(0.8, 0.76)){
  y.H <- array_tree(y.H, margin = 3)
  
  est_harvest <- y.H %>% 
    map(~ as.data.frame(.x) %>%
          slice_head(n = 3) %>%
          pivot_longer(1:2, names_to = "bag", values_to = "N") %>%
          group_by(bag) %>%
          summarise(avg = mean(N)) %>%
          arrange(desc(bag)) %>% ungroup() %>%
          mutate(hest = avg / (harvest.rate * report.rate)) %>%
          select(-avg) %>%
          arrange(desc(bag)) %>%
          split(.$bag) %>%
          map(~ .x %>% select(-bag)))
  
  
  yearly_props <- list("antlered" =
                         data.frame("antlered" = apply(abind(array_tree(y.Af, margin = 3) %>% 
                                                               replace(0, 1) %>%
                                                               map(~ as.data.frame(.x) %>% 
                                                                     mutate(across(where(is.numeric), 
                                                                                   prop.table))), 
                                                             along = 3), 1, mean)), 
                       "antlerless" =
                         data.frame("antlerless" = apply(abind(array_tree(y.Af, margin = 3) %>% 
                                                                 replace(0, 1) %>%
                                                                 map(~ as.data.frame(.x) %>% 
                                                                       mutate(across(where(is.numeric), 
                                                                                     prop.table))), 
                                                               along = 3), 1, mean)))
  
  
  init_pop <- abind::abind(map(est_harvest, 
                               ~ as.data.frame(do.call(cbind, map2(.x, yearly_props, 
                                                                   ~ round(mapply(`*`, .y, .x), 0)))) %>%
                                 select("antlerless", "antlered")), along = 3)
  
  for(dmu in 1:dim(init_pop)[3]){
    init_pop[1,,dmu] <- round(init_pop[2,,dmu] / (0.6 * (1 - harvest.rate)), 0)
  }
  
  init_pop <- log(init_pop)
  
  return(init_pop)
  
}


# Function to estimate/calculate initial population size for deer based on 
# harvest data (total and sexed/aged data), reporting rate, and harvest rate. 
init_pop_values <- function(y.H, y.Af, y.Am, harvest.rate = c(0.09, 0.20), 
                            report.rate = c(0.8, 0.76)){
  y.H <- array_tree(y.H, margin = 3)
  
  est_harvest <- y.H %>% 
    map(~ as.data.frame(t(.x)) %>%
          rownames_to_column(var = "bag") %>% 
          mutate(across(where(is.numeric), ~ .x / (harvest.rate * report.rate))) %>%
          split(.$bag) %>%
          map(~ .x %>% select(-bag)))
  
  
  yearly_props <- list("antlered" =
                         data.frame("antlered" = apply(abind(array_tree(y.Af, margin = 3) %>% 
                                                               replace(0, 1) %>%
                                                               map(~ as.data.frame(.x) %>% 
                                                                     mutate(across(where(is.numeric), 
                                                                                   prop.table))), 
                                                             along = 3), 1, mean)), 
                       "antlerless" =
                         data.frame("antlerless" = apply(abind(array_tree(y.Af, margin = 3) %>% 
                                                                 replace(0, 1) %>%
                                                                 map(~ as.data.frame(.x) %>% 
                                                                       mutate(across(where(is.numeric), 
                                                                                     prop.table))), 
                                                               along = 3), 1, mean)))
  
  
  init_pop <- abind::abind(
    map(est_harvest, 
        ~ aperm(apply(abind::abind(map2(.x, yearly_props, 
                                        ~ round(mapply(`*`, .y, .x, USE.NAMES = FALSE), 0)), 
                                   along = 3), 1:2, function(x) sort(x, decreasing = TRUE)), c(2,3,1))),
    along = 4)
  
  for(dmu in 1:dim(init_pop)[4]){
    for(bag in 1:dim(init_pop)[3]){
      for(yr in 1:(dim(init_pop)[2] - 1)){
        h <- harvest.rate[bag]
        init_pop[1,yr,bag,dmu] <- round(init_pop[2,yr+1,bag,dmu] / (0.6 * (1 - h)), 0)
      }
      init_pop[3,,bag,dmu] <- NA
      init_pop[,1,bag,dmu] <- NA
    }
  }
  
  return(init_pop)
  
}

# Wrapper function to generate priors for survival, productivity, and reporting 
# rate.
generate_priors <- function(values, dist.type, ...) {
  
  if(dist.type == "beta") {
    rates <- t(apply(values, 1, 
                     FUN = function(x) unlist(beta.MoM.fcn(mean = x[1], 
                                                           var = x[2]))))
    
  } else if(dist.type == "gamma") {
    rates <- t(apply(values, 1, 
                     FUN = function(x) unlist(gamma.MoM.fcn(mean = x[1], 
                                                            var = x[2]))))
  } else {
    warning("Please provide dist.type of either 'beta', or 'gamma'")
  }
  
  return(rates)
}

# MoM for beta distribution
beta.MoM.fcn <- function (mean, var, plot = FALSE){
  alpha <-mean * ((mean * (1 - mean) / var) - 1)
  beta <-(1 - mean) * ((mean * (1 - mean) / var) - 1)
  y <-rbeta(1:1000, alpha, beta)
  if(plot == TRUE){
    beta.plot <-plot(y, dbeta(y, alpha, beta), ylab = "Frequency", xlab = "Rate")
  }
  return(c(list(alpha = alpha, beta = beta)))
}


# MoM for gamma distribution
gamma.MoM.fcn <- function(mean, var, plot = FALSE) {
  shape <- (mean ^ 2) / var 
  rate <- mean / var
  y <- rgamma(1:1000, shape = shape, rate = rate)
  if(plot == TRUE){
    gamma.plot <- plot(y, dgamma(y, shape, rate), ylab = "Frequency", xlab = "Rate")
  }
  return(c(list(shape = shape, rate = rate)))
}
