##############################################################################L

#   Functions for TN Abundance Estimation Analysis



#   Allison C Keever
#   TNTech University
#   Last updated: 5/25/23

##############################################################################L


# Function to run model in JAGS. Useful if running many different models with 
# different data
run_model <- function(script, ni, nb, nt, nc, win.data, params, inits) {
  # Run the model
  out <- jags(win.data, inits, params, script, 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
              jags.module = c("glm", "dic"))
  
  if(any(out$BUGSoutput$summary[,"Rhat"] > 1.05)){
    out <- autojags(out, n.iter = 200000, n.thin = nt, n.update = 5)
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
  
  init_pop <- init_pop
  
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