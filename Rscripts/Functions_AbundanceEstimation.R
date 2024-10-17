##############################################################################L

#   Functions for TN Abundance Estimation Analysis



#   Allison C Keever
#   TNTech University
#   Last updated: 10/17/24

##############################################################################L


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
  
  
  yearly_props <- list(
    "antlered" = data.frame(
      "antlered" = apply(abind(array_tree(y.Am, margin = 3) %>% 
                                 replace(0, 1) %>%
                                 map(~ as.data.frame(.x) %>% 
                                       mutate(across(where(is.numeric), 
                                                     prop.table))), 
                               along = 3), 1, mean)), 
    "antlerless" = data.frame(
      "antlerless" = apply(abind(array_tree(y.Af, margin = 3) %>% 
                                   replace(0, 1) %>%
                                   map(~ as.data.frame(.x) %>% 
                                         mutate(across(where(is.numeric), 
                                                       prop.table))), 
                                 along = 3), 1, mean)))
  
  
  init_pop <- abind::abind(map(est_harvest, 
                               ~ as.data.frame(
                                 do.call(
                                   cbind,
                                   map2(.x, 
                                        yearly_props,
                                        ~ round(mapply(`*`, .y, .x), 0)) %>%
                                     assign_in(list(1), 
                                               rbind(pluck(., 2, 4), 
                                                     pluck(., 1)))  %>% 
                                     assign_in(list(2), 
                                               t(rbind(pluck(., 2)[-4,]))))) %>% 
                                 rename("antlerless" = V2) %>%
                                 select("antlerless", "antlered")), along = 3)
  
  for(dmu in 1:dim(init_pop)[3]){
    init_pop[1,,dmu] <- round(init_pop[2,,dmu] / (0.6 * (1 - harvest.rate)), 0)
  }
  
  init_pop <- init_pop * 0.88
  
  return(init_pop)
  
}


# Function to estimate/calculate initial population size for deer based on 
# harvest data (total and sexed/aged data), reporting rate, and harvest rate. 
init_pop_values <- function(y.H, y.Af, y.Am, harvest.rate = c(0.15, 0.30), 
                            report.rate = c(0.82, 0.78)){
  y.H <- array_tree(y.H, margin = 3)
  
  est_harvest <- y.H %>% 
    map(~ as.data.frame(t(.x)) %>%
          rownames_to_column(var = "bag") %>% 
          mutate(across(where(is.numeric), ~ .x / (harvest.rate * report.rate))) %>%
          split(.$bag) %>%
          map(~ .x %>% select(-bag)))
  
  
  yearly_props <- list("antlered" =
                         data.frame("antlered" = 
                                      apply(abind(
                                        array_tree(y.Am, margin = 3) %>% 
                                          replace(0, 1) %>%
                                          map(~ as.data.frame(.x) %>% 
                                                mutate(across(where(is.numeric), 
                                                              prop.table))), 
                                        along = 3), 1, mean)), 
                       "antlerless" =
                         data.frame("antlerless" = 
                                      apply(abind(
                                        array_tree(y.Af, margin = 3) %>% 
                                          replace(0, 1) %>%
                                          map(~ as.data.frame(.x) %>% 
                                                mutate(across(where(is.numeric), 
                                                              prop.table))), 
                                        along = 3), 1, mean)))
  
  
  init_pop <- abind::abind(
    map(est_harvest, 
        ~ aperm(apply(abind::abind(
          map2(.x, 
               yearly_props,
               ~ round(mapply(`*`, .y, .x, USE.NAMES = FALSE), 0)) %>%
            assign_in(list(1), 
                      rbind(pluck(., 2)[4,], 
                            pluck(., 1)))  %>% 
            assign_in(list(2), 
                      rbind(pluck(., 2)[-4,])),
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

