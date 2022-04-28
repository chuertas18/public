rm(list=ls())
###################################################################################################
# Paracou local allometry model process
###################################################################################################
library(brms)
library(data.table)

###################################################################################################
# 1/ Databases
###################################################################################################
# Establish working directory
setwd("Y:/users/ClaudiaHuertas/AllometryData/Database_unified_ALS_FTH_James_v14") 

# Alometry database
# data_allo <-fread("Database_unified_ALS_FTH_James_v14_5.csv")
###################################################################################################
# Model 1/ Single allometry for the entire site 
###################################################################################################
# Model 1 - Power law
model1_pw = brm(
  bf(H~alpha*dbh^beta, ## Power law
     alpha ~ 1,
     beta ~ 1,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 1000, 
  warmup = 500, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


# Model 1 - Michaelis Menten
model1_MM = brm(
  bf(H~(alpha*dbh)/(beta+dbh), ## 
     alpha ~ 1,
     beta ~ 1,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 1000, 
  warmup = 500, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


# Model 1 - Weillbull 2 parameters
model1_w2 = brm(
  bf(H~alpha*(1-exp(-dbh/beta)), ## Weillbull two parameters
     alpha ~ 1,
     beta ~ 1,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 1000, 
  warmup = 500, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)

###################################################################################################
# 3/ Second model (model2) includes canopy height in a local neighborhood (HC)
###################################################################################################
# Model 2 - Power law
model2_pw = brm(
  bf(H~alpha*dbh^beta, ## Power law
     alpha ~ 1 + HC,
     beta ~ 1,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 2000, 
  warmup = 1000, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


# Model 2 - Michaelis Menten
model2_MM = brm(
  bf(H~(alpha*dbh)/(beta+dbh), 
     alpha ~ 1 + HC,
     beta ~ 1,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 2000, 
  warmup = 1000, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


# Model 2 - Weillbull 2 parameters
model2_w2 = brm(
  bf(H~alpha*(1-exp(-dbh/beta)), ## Weillbull two parameters
     alpha ~ 1 + HC,
     beta ~ 1,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 2000, 
  warmup = 1000, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


###################################################################################################
# 3/ model 2b species
###################################################################################################
# Model 2b- Power law
model2b_pw_nsing= brm(
  bf(H~alpha*dbh^beta,
     alpha ~ 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
)



# Model 2b Michaelis-Menten
model2b_MM_nsing = brm(
  bf(H~(alpha*dbh)/(beta+dbh),
     alpha ~ 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  cores = 7,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
)



# Model 2 Weilbull
model2b_w2_nsing = brm(
  bf(H~alpha*(1-exp(-dbh/beta)),
     alpha ~ 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  cores = 7,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
  # backend = "cmdstanr",
  # threads = threading(4)
)

###################################################################################################
# 4/ Third model HC + species
###################################################################################################
# Model 3- Power law
model3_pw_nsing= brm(
  bf(H~alpha*dbh^beta,
     alpha ~ HC + 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
)

saveRDS(model3_pw_nsing,file = "D:/temp/bayesian/model3_pw_nsing")

# Model 3 Michaelis-Menten
# HC + species(random effects
model3_MM_nsing = brm(
  bf(H~(alpha*dbh)/(beta+dbh),
     alpha ~ HC + 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
)


model3_w2_nsing = brm(
  bf(H~alpha*(1-exp(-dbh/beta)),
     alpha ~ HC + 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_allo,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
)




