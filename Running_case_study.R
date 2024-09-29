# Preparation
setwd("~")
source("Modelling_function.R")

# Parameter estimation
## Age-related model
### Parallel sampling for burn-in
source("Initial_setting.R")
para_t <- init(1)
for(i in 1:3){
  para_t <- paral_sample(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length)
}

### Posterior sampling
K <- 40
sfInit(parallel = TRUE, cpus = K) 
sfSource("Modelling_function.R")
sfSource("Initial_setting.R")
sfExport("para_t")  
Result <- sfLapply(1:K, Poster_samp,
                   para_t = para_t, dat = dat, new_tranf_ob = new_tranf_ob, 
                   contact_para = contact_para, N = N,
                   time_length = time_length, age_length = age_length, K = K)
sfStop()

a <- 4000 / K
para <- lapply(1:(length(Result) * a), function(i){
  Result[[(ceiling(i / a))]][[(i - a * (ceiling(i / a) - 1))]]
})

save(para, file = paste0("Sample_Para", ".rda"), version = 2)

## Basic Model Comparison
### Parallel sampling for burn-in
source("Initial_setting_basic.R")
para_t <- init(1)
for(i in 1:3){
  para_t <- paral_sample(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length)
}

### Posterior sampling
K <- 40
sfInit(parallel = TRUE, cpus = K) 
sfSource("Modelling_function.R")
sfSource("Initial_setting_basic.R")
sfExport("para_t")  
Result <- sfLapply(1:K, Poster_samp,
                   para_t = para_t, dat = dat, new_tranf_ob = new_tranf_ob, 
                   contact_para = contact_para, N = N,
                   time_length = time_length, age_length = age_length, K = K)
sfStop()

a <- 4000 / K
para <- lapply(1:(length(Result) * a), function(i){
  Result[[(ceiling(i / a))]][[(i - a * (ceiling(i / a) - 1))]]
})

save(para, file = paste0("Sample_Para_2", ".rda"), version = 2)

## Basic Model 2 Comparison
### Parallel sampling for burn-in
source("Initial_setting_basic_2.R")
para_t <- init(1)
for(i in 1:3){
  para_t <- paral_sample(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length)
}

### Posterior sampling
K <- 40
sfInit(parallel = TRUE, cpus = K) 
sfSource("Modelling_function.R")
sfSource("Initial_setting_basic_2.R")
sfExport("para_t")  
Result <- sfLapply(1:K, Poster_samp,
                   para_t = para_t, dat = dat, new_tranf_ob = new_tranf_ob, 
                   contact_para = contact_para, N = N,
                   time_length = time_length, age_length = age_length, K = K)
sfStop()

a <- 4000 / K
para <- lapply(1:(length(Result) * a), function(i){
  Result[[(ceiling(i / a))]][[(i - a * (ceiling(i / a) - 1))]]
})

save(para, file = paste0("Sample_Para_3", ".rda"), version = 2)

# Model selection
source("Modelling_function.R")

load("Sample_Para.rda")
source("Initial_setting.R")
Model_select(para, dat, new_tranf_ob)

load("Sample_Para_2.rda")
source("initial_setting_basic.R")
Model_select(para, dat, new_tranf_ob)

load("Sample_Para_3.rda")
source("initial_setting_basic_2.R")
Model_select(para, dat, new_tranf_ob)

