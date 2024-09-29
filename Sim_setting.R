## Preparation
source("Modelling_function_sim.R")
source("Initial_setting_sim_age.R")

set.seed(138)
## parameter matrix
ratio_as <- rbeta(age_length, 2, 2) 
remove <- remove_para[1,] / remove_para[2,]
contact_4 <- c(0, 21, 0.5, 1)
period <- c(latent, bef_onset, self_heal)

remove_a <- 14
ratio_beta <- rep(0.5, age_length)
v <- rep(runif(1, 0, 100), 2)

tr_beta <- 0.4 / (exp(- bs(1:age_length, degree = 4, Boundary.knots = c(0, 8), intercept = F) %*% rnorm(4, 0, 1)) + 1) 
tr_beta <- cbind(tr_beta, tr_beta * ratio_beta)
tr_beta
contact <- contact_func(contact_4[1], contact_4[2], contact_4[3], contact_4[4], contact_para, time_length)
v <- rep(0.01, 2)

Y <- array(0, c(age_length, time_length, comp_length))
Y[,1,] <- dat
for(t in 1:(time_length - 1)){
  Y[,t+1,] <- One_step(t, Y[,t,], period, contact[[t]], tr_beta, N, Sus,
                       remove, remove_a, ratio_as)
}
fit_new_tranf <- new_tranf(Y)
true <- list(tr_beta[,1], tr_beta[,2], ratio_as, rep(remove_a, age_length), contact_4, v)
Ture_value <- list(New_num = fit_new_tranf, 
                   dis_par = v, 
                   init = c(dat),
                   con_par = contact_para, 
                   population = N,
                   Other_par = true)
colSums(fit_new_tranf[[1]])
colSums(fit_new_tranf[[2]])

save(Ture_value, file = paste0("Ture_value_sim", ".rda"), version = 2)
