load("Data/All_dat.rda")

contact_para <- dat_list$`Contact matrix`
Sus <- dat_list$Susceptbility
N <- dat_list$Population
dat <- dat_list$`Initial value`
Y <- dat_list$`Observed data`
remove_para <- dat_list$`Transmission period`

age <- c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60+")
start_date <- as.Date("2020-01-08")
end_date <- as.Date("2020-02-22")

latent <-  3
bef_onset <-  3
self_heal <- 17

comp_length <-  8
time_length <-  47
age_length <- 7

new_tranf_ob <- new_tranf(Y)
dat <- c(dat)
new_tranf_ob <- cbind(c(new_tranf_ob[[1]]), c(new_tranf_ob[[2]]))

