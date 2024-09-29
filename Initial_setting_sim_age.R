load("Data/All_dat.rda")

contact_para <- dat_list$`Contact matrix`
Sus <- dat_list$Susceptbility
N <- dat_list$Population
dat <- dat_list$`Initial value`
remove_para <- dat_list$`Transmission period`

age <- c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60+")
start_date <- as.Date("2020-01-08")
end_date <- as.Date("2020-02-22")

latent <-  3
bef_onset <-  3
self_heal <-  17

comp_length <-  8
time_length <-  60
age_length <- 7

dat <- c(dat)

