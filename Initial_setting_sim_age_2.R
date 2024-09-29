load("Data/All_dat.rda")

latent <-  3
bef_onset <-  3
self_heal <-  17

comp_length <-  8
time_length <-  60
age_length <- 7

contact_para <- dat_list$`Contact matrix`

Sus <- dat_list$Susceptbility

N <- dat_list$Population

dat <- dat_list$`Initial value`
remove_para <- dat_list$`Transmission period`

age <- c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60+")
start_date <- as.Date("2020-01-08")
end_date <- as.Date("2020-02-22")

dat <- c(dat)

gibbs <- function(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length){ 
  
  xi <- para_t[[1]]
  tr_beta <- para_t[[2]]
  ratio_as <-  para_t[[3]]
  remove <- para_t[[4]]
  remove_a <- para_t[[5]]
  period <- para_t[[6]]
  Sus <- para_t[[7]]
  contact_4 <- para_t[[8]]
  v <- para_t[[9]]
  cov_beta <- para_t[[10]]
  sigma <- para_t[[11]]
  max_l <- para_t[[12]]
  tol_lik <- para_t[[13]]
  contact <- para_t[[14]]
  h <- para_t[[15]] 
  mean_t <- para_t[[16]]
  
  para_matrix <- matrix(0, age_length * 8, age_length * 8)
  diag(para_matrix) <- c(rep(0, age_length), rep(- 1 / period[1], age_length),
                         rep(- 1 / period[2], age_length), - 1 / remove,
                         - (1 / period[3] + 1 / remove_a), rep(0, 3 * age_length))
  diag(para_matrix[15:21,8:14]) <- 1 / period[1] * (1 - ratio_as)
  diag(para_matrix[22:28,15:21]) <- rep(1 / period[2], age_length)
  diag(para_matrix[29:35,8:14]) <- 1 / period[1] * ratio_as
  diag(para_matrix[36:42,22:28]) <- 1 / remove
  diag(para_matrix[50:56,29:35]) <- 1 / period[3]
  diag(para_matrix[43:49,29:35]) <- 1 / remove_a
  para_matrix <- lapply(1:(time_length-1), function(t){
    para_matrix[8:14,15:21] <- para_matrix[8:14,22:28] <- t(t((N * Sus) * contact[[t]]) * (tr_beta[,1] / N))
    para_matrix[8:14,29:35] <- t(t(para_matrix[8:14,15:21]) * tr_beta[,2])
    return(para_matrix)
  })
  
  tr_beta_c <- log(tr_beta[,1]) - log(1 - tr_beta[,1])
  ratio_as_c <- log(ratio_as) - log(1 - ratio_as)
  para_age <- cbind(tr_beta_c, ratio_as_c)
  h <- h + 1
  for(hh in 1:7){
    cov_beta[[hh]] <- (h - 1) / h * cov_beta[[hh]] + 5.76 / length(para_age[hh,]) / h * 
      (h * tcrossprod(mean_t[[hh]][[1]]) - (h + 1) * tcrossprod(mean_t[[hh]][[2]]) +
         tcrossprod(c(para_age[hh,])) + diag(rep(10 ^ (-10), 2)))
    if(tr_beta[1,2] <= 1){
      min_T <- max(tr_beta[-hh,1]) /  tr_beta[1,2] / 30
      min_T <- log(min_T) - log(1 - min_T)
      max_T <- min(min(tr_beta[-hh,1]) * 30 * tr_beta[1,2], 1)
      max_T <- log(max_T) - log(1 - max_T)
    }else{
      min_T <- max(tr_beta[-hh,1]) *  tr_beta[1,2] / 30
      min_T <- log(min_T) - log(1 - min_T)
      max_T <- min(min(tr_beta[-hh,1]) * 30 / tr_beta[1,2], 1)
      max_T <- log(max_T) - log(1 - max_T)
    }
    para_age_t <- as.vector(rtmvnorm(1, para_age[hh,], cov_beta[[hh]], lb = c(min_T, - Inf), ub = c(max_T, Inf)))
    para_age_ct <- 1 / (1 + exp(- para_age_t))
    tr_beta_ct <- tr_beta_c
    tr_beta_ct[hh] <- c(para_age_t[1])
    tr_beta_t <- 1 / (1 + exp(- tr_beta_ct))
    tr_beta_t <- cbind(tr_beta_t, tr_beta[,2])
    para_matrix_t <- lapply(1:(time_length-1), function(t){
      para_matrix[[t]][8:14,14+hh] <- para_matrix[[t]][8:14,21+hh] <- N * Sus * contact[[t]][,hh] * para_age_ct[1] / N[hh]
      para_matrix[[t]][8:14,28+hh] <- para_matrix[[t]][8:14,14+hh] * tr_beta[hh,2]
      para_matrix[[t]][14+hh,7+hh] <- 1 / period[1] * (1 - para_age_ct[2])
      para_matrix[[t]][28+hh,7+hh] <- 1 / period[1] * para_age_ct[2]
      return(para_matrix[[t]])
    })
    l_t <- l(dat, new_tranf_ob, N, para_matrix_t, v)
    r <- l_t - max_l 
      # sum(dnorm(diff(tr_beta_ct, lag = 1, differences = 1), 0, 1 / sqrt(sigma[1]), log = T)) - 
      # sum(dnorm(diff(tr_beta_c, lag = 1, differences = 1), 0, 1 / sqrt(sigma[1]), log = T)) 
    U <- log(runif(1))
    if(U < r){
      para_age[hh,] <- para_age_t
      tr_beta_c <- tr_beta_ct
      tr_beta <- tr_beta_t
      para_matrix <- para_matrix_t
      max_l <- l_t
    }
    mean_t[[hh]][[1]] <- mean_t[[hh]][[2]]
    mean_t[[hh]][[2]] <- as.vector((h + 1) / (h + 2) * mean_t[[hh]][[1]] + c(para_age[hh,]) / (h + 2)) 
  }
  
  ratio_as <- 1 / (1 + exp(- para_age[,2]))
  sigma[1] <- rgamma(1, 1 + 0.5 * (age_length - 1), 1 / 100 + 0.5 * sum(diff(tr_beta_c, lag = 1, differences = 1) ^ 2))
  
  ratio <- max(tr_beta[,1]) / min(tr_beta[,1])
  for(rep in 1:3){
    if(runif(1) < 0.5){
      a <- 1.1
    }else{
      a <- 1.01
    }
    tr_beta_t <- rep(rlnormTrunc(1, log(tr_beta[1,2]), log(a), min = ratio / 30, max = 30 / ratio), age_length)
    para_matrix_t <- lapply(1:(time_length-1), function(t){
      para_matrix[[t]][8:14,29:35] <- t(t(para_matrix[[t]][8:14,15:21]) * tr_beta_t)
      return(para_matrix[[t]])
    })
    l_t <- l(dat, new_tranf_ob, N, para_matrix_t, v)
    r <- l_t - max_l +
      log(dlnormTrunc(tr_beta[1,2], meanlog = log(tr_beta_t[1]), sdlog = log(a), min = ratio / 30, max = 30 / ratio)) -
      log(dlnormTrunc(tr_beta_t[1], meanlog = log(tr_beta[1,2]), sdlog = log(a), min = ratio / 30, max = 30 / ratio)) 
    U <- log(runif(1))
    if((U < r)){
      tr_beta[,2] <- tr_beta_t
      para_matrix <- para_matrix_t
      max_l <- l_t
    }
  }
  
  if(runif(1) < 0.5){
    a <- 1.1
  }else{
    a <- 1.01
  }
  remove_a_i <- rep(rlnormTrunc(1, meanlog = log(remove_a[1]), sdlog = log(a), min = 0, max = 60), age_length)
  para_matrix_t <- lapply(1:(time_length-1), function(t){
    diag(para_matrix[[t]][29:35,29:35]) <- - (1 / period[3] + 1 / remove_a_i)
    diag(para_matrix[[t]][43:49,29:35]) <- 1 / remove_a_i
    return(para_matrix[[t]])
  })
  l_t <- l(dat, new_tranf_ob, N, para_matrix_t, v)
  r <- l_t - max_l +
    log(dlnormTrunc(remove_a[1], meanlog = log(remove_a_i[1]), sdlog = log(a), min = 0, max = 60)) -
    log(dlnormTrunc(remove_a_i[1], meanlog = log(remove_a[1]), sdlog = log(a), min = 0, max = 60))
  U <- log(runif(1))
  if(U < r){
    remove_a <- remove_a_i
    para_matrix <- para_matrix_t
    max_l <- l_t
  }
  
  if(runif(1) < 0.5){
    a <- 1.1
  }else{
    a <- 1.01
  }
  v_t <- rep(rlnormTrunc(1, meanlog = log(v[1]), sdlog = log(a), min = 0, max = 100), 2)
  l_t <- l(dat, new_tranf_ob, N, para_matrix, v_t)
  r <- l_t - max_l +
    log(dlnormTrunc(v[1], meanlog = log(v_t[1]), sdlog = log(a), min = 0, max = 100)) -
    log(dlnormTrunc(v_t[1], meanlog = log(v[1]), sdlog = log(a), min = 0, max = 100)) 
  U <- log(runif(1))
  if(U < r){
    v <- v_t
    max_l <- l_t
  }
  
  for(i in 2){
    if(runif(1) < 0.5){
      a <- 1.1
    }else{
      a <- 1.01
    }
    contact_4_t <- contact_4
    contact_4_t[i] <- rlnormTrunc(1, log(contact_4[i]), log(a), min = 0, max = 30)
    contact_t <- contact_func(contact_4_t[1], contact_4_t[2],
                              contact_4_t[3], contact_4_t[4], contact_para, time_length)
    para_matrix_t <- lapply(1:(time_length-1), function(t){
      para_matrix[[t]][8:14,15:21] <- para_matrix[[t]][8:14,22:28] <- t(t((N * Sus) * contact_t[[t]]) * (tr_beta[,1] / N))
      para_matrix[[t]][8:14,29:35] <- t(t(para_matrix[[t]][8:14,15:21]) * tr_beta[,2])
      return(para_matrix[[t]])
    })
    l_t <- l(dat, new_tranf_ob, N, para_matrix_t, v)
    r <- l_t - max_l +
      log(dlnormTrunc(contact_4[i], log(contact_4_t[i]), log(a), min = 0, max = 30)) -
      log(dlnormTrunc(contact_4_t[i], log(contact_4[i]), log(a), min = 0, max = 30))
    U <- log(runif(1))
    if(U < r){
      contact_4 <- contact_4_t
      contact <- contact_t
      para_matrix <- para_matrix_t
      max_l <- l_t
    }
  }
  
  for(i in 3){
    if(runif(1) < 0.5){
      a <- 1.1
    }else{
      a <- 1.01
    }
    contact_4_t <- contact_4
    contact_4_t[i] <- rlnormTrunc(1, log(contact_4[i]), log(a), min = 0, max = 2)
    contact_t <- contact_func(contact_4_t[1], contact_4_t[2],
                              contact_4_t[3], contact_4_t[4], contact_para, time_length)
    para_matrix_t <- lapply(1:(time_length-1), function(t){
      para_matrix[[t]][8:14,15:21] <- para_matrix[[t]][8:14,22:28] <- t(t((N * Sus) * contact_t[[t]]) * (tr_beta[,1] / N))
      para_matrix[[t]][8:14,29:35] <- t(t(para_matrix[[t]][8:14,15:21]) * tr_beta[,2])
      return(para_matrix[[t]])
    })
    l_t <- l(dat, new_tranf_ob, N, para_matrix_t, v)
    r <- l_t - max_l +
      log(dlnormTrunc(contact_4[i], log(contact_4_t[i]), log(a), min = 0, max = 2)) -
      log(dlnormTrunc(contact_4_t[i], log(contact_4[i]), log(a), min = 0, max = 2))  
    U <- log(runif(1))
    if(U < r){
      contact_4 <- contact_4_t
      contact <- contact_t
      para_matrix <- para_matrix_t
      max_l <- l_t
    }
  }
  
  tol_lik_t <- max_l + 
    sum(dnorm(diff(tr_beta_c, lag = 1, differences = 1), 0, 1 / sqrt(sigma[1]), log = T)) +
    sum(dgamma(sigma[1], 1, 1 / 100, log = T)) 
  
  return(
    list(
      xi,
      tr_beta,
      ratio_as,
      remove,
      remove_a,
      period,
      Sus,
      contact_4,
      v,
      cov_beta,
      sigma,
      max_l,
      c(tol_lik, tol_lik_t),
      contact,
      h,
      mean_t
    )
  )
}
