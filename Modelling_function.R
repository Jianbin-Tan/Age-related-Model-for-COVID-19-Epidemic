# Package
library(deSolve)
library(extraDistr)
library(EnvStats)
library(truncdist)
library(tidyr)
library(ggplot2)
library(reshape2)
library(splines)
library(forcats)
library(dplyr)
library(viridis)
library(ggsci)
library(tidybayes)
library(coda)
library(ggalluvial)
library(scales)
library(xtable)
library(latex2exp)
library(snowfall)
library(expm)
library(mvtnorm)
library(TruncatedNormal)

# Function
## Observed data in likelihood
new_tranf <- function(Y){
  a <- apply(Y[,-1,c(4,6)] - Y[,-time_length,c(4,6)], c(1,2), sum) # New onset symptomatic
  c <- Y[,-1,7] - Y[,-time_length,7] # New comfirmed asymptomatic
  return(list(a,c))
}

## Contact function 
f_alp <- function(k, d, m, cont, down, up){
  o_contact <- cont[2] * down
  (cont[1] * up - o_contact) / (1 + exp(2 * log(99) / m * (k - d - m / 2))) + o_contact
}

contact_func <- function(d, m, down, up, contact_para, time_length){
  lapply(1:(time_length-1), function(k){
    apply(contact_para, c(2,3), f_alp, k = k, d = d, m = m, down = down, up = up)
  })
}

## ODE equation
eqn <- function(time, class, parms, period, contact, tr_beta, N, Sus,
                remove, remove_a, ratio){
  
  class <- matrix(class, age_length)
  
  a <- class[,2] / period[1] # latent
  b <- class[,3] / period[2] # before onset
  
  dS <- - contact %*% (rowSums(class[,3:5] * tr_beta[,c(1,1,2)]) / N) * (class[,1] * Sus)
  
  dRs <- class[,4] / remove
  dRa <- class[,5] / remove_a
  dRua <- class[,5] / period[3] # self healing
  
  dE <- - dS - a
  dIp <- a *  (1 - ratio) - b
  dIs <- b - dRs
  dIa <- a * ratio - dRa - dRua
  
  return(list(as.vector(cbind(dS, dE, dIp, dIs, dIa, dRs, dRa, dRua))))
}

## One step ODE
One_step <- function(t, dat, period, contact, tr_beta, N, Sus,
                     remove, remove_a, ratio){
  out <- as.numeric(ode(y = c(dat), times = (t-1):t, eqn, parms = NUll,
                        period = period, contact = contact,
                        tr_beta = tr_beta, N = N, Sus = Sus,
                        remove = remove, remove_a = remove_a,
                        ratio = ratio)[2,-1])
  return(matrix(out, age_length))
}

### Likelihood
l <- function(dat, new_tranf_ob, N, para_matrix, v){
  fit_new_tranf_1 <- vector()
  fit_new_tranf_2 <- vector()
  para_matrix <- lapply(1:(time_length-1), function(t) {expm(para_matrix[[t]], method = "Higham08")})
  for(t in 1:(time_length - 1)){
    dat_t <- para_matrix[[t]] %*% dat
    fit_new_tranf_1 <- c(fit_new_tranf_1, dat_t[22:28] + dat_t[36:42] - dat[22:28] - dat[36:42])
    fit_new_tranf_2 <- c(fit_new_tranf_2, dat_t[43:49] - dat[43:49])
    dat <- dat_t
  }
  logLik <- sum(dnbinom(new_tranf_ob[,1], mu = fit_new_tranf_1, size = 1 / v[1], log = T)) + 
    sum(dnbinom(new_tranf_ob[,2], mu = fit_new_tranf_2, size = 1 / v[2], log = T))
  
  if(is.finite(logLik) == F){
    logLik <- -Inf
  }
  return(logLik)
}

## Gibbs sampling
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
  para_matrix <- lapply(1:(time_length-1), function(t){
    if(t >= 21){
      diag(para_matrix[43:49,29:35]) <- 1 / remove_a
    }else{
      diag(para_matrix[29:35,29:35]) <- - 1 / period[3]
    }
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
    r <- l_t - max_l + 
      sum(dnorm(diff(tr_beta_ct, lag = 1, differences = 1), 0, 1 / sqrt(sigma[1]), log = T)) - 
      sum(dnorm(diff(tr_beta_c, lag = 1, differences = 1), 0, 1 / sqrt(sigma[1]), log = T)) 
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
    if(t >= 21){
      diag(para_matrix[[t]][29:35,29:35]) <- - (1 / period[3] + 1 / remove_a_i)
      diag(para_matrix[[t]][43:49,29:35]) <- 1 / remove_a_i
    }
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
    (- sum(abs(diff(tr_beta_c, lag = 1, differences = 1)) ^ 2) * sigma / 2) +
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
      mean_t,
      tol_lik_t
    )
  )
}

## Genearte initial value
init <- function(it){
  xi <- rep(1, 2)
  tr_beta <- cbind(rep(0.2, age_length), rep(0.5, age_length))
  ratio_as <- rep(0.5, age_length)
  remove <- remove_para[1,] / remove_para[2,]
  remove_a <- rep(14, age_length)
  period <- c(latent, bef_onset, self_heal)
  Sus <- Sus
  contact_4 <- c(0, 14, 1, 1)
  v <- c(0.01, 0.01)
  cov_beta <- lapply(1:age_length, function(i){
    diag(rep(10^(-2), 2))
  })
  mean_t <- lapply(1:age_length, function(i){
    list(c(0, 0),
         c(0, 0))
  })
  
  sigma <- c(30)
  max_l <- -10 ^ (10)
  tol_lik <- -10 ^ (10)
  contact <- contact_func(contact_4[1], contact_4[2], contact_4[3], contact_4[4], contact_para, time_length)
  
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
      tol_lik,
      contact, 
      1,
      mean_t
    )
  )
}

## Convergence test
paral_sample <- function(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length){
  
  P <- 0
  while(P != 1){
    para_t[[13]] <- para_t[[13]][length(para_t[[13]])]
    para_t[[15]] <- 1
    para_t[[16]] <- lapply(1:age_length, function(i){
      list(c(log(para_t[[2]][i,1]) - log(1 - para_t[[2]][i,1]),
             log(para_t[[3]][i]) - log(1 - para_t[[3]][i])),
           c(log(para_t[[2]][i,1]) - log(1 - para_t[[2]][i,1]),
             log(para_t[[3]][i]) - log(1 - para_t[[3]][i])))
    })
    
    para <- list()
    para[[1]] <- para_t
    for(h in 1:199){
      para[[h+1]] <- gibbs(para[[h]], dat, new_tranf_ob, contact_para, N, time_length, age_length)
    }
    
    P <- which.max(para[[200]][[13]])
    para_t <- para[[P]]
  }
  para_t <- para[[1]]
  
  return(para_t)
}

## Posterior sampling
Poster_samp <- function(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length, K, seed){
  para <- list()
  para[[1]] <- para_t
  G <- 4000 / K * 5
  for(h in 2:G){
    para[[h]] <- gibbs(para[[h-1]], dat, new_tranf_ob, contact_para, N, time_length, age_length)
  }
  para <- lapply(seq(5, G, 5), function(i){para[[i]]})
  return(para)
}

## Functions for prediction
pred_func <- function(para, dat, mark){
  sapply(mark, function(k){
    para_t <- para[[k]]
    xi <- para_t[[1]]
    tr_beta <- para_t[[2]]
    ratio_as <-  para_t[[3]]
    remove <- para_t[[4]]
    remove_a <- para_t[[5]]
    period <- para_t[[6]]
    Sus <- para_t[[7]]
    contact_4 <- para_t[[8]]
    v <- para_t[[9]]
    
    contact <- contact_func(contact_4[1], contact_4[2],
                            contact_4[3], contact_4[4], contact_para, time_length + pred_day)
    para_matrix <- matrix(0, age_length * 8, age_length * 8)
    diag(para_matrix) <- c(rep(0, age_length), rep(- 1 / period[1], age_length),
                           rep(- 1 / period[2], age_length), - 1 / remove,
                           - (1 / period[3] + 1 / remove_a), rep(0, 3 * age_length))
    diag(para_matrix[15:21,8:14]) <- 1 / period[1] * (1 - ratio_as)
    diag(para_matrix[22:28,15:21]) <- rep(1 / period[2], age_length)
    diag(para_matrix[29:35,8:14]) <- 1 / period[1] * ratio_as
    diag(para_matrix[36:42,22:28]) <- 1 / remove
    diag(para_matrix[50:56,29:35]) <- 1 / period[3]
    para_matrix <- lapply(1:(time_length-1 + pred_day), function(t){
      if(t >= 21){
        diag(para_matrix[43:49,29:35]) <- 1 / remove_a
      }else{
        diag(para_matrix[29:35,29:35]) <- - 1 / period[3]
      }
      para_matrix[8:14,15:21] <- para_matrix[8:14,22:28] <- t(t((N * Sus) * contact[[t]]) * (tr_beta[,1] / N))
      para_matrix[8:14,29:35] <- t(t(para_matrix[8:14,15:21]) * tr_beta[,2])
      return(para_matrix)
    })
    
    fit_new_tranf_1 <- vector()
    fit_new_tranf_2 <- vector()
    for(t in 1:(time_length - 1)){
      dat_t <- expm(para_matrix[[t]]) %*% dat
      fit_new_tranf_1 <- c(fit_new_tranf_1, dat_t[22:28] + dat_t[36:42] - dat[22:28] - dat[36:42])
      fit_new_tranf_2 <- c(fit_new_tranf_2, dat_t[43:49] - dat[43:49])
      dat <- dat_t
    }
    return(cbind(fit_new_tranf_1, fit_new_tranf_2))
  }, simplify = "array")
}

fit_tranf <- function(new_tranf_ob_age){
  apply(new_tranf_ob_age, c(2), sum)
}

## Model selection
Model_select <- function(para, dat, new_tranf_ob){
  para_median <- lapply(c(1:9, 11), function(k){
    sapply(1:length(para[[1]][[k]]), function(m){
      median(sapply(1:length(para), function(i){
        c(para[[i]][[k]])[m]
      }))
    })
  })
  
  para_median[[2]] <- matrix(para_median[[2]], age_length)
  tr_beta <- para_median[[2]]
  ratio_as <-  para_median[[3]]
  remove <- para_median[[4]]
  remove_a <- para_median[[5]]
  period <- para_median[[6]]
  Sus <- para_median[[7]]
  contact_4 <- para_median[[8]]
  v <- para_median[[9]]
  sigma <- para_median[[10]]
  
  contact <- contact_func(contact_4[1], contact_4[2],
                          contact_4[3], contact_4[4], contact_para, time_length)
  para_matrix <- matrix(0, age_length * 8, age_length * 8)
  diag(para_matrix) <- c(rep(0, age_length), rep(- 1 / period[1], age_length),
                         rep(- 1 / period[2], age_length), - 1 / remove,
                         - (1 / period[3] + 1 / remove_a), rep(0, 3 * age_length))
  diag(para_matrix[15:21,8:14]) <- 1 / period[1] * (1 - ratio_as)
  diag(para_matrix[22:28,15:21]) <- rep(1 / period[2], age_length)
  diag(para_matrix[29:35,8:14]) <- 1 / period[1] * ratio_as
  diag(para_matrix[36:42,22:28]) <- 1 / remove
  diag(para_matrix[50:56,29:35]) <- 1 / period[3]
  para_matrix <- lapply(1:(time_length-1), function(t){
    if(t >= 21){
      diag(para_matrix[43:49,29:35]) <- 1 / remove_a
    }else{
      diag(para_matrix[29:35,29:35]) <- - 1 / period[3]
    }
    para_matrix[8:14,15:21] <- para_matrix[8:14,22:28] <- t(t((N * Sus) * contact[[t]]) * (tr_beta[,1] / N))
    para_matrix[8:14,29:35] <- t(t(para_matrix[8:14,15:21]) * tr_beta[,2])
    return(para_matrix)
  })
  
  tr_beta_c <- log(tr_beta[,1]) - log(1 - tr_beta[,1])
  if(sum(diff(tr_beta_c, lag = 1, differences = 1) ^ 2) != 0){
    lik_tr <- sum(dnorm(diff(tr_beta_c, lag = 1, differences = 1), 0, 1 / sqrt(sigma), log = T))
  }else{
    lik_tr <- 0
  }
  a <- - 2 * (l(dat, new_tranf_ob, N, para_matrix, v) + lik_tr)
  
  b <- mean(sapply(1:length(para), function(k){
    if(sum(diff(para[[k]][[2]][,1], lag = 1, differences = 1) ^ 2) != 0){
      tr_beta_c <- log(para[[k]][[2]][,1]) - log(1 - para[[k]][[2]][,1])
      lik_tr <- sum(dnorm(diff(tr_beta_c, lag = 1, differences = 1), 0, 1 / sqrt(para[[k]][[11]]), log = T))
    }else{
      lik_tr <- 0
    }
    return(- 2 * (para[[k]][[12]] + lik_tr))
  })) 
  
  df <- b - a
  
  return(c(b, df, b + df, b + 2 * df))
}

## Simulated number
sim_cont <- function(age_length, time_length, comp_length, dat, period,
                     contact, tr_beta, N, Sus, remove, remove_a, ratio_as){
  Y <- array(0, c(age_length, time_length, comp_length))
  Y[,1,1:comp_length] <- matrix(dat, age_length)
  for(t in 1:(time_length - 1)){
    if(t < 21){
      Y[,t+1,] <- One_step(t, Y[,t,], period, contact[[t]], tr_beta, N, Sus,
                           remove, rep(Inf, age_length), ratio_as)
    }else{
      Y[,t+1,] <- One_step(t, Y[,t,], period, contact[[t]], tr_beta, N, Sus,
                           remove, remove_a, ratio_as)
    }
  }
  return(sum(Y[,time_length,-1]))
}

