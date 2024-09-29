# Preparation
## Set your address
setwd("~")

## Load the true values of parameters for simulations
load("Ture_value_sim.rda")

## Input some R packages
library(snowfall)

# Simulation
## Simulations for age-related model
K <- 45
sfInit(parallel = TRUE, cpus = K)
sfSource("Modelling_function_sim.R")
sfSource("Initial_setting_sim_age.R")
sfExport("Ture_value")  
Result <- sfLapply(floor(seq(200, 70000, length.out = 100)), 
                   sim_function, Ture_value = Ture_value)
sfStop()

save(Result, file = paste0("Result_sim_1", ".rda"), version = 2)

new_tranf_ob_tol <- lapply(1:length(Result), function(k){
  Result[[k]][[2]]
})

## Simulations for basic model
sfInit(parallel = TRUE, cpus = K)
sfSource("Modelling_function_sim.R")
sfSource("Initial_setting_sim_basic.R")
sfExport("Ture_value")  
sfExport("new_tranf_ob_tol") 
Result_2 <- sfLapply(1:100, sim_function_2, Ture_value = Ture_value, 
                     new_tranf_ob_tol = new_tranf_ob_tol)
sfStop()

save(Result_2, file = paste0("Result_sim_2", ".rda"), version = 2)

## Simulations for basic model 2
sfInit(parallel = TRUE, cpus = K)
sfSource("Modelling_function_sim.R")
sfSource("Initial_setting_sim_basic_2.R")
sfExport("Ture_value")  
sfExport("new_tranf_ob_tol") 
Result_3 <- sfLapply(1:100, sim_function_2, Ture_value = Ture_value, 
                     new_tranf_ob_tol = new_tranf_ob_tol)
sfStop()

save(Result_3, file = paste0("Result_sim_3", ".rda"), version = 2)

## Comparison
true <- Ture_value$Other_par
true[[2]] <- true[[2]] / true[[1]]

for(k in 1:length(Result)){
  for(j in 1:length(Result[[1]][[1]])){
    Result[[k]][[1]][[j]][[1]] <- Result[[k]][[1]][[j]][[2]][,2]
    Result[[k]][[1]][[j]][[2]] <- Result[[k]][[1]][[j]][[2]][,1]
  }
}

re_1 <- lapply(c(2, 1, 3, 5, 8, 9), function(m){
  matrix(sapply(1:length(Result), function(k){
    sapply(1:length(Result[[k]][[1]][[1]][[m]]), function(i){
      median(sapply(1:length(Result[[k]][[1]]), function(j){
        Result[[k]][[1]][[j]][[m]][i]
      }))
    })
  }), nrow = length(Result[[1]][[1]][[1]][[m]]))
})

re_1 <- lapply(1:length(re_1), function(m){
  sapply(1:nrow(re_1[[m]]), function(k){
    a <- abs((mean(re_1[[m]][k,] - true[[m]][k])) / true[[m]][k] * 100)
    b <- sqrt(mean((re_1[[m]][k,] - mean(re_1[[m]][k,])) ^ 2)) / true[[m]][k] * 100
    c <- sqrt(mean((re_1[[m]][k,] - true[[m]][k]) ^ 2)) / true[[m]][k] * 100
    return(round(c(a, b, c), 2))
  })
})

for(k in 1:length(Result_2)){
  for(j in 1:length(Result_2[[1]])){
    Result_2[[k]][[j]][[1]] <- Result_2[[k]][[j]][[2]][,2]
    Result_2[[k]][[j]][[2]] <- Result_2[[k]][[j]][[2]][,1]
  }
}

re_2 <- lapply(c(2, 1, 3, 5, 8, 9), function(m){
  matrix(sapply(1:length(Result_2), function(k){
    sapply(1:length(Result_2[[k]][[1]][[m]]), function(i){
      median(sapply(1:length(Result_2[[k]]), function(j){
        Result_2[[k]][[j]][[m]][i]
      }))
    })
  }), nrow = length(Result_2[[1]][[1]][[m]]))
})

re_2 <- lapply(1:length(re_2), function(m){
  sapply(1:nrow(re_2[[m]]), function(k){
    a <- abs((mean(re_2[[m]][k,] - true[[m]][k])) / true[[m]][k] * 100)
    b <- sqrt(mean((re_2[[m]][k,] - mean(re_2[[m]][k,])) ^ 2)) / true[[m]][k] * 100
    c <- sqrt(mean((re_2[[m]][k,] - true[[m]][k]) ^ 2)) / true[[m]][k] * 100
    return(round(c(a, b, c), 2))
  })
})

for(k in 1:length(Result_3)){
  for(j in 1:length(Result_3[[1]])){
    Result_3[[k]][[j]][[1]] <- Result_3[[k]][[j]][[2]][,2]
    Result_3[[k]][[j]][[2]] <- Result_3[[k]][[j]][[2]][,1]
  }
}

re_3 <- lapply(c(2, 1, 3, 5, 8, 9), function(m){
  matrix(sapply(1:length(Result_3), function(k){
    sapply(1:length(Result_3[[k]][[1]][[m]]), function(i){
      median(sapply(1:length(Result_3[[k]]), function(j){
        Result_3[[k]][[j]][[m]][i]
      }))
    })
  }), nrow = length(Result_3[[1]][[1]][[m]]))
})

re_3 <- lapply(1:length(re_3), function(m){
  sapply(1:nrow(re_3[[m]]), function(k){
    a <- abs((mean(re_3[[m]][k,] - true[[m]][k])) / true[[m]][k] * 100)
    b <- sqrt(mean((re_3[[m]][k,] - mean(re_3[[m]][k,])) ^ 2)) / true[[m]][k] * 100
    c <- sqrt(mean((re_3[[m]][k,] - true[[m]][k]) ^ 2)) / true[[m]][k] * 100
    return(round(c(a, b, c), 2))
  })
})

dat_latex <- rbind(cbind(t(re_1[[3]]), t(re_2[[3]]), t(re_3[[3]])),
                   cbind(t(re_1[[1]]), t(re_2[[1]]), t(re_3[[1]])),
                   c(t(re_1[[2]])[1,], t(re_2[[2]])[1,], t(re_3[[2]])[1,]),
                   c(t(re_1[[4]])[1,], t(re_2[[4]])[1,], t(re_3[[4]])[1,]),
                   cbind(t(re_1[[5]])[2:3,], t(re_2[[5]])[2:3,], t(re_3[[5]])[2:3,])
                   )
xtable(dat_latex) # Table 1

## More comparison for the age-fused procedure
sfInit(parallel = TRUE, cpus = K)
sfSource("Modelling_function_sim.R")
sfSource("Initial_setting_sim_age_2.R")
sfExport("Ture_value")  
sfExport("new_tranf_ob_tol") 
Result_4 <- sfLapply(1:100, sim_function_2, Ture_value = Ture_value, 
                     new_tranf_ob_tol = new_tranf_ob_tol)
sfStop()

save(Result_4, file = paste0("Result_sim_4", ".rda"), version = 2)

for(k in 1:length(Result_4)){
  for(j in 1:length(Result_4[[1]])){
    Result_4[[k]][[j]][[1]] <- Result_4[[k]][[j]][[2]][,2]
    Result_4[[k]][[j]][[2]] <- Result_4[[k]][[j]][[2]][,1]
  }
}

re_4 <- lapply(2, function(m){
  matrix(sapply(1:length(Result_4), function(k){
    sapply(1:length(Result_4[[k]][[1]][[m]]), function(i){
      median(sapply(1:length(Result_4[[k]]), function(j){
        Result_4[[k]][[j]][[m]][i]
      }))
    })
  }), nrow = length(Result_4[[1]][[1]][[m]]))
})

re_1 <- lapply(2, function(m){
  matrix(sapply(1:length(Result), function(k){
    sapply(1:length(Result[[k]][[1]][[1]][[m]]), function(i){
      median(sapply(1:length(Result[[k]][[1]]), function(j){
        Result[[k]][[1]][[j]][[m]][i]
      }))
    })
  }), nrow = length(Result[[1]][[1]][[1]][[m]]))
})

## Figure 3
Ture_age_dis <- Ture_value$Other_par[[1]]
fit_age_dis_sub_1 <- cbind(c((re_1[[1]])), rep(1, 7 * 100), rep(1:7, 100))
fit_age_dis_sub_2 <- cbind(c((re_4[[1]])), rep(2, 7 * 100), rep(1:7, 100))

dat_plot <- as.data.frame(rbind(fit_age_dis_sub_1, fit_age_dis_sub_2))
colnames(dat_plot) <- c("y", "class", "age")
dat_plot$class <- factor(dat_plot$class, levels = c("1", "2"),
                         labels = c("With shrinkage priors", "Without shrinkage priors"))
dat_plot$age <- factor(dat_plot$age, levels = c("1", "2", "3", "4", "5", "6", "7"),
                       labels = age)
library(ggsci)
set.seed(100)
ggplot() +
  geom_jitter(aes(x = dat_plot$age, y = dat_plot$y, color = dat_plot$class), position = position_jitterdodge(
    jitter.width = NULL,
    jitter.height = 0,
    dodge.width = 1,
    seed = NA
  )) +
  geom_point(aes(x = rep(age, 100), y = rep(Ture_age_dis, 100), fill = "True value"), size = 3) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Age", y = "Numbers", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1),
         color = T) +
  scale_color_manual(values = pal_jco()(2))
ggsave(paste0("Plot/age_sim.pdf"), width = 6.5, height = 4, dpi = 300)

