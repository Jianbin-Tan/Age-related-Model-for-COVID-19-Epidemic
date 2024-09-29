# Preparation
setwd("~")
source("Modelling_function.R")

## Sensitivity analysis for susceptibility
for(label in 1:3){
  source("Initial_setting.R")
  para_t <- init(1)
  para_t[[7]] <- para_t[[7]] + rnorm(age_length, 0, 0.1)
  para_t[[7]][para_t[[7]] < 0] <- 0
  for(i in 1:3){
    para_t <- paral_sample(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length)
  }
  
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
  
  save(para, file = paste0("Sample_Para_sen_", label, ".rda"), version = 2)
}

### Plot
result <- list()
for(label in 1:3){
  load(paste0("Sample_Para_sen_", label, ".rda"))
  result[[label]] <- para
}

x <- 1:7
mark <- 1:length(result[[1]])
tran <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1]
  })
})

tran_2 <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1] * result[[k]][[i]][[2]][,2]
  })
})

data_plot <- lapply(1:3, function(k){
  sym_Susceptibility <- tran[[k]]
  asym_Susceptibility <- tran_2[[k]]
  
  rownames(sym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  sym_Susceptibility <- mcmc(t(sym_Susceptibility))
  rownames(asym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  asym_Susceptibility <- mcmc(t(asym_Susceptibility))
  
  name_i <- factor(c(" Symptomatic", 'Asymptomatic'), ordered = F)
  a1 <- as.data.frame(spread_draws(sym_Susceptibility, age[i]))
  a1$class <- rep(name_i[1], nrow(a1))
  a2 <- as.data.frame(spread_draws(asym_Susceptibility, age[i]))
  a2$class <- rep(name_i[2], nrow(a2))
  dat_plot_1 <- rbind(a1, a2)[,c(1,2,6)]
  dat_plot_1$i <- as.factor(dat_plot_1$i)
  
  dat_plot_1$Case <- rep(paste0("Case ", k), nrow(a1))
  return(dat_plot_1)
})

data_plot <- rbind(data_plot[[1]], data_plot[[2]], data_plot[[3]])

ggplot(data_plot, aes(x = i, y = age, color = class)) +
  geom_boxplot(outlier.alpha = 0) +
  facet_grid(~Case) +
  scale_x_discrete(breaks = 1:7, label = age) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 12),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Age", y = "Transmissibility", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1),
         color = T) + 
  scale_color_manual(values = c("#de2d26", "#2b8cbe")) +
  scale_y_continuous(breaks = seq(0, 0.5, length.out = 6), limits = c(0, 0.5))

ggsave(paste0("Plot/Sen_1.pdf"), width = 10, height = 4, dpi = 300)


## Sensitivity analysis for contact number
for(label in 4:6){
  source("Initial_setting.R")
  contact_para_t <- contact_para + array(rnorm(length(contact_para), 0, 0.1), dim = dim(contact_para))
  contact_para_t[contact_para_t < 0] <- 0
  para_t <- init(1)
  for(i in 1:3){
    para_t <- paral_sample(para_t, dat, new_tranf_ob, contact_para_t, N, time_length, age_length)
  }
  
  K <- 40
  sfInit(parallel = TRUE, cpus = K) 
  sfSource("Modelling_function.R")
  sfSource("Initial_setting.R")
  sfExport("para_t", "contact_para_t")  
  Result <- sfLapply(1:K, Poster_samp,
                     para_t = para_t, dat = dat, new_tranf_ob = new_tranf_ob, 
                     contact_para = contact_para_t, N = N,
                     time_length = time_length, age_length = age_length, K = K)
  sfStop()
  
  a <- 4000 / K
  para <- lapply(1:(length(Result) * a), function(i){
    Result[[(ceiling(i / a))]][[(i - a * (ceiling(i / a) - 1))]]
  })
  
  save(para, file = paste0("Sample_Para_sen_", label, ".rda"), version = 2)
}

### Plot
result <- list()
for(label in 4:6){
  load(paste0("Sample_Para_sen_", label, ".rda"))
  result[[label - 3]] <- para
}

x <- 1:7
mark <- 1:length(result[[1]])
tran <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1]
  })
})

tran_2 <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1] * result[[k]][[i]][[2]][,2]
  })
})

data_plot <- lapply(1:3, function(k){
  sym_Susceptibility <- tran[[k]]
  asym_Susceptibility <- tran_2[[k]]
  
  rownames(sym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  sym_Susceptibility <- mcmc(t(sym_Susceptibility))
  rownames(asym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  asym_Susceptibility <- mcmc(t(asym_Susceptibility))
  
  name_i <- factor(c(" Symptomatic", 'Asymptomatic'), ordered = F)
  a1 <- as.data.frame(spread_draws(sym_Susceptibility, age[i]))
  a1$class <- rep(name_i[1], nrow(a1))
  a2 <- as.data.frame(spread_draws(asym_Susceptibility, age[i]))
  a2$class <- rep(name_i[2], nrow(a2))
  dat_plot_1 <- rbind(a1, a2)[,c(1,2,6)]
  dat_plot_1$i <- as.factor(dat_plot_1$i)
  
  dat_plot_1$Case <- rep(paste0("Case ", k), nrow(a1))
  return(dat_plot_1)
})

data_plot <- rbind(data_plot[[1]], data_plot[[2]], data_plot[[3]])

ggplot(data_plot, aes(x = i, y = age, color = class)) +
  geom_boxplot(outlier.alpha = 0) +
  facet_grid(~Case) +
  scale_x_discrete(breaks = 1:7, label = age) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Age", y = "Transmissibility", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1),
         color = T) + 
  scale_color_manual(values = c("#de2d26", "#2b8cbe")) +
  scale_y_continuous(breaks = seq(0, 0.5, length.out = 6), limits = c(0, 0.5))

ggsave(paste0("Plot/Sen_2.pdf"), width = 10, height = 4, dpi = 300)

## Sensitivity analysis for observed new numbers
for(label in 7:9){
  source("Initial_setting.R")
  new_tranf_ob_t <- matrix(rpois(length(new_tranf_ob ),new_tranf_ob), nrow = nrow(new_tranf_ob))
  new_tranf_ob_t[new_tranf_ob_t < 0] <- 0
  para_t <- init(1)
  for(i in 1:3){
    para_t <- paral_sample(para_t, dat, new_tranf_ob_t, contact_para, N, time_length, age_length)
  }
  
  K <- 40
  sfInit(parallel = TRUE, cpus = K) 
  sfSource("Modelling_function.R")
  sfSource("Initial_setting.R")
  sfExport("para_t", "new_tranf_ob_t")  
  Result <- sfLapply(1:K, Poster_samp,
                     para_t = para_t, dat = dat, new_tranf_ob = new_tranf_ob_t, 
                     contact_para = contact_para, N = N,
                     time_length = time_length, age_length = age_length, K = K)
  sfStop()
  
  a <- 4000 / K
  para <- lapply(1:(length(Result) * a), function(i){
    Result[[(ceiling(i / a))]][[(i - a * (ceiling(i / a) - 1))]]
  })
  
  save(para, file = paste0("Sample_Para_sen_", label, ".rda"), version = 2)
}

### Plot
result <- list()
for(label in 7:9){
  load(paste0("Sample_Para_sen_", label, ".rda"))
  result[[label - 6]] <- para
}

x <- 1:7
mark <- 1:length(result[[1]])
tran <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1]
  })
})

tran_2 <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1] * result[[k]][[i]][[2]][,2]
  })
})

data_plot <- lapply(1:3, function(k){
  sym_Susceptibility <- tran[[k]]
  asym_Susceptibility <- tran_2[[k]]
  
  rownames(sym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  sym_Susceptibility <- mcmc(t(sym_Susceptibility))
  rownames(asym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  asym_Susceptibility <- mcmc(t(asym_Susceptibility))
  
  name_i <- factor(c(" Symptomatic", 'Asymptomatic'), ordered = F)
  a1 <- as.data.frame(spread_draws(sym_Susceptibility, age[i]))
  a1$class <- rep(name_i[1], nrow(a1))
  a2 <- as.data.frame(spread_draws(asym_Susceptibility, age[i]))
  a2$class <- rep(name_i[2], nrow(a2))
  dat_plot_1 <- rbind(a1, a2)[,c(1,2,6)]
  dat_plot_1$i <- as.factor(dat_plot_1$i)
  
  dat_plot_1$Case <- rep(paste0("Case ", k), nrow(a1))
  return(dat_plot_1)
})

data_plot <- rbind(data_plot[[1]], data_plot[[2]], data_plot[[3]])


ggplot(data_plot, aes(x = i, y = age, color = class)) +
  geom_boxplot(outlier.alpha = 0) +
  facet_grid(~Case) +
  scale_x_discrete(breaks = 1:7, label = age) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 12),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Age", y = "Transmissibility", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1),
         color = T) + 
  scale_color_manual(values = c("#de2d26", "#2b8cbe")) +
  scale_y_continuous(breaks = seq(0, 0.5, length.out = 6), limits = c(0, 0.5))

ggsave(paste0("Plot/Sen_3.pdf"), width = 10, height = 4, dpi = 300)

## Sensitivity analysis for periods
for(label in 10:12){
  source("Initial_setting.R")
  para_t <- init(1)
  para_t[[6]] <- rlnorm(3, log(para_t[[6]]), 0.1)
  for(i in 1:3){
    para_t <- paral_sample(para_t, dat, new_tranf_ob, contact_para, N, time_length, age_length)
  }
  
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
  
  save(para, file = paste0("Sample_Para_sen_", label, ".rda"), version = 2)
}

### Plot
result <- list()
for(label in 10:12){
  load(paste0("Sample_Para_sen_", label, ".rda"))
  result[[label-9]] <- para
}

x <- 1:7
mark <- 1:length(result[[1]])
tran <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1]
  })
})

tran_2 <- lapply(1:3, function(k){
  sapply(mark, function(i){
    result[[k]][[i]][[2]][,1] * result[[k]][[i]][[2]][,2]
  })
})

data_plot <- lapply(1:3, function(k){
  sym_Susceptibility <- tran[[k]]
  asym_Susceptibility <- tran_2[[k]]
  
  rownames(sym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  sym_Susceptibility <- mcmc(t(sym_Susceptibility))
  rownames(asym_Susceptibility) <- sapply(1:age_length, function(i){paste0("age", "[", i, "]")})
  asym_Susceptibility <- mcmc(t(asym_Susceptibility))
  
  name_i <- factor(c(" Symptomatic", 'Asymptomatic'), ordered = F)
  a1 <- as.data.frame(spread_draws(sym_Susceptibility, age[i]))
  a1$class <- rep(name_i[1], nrow(a1))
  a2 <- as.data.frame(spread_draws(asym_Susceptibility, age[i]))
  a2$class <- rep(name_i[2], nrow(a2))
  dat_plot_1 <- rbind(a1, a2)[,c(1,2,6)]
  dat_plot_1$i <- as.factor(dat_plot_1$i)
  
  dat_plot_1$Case <- rep(paste0("Case ", k), nrow(a1))
  return(dat_plot_1)
})

data_plot <- rbind(data_plot[[1]], data_plot[[2]], data_plot[[3]])

ggplot(data_plot, aes(x = i, y = age, color = class)) +
  geom_boxplot(outlier.alpha = 0) +
  facet_grid(~Case) +
  scale_x_discrete(breaks = 1:7, label = age) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 12),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Age", y = "Transmissibility", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1),
         color = T) + 
  scale_color_manual(values = c("#de2d26", "#2b8cbe")) +
  scale_y_continuous(breaks = seq(0, 0.5, length.out = 6), limits = c(0, 0.5))

ggsave(paste0("Plot/Sen_4.pdf"), width = 10, height = 4, dpi = 300)
