## Preparation
setwd("~")
source("Modelling_function.R")
source("Initial_setting.R")
load("Sample_Para.rda")
mark <- 1:length(para)

## New cases plot
a <- matrix(new_tranf_ob[,1], age_length)
rownames(a) <- 1:7
colnames(a) <- 1:(time_length-1)
a <- melt(a)
a <- cbind(a, rep(1, nrow(a)))

b <- matrix(new_tranf_ob[,2], age_length)
rownames(b) <- 1:7
colnames(b) <- 1:(time_length-1)
b <- melt(b)
b <- cbind(b, rep(2, nrow(b)))

dat_fit <- rbind(as.matrix(a), as.matrix(b))
dat_fit <- as.data.frame(dat_fit)
colnames(dat_fit) <- c("age", "time", "number", "mark")
dat_fit$mark <- factor(
  dat_fit$mark,
  labels = c("New onset symptomatic cases",
             "Newly comfirmed asymptomatic cases"),
  levels = 1:2
)
dat_fit$age <- factor(
  dat_fit$age,
  labels = age,
  levels = 1:7
)

ggplot(dat_fit) + 
  geom_bar(aes(x = time, y = number, fill = age), stat = "identity") + 
  facet_wrap(~mark, nrow = 2) + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 15),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Time", y = "Numbers", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1),
         color = T) +
  scale_x_continuous(breaks = seq(1, time_length - 1, 5), 
                     label = format(seq.Date(from = start_date, by = "day", length.out = time_length), "%m/%d")[seq(1, time_length - 1, 5)]) + 
  scale_fill_manual(values = c('#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd'))

ggsave(paste0("Plot/News.pdf"), width = 8, height = 6, dpi = 300) # Figure 1

## Contact function
contact <- contact_func(0, 21, 1, 1, contact_para, time_length)
dat_contact <- sapply(1:(time_length-1), function(t){
  c(contact[[t]])
})
rownames(dat_contact) <- 1:49
colnames(dat_contact) <- 1:(time_length -1)
dat_contact <- melt(dat_contact)
colnames(dat_contact) <- c("age", "time", "Numbers")
dat_contact$age <- factor(dat_contact$age)

ggplot(dat_contact) + 
  geom_line(aes(x = time, y = Numbers, color = age), size = 0.8) +
  geom_vline(xintercept = 21, linetype = 2, size = 1, alpha = 0.8, color = "#999999") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Time", y = "Numbers", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(),
         fill = guide_legend(nrow = 1, order = 2),
         # shape = guide_legend(nrow = 1),
         size = guide_legend(nrow = 1)) +
  theme(legend.position = "none",
        text = element_text(size = 15)) +
  scale_x_continuous(breaks = seq(1, time_length, 5), 
                     label = format(seq.Date(from = start_date, by = "day", length.out = time_length), "%m/%d")[seq(1, time_length, 5)]) 

ggsave(paste0("Plot/Contact.pdf"), width = 8, height = 6, dpi = 300) # Figure 1 in the Supplementary Material

## Transmissibility
x <- 1:7
sym_Susceptibility <- sapply(mark, function(i){
  para[[i]][[2]][,1]
})

asym_Susceptibility <- sapply(mark, function(i){
  para[[i]][[2]][,1] * para[[i]][[2]][,2]
})

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

ggplot(dat_plot_1, aes(x = i, y = age, color = class)) +
  geom_boxplot(outlier.alpha = 0) +
  scale_x_discrete(breaks = 1:7, label = age) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 10),
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
ggsave(paste0("Plot/Trans.pdf"), width = 4.5, height = 3, dpi = 300) # Figure 3

## Transmission contribution
para_median <- lapply(1:9, function(k){
  sapply(1:length(para[[1]][[k]]), function(m){
    median(sapply(1:length(para), function(i){
      c(para[[i]][[k]])[m]
    }))
  })
})

para_median[[2]] <- matrix(para_median[[2]], age_length)
tr_beta <- para_median[[2]]
tr_beta[,2] <- tr_beta[,2] * tr_beta[,1]
ratio_as <-  para_median[[3]]
remove <- para_median[[4]]
remove_a <- para_median[[5]]
period <- para_median[[6]]
Sus <- para_median[[7]]
contact_4 <- para_median[[8]]

contact <- contact_func(contact_4[1], contact_4[2],
                        contact_4[3], contact_4[4], contact_para, time_length)

base_num <- sim_cont(age_length, time_length, comp_length, dat, period,
                     contact, tr_beta, N, Sus, remove, remove_a, ratio_as)
age_dat <- matrix(0, age_length, 2)
for(i in 1:7){
  for(j in 1:2){
    tr_beta_sim <- tr_beta
    tr_beta_sim[i,j] <- tr_beta_sim[i,j] / 2
    age_dat[i,j] <- (base_num - sim_cont(age_length, time_length, comp_length, dat, period,
                                         contact, tr_beta_sim, N, Sus, remove, remove_a, ratio_as)) / base_num
  }
}

colnames(age_dat) <- c("Symptomatic"," Asymptomatic")
rownames(age_dat) <- age

rowSums(age_dat)
colMeans(age_dat)
xtable(round(cbind(age_dat, rowMeans(age_dat)) * 100, 2)) # Table 3 in Supplementary Material

age_dat <- melt(age_dat)
colnames(age_dat) <- c("age", "mark", "Proportion")

ggplot(age_dat) +
  geom_bar(aes(x = age, y = Proportion, fill = mark),
           stat = "identity", position = position_dodge(), width = 0.9, linetype = 1) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 20),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "", y = "Proportion", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1, order = 2),
         shape = "none",
         size = guide_legend(nrow = 1),
         color = T) +
  scale_y_continuous(breaks = seq(0, 0.5, length.out = 6), limits = c(0, 0.6), 
                     label = c("0", "10%", "20%", "30%", "40%", "50%")) + 
  scale_fill_manual(values = c("#de2d26", "#2b8cbe"))

ggsave(paste0("Plot/Prop_contri.pdf"), width = 7, height = 4, dpi = 300) # Figure 4(A)

## Symptomatic transition prabability
dat_plot_2 <- sapply(mark, function(k){
  1 - para[[k]][[3]]
})

dat_plot_2 <- sapply(1:age_length, function(i){
  quantile(dat_plot_2[i,], c(0.025, 0.5, 0.975))
})

dat_plot_2 <- data.frame(t(dat_plot_2), rep("Symptomatic transition rate", age_length))
colnames(dat_plot_2) <- c("low", "median", "upper", "Symptomatic transition rate")

ggplot(dat_plot_2) +
  geom_bar(aes(x = x, y = median, fill = "Symptomatic transition probability"),
           stat = "identity", position = position_dodge(), width = 0.5, linetype = 1) +
  geom_errorbar(aes(x = x, ymax = upper, ymin = low),
                stat = "identity", position = position_dodge(0.9), width = 0.1, size = 0.5) +
  scale_x_continuous(breaks = 1:7, label = age) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 20),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "", y = "Proportion", 
       title = "", 
       colour = "", fill = "") + 
  guides(linetype = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1, order = 2),
         shape = "none",
         size = guide_legend(nrow = 1),
         color = T) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.2), 
                     label = c("0", "20%", "40%", "60%", "80%", "100%")) + 
  scale_fill_manual(values = "#fec44f")

ggsave(paste0("Plot/Prop_asym.pdf"), width = 7, height = 4, dpi = 300) # Figure 4(B)

