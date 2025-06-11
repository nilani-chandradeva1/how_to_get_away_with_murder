require(tidyverse)

df_1_emerge_T <- readRDS("2.output/model_results_base_df_constant_emergence_TRUE.rds") #baseline constant emergence T
df_2_emerge_T <- readRDS("2.output/model_results_df_constant_emergence_TRUE.rds") #int constant emergence T

df_3_emerge_F <- readRDS("2.output/model_results_base_df_constant_emergence_FALSE.rds") #baseline constant emergence F
df_4_emerge_F <- readRDS("2.output/model_results_df_constant_emergence_FALSE.rds") #int constant emergence F

int_scenarios <- rbind(df_2_emerge_T, df_4_emerge_F)

int_scenarios_main <- int_scenarios %>%
  filter(delta_D == 1000 & M0 == 2000)

mosq_killed_plot <- ggplot(int_scenarios_main, aes(x = t, y = D, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.8, 0.5))+
  labs(col = "Duration of killing (days)", lty = "Constant adult emergence")+
  ylab("Number of mosquitoes killed by \n intervention")

prevalence_plot <- ggplot(int_scenarios_main, aes(x = t, y = (I_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Prevalence(%) in humans")+
  guides(col = "none", linetype = "none")+
  ylim(21, 24)

mosq_pop_plot <- ggplot(int_scenarios_main, aes(x = t, y = M, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Mosquito population size")+
  guides(col = "none", linetype = "none")

mosq_Iv_plot <- ggplot(int_scenarios_main, aes(x = t, y = (I_v/M)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Prevalence (%) in mosquitoes")+
  guides(col = "none", linetype = "none")+
  ylim(15, 21)

figure_dynamics <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot, 
                                      mosq_Iv_plot, prevalence_plot, 
                                      labels = c("A", "B", "C", "D"))
#then epi impact plot for:
#different m0 (vector to host ratio - defines endemicity)
#different proportion of total vectors killed
#each faceted by the time of the measurement

ggplot(df_1_emerge_T, aes(x = t, y = I_h/N, col = as.factor(m0)))+
  geom_line()

df_1_emerge_T_summary <- df_1_emerge_T %>%
  group_by(m0) %>%
  summarise(mean_prev = mean(I_h/N))

df_2_emerge_T_summary <- df_2_emerge_T %>%
  filter(between(t, 100, 100+30)) %>% ##timing for measuring difference. Mean over whole period - intervention gets washed out
  group_by(delta_t, delta_D, M0, m0) %>%
  summarise(mean_prev = mean(I_h/N))
