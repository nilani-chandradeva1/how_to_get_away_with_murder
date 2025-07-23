require(tidyverse)

df_1_emerge_T <- readRDS("2.output/model_results_base_df_constant_emergence_TRUE.rds")  #baseline constant emergence T
df_2_emerge_T <- readRDS("2.output/model_results_df_constant_emergence_TRUE.rds") #int constant emergence T

df_3_emerge_F <- readRDS("2.output/model_results_base_df_logistic_growth.rds") #baseline constant emergence F
df_4_emerge_F <- readRDS("2.output/model_results_df_logistic_growth.rds") #int constant emergence F
 
df_1_emerge_T <- df_1_emerge_T %>%
  select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D) %>%
  mutate(C0_daily = c(C0[1], diff(C0))) %>%
  ungroup()

df_2_emerge_T <- df_2_emerge_T %>%
  select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed) %>%
  mutate(C_daily = c(C[1], diff(C))) %>%
  ungroup()

df_3_emerge_F <- df_3_emerge_F %>%
  #select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D) %>%
  mutate(C0_daily = c(C0[1], diff(C0))) %>%
  ungroup()

df_4_emerge_F <- df_4_emerge_F %>%
  #select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed) %>%
  mutate(C_daily = c(C[1], diff(C))) %>%
  ungroup()

int_scenarios <- rbind(df_2_emerge_T, df_4_emerge_F)



#add on information for baseline or intervention, so can show baseline line
#convert incidence to daily so can plot that too. 

delta_D_vec <- unique(int_scenarios$delta_D)
M0_vec <- unique(int_scenarios$M0)

int_scenarios_main <- int_scenarios %>%
  filter(delta_D == delta_D_vec[2] & M0 == M0_vec[2])



base_scenarios <- rbind(df_1_emerge_T, df_3_emerge_F) %>%
  rename(C = C0, 
         C_daily = C0_daily)

base_scenarios_main <- base_scenarios %>%
  filter(M0 == M0_vec[2])

df_all_main <- rbind(int_scenarios_main, base_scenarios_main)
unique(df_all_main$delta_t)


start_int <- 100



mosq_killed_plot <- ggplot(df_all_main, aes(x = t-start_int, y = D, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.8, 0.5))+
  labs(col = "Duration of killing (days)")+
  ylab("Number of mosquitoes killed by \n intervention")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-99,400)+
  xlab("Time since intervention started (days)")+
  ylim(0,2000)

daily_inc_plot <- ggplot(df_all_main, aes(x = t-start_int, y = C_daily, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.8, 0.5))+
  guides(col = "none", linetype = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Number of mosquitoes killed by \n intervention")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-99,400)+
  xlab("Time since intervention started (days)")

prevalence_plot <- ggplot(df_all_main, aes(x = t-start_int, y = (I_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Prevalence(%) in humans")+
  guides(col = "none", linetype = "none")+
  ylim(0, 24)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                       values = c("solid", "dotdash"))+
  xlim(-99,400)+
  xlab("Time since intervention started (days)")


mosq_pop_plot <- ggplot(df_all_main, aes(x = t-start_int, y = M, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Mosquito population size")+
  guides(col = "none", linetype = "none")+
  ylim(0,2000)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-99,400)+
  xlab("Time since intervention started (days)")

mosq_Iv_plot <- ggplot(df_all_main, aes(x = t-start_int, y = (I_v/M)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Prevalence (%) in mosquitoes")+
  guides(col = "none", linetype = "none")+
  ylim(0, 21)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-99,400)+
  xlab("Time since intervention started (days)")

figure_dynamics <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot, 
                                      daily_inc_plot, prevalence_plot, 
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
