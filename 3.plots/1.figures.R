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

scenario_pals <- c('#1b9e77','#d95f02','#7570b3') #colour for each product
scenario_pals2 <- c('#e7298a', scenario_pals) #baseline colour

mosq_killed_plot <- ggplot(df_all_main, aes(x = t-start_int, y = D, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme( 
        text = element_text(size = 14))+
  guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Number of mosquitoes killed by \n intervention")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  #xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  #ylim(0,2000)+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 300), ylim = c(0,2000))

inf_mosq <- ggplot(df_all_main, aes(x = t-start_int, y = (I_v/M)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme( 
    text = element_text(size = 14))+
  guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Infectious vectors (%)")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  #xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  #ylim(0,2000)+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 300), ylim = c(0,25))

susceptible_people <- ggplot(df_all_main, aes(x = t-start_int, y = (S_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme( 
    text = element_text(size = 14))+
  guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Susceptible humans (%)")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  #xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  #ylim(0,2000)+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 300), ylim = c(0,10))

daily_inc_plot <- ggplot(df_all_main, aes(x = t-start_int, y = C_daily, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.8, 0.5), 
        text = element_text(size = 14))+
  guides(col = "none", linetype = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Daily incidence")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  #xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 300))

prevalence_plot <- ggplot(df_all_main, aes(x = t-start_int, y = (I_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme( 
        text = element_text(size = 14))+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Prevalence(%) in humans")+
  #guides(col = "none", linetype = "none")+
  ylim(0, 24)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                       values = c("solid", "dotdash"))+
  #xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  #labs(col = "Time to complete MDA")+
  theme(legend.position = c(0.5, 0.5))+
  scale_colour_manual(values = scenario_pals2, labels = c("Baseline (no intervention)", "10 days", "30 days", "90 days"), 
                      name = "Time taken to kill target mosquitoes")+
  coord_cartesian(xlim = c(-10, 300))


mosq_pop_plot <- ggplot(df_all_main, aes(x = t-start_int, y = M, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(
        text = element_text(size = 14))+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Mosquito population size")+
  guides(col = "none", linetype = "none")+
  #ylim(0,2000)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 300), ylim = c(0,2000))


df_all_main %>%
  group_by(delta_t, constant_emergence) %>%
  summarise(min_M = min(M))

#mosq_Iv_plot <- ggplot(df_all_main, aes(x = t-start_int, y = (I_v/M)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
#  geom_line(linewidth = 0.9)+
#  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
#  theme_bw()+
#  #theme(legend.position = c(0.7, 0.3))+
#  ylab("Prevalence (%) in mosquitoes")+
#  guides(col = "none", linetype = "none")+
#  ylim(0, 21)+
#  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
#                        values = c("solid", "dotdash"))+
#  xlim(-10,300)+
#  xlab("Time since intervention started (days)")+
#  scale_colour_manual(values = scenario_pals2)

Re_t_plot <- ggplot(df_all_main, aes(x = t-start_int, y = Re_t, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Effective reproduction number (Re,t)")+
  guides(col = "none", linetype = "none")+
  ylim(0, 2)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)

figure_dynamics <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot, 
                                      daily_inc_plot, prevalence_plot, 
                                      Re_t_plot,
                                      align = "v",
                                      labels = c("A", "B", "C", "D"))

sens_dynamics_epi <- cowplot::plot_grid(inf_mosq, susceptible_people, 
                                        labels = c("A", "B"), 
                                        nrow = 2, 
                                        align = "v")

ggsave(sens_dynamics_epi, file = "3.plots/simple_model_plot_epi.pdf")
#then epi impact plot for:
#different m0 (vector to host ratio - defines endemicity)
#different proportion of total vectors killed
#each faceted by the time of the measurement

#epi impact

#measure impact at different time periods

eqm_point <- 350

#assign time periods 
time_ranges <- tibble(
  time_period = c("10d", "30d", "90d", "250d"), #number of days since start
  start = c(100, 100, 100, 100), #start times for measurements
  end = c(110, 130, 190, eqm_point) #end times for measurements
)


model_base_long <- base_scenarios_main %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)


model_base_epi <- model_base_long %>%
  group_by(m0, M0, time_period, constant_emergence) %>% #for each m0, M0 and time period, sum total cases
  summarise(tot_cases_baseline = sum(C_daily), .groups = "drop")

model_int_long <- int_scenarios_main %>%
  select(t,delta_t, m0, M0, delta_D, prop_killed, C_daily, constant_emergence) %>%
  filter(prop_killed == 0.95) %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)

model_int_summary <- model_int_long %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed, time_period, constant_emergence) %>%
  summarise(tot_cases_int = sum(C_daily), .groups = "drop")

#compare difference in incidence
summary_impact <- left_join(model_int_summary, model_base_epi) %>% #at each time period, for each scenario, what is rel diff in prev
  mutate(abs_diff_cases = tot_cases_baseline-tot_cases_int, 
         rel_diff_cases = ((tot_cases_baseline-tot_cases_int)/tot_cases_baseline)*100, 
         constant_emergence = case_when(constant_emergence == TRUE ~ "Constant emergence", 
                                        constant_emergence == FALSE ~ "Carrying capacity"))


#need to work on formatting of plot here and cross-check with the figure showed T & E last week
#make axis suitable for log scale

impact_plot <- ggplot(summary_impact, aes(x = factor(time_period, levels = c("10d", "30d", "90d", "250d")), y = log(rel_diff_cases+1), fill = as.factor(delta_t)))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(constant_emergence))+
  theme_minimal()+
  #ylim(0,4)+
  scale_y_continuous(limits = c(0,4), labels = c(0, 10, 20, 30, 40))+
  ylab("Efficacy (%)")+
  xlab("Time period over which incidence measured since intervention start")+
  labs(fill = "Time to complete MDA (days)")+
  theme(legend.position = c(0.8, 0.9))+
  scale_fill_manual(values = scenario_pals)+
  guides(fill = "none")

library(ggpattern)

impact_plot_main <- ggplot(summary_impact, 
                           aes(x = factor(time_period, levels = c("10d", "30d", "90d", "250d")),
                               y = log(rel_diff_cases + 1),
                               fill = as.factor(delta_t),
                               pattern = as.factor(constant_emergence))) +
  geom_bar_pattern(
    stat = "identity", 
    position = position_dodge(),
    colour = "black",                # Border of bars
    pattern_colour = "black",        # Pattern line color
    pattern_fill = NA,               # Transparent so bar fill shows
    pattern_density = 0.4,           
    pattern_spacing = 0.02,          
    pattern_key_scale_factor = 0.5
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 4), labels = c(0, 10, 20, 30, 40)) +
  ylab("Efficacy (%)") +
  xlab("Time period over which incidence measured since intervention start") +
  labs(fill = "Time to complete MDA (days)",
       pattern = "Adult emergence") +
  theme(legend.position = c(0.8, 0.9), 
        text = element_text(size = 14)) +
  scale_fill_manual(values = scenario_pals) +
  scale_pattern_manual(values = c("Carrying capacity" = "none",
                                  "Constant emergence" = "stripe"))+
  guides(pattern = guide_legend(
    override.aes = list(fill = "white"), # Force fill color in legend
  ), 
  #fill = guide_legend(override.aes = list(pattern = "none")), 
  fill = "none")

summary_impact %>%
  filter(time_period == "250d") %>%
  mutate(rel_rel_cases = case_when(constant_emergence == "Carrying capacity" ~ 0.210, 
                                  TRUE ~ 0.290), 
         abs_rel_diff = round(rel_diff_cases - rel_rel_cases, 2))



figure_dynamics <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot, 
                                      daily_inc_plot, prevalence_plot, 
                                     
                                      align = "v",
                                      labels = c("A", "B", "C", "D"))

figure_dynamic_impact <- cowplot::plot_grid(figure_dynamics, impact_plot_main, 
                                            labels = c("", "E"))

ggsave(figure_dynamic_impact, file = "3.plots/simple_model_plot.png")      
ggsave(figure_dynamic_impact, file = "3.plots/simple_model_plot.pdf") 

#sensitivity for higher V-H ratio but same proportion killed
model_base_long2 <- base_scenarios %>%
  filter(M0 == M0_vec[1]) %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)


model_base_epi2 <- model_base_long2 %>%
  group_by(m0, M0, time_period, constant_emergence) %>% #for each m0, M0 and time period, sum total cases
  summarise(tot_cases_baseline = sum(C_daily), .groups = "drop")

model_int_long2 <- int_scenarios %>%
  select(t,delta_t, m0, M0, delta_D, prop_killed, C_daily, constant_emergence) %>%
  filter(prop_killed == 0.95 & M0 == M0_vec[1]) %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)

model_int_summary2 <- model_int_long2 %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed, time_period, constant_emergence) %>%
  summarise(tot_cases_int = sum(C_daily), .groups = "drop")

#compare difference in incidence
summary_impact2 <- left_join(model_int_summary2, model_base_epi2) %>% #at each time period, for each scenario, what is rel diff in prev
  mutate(abs_diff_cases = tot_cases_baseline-tot_cases_int, 
         rel_diff_cases = ((tot_cases_baseline-tot_cases_int)/tot_cases_baseline)*100, 
         constant_emergence = case_when(constant_emergence == TRUE ~ "Constant emergence", 
                                        constant_emergence == FALSE ~ "Carrying capacity"))



impact_plot_sens <- ggplot(summary_impact2, 
       aes(x = factor(time_period, levels = c("10d", "30d", "90d", "250d")),
           #y = sqrt(rel_diff_cases),
          y = rel_diff_cases,
           fill = as.factor(delta_t),
           pattern = as.factor(constant_emergence))) +
  geom_bar_pattern(
    stat = "identity", 
    position = position_dodge(),
    colour = "black",                # Border of bars
    pattern_colour = "black",        # Pattern line color
    pattern_fill = NA,               # Transparent so bar fill shows
    pattern_density = 0.4,           
    pattern_spacing = 0.02,          
    pattern_key_scale_factor = 0.5
  ) +
  theme_minimal() +
  scale_y_sqrt(limits = c(0, 2.5), 
               breaks = c(0, 0.5, 1, 1.5, 2.0, 2.5), #original scale breaks
                     labels = c(0, 0.5, 1, 1.5, 2.0, 2.5)) +
  ylab("Efficacy (%)") +
  #ylim(0,2.25)+
  xlab("Time period over which incidence measured since intervention start") +
  labs(fill = "Time to complete MDA (days)",
       pattern = "Adult emergence") +
  theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = scenario_pals) +
  scale_pattern_manual(values = c("Carrying capacity" = "none",
                                  "Constant emergence" = "stripe"))+
  guides(pattern = guide_legend(
    override.aes = list(fill = "white"), # Force fill color in legend
  ), 
  fill = guide_legend(override.aes = list(pattern = "none")))



int_scenarios_sens <- int_scenarios %>%
  filter(prop_killed == 0.95 & M0 == M0_vec[1])



base_scenarios_sens <- base_scenarios %>%
  filter(M0 == M0_vec[1])

ggplot(base_scenarios_sens, aes(x = t, y = M, col = as.factor(constant_emergence)))+
  geom_line()


unique(base_scenarios_sens$delta_t)
unique(int_scenarios_sens$delta_t)


df_all_main_sens <- rbind(int_scenarios_sens, base_scenarios_sens)

#sensitivity plot



sens_mosq_killed_plot <- ggplot(df_all_main_sens, aes(x = t-start_int, y = D, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Number of mosquitoes killed by \n intervention")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)+
  ylim(0, 2e5)

sens_daily_inc_plot <- ggplot(df_all_main_sens, aes(x = t-start_int, y = C_daily, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.8, 0.5))+
  guides(col = "none", linetype = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Daily incidence")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)

sens_prevalence_plot <- ggplot(df_all_main_sens, aes(x = t-start_int, y = (I_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Prevalence(%) in humans")+
  #guides(col = "none", linetype = "none")+
  #ylim(0, 24)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  #labs(col = "Time to complete MDA")+
  theme(legend.position = c(0.5, 0.5))+
  scale_colour_manual(values = scenario_pals2, labels = c("Baseline (no intervention)", "10 days", "30 days", "90 days"), 
                      name = "Time taken to kill target mosquitoes")+
  ylim(0, 24)


sens_mosq_pop_plot <- ggplot(df_all_main_sens, aes(x = t-start_int, y = M, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Mosquito population size")+
  guides(col = "none", linetype = "none")+
  ylim(0,2e5)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)

sens_Re_t_plot <- ggplot(df_all_main_sens, aes(x = t-start_int, y = Re_t, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  #theme(legend.position = c(0.7, 0.3))+
  ylab("Effective reproduction number (Re,t)")+
  guides(col = "none", linetype = "none")+
  ylim(0, 2)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Carrying capacity", "Constant emergence"), 
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)

figure_dynamics_sens <- cowplot::plot_grid(sens_mosq_killed_plot, sens_mosq_pop_plot, 
                                      sens_daily_inc_plot, sens_prevalence_plot, 
                                      #sens_Re_t_plot,
                                      align = "v",
                                      labels = c("A", "B", "C", "D"))
figure_dynamics_sens_impact <- cowplot::plot_grid(figure_dynamics_sens, 
                                                 impact_plot_sens)
