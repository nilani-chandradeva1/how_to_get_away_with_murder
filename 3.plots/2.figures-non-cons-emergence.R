model_results_df <- readRDS("2.output/model_results_df_constant_emergence_FALSE.rds")
model_results_base_df <- readRDS("2.output/model_results_base_df_constant_emergence_FALSE.rds")
constant_emergence <- 1 #if 0, phi = mu*M0, if 1, phi = mu*M
# Define model parameters here
params_base <- list(
  N0 = 1000,              # Initial human population (should remain constant)
  
  #M0 = 2000,              # Initial mosquito population. Change to change V/H ratio and endemicity
  M0 = c(200000, 2000),
  
  mu_h = 1 / (70 * 365),  # Human natural death rate
  ## (WHO Life Tables; 70y life expectancy)
  #
  sigma_h = 1 / 12,       # Human incubation rate
  ## (~12 days incubation; Smith et al. 2012)
  #
  gamma_h = 1 / 60,       # Human recovery rate
  ## (~60 days infectious; Griffin et al. 2010)
  #
  omega_h = 1 / 180,      # Human loss of immunity rate
  ## (~180 days immune; White et al. 2014)
  
  #beta_hv = 0.3,          # Daily transmission prob: human to mosquito
  ## (Smith et al. 2012)
  #a = 1/3,                #1 bite every 3 days (Le Menach)
  #b_hv = 0.05,
  #b_vh = 0.2,
  beta_hv = 0.3,
  beta_vh = 0.3,
  
  mu_v = c(0.1,0.03),              # Mosquito natural death rate or longer-living vector e.g a blackfly
  ## (~10 days lifespan; Lines et al. 1987)
  #
  sigma_v = 1 / 10,       # Mosquito incubation rate (EIP)
  ## (~10 days at 25°C; Guerra et al. 2010)
  #
  
  #beta_vh = 0.3,          # Daily transmission prob: mosquito to human
  ## (Smith et al. 2012)
  
  tau = 100,              # Start time of intervention (days)
  #delta_D = 1000,          # Target mosquitoes killed
  delta_D = c(500, 1000),
  #times = seq(0, 500, by = 0.1), # All time points inclusive of pre-intervention
  delta_t_vec = c(10, 30, 90) # Intervention durations (days)
  #delta_t_vec = 10
)
params_base$m0 <- with(params_base, {M0 / N0}) # Initial mosquito to human ratio
#params_base$beta_hv <- with(params_base, {b_hv*a})
#params_base$beta_vh <- with(params_base, {b_vh*a})
params_base$constant_emergence <- constant_emergence

#plots for analysis####

#first look at killing 1000 out of 2000 mosquitoes (50%) at 2:1 vector to host ratio (2000 mosq to 1000 people)
mu_v_vec <- unique(params_base$mu_v)
M0_vec <- unique(params_base$M0)
delta_D_vec <- unique(params_base$delta_D)
m0_vec <- unique(params_base$m0)

model_results_base_main <- model_results_base_df %>%
  filter(m0 == m0_vec[2] & mu_v == mu_v_vec[1] & M0 == M0_vec[2])

model_results_main <- model_results_df %>%
  filter(delta_D == delta_D_vec[2] & M0 == M0_vec[2] & mu_v == mu_v_vec[1] & m0 == m0_vec[2])

model_results_all_main <- model_results_main %>%
  left_join(model_results_base_main[,c("t", "C0", "mu_v", "M0", "m0")]) %>% # let it auto-select the columns to join by
  mutate(delta_C = C0-C, 
         rel_delta_C = ((C0-C)/C0)*100)

#take measurements of prevalence averted at day 10, day 30 and 90
scenario_pals <- c('#1b9e77','#d95f02','#7570b3') #colour for each product

#dynamics plots, split across two figures

model_results_all_main %>%
  group_by(delta_t) %>%
  summarise(mean_mosq = mean(M), 
            psi = unique(psi))

mosq_killed_dynamics <- ggplot(model_results_all_main, aes(x = t - params_base$tau, y = D, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = scenario_pals[2])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.7, 0.8))+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  ylab("Number of mosquitoes killed by \n intervention")+
  xlab("Time since intervention started (days)")+
  ylim(0, 2001)

prev_dynamics <- ggplot(model_results_all_main, aes(x = t-params_base$tau, y = (I_h/N)*100, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  #geom_line(data = model_results_base_df, aes(x = t, y = I_h/N, col = "baseline"), linetype = "dashed")+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = scenario_pals[2])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  guides(col = "none")+
  ylab("Prevalence (%)")+
  xlab("Time since intervention started (days)")+
  ylim(0, 24)

delta_C_dynamics <- ggplot(model_results_all_main, aes(x = t-params_base$tau, y = delta_C, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = scenario_pals[2])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  theme(legend.position = c(0.7, 0.2))+
  xlab("Time since intervention started (days)")+
  ylab("Difference in number of cases \n compared to baseline")

mosq_dynamics <- ggplot(model_results_all_main, aes(x = t- params_base$tau, y = M, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col =  scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col =  scenario_pals[2])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col =  scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  guides(col = "none")+
  xlab("Time since intervention started (days)")+
  ylab("Mosquito population size")+
  ylim(0, 2001)


Re_t_plot <- ggplot(model_results_all_main, aes(x = t-params_base$tau, y = Re_t, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = scenario_pals[2])+
  geom_vline(xintercept = 100 - params_base$tau+ params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  guides(col = "none")+
  xlab("Time since intervention started (days)")+
  ylab("Re(t)")+
  ylim(0, 2)

R0_t_plot <- ggplot(model_results_all_main, aes(x = t-params_base$tau, y = R0_t, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = scenario_pals[2])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  guides(col = "none")+
  xlab("Time since intervention started (days)")+
  ylab("R0(t)")+
  ylim(0, 60)

fig_1 <- cowplot::plot_grid(mosq_killed_dynamics, mosq_dynamics, 
                            
                            labels = c("A", "B"))

fig_2 <- cowplot:: plot_grid(prev_dynamics,delta_C_dynamics,
                             Re_t_plot, R0_t_plot, 
                             labels = c("A", "B", "C", "D"))
ggsave(fig_1, file = "3.plots/non_cons_emerg_pop_dynamics.svg")
ggsave(fig_2, file = "3.plots/non_cons_emerg_trans_dynamics.svg")

#which has a greater epi impact?

#checking dynamics
ggplot(model_results_base_df, aes(x =t, y = I_h/N, col = as.factor(label)))+
  geom_line()+
  ylim(0,1)

ggplot(model_results_base_df, aes(x =t, y = C0, col = as.factor(label)))+
  geom_line()

ggplot(model_results_df, aes(x =t, y = I_h/N, col = as.factor(label)))+
  geom_line()

ggplot(model_results_df, aes(x =t, y = C, col = as.factor(label)))+
  geom_line()

#summary stats for baseline scenario (nothing changes over time, so don't need to filter by time here)

model_base_epi <- model_results_base_df %>%
  group_by(m0, mu_v, M0) %>%
  summarise(mean_prev_baseline = mean(I_h/N), 
            mean_R0_t_baseline = mean(R0_t), 
            mean_Re_t_baseline = mean(Re_t), 
            mean_C0_baseline = mean(C0)) 

#summary stats for each intervention scenario, depending on days since int that we take the measurement
model_all_summary_d10 <- model_results_df %>%
  filter(between(t, params_base$tau, params_base$tau+params_base$delta_t_vec[1])) %>%
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t), 
    mean_C_int = mean(C)) %>%
  mutate(measurement_t = params_base$delta_t_vec[1])

model_all_summary_d30 <- model_results_df%>%
  filter(between(t, params_base$tau, params_base$tau+params_base$delta_t_vec[2])) %>%
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t), 
    mean_C_int = mean(C))%>%
  mutate(measurement_t = params_base$delta_t_vec[2])

model_all_summary_d90 <- model_results_df%>%
  filter(between(t, params_base$tau, params_base$tau+params_base$delta_t_vec[3])) %>%
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t), 
    mean_C_int = mean(C))%>%
  mutate(measurement_t = params_base$delta_t_vec[3])

model_all_summary_d150 <- model_results_df %>%
  filter(between(t, params_base$tau, 250)) %>% #all back to eqm now in mosq pop
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t), 
    mean_C_int = mean(C))%>%
  mutate(measurement_t = 150)


model_all_summary <- do.call("rbind", list(model_all_summary_d10, model_all_summary_d30, model_all_summary_d90, model_all_summary_d150))

#compare diff in prevalence
summary_impact <- left_join(model_all_summary, model_base_epi) %>%
  mutate(abs_diff_prev = mean_prev_baseline - mean_prev_int, 
         rel_diff_prev = ((mean_prev_baseline - mean_prev_int)/mean_prev_baseline)*100, 
         
         abs_diff_C = mean_C0_baseline - mean_C_int, 
         rel_diff_C = ((mean_C0_baseline - mean_C_int)/mean_C0_baseline)*100, 
         
         abs_diff_Re = mean_Re_t_baseline - mean_Re_t_int, 
         rel_diff_Re = ((mean_Re_t_baseline - mean_Re_t_int)/mean_Re_t_baseline)*100, 
         
         abs_diff_R0 = mean_R0_t_baseline - mean_R0_t_int, 
         rel_diff_R0 = ((mean_R0_t_baseline - mean_R0_t_int)/ mean_R0_t_baseline)*100)

scenario_pals2 <- c(scenario_pals, '#e7298a')

rel_diff_prev_plot <- summary_impact %>%
  filter(mu_v == mu_v_vec[1] & delta_D == delta_D_vec[2] & m0 ==m0_vec[2] & M0 == M0_vec[2]) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_prev, fill = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  xlab("Time of measurement (days) after intervention on")+
  ylab("Relative difference (%) in prevalence compared to baseline")+
  theme_bw()+
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)+
  guides(fill = "none")


rel_diff_Re <- summary_impact %>%
  filter(mu_v == mu_v_vec[1] & delta_D == delta_D_vec[2] & m0 ==m0_vec[2] & M0 == M0_vec[2]) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_Re, fill = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  xlab("Time of measurement (days) after intervention on")+
  ylab("Relative reduction (%) in Rt compared to baseline")+
  theme_bw()+
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)+
  theme(legend.position = c(0.7, 0.8))

rel_diff_C <- summary_impact %>%
  filter(mu_v == mu_v_vec[1] & delta_D == delta_D_vec[2] & m0 ==m0_vec[2] & M0 == M0_vec[2]) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_C, fill = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  xlab("Time of measurement (days) after intervention on")+
  ylab("Relative difference (%) in incidence compared to baseline")+
  theme_bw()+
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)+
  guides(fill = "none")


fig_3 <- cowplot::plot_grid(rel_diff_prev_plot, rel_diff_C, rel_diff_Re,
                            labels = c("A", "B", "C"), ncol = 3)

ggsave(fig_3, file = "3.plots/non_cons_emerg_epi_impact.svg")

#generic facet wrap

summary_impact %>%
  ggplot()+
  aes(x = as.factor(delta_t), y = abs_diff_prev, fill = as.factor(measurement_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(mu_v, delta_D, m0), labeller = label_both)+
  xlab("Duration of killing")+
  ylab("Absolute difference in prevalence compared to baseline")+
  theme_bw()+
  scale_fill_manual(name = "Time of measurement", labels = c("10 days after int on", 
                                                             "30 days after int on", 
                                                             "90 days after int on", 
                                                             "When mosq pop back to eqm"), 
                    values = scenario_pals2)

#another way to show this. 
low_endemicity_facet <- summary_impact %>%
  filter(m0 == 1) %>%
  ggplot()+
  aes(x = as.factor(delta_t), y = rel_diff_prev, fill = as.factor(measurement_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(mu_v, delta_D, m0), labeller = label_both)+
  xlab("Duration of killing")+
  ylab("Relative difference in prevalence compared to baseline (%)")+
  theme_bw()+
  scale_fill_manual(name = "Time of measurement", labels = c("10 days after int on", 
                                                             "30 days after int on", 
                                                             "90 days after int on", 
                                                             "When mosq pop back to eqm"), 
                    values = scenario_pals2)

#think it is clearer this way, emphasis is on how measured impact varies depending on time of measurement. 
low_endemicity_facet_prev <- summary_impact %>%
  filter(m0 == 1) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_prev, fill = as.factor(delta_t), 
      col = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(mu_v, delta_D, m0), labeller = label_both)+
  xlab("Time of measurement (days) after intervention on")+
  ylab("Relative difference in prevalence compared to baseline (%)")+
  theme_bw() +
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)+
  scale_colour_manual(name = "Intervention's duration of killing (days)", 
                      values = scenario_pals)+
  guides(col = "none")+
  theme(legend.position = c(0.1, 0.8))

ggsave(low_endemicity_facet_prev, file = "3.plots/cons_emerg_low_endem_facet_prev.svg")



low_endemicity_facet_Re <- summary_impact %>%
  filter(m0 == 1) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_Re, fill = as.factor(delta_t), 
      col = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(mu_v, delta_D, m0), labeller = label_both)+
  xlab("Time of measurement (days) after intervention on")+
  ylab("Relative difference in Re compared to baseline")+
  theme_bw() +
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)+
  scale_colour_manual(name = "Intervention's duration of killing (days)", 
                      values = scenario_pals)+
  guides(col = "none")+
  theme(legend.position = c(0.3, 0.1))

ggsave(low_endemicity_facet_Re, file = "3.plots/cons_emerg_low_endem_facet_Re.svg")

high_endemicity_facet_prev <- summary_impact %>%
  filter(m0 == 200) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_prev, fill = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(mu_v, delta_D, m0), labeller = label_both)+
  xlab("Duration of killing")+
  ylab("Relative difference in prevalence compared to baseline (%)")+
  theme_bw()+
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)

high_endemicity_facet_Re <- summary_impact %>%
  filter(m0 == 200) %>%
  ggplot()+
  aes(x = as.factor(measurement_t), y = rel_diff_Re, fill = as.factor(delta_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(mu_v, delta_D, m0), labeller = label_both)+
  xlab("Duration of killing")+
  ylab("Relative difference in Re compared to baseline (%)")+
  theme_bw()+
  scale_fill_manual(name = "Intervention's duration of killing (days)", 
                    values = scenario_pals)



#dynamics for longer-living vector
model_results_SM <- model_results_df %>%
  filter(delta_D == delta_D_vec[2] & M0 == M0_vec[2] & mu_v == mu_v_vec[2] & m0 == m0_vec[2])

model_results_all_SM <- model_results_SM %>%
  left_join(model_results_base_main[,c("t", "C0", "mu_v", "M0", "m0")]) %>% # let it auto-select the columns to join by
  mutate(delta_C = C0-C, 
         rel_delta_C = ((C0-C)/C0)*100)

Re_t_plot_SM <- ggplot(model_results_all_SM, aes(x = t-params_base$tau, y = Re_t, col = as.factor(delta_t)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = scenario_pals[1])+
  geom_vline(xintercept = 100- params_base$tau + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = scenario_pals[2])+
  geom_vline(xintercept = 100 - params_base$tau+ params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = scenario_pals[3])+
  geom_vline(xintercept = 250- params_base$tau, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  scale_colour_manual(name = "Duration of killing (days)", values = c(scenario_pals[1], 
                                                                      scenario_pals[2], 
                                                                      scenario_pals[3]))+
  guides(col = "none")+
  xlab("Time since intervention started (days)")+
  ylab("Re(t)")


#how R0 changes over time depending on mortality rate
model_results_df %>%
  filter(delta_D == delta_D_vec[2] & M0 == M0_vec[2] & m0 == m0_vec[2]) %>%
  ggplot()+
  aes(x = t-params_base$tau, y = Re_t, col = as.factor(psi), lty = as.factor(mu_v))+
  geom_line(linewidth = 1.1)+
  facet_wrap(vars(delta_t), labeller = label_both)+
  theme_bw()+
  xlab("Time since intervation started (days)")
#with a higher mortality rate, the rebound in Rt is much higher

R0_calc <- function(mu_v){
  m <- 2
  beta_hv = 0.3
  beta_vh = 0.3
  sigma_h = 1 / 12
  sigma_v = 1 / 10
  mu_h = 1 / (70 * 365)
  gamma_h = 1 / 60
  R0 <- (beta_hv*beta_vh*sigma_h*sigma_v)/((sigma_h + mu_h) * (gamma_h+mu_h) * (sigma_v + mu_v) * mu_v)
  R0_t <- R0*m
  df <- data.frame(mu_v, R0_t)
  return(df)
}

mu_v_input <- list(seq(0.10, 0.9, by = 0.1))
R0_out <- as.data.frame(lapply(mu_v_input, R0_calc))
ggplot(R0_out, aes(x = mu_v, y = R0_t))+
  geom_point()+
  geom_line()+
  theme_bw()+
  xlim(0.1, 0.9)