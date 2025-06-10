require(odin)
require(tidyverse)

constant_emergence <- 0 #if 0, phi = mu*M0, if 1, phi = mu*M
odin::can_compile()

# Define model parameters here
params_base <- list(
  N0 = 1000,              # Initial human population (should remain constant)
  
  #M0 = 2000,              # Initial mosquito population. Change to change V/H ratio and endemicity
  M0 = c(1000,2000),
  
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

# --- BASELINE PSI CALCULATION ---
# This calculates the value of psi required for equal birth and death rates.
# Note, to account for limits in machine precision, a lower bound is placed on
# 1 - delta_D_target/M0 = epsilon, where epsilon is the smallest non-zero
# normalised floating-point number. Also applied when delta_D_target >= M0.
baseline_psi <- function(Delta_t, delta_D_target, M0) {
  log_value <- max((M0 - delta_D_target) / M0, .Machine$double.xmin)
  -log(log_value) / Delta_t # Return value of psi
}

# --- FITTING FUNCTION FOR PSI ---
# This section calculates psi for a given a constant emergence rate.

# This is the objective function to minimise
psi_objective_function <- function(psi, Delta_t, delta_D_target, mu_v, M0) {
  term1 <- mu_v * Delta_t
  term2 <- (psi / (psi + mu_v)) * (1 - exp(-(psi + mu_v) * Delta_t))
  D_pred <- (psi * M0 / (psi + mu_v)) * (term1 + term2)
  (delta_D_target - D_pred)^2
}




# This fits the value of psi that minimises the objective function
fit_psi <- function(Delta_t, delta_D_target, mu_v, M0) {
  
  # Use the value of gamma under equal birth and death rates as initial guess
  guess <- baseline_psi(Delta_t, delta_D_target, M0)
  
  # Minimise the objective function
  fit <- optim(
    par = guess,
    fn = psi_objective_function,
    method = "L-BFGS-B",
    lower = 1e-6,
    upper = 100,
    control = list(factr = 1e2, pgtol = 1e-8, maxit = 1000),
    Delta_t = Delta_t,
    delta_D_target = delta_D_target,
    mu_v = mu_v,
    M0 = M0
  )
  
  # Return the estimate of psi
  fit$par
  
}

grid <- expand.grid(delta_t = params_base$delta_t_vec, 
                    delta_D = params_base$delta_D, 
                    M0 = params_base$M0, 
                    mu_v = params_base$mu_v)

#psi_vec <- numeric(length(params_base$delta_t_vec))  # store
psi_vec <- numeric(nrow(grid))

for (i in seq_len(nrow(grid))) {
  Delta_t <- grid$delta_t[i]
  Delta_D <- grid$delta_D[i]
  M0 <- grid$M0[i]
  mu_v <- grid$mu_v[i]
  
  if (constant_emergence == 0) {
    psi_fit <- fit_psi(Delta_t, Delta_D, mu_v, M0)
  } else {
    psi_fit <- baseline_psi(Delta_t, Delta_D, M0)
  }
  
  psi_vec[i] <- psi_fit  # store the result
}

grid$psi <- with(grid, {psi_vec})
#params_base$psi <- with(params_base, {psi_vec})


endemic_eqm <- function(params_base){
  with(params_base, {
    
    # Calculate infectious humans at equilibrium
    K_h <- 1 + (gamma_h + mu_h) / sigma_h + gamma_h / (omega_h + mu_h)
    K_hv <- (mu_h * K_h + (omega_h * gamma_h) / (omega_h + mu_h)) *
      (sigma_v + mu_v)
    
    numerator <- beta_hv * sigma_v * m0 - (K_hv * mu_v) / beta_vh
    denominator <- beta_hv * sigma_v * m0 * K_h + K_hv
    
    I_h_star <- N0 * numerator / denominator
    
    
    # Calculate other human variables
    E_h_star <- (gamma_h + mu_h) / sigma_h * I_h_star
    R_h_star <- gamma_h / (omega_h + mu_h) * I_h_star
    S_h_star <- N0 - (E_h_star + I_h_star + R_h_star)
    
    # Calculate mosquito variables
    S_v_star <- (mu_v * M0) / (mu_v + beta_vh * I_h_star / N0)
    E_v_star <- (beta_vh * S_v_star * I_h_star) / (N0 * (sigma_v + mu_v))
    I_v_star <- (sigma_v / mu_v) * E_v_star
    
    # Return vector of variables at the endemic equilibrium
    c(
      init_S_h = S_h_star, init_E_h = E_h_star, init_I_h = I_h_star, init_R_h = R_h_star, init_C = 0,
      init_S_v = S_v_star, init_E_v = E_v_star, init_I_v = I_v_star, init_D = 0
    )
    
  })
}



malaria_model <- odin::odin({
  #human model 
  deriv(S_h) <-  mu_h * N - beta_hv * S_h * I_v / N  + omega_h * R_h - mu_h * S_h
  deriv(E_h) <- beta_hv * S_h * I_v / N - sigma_h * E_h - mu_h * E_h
  deriv(I_h) <- sigma_h * E_h - gamma_h * I_h - mu_h * I_h
  deriv(R_h) <- gamma_h * I_h - omega_h * R_h - mu_h * R_h
  deriv(C) <-  beta_hv * S_h * I_v / N #this COLLECTS what ends up in Eh 
  
  #mosquito model 
  deriv(S_v) <- if (t >= tau_real && t <= (tau_real + delta_t)) phi - beta_vh * S_v * I_h/N - mu_v*S_v - psi*S_v else phi - beta_vh * S_v * I_h/N - mu_v*S_v
  deriv(E_v) <- if (t >= tau_real && t <= (tau_real + delta_t)) beta_vh * S_v * I_h / N - sigma_v * E_v - mu_v*E_v - psi*E_v else beta_vh * S_v * I_h / N - sigma_v * E_v - mu_v*E_v
  deriv(I_v) <- if (t >= tau_real && t <= (tau_real + delta_t)) sigma_v * E_v - mu_v*I_v - psi*I_v else sigma_v * E_v - mu_v*I_v 
  deriv(D) <-   if (t >= tau_real && t <= (tau_real + delta_t)) psi*S_v + psi*E_v + psi*I_v else 0
  
  
  
  init_S_h <- user()
  init_E_h <- user()
  init_I_h <- user()
  init_R_h <- user()
  init_C <- user()
  
  init_S_v <- user()
  init_E_v <- user()
  init_I_v <- user()
  init_D <- user()
  
  initial(S_h) <- init_S_h
  initial(E_h) <- init_E_h
  initial(I_h) <- init_I_h
  initial(R_h) <- init_R_h
  initial(C) <- init_C
  
  N <- S_h + E_h + I_h + R_h
  N0 <- user()
  
  initial(S_v) <- init_S_v
  initial(E_v) <- init_E_v
  initial(I_v) <- init_I_v
  initial(D) <- init_D
  
  M <- S_v + E_v + I_v
  m <- M/N
  s_h <- S_h/N
  s_v <- S_v/M
  M0 <- user()
  
  #reproduction number
  R0 <- (beta_hv*beta_vh*sigma_h*sigma_v)/((sigma_h + mu_h) * (gamma_h+mu_h) * (sigma_v + mu_v) * mu_v)
  R0_t <- R0*m #R0 at time t, given mosquito density at time t
  Re_t <- R0_t*s_h*s_v #effective reproduction number
  
  #rates
  mu_h <- user()
  beta_hv <- user()
  omega_h <- user()
  sigma_h <- user()
  gamma_h <- user()
  constant_emergence <- user()
  phi <- if(constant_emergence == 0) mu_v*M0 else mu_v*M
  beta_vh <- user()
  mu_v <- user()
  sigma_v <- user()
  psi <- user()
  delta_t <- user()
  n_times <- user()
  
  int_on[] <- user()
  dim(int_on) <- n_times
  
  tau_real <- interpolate(ttt, int_on, "constant")
  ttt[] <- user()
  dim(ttt) <- n_times
  
  #define outputs
  #output(prev) <- I_h/N
  output(R0_t) <- R0_t
  output(Re_t) <- Re_t
})


#mod <- malaria_model$new(user = params_base)


#param_grid <- tibble(
#  delta_t = params_base$delta_t_vec, 
#  psi = params_base$psi
#)

param_grid <- tibble(
  delta_D = grid$delta_D, 
  delta_t = grid$delta_t, 
  psi = grid$psi, 
  M0 = grid$M0, 
  mu_v = grid$mu_v
)

model_results <- vector("list", nrow(param_grid))

time_period <- 500
#psi_vec <- psi_vec

for (i in seq_len(nrow(param_grid))){
  in_delta_D <- param_grid$delta_D[i]
  in_delta_t <- param_grid$delta_t[i]
  in_psi <- param_grid$psi[i]
  in_M0 <- param_grid$M0[i]
  in_mu_v <- param_grid$mu_v[i]
  
  in_params <- params_base
  in_params$delta_t <- in_delta_t
  in_params$delta_D <- in_delta_D
  in_params$M0 <- in_M0
  in_params$mu_v <- in_mu_v
  #in_params$psi <- in_psi
  in_params$psi <- psi_vec[i]
  
  #update the mosquito-human ratio 
  in_params$m0 <- in_M0/in_params$N0
  
  #calc eqm for the parameter sets
  state <- endemic_eqm(in_params)
  in_params <- append(as.list(state), in_params)
  
  
  in_params$ttt <- seq(0, time_period, by = 1)
  in_params$n_times <- length(in_params$ttt)
  in_params$int_on <- rep(in_params$tau, in_params$n_times)
  
  #run model 
  mod <- malaria_model$new(user = in_params)
  print(glue::glue("Running delta_t = {in_delta_t}, psi = {in_psi}"))
  out <- mod$run(in_params$ttt, atol = 1e-6, rtol = 1e-6)
  
  #store output 
  model_results[[i]] <- cbind(as.data.frame(out), delta_t = in_delta_t, psi = in_psi, 
                              delta_D = in_delta_D, M0 = in_M0, m0 = in_params$m0, 
                              mu_v = in_mu_v)
}

model_results_df <- dplyr::bind_rows(model_results)
unique(model_results_df$delta_D)
unique(model_results_df$M0)
unique(model_results_df$mu_v)

#also run a baseline scenario: no interventions 

grid_baseline <- expand.grid(
                    M0 = params_base$M0, 
                    mu_v = params_base$mu_v)

param_grid_base <- tibble(
  delta_t = 0, 
  psi = 0, 
  delta_D = 0,
  M0 = grid_baseline$M0, 
  mu_v = grid_baseline$mu_v
)

model_results_base <- vector("list", nrow(param_grid_base))

time_period <- 500
#psi_vec <- psi_vec

for (i in seq_len(nrow(param_grid_base))){
  in_delta_D <- param_grid_base$delta_D[i]
  in_delta_t <- param_grid_base$delta_t[i]
  in_psi <- param_grid_base$psi[i]
  in_M0 <- param_grid_base$M0[i]
  in_mu_v <- param_grid_base$mu_v[i]
  
  in_params <- params_base
  in_params$delta_t <- in_delta_t
  in_params$delta_D <- in_delta_D
  in_params$M0 <- in_M0
  in_params$mu_v <- in_mu_v
  #in_params$psi <- in_psi
  in_params$psi <- psi_vec[i]
  
  #update the mosquito-human ratio 
  in_params$m0 <- in_M0/in_params$N0
  
  #calc eqm for the parameter sets
  state <- endemic_eqm(in_params)
  in_params <- append(as.list(state), in_params)
  
  
  in_params$ttt <- seq(0, time_period, by = 1)
  in_params$n_times <- length(in_params$ttt)
  in_params$int_on <- rep(in_params$tau, in_params$n_times)
  
  #run model 
  mod <- malaria_model$new(user = in_params)
  print(glue::glue("Running delta_t = {in_delta_t}, psi = {in_psi}"))
  out <- mod$run(in_params$ttt, atol = 1e-6, rtol = 1e-6)
  
  #store output 
  model_results_base[[i]] <- cbind(as.data.frame(out), delta_t = in_delta_t, psi = in_psi, 
                                   delta_D = in_delta_D, M0 = in_M0, m0 = in_params$m0, mu_v = in_params$mu_v)
}

model_results_base_df <- dplyr::bind_rows(model_results_base)
unique(model_results_base_df$m0)
model_results_base_df <- model_results_base_df %>%
  mutate(M = S_v + E_v + I_v, 
         N = S_h+E_h+I_h+R_h, 
         label =paste0("baseline: m0 = ", m0, 
                       ", mu_v = ",mu_v)) %>%
  rename(C0 = C)

unique(model_results_base_df$label)
unique(model_results_df$M0)

model_results_df <- model_results_df %>%
  mutate(N = S_h+E_h+I_h+R_h, 
         M = S_v+E_v+I_v, 
         label = paste0("psi = ", signif(psi, 3), 
                        ", delta_t = ", delta_t, 
                        ", D = ", delta_D, 
                        ", M0 = ", M0, 
                        ", m0=",m0, 
                        ", mu_v= ", mu_v))
unique(model_results_df$label)

model_results_df_all <- model_results_df %>%
  left_join(model_results_base_df[,c("t", "C0")], by = c("t")) %>%
  mutate(delta_C = C0-C, #baseline minus int
         rel_delta_C = ((C0-C)/C0)*100)

range(model_results_df_all$rel_delta_C, na.rm = TRUE) 
range(model_results_df_all$delta_C, na.rm = TRUE)

ggplot(model_results_base_df, aes(x = t, y = C0))+
  facet_wrap(vars(label))+
  geom_line()


ggplot(model_results_df, aes(x = t, y = I_h/N))+
  facet_wrap(vars(label))+
  geom_line()

ggplot(model_results_df, aes(x = t, y = I_v))+
  facet_wrap(vars(label))+
  geom_line()

ggplot(model_results_df, aes(x = t, y = D))+
  facet_wrap(vars(label))+
  geom_line()

ggplot(model_results_df, aes(x = t, y = C))+
  facet_wrap(vars(label))+
  geom_line()

ggplot(model_results_base_df, aes(x = t, y = I_h/N))+
  facet_wrap(vars(label))+
  geom_line()

model_results_base_df %>%
  group_by(m0, mu_v) %>%
  summarise(mean_prev = mean(I_h/N))



model_results_df_all <- model_results_df_all %>%
  mutate(prop_killed = delta_D/M0)
unique(model_results_df_all$prop_killed)

model_results_df <- model_results_df %>%
  mutate(constant_emergence = TRUE)

saveRDS(model_results_df, file = "2.output/model_results_df_constant_emergence_TRUE.rds")

model_results_base_df <- model_results_base_df %>%
  mutate(constant_emergence = TRUE)

saveRDS(model_results_base_df, file = "2.output/model_results_base_df_constant_emergence_TRUE.rds")




#plots for analysis####

#first look at killing 1000 out of 2000 mosquitoes (50%) at 2:1 vector to host ratio (2000 mosq to 1000 people)
death_rates <- unique(params_base$mu_v)
M0_vec <- unique(params_base$M0)
delta_D_vec <- unique(params_base$delta_D)
m0_vec <- unique(params_base$m0)

model_results_base_main <- model_results_base_df %>%
  filter(m0 == m0_vec & mu_v == mu_v_vec[1] & M0 == M0_vec[2])

model_results_main <- model_results_df %>%
  filter(delta_D == delta_D_vec[2] & M0 == M0_vec[2] & mu_v == mu_v_vec[1])

model_results_all_main <- model_results_main %>%
  left_join(model_results_base_main[,c("t", "C0", "mu_v", "M0")]) %>% # let it auto-select the columns to join by
  mutate(delta_C = C0-C, 
         rel_delta_C = ((C0-C)/C0)*100)

#take measurements of prevalence averted at day 10, day 30 and 90

mosq_killed_dynamics <- ggplot(model_results_all_main, aes(x = t, y = D, col = as.factor(label)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = "grey")+
  theme_bw()+
  guides(col = "none")+
  labs(col = "Scenario")+
  ylab("Number of mosquitoes killed by \n intervention")+
  ylim(0,2000)

prev_dynamics <- ggplot(model_results_all_main, aes(x = t, y = (I_h/N)*100, col = as.factor(label)))+
  geom_line(linewidth = 1.1)+
  #geom_line(data = model_results_base_df, aes(x = t, y = I_h/N, col = "baseline"), linetype = "dashed")+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = "grey")+
  theme_bw()+
  guides(col = "none")+
  ylab("Prevalence (%)")+
  ylim(0, 24)

delta_C_dynamics <- ggplot(model_results_all_main, aes(x = t, y = delta_C, col = as.factor(label)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = "grey")+
  theme_bw()+
  guides(col = "none")+
  ylab("Difference in number of cases \n compared to baseline")

mosq_dynamics <- ggplot(model_results_all_main, aes(x = t, y = M, col = as.factor(label)))+
  geom_line(linewidth = 1.1)+
  geom_vline(xintercept = 100, linetype = "dashed", linewidth = 1.1)+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[1], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[2], linetype = "dotted", linewidth = 1.1, col = "grey")+
  geom_vline(xintercept = 100 + params_base$delta_t_vec[3], linetype = "dotted", linewidth = 1.1, col = "grey")+
  theme_bw()+
  theme(legend.position = c(0.6, 0.3))+
  ylab("Mosquito population size")+
  ylim(0, 2000)

fig_1 <- cowplot::plot_grid(mosq_killed_dynamics, mosq_dynamics, 
                            prev_dynamics, delta_C_dynamics, 
                            labels = c("A", "B", "C", "D"))


#which has a greater epi impact?

ggplot(model_results_base_df, aes(x =t, y = I_h/N, col = as.factor(label)))+
  geom_line()

ggplot(model_results_df_all, aes(x =t, y = I_h/N, col = as.factor(label)))+
  geom_line()

#summary stats for baseline scenario (nothing changes over time, so don't need to filter by time here)

model_base_epi <- model_results_base_df %>%
  group_by(m0, mu_v, M0) %>%
  summarise(mean_prev_baseline = mean(I_h/N), 
          mean_R0_t_baseline = mean(R0_t), 
          mean_Re_t_baseline = mean(Re_t)) 


model_all_summary_d10 <- model_results_df_all %>%
  filter(between(t, params_base$tau, params_base$tau+params_base$delta_t_vec[1])) %>%
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t)) %>%
  mutate(measurement_t = params_base$delta_t_vec[1])

model_all_summary_d30 <- model_results_df_all %>%
  filter(between(t, params_base$tau, params_base$tau+params_base$delta_t_vec[2])) %>%
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t))%>%
  mutate(measurement_t = params_base$delta_t_vec[2])

model_all_summary_d90 <- model_results_df_all %>%
  filter(between(t, params_base$tau, params_base$tau+params_base$delta_t_vec[3])) %>%
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t))%>%
  mutate(measurement_t = params_base$delta_t_vec[3])

model_all_summary_d300 <- model_results_df_all %>%
  filter(between(t, params_base$tau, 300)) %>% #all back to eqm now
  group_by(delta_t, delta_D, m0, mu_v, M0) %>%
  summarise(
    mean_prev_int = mean(I_h/N),
    mean_R0_t_int = mean(R0_t), 
    mean_Re_t_int = mean(Re_t))%>%
  mutate(measurement_t = 300)


model_all_summary <- do.call("rbind", list(model_all_summary_d10, model_all_summary_d30, model_all_summary_d90, model_all_summary_d300))

#compare diff in prevalence
summary_impact <- left_join(model_all_summary, model_base_epi, by = c("m0", "M0", "mu_v")) %>%
  mutate(abs_diff_prev = mean_prev_baseline - mean_prev_int, 
         rel_diff_prev = ((mean_prev_baseline - mean_prev_int)/mean_prev_baseline)*100)

rel_diff_prev_plot <- summary_impact %>%
  filter(mu_v == 0.10 & delta_D == 1000 & m0 ==2 & M0 == 2000) %>%
  ggplot()+
  aes(x = as.factor(delta_t), y = rel_diff_prev, fill = as.factor(measurement_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  xlab("Duration of killing")+
  ylab("Absolute difference in prevalence compared to baseline")

abs_diff_prev_plot <- summary_impact %>%
  filter(mu_v == 0.10 & delta_D == 1000 & m0 ==2 & M0 == 2000) %>%
  ggplot()+
  aes(x = as.factor(delta_t), y = abs_diff_prev, fill = as.factor(measurement_t))+
  geom_bar(stat = "identity", position = position_dodge())+
  xlab("Duration of killing")+
  ylab("Relative difference in prevalence compared to baseline (%)")
  
