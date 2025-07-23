library(odin)
library(tidyverse)

# ------------------------------------------------
# 1. Flag for emergence type (0 = constant, 1 = logistic)
# ------------------------------------------------
constant_emergence <- 1

# ------------------------------------------------
# 2. Base parameters (named list)
# ------------------------------------------------
params_base <- list(
  N0 = 1000,                # Initial human population
  M0 = 2000,                # Initial mosquito population
  mu_h = 1 / (70 * 365),    # Human death rate
  sigma_h = 1 / 12,         # Human incubation rate
  gamma_h = 1 / 60,         # Human recovery rate
  omega_h = 1 / 180,        # Loss of immunity rate
  beta_hv = 0.3,            # Transmission prob. human to vector
  beta_vh = 0.3,            # Transmission prob. vector to human
  mu_v = 0.1,               # Mosquito death rate
  sigma_v = 1 / 10,         # Mosquito incubation rate (EIP)
  tau = 100,                  # Intervention start time (for fitting psi)
  delta_D = 1900,           # Target mosquitoes killed
  constant_emergence = constant_emergence,
  psi = 0.1                 # Initial guess for psi
)

# Mosquito-to-human ratio
params_base$m0 <- params_base$M0 / params_base$N0

# ------------------------------------------------
# 3. Endemic equilibrium for initial conditions
# ------------------------------------------------
endemic_eqm <- function(params) {
  with(params, {
    K_h <- 1 + (gamma_h + mu_h) / sigma_h + gamma_h / (omega_h + mu_h)
    K_hv <- (mu_h * K_h + (omega_h * gamma_h) / (omega_h + mu_h)) *
      (sigma_v + mu_v)
    numerator <- beta_hv * sigma_v * m0 - (K_hv * mu_v) / beta_vh
    denominator <- beta_hv * sigma_v * m0 * K_h + K_hv
    I_h_star <- N0 * numerator / denominator
    E_h_star <- (gamma_h + mu_h) / sigma_h * I_h_star
    R_h_star <- gamma_h / (omega_h + mu_h) * I_h_star
    S_h_star <- N0 - (E_h_star + I_h_star + R_h_star)
    S_v_star <- (mu_v * M0) / (mu_v + beta_vh * I_h_star / N0)
    E_v_star <- (beta_vh * S_v_star * I_h_star) / (N0 * (sigma_v + mu_v))
    I_v_star <- (sigma_v / mu_v) * E_v_star
    
    c(
      init_S_h = S_h_star, init_E_h = E_h_star, init_I_h = I_h_star,
      init_R_h = R_h_star, init_C = 0,
      init_S_v = S_v_star, init_E_v = E_v_star, init_I_v = I_v_star,
      init_D = 0
    )
  })
}

# ------------------------------------------------
# 4. Odin model
# ------------------------------------------------
malaria_model <- odin::odin({
  # Human
  deriv(S_h) <- mu_h * N - beta_hv * S_h * I_v / N + omega_h * R_h - mu_h * S_h
  deriv(E_h) <- beta_hv * S_h * I_v / N - sigma_h * E_h - mu_h * E_h
  deriv(I_h) <- sigma_h * E_h - gamma_h * I_h - mu_h * I_h
  deriv(R_h) <- gamma_h * I_h - omega_h * R_h - mu_h * R_h
  deriv(C)   <- beta_hv * S_h * I_v / N
  
  # Mosquito
  deriv(S_v) <- if (t >= tau_real && t <= (tau_real + delta_t))
    phi - beta_vh * S_v * I_h/N - mu_v*S_v - psi*S_v
  else
    phi - beta_vh * S_v * I_h/N - mu_v*S_v
  deriv(E_v) <- if (t >= tau_real && t <= (tau_real + delta_t))
    beta_vh * S_v * I_h/N - sigma_v*E_v - mu_v*E_v - psi*E_v
  else
    beta_vh * S_v * I_h/N - sigma_v*E_v - mu_v*E_v
  deriv(I_v) <- if (t >= tau_real && t <= (tau_real + delta_t))
    sigma_v*E_v - mu_v*I_v - psi*I_v
  else
    sigma_v*E_v - mu_v*I_v
  deriv(D)   <- if (t >= tau_real && t <= (tau_real + delta_t))
    psi*(S_v + E_v + I_v)
  else 0
  
  # Initial conditions
  init_S_h <- user(); initial(S_h) <- init_S_h
  init_E_h <- user(); initial(E_h) <- init_E_h
  init_I_h <- user(); initial(I_h) <- init_I_h
  init_R_h <- user(); initial(R_h) <- init_R_h
  init_C   <- user(); initial(C)   <- init_C
  
  init_S_v <- user(); initial(S_v) <- init_S_v
  init_E_v <- user(); initial(E_v) <- init_E_v
  init_I_v <- user(); initial(I_v) <- init_I_v
  init_D   <- user(); initial(D)   <- init_D
  
  N  <- S_h + E_h + I_h + R_h
  M  <- S_v + E_v + I_v
  m  <- M/N
  s_h <- S_h/N
  s_v <- S_v/M
  
  # Parameters
  N0 <- user()
  M0 <- user()
  mu_h <- user()
  beta_hv <- user()
  omega_h <- user()
  sigma_h <- user()
  gamma_h <- user()
  constant_emergence <- user()
  phi <- if(constant_emergence == 0) mu_v*M0 else mu_v*M0 + r*M*(1-M/M0)
  r <- mu_v
  beta_vh <- user()
  mu_v <- user()
  sigma_v <- user()
  psi <- user()
  delta_t <- user()
  n_times <- user()
  
  int_on[] <- user()
  dim(int_on) <- n_times
  
  ttt[] <- user()
  dim(ttt) <- n_times
  tau_real <- interpolate(ttt, int_on, "constant")
  
  # reproduction number
  R0 <- (beta_hv*beta_vh*sigma_h*sigma_v)/((sigma_h + mu_h) * (gamma_h+mu_h) * (sigma_v + mu_v) * mu_v)
  R0_t <- R0*m #R0 at time t, given mosquito density at time t
  Re_t <- R0_t*s_h*s_v #effective reproduction number
  
  #outputs
  output(R0_t) <- R0_t
  output(Re_t) <- Re_t
  
})

# ------------------------------------------------
# 5. psi fitting (logistic emergence)
# ------------------------------------------------
psi_objective_logistic <- function(psi, Delta_t, delta_D_target, params_base) {
  # Copy and update parameters
  params <- params_base
  params$psi <- psi
  params$delta_t <- Delta_t
  params$n_times <- 1
  params$ttt <- 0
  params$int_on <- 0
  
  # Add initial conditions (required user params)
  eqm <- endemic_eqm(params)
  params <- c(params, as.list(eqm))
  
  # Now create the model
  mod <- malaria_model$new(user = params)
  out <- mod$run(seq(0, Delta_t, by = 0.1))
  
  D_pred <- tail(out[,"D"], 1)
  (delta_D_target - D_pred)^2
}


fit_psi_logistic <- function(Delta_t, delta_D_target, params_base) {
  optim(
    par = 0.1,
    fn = psi_objective_logistic,
    method = "L-BFGS-B",
    lower = 1e-6,
    upper = 10,
    Delta_t = Delta_t,
    delta_D_target = delta_D_target,
    params_base = params_base
  )$par
}

#test block####

test_Delta_t <- 30       # 30-day intervention
test_delta_D <- 1900     # target mosquito deaths
test_params <- params_base
test_params$M0 <- 2000   # test mosquito population
test_params$m0 <- test_params$M0 / test_params$N0

psi_test <- fit_psi_logistic(
  Delta_t = test_Delta_t,
  delta_D_target = test_delta_D,
  params_base = test_params
)

cat(sprintf("Test run: Delta_t=%d days, delta_D=%d, psi=%.4f\n",
            test_Delta_t, test_delta_D, psi_test))

#####

# ------------------------------------------------
# 6. Parameter grid and loop
# ------------------------------------------------
grid <- expand.grid(
  delta_t = c(10, 30, 90),
  delta_D = c(190000, 1900),
  M0 = c(200000, 2000),
  mu_v = params_base$mu_v, 
  N0 = params_base$N0, 
  tau = params_base$tau
)

write_rds(grid,"2.output_params_base_logistic.rds")

psi_vec <- numeric(nrow(grid))

for (i in seq_len(nrow(grid))) {
  Delta_t <- grid$delta_t[i]
  Delta_D <- grid$delta_D[i]
  M0 <- grid$M0[i]
  mu_v <- grid$mu_v[i]
  
  params_iter <- params_base
  params_iter$M0 <- M0
  params_iter$mu_v <- mu_v
  params_iter$m0 <- M0 / params_iter$N0
  
  psi_fit <- fit_psi_logistic(Delta_t, Delta_D, params_iter)
  psi_vec[i] <- psi_fit
  cat(sprintf("Run %d: delta_t=%.1f, delta_D=%.1f, M0=%.0f, psi=%.4f\n",
              i, Delta_t, Delta_D, M0, psi_fit))
}

# ------------------------------------------------
# 7. Pass into the model 
# ------------------------------------------------

param_grid <- tibble(
  delta_D = grid$delta_D, 
  delta_t = grid$delta_t, 
  psi = psi_vec,            # use fitted psi values
  M0 = grid$M0, 
  mu_v = grid$mu_v
)

model_results <- vector("list", nrow(param_grid))

time_period <- 500  # total simulation time

for (i in seq_len(nrow(param_grid))) {
  in_delta_D <- param_grid$delta_D[i]
  in_delta_t <- param_grid$delta_t[i]
  in_psi <- param_grid$psi[i]
  in_M0 <- param_grid$M0[i]
  in_mu_v <- param_grid$mu_v[i]
  
  # Base parameter list
  in_params <- params_base
  in_params$delta_t <- in_delta_t
  in_params$delta_D <- in_delta_D
  in_params$M0 <- in_M0
  in_params$mu_v <- in_mu_v
  in_params$psi <- in_psi
  
  # Update mosquito-human ratio
  in_params$m0 <- in_M0 / in_params$N0
  
  # Compute endemic equilibrium state
  state <- endemic_eqm(in_params)
  in_params <- c(as.list(state), in_params)
  
  # Time-related inputs
  in_params$ttt <- seq(0, time_period, by = 1)
  in_params$n_times <- length(in_params$ttt)
  in_params$int_on <- rep(in_params$tau, in_params$n_times)
  
  # Ensure all required parameters are present
  required_params <- c(
    "init_S_h","init_E_h","init_I_h","init_R_h","init_C",
    "init_S_v","init_E_v","init_I_v","init_D",
    "N0","M0","mu_h","beta_hv","omega_h","sigma_h","gamma_h",
    "constant_emergence","beta_vh","mu_v","sigma_v","psi","delta_t",
    "n_times","int_on","ttt"
  )
  
  missing <- setdiff(required_params, names(in_params))
  if (length(missing) > 0) {
    stop(sprintf("Missing parameters: %s", paste(missing, collapse = ", ")))
  }
  
  # Run the model
  mod <- malaria_model$new(user = in_params)
  print(glue::glue("Running delta_t = {in_delta_t}, psi = {in_psi}"))
  
  out <- mod$run(in_params$ttt, atol = 1e-6, rtol = 1e-6)
  
  # Store output
  model_results[[i]] <- cbind(
    as.data.frame(out),
    delta_t = in_delta_t, psi = in_psi, 
    delta_D = in_delta_D, M0 = in_M0, m0 = in_params$m0, 
    mu_v = in_mu_v
  )
}

model_results_df <- dplyr::bind_rows(model_results)

#save output
model_results_df <- model_results_df %>%
  mutate(N = S_h+E_h+I_h+R_h, 
         M = S_v+E_v+I_v,
         constant_emergence = FALSE) #we use the logistic growth model 

saveRDS(model_results_df, file = "2.output/model_results_df_logistic_growth.rds")

#some plots to check all is well#

model_results_main <- model_results_df %>%
  filter(delta_D == 1900 &  M0 == 2000)

model_results_main %>%
  group_by(delta_t) %>%
  summarise(max_D = max(D))

ggplot(model_results_main, aes(x = t, y = D, col = as.factor(delta_t)))+
  geom_line()

ggplot(model_results_main, aes(x = t, y = M, col = as.factor(delta_t)))+
  geom_line()+
  ylim(0,2000)

ggplot(model_results_main, aes(x = t, y = C, col = as.factor(delta_t)))+
  geom_line()+
  ylim(0,2000)

ggplot(model_results_main, aes(x = t, y = I_h/N, col = as.factor(delta_t)))+
  geom_line()+
  ylim(0, 0.3)

ggplot(model_results_main, aes(x = t, y = I_v/M, col = as.factor(delta_t)))+
  geom_line()+
  ylim(0, 0.3)

#run the baseline scenario

# ------------------------------------------------
# 8. Pass into the model for baseline scenario
# ------------------------------------------------

grid_baseline <- expand.grid(M0 = params_base$M0, 
                             mu_v = params_base$mu_v)

param_grid_base <- tibble(
  delta_D = 0, 
  delta_t = 0, 
  psi =0,            # no intervention
  M0 = grid_baseline$M0, 
  mu_v = grid_baseline$mu_v
)

model_results_base <- vector("list", nrow(param_grid_base))

time_period <- 500  # total simulation time

for (i in seq_len(nrow(param_grid_base))) {
  in_delta_D <- param_grid_base$delta_D[i]
  in_delta_t <- param_grid_base$delta_t[i]
  in_psi <- param_grid_base$psi[i]
  in_M0 <- param_grid_base$M0[i]
  in_mu_v <- param_grid_base$mu_v[i]
  
  # Base parameter list
  in_params <- params_base
  in_params$delta_t <- in_delta_t
  in_params$delta_D <- in_delta_D
  in_params$M0 <- in_M0
  in_params$mu_v <- in_mu_v
  in_params$psi <- in_psi
  
  # Update mosquito-human ratio
  in_params$m0 <- in_M0 / in_params$N0
  
  # Compute endemic equilibrium state
  state <- endemic_eqm(in_params)
  in_params <- c(as.list(state), in_params)
  
  # Time-related inputs
  in_params$ttt <- seq(0, time_period, by = 1)
  in_params$n_times <- length(in_params$ttt)
  in_params$int_on <- rep(in_params$tau, in_params$n_times)
  
  # Ensure all required parameters are present
  required_params <- c(
    "init_S_h","init_E_h","init_I_h","init_R_h","init_C",
    "init_S_v","init_E_v","init_I_v","init_D",
    "N0","M0","mu_h","beta_hv","omega_h","sigma_h","gamma_h",
    "constant_emergence","beta_vh","mu_v","sigma_v","psi","delta_t",
    "n_times","int_on","ttt"
  )
  
  missing <- setdiff(required_params, names(in_params))
  if (length(missing) > 0) {
    stop(sprintf("Missing parameters: %s", paste(missing, collapse = ", ")))
  }
  
  # Run the model
  mod <- malaria_model$new(user = in_params)
  print(glue::glue("Running delta_t = {in_delta_t}, psi = {in_psi}"))
  
  out <- mod$run(in_params$ttt, atol = 1e-6, rtol = 1e-6)
  
  # Store output
  model_results_base[[i]] <- cbind(
    as.data.frame(out),
    delta_t = in_delta_t, psi = in_psi, 
    delta_D = in_delta_D, M0 = in_M0, m0 = in_params$m0, 
    mu_v = in_mu_v
  )
}

model_results_base_df <- dplyr::bind_rows(model_results_base)

#save output
model_results_base_df <- model_results_base_df %>%
  rename(C0 = C)  %>%
  mutate(N = S_h+E_h+I_h+R_h, 
         M = S_v+E_v+I_v,
         constant_emergence = FALSE) #we use the logistic growth model 

saveRDS(model_results_base_df, file = "2.output/model_results_base_df_logistic_growth.rds")

