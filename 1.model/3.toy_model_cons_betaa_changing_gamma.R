#a simple toy model using odin to determine whether to kill "100 mosquitoes in a day or 1 mosquito a day for 100 days"
#constant betaa and changing gamma
require(tidyverse)
require(odin)

model1 <- odin::odin({
  #differential eqs
  deriv(M) <- if (t >= int_on_real && t < (int_on_real + dur)) betaa - mu*M - gamma*M else betaa - mu*M
  deriv(int_dead) <- if (t >= int_on_real &&  t < (int_on_real + dur)) gamma*M else 0
  deriv(int_natural_dead) <-  if (t >= int_on_real &&  t < (int_on_real + dur)) mu*M else 0
  
  #states
  initial(M) <- M0
  initial(int_dead) <- 0
  initial(int_natural_dead) <- 0
  
  M0 <- 1000
  
  #rate parameter values
  betaa <- M0*mu #constant emergence rate
  mu <- 0.132
  AG_out <-user() #this is the rates calculated frm Andrew's work 
  mosq_to_kill <- 100 #need to check unit on this, i.e. over what duration? talk to TC. Can see max 200 mosq killed by int if this is 100
  gammaT <- AG_out*mosq_to_kill #this needs checking against AG soln
  gamma <- gammaT/M 
  
  #define times 
  dur <- user() #this is dependent on the endectocide
  n_times <- user()
  
  int_on[] <- user()
  dim(int_on) <- n_times
  
  int_on_real <- interpolate(ttt, int_on, "constant")
  ttt[] <- user()
  dim(ttt) <- n_times
  
  #define outputs: can't have same as variable name
  output(M_out) <- M
  output(M0_out) <- M0
  output(int_dead_out) <- int_dead
  output(int_natural_dead_out) <- int_natural_dead
  output(betaa_out) <- betaa
  output(gamma_out) <- gamma
  output(mu_out) <- mu
  output(int_on_real_out) <- int_on_real
  output(dur_out) <- dur
})

time_period <- 365*3
step_size <- 1
ttt <- seq(0, time_period, step_size)
n_times <- length(ttt)

dur1 <- 23
dur2 <- 20
dur3 <- 26
dur4 <- 40


#check int_on matches this length 
int_on1 <- rep(time_period, n_times)
int_on1[ttt >= 300 & ttt < (300 + dur1)] <- rep(300,dur1)
endec1 <- model1$new(dur = dur1, n_times = n_times, ttt = ttt, int_on = int_on1, AG_out = 0.09)

#then run the model over the time period 
tt_endec1 <- seq(0, time_period, length.out = time_period)
run_endec1 <- endec1$run(tt_endec1)

df_endec1 <- as.data.frame(run_endec1) 
ggplot(df_endec1, aes(x = t, y = M_out))+
  geom_line()

ggplot(df_endec1, aes(x = t, y = int_dead_out))+
  geom_line()

gamma_estimator <- function(int_rateA, durationA, durationB){
  int_rateB <- (int_rateA*durationA)/durationB
  return(int_rateB)
}

endec1_int_rate = 0.09

#for endec2 
endec2_int_rate <- gamma_estimator(int_rateA = endec1_int_rate, durationA = dur1, 
                                   durationB = dur2)
endec2_int_rate #0.1035

#for endec3 
endec3_int_rate <- gamma_estimator(int_rateA =  endec1_int_rate, durationA = dur1, 
                                   durationB = dur3)
#0.07961538

#endec4 
endec4_int_rate <- gamma_estimator(int_rateA =  endec1_int_rate, durationA = dur1, 
                                   durationB = dur4)
#0.05175

#then run for the different products

#endec2###
int_on2 <- rep(time_period, n_times)
int_on2[ttt >= 300 & ttt < (300 + dur2)] <- rep(300,dur2)
endec2 <- model1$new(dur = dur2, n_times = n_times, ttt = ttt, int_on = int_on2, AG_out = endec2_int_rate)

#then run the model over the time period 
tt_endec2 <- seq(0, time_period, length.out = time_period)
run_endec2 <- endec2$run(tt_endec2)

df_endec2 <- as.data.frame(run_endec2)

#endec3###
int_on3 <- rep(time_period, n_times)
int_on3[ttt >= 300 & ttt < (300 + dur3)] <- rep(300,dur3)
endec3 <- model1$new(dur = dur3, n_times = n_times, ttt = ttt, int_on = int_on3, AG_out = endec3_int_rate)

#then run the model over the time period 
tt_endec3 <- seq(0, time_period, length.out = time_period)
run_endec3 <- endec3$run(tt_endec3)

df_endec3 <- as.data.frame(run_endec3)

#endec4###
int_on4 <- rep(time_period, n_times)
int_on4[ttt >= 300 & ttt < (300 + dur4)] <- rep(300,dur4)
endec4 <- model1$new(dur = dur4, n_times = n_times, ttt = ttt, int_on = int_on4, AG_out = endec4_int_rate)

#then run the model over the time period 
tt_endec4 <- seq(0, time_period, length.out = time_period)
run_endec4 <- endec4$run(tt_endec4)

df_endec4 <- as.data.frame(run_endec4)

#combine the datasets
all_products <- do.call("rbind", list(df_endec1, df_endec2, df_endec3, df_endec4))

mosq_dead <- ggplot(all_products, aes(x = t, y = int_dead_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Number of \n intervention-killed mosquitoes")+
  labs(col = "Intervention duration")+
  ggtitle("Constant betaa, varying gamma")+
  xlim(250, 420)

mosq_all <- ggplot(all_products, aes(x = t, y = M_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Total number of mosquitoes")+
  labs(col = "Intervention duration")+
  xlim(250, 420)

betaa_plot <- ggplot(all_products, aes(x = t, y = betaa_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Number of mosquitoes \n emerging per day")+
  labs(col = "Intervention duration")+
  xlim(250, 420)

cowplot::plot_grid(mosq_dead, mosq_all, betaa_plot, nrow = 3)


#baseline model 
endec0 <- model1$new(dur = dur1, n_times = n_times, ttt = ttt, int_on = int_on1, AG_out = 0)
tt_endec0 <- seq(0, time_period, length.out = time_period)
run_endec0 <- endec0$run(tt_endec0)
df_endec0 <- as.data.frame(run_endec0) %>%
  mutate(dur_out = 0)

ggplot(df_endec0, aes(x = t, y = int_dead_out))+
  geom_line()
ggplot(df_endec0, aes(x = t, y = M_out))+
  geom_line()
ggplot(df_endec0, aes(x = t, y = betaa_out))+
  geom_line()

#combine the datasets
all_products <- do.call("rbind", list(df_endec0, df_endec1, df_endec2, df_endec3, df_endec4))

mosq_dead <- ggplot(all_products, aes(x = t, y = int_dead_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Number of \n intervention-killed mosquitoes")+
  labs(col = "Intervention duration")+
  ggtitle("Constant betaa, varying gamma")+
  xlim(290, 420)

mosq_all <- ggplot(all_products, aes(x = t, y = M_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Total number of mosquitoes")+
  labs(col = "Intervention duration")+
  xlim(290, 420)

betaa_plot <- ggplot(all_products, aes(x = t, y = betaa_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Number of mosquitoes \n emerging per day")+
  labs(col = "Intervention duration")+
  xlim(290, 420)

death_rate_plot <- ggplot(all_products, aes(x = t, y = gamma_out, col = as.factor(dur_out)))+
  geom_line(size = 1.5, alpha = 0.6)+
  theme_bw()+
  ylab("Int death rate")+
  labs(col = "Intervention duration")+
  xlim(290, 420)

cowplot::plot_grid(mosq_dead, mosq_all, betaa_plot,death_rate_plot, nrow = 4)

#measure impact in a i) 4-month period from int turned on ii) for duration of longest product and
#iii) duration of each product 

mean_M_out_4_months <- all_products %>%
  group_by(dur_out) %>%
  filter(between(t, 300, 300+(30*4))) %>%
  summarise(mean_M_out = mean(M_out))

mean_M_out_4_months_wide <- mean_M_out_4_months %>%
  pivot_wider(names_from = dur_out, values_from = mean_M_out)

impact_endec1_4m <- ((mean_M_out_4_months_wide$`0` - mean_M_out_4_months_wide$`23`)/mean_M_out_4_months_wide$`0`) *100
impact_endec3_4m <- ((mean_M_out_4_months_wide$`0` - mean_M_out_4_months_wide$`26`)/mean_M_out_4_months_wide$`0`) *100
impact_endec2_4m <- ((mean_M_out_4_months_wide$`0` - mean_M_out_4_months_wide$`20`)/mean_M_out_4_months_wide$`0`) *100
impact_endec4_4m <- ((mean_M_out_4_months_wide$`0` - mean_M_out_4_months_wide$`40`)/mean_M_out_4_months_wide$`0`) *100


measurement_method <- c("4 months", "Dur of longest product (40 days)", "for each product's killing duration")
results_df <- data.frame(measurement_method)

results_df <- results_df %>%
  mutate(
    prod_endec1 = case_when(
      measurement_method == "4 months" ~ impact_endec1_4m,
      TRUE ~ NA
    ),
    prod_endec2 = case_when(
      measurement_method == "4 months" ~ impact_endec2_4m,
      TRUE ~ NA
    ),
    prod_endec3 = case_when(
      measurement_method == "4 months" ~ impact_endec3_4m,
      TRUE ~ NA
    ),
    prod_endec4 = case_when(
      measurement_method == "4 months" ~ impact_endec4_4m,
      TRUE ~ NA
    )
  )


mean_M_out_dur4 <- all_products %>%
  group_by(dur_out) %>%
  filter(between(t, 300, 300+(dur4))) %>%
  summarise(mean_M_out = mean(M_out))

mean_M_out_dur4_wide <- mean_M_out_dur4 %>%
  pivot_wider(names_from = dur_out, values_from = mean_M_out)

impact_endec1_dur4 <- ((mean_M_out_dur4_wide$`0` - mean_M_out_dur4_wide$`23`)/mean_M_out_dur4_wide$`0` )*100
impact_endec3_dur4 <- ((mean_M_out_dur4_wide$`0` - mean_M_out_dur4_wide$`26`)/mean_M_out_dur4_wide$`0` )*100
impact_endec2_dur4 <- ((mean_M_out_dur4_wide$`0` - mean_M_out_dur4_wide$`20`)/mean_M_out_dur4_wide$`0` )*100
impact_endec4_dur4 <- ((mean_M_out_dur4_wide$`0` - mean_M_out_dur4_wide$`40`)/mean_M_out_dur4_wide$`0` )*100

results_df <- results_df %>%
  mutate(
    prod_endec1 = case_when(
      measurement_method == measurement_method[2] ~ impact_endec1_dur4,
      TRUE ~ prod_endec1
    ),
    prod_endec2 = case_when(
      measurement_method == measurement_method[2] ~ impact_endec2_dur4,
      TRUE ~ prod_endec2
    ),
    prod_endec3 = case_when(
      measurement_method == measurement_method[2] ~ impact_endec3_dur4,
      TRUE ~ prod_endec3
    ),
    prod_endec4 = case_when(
      measurement_method == measurement_method[2] ~ impact_endec4_dur4,
      TRUE ~ prod_endec4
    )
  )

#for dur of each product

endec0_dur1 <- df_endec0 %>%
  filter(between(t, 300, 300+dur1)) %>%
  summarise(mean_M_out = mean(M_out))

endec0_dur2 <- df_endec0 %>%
  filter(between(t, 300, 300+dur2)) %>%
  summarise(mean_M_out = mean(M_out))

endec0_dur3 <- df_endec0 %>%
  filter(between(t, 300, 300+dur3)) %>%
  summarise(mean_M_out = mean(M_out))

endec0_dur4 <- df_endec0 %>%
  filter(between(t, 300, 300+dur4)) %>%
  summarise(mean_M_out = mean(M_out))

avert_endec1 <- df_endec1 %>%
  filter(between(t, 300, 300+dur1)) %>%
  summarise(mean_M_out = mean(M_out))
avert_endec2 <- df_endec2 %>%
  filter(between(t, 300, 300+dur2)) %>%
  summarise(mean_M_out = mean(M_out))
avert_endec3 <- df_endec3 %>%
  filter(between(t, 300, 300+dur3)) %>%
  summarise(mean_M_out = mean(M_out))
avert_endec4 <- df_endec4 %>%
  filter(between(t, 300, 300+dur4)) %>%
  summarise(mean_M_out = mean(M_out))

impact_endec1_dur1 <- ((endec0_dur1-avert_endec1)/endec0_dur1)*100
impact_endec2_dur2 <- ((endec0_dur2-avert_endec2)/endec0_dur2)*100
impact_endec3_dur3 <- ((endec0_dur3-avert_endec3)/endec0_dur3)*100
impact_endec4_dur4 <- ((endec0_dur4-avert_endec4)/endec0_dur4)*100

results_df <- results_df %>%
  mutate(
    prod_endec1 = case_when(
      measurement_method == measurement_method[3] ~ impact_endec1_dur1[[1]],
      TRUE ~ prod_endec1
    ),
    prod_endec2 = case_when(
      measurement_method == measurement_method[3] ~ impact_endec2_dur2[[1]],
      TRUE ~ prod_endec2
    ),
    prod_endec3 = case_when(
      measurement_method == measurement_method[3] ~ impact_endec3_dur3[[1]],
      TRUE ~ prod_endec3
    ),
    prod_endec4 = case_when(
      measurement_method == measurement_method[3] ~ impact_endec4_dur4[[1]],
      TRUE ~ prod_endec4
    )
  )

results_df_method3 <- results_df %>%
  rename(Killing_duration_23d = prod_endec1, 
         Killing_duration_20d = prod_endec2, 
         Killing_duration_26d = prod_endec3, 
         Killing_duration_40d = prod_endec4,
         percent_M_averted_over_time = measurement_method)
