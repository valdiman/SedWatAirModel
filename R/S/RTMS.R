# Code to model PCB 19 in laboratory experiments
# using sediment from Altavista, VI. Passive measurements
# of PCB 19 in the water and the air phases are predicted and
# linked to the water and air concentrations from the passive
# samplers.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
  library(gridExtra)
}

# Read data ---------------------------------------------------------------
{
  obs.data <- read.csv("Data/AVL_S_data_long.csv", check.names = FALSE)
  # Physical chemical properties
  pcp.data <- read.csv("Data/AVL_S_PCP.csv")
}

# Extract individual PCBi -------------------------------------------------
{
  pcb.name <- "PCB19"
  obs.data.pcbi <- obs.data[, c("sampler", "time", pcb.name)]
  pcp.data.pcbi <- pcp.data[pcp.data$congener == pcb.name, ]
}

# Organize data -----------------------------------------------------------
{
  spme <- obs.data.pcbi %>%
    filter(sampler == "SPME")
  puf <- obs.data.pcbi %>%
    filter(sampler == "PUF") %>%
    mutate(across(starts_with("PCB"),  ~ . * 10))
  sed <- obs.data.pcbi %>%
    filter(sampler == "sed")
}

# Reactive transport function ---------------------------------------------
rtm.PCB19 = function(t, state, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- pcp.data.pcbi$MW # g/mol PCB 19 molecular weight
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- pcp.data.pcbi$Kaw # PCB 19 dimensionless Henry's law constant @ 25 C
  dUaw <- pcp.data.pcbi$dUaw # internal energy for the transfer of air-water for PCB 19 (J/mol)
  Kaw.t <- Kaw*exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  Kow <- pcp.data.pcbi$Kow # PCB 19 octanol-water equilibrium partition coefficient (low value!!)
  dUow <-  pcp.data.pcbi$dUow # internal energy for the transfer of octanol-water for PCB 19 (J/mol)
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Koa <- pcp.data.pcbi$Koa # PCB 19 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Apuf <- 7.07 # cm2
  Vpuf <- 29 # cm3 volume of PUF
  d <- 0.0213*100^3 # g/m3 density of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) # PCB 19-PUF equilibrium partition coefficient [m3/g]
  Kpuf <- Kpuf * d # [La/Lpuf]
  
  # SPME fiber constants
  Aspme <- 0.138 # cm2/cm SPME area
  Vspme <- 0.000000069 * 1000 # cm3/cm SPME volume/area
  Lspme <- 1 # cm SPME length normalization to 1 cm
  Kspme <- 10^(1.06 * log10(Kow.t) - 1.16) # PCB 19-SPME equilibrium partition coefficient
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 19's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 19's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 19
  
  # kaw calculations (air-water mass transfer coefficient)
  # i) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 19. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 19
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 19, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioremediation rate
  kb <- parms$kb
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  
  # Passive sampler rates
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  ro <- parms$ro # cm/d sampling rate for PUF
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  Cspme <- state[3]
  Ca <- state[4]
  Cpuf <- state[5]
  
  dCsdt <- - f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw
  dCwdt <- - ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) - 
    ko * Aspme * Lspme / Vw * (Cw - Cspme / Kspme) -
    kb * Cw # [ng/L]
  dCspmedt <- ko * Aspme / Vspme * (Cw - Cspme / Kspme) # Cw = [ng/L], Cspme = [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf) # Ca = [ng/L]
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf) # Ca = [ng/L], Cpuf = [ng/L]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dCspmedt, dCadt, dCpufdt)))
}

# Initial conditions and run function
{
  # Estimating Cs0
  Ct <- sed$PCB19 # ng/g PCB 19 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = Cs0, Cw = 0, Cspme = 0, Ca = 0, Cpuf = 0)
parms <- list(ro = 540.409, ko = 10, kdf = 1.9, kds = 0.001, f = 0.8,
              ka = 150, kb = 0) # Input
t.1 <- unique(spme$time)
# Run the ODE function without specifying parms
model.result <- ode(y = cinit, times = t.1, func = rtm.PCB19, parms = parms)
head(model.result)

{
  # Transform Cf and Cpuf to mass/cm and mass/puf
  model.result <- as.data.frame(model.result)
  
  # Calculate Mf and Mpuf based on volumes
  ms <- 10 # [g]
  M <- 0.1 # kg/L solid-water ratio
  Vw <- 100 # [cm3]
  Va <- 125 # [cm3]
  Vspme <- 0.000000069 # L/cm SPME
  Vpuf <- 29 # cm3 volume of PUF
  model.result$msed <- model.result$Cs * ms / (M * 1000)
  model.result$mW <- model.result$Cw * Vw / 1000
  model.result$mspme <- model.result$Cspme * Vspme  # [ng/cm]
  model.result$mA <-  model.result$Ca * Va / 1000
  model.result$mpuf <- model.result$Cpuf * Vpuf / 1000 # [ng/puf]
  model.result$Mt <- model.result$msed + model.result$mW + model.result$Cspme * Vspme +
    model.result$mA + model.result$Cpuf * Vpuf / 1000
  
  # Filter out SPME and PUF rows
  spme <- obs.data.pcbi[obs.data.pcbi$sampler == "SPME", c("time", "PCB19")]
  puf  <- obs.data.pcbi[obs.data.pcbi$sampler == "PUF",  c("time", "PCB19")]
  # Combine them by row
  obs.data.pcbi.2 <- data.frame(
    time = spme$time,
    SPME = spme$PCB19,
    PUF  = puf$PCB19 * 10
  )
  
  # Convert model results to tibble and select relevant columns
  model_results <- as_tibble(model.result) %>%
    mutate(time = as.numeric(time)) %>%
    select(time, mspme, mpuf)
  
  # Merge model results with observed data
  comparison_data <- model_results %>%
    left_join(obs.data.pcbi.2, by = "time")
  
  # Calculate the averages of mf and mpuf within each group (e.g., per time)
  grouped_comparison <- comparison_data %>%
    group_by(time) %>%  # Adjust the grouping variable if needed
    summarise(
      avg_mspme_model = mean(mspme, na.rm = TRUE),
      avg_mspme_observed = mean(SPME, na.rm = TRUE),
      avg_mpuf_model = mean(mpuf, na.rm = TRUE),
      avg_mpuf_observed = mean(PUF, na.rm = TRUE)
    )
  
  # Define function to calculate R-squared, handling NA values
  m_r2 <- function(predicted, observed) {
    # Remove NA values from both predicted and observed
    valid_indices <- complete.cases(predicted, observed)
    predicted <- predicted[valid_indices]
    observed <- observed[valid_indices]
    
    # Calculate R-squared
    ss_total <- sum((observed - mean(observed))^2)
    ss_residual <- sum((observed - predicted)^2)
    1 - (ss_residual / ss_total)
  }
  
  # Calculate R-squared values based on grouped average data
  mspme_r2_value <- m_r2(grouped_comparison$avg_mspme_model, grouped_comparison$avg_mspme_observed)
  mpuf_r2_value <- m_r2(grouped_comparison$avg_mpuf_model, grouped_comparison$avg_mpuf_observed)
  
  # Print R-squared values
  print(paste("R-squared for mspme (average): ", mspme_r2_value))
  print(paste("R-squared for mpuf (average): ", mpuf_r2_value))
  
  # Plot
  # Run the model with the new time sequence
  cinit <- c(Cs = Cs0, Cw = 0, Cspme = 0, Ca = 0, Cpuf = 0)
  t_daily <- seq(0, 40, by = 1)  # Adjust according to your needs
  out_daily <- ode(y = cinit, times = t_daily, func = rtm.PCB19, parms = parms)
  head(out_daily)
  
  # Transform Cf and Cpuf to mass/cm and mass/puf
  out.daily <- as.data.frame(out_daily)
  colnames(out.daily) <- c("time", "Cs", "Cw", "Cspme", "Ca", "Cpuf")
  
  # Calculate Mf and Mpuf based on volumes
  out.daily$mspme <- out.daily$Cspme * Vspme # [ng]
  out.daily$mpuf <- out.daily$Cpuf * Vpuf / 1000  # [ng] Check this 10!!
  
  # Convert model results to tibble and ensure numeric values
  model_results_daily_clean <- as_tibble(out.daily) %>%
    mutate(across(c(mspme, mpuf, Cw, time), as.numeric)) %>%  # Ensure all relevant columns are numeric
    select(time, mspme, mpuf)
  
  # Export data
  #write.csv(model_results_daily_clean, file = "Output/Data/RTM/S/AVL/PCB19AVLSControlFV.csv")
  
  # Prepare model data for plotting
  model_data_long <- model_results_daily_clean %>%
    pivot_longer(cols = c(mspme, mpuf), 
                 names_to = "variable", 
                 values_to = "model_value") %>%
    mutate(type = "Model")
  
  # Clean observed data and prepare for plotting
  observed_data_clean <- obs.data.pcbi.2 %>%
    pivot_longer(cols = c(SPME, PUF), 
                 names_to = "variable", 
                 values_to = "observed_value") %>%
    mutate(variable = recode(variable, 
                             "SPME" = "mspme", 
                             "PUF" = "mpuf"),
           type = "Observed")
  
  plot_data_daily <- bind_rows(
    model_data_long %>%
      rename(value = model_value) %>%
      mutate(type = "Model"),
    observed_data_clean %>%
      rename(value = observed_value) %>%
      mutate(type = "Observed")
  )
  
  # Plot mf
  p_spme <- ggplot(plot_data_daily %>% filter(variable == "mspme"), aes(x = time)) +
    geom_line(data = . %>% filter(type == "Model"),
              aes(y = value, color = "Model"), linewidth = 1) +
    geom_point(data = . %>% filter(type == "Observed"),
               aes(y = value, color = "Observed"), size = 2) +
    labs(x = "Time", y = "mf [ng/cm]") +
    scale_color_manual(values = c("Model" = "blue", "Observed" = "red")) +
    theme_bw() +
    theme(legend.title = element_blank())
  
  # Plot mpuf
  p_puf <- ggplot(plot_data_daily %>% filter(variable == "mpuf"), aes(x = time)) +
    geom_line(data = . %>% filter(type == "Model"),
              aes(y = value, color = "Model"), linewidth = 1) +
    geom_point(data = . %>% filter(type == "Observed"),
               aes(y = value, color = "Observed"), size = 2) +
    labs(x = "Time", y = "mpuf [ng/PUF]") +
    scale_color_manual(values = c("Model" = "blue", "Observed" = "red")) +
    theme_bw() +
    theme(legend.title = element_blank())
  
}

# Arrange plots side by side
p.19 <- grid.arrange(p_spme, p_puf, ncol = 2)

# Save plot in folder
#ggsave("Output/Plots/RTM/S/AVL/PCB19ALV_S_ControlFV.png", plot = p.19, width = 15,
#       height = 5, dpi = 500)


# Recalculate msed from mass balance
# Code to model PCB 19 in laboratory experiments
# using sediment from Altavista, VI. Passive measurements
# of PCB 19 in the water and the air phases are predicted and
# linked to the water and air concentrations from the passive
# samplers.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")
install.packages("FME")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
  library(gridExtra)
  library(FME)
}

# Read data ---------------------------------------------------------------
{
  obs.data <- read.csv("Data/AVL_S_data_long.csv", check.names = FALSE)
  # Physical chemical properties
  pcp.data <- read.csv("Data/AVL_S_PCP.csv")
}

# Extract individual PCBi -------------------------------------------------
{
  pcb.name <- "PCB19"
  obs.data.pcbi <- obs.data[, c("sampler", "time", pcb.name)]
  pcp.data.pcbi <- pcp.data[pcp.data$congener == pcb.name, ]
}

# Organize data -----------------------------------------------------------
{
  spme <- obs.data.pcbi %>%
    filter(sampler == "SPME")
  puf <- obs.data.pcbi %>%
    filter(sampler == "PUF") %>%
    mutate(across(starts_with("PCB"),  ~ . * 10))
  sed <- obs.data.pcbi %>%
    filter(sampler == "sed")
}

# Reactive transport function ---------------------------------------------
rtm.PCB19 = function(t, state, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- pcp.data.pcbi$MW # g/mol PCB 19 molecular weight
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- pcp.data.pcbi$Kaw # PCB 19 dimensionless Henry's law constant @ 25 C
  dUaw <- pcp.data.pcbi$dUaw # internal energy for the transfer of air-water for PCB 19 (J/mol)
  Kaw.t <- Kaw*exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  Kow <- pcp.data.pcbi$Kow # PCB 19 octanol-water equilibrium partition coefficient (low value!!)
  dUow <-  pcp.data.pcbi$dUow # internal energy for the transfer of octanol-water for PCB 19 (J/mol)
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Koa <- pcp.data.pcbi$Koa # PCB 19 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Apuf <- 7.07 # cm2
  Vpuf <- 29 # cm3 volume of PUF
  d <- 0.0213*100^3 # g/m3 density of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) # PCB 19-PUF equilibrium partition coefficient [m3/g]
  Kpuf <- Kpuf * d # [La/Lpuf]
  
  # SPME fiber constants
  Aspme <- 0.138 # cm2/cm SPME area
  Vspme <- 0.000000069 * 1000 # cm3/cm SPME volume/area
  Lspme <- 1 # cm SPME length normalization to 1 cm
  Kspme <- 10^(1.06 * log10(Kow.t) - 1.16) # PCB 19-SPME equilibrium partition coefficient
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 19's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 19's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 19
  
  # kaw calculations (air-water mass transfer coefficient)
  # i) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 19. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 19
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 19, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioremediation rate
  kb <- parms$kb
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  
  # Passive sampler rates
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  ro <- parms$ro # cm/d sampling rate for PUF
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  Cspme <- state[3]
  Ca <- state[4]
  Cpuf <- state[5]
  
  dCsdt <- - f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw
  dCwdt <- - ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) - 
    ko * Aspme * Lspme / Vw * (Cw - Cspme / Kspme) -
    kb * Cw # [ng/L]
  dCspmedt <- ko * Aspme / Vspme * (Cw - Cspme / Kspme) # Cw = [ng/L], Cspme = [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf) # Ca = [ng/L]
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf) # Ca = [ng/L], Cpuf = [ng/L]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dCspmedt, dCadt, dCpufdt)))
}

Ct <- obs.data.pcbi %>% filter(sampler == "sed") %>% pull(PCB19) # ng/g
M <- 0.1   # kg/L
Cs0 <- Ct * M * 1000  # [ng/L]
cinit <- c(Cs = Cs0, Cw = 0, Cspme = 0, Ca = 0, Cpuf = 0)

# 2. Prepare observed data (do this FIRST)
obs.data.pcbi.2 <- obs.data.pcbi %>%
  filter(sampler %in% c("SPME", "PUF")) %>%
  # First check for duplicates
  group_by(time, sampler) %>%
  summarise(PCB19 = mean(PCB19, na.rm = TRUE),  # Average duplicates
            .groups = "drop") %>%
  # Then pivot
  tidyr::pivot_wider(
    names_from = sampler, 
    values_from = PCB19,
    values_fn = mean  # Explicitly handle duplicates
  ) %>%
  mutate(
    PUF = as.numeric(PUF) * 10,  # Convert to numeric before multiplication
    SPME = as.numeric(SPME)
  )

# 3. Define time points from data
t.1 <- unique(obs.data.pcbi.2$time)

# 4. Define cost function WITH data in scope
cost_func <- function(parms) {
  parms_list <- list(
    ro = parms["ro"],
    ko = parms["ko"],
    kdf = parms["kdf"],
    kds = parms["kds"],
    f = parms["f"],
    ka = parms["ka"],
    kb = 0
  )
  
  out <- ode(y = cinit, times = t.1, func = rtm.PCB19, parms = parms_list)
  out_df <- as.data.frame(out)
  
  # Calculate predicted masses
  pred_spme <- out_df$Cspme * 6.9e-8  # [ng/cm]
  pred_puf <- out_df$Cpuf * 29 / 1000  # [ng/puf]
  
  # Get matching observed values
  obs <- obs.data.pcbi.2 %>% 
    filter(time %in% out_df$time)
  
  # Return residuals
  c(pred_spme - obs$SPME, pred_puf - obs$PUF)
}

# 5. Parameter fitting
par_guess <- c(ro = 540, ko = 10, kdf = 1.9, kds = 0.001, f = 0.8, ka = 150)
lower_bounds <- c(ro = 1, ko = 1e-3, kdf = 1e-3, kds = 1e-6, f = 0, ka = 1)
upper_bounds <- c(ro = 1e4, ko = 1e3, kdf = 10, kds = 1, f = 1, ka = 1e3)

fit <- modFit(f = cost_func, p = par_guess, 
              lower = lower_bounds, upper = upper_bounds)

# 6. Extract best parameters
best_parms <- as.list(fit$par)
best_parms$kb <- 0

# 7. Run final model
final_out <- ode(y = cinit, times = t.1, func = rtm.PCB19, parms = best_parms)
final_out_df <- as.data.frame(final_out) %>%
  mutate(
    mspme = Cspme * 6.9e-8,
    mpuf = Cpuf * 29 / 1000,
    Mt = (Cs * 10 / 100) + (Cw * 0.1) + mspme + (Ca * 0.125) + mpuf
  )

# 8. Visualization
ggplot(final_out_df) +
  geom_line(aes(time, mspme, color = "Predicted SPME")) +
  geom_point(data = obs.data.pcbi.2, aes(time, SPME, color = "Observed SPME")) +
  geom_line(aes(time, mpuf, color = "Predicted PUF")) +
  geom_point(data = obs.data.pcbi.2, aes(time, PUF, color = "Observed PUF")) +
  labs(title = "Model Fit", x = "Time", y = "Mass") +
  scale_color_manual(values = c(
    "Observed SPME" = "red", "Predicted SPME" = "pink",
    "Observed PUF" = "blue", "Predicted PUF" = "lightblue"
  ))




