# Code to model PCB 32 in laboratory experiments
# using sediment from NBH (site INT 222). Passive measurements
# of PCB 32 in the water and the air phases are predicted and
# linked to the water and air concentrations from the passive
# samplers. Control experiment, no biochar, no LB400

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
  exp.data <- read.csv("Data/uncoated_biochar_V2.csv")
  # Select individual congener from datasets
  pcb.ind <- "PCB_32"
  # Extract relevant columns
  pcbi <- exp.data[, c("Sample_medium", "Experiment", "percent_biochar",
                       "Group", "time", "Replicate", pcb.ind)]
}

# Organize data -----------------------------------------------------------
{
  # Using time series experiments
  # Pull congener-specific data from the dataset without averaging
  # Select SPME control samples
  pcbi.spme.control.5 <- pcbi %>%
    filter(Sample_medium == "SPME", Experiment == "biochar_timeseries",
           Group == "Control", percent_biochar == 5.0) %>%
    rename("mf_control_5" = PCB_32)
  
  # Select PUF control samples
  pcbi.puf.control.5 <- pcbi %>%
    filter(Sample_medium == "PUF", Experiment == "biochar_timeseries",
           Group == "Control", percent_biochar == 5.0) %>%
    rename("mpuf_control_5" = PCB_32)
  
  # Combine the mf and mpuf data for Control
  pcb_combined_control_5 <- cbind(
    pcbi.spme.control.5 %>%
      select(time, mf_control_5),
    pcbi.puf.control.5 %>%
      select(mpuf_control_5)
  )
  # Add a row for time = 0
  pcb_combined_control_5 <- rbind(
    data.frame(time = 0, mf_control_5 = 0, mpuf_control_5 = 0),
    pcb_combined_control_5
  )
}

# Reactive transport function ---------------------------------------------
# ---- ODE function (minimal, no defensive checks) ----
rtm.PCB32 <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # --- geometry / fixed geometry values (cm3, cm2)
    Vw  <- 100    # cm3 water volume
    Vpw <- 4      # cm3 porewater volume
    Va  <- 125    # cm3 headspace
    Aaw <- 20     # cm2 air-water area
    Aws <- 30     # cm2 sediment-water area
    
    # --- sediment mass (g) and compute Vs (cm3) from porosity + density
    ms <- 10 # g
    n  <- 0.42
    ds <- 1540        # g/L
    M  <- ds * (1 - n) / n           # g solids per L porewater
    Vs <- ms / M * 1000        # cm3 porewater associated with ms_local
    
    # --- congener-specific constants (as you provided)
    MH2O <- 18.0152; MCO2 <- 44.0094; MW.pcb <- 257.532
    R <- 8.3144
    Tst.1 <- 273.15 + 25; Tw.1 <- 273.15 + 20
    
    Kaw <- 0.016011984; dUaw <- 52.59022
    Kaw.t <- Kaw * exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
    Kow <- 10^(5.44); dUow <- -22.88894
    Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 - 1 / Tst.1))
    Koa <- 10^(7.709482239)
    
    # PUF & SPME
    Apuf <- 7.07
    Vpuf <- 29
    d <- 0.0213 * 100^3
    Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774)
    Kpuf <- Kpuf * d
    Af <- 0.138
    Vf <- 0.000000069 * 1000 # [cm3/cm SPME]
    L  <- 1
    Vf_tot <- Vf * L
    Kf <- 10^(1.06 * log10(Kow.t) - 1.16)
    
    # Diffusion / transfer coefficients
    D.water.air <- 0.2743615
    D.co2.w <- 1.67606E-05
    D.pcb.air <- D.water.air * (MW.pcb/MH2O)^(-0.5)
    D.pcb.water <- D.co2.w * (MW.pcb/MCO2)^(-0.5)
    v.H2O <- 0.010072884
    V.water.air <- 0.003
    V.co2.w <- 4.1*10^-2
    SC.pcb.w <- v.H2O / D.pcb.water
    
    bl <- 0.21
    kpw <- D.pcb.water * 60 * 60 * 24 / bl   # cm/day
    
    Kaw.a <- V.water.air * (D.pcb.air/D.water.air)^(0.67)
    Kaw.w <- V.co2.w * (SC.pcb.w/600)^(-0.5)
    kaw.o <- (1 / (Kaw.a * Kaw.t) + (1 / Kaw.w))^-1
    kaw.o <- kaw.o * 100 * 60 * 60 * 24   # cm/day
    
    ksed <- 6.255 # [1/d] from optimization code
    
    # sampler rates / kb from parms
    ko <- parms$ko
    ro <- parms$ro
    kb <- parms$kb
    Kd <- parms$Kd
    
    # --- state variables (order must match cinit)
    Cs   <- state[1]   # ng/g (solid)
    # [ng/L] -> [ng/cm3]
    Cpw  <- state[2] / 1000
    Cw   <- state[3] / 1000
    Cf   <- state[4] / 1000
    Ca   <- state[5] / 1000
    Cpuf <- state[6] / 1000
    
    # --- current Cs -> porewater-equivalent [ng/L] -> [ng/cm3]
    Cs_pw_eq <- Cs / Kd
    
    # --- ODEs (mass-consistent)
    dCsdt  <- - ksed * Vs / ms * (Cs_pw_eq - Cpw)
    
    dCpwdt <-   ksed * Vs / Vpw * (Cs_pw_eq - Cpw) -
      kpw * Aws / Vpw * (Cpw - Cw) -
      kb * Cpw
    
    dCwdt <- kpw * Aws / Vw * (Cpw - Cw) -
      kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
      ko * Af * L / Vw * (Cw - Cf / Kf)
    
    dCfdt <- ko * Af / Vf_tot * (Cw - Cf / Kf)
    
    dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
      ro * Apuf / Va * (Ca - Cpuf / Kpuf)
    
    dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf)
    
    # Convert back to ng/L/day
    return(list(c(dCsdt,
                  dCpwdt * 1000,
                  dCwdt * 1000,
                  dCfdt * 1000,
                  dCadt * 1000,
                  dCpufdt * 1000)))
  })
}

# ---- Outside: compute Kd, initial conditions, parms, run ----
Ct <- 520   # ng/g sediment measured
E <- 1.74; S <- 1.35; A <- 0; B <- 0.17; V <- 1.6914

logKoc <- 1.1 * E - 0.72 * S + 0.15 * A - 1.98 * B + 2.28 * V + 0.14
Koc <- 10^(logKoc)
foc <- 0.03
# Add PCB sorption to biochar
Kbc <- 10^(3.8) # [Lw/KgBC] From Dong et al 2025
fbc <- 0.05 # 5% of total sediment
Kd  <- (Koc * foc  + Kbc * fbc) / (1 + fbc) # L/kg sediment

Cpw0 <- Ct * 1000 / Kd   # ng/L

cinit <- c(Cs = Ct, Cpw = Cpw0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)

parms <- list(ro = 420, ko = 3, kb = 0, Kd = Kd)

t.1 <- unique(pcb_combined_control_5$time)
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB32, parms = parms)

# ---- post-process masses (ng) ----
df.1 <- as.data.frame(out.1)
colnames(df.1) <- c("time","Cs","Cpw","Cw","Cf","Ca","Cpuf")

msed_g  <- 10
Vpw_cm3 <- 4
Vw_cm3  <- 100
Va_cm3  <- 125
Vpuf_cm3 <- 29

Vpw_L  <- Vpw_cm3 / 1000     # L
Vw_L   <- Vw_cm3  / 1000
Va_L   <- Va_cm3  / 1000
Vpuf_L <- Vpuf_cm3 / 1000

# SPME: you supplied Vf as cm3 per cm of fiber
Vf_cm3_per_cm <- 0.000000069 * 1000  # cm3 per cm (keep your original source)
fiber_length_cm <- 1                  # exposed fiber length in cm (you used L=1)
Vf_cm3_total <- Vf_cm3_per_cm * fiber_length_cm
Vf_L <- Vf_cm3_total / 1000           # L

# ---- compute compartment masses (ng) ----
# Cs is ng/g, multiply by sediment mass (g) -> ng
df.1$ms   <- df.1$Cs * msed_g

# aqueous concentrations are ng/L: multiply by corresponding volume (L) -> ng
df.1$mpw  <- df.1$Cpw  * Vpw_L     # porewater mass (ng)
df.1$mw   <- df.1$Cw   * Vw_L      # water column mass (ng)
df.1$mf   <- df.1$Cf   * Vf_L      # SPME total mass (ng)
df.1$ma   <- df.1$Ca   * Va_L      # air mass (ng)
df.1$mpuf <- df.1$Cpuf * Vpuf_L    # PUF mass (ng)

# ---- total mass & fractions (guard against zero total) ----
df.1$mt <- df.1$ms + df.1$mpw + df.1$mw + df.1$mf + df.1$ma + df.1$mpuf

# Ensure observed data is in a tibble
observed_data <- as_tibble(pcb_combined_control_5) %>%
  select(time, mf_control_5, mpuf_control_5)

# Convert model results to tibble and select relevant columns
model_results <- as_tibble(df.1) %>%
  mutate(time = as.numeric(time)) %>%
  select(time, mf, mpuf)

# Merge model results with observed data
comparison_data <- model_results %>%
  left_join(observed_data, by = "time")

# Calculate the averages of mf and mpuf within each group (e.g., per time)
grouped_comparison <- comparison_data %>%
  group_by(time) %>%  # Adjust the grouping variable if needed
  summarise(
    avg_mf_model = mean(mf, na.rm = TRUE),
    avg_mf_observed = mean(mf_control_5, na.rm = TRUE),
    avg_mpuf_model = mean(mpuf, na.rm = TRUE),
    avg_mpuf_observed = mean(mpuf_control_5, na.rm = TRUE)
  )

# Define function to calculate R-squared, handling NA values
mf_r2 <- function(predicted, observed) {
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
mf_r2_value <- mf_r2(grouped_comparison$avg_mf_model, grouped_comparison$avg_mf_observed)
mpuf_r2_value <- mf_r2(grouped_comparison$avg_mpuf_model, grouped_comparison$avg_mpuf_observed)

# Print R-squared values
print(paste("R-squared for mf (average): ", mf_r2_value))
print(paste("R-squared for mpuf (average): ", mpuf_r2_value))

# Plot
{  
  # Run the model with the new time sequence
  cinit <- c(Cs = Ct, Cpw = Cpw0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
  t_daily <- seq(0, 130, by = 1)  # Adjust according to your needs
  out_daily <- ode(y = cinit, times = t_daily, func = rtm.PCB32,
                   parms = parms)
  head(out_daily)
  
  # Transform Cf and Cpuf to mass/cm and mass/puf
  out.daily <- as.data.frame(out_daily)
  colnames(out.daily) <- c("time", "Cs", "Cpw", "Cw", "Cf", "Ca", "Cpuf")
  
  # Calculate Mf and Mpuf based on volumes
  out.daily$mf <- out.daily$Cf * Vf_L # [ng]
  out.daily$mpuf <- out.daily$Cpuf * Vpuf_L  # [ng]
  
  # Convert model results to tibble and ensure numeric values
  model_results_daily_clean <- as_tibble(out.daily) %>%
    mutate(across(c(mf, mpuf, Cw, time), as.numeric)) %>%  # Ensure all relevant columns are numeric
    select(time, mf, mpuf)
  
  # Export data
  #write.csv(model_results_daily_clean, file = "Output/Data/RTM/PCB32Control.csv")
  
  # Prepare model data for plotting
  model_data_long <- model_results_daily_clean %>%
    pivot_longer(cols = c(mf, mpuf), 
                 names_to = "variable", 
                 values_to = "model_value") %>%
    mutate(type = "Model")
  
  # Clean observed data and prepare for plotting
  observed_data_clean <- observed_data %>%
    pivot_longer(cols = c(mf_control_5, mpuf_control_5), 
                 names_to = "variable", 
                 values_to = "observed_value") %>%
    mutate(variable = recode(variable, 
                             "mf_control_5" = "mf", 
                             "mpuf_control_5" = "mpuf"),
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
  p_mf <- ggplot(plot_data_daily %>% filter(variable == "mf"), aes(x = time)) +
    geom_line(data = . %>% filter(type == "Model"),
              aes(y = value, color = "Model"), linewidth = 1) +
    geom_point(data = . %>% filter(type == "Observed"),
               aes(y = value, color = "Observed"), size = 2) +
    labs(x = "Time", y = "mf [ng/cm]") +
    scale_color_manual(values = c("Model" = "blue", "Observed" = "red")) +
    theme_bw() +
    theme(legend.title = element_blank())
  
  # Plot mpuf
  p_mpuf <- ggplot(plot_data_daily %>% filter(variable == "mpuf"), aes(x = time)) +
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
p.32 <- grid.arrange(p_mf, p_mpuf, ncol = 2)

# Save plot in folder
ggsave("Output/Plots/RTM/PCB32Control.png", plot = p.32, width = 6,
       height = 5, dpi = 500)

