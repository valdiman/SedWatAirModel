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
rtm.PCBi = function(t, state, parms){
  
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

Ct <- obs.data.pcbi %>% filter(sampler == "sed") %>% pull(!!sym(pcb.name)) # ng/g
M <- 0.1   # kg/L
Cs0 <- Ct * M * 1000  # [ng/L]
cinit <- c(Cs = Cs0, Cw = 0, Cspme = 0, Ca = 0, Cpuf = 0)

# 2. Prepare observed data (do this FIRST)
obs.data.pcbi.2 <- obs.data.pcbi %>%
  filter(sampler %in% c("SPME", "PUF")) %>%
  group_by(time, sampler) %>%
  summarise(!!sym(pcb.name) := mean(!!sym(pcb.name), na.rm = TRUE),  # Dynamic mean
            .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sampler, 
    values_from = all_of(pcb.name),  # Dynamic column
    values_fn = mean
  ) %>%
  mutate(
    PUF = as.numeric(PUF) * 10,
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
  
  out <- ode(y = cinit, times = t.1, func = rtm.PCBi, parms = parms_list)
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

# See parameters
summary(fit)

# 6. Extract best parameters
best_parms <- as.list(fit$par)
best_parms$kb <- 0

# 7. Run final model
final_out <- ode(y = cinit, times = t.1, func = rtm.PCBi, parms = best_parms)
final_out_df <- as.data.frame(final_out) %>%
  mutate(
    mspme = Cspme * 6.9e-8,
    mpuf = Cpuf * 29 / 1000,
    Mt = (Cs * 10 / 100) + (Cw * 0.1) + mspme + (Ca * 0.125) + mpuf
  )

# 8. Visualization
# Create daily time points (adjust max time as needed)
fine_time <- seq(from = min(t.1), to = max(t.1), by = 0.5)  # 0.5-day intervals

# Run model with fine time resolution
fine_out <- ode(y = cinit, times = fine_time, func = rtm.PCBi, parms = best_parms)
fine_out_df <- as.data.frame(fine_out) %>%
  mutate(
    mspme = Cspme * 6.9e-8,  # [ng/cm]
    mpuf = Cpuf * 29 / 1000,   # [ng/puf]
    Mt = (Cs * 0.01) + (Cw * 0.1) +
      mspme + (Ca * 0.125) + mpuf
  )

# SPME Plot
spme_plot <- ggplot() +
  geom_line(data = fine_out_df, 
            aes(time, mspme, color = "Predicted"), 
            linewidth = 0.8) +
  geom_point(data = obs.data.pcbi.2, 
             aes(time, SPME, color = "Observed"),
             size = 3, shape = 19) +
  labs(title = "SPME: Model Fit", x = "Time (days)", y = "Mass (ng/cm)") +
  scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(spme_plot)

# PUF Plot
puf_plot <- ggplot() +
  geom_line(data = fine_out_df, 
            aes(time, mpuf, color = "Predicted"), 
            linewidth = 0.8) +
  geom_point(data = obs.data.pcbi.2, 
             aes(time, PUF, color = "Observed"),
             size = 3, shape = 17) +  # Triangle points
  labs(title = "PUF: Model Fit", x = "Time (days)", y = "Mass (ng/puf)") +
  scale_color_manual(values = c("Observed" = "green4", "Predicted" = "orange2")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(puf_plot)

# Calculate goodness-of-fit metrics ----------------------------------------

# 1. Prepare observed and predicted data
obs_pred <- obs.data.pcbi.2 %>%
  left_join(as.data.frame(final_out), by = "time") %>%
  mutate(
    pred_spme = Cspme * 6.9e-8,
    pred_puf = Cpuf * 29 / 1000
  )

# 2. Calculate metrics for SPME
spme_metrics <- obs_pred %>%
  filter(!is.na(SPME)) %>%
  summarise(
    RMSE_spme = sqrt(mean((SPME - pred_spme)^2)),
    MAE_spme = mean(abs(SPME - pred_spme)),
    R2_spme = cor(SPME, pred_spme)^2,
    NSE_spme = 1 - sum((SPME - pred_spme)^2)/sum((SPME - mean(SPME))^2)
  )

# 3. Calculate metrics for PUF
puf_metrics <- obs_pred %>%
  filter(!is.na(PUF)) %>%
  summarise(
    RMSE_puf = sqrt(mean((PUF - pred_puf)^2)),
    MAE_puf = mean(abs(PUF - pred_puf)),
    R2_puf = cor(PUF, pred_puf)^2,
    NSE_puf = 1 - sum((PUF - pred_puf)^2)/sum((PUF - mean(PUF))^2)
  )

# 4. Print metrics
cat("\nSPME Metrics:\n")
print(spme_metrics)
cat("\nPUF Metrics:\n")
print(puf_metrics)

# 5. Add metrics to plots
spme_plot <- spme_plot +
  annotate("text", x = max(t.1)*0.7, y = max(obs_pred$SPME, na.rm=TRUE)*0.9,
           label = paste0("R² = ", round(spme_metrics$R2_spme, 3), "\n",
                          "RMSE = ", round(spme_metrics$RMSE_spme, 3)),
           hjust = 0, size = 4)

puf_plot <- puf_plot +
  annotate("text", x = max(t.1)*0.7, y = max(obs_pred$PUF, na.rm=TRUE)*0.9,
           label = paste0("R² = ", round(puf_metrics$R2_puf, 3), "\n",
                          "RMSE = ", round(puf_metrics$RMSE_puf, 3)),
           hjust = 0, size = 4)

# Re-print plots
print(spme_plot)
print(puf_plot)


