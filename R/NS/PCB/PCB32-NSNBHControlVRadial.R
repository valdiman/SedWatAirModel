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
exp.data <- read.csv("Data/01_NS_SPME_PUF.csv")
pcb.ind <- "PCB_32"
pcbi <- exp.data[, c("ID", "Group", "time", "Sample_medium", pcb.ind)]

# Organize data -----------------------------------------------------------
pcbi.spme.control <- pcbi %>%
  filter(ID == "NBH_NS", Group == "Control",
         Sample_medium == "SPME") %>%
  rename(mf_control = !!sym(pcb.ind))

pcbi.puf.control <- pcbi %>%
  filter(ID == "NBH_NS", Group == "Control",
         Sample_medium == "PUF") %>%
  rename(mpuf_control = !!sym(pcb.ind))

pcb_combined_control <- cbind(
  pcbi.spme.control %>% select(time, mf_control),
  pcbi.puf.control %>% select(mpuf_control)
)

# Add a row for time = 0
pcb_combined_control <- rbind(
  data.frame(time = 0, mf_control = 0, mpuf_control = 0),
  pcb_combined_control
)

# Fixed phys-chem and geometry / precomputed ------------------------------
pc <- read.csv("Data/04_PCP.csv", stringsAsFactors = FALSE)
pc_row <- pc[pc$congener == pcb.ind, ]

MW.pcb <- pc_row$MW
Kow    <- pc_row$Kow
dUow   <- pc_row$dUow
Kaw    <- pc_row$Kaw
dUaw   <- pc_row$dUaw
Koa    <- pc_row$Koa
E <- pc_row$E; S <- pc_row$S; A <- pc_row$A; B <- pc_row$B; V <- pc_row$V

# geometry / fixed values
Vw_cm3   <- 100    # cm3 water volume
Vpw_cm3  <- 4      # cm3 porewater volume
Va_cm3   <- 125    # cm3 headspace
Vpuf_cm3 <- 29     # cm3 PUF
Aaw <- 20     # cm2 air-water area
Aws <- 30     # cm2 sediment-water area
ms_g <- 10    # g sediment in the experimental sediment layer

# derived volumes in liters
Vw_L   <- Vw_cm3   / 1000
Vpw_L  <- Vpw_cm3  / 1000
Va_L   <- Va_cm3   / 1000
Vpuf_L <- Vpuf_cm3 / 1000

# compute Vs (cm3 porewater associated with ms_g)
n  <- 0.42
ds <- 1540      # g / L (sediment dry density)
M  <- ds * (1 - n) / n     # g solids per L porewater
Vs <- ms_g / M * 1000      # cm3 porewater associated with ms_g

# compute Kd once (Koc model) and document units
logKoc <- 1.1 * E - 0.72 * S + 0.15 * A - 1.98 * B + 2.28 * V + 0.14
Koc <- 10^(logKoc)
foc <- 0.03
Kd  <- Koc * foc   # L/kg sediment

# Temperature-corrected partitioning coefficients
MH2O <- 18.0152; MCO2 <- 44.0094; R <- 8.3144
Tst.1 <- 273.15 + 25; Tw.1 <- 273.15 + 20
Kaw.t <- Kaw * exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 - 1 / Tst.1))

# diffusion / transfer (documented magic numbers)
D.water.air <- 0.2743615
D.co2.w <- 1.67606E-05
D.pcb.air <- D.water.air * (MW.pcb/MH2O)^(-0.5)
D.pcb.water <- D.co2.w * (MW.pcb/MCO2)^(-0.5)
v.H2O <- 0.010072884
V.water.air <- 0.003
V.co2.w <- 4.1*10^-2
SC.pcb.w <- v.H2O / D.pcb.water
bl <- 0.1
kpw <- D.pcb.water * 60 * 60 * 24 / bl   # cm/day

# Kaw combined (air+water resistances) then converted to cm/day
Kaw.a <- V.water.air * (D.pcb.air/D.water.air)^(0.67)
Kaw.w <- V.co2.w * (SC.pcb.w/600)^(-0.5)
kaw.o <- (1 / (Kaw.a * Kaw.t) + (1 / Kaw.w))^-1
kaw.o <- kaw.o * 100 * 60 * 60 * 24   # cm/day

# PUF & SPME derived
Apuf <- 7.07
Vpuf <- Vpuf_cm3
dpuf <- 0.0213 * 100^3
Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) * dpuf

Af <- 0.138
Vf_cm3_per_cm <- 0.000000069 * 1000
fiber_length_cm <- 1
Vf_cm3_total <- Vf_cm3_per_cm * fiber_length_cm
Vf_L <- Vf_cm3_total / 1000
L  <- 1
Vf_tot <- Vf_cm3_total * L
Kf <- 10^(1.06 * log10(Kow.t) - 1.16)

# -------------------------
# 3) Desorption_single: Werner single-domain truncated series
desorption_single <- function(t_days, D_r2, mtot_ng, Nterms = 200) {
  coef <- 6 / (pi^2)
  n <- seq_len(Nterms)
  tvec <- as.numeric(t_days)
  mfrac <- numeric(length(tvec))
  dm_dt <- numeric(length(tvec))
  for (i in seq_along(tvec)) {
    t <- tvec[i]
    exp_term <- exp(- (D_r2) * (n^2) * (pi^2) * t)
    sum_term <- sum((1 / (n^2)) * exp_term)
    mfrac[i] <- 1 - coef * sum_term
    dm_dt[i] <- mtot_ng * ( coef * sum( D_r2 * (pi^2) * exp_term ) )
  }
  return(list(mfrac = mfrac, dm_dt = dm_dt))
}

# ---------- assemble parms ----------
parms <- list(
  # fitted/used rates (set directly)
  ro = 33.2319,        # fitted ro (cm/day)
  ko = 3/50, kb = 0,
  # fixed params & precomputed
  Kd = Kd, MW.pcb = MW.pcb,
  Vw = Vw_cm3, Vpw = Vpw_cm3, Va = Va_cm3, Aws = Aws, Aaw = Aaw,
  ms_g = ms_g, Vs = Vs,
  Kaw.t = Kaw.t, Kow.t = Kow.t,
  kpw = kpw, kaw.o = kaw.o,
  Apuf = Apuf, Vpuf = Vpuf, Kpuf = Kpuf,
  Af = Af, Vf_tot = Vf_tot, Kf = Kf,
  L = L,
  # desorption-specific defaults (use fitted D_r2)
  D_r2 = 8.04777e-06,   # fitted D/r^2 (day^-1)
  M_sed_init = NA,      # will set below from Cs_init * ms_g
  Nterms_des = 200
)

# ---- initial total mass and initial conditions ----
bulk_conc <- read.csv("Data/03_NBH_SedimentPCB.csv", stringsAsFactors = FALSE)
Ct <- mean(bulk_conc[[pcb.ind]])   # ng/g sediment mean measured

# Total mass (ng) in sediment layer (use ms_g in g)
M_sed_init <- Ct * ms_g     # ng (ms_g in g; Ct in ng/g)
parms$M_sed_init <- M_sed_init

# Solve for Cs_init so that M_sed_init = Cs*ms_g + dissolved mass in porewater
Cs_init <- M_sed_init / (ms_g + (1000 * Vpw_L) / Kd)  # ng/g
Cpw_init <- Cs_init * 1000 / Kd                      # ng/L

# include Cf initial (SPME) as 0 explicitly; order must match ODE: Cs, Cpw, Cw, Cf, Ca, Cpuf
cinit <- c(Cs = Cs_init, Cpw = Cpw_init, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)

# ---- ODE uses parms only, minimal inline computation ----
rtm.PCB <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # State (input units)
    # Cs: ng/g (state[1])
    # Cpw, Cw, Cf, Ca, Cpuf: ng/L (state[2:6])
    # Convert aqueous/gas to ng/cm3 for flux calc (1 L = 1000 cm3 -> divide by 1000)
    Cs   <- state[1]                 # ng/g
    Cpw_cm3  <- state[2] / 1000      # ng/cm3
    Cw_cm3   <- state[3] / 1000
    Cf_cm3   <- state[4] / 1000
    Ca_cm3   <- state[5] / 1000
    Cpuf_cm3 <- state[6] / 1000
    
    # ---- use desorption_single to supply dm_dt from sediment ----
    # desorption_single returns dm_dt in ng/day (mass released from solids)
    des <- desorption_single(t_days = t, D_r2 = parms$D_r2, mtot_ng = parms$M_sed_init, Nterms = parms$Nterms_des)
    dm_dt_ng_per_day <- des$dm_dt[1]
    
    # dCs/dt now directly from dm_dt (ng/g/day)
    dCsdt  <- - dm_dt_ng_per_day / ms_g
    
    # dCpw/dt receives the dm_dt divided into the porewater volume (ng/cm3/day)
    dCpwdt <-   (dm_dt_ng_per_day / Vpw) -
      kpw * Aws / Vpw * (Cpw_cm3 - Cw_cm3) -
      kb * Cpw_cm3
    
    # rest unchanged
    dCwdt <- kpw * Aws / Vw * (Cpw_cm3 - Cw_cm3) -
      kaw.o * Aaw / Vw * (Cw_cm3 - Ca_cm3 / Kaw.t) -
      ko * Af * L / Vw * (Cw_cm3 - Cf_cm3 / Kf)
    
    dCfdt <- ko * Af / Vf_tot * (Cw_cm3 - Cf_cm3 / Kf)
    
    dCadt <- kaw.o * Aaw / Va * (Cw_cm3 - Ca_cm3 / Kaw.t) -
      ro * Apuf / Va * (Ca_cm3 - Cpuf_cm3 / Kpuf)
    
    dCpufdt <- ro * Apuf / Vpuf * (Ca_cm3 - Cpuf_cm3 / Kpuf)
    
    # return derivatives in same units as state
    return(list(c(
      dCsdt,
      dCpwdt * 1000,
      dCwdt * 1000,
      dCfdt * 1000,
      dCadt * 1000,
      dCpufdt * 1000
    )))
  })
}

# ---- run ----
t.1 <- sort(unique(pcb_combined_control$time))
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB, parms = parms, rtol = 1e-6, atol = 1e-8)

# ---- post-process masses (ng) ----
df.1 <- as.data.frame(out.1)
colnames(df.1) <- c("time","Cs","Cpw","Cw","Cf","Ca","Cpuf")

# convert L volumes (we already have Vw_L etc.)
# compute compartment masses (ng)
df.1$ms   <- df.1$Cs * ms_g
df.1$mpw  <- df.1$Cpw  * Vpw_L     # porewater mass (ng)
df.1$mw   <- df.1$Cw   * Vw_L      # water column mass (ng)
df.1$mf   <- df.1$Cf   * Vf_L      # SPME total mass (ng)
df.1$ma   <- df.1$Ca   * Va_L      # air mass (ng)
df.1$mpuf <- df.1$Cpuf * Vpuf_L    # PUF mass (ng)

# ---- total mass & fractions (guard against zero total) ----
df.1$mt <- df.1$ms + df.1$mpw + df.1$mw + df.1$mf + df.1$ma + df.1$mpuf

print(df.1)

# Comparison and plot
{
  # Ensure observed data is in a tibble
  observed_data <- as_tibble(pcb_combined_control) %>%
    select(time, mf_control, mpuf_control)
  
  # Convert model results to tibble and select relevant columns
  model_results <- as_tibble(df.1) %>%
    mutate(time = as.numeric(time),
           mf   = Cf   * Vf_L,
           mpuf = Cpuf * Vpuf_L) %>%
    select(time, mf, mpuf)
  
  # Merge model results with observed data
  comparison_data <- model_results %>%
    left_join(observed_data, by = "time")
  
  # Calculate the averages of mf and mpuf within each group (e.g., per time)
  grouped_comparison <- comparison_data %>%
    group_by(time) %>%  # Adjust the grouping variable if needed
    summarise(
      avg_mf_model = mean(mf, na.rm = TRUE),
      avg_mf_observed = mean(mf_control, na.rm = TRUE),
      avg_mpuf_model = mean(mpuf, na.rm = TRUE),
      avg_mpuf_observed = mean(mpuf_control, na.rm = TRUE)
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
  # Run the model with the new time sequence
  cinit <- c(Cs = Cs_init, Cpw = Cpw_init, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
  t_daily <- seq(0, 80, by = 1)  # Adjust according to your needs
  out_daily <- ode(y = cinit, times = t_daily, func = rtm.PCB,
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
  #write.csv(model_results_daily_clean, file = "Output/Data/NS/PCB32ControlNBH.csv")
  
  # Prepare model data for plotting
  model_data_long <- model_results_daily_clean %>%
    pivot_longer(cols = c(mf, mpuf), 
                 names_to = "variable", 
                 values_to = "model_value") %>%
    mutate(type = "Model")
  
  # Clean observed data and prepare for plotting
  observed_data_clean <- observed_data %>%
    pivot_longer(cols = c(mf_control, mpuf_control), 
                 names_to = "variable", 
                 values_to = "observed_value") %>%
    mutate(variable = recode(variable, 
                             "mf_control" = "mf", 
                             "mpuf_control" = "mpuf"),
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
plot <- grid.arrange(p_mf, p_mpuf, ncol = 2)

# Save plot in folder
ggsave("Output/Plot/NS/PCB32ControlNBH.png", plot = plot, width = 6,
       height = 5, dpi = 500)

