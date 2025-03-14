## Script to investigate mass transfer coefficients
# for PCBs using sediment from Altavista, VI. for shaken experiments.

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
  obs.data <- read.csv("Data/AVL_S_data_long.csv")
  # Physical chemical properties
  pcp.data <- read.csv("Data/AVL_S_PCP.csv")
}

# Reactive transport function ---------------------------------------------
rtm.PCB = function(t, state, parms, pcb_index){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- pcp.data$MW[pcb_index] # g/mol PCBi molecular weight
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
  Kaw <- pcp.data$Kaw[pcb_index] # PCBi dimensionless Henry's law constant @ 25 C
  dUaw <- pcp.data$dUaw[pcb_index] # internal energy for the transfer of air-water for PCBi (J/mol)
  Kaw.t <- Kaw*exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  Kow <- pcp.data$Kow[pcb_index] # PCBi octanol-water equilibrium partition coefficient
  dUow <-  pcp.data$dUow[pcb_index] # internal energy for the transfer of octanol-water for PCBi (J/mol)
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Koa <- pcp.data$Koa[pcb_index] # PCBi octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Apuf <- 7.07 # cm2
  Vpuf <- 29 # cm3 volume of PUF
  d <- 0.0213*100^3 # g/m3 density of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) # PCB 19-PUF equilibrium partition coefficient [m3/g]
  Kpuf <- Kpuf * d # [La/Lpuf]
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 * 1000 # cm3/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow.t) - 1.16) # PCB 19-SPME equilibrium partition coefficient
  
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
  Cf <- state[3]
  Ca <- state[4]
  Cpuf <- state[5]
  
  dCsdt <- - f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw
  dCwdt <- - ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) - 
    ko * Af * L / Vw * (Cw - Cf / Kf) # [ng/L]
  dCfdt <- ko * Af / Vf * (Cw - Cf / Kf) # Cw = [ng/L], Cf = [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf) # Ca = [ng/L]
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf) # Ca = [ng/L], Cpuf = [ng/L]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dCfdt, dCadt, dCpufdt)))
}

# Define the number of PCBs to optimize
n <- 2  # Now testing for PCB1 and PCB2
parms_list <- vector("list", n)  # Initialize list for parameters
optimization_results <- vector("list", n)  # Initialize list for storing results

# Define initial parameters
init_parms <- c(
  ro = 500,   # Initial guess for ro
  ko = 10,    # Initial guess for ko
  kdf = 2,    # Initial guess for kdf
  kds = 0.5,  # Initial guess for kds
  f = 0.5,    # Initial guess for f
  ka = 100    # Initial guess for ka
)

# Define the objective function
objective_function <- function(parms, pcb_index, obs.data) {
  
  # Update parameter list for current PCB
  parms_list[[pcb_index]] <<- list(
    ro = parms[1],
    ko = parms[2],
    kdf = parms[3],
    kds = parms[4],
    f = parms[5],
    ka = parms[6]
  )
  
  # Constants and initial conditions
  ms <- 10      # [g]
  M <- 0.1      # kg/L solid-water ratio
  Vw <- 100     # [cm3]
  Va <- 125     # [cm3]
  Vf <- 0.000000069  # L/cm SPME
  Vpuf <- 29    # cm3 volume of PUF
  
  # Extract initial concentrations and time range
  Ct <- obs.data[25, paste0("PCB", pcb_index)]  # PCB concentration at t=25
  Cs0 <- Ct * M * 1000  # Convert to ng/L
  cinit <- c(Cs = Cs0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
  
  # Extract time range for ODE
  time_range <- obs.data$time[1:12]
  
  # Generate a more dynamic mock model for testing (using parms)
  modeled_mf <- parms[5] * sin(time_range / 3)   # Sinusoidal model
  modeled_mpuf <- parms[6] * cos(time_range / 3)  # Cosine model
  
  # Compute observed means
  observed_means <- obs.data %>%
    filter(sample != "sed") %>%
    group_by(time) %>%
    summarise(
      mf_PCB = mean(get(paste0("PCB", pcb_index)), na.rm = TRUE),
      mpuf_PCB = mean(get(paste0("PCB", pcb_index)), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Extract observed values
  observed_mf <- observed_means$mf_PCB  
  observed_mpuf <- observed_means$mpuf_PCB  
  
  # Compute error
  error_mf <- sum((observed_mf - modeled_mf)^2, na.rm = TRUE)
  error_mpuf <- sum((observed_mpuf - modeled_mpuf)^2, na.rm = TRUE)
  total_error <- error_mf + error_mpuf
  
  return(total_error)
}

# Run optimization for each PCB
for (i in 1:n) {
  result <- optim(
    par = init_parms,  
    fn = function(parms) objective_function(parms, i, obs.data),  
    method = "L-BFGS-B",
    lower = c(100, 1, 0.5, 0, 0.1, 50),
    upper = c(1000, 50, 5, 1, 1, 500)
  )
  
  # Store results
  optimization_results[[i]] <- result
}

# Print all optimization results
for (i in 1:n) {
  print(paste("Optimization result for PCB", i, ":"))
  print(optimization_results[[i]])
}
