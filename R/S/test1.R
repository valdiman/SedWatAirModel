

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


rtm.PCB = function(t, state, parms){

Vw <- 100 # cm3 water volume
Va <- 125 # cm3 headspace volumne
Aaw <- 20 # cm2

Af <- 0.138 # cm2/cm SPME area
Vf <- 0.000000069 * 1000 # cm3/cm SPME volume/area
L <- 1 # cm SPME length normalization to 1 cm

Apuf <- 7.07 # cm2
Vpuf <- 29 # cm3 volume of PUF

kaw.o <- 101.6757 #[cm/d]
Kaw.t <- 0.008619 # 
Kf <- 4328.106
Kpuf <- 161815.9 # [La/Lpuf]

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

# Initial conditions and run function
{
  # Estimating Cs0
  Ct <- 259.8342356 # ng/g PCB 19 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = Cs0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
parms <- list(ro = 540.409, ko = 10, kdf = 1.9, kds = 0.001, f = 0.8,
              ka = 150) # Input
t <- c(0, 3, 11, 16, 35)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t, func = rtm.PCB, parms = parms)
head(out.1)

# Transform Cf and Cpuf to mass/cm and mass/puf
out.1 <- as.data.frame(out.1)
colnames(out.1) <- c("time", "Cs", "Cw", "Cf", "Ca", "Cpuf")
  
# Calculate Mf and Mpuf based on volumes
ms <- 10 # [g]
M <- 0.1 # kg/L solid-water ratio
Vw <- 100 # [cm3]
Va <- 125 # [cm3]
Vf <- 0.000000069 # L/cm SPME
Vpuf <- 29 # cm3 volume of PUF
out.1$mf <- out.1$Cf * Vf  # [ng/cm]
out.1$mpuf <- out.1$Cpuf * Vpuf / 1000 # [ng/puf]
out.1$Mt <- out.1$Cs * ms / (M * 1000) + out.1$Cw * Vw / 1000 + out.1$Cf * Vf +
  out.1$Ca * Va / 1000 + out.1$Cpuf * Vpuf / 1000


