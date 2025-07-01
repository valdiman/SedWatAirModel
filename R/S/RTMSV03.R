library(FME)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Read data ---------------------------------------------------------------
obs.data <- read.csv("Data/AVL_S_data_long.csv", check.names = FALSE)
pcp.data <- read.csv("Data/AVL_S_PCP.csv")

# --- Extract individual PCBi -------------------------------------------------
pcb.name <- "PCB19"
obs.data.pcbi <- obs.data[, c("sampler", "time", pcb.name)]
pcp.data.pcbi <- pcp.data[pcp.data$congener == pcb.name, ]

# --- Organize data -----------------------------------------------------------
spme <- obs.data.pcbi %>% filter(sampler == "SPME")
puf <- obs.data.pcbi %>% filter(sampler == "PUF") %>% mutate(across(starts_with("PCB"), ~ . * 10))
sed <- obs.data.pcbi %>% filter(sampler == "sed")

# Prepare observed data in long format ------------------------------------------------
# Combine all samplers to one dataframe, rename column to 'value', and add 'var' column matching model output names
obs <- bind_rows(
  spme %>% select(time, value = !!sym(pcb.name)) %>% mutate(var = "Cspme"),
  puf %>% select(time, value = !!sym(pcb.name)) %>% mutate(var = "Cpuf"),
  sed %>% select(time, value = !!sym(pcb.name)) %>% mutate(var = "Cs")
) %>% arrange(var, time)

# --- Define the reactive transport model -----------------------------------
rtm.PCBi <- function(t, state, parms) {
  # Constants
  MH2O <- 18.0152
  MCO2 <- 44.0094
  MW.pcb <- pcp.data.pcbi$MW
  R <- 8.3144
  Tst <- 25
  Tst.1 <- 273.15 + Tst
  Tw <- 20
  Tw.1 <- 273.15 + Tw
  
  # Volumes and areas (cm3, cm2)
  Vw <- 100
  Va <- 125
  Aaw <- 20
  Aws <- 30
  Vs <- 10
  
  # PCB-specific constants and temperature adjustments
  Kaw <- pcp.data.pcbi$Kaw
  dUaw <- pcp.data.pcbi$dUaw
  Kaw.t <- Kaw * exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  
  Kow <- pcp.data.pcbi$Kow
  dUow <- pcp.data.pcbi$dUow
  Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 - 1 / Tst.1))
  
  Koa <- pcp.data.pcbi$Koa
  
  # PUF constants
  Apuf <- 7.07
  Vpuf <- 29
  d <- 0.0213 * 100^3
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774)
  Kpuf <- Kpuf * d
  
  # SPME fiber constants
  Aspme <- 0.138
  Vspme <- 0.000000069 * 1000
  Lspme <- 1
  Kspme <- 10^(1.06 * log10(Kow.t) - 1.16)
  
  # Air & water physical conditions
  D.water.air <- 0.2743615
  D.co2.w <- 1.67606E-05
  D.pcb.air <- D.water.air * (MW.pcb / MH2O)^(-0.5)
  D.pcb.water <- D.co2.w * (MW.pcb / MCO2)^(-0.5)
  
  v.H2O <- 0.010072884
  V.water.air <- 0.003
  V.co2.w <- 4.1e-2
  SC.pcb.w <- v.H2O / D.pcb.water
  
  # Air-water mass transfer coefficients
  Kaw.a <- V.water.air * (D.pcb.air / D.water.air)^(0.67)
  Kaw.w <- V.co2.w * (SC.pcb.w / 600)^(-0.5)
  kaw.o <- (1 / (Kaw.a * Kaw.t) + (1 / Kaw.w))^-1
  kaw.o <- kaw.o * 100 * 60 * 60 * 24
  
  # Parameters
  ka <- parms$ka
  kdf <- parms$kdf
  kds <- parms$kds
  f <- parms$f
  kb <- parms$kb
  ko <- parms$ko
  ro <- parms$ro
  Mt <- parms$Mt
  
  # State variables
  Cw <- state[1]
  Cspme <- state[2]
  Ca <- state[3]
  Cpuf <- state[4]
  
  # Compute Cs from mass balance
  Cs <- (Mt - Cw * Vw - Cspme * Vspme - Ca * Va - Cpuf * Vpuf) / Vs
  
  # Derivatives
  dCwdt <- -ka * Cw +
    f * kdf * Cs + (1 - f) * kds * Cs -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    ko * Aspme * Lspme / Vw * (Cw - Cspme / Kspme) -
    kb * Cw
  
  dCspmedt <- ko * Aspme / Vspme * (Cw - Cspme / Kspme)
  
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf)
  
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf)
  
  return(list(c(dCwdt, dCspmedt, dCadt, dCpufdt), Cs = Cs))
}

# --- Initial total mass Mt from sediment concentration -------------------
Ct <- sed %>% filter(time == 0) %>% pull(!!sym(pcb.name))  # ng/g sediment at time zero
M <- 0.1  # kg/L sediment density (or sediment water content?)
Vs <- 10  # cm3 sediment volume in reactor (consistent with model)
Cs0 <- Ct * M * 1000   # ng/L converted from ng/g * kg/L * 1000 (L/m3)
Mt <- Cs0 * Vs         # total ng mass

# Initial state vector (no initial concentrations except sediment Cs calculated inside model)
state <- c(Cw = 0, Cspme = 0, Ca = 0, Cpuf = 0)

# Time sequence for model simulation
time <- sort(unique(obs$time))

# --- Initial parameter guesses for optimization ---------------------------
parms <- list(
  ka = 0.2,
  kdf = 0.1,
  kds = 0.05,
  f = 0.3,
  ko = 0.02,
  ro = 0.015,
  kb = 0.01,
  Mt = Mt
)

# --- Cost function for parameter fitting -----------------------------------
costFunc <- function(p_vec) {
  p <- as.list(p_vec)
  p$Mt <- parms$Mt
  
  sim <- ode(y = state, times = time, func = rtm.PCBi, parms = p)
  sim <- as.data.frame(sim) %>%
    pivot_longer(cols = -time, names_to = "var", values_to = "model")
  
  obs.joined <- left_join(obs, sim, by = c("time", "var")) %>%
    filter(!is.na(value), !is.na(model))
  
  residuals <- obs.joined$value - obs.joined$model
  return(residuals)
}


# --- Prepare numeric parameter vector for optimization (excluding Mt) ------
p_fit <- unlist(parms[names(parms) != "Mt"])

# --- Run optimization ------------------------------------------------------
fit <- modFit(f = costFunc, p = p_fit)

# --- View optimized parameters ---------------------------------------------
print(fit$par)

# --- Run final simulation with optimized parameters ------------------------
p_opt <- as.list(fit$par)
p_opt$Mt <- parms$Mt

out <- ode(y = state, times = time, func = rtm.PCBi, parms = p_opt) %>%
  as.data.frame() %>%
  pivot_longer(cols = -time, names_to = "var", values_to = "model")

# --- Plot observed vs modeled -----------------------------------------------
library(ggplot2)

ggplot() +
  geom_point(data = obs, aes(time, value, color = var), shape = 1) +
  geom_line(data = out, aes(time, model, color = var)) +
  labs(title = paste0("PCB Model Fit for ", pcb.name),
       x = "Time", y = "Concentration (ng/L or ng/g)") +
  theme_minimal()
