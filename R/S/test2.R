
library(deSolve)

library(dplyr)

puf_summary <- puf %>%
  group_by(time) %>%
  summarize(mean_PUF = mean(PCB19))


library(deSolve)
library(dplyr)

# Your model function, params as named vector
sed_puf_model <- function(t, state, params) {
  Cs <- state[1]
  Cpuf <- state[2]
  
  k_sed_to_puf <- params["k_sed_to_puf"]
  k_puf_to_sed <- 0  # irreversible
  
  dCsdt <- - k_sed_to_puf * Cs
  dCpufdt <- k_sed_to_puf * Cs
  
  list(c(dCsdt, dCpufdt))
}

# Objective function using interpolation for predicted PUF at observed times
objective_fn <- function(params, times, observed_puf, Cs0, Cpuf0) {
  # Remove duplicate times for ODE
  times_unique <- sort(unique(times))
  
  state <- c(Cs = Cs0, Cpuf = Cpuf0)
  out <- ode(y = state, times = times_unique, func = sed_puf_model, parms = params)
  out_df <- as.data.frame(out)
  
  # Interpolate predicted Cpuf at original observed times
  pred_puf <- approx(out_df$time, out_df$Cpuf, xout = times)$y
  
  # Sum of squared errors
  ssq <- sum((observed_puf - pred_puf)^2)
  return(ssq)
}

# Summarize PUF data (average replicates at same times)
puf_summary <- puf %>%
  group_by(time) %>%
  summarize(mean_PUF = mean(PCB19))

# Initial conditions
Ct <- sed$PCB19 # your sediment PCB19 conc (ng/g)
M <- 0.1        # solid-water ratio (kg/L)
Cs0 <- Ct * M * 1000 # Convert to ng/L
Cpuf0 <- 0

times <- puf_summary$time
observed_puf <- puf_summary$mean_PUF

# Initial guess for rate constant
init_params <- c(k_sed_to_puf = 0.5)

# Fit with bounds (you can adjust upper bound)
fit <- optim(par = init_params, fn = objective_fn, times = times, observed_puf = observed_puf, Cs0 = Cs0, Cpuf0 = Cpuf0,
             method = "L-BFGS-B", lower = 0, upper = 10)

print(fit$par)  # Estimated k_sed_to_puf

# Plotting model fit
params_fit <- fit$par
state <- c(Cs = Cs0, Cpuf = Cpuf0)
times_fine <- seq(min(times), max(times), length.out = 100)
out <- ode(y = state, times = times_fine, func = sed_puf_model, parms = params_fit)
out_df <- as.data.frame(out)

plot(times, observed_puf, pch = 19, col = "red", xlab = "Time", ylab = "PUF PCB19 concentration (ng/L)")
lines(out_df$time, out_df$Cpuf, col = "blue", lwd = 2)
legend("topleft", legend = c("Observed PUF", "Model Prediction"), col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1))
