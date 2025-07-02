# === 1. Packages & Libraries ===
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

# === 2. Data Input ===
obs.data <- read.csv("Data/AVL_S_data_long.csv", check.names = FALSE)
pcp.data <- read.csv("Data/AVL_S_PCP.csv")

# === 3. Select Target PCB Congener ===
pcb.name <- "PCB19"
obs.data.pcbi <- obs.data[, c("sampler", "time", pcb.name)]
pcp.data.pcbi <- pcp.data[pcp.data$congener == pcb.name, ]

# === 4. Organize Observations ===
spme <- obs.data.pcbi %>% filter(sampler == "SPME")
puf  <- obs.data.pcbi %>% filter(sampler == "PUF") %>% mutate(across(starts_with("PCB"), ~ . * 10))
sed  <- obs.data.pcbi %>% filter(sampler == "sed")

# === 5. Reactive Transport Model Function ===
rtm.PCBi <- function(t, state, parms) {
  # Constants
  MH2O <- 18.0152; MCO2 <- 44.0094; R <- 8.3144
  Tst <- 25; Tst.K <- 273.15 + Tst
  Tw  <- 20; Tw.K  <- 273.15 + Tw
  
  # Bioreactor Geometry
  Vw <- 100; Va <- 125; Aaw <- 20; Aws <- 30
  
  # Congener-Specific Properties
  MW.pcb <- pcp.data.pcbi$MW
  Kaw <- pcp.data.pcbi$Kaw
  dUaw <- pcp.data.pcbi$dUaw
  Kaw.t <- Kaw * exp(-dUaw / R * (1 / Tw.K - 1 / Tst.K)) * Tw.K / Tst.K
  
  Kow <- pcp.data.pcbi$Kow
  dUow <- pcp.data.pcbi$dUow
  Kow.t <- Kow * exp(-dUow / R * (1 / Tw.K - 1 / Tst.K))
  
  Koa <- pcp.data.pcbi$Koa
  
  # PUF Constants
  Apuf <- 7.07; Vpuf <- 29
  d <- 0.0213 * 100^3
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) * d
  
  # SPME Constants
  Aspme <- 0.138; Vspme <- 0.000000069 * 1000; Lspme <- 1
  Kspme <- 10^(1.06 * log10(Kow.t) - 1.16)
  
  # Physical Conditions
  D.water.air <- 0.2743615
  D.co2.w <- 1.67606E-05
  D.pcb.air <- D.water.air * (MW.pcb / MH2O)^(-0.5)
  D.pcb.water <- D.co2.w * (MW.pcb / MCO2)^(-0.5)
  v.H2O <- 0.010072884
  V.water.air <- 0.003
  V.co2.w <- 4.1e-2
  SC.pcb.w <- v.H2O / D.pcb.water
  
  # Mass Transfer Coefficients
  Kaw.a <- V.water.air * (D.pcb.air / D.water.air)^0.67
  Kaw.w <- V.co2.w * (SC.pcb.w / 600)^(-0.5)
  kaw.o <- (1 / (Kaw.a * Kaw.t) + 1 / Kaw.w)^(-1) * 100 * 60 * 60 * 24
  
  # Parameters
  with(parms, {
    Cs <- state[1]; Cw <- state[2]; Cspme <- state[3]; Ca <- state[4]; Cpuf <- state[5]
    
    dCsdt    <- -f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw
    dCwdt    <- -ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
      kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
      ko * Aspme * Lspme / Vw * (Cw - Cspme / Kspme) -
      kb * Cw
    dCspmedt <- ko * Aspme / Vspme * (Cw - Cspme / Kspme)
    dCadt    <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) - ro * Apuf / Va * (Ca - Cpuf / Kpuf)
    dCpufdt  <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf)
    
    list(c(dCsdt, dCwdt, dCspmedt, dCadt, dCpufdt))
  })
}

# === 6. Initial Conditions ===
Ct <- sed %>% pull(!!sym(pcb.name)) # [ng/g]
M <- 0.1 # [kg/L]
Cs0 <- Ct * M * 1000 # [ng/L]
cinit <- c(Cs = Cs0, Cw = 0, Cspme = 0, Ca = 0, Cpuf = 0)

# === 7. Format Observed Data for Fitting ===
obs.data.pcbi.2 <- obs.data.pcbi %>%
  filter(sampler %in% c("SPME", "PUF")) %>%
  group_by(time, sampler) %>%
  summarise(!!sym(pcb.name) := mean(!!sym(pcb.name), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sampler, values_from = all_of(pcb.name)) %>%
  mutate(PUF = as.numeric(PUF) * 10, SPME = as.numeric(SPME))

# === 8. Define Time Points ===
t.1 <- unique(obs.data.pcbi.2$time)

# === 9. Cost Function ===
cost_func <- function(parms) {
  out <- ode(y = cinit, times = t.1, func = rtm.PCBi, parms = list(
    ro = parms["ro"], ko = parms["ko"], kdf = parms["kdf"],
    kds = parms["kds"], f = parms["f"], ka = parms["ka"], kb = 0
  ))
  out_df <- as.data.frame(out)
  
  pred_spme <- out_df$Cspme * 6.9e-8
  pred_puf  <- out_df$Cpuf * 29 / 1000
  
  obs <- obs.data.pcbi.2 %>% filter(time %in% out_df$time)
  c(pred_spme - obs$SPME, pred_puf - obs$PUF)
}

# === 10. Fit Model ===
par_guess <- c(ro = 540, ko = 10, kdf = 1.9, kds = 0.001, f = 0.8, ka = 150)
bounds <- list(lower = c(ro = 1, ko = 1e-3, kdf = 1e-3, kds = 1e-6, f = 0, ka = 1),
               upper = c(ro = 1e4, ko = 1e3, kdf = 10, kds = 1, f = 1, ka = 1e3))

fit <- modFit(f = cost_func, p = par_guess, lower = bounds$lower, upper = bounds$upper)
summary(fit)

# === 11. Final Model Run ===
best_parms <- as.list(fit$par); best_parms$kb <- 0
final_out_df <- ode(y = cinit, times = t.1, func = rtm.PCBi, parms = best_parms) %>%
  as.data.frame() %>%
  mutate(
    mspme = Cspme * 6.9e-8,
    mpuf = Cpuf * 29 / 1000,
    Mt = (Cs * 10 / (0.1 * 1000)) + (Cw * 0.1) + mspme + (Ca * 0.125) + mpuf
  )

# === 12. Fine Prediction ===
fine_time <- seq(min(t.1), max(t.1), by = 0.5)
fine_out_df <- ode(y = cinit, times = fine_time, func = rtm.PCBi, parms = best_parms) %>%
  as.data.frame() %>%
  mutate(
    mspme = Cspme * 6.9e-8,
    mpuf = Cpuf * 29 / 1000
  )

# === 13. Model Fit Evaluation ===
obs_pred <- obs.data.pcbi.2 %>%
  left_join(final_out_df, by = "time") %>%
  mutate(
    pred_spme = Cspme * 6.9e-8,
    pred_puf = Cpuf * 29 / 1000
  )

metrics <- function(obs, pred) {
  data.frame(
    RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    MAE = mean(abs(obs - pred), na.rm = TRUE),
    R2 = cor(obs, pred, use = "complete.obs")^2,
    NSE = 1 - sum((obs - pred)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  )
}

spme_metrics <- metrics(obs_pred$SPME, obs_pred$pred_spme)
puf_metrics  <- metrics(obs_pred$PUF, obs_pred$pred_puf)

# === 14. Plots ===
plot_func <- function(data, obs_col, pred_col, title, ylab, col_obs, col_pred, shape) {
  ggplot() +
    geom_line(data = fine_out_df, aes(time, !!sym(pred_col), color = "Predicted"), linewidth = 0.8) +
    geom_point(data = obs.data.pcbi.2, aes(time, !!sym(obs_col), color = "Observed"), size = 3, shape = shape) +
    labs(title = title, x = "Time (days)", y = ylab) +
    scale_color_manual(values = c("Observed" = col_obs, "Predicted" = col_pred)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
}

spme_plot <- plot_func(obs_pred, "SPME", "mspme", "SPME: Model Fit", "Mass (ng/cm)", "red", "blue", 19) +
  annotate("text", x = max(t.1)*0.7, y = max(obs_pred$SPME, na.rm = TRUE)*0.9,
           label = paste0("R² = ", round(spme_metrics$R2, 3), "\nRMSE = ", round(spme_metrics$RMSE, 3))
           hjust = 0, size = 4)

puf_plot <- plot_func(obs_pred, "PUF", "mpuf", "PUF: Model Fit", "Mass (ng/puf)", "green4", "orange2", 17) +
  annotate("text", x = max(t.1)*0.7, y = max(obs_pred$PUF, na.rm = TRUE)*0.9,
           label = paste0("R² = ", round(puf_metrics$R2, 4), "\nRMSE = ", round(puf_metrics$RMSE, 3)),
           hjust = 0, size = 4)

# === 15. Display Plots ===
print(spme_plot)
print(puf_plot)
