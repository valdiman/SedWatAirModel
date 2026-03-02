# full_model_avg_obs_with_ko.R
library(deSolve)
library(dplyr)

# -------------------------
# User inputs / experiment constants
# -------------------------
Vw_L   <- 0.100   # L
# Vpw_L will be computed from sediment geometry
Va_L   <- 0.125   # L
Vpuf_L <- 0.029   # L

Aaw  <- 20    # cm2 air-water area
Aws  <- 30    # cm2 sediment-water area
Apuf <- 7.07  # cm2 PUF area

pcb.ind <- "PCB_32"

# -------------------------
# Read & organize data
# -------------------------
exp.data <- read.csv("Data/01_NS_SPME_PUF.csv", stringsAsFactors = FALSE, check.names = FALSE)
pcbi <- exp.data[, c("ID", "Group", "time", "Sample_medium", pcb.ind)]

pcbi.spme.control <- pcbi %>%
  filter(ID == "NBH_NS", Group == "Control", Sample_medium == "SPME") %>%
  rename(mf = !!sym(pcb.ind))

pcbi.puf.control <- pcbi %>%
  filter(ID == "NBH_NS", Group == "Control", Sample_medium == "PUF") %>%
  rename(mpuf = !!sym(pcb.ind))

# -------------------------
# Aggregate replicates by time (means) and combine
# -------------------------
pcbi.spme.mean <- pcbi.spme.control %>%
  group_by(time) %>%
  summarise(mf_control = mean(mf, na.rm = TRUE),
            mf_n = sum(!is.na(mf)),
            mf_sd = ifelse(mf_n > 1, sd(mf, na.rm = TRUE), NA_real_),
            .groups = "drop")

pcbi.puf.mean <- pcbi.puf.control %>%
  group_by(time) %>%
  summarise(mpuf_control = mean(mpuf, na.rm = TRUE),
            mpuf_n = sum(!is.na(mpuf)),
            mpuf_sd = ifelse(mpuf_n > 1, sd(mpuf, na.rm = TRUE), NA_real_),
            .groups = "drop")

pcb_combined_control_avg <- full_join(pcbi.spme.mean, pcbi.puf.mean, by = "time") %>%
  arrange(time)

if (!any(pcb_combined_control_avg$time == 0)) {
  pcb_combined_control_avg <- bind_rows(
    data.frame(time = 0, mf_control = 0, mf_n = 1, mf_sd = NA_real_,
               mpuf_control = 0, mpuf_n = 1, mpuf_sd = NA_real_),
    pcb_combined_control_avg
  ) %>% arrange(time)
} else {
  pcb_combined_control_avg <- pcb_combined_control_avg %>%
    mutate(mf_control  = ifelse(time == 0 & is.na(mf_control), 0, mf_control),
           mpuf_control = ifelse(time == 0 & is.na(mpuf_control), 0, mpuf_control))
}

obs_times <- sort(unique(pcb_combined_control_avg$time))

# -------------------------
# Congener properties & partitioning
# -------------------------
pc <- read.csv("Data/05_PCP.csv", stringsAsFactors = FALSE, check.names = FALSE)
pc_row <- pc[pc$congener == pcb.ind, ]
if (nrow(pc_row) == 0) stop("Congener not found in Data/05_PCP.csv: ", pcb.ind)

MW.pcb <- pc_row$MW
Kow <- pc_row$Kow
dUow <- pc_row$dUow
Kaw <- pc_row$Kaw
dUaw <- pc_row$dUaw
Koa <- pc_row$Koa
E <- pc_row$E; S <- pc_row$S; A <- pc_row$A; B <- pc_row$B; V <- pc_row$V

MH2O <- 18.0152; MCO2 <- 44.0094; R <- 8.3144
Tst.1 <- 273.15 + 25; Tw.1 <- 273.15 + 20

Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 - 1 / Tst.1))
Kaw.t <- Kaw * exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1

# -------------------------
# PUF & SPME params (use your Af value)
# -------------------------
Vpuf_cm3 <- Vpuf_L * 1000
d_puf <- 0.0213 * 100^3
Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) * d_puf

Af <- 0.007226      # <- keep the Af number you set (units cm2, for the full fiber length if L=1)
Vf_cm3_per_cm <- 0.000000069 * 1000
L <- 1
Vf_tot <- Vf_cm3_per_cm * L   # cm3
Kf <- 10^(1.06 * log10(Kow.t) - 1.16)

# -------------------------
# Diffusion / transfer precompute
# -------------------------
D.water.air <- 0.2743615
D.co2.w <- 1.67606E-05
D.pcb.air <- D.water.air * (MW.pcb / MH2O)^(-0.5)
D.pcb.water <- D.co2.w * (MW.pcb / MCO2)^(-0.5)
v.H2O <- 0.010072884
V.water.air <- 0.003
V.co2.w <- 4.1e-2
SC.pcb.w <- v.H2O / D.pcb.water

bl <- 0.1
kpw <- D.pcb.water * 60 * 60 * 24 / bl

Kaw.a <- V.water.air * (D.pcb.air / D.water.air)^(0.67)
Kaw.w <- V.co2.w * (SC.pcb.w / 600)^(-0.5)
kaw.o <- (1 / (Kaw.a * Kaw.t) + (1 / Kaw.w))^-1
kaw <- kaw.o * 100 * 60 * 60 * 24

# Kd via Abraham Koc
logKoc <- 1.1 * E - 0.72 * S + 0.15 * A - 1.98 * B + 2.28 * V + 0.14
Koc <- 10^(logKoc)
foc <- 0.03
Kd <- Koc * foc   # L/kg

# -------------------------
# Sediment geometry -> Vpw
# -------------------------
n  <- 0.42
ds <- 1540      # g / L
ms_g <- 10      # g sediment in vial
M  <- ds * (1 - n) / n
Vpw_cm3 <- ms_g / M * 1000
Vpw_L <- Vpw_cm3 / 1000

cat("Computed Vpw_cm3 =", signif(Vpw_cm3,6), "cm3; Vpw_L =", signif(Vpw_L,6), "L\n")

# -------------------------
# Werner desorption single-domain
# -------------------------
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

# -------------------------
# Build parms_base (include Kf and correct Vf name)
# -------------------------
bulk_conc <- read.csv("Data/03_NBH_SedimentPCB.csv", stringsAsFactors = FALSE, check.names = FALSE)
if (!(pcb.ind %in% colnames(bulk_conc))) stop("Congener not found in 03_NBH_SedimentPCB.csv: ", pcb.ind)
Ct <- mean(bulk_conc[[pcb.ind]], na.rm = TRUE)

parms_base <- list(
  ms   = ms_g,
  Kd   = Kd,
  kpw  = kpw,
  kaw  = kaw,
  kb   = 0,
  ro   = 400,
  ko   = 3,
  Kpuf = Kpuf,
  Kaw.t= Kaw.t,
  Af   = Af,
  L    = L,
  Vf_tot = Vf_tot,
  Apuf = Apuf,
  Aaw  = Aaw,
  Aws  = Aws,
  Vpw  = Vpw_cm3,
  Vpw_L = Vpw_L,
  Vw_L = Vw_L,
  Va_L = Va_L,
  Vpuf_L = Vpuf_L,
  Vw  = Vw_L * 1000,
  Va  = Va_L * 1000,
  Vpuf = Vpuf_cm3,
  Vf_tot_cm3 = Vf_tot,
  Kf   = Kf,
  D_r2 = 3.5e-6,
  M_sed_init = NA,
  Nterms_des = 200
)

M_sed_init <- Ct * parms_base$ms
parms_base$M_sed_init <- M_sed_init

# initial conditions
Cs_init <- Ct
Cpw_init <- 0
cinit_pcb <- c(Cs = Cs_init, Cpw = Cpw_init, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)

cat("Initial Cs (ng/g):", signif(Cs_init,6), "  initial Cpw (ng/L):", signif(Cpw_init,6), "\n")
cat("M_sed_init (ng):", signif(M_sed_init,6), "\n")

# -------------------------
# Robust ODE function & run helper
# -------------------------
rtm.PCB_safe <- function(t, state, parms) {
  if (length(state) != 6) stop("rtm.PCB_safe: state must be length 6")
  required <- c("Vpw", "Vw", "Va", "Vpuf", "Vf_tot_cm3", "Af", "Apuf", "Aws", "Aaw",
                "L", "kpw", "kaw", "kb", "ko", "Kaw.t", "Kf", "Kpuf", "ro", "ms", "M_sed_init", "D_r2", "Nterms_des")
  miss <- setdiff(required, names(parms))
  if (length(miss) > 0) stop("rtm.PCB_safe: missing parms: ", paste(miss, collapse = ", "))
  
  Vpw_cm3 <- as.numeric(parms$Vpw)
  Vw_cm3  <- as.numeric(parms$Vw)
  Va_cm3  <- as.numeric(parms$Va)
  Vpuf_cm3<- as.numeric(parms$Vpuf)
  Vf_tot  <- as.numeric(parms$Vf_tot_cm3)
  
  Af_val  <- as.numeric(parms$Af)
  Apuf_val<- as.numeric(parms$Apuf)
  Aws_val <- as.numeric(parms$Aws)
  Aaw_val <- as.numeric(parms$Aaw)
  L_val   <- as.numeric(parms$L)
  
  kpw_val <- as.numeric(parms$kpw)
  kaw_val <- as.numeric(parms$kaw)
  kb_val  <- as.numeric(parms$kb)
  ko_val  <- as.numeric(parms$ko)
  Kawt    <- as.numeric(parms$Kaw.t)
  Kf_val  <- as.numeric(parms$Kf)
  Kpuf_val<- as.numeric(parms$Kpuf)
  ro_val  <- as.numeric(parms$ro)
  ms_g    <- as.numeric(parms$ms)
  
  Cs   <- as.numeric(state[1])
  Cpw  <- as.numeric(state[2])
  Cw   <- as.numeric(state[3])
  Cf   <- as.numeric(state[4])
  Ca   <- as.numeric(state[5])
  Cpuf <- as.numeric(state[6])
  
  Cpw_cm3  <- Cpw  / 1000
  Cw_cm3   <- Cw   / 1000
  Cf_cm3   <- Cf   / 1000
  Ca_cm3   <- Ca   / 1000
  Cpuf_cm3 <- Cpuf / 1000
  
  des <- desorption_single(t_days = t, D_r2 = parms$D_r2, mtot_ng = parms$M_sed_init, Nterms = parms$Nterms_des)
  dm_dt_ng_per_day <- as.numeric(des$dm_dt[1])
  
  dCsdt <- - dm_dt_ng_per_day / ms_g
  
  dCpwdt <- (dm_dt_ng_per_day / Vpw_cm3) -
    kpw_val * Aws_val / Vpw_cm3 * (Cpw_cm3 - Cw_cm3) -
    kb_val * Cpw_cm3
  
  dCwdt <- kpw_val * Aws_val / Vw_cm3 * (Cpw_cm3 - Cw_cm3) -
    kaw_val * Aaw_val / Vw_cm3 * (Cw_cm3 - Ca_cm3 / Kawt) -
    ko_val * Af_val * L_val / Vw_cm3 * (Cw_cm3 - Cf_cm3 / Kf_val)
  
  dCfdt <- ko_val * Af_val / Vf_tot * (Cw_cm3 - Cf_cm3 / Kf_val)
  
  dCadt <- kaw_val * Aaw_val / Va_cm3 * (Cw_cm3 - Ca_cm3 / Kawt) -
    ro_val * Apuf_val / Va_cm3 * (Ca_cm3 - Cpuf_cm3 / Kpuf_val)
  
  dCpufdt <- ro_val * Apuf_val / Vpuf_cm3 * (Ca_cm3 - Cpuf_cm3 / Kpuf_val)
  
  derivs <- c(
    as.numeric(dCsdt),
    as.numeric(dCpwdt * 1000),
    as.numeric(dCwdt * 1000),
    as.numeric(dCfdt * 1000),
    as.numeric(dCadt * 1000),
    as.numeric(dCpufdt * 1000)
  )
  
  if (length(derivs) != 6) stop("rtm.PCB_safe: derivative vector length != 6")
  if (!is.numeric(derivs)) stop("rtm.PCB_safe: derivatives not numeric")
  
  return(list(derivs))
}

run_model_PCB_safe <- function(pars, times_run = obs_times, cinit = cinit_pcb) {
  if (length(cinit) != 6) stop("run_model_PCB_safe: cinit must be length 6")
  need <- c("Vpw","Vw","Va","Vpuf","Vf_tot_cm3","Af","Apuf","Aws","Aaw","L",
            "kpw","kaw","kb","ko","Kaw.t","Kf","Kpuf","ro","ms","M_sed_init","D_r2","Nterms_des")
  miss2 <- setdiff(need, names(pars))
  if (length(miss2) > 0) stop("run_model_PCB_safe: missing parms: ", paste(miss2, collapse = ", "))
  
  out <- ode(y = cinit, times = times_run, func = rtm.PCB_safe, parms = pars, atol = 1e-8, rtol = 1e-8)
  df <- as.data.frame(out)
  colnames(df) <- c("time","Cs","Cpw","Cw","Cf","Ca","Cpuf")
  df$mf_model <- (df$Cf / 1000) * pars$Vf_tot_cm3
  df$mpuf_model <- (df$Cpuf / 1000) * pars$Vpuf
  return(df)
}

# -------------------------
# Fit D_r2, ro, ko to averaged observations (log1p SSE)
# -------------------------
parms_PCB <- parms_base

obj_fit_Dr2_ro_ko_avg <- function(par_log) {
  D_r2_try <- exp(par_log[1])
  ro_try   <- exp(par_log[2])
  ko_try   <- exp(par_log[3])
  
  p <- parms_PCB
  p$D_r2 <- D_r2_try
  p$ro   <- ro_try
  p$ko   <- ko_try
  
  out <- tryCatch(run_model_PCB_safe(p, times_run = obs_times, cinit = cinit_pcb), error = function(e) NULL)
  if (is.null(out)) return(1e30)
  
  cmp <- out %>%
    select(time, mf_model, mpuf_model) %>%
    right_join(pcb_combined_control_avg %>% select(time, mf_control, mpuf_control), by = "time")
  if (nrow(cmp) == 0) return(1e30)
  
  sse_mf   <- sum((log1p(cmp$mf_model)  - log1p(cmp$mf_control))^2, na.rm = TRUE)
  sse_mpuf <- sum((log1p(cmp$mpuf_model) - log1p(cmp$mpuf_control))^2, na.rm = TRUE)
  sse <- sse_mf + sse_mpuf
  return(sse)
}

init3 <- c(log(parms_PCB$D_r2), log(parms_PCB$ro), log(parms_PCB$ko))
lower3 <- c(log(1e-12), log(1), log(1e-3))
upper3 <- c(log(1e-1), log(1e4), log(1e4))

res3_avg <- optim(par = init3, fn = obj_fit_Dr2_ro_ko_avg, method = "L-BFGS-B",
                  lower = lower3, upper = upper3, control = list(maxit = 2000))

par_hat <- exp(res3_avg$par)
D_r2_hat_new <- par_hat[1]
ro_hat_new   <- par_hat[2]
ko_hat_new   <- par_hat[3]

cat("Joint fitted parameters (averaged obs):\n")
cat(" D/r^2 (day^-1):", signif(D_r2_hat_new,6), "\n")
cat(" ro (cm/day):    ", signif(ro_hat_new,6), "\n")
cat(" ko (cm/day):    ", signif(ko_hat_new,6), "\n")
cat(" SSE (log-space):", signif(res3_avg$value,6), "\n\n")

# -------------------------
# Final run & comparison
# -------------------------
parms_PCB$D_r2 <- D_r2_hat_new
parms_PCB$ro   <- ro_hat_new
parms_PCB$ko   <- ko_hat_new

df_best_pcb <- run_model_PCB_safe(parms_PCB, times_run = obs_times, cinit = cinit_pcb)

comparison_pcb <- df_best_pcb %>%
  select(time, mf_model, mpuf_model) %>%
  right_join(pcb_combined_control_avg %>% select(time, mf_control, mpuf_control), by = "time") %>%
  arrange(time)

cat("Model vs Observed (mf, mpuf in ng) using averaged obs:\n")
print(comparison_pcb)

# mass-balance diagnostics
df_diag_pcb <- df_best_pcb %>%
  mutate(ms_mass = Cs * parms_PCB$ms,
         mpw = Cpw * parms_PCB$Vpw / 1000,
         mw  = Cw  * parms_PCB$Vw_L,
         ma  = Ca  * parms_PCB$Va_L,
         mpuf = mpuf_model,
         mf   = mf_model) %>%
  mutate(mt = ms_mass + mpw + mw + ma + mpuf + mf) %>%
  select(time, Cs, ms_mass, mpw, mw, ma, mf, mpuf, mt)

cat("\nMass-balance diagnostics (ng):\n")
print(df_diag_pcb)