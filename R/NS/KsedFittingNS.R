# Fit sediment desorption rate (ksed) for PCB transport model.
# NS = non-shaken

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("deSolve")
install.packages("tibble")

# Load libraries
{
  library(deSolve)
  library(dplyr)
  library(tibble)
}

# -------------------------
# Volumes
Vw_L   <- 0.100   # L (100 cm3)
Va_L   <- 0.125   # L (125 cm3)
Vpuf_L <- 0.029   # L (29 cm3)

# Areas and sampler geometry
Aaw  <- 20    # cm2 air-water area
Aws  <- 30    # cm2 sediment-water area
Apuf <- 7.07  # cm2 (PUF area, whichever unit you've used before)

# Non-shaken data
exp.ns.data <- read.csv("Data/01_NS_SPME_PUF.csv", stringsAsFactors = FALSE,
                        check.names = FALSE)

pcb.ind <- "PCB_18+30"

pcbi <- exp.ns.data[, c("ID", "Sample_medium", "Group",
                        "time", "Replicate", pcb.ind)]

# PUF control samples (mpuf)
pcbi.puf.control <- pcbi %>%
  filter(ID == "NBH_NS", Sample_medium == "PUF",
         Group == "Control") %>%
  rename(mpuf_control = !!sym(pcb.ind)) %>%
  select(time, mpuf_control)

# Add t = 0
pcbi_control_0 <- rbind(
  data.frame(time = 0, mpuf_control = 0),
  pcbi.puf.control
)

# Observed PUF by time (mean across replicates)
obs_df <- as_tibble(pcbi_control_0) %>%
  group_by(time) %>%
  summarise(mpuf_obs = mean(mpuf_control, na.rm = TRUE), .groups = "drop") %>%
  arrange(time) %>%
  filter(!is.na(mpuf_obs))

# -------------------------
# 2) Simplified ODE model
#    states: Cs (ng/g), Cpw (ng/L), Cw (ng/L), Ca (ng/L), Cpuf (ng/L)
# -------------------------
rtm_simple <- function(t, state, parms) {
  # ensure state is numeric vector of length 5
  if (length(state) < 5) stop("rtm_simple: state vector must have length 5 (Cs, Cpw, Cw, Ca, Cpuf).")
  # extract states by index to avoid name-mismatch problems
  Cs   <- as.numeric(state[1])   # ng/g
  Cpw  <- as.numeric(state[2])   # ng/L
  Cw   <- as.numeric(state[3])   # ng/L
  Ca   <- as.numeric(state[4])   # ng/L
  Cpuf <- as.numeric(state[5])   # ng/L
  
  # quick checks that parms contain required entries
  required_parms <- c("Vpw_cm3","Vw_L","Va_L","Vpuf_L","Aaw","Aws","Apuf",
                      "Kd","ms","ksed","kpw","kaw","kb","ro","Kpuf","Kaw.t",
                      "Af","Vf_tot")
  missing_p <- setdiff(required_parms, names(parms))
  if (length(missing_p) > 0) {
    stop("rtm_simple: missing parms: ", paste(missing_p, collapse = ", "))
  }
  
  # volumes / conversions
  Vpw_cm3  <- as.numeric(parms$Vpw_cm3)
  Vw_cm3   <- as.numeric(parms$Vw_L)  * 1000
  Va_cm3   <- as.numeric(parms$Va_L)  * 1000
  Vpuf_cm3 <- as.numeric(parms$Vpuf_L) * 1000
  
  # convert to ng/cm3 for flux calcs
  Cpw_cm3  <- Cpw  / 1000
  Cw_cm3   <- Cw   / 1000
  Ca_cm3   <- Ca   / 1000
  Cpuf_cm3 <- Cpuf / 1000
  
  # convert Cs (ng/g) to porewater-equivalent in ng/cm3 using Kd (L/kg)
  # Cs_pw_eq_cm3 = Cs / Kd   (ng/cm3)
  Cs_pw_eq_cm3 <- Cs / parms$Kd
  
  # parameters
  ksed <- as.numeric(parms$ksed)  # 1/day
  kpw  <- as.numeric(parms$kpw)   # cm/day
  kaw  <- as.numeric(parms$kaw)   # cm/day
  kb   <- as.numeric(parms$kb)    # 1/day
  ro   <- as.numeric(parms$ro)    # 1/day
  Kpuf <- as.numeric(parms$Kpuf)
  Kaw.t<- as.numeric(parms$Kaw.t)
  Af   <- as.numeric(parms$Af)
  Vf_tot <- as.numeric(parms$Vf_tot)
  Apuf <- as.numeric(parms$Apuf)
  Aws  <- as.numeric(parms$Aws)
  Aaw  <- as.numeric(parms$Aaw)
  ms_g <- as.numeric(parms$ms)
  
  # ---- ODEs (use ng/cm3 and cm3 volumes) ----
  # dCs/dt (ng/g/day)  (note: Vpw_cm3 used to convert to same basis)
  dCsdt <- - ksed * Vpw_cm3 / ms_g * (Cs_pw_eq_cm3 - Cpw_cm3)
  
  # dCpw/dt (ng/cm3/day)  (exchange with sediment porewater and water)
  dCpwdt <-  ksed * (Cs_pw_eq_cm3 - Cpw_cm3) -    # simplified as discussed (Vpw_cm3 cancels)
    kpw * Aws / Vpw_cm3 * (Cpw_cm3 - Cw_cm3) -
    kb * Cpw_cm3
  
  # dCw/dt (ng/cm3/day)
  dCwdt <- kpw * Aws / Vw_cm3 * (Cpw_cm3 - Cw_cm3) -
    kaw * Aaw / Vw_cm3 * (Cw_cm3 - Ca_cm3 / Kaw.t)
  
  # dCa/dt (ng/cm3/day)
  dCadt <- kaw * Aaw / Va_cm3 * (Cw_cm3 - Ca_cm3 / Kaw.t) -
    ro * Apuf / Va_cm3 * (Ca_cm3 - Cpuf_cm3 / Kpuf)
  
  # dCpuf/dt (ng/cm3/day)
  dCpufdt <- ro * Apuf / Vpuf_cm3 * (Ca_cm3 - Cpuf_cm3 / Kpuf)
  
  # convert aqueous derivatives back to ng/L/day (multiply by 1000)
  out <- c(
    as.numeric(dCsdt),
    as.numeric(dCpwdt * 1000),
    as.numeric(dCwdt  * 1000),
    as.numeric(dCadt  * 1000),
    as.numeric(dCpufdt * 1000)
  )
  
  # final guard: ensure length is 5 and numeric
  if (length(out) != 5) stop("rtm_simple: derivative vector length != 5")
  if (!is.numeric(out)) stop("rtm_simple: derivatives not numeric")
  
  return(list(out))
}

# -------------------------
# 3) Congener-specific chemistry & transport precompute (PCB_32)
# -------------------------
pc <- read.csv("Data/05_PCP.csv", stringsAsFactors = FALSE)
pc_row <- pc[pc$congener == pcb.ind, ]

MW.pcb <- pc_row$MW
Kow <- pc_row$Kow
dUow <- pc_row$dUow
Kaw <- pc_row$Kaw
dUaw <- pc_row$dUaw
Koa <- pc_row$Koa
E <- pc_row$E; S <- pc_row$S; A <- pc_row$A; B <- pc_row$B; V <- pc_row$V

# physical constants & temps
MH2O <- 18.0152; MCO2 <- 44.0094; R <- 8.3144
Tst.1 <- 273.15 + 25; Tw.1 <- 273.15 + 20

# temp-corrected partition coefficients
Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 - 1 / Tst.1))
Kaw.t <- Kaw * exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1

# PUF & SPME parameters
Vpuf_cm3 <- Vpuf_L * 1000
d_puf <- 0.0213 * 100^3
Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) * d_puf

Af <- 0.138
Vf_cm3_per_cm <- 0.000000069 * 1000
L <- 1
Vf_tot <- Vf_cm3_per_cm * L
Kf <- 10^(1.06 * log10(Kow.t) - 1.16)

# Diffusion / transfer coefficients
D.water.air <- 0.2743615
D.co2.w <- 1.67606E-05
D.pcb.air <- D.water.air * (MW.pcb / MH2O)^(-0.5)
D.pcb.water <- D.co2.w * (MW.pcb / MCO2)^(-0.5)
v.H2O <- 0.010072884
V.water.air <- 0.003
V.co2.w <- 4.1e-2
SC.pcb.w <- v.H2O / D.pcb.water

# kpw (porewater-water) [cm/day]
bl <- 0.1 # it can be modified
kpw <- D.pcb.water * 60 * 60 * 24 / bl

# kaw (air-water) [cm/day]
Kaw.a <- V.water.air * (D.pcb.air / D.water.air)^(0.67)
Kaw.w <- V.co2.w * (SC.pcb.w / 600)^(-0.5)
kaw.o <- (1 / (Kaw.a * Kaw.t) + (1 / Kaw.w))^-1
kaw <- kaw.o * 100 * 60 * 60 * 24

# ksed initial (Koelmans regression)
logksed <- -0.832 * log10(Kow.t) + 1.34
ksed_init <- 10^(logksed)

# Kd via Abraham Koc
logKoc <- 1.1 * E - 0.72 * S + 0.15 * A - 1.98 * B + 2.28 * V + 0.14
Koc <- 10^(logKoc)
foc <- 0.03
Kd <- Koc * foc   # L/kg

# -------------------------
# 4) Fixed parameters
# -------------------------
# compute Vs in cm3 associated with ms (same formula you used earlier)
n  <- 0.42
ds <- 1540      # g / L
ms_g <- 10      # g sediment in vial
M  <- ds * (1 - n) / n
Vpw_cm3 <- ms_g / M * 1000
Vpw_L <- Vpw_cm3 / 1000

parms_base <- list(
  ms   = ms_g,           # g sediment in vial
  Kd   = Kd,             # L/kg
  ksed = ksed_init,      # initial guess (1/day)
  kpw  = kpw,            # cm/day
  kaw  = kaw,            # cm/day
  kb   = 0,              # 1/day (biodeg)
  ro   = 400,            # 1/day (PUF-side exchange)
  Kpuf = Kpuf,
  Kaw.t= Kaw.t,
  Af   = Af,
  L    = L,
  Vf_tot = Vf_tot,
  Apuf = Apuf,
  Aaw  = Aaw,
  Aws  = Aws,
  Vpw_cm3 = Vpw_cm3,
  Vw_L = Vw_L,
  Vpw_L = Vpw_L,
  Va_L = Va_L,
  Vpuf_L = Vpuf_L
)

# -------------------------
# 5) Initial conditions & times (mass-consistent)
# -------------------------
bulk_conc <- read.csv("Data/03_NBH_SedimentPCB.csv",
                      stringsAsFactors = FALSE,
                      check.names = FALSE)
Ct <- mean(bulk_conc[[pcb.ind]])   # ng/g

# Total mass in sediment layer (ng)
M_sed_init <- Ct * parms_base$ms

# Solve Cs_init so that M_sed_init = ms*Cs + dissolved mass in porewater
# dissolved mass = Cpw (ng/L) * Vpw_L (L), and Cpw = Cs*1000 / Kd (ng/L)
# -> Cs_init = M_sed_init / (ms + (1000 * Vpw_L)/Kd)
Cs_init <- M_sed_init / (parms_base$ms + (1000 * parms_base$Vpw_L) / parms_base$Kd)
Cpw_init <- Cs_init * 1000 / parms_base$Kd

cinit <- c(Cs = Cs_init, Cpw = Cpw_init, Cw = 0, Ca = 0, Cpuf = 0)
times <- sort(unique(c(0, obs_df$time)))

# -------------------------
# 6) Helper: run model and compute mpuf (ng)
# -------------------------
run_model <- function(ksed_try, parms_template = parms_base, times_run = times) {
  parms <- parms_template
  parms$ksed <- ksed_try
  out <- ode(y = cinit, times = times_run, func = rtm_simple, parms = parms)
  df <- as.data.frame(out)
  colnames(df) <- c("time","Cs","Cpw","Cw","Ca","Cpuf")
  df$mpuf <- df$Cpuf * parms$Vpuf_L  # ng
  return(df)
}

# -------------------------
# 7) Objective: SSE between modeled mpuf and observed mpuf_obs at obs times
# -------------------------
obj_fn <- function(log_ksed) {
  ksed_try <- exp(log_ksed)
  df <- run_model(ksed_try)
  cmp <- df %>% select(time, mpuf) %>% right_join(obs_df, by = "time")
  if(nrow(cmp) == 0) return(1e30)
  sse <- sum((cmp$mpuf - cmp$mpuf_obs)^2, na.rm = TRUE)
  return(sse)
}

# -------------------------
# 8) Optimize (log-space)
# -------------------------
init_guess <- parms_base$ksed
res <- optim(par = log(init_guess), fn = obj_fn,
             method = "L-BFGS-B", lower = log(1e-12), upper = log(10),
             control = list(maxit = 500))
ksed_best <- exp(res$par)
n_obs <- nrow(obs_df)
RMSE <- sqrt(res$value / max(1, n_obs))
cat("Best-fit ksed (1/day):", signif(ksed_best,4), "\nRMSE:", signif(RMSE,4), "\n\n")

# -------------------------
# 9) Final run & comparison
# -------------------------
df_best <- run_model(ksed_best)
comparison <- df_best %>% select(time, mpuf_model = mpuf) %>% right_join(obs_df, by = "time")
cat("Model vs Observed (mpuf in ng):\n")
print(comparison)

# -------------------------
# 10) Mass-balance diagnostics (ng)
# -------------------------
df_diag <- df_best %>%
  mutate(ms_mass = Cs * parms_base$ms,          # sorbed mass on solids (ng)
         mpw = Cpw * parms_base$Vpw_L,          # porewater mass in sediment (ng)
         mw  = Cw  * parms_base$Vw_L,
         ma  = Ca  * parms_base$Va_L,
         mpuf = mpuf) %>%
  mutate(mt = ms_mass + mpw + mw + ma + mpuf) %>%
  filter(time %in% obs_df$time) %>%
  select(time, Cs, ms_mass, mpw, mw, ma, mpuf, mt) %>%
  arrange(time)

cat("\nMass-balance diagnostics at observation times (ng):\n")
print(df_diag)

# --- Diagnostics: Deff vs molecular diffusion for various radii ----
# Requires that you already have ksed_best (1/day) and parms_base$ksed (Koelmans) in workspace.
# If those names differ, change them below.

# inputs (use your values or adjust)
ksed_fit_day <- ksed_best           # fitted ksed [1/day]
ksed_koel_day <- parms_base$ksed    # Koelmans ksed [1/day], fallback

# compute (approx) molecular diffusion in water (m^2/s)
D_co2_w <- 1.67606e-05              # cm^2/s (your constant)
D_pcb_water_cm2_s <- D_co2_w * (MW.pcb / 44.0094)^(-0.5)
D_mol_m2_s <- D_pcb_water_cm2_s * 1e-4

# radii to test (in meters) -- change/add values you like (mm -> m)
r_mm <- c(0.5, 1, 2, 3, 5, 10)   # mm
r_m  <- r_mm / 1000

# helper conversions
sec_per_day <- 86400

# build table
out <- data.frame(
  r_mm = r_mm,
  r_m = r_m
)

out$Deff_fit_m2_s    <- (ksed_fit_day / sec_per_day) * (out$r_m^2)
out$Deff_fit_m2_day  <- ksed_fit_day * (out$r_m^2)
out$ratio_fit_to_Dmol<- out$Deff_fit_m2_s / D_mol_m2_s

out$Deff_koel_m2_s    <- (ksed_koel_day / sec_per_day) * (out$r_m^2)
out$Deff_koel_m2_day  <- ksed_koel_day * (out$r_m^2)
out$ratio_koel_to_Dmol<- out$Deff_koel_m2_s / D_mol_m2_s

# ksed_max that keeps Deff <= D_mol for that r: ksed_max = D_mol / r^2 (in 1/s), converted to 1/day
out$ksed_max_1_day <- (D_mol_m2_s / (out$r_m^2)) * sec_per_day

# round/format for display
out_print <- out %>%
  mutate(
    Deff_fit_m2_s = signif(Deff_fit_m2_s, 4),
    Deff_fit_m2_day = signif(Deff_fit_m2_day, 4),
    ratio_fit_to_Dmol = signif(ratio_fit_to_Dmol, 4),
    Deff_koel_m2_s = signif(Deff_koel_m2_s, 4),
    Deff_koel_m2_day = signif(Deff_koel_m2_day, 4),
    ratio_koel_to_Dmol = signif(ratio_koel_to_Dmol, 4),
    ksed_max_1_day = signif(ksed_max_1_day, 4)
  )

cat("Molecular diffusion (approx) D_mol =", format(D_mol_m2_s, scientific=TRUE, digits=3), "m^2/s\n\n")
print(out_print)

