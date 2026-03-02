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
Apuf <- 7.07  # cm2 (PUF area)

# -------------------------
# 1) Read & prepare observed data (Non-shaken)
exp.ns.data <- read.csv("Data/01_NS_SPME_PUF.csv", stringsAsFactors = FALSE, check.names = FALSE)

pcb.ind <- "PCB_32"

pcbi <- exp.ns.data[, c("ID", "Sample_medium", "Group", "time", "Replicate", pcb.ind)]

# PUF control samples (mpuf) for NBH_NS
pcbi.puf.control <- pcbi %>%
  filter(ID == "NBH_NS", Sample_medium == "PUF", Group == "Control") %>%
  rename(mpuf_control = !!sym(pcb.ind)) %>%
  select(time, mpuf_control)

# Add t = 0
pcbi_control_0 <- bind_rows(data.frame(time = 0, mpuf_control = 0), pcbi.puf.control) %>%
  arrange(time)

# Observed PUF by time (mean across replicates)
obs_df <- as_tibble(pcbi_control_0) %>%
  group_by(time) %>%
  summarise(mpuf_obs = mean(mpuf_control, na.rm = TRUE), .groups = "drop") %>%
  arrange(time) %>%
  filter(!is.na(mpuf_obs))

# times vector for static experiment
times <- sort(unique(c(0, obs_df$time)))

# -------------------------
# 2) Congener chemistry & physical constants (load PCP table)
# -------------------------
pc <- read.csv("Data/05_PCP.csv", stringsAsFactors = FALSE)
pc_row <- pc[pc$congener == pcb.ind, ]

if(nrow(pc_row) == 0) stop("Congener not found in Data/05_PCP.csv: ", pcb.ind)

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

# -------------------------
# 3) PUF & SPME parameters (derived)
# -------------------------
Vpuf_cm3 <- Vpuf_L * 1000
d_puf <- 0.0213 * 100^3
Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) * d_puf

# SPME fiber params (Af, Vf_tot, Kf)
Af <- 0.138
Vf_cm3_per_cm <- 0.000000069 * 1000
L <- 1
Vf_tot <- Vf_cm3_per_cm * L
Kf <- 10^(1.06 * log10(Kow.t) - 1.16)

# -------------------------
# 4) Diffusion / transfer coefficients (precompute)
# -------------------------
D.water.air <- 0.2743615
D.co2.w <- 1.67606E-05
D.pcb.air <- D.water.air * (MW.pcb / MH2O)^(-0.5)
D.pcb.water <- D.co2.w * (MW.pcb / MCO2)^(-0.5)
v.H2O <- 0.010072884
V.water.air <- 0.003
V.co2.w <- 4.1e-2
SC.pcb.w <- v.H2O / D.pcb.water

# kpw (porewater-water) [cm/day] using boundary layer bl (cm)
bl <- 0.1
kpw <- D.pcb.water * 60 * 60 * 24 / bl

# kaw (air-water) [cm/day]
Kaw.a <- V.water.air * (D.pcb.air / D.water.air)^(0.67)
Kaw.w <- V.co2.w * (SC.pcb.w / 600)^(-0.5)
kaw.o <- (1 / (Kaw.a * Kaw.t) + (1 / Kaw.w))^-1
kaw <- kaw.o * 100 * 60 * 60 * 24

# -------------------------
# 5) Partitioning: Kd via Abraham Koc
# -------------------------
logKoc <- 1.1 * E - 0.72 * S + 0.15 * A - 1.98 * B + 2.28 * V + 0.14
Koc <- 10^(logKoc)
foc <- 0.03
Kd <- Koc * foc   # L/kg

# -------------------------
# 6) Sediment geometry -> compute porewater volume (Vpw) once and for all
# -------------------------
n  <- 0.42
ds <- 1540      # g / L (particle density)
ms_g <- 10      # g sediment in vial

M  <- ds * (1 - n) / n
Vpw_cm3 <- ms_g / M * 1000   # porewater volume connected to ms (cm3)
Vpw_L   <- Vpw_cm3 / 1000

# -------------------------
# 7) Single-domain desorption (Werner truncated series)
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
# 8) Build parameter list (canonical fields; Vpw included)
# -------------------------
# compute initial total mass in sediment (ng) from measured bulk conc (use data file you have)
bulk_conc <- read.csv("Data/03_NBH_SedimentPCB.csv", stringsAsFactors = FALSE, check.names = FALSE)
if (!(pcb.ind %in% colnames(bulk_conc))) stop("Congener not found in 03_NBH_SedimentPCB.csv: ", pcb.ind)
Ct <- mean(bulk_conc[[pcb.ind]], na.rm = TRUE)   # ng/g

parms_base <- list(
  ms   = ms_g,
  Kd   = Kd,
  kpw  = kpw,
  kaw  = kaw,
  kb   = 0,
  ro   = 400,            # initial placeholder; will be fitted unless you set otherwise
  Kpuf = Kpuf,
  Kaw.t= Kaw.t,
  Af   = Af,
  L    = L,
  Vf_tot = Vf_tot,
  Apuf = Apuf,
  Aaw  = Aaw,
  Aws  = Aws,
  Vpw_cm3 = Vpw_cm3,
  Vpw_L   = Vpw_L,
  Vw_L = Vw_L,
  Va_L = Va_L,
  Vpuf_L = Vpuf_L,
  # desorption defaults
  D_r2 = 3.5e-6,         # initial guess for D_a / r^2 (day^-1)
  M_sed_init = NA,       # computed below
  Nterms_des = 200
)

# -------------------------
# 9) Initial conditions (mass-consistent Cs_init recommended)
# -------------------------
# Total mass in sediment layer (ng) if you want to use Ct (bulk measured)
M_sed_init <- Ct * parms_base$ms   # ng (ms in g; Ct in ng/g)
# Solve Cs_init so that M_sed_init = ms*Cs + dissolved mass in porewater
# dissolved mass = Cpw (ng/L) * Vpw_L (L), and Cpw = Cs*1000 / Kd (ng/L)
# -> Cs_init = M_sed_init / (ms + (1000 * Vpw_L)/Kd)
Cs_init <- M_sed_init / (parms_base$ms + (1000 * parms_base$Vpw_L) / parms_base$Kd)
Cpw_init <- Cs_init * 1000 / parms_base$Kd

parms_base$M_sed_init <- M_sed_init

cinit <- c(Cs = Cs_init, Cpw = Cpw_init, Cw = 0, Ca = 0, Cpuf = 0)

# -------------------------
# 10) ODE: rtm_simple using desorption_single as source
# -------------------------
rtm_simple <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # volumes (canonical)
    Vpw_cm3  <- parms$Vpw_cm3
    Vw_cm3   <- parms$Vw_L  * 1000
    Va_cm3   <- parms$Va_L  * 1000
    Vpuf_cm3 <- parms$Vpuf_L * 1000
    
    # states (ensure numeric)
    Cs   <- as.numeric(state[1])   # ng/g
    Cpw  <- as.numeric(state[2])   # ng/L
    Cw   <- as.numeric(state[3])   # ng/L
    Ca   <- as.numeric(state[4])   # ng/L
    Cpuf <- as.numeric(state[5])   # ng/L
    
    # convert to ng/cm3 for flux calcs
    Cpw_cm3  <- Cpw  / 1000
    Cw_cm3   <- Cw   / 1000
    Ca_cm3   <- Ca   / 1000
    Cpuf_cm3 <- Cpuf / 1000
    
    # convert Cs (ng/g) to porewater-equivalent (ng/cm3)
    Cs_pw_eq_cm3 <- Cs / parms$Kd
    
    # params
    kpw  <- parms$kpw
    kaw  <- parms$kaw
    kb   <- parms$kb
    ro   <- parms$ro
    Kpuf <- parms$Kpuf
    Kaw.t<- parms$Kaw.t
    Apuf <- parms$Apuf
    Aws  <- parms$Aws
    Aaw  <- parms$Aaw
    ms_g <- parms$ms
    
    # desorption source: dm_dt in ng/day
    des <- desorption_single(t_days = t, D_r2 = parms$D_r2, mtot_ng = parms$M_sed_init, Nterms = parms$Nterms_des)
    dm_dt_ng_per_day <- des$dm_dt[1]
    
    # ODEs inside ng/cm3/day (except Cs ng/g/day)
    dCsdt <- - dm_dt_ng_per_day / ms_g
    
    dCpwdt <- (dm_dt_ng_per_day / Vpw_cm3) -
      kpw * Aws / Vpw_cm3 * (Cpw_cm3 - Cw_cm3) -
      kb * Cpw_cm3
    
    dCwdt <- kpw * Aws / Vw_cm3 * (Cpw_cm3 - Cw_cm3) -
      kaw * Aaw / Vw_cm3 * (Cw_cm3 - Ca_cm3 / Kaw.t)
    
    dCadt <- kaw * Aaw / Va_cm3 * (Cw_cm3 - Ca_cm3 / Kaw.t) -
      ro * Apuf / Va_cm3 * (Ca_cm3 - Cpuf_cm3 / Kpuf)
    
    dCpufdt <- ro * Apuf / Vpuf_cm3 * (Ca_cm3 - Cpuf_cm3 / Kpuf)
    
    # convert back to state units: aqueous/gas ng/cm3/day -> ng/L/day multiply 1000
    return(list(c(
      as.numeric(dCsdt),
      as.numeric(dCpwdt * 1000),
      as.numeric(dCwdt  * 1000),
      as.numeric(dCadt  * 1000),
      as.numeric(dCpufdt * 1000)
    )))
  })
}

# -------------------------
# 11) Helper run_model(pars)
# -------------------------
run_model <- function(pars, times_run = times) {
  out <- ode(y = cinit, times = times_run, func = rtm_simple, parms = pars, atol = 1e-8, rtol = 1e-8)
  df <- as.data.frame(out)
  colnames(df) <- c("time","Cs","Cpw","Cw","Ca","Cpuf")
  df$mpuf <- df$Cpuf * pars$Vpuf_L  # ng
  return(df)
}

# -------------------------
# 12) Objective: fit D_r2 and ro simultaneously (log-space)
# -------------------------
obj_fit_Dr2_ro <- function(par_log) {
  D_r2_try <- exp(par_log[1])
  ro_try   <- exp(par_log[2])
  
  p <- parms_base
  p$D_r2 <- D_r2_try
  p$ro   <- ro_try
  
  out <- tryCatch({
    run_model(p, times_run = times)
  }, error = function(e) {
    message("ODE error: ", e$message)
    return(NULL)
  })
  
  if (is.null(out)) return(1e30)
  
  cmp <- out %>% select(time, mpuf) %>% right_join(obs_df, by = "time")
  if (nrow(cmp) == 0) return(1e30)
  
  # log-space SSE (robust to zeros)
  sse <- sum((log1p(cmp$mpuf) - log1p(cmp$mpuf_obs))^2, na.rm = TRUE)
  return(sse)
}

# initial guesses (log-space)
init2 <- c(log(3.5e-6), log(300))   # log(D_r2_init), log(ro_init)

# bounds: D_r2 between 1e-12..1e-1 day^-1, ro between 1..1e4 cm/day
lower2 <- c(log(1e-12), log(1))
upper2 <- c(log(1e-1), log(1e4))

# -------------------------
# 13) Fit (this may take some time)
# -------------------------
res2 <- optim(par = init2, fn = obj_fit_Dr2_ro, method = "L-BFGS-B",
              lower = lower2, upper = upper2,
              control = list(maxit = 2000))

pars_fit <- exp(res2$par)
D_r2_hat <- pars_fit[1]
ro_hat <- pars_fit[2]

cat("Fitted D/r^2 (day^-1):", signif(D_r2_hat,6), "\n")
cat("Fitted ro (cm/day):", signif(ro_hat,6), "\n")
cat("SSE (log-space):", signif(res2$value,6), "\n\n")

# -------------------------
# 14) Final run & comparisons
# -------------------------
p_final <- parms_base
p_final$D_r2 <- D_r2_hat
p_final$ro   <- ro_hat

df_best <- run_model(p_final, times_run = times)
comparison <- df_best %>% select(time, mpuf_model = mpuf) %>% right_join(obs_df, by = "time")
cat("Model vs Observed (mpuf in ng):\n")
print(comparison)

# -------------------------
# 15) Mass-balance diagnostics (ng)
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
