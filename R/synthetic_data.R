###############################################################################
# synthetic_data.R
# Hourly synthetic DO/T records with KNOWN metabolism, internal waves,
# non-coherent physical noise and (surface case) air-lake exchange and
# diel heating. Two cases:
#   "metalimnion": sensor below the mixed layer, moderate IW, no gas exchange
#   "metalimnion_strongIW": same, IW-induced DO variance >> diel biological signal
#   "surface"    : sensor in the mixed layer, weak IW, gas exchange,
#                  diel heating (T' correlated with the metabolic cycle)
###############################################################################

make_synthetic_lake <- function(case = c("metalimnion", "metalimnion_strongIW", "surface"), days = 80,
                                start = "2019-06-01", lat = 46.5, lon = 6.66,
                                elev = 372, kd = 0.3, seed = 42) {
  case <- match.arg(case)
  set.seed(seed)
  time <- seq(as.POSIXct(start, tz = "UTC"), by = 3600, length.out = days * 24)
  nt <- length(time); th <- (seq_len(nt) - 1)           # hours since start
  dnum <- th / 24

  ar1 <- function(n, tau_h, sd) {                       # red noise, unit-variance scaled
    phi <- exp(-1 / tau_h); e <- rnorm(n); x <- numeric(n); x[1] <- e[1]
    for (i in 2:n) x[i] <- phi * x[i - 1] + sqrt(1 - phi^2) * e[i]
    sd * x
  }

  # --- light: clear-sky sine between sunrise and sunset (solar time) x clouds
  DL   <- day_length_h(time, lat)
  hsol <- (as.numeric(format(time, "%H")) + lon / 15) %% 24
  sunrise <- 12 - DL / 2
  cloud <- rep(runif(days, 0.35, 1), each = 24)
  I0 <- pmax(0, 850 * cloud * sin(pi * (hsol - sunrise) / DL))
  I0[hsol < sunrise | hsol > sunrise + DL] <- 0

  p <- switch(case,
    metalimnion = list(z = 10, zmix = 5, T0 = 14, dTdz = -1.0, dDOdz = -5,
                       zeta_sd = 0.6, seiche = c(8, 0.5, 22, 0.3), heat = 0,
                       Pmax = function(d) 1.0 + 0.5 * sin(2 * pi * d / 40),
                       alpha = 0.02, R20 = 0.25, noise_sd = 3, DO0 = 300),
    metalimnion_strongIW = list(z = 10, zmix = 5, T0 = 14, dTdz = -1.0, dDOdz = -8,
                       zeta_sd = 1.5, seiche = c(8, 1.0, 22, 0.6), heat = 0,
                       Pmax = function(d) 1.0 + 0.5 * sin(2 * pi * d / 40),
                       alpha = 0.02, R20 = 0.25, noise_sd = 3, DO0 = 300),
    surface     = list(z = 1, zmix = 6, T0 = 20, dTdz = -0.1, dDOdz = -1,
                       zeta_sd = 0.3, seiche = c(8, 0.2, 22, 0.1), heat = 0.8,
                       Pmax = function(d) 2.5 + 1.0 * sin(2 * pi * d / 40),
                       alpha = 0.01, R20 = 1.0, noise_sd = 2, DO0 = 290))

  # --- vertical displacement of isotherms/iso-oxygen surfaces (m)
  env <- 0.6 + 0.4 * sin(2 * pi * dnum / 9)^2           # "wind events" modulate seiches
  zeta <- ar1(nt, 12, p$zeta_sd) +
          env * p$seiche[2] * sin(2 * pi * th / p$seiche[1]) +
          env * p$seiche[4] * sin(2 * pi * th / p$seiche[3] + 1)

  # --- temperature: background + diel heating (surface) + IW displacement
  T_bg <- p$T0 + 0.03 * dnum + p$heat * sin(2 * pi * (hsol - 10) / 24)
  T_obs <- T_bg + p$dTdz * zeta + rnorm(nt, 0, 0.01)

  # --- true metabolism (hourly, mmol m-3 h-1)
  Iz  <- I0 * exp(-kd * p$z)
  Pm  <- p$Pmax(dnum)
  gpp <- Pm * tanh(p$alpha * Iz / Pm)
  r   <- p$R20 * 1.07^(T_bg - 20)                        # respiration at in-situ T

  # --- biological DO (and gas exchange in the surface case), Euler, 6-min steps
  wind <- pmax(0.5, 3 + 1.5 * sin(2 * pi * (hsol - 14) / 24) + ar1(nt, 24, 1.5))
  kg   <- if (case == "surface") k_gas_cole(wind, T_bg) else rep(0, nt)
  DO_bio <- numeric(nt); DO_bio[1] <- p$DO0; sub <- 10
  for (i in 2:nt) {
    x <- DO_bio[i - 1]
    for (s in seq_len(sub)) {
      Fg <- kg[i - 1] * (x - o2_sat(T_bg[i - 1], elev))
      x  <- x + (gpp[i - 1] - r[i - 1] - Fg / p$zmix) / sub
    }
    DO_bio[i] <- x
  }

  # --- observed DO: biology + IW dislocation (non-linear) + non-coherent noise
  DO_obs <- DO_bio + p$dDOdz * zeta + 0.3 * p$dDOdz / 8 * zeta^2 +
            ar1(nt, 6, p$noise_sd) + rnorm(nt, 0, 0.5)

  data.table(time = time, DO = DO_obs, T = T_obs, I0 = I0, wind = wind,
             zmix = p$zmix, gpp_true = gpp, r_true = r, DO_bio = DO_bio,
             zeta = zeta, z_sensor = p$z)
}

truth_daily <- function(syn) {
  syn[, .(GPP_true = sum(gpp_true), R_true = sum(r_true),
          NEP_true = sum(gpp_true - r_true)), by = .(day = as.Date(time, tz = "UTC"))]
}
