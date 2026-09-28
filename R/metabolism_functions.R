###############################################################################
# metabolism_functions.R
# Shared helpers for the diel oxygen methods of Fernández Castro et al. (2021)
#   - O2 saturation, Schmidt number, piston velocity, air-lake flux
#   - day length (for t_night)
#   - Gaussian smoother (the "<.>" low-passed signals of the original code)
#   - smoothed DO budget -> NEP (used by the spectral method)
#   - time-domain diel method: DO_clean -> hourly nep -> Eq. S7 fit -> GPP, R, NEP
#
# Units throughout: DO in mmol m-3, T in degC, time step in hours,
# hourly rates in mmol m-3 h-1, daily rates in mmol m-3 d-1,
# k_gas in m h-1, F_gas in mmol m-2 h-1 (POSITIVE = O2 leaving the lake).
###############################################################################

suppressPackageStartupMessages(library(data.table))

## ---------------------------------------------------------------- checks ----
stop_if_missing <- function(dt, cols, where = "") {
  miss <- setdiff(cols, names(dt))
  if (length(miss))
    stop(where, ": missing column(s) ", paste(miss, collapse = ", "),
         "\n  Available: ", paste(names(dt), collapse = ", "), call. = FALSE)
}

# Put a series on a regular grid (default hourly), averaging within bins.
# The spectral method assumes 24 samples per day (Welch nperseg = 128 h, etc.).
regularize <- function(dt, step_min = 60, tz = "UTC") {
  stop_if_missing(dt, "time", "regularize")
  dt <- as.data.table(dt)
  num <- names(dt)[vapply(dt, is.numeric, logical(1))]
  step_s <- step_min * 60
  dt[, tbin := as.POSIXct(round(as.numeric(time) / step_s) * step_s,
                          origin = "1970-01-01", tz = tz)]
  out <- dt[, lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)),
            by = tbin, .SDcols = num]
  grid <- data.table(tbin = seq(min(out$tbin), max(out$tbin), by = step_s))
  out <- out[grid, on = "tbin"]
  setnames(out, "tbin", "time")
  out[]
}

## ------------------------------------------------------ O2 & gas exchange ----
# O2 saturation (Garcia & Gordon 1992, Benson & Krause fit), S = 0,
# with a barometric correction for elevation. Returns mmol m-3.
o2_sat <- function(T, elev = 0) {
  Ts <- log((298.15 - T) / (273.15 + T))
  lnC <- 5.80871 + 3.20291 * Ts + 4.17887 * Ts^2 + 5.10006 * Ts^3 -
    9.86643e-2 * Ts^4 + 3.80369 * Ts^5                      # umol kg-1 at 1 atm
  rho <- 1000 * (1 - (T + 288.9414) / (508929.2 * (T + 68.12963)) * (T - 3.9863)^2)
  P   <- (1 - 2.25577e-5 * elev)^5.25588                         # atm
  exp(lnC) * rho / 1000 * P
}

# Schmidt number of O2 in freshwater (Wanninkhof 1992)
schmidt_o2 <- function(T) 1800.6 - 120.10 * T + 3.7818 * T^2 - 0.047608 * T^3

# k600 from wind (Cole & Caraco 1998), cm h-1  ->  k_O2 in m h-1
k_gas_cole <- function(u10, T, n = -0.5) {
  k600 <- 2.07 + 0.215 * u10^1.7
  k600 * (schmidt_o2(T) / 600)^n / 100
}

# Air-lake flux, mmol m-2 h-1, positive = outgassing
gas_flux <- function(DO, T, k_gas, elev = 0) k_gas * (DO - o2_sat(T, elev))

## ------------------------------------------------------------------ sun ----
# Astronomical day length (h), with -0.833 deg refraction/disc correction
day_length_h <- function(time, lat) {
  doy  <- as.numeric(format(time, "%j"))
  decl <- 23.44 * pi / 180 * sin(2 * pi * (284 + doy) / 365)
  phi  <- lat * pi / 180
  cosw <- (sin(-0.833 * pi / 180) - sin(phi) * sin(decl)) / (cos(phi) * cos(decl))
  2 * acos(pmin(pmax(cosw, -1), 1)) * 180 / pi / 15
}

## ------------------------------------------------------------ smoothing ----
# Gaussian smoother used in DO_budget.py:
#   w = exp(-(dt)^2 / (0.5*hw)^2) within +/- hw  (hw = 30 d in the original)
gauss_smooth <- function(x, time, hw_days = 30) {
  td <- as.numeric(time) / 86400
  out <- rep(NA_real_, length(x))
  ok  <- is.finite(x)
  for (i in seq_along(x)) {
    j <- which(ok & abs(td - td[i]) <= hw_days)
    if (!length(j)) next
    w <- exp(-((td[j] - td[i])^2) / (0.5 * hw_days)^2)
    out[i] <- sum(w * x[j]) / sum(w)
  }
  out
}

# First centred differences (edges one-sided), dy/dx on finite points only
centered_diff <- function(x, y) {
  out <- rep(NA_real_, length(y))
  ok <- which(is.finite(x) & is.finite(y))
  if (length(ok) < 2) return(out)
  x0 <- x[ok]; y0 <- y[ok]; n <- length(ok)
  d <- numeric(n)
  d[1] <- (y0[2] - y0[1]) / (x0[2] - x0[1])
  d[n] <- (y0[n] - y0[n - 1]) / (x0[n] - x0[n - 1])
  if (n > 2) d[2:(n - 1)] <- (y0[3:n] - y0[1:(n - 2)]) / (x0[3:n] - x0[1:(n - 2)])
  out[ok] <- d
  out
}

## ------------------------------------------------ physical source terms ----
# Adds F_gas (mmol m-2 h-1) and gas_term = F_gas/zmix for sensors in the mixed
# layer (mmol m-3 h-1), plus a transport term (default 0).
#   dDO/dt = nep - gas_term + transport   ->   nep = dDO/dt + gas_term - transport
add_physics_terms <- function(dt, cfg) {
  dt <- copy(dt)
  if (!"transport" %in% names(dt)) dt[, transport := 0]
  if (!"zmix" %in% names(dt)) dt[, zmix := cfg$zmix_default]
  if (!"k_gas" %in% names(dt)) {
    if ("wind" %in% names(dt)) dt[, k_gas := k_gas_cole(wind, T)]
    else dt[, k_gas := NA_real_]
  }
  dt[, in_ml := cfg$z_sensor < zmix]
  if (any(dt$in_ml, na.rm = TRUE) && all(is.na(dt$k_gas)))
    stop("Sensor is in the mixed layer but neither 'k_gas' nor 'wind' is provided.", call. = FALSE)
  dt[, F_gas := fifelse(in_ml, gas_flux(DO, T, k_gas, cfg$elev), 0)]
  dt[, gas_term := fifelse(in_ml & is.finite(zmix) & zmix > 0, F_gas / zmix, 0)]
  dt[]
}

## ------------------------------------------- NEP from the smoothed budget ----
# Long-term (fortnight-scale) NEP used by the spectral method, mmol m-3 d-1:
#   NEP = d<DO>/dt + <F_gas/zmix> - <transport>, where <.> = gauss_smooth(30 d)
nep_budget_smooth <- function(dt, hw_days = 30) {
  sDO  <- gauss_smooth(dt$DO, dt$time, hw_days)
  sgas <- gauss_smooth(dt$gas_term, dt$time, hw_days) * 24
  str  <- gauss_smooth(dt$transport, dt$time, hw_days) * 24
  dsDO <- centered_diff(as.numeric(dt$time) / 86400, sDO)
  list(sDO = sDO, NEP = dsDO + sgas - str)
}

## ------------------------------------ TIME-DOMAIN diel method (Text S2) ----
# Eq. S7: nep_i = Pmax * tanh(alpha * Iz_i / Pmax) - R20 * theta^(T_i - 20)
fit_eq7 <- function(nep, Iz, T, theta = 1.07) {
  ok <- is.finite(nep) & is.finite(Iz) & is.finite(T)
  if (sum(ok) < 8) return(NULL)
  nep <- nep[ok]; Iz <- Iz[ok]; T <- T[ok]
  pred <- function(p) p[1] * tanh(p[2] * Iz / p[1]) - p[3] * theta^(T - 20)
  sse  <- function(p) sum((nep - pred(p))^2)
  night <- Iz < 1
  R0 <- if (any(night)) max(-mean(nep[night]), 1e-3) else max(-min(nep), 1e-3)
  P0 <- max(max(nep) + R0, 1e-2)
  a0 <- P0 / max(stats::median(Iz[Iz > 1]), 1)
  best <- NULL
  for (s in list(c(P0, a0, R0), c(2 * P0, a0 / 2, R0), c(P0 / 2, 2 * a0, R0))) {
    f <- try(optim(s, sse, method = "L-BFGS-B",
                   lower = c(1e-6, -Inf, -Inf), upper = c(100 * P0, Inf, Inf)),
             silent = TRUE)
    if (!inherits(f, "try-error") && (is.null(best) || f$value < best$value)) best <- f
  }
  # L-BFGS-B code 52 (line-search warning) usually sits at the optimum: accept
  if (is.null(best) || !is.finite(best$value) || best$convergence == 1) return(NULL)
  p <- best$par
  list(Pmax = p[1], alpha = p[2], R20 = p[3],
       r2 = 1 - best$value / sum((nep - mean(nep))^2),
       pred = pred, par = p)
}

# dt: regular hourly data with time, DO, T, I0 (+ zmix, k_gas|wind, transport)
# td_env: environment where iw_correction_time_domain.R was sourced
td_metabolism <- function(dt, cfg, td_env) {
  td_cfg <- modifyList(td_env$cfg, list(make_plots = FALSE, tz = cfg$tz,
                                        depth_col = NULL, input_file = NULL,   # single sensor here
                                        day_start_hour = cfg$day_start_hour,
                                        min_coverage = cfg$min_coverage))
  iw <- td_env$correct_internal_waves(
    td_env$prepare_input(dt[, .(time, DO, T)], td_cfg), td_cfg)
  x <- merge(dt, iw$data[, .(time, day, DO_clean)], by = "time", all.x = TRUE)
  x <- merge(x, iw$daily[, .(day, r2_IW = r2, status_IW = status, diel_flag)],
             by = "day", all.x = TRUE)
  x[, Iz := I0 * exp(-cfg$kd * cfg$z_sensor)]
  x[, th := as.numeric(time) / 3600]
  # hourly nep (mmol m-3 h-1), derivative taken within each day only
  x[, dDOdt := centered_diff(th, DO_clean), by = day]
  x[, nep := dDOdt + gas_term - transport]

  daily <- x[, {
    NEP_obs <- mean(nep, na.rm = TRUE) * 24
    f <- if (status_IW[1] == "ok") fit_eq7(nep, Iz, T, cfg$theta) else NULL
    GPP <- R <- NEP <- NA_real_; flag <- "ok"; r2_fit <- NA_real_
    if (is.null(f)) {
      flag <- if (status_IW[1] != "ok") paste0("iw_", status_IW[1]) else "fit_failed"
    } else {
      gpp_i <- f$par[1] * tanh(f$par[2] * Iz / f$par[1])
      r_i   <- f$par[3] * cfg$theta^(T - 20)
      GPP <- mean(gpp_i, na.rm = TRUE) * 24
      R   <- mean(r_i,   na.rm = TRUE) * 24
      NEP <- GPP - R; r2_fit <- f$r2
      if (is.finite(r2_IW[1]) && r2_IW[1] > cfg$r2_IW_max) flag <- "physics_dominated"
      else if (GPP < 0 || R < 0) flag <- "negative_rate"
    }
    if (flag != "ok") {            # Text S2 fallback
      NEP <- NEP_obs
      if (is.finite(NEP)) { GPP <- max(NEP, 0); R <- max(-NEP, 0) }
    }
    list(GPP = GPP, R = R, NEP = NEP, NEP_obs = NEP_obs, r2_fit = r2_fit,
         r2_IW = r2_IW[1], diel_flag = diel_flag[1], flag = flag)
  }, by = day]
  list(hourly = x[], daily = daily[])
}
