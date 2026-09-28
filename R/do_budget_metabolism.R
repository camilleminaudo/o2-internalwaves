###############################################################################
# do_budget_metabolism.R
# Oxygen budget at each sensor depth from high-frequency DO, T and PAR:
#
#   dDO/dt = NEP  +  gas  +  diffusion  +  residual
#
#   NEP = GPP - R     Eq. S7 (Fernández Castro et al. 2021):
#                     nep_i = Pmax tanh(alpha PAR_i / Pmax) - R20 theta^(T_i - 20)
#   gas               air-lake exchange, mixed-layer sensors only: -F_gas / z_th
#                     F_gas = k_O2 (DO - DO_sat), k from wind (Cole & Caraco 1998
#                     or Vachon & Prairie 2013), positive = outgassing
#   diffusion         Kz dDO/dz across the thermocline and between sensors;
#                     Kz from Hondzo & Stefan (1993): Kz = a (N^2)^-0.43
#   residual          everything else (entrainment, lateral advection, ...)
#
# DO is first corrected for internal waves with iw_correction_time_domain.R
# (time-domain method, near-surface rule included).
#
# NEW: SPECTRAL SCALE SEPARATION (run_spectral = TRUE), Supplement Text S1,
# using spectral_iw_correction.R (same code as the Python version):
#   R   = 24 * 2 A24 / t_night,  A24 from the diel (18-28 h) DO variance after
#         IW removal (co-spectrum with T), baseline removal and phase check
#   NEP = d<DO>/dt - <gas> - <diffusion>   (<.> = Gaussian low-pass, 30 d)
#   GPP = R + NEP
# Needs hourly data (aggregated here) and >= 257 h of record. The spectral NEP
# is a fortnight-scale quantity that includes all unknown physics; the daily
# 'departure' = daily nep_obs - spectral NEP is reported for comparison with
# the time-domain residual, but it also contains day-to-day biology.
#
# INPUT (wide csv): time, DO_<z>, T_<z>, PAR_<z>, wind, z_th
#   <z> = sensor depth in m (e.g. DO_1, DO_9.5). PAR can be measured at fewer
#   depths than DO: missing depths are extrapolated with the attenuation
#   coefficient kd fitted on the PAR sensors (or cfg$kd_default).
#   z_th = thermocline / mixed-layer depth (m).
#
# IMPORTANT (identifiability): with fit_days = 1 (paper), GPP, R and the
# residual come from the same day of data. R is nearly constant within a day,
# so any physical flux that persists over the day is absorbed into R (or
# GPP) and the daily residual is ~0 by construction. With fit_days > 1 the
# biological parameters (Pmax, alpha, R20) are shared over fit_days days and
# day-to-day physical events show up in the residual instead. The residual
# then also contains any day-to-day change in the biological parameters.
#
# Units: internally mmol O2 m-3, hours. Outputs in mmol m-3 d-1 (areal:
# mmol m-2 d-1) or mg L-1 d-1 (areal: g m-2 d-1) with out_units = "mg".
###############################################################################

## ---------------------------------------------------------------- CONFIG ----
cfg <- list(
  input_file      = "C:/Projects/myGit/o2-internalwaves/data/Gergal_example_wide.csv",                  # wide csv; NULL -> synthetic test with known truth
  iw_script       = "iw_correction_time_domain.R",
  out_dir         = "results_budget",
  time_col        = "time",
  time_format     = NULL,                  # NULL = auto (ISO, "Y-m-d H:M:S", epoch s)
  tz              = "GMT",
  do_prefix       = "DO", t_prefix = "T", par_prefix = "PAR",
  wind_col        = "wind",
  zth_col         = "z_th",
  do_units        = "mg/L",                # "mg/L" or "mmol/m3"
  out_units       = "mmol",                # "mmol" (mmol m-3 d-1) or "mg" (mg L-1 d-1)
  elev            = 0,                     # m a.s.l., for O2 saturation
  lake_area_km2   = 1,                     # Hondzo & Stefan Kz and Vachon k600
  z_bottom        = NULL,                  # bottom of the deepest layer (m); NULL = deepest sensor + half spacing
  wind_height     = 10,                    # m, height of the anemometer
  k600_model      = "cole",                # "cole" | "vachon"
  kz_min          = 1.4e-7,                # m2 s-1 (molecular heat diffusivity)
  kz_max          = 1e-3,                  # m2 s-1
  kz_const        = NULL,                  # m2 s-1; if set, replaces Hondzo & Stefan
  profile_smooth_h = 24,                   # running mean for gradients (removes IW noise)
  theta           = 1.07,
  par_min         = 1,                     # umol m-2 s-1; below = dark
  kd_default      = 0.5,                   # m-1, used if < 2 PAR sensors
  zeu_frac        = 0.01,                  # euphotic depth = depth of 1% of PAR at the top PAR sensor
  fit_days        = 3,                     # days sharing Pmax, alpha, R20. 1 = paper, but then the
                                           # residual is ~0 by construction (see header)
  r2_IW_max       = 0.75,                  # Text S2: day rejected if IW regression R2 > this (NULL = off)
  day_start_hour  = 0,
  min_coverage    = 0.8,
  iw              = list(surface_rule = "hybrid", surface_depth_max = 3),   # passed to the IW script

  ## spectral scale separation ------------------------------------
  run_spectral    = TRUE,
  spectral_script = "spectral_iw_correction.R",
  spectral        = list(band_integration = "trapz"),   # overrides of spectral_defaults ("rect": less biased)
  spectral_iw_rule = "follow_td",          # "all": step i at every depth (paper);
                                           # "follow_td": no step i at near-surface depths that the
                                           # time-domain hybrid rule leaves uncorrected on most days
  spec_smooth_hw_days = 30,                # low-pass half-width for <DO>, <T> and the NEP budget
  make_plots      = TRUE
)

suppressPackageStartupMessages(library(data.table))

## ------------------------------------------------------------ INPUT -------
stop_if_missing <- function(dt, cols, where = "") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "),
                         "\n  Available: ", paste(names(dt), collapse = ", "), call. = FALSE)
}

parse_time <- function(x, cfg) {
  if (inherits(x, "POSIXct")) return(as.POSIXct(format(x, tz = cfg$tz, usetz = FALSE), tz = cfg$tz))
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = cfg$tz))
  t <- if (!is.null(cfg$time_format)) as.POSIXct(x, format = cfg$time_format, tz = cfg$tz)
       else as.POSIXct(x, tz = cfg$tz, tryFormats = c("%Y-%m-%dT%H:%M:%OSZ", "%Y-%m-%d %H:%M:%OS",
                                                       "%Y-%m-%d %H:%M", "%d/%m/%Y %H:%M"))
  bad <- which(is.na(t) & !is.na(x))
  if (length(bad)) stop(length(bad), " time value(s) not parsed (e.g. '", x[bad[1]],
                        "'). Set cfg$time_format.", call. = FALSE)
  t
}

read_wide <- function(df, cfg) {
  dt <- as.data.table(df)
  stop_if_missing(dt, c(cfg$time_col, cfg$wind_col, cfg$zth_col))
  setnames(dt, c(cfg$time_col, cfg$wind_col, cfg$zth_col), c("time", "wind", "z_th"))
  dt[, time := parse_time(time, cfg)]
  pat <- function(p) sprintf("^%s_([0-9]+(\\.[0-9]+)?)$", p)
  to_long <- function(prefix, name) {
    cols <- grep(pat(prefix), names(dt), value = TRUE)
    if (!length(cols)) return(NULL)
    m <- melt(dt[, c("time", cols), with = FALSE], id.vars = "time", variable.name = "col",
              value.name = name, variable.factor = FALSE)
    m[, depth := as.numeric(sub(pat(prefix), "\\1", col))][, col := NULL]
    m[, (name) := as.numeric(get(name))]
    m
  }
  DOl <- to_long(cfg$do_prefix, "DO"); Tl <- to_long(cfg$t_prefix, "T"); Pl <- to_long(cfg$par_prefix, "PAR")
  if (is.null(DOl)) stop("No ", cfg$do_prefix, "_<depth> columns found.", call. = FALSE)
  if (is.null(Tl))  stop("No ", cfg$t_prefix,  "_<depth> columns found.", call. = FALSE)
  if (is.null(Pl))  stop("No ", cfg$par_prefix, "_<depth> columns found (needed for GPP).", call. = FALSE)
  noT <- setdiff(unique(DOl$depth), unique(Tl$depth))
  if (length(noT)) stop("DO sensor(s) without temperature at depth(s): ", paste(noT, collapse = ", "), call. = FALSE)
  if (length(unique(DOl$depth)) < 2)
    warning("Only one DO depth: no diffusive flux can be estimated.", call. = FALSE)
  if (all(is.na(dt$z_th))) stop("z_th is empty.", call. = FALSE)
  if (cfg$do_units == "mg/L") DOl[, DO := DO * 1000 / 32]           # -> mmol m-3
  else if (cfg$do_units != "mmol/m3") stop("cfg$do_units must be 'mg/L' or 'mmol/m3'.", call. = FALSE)
  long <- merge(DOl, Tl, by = c("time", "depth"), all.x = TRUE)
  setorder(long, depth, time)
  list(long = long, par = Pl[is.finite(PAR)], met = dt[, .(time, wind = as.numeric(wind), z_th = as.numeric(z_th))])
}

## --------------------------------------------------- PHYSICAL HELPERS -----
o2_sat <- function(T, elev = 0) {                      # mmol m-3, Garcia & Gordon (1992), S = 0
  Ts <- log((298.15 - T) / (273.15 + T))
  lnC <- 5.80871 + 3.20291 * Ts + 4.17887 * Ts^2 + 5.10006 * Ts^3 - 9.86643e-2 * Ts^4 + 3.80369 * Ts^5
  exp(lnC) * rho_w(T) / 1000 * (1 - 2.25577e-5 * elev)^5.25588
}
rho_w <- function(T) 1000 * (1 - (T + 288.9414) / (508929.2 * (T + 68.12963)) * (T - 3.9863)^2)
schmidt_o2 <- function(T) 1800.6 - 120.10 * T + 3.7818 * T^2 - 0.047608 * T^3

u10 <- function(u, h) u * (1 + sqrt(0.00114) / 0.41 * log(10 / h))   # neutral log profile
k_o2 <- function(wind, T, cfg) {                        # m h-1
  U <- pmax(u10(wind, cfg$wind_height), 0)
  k600 <- switch(cfg$k600_model,
                 cole   = 2.07 + 0.215 * U^1.7,
                 vachon = 2.51 + 1.48 * U + 0.39 * U * log10(cfg$lake_area_km2),
                 stop("cfg$k600_model must be 'cole' or 'vachon'.", call. = FALSE))
  n <- ifelse(U < 3.7, -2/3, -1/2)
  k600 * (schmidt_o2(T) / 600)^n / 100
}

kz_hs <- function(N2, cfg) {                            # m2 h-1
  if (!is.null(cfg$kz_const)) return(rep(cfg$kz_const * 3600, length(N2)))
  a <- 8.17e-4 * cfg$lake_area_km2^0.56                 # cm2 s-1
  kz <- a * pmax(N2, 1e-7)^(-0.43) * 1e-4               # m2 s-1
  pmin(pmax(kz, cfg$kz_min), cfg$kz_max) * 3600
}

run_mean <- function(x, n) {                            # centred, NA-tolerant
  if (n <= 1) return(x)
  w <- rep(1, n); ok <- is.finite(x)
  num <- stats::filter(ifelse(ok, x, 0), w, sides = 2)
  den <- stats::filter(as.numeric(ok), w, sides = 2)
  r <- as.numeric(num / den)
  r[is.na(r)] <- x[is.na(r)]                            # edges: keep raw
  r
}

centered_diff <- function(x, y) {
  out <- rep(NA_real_, length(y)); ok <- which(is.finite(x) & is.finite(y))
  if (length(ok) < 2) return(out)
  x0 <- x[ok]; y0 <- y[ok]; n <- length(ok); d <- numeric(n)
  d[1] <- (y0[2] - y0[1]) / (x0[2] - x0[1]); d[n] <- (y0[n] - y0[n - 1]) / (x0[n] - x0[n - 1])
  if (n > 2) d[2:(n - 1)] <- (y0[3:n] - y0[1:(n - 2)]) / (x0[3:n] - x0[1:(n - 2)])
  out[ok] <- d; out
}

## ------------------------------------------------------------ PAR ---------
# PAR at each DO depth: measured if available, otherwise extrapolated from the
# nearest PAR sensor with the daily kd fitted on all PAR sensors.
par_at_depths <- function(par, depths, times, cfg) {
  par <- copy(par)[time %in% times]
  par[, day := day_of(time, cfg)]
  kdt <- par[PAR > cfg$par_min, if (.N >= 2 && var(depth) > 0)
               .(kd = -cov(depth, log(PAR)) / var(depth)) else .(kd = NA_real_), by = time]
  kdt[, day := day_of(time, cfg)]
  kd_day <- kdt[is.finite(kd), .(kd = median(kd)), by = day]
  days <- unique(day_of(times, cfg))
  kd_day <- merge(data.table(day = days), kd_day, by = "day", all.x = TRUE)
  if (any(is.na(kd_day$kd))) {
    if (all(is.na(kd_day$kd))) message("kd: fewer than 2 PAR sensors with light; using kd_default = ", cfg$kd_default)
    kd_day[is.na(kd), kd := cfg$kd_default]
  }
  pdep <- sort(unique(par$depth))
  out <- rbindlist(lapply(depths, function(z) {
    zr <- pdep[which.min(abs(pdep - z))]
    p  <- par[depth == zr, .(time, PAR_ref = PAR)]
    d  <- data.table(time = times, depth = z)[p, on = "time", PAR_ref := i.PAR_ref]
    d[, day := day_of(time, cfg)][kd_day, on = "day", kd := i.kd]
    d[, PAR := if (z == zr) PAR_ref else PAR_ref * exp(-kd * (z - zr))]
    d[, .(time, depth, PAR)]
  }))
  top <- min(pdep)
  kd_day[, z_eu := top + log(1 / cfg$zeu_frac) / kd]
  list(par = out, kd_day = kd_day)
}

day_of <- function(t, cfg) as.Date(as.POSIXct(as.numeric(t) - cfg$day_start_hour * 3600,
                                             origin = "1970-01-01", tz = cfg$tz), tz = cfg$tz)

## ------------------------------------------------------ PHYSICAL FLUXES ---
# Returns, per time and depth: gas and diffusion contributions to dDO/dt
# (mmol m-3 h-1), layer thickness h (m), Kz and fluxes at interfaces.
physical_terms <- function(hr, met, cfg) {
  depths <- sort(unique(hr$depth)); nd <- length(depths)
  W <- function(v) as.matrix(dcast(hr, time ~ depth, value.var = v)[, -1])
  times <- sort(unique(hr$time)); nt <- length(times)
  step_h <- median(diff(as.numeric(times))) / 3600
  ns <- max(1, round(cfg$profile_smooth_h / step_h))
  Ts  <- apply(W("T"), 2, run_mean, n = ns)
  DOs <- apply(W("DO_clean"), 2, run_mean, n = ns)
  DOr <- W("DO"); Tr <- W("T")
  zth <- met$z_th
  inML <- outer(zth, depths, ">")
  z_bot <- if (!is.null(cfg$z_bottom)) cfg$z_bottom else
           if (nd > 1) depths[nd] + (depths[nd] - depths[nd - 1]) / 2 else depths[nd] + 1

  # interfaces between consecutive sensors
  Fint <- Kzint <- matrix(0, nt, max(nd - 1, 1)); is_th <- matrix(FALSE, nt, max(nd - 1, 1))
  if (nd > 1) for (k in 1:(nd - 1)) {
    a <- k; b <- k + 1; dz <- depths[b] - depths[a]
    N2 <- 9.81 / 1000 * (rho_w(Ts[, b]) - rho_w(Ts[, a])) / dz
    Kz <- kz_hs(N2, cfg)
    Fd <- Kz * (DOs[, a] - DOs[, b]) / dz              # downward diffusive flux, mmol m-2 h-1
    both_ml <- inML[, a] & inML[, b]
    Fd[both_ml] <- 0; Kz[both_ml] <- NA
    Fint[, k] <- Fd; Kzint[, k] <- Kz; is_th[, k] <- inML[, a] & !inML[, b]
  }
  Fint[!is.finite(Fint)] <- NA

  # layer thickness: ML sensors share the ML box; below-ML sensors: interface to interface
  h <- matrix(NA_real_, nt, nd); diffc <- matrix(0, nt, nd)
  n_ml <- rowSums(inML)
  F_th <- if (nd > 1) rowSums(Fint * is_th, na.rm = TRUE) else rep(0, nt)
  for (b in seq_len(nd)) {
    ml <- inML[, b]
    top <- if (b == 1) zth else ifelse(inML[, b - 1], zth, (depths[b - 1] + depths[b]) / 2)
    bot <- if (b == nd) rep(z_bot, nt) else rep((depths[b] + depths[b + 1]) / 2, nt)
    hb <- ifelse(ml, zth / pmax(n_ml, 1), bot - top)
    hb[hb <= 0] <- NA
    h[, b] <- hb
    up <- if (b > 1) Fint[, b - 1] else rep(0, nt)          # flux entering from above
    dn <- if (b < nd) Fint[, b] else rep(0, nt)             # flux leaving downward
    below <- (up - dn) / hb
    diffc[, b] <- ifelse(ml, -F_th / zth, below)
  }
  # gas exchange with the shallowest sensor (must be in the mixed layer)
  wind <- met$wind
  kO2 <- k_o2(wind, Tr[, 1], cfg)
  Fgas <- kO2 * (DOr[, 1] - o2_sat(Tr[, 1], cfg$elev))    # mmol m-2 h-1, + = outgassing
  Fgas[!inML[, 1]] <- NA                                   # surface sensor below z_th: ML not sampled
  gasc <- ifelse(inML, -matrix(Fgas, nt, nd) / zth, 0)
  gasc[is.na(gasc)] <- 0
  if (any(!inML[, 1], na.rm = TRUE))
    message(sum(!inML[, 1], na.rm = TRUE), " time steps with z_th above the shallowest sensor: ",
            "gas exchange not attributed to any sensor then.")

  long <- function(M, name) {
    d <- as.data.table(M); setnames(d, as.character(depths)); d[, time := times]
    m <- melt(d, id.vars = "time", variable.name = "depth", value.name = name, variable.factor = FALSE)
    m[, depth := as.numeric(depth)]
  }
  out <- Reduce(function(x, y) merge(x, y, by = c("time", "depth")),
                list(long(gasc, "gas"), long(diffc, "diffusion"), long(h, "h"), long(inML * 1, "in_ml")))
  kz_out <- if (nd > 1) {
    d <- as.data.table(Kzint / 3600)
    setnames(d, sprintf("Kz_%g_%g", depths[-nd], depths[-1])); cbind(time = times, d)
  } else NULL
  list(terms = out, Fgas = data.table(time = times, F_gas = Fgas, k_O2 = kO2), kz = kz_out, z_bottom = z_bot)
}

## ------------------------------------------------------- METABOLISM FIT ---
fit_eq7 <- function(nep, PAR, T, theta, par_min) {
  ok <- is.finite(nep) & is.finite(PAR) & is.finite(T)
  if (sum(ok) < 8) return(NULL)
  nep <- nep[ok]; PAR <- PAR[ok]; T <- T[ok]; g <- theta^(T - 20)
  if (max(PAR) < par_min) {                                # aphotic: respiration only
    R20 <- -sum(nep * g) / sum(g^2)
    return(list(par = c(Pmax = 0, alpha = 0, R20 = R20), aphotic = TRUE,
                r2 = 1 - sum((nep + R20 * g)^2) / sum((nep - mean(nep))^2)))
  }
  pred <- function(p) p[1] * tanh(p[2] * PAR / p[1]) - p[3] * g
  sse <- function(p) sum((nep - pred(p))^2)
  dark <- PAR < par_min
  R0 <- if (any(dark)) max(-mean(nep[dark] / g[dark]), 1e-3) else max(-min(nep), 1e-3)
  P0 <- max(max(nep) + R0, 1e-2); a0 <- P0 / max(median(PAR[!dark]), 1)
  best <- NULL
  for (s in list(c(P0, a0, R0), c(2 * P0, a0 / 2, R0), c(P0 / 2, 2 * a0, R0))) {
    f <- try(optim(s, sse, method = "L-BFGS-B", lower = c(1e-6, -Inf, -Inf), upper = c(100 * P0, Inf, Inf)),
             silent = TRUE)
    if (!inherits(f, "try-error") && is.finite(f$value) && (is.null(best) || f$value < best$value)) best <- f
  }
  if (is.null(best) || best$convergence == 1) return(NULL)
  list(par = setNames(best$par, c("Pmax", "alpha", "R20")), aphotic = FALSE,
       r2 = 1 - best$value / sum((nep - mean(nep))^2))
}

## ------------------------------------------------------------ PIPELINE ----
run_budget <- function(df, cfg, iw_env, sp_env = NULL) {
  inp <- read_wide(df, cfg)

  # regular grid
  step <- median(diff(sort(unique(as.numeric(inp$long$time)))))
  grid <- seq(min(inp$long$time), max(inp$long$time), by = step)
  depths <- sort(unique(inp$long$depth))
  hr <- CJ(time = grid, depth = depths)[inp$long, on = .(time, depth), `:=`(DO = i.DO, T = i.T)]
  met <- data.table(time = grid)[inp$met, on = "time", `:=`(wind = i.wind, z_th = i.z_th)]
  for (v in c("wind", "z_th")) {                           # fill short gaps linearly
    ok <- is.finite(met[[v]])
    if (sum(ok) < 2) stop("Column ", v, " has fewer than 2 values.", call. = FALSE)
    met[, (v) := approx(as.numeric(time)[ok], get(v)[ok], as.numeric(time), rule = 2)$y]
  }

  # 1) internal-wave correction (current version of the time-domain script)
  iw_cfg <- modifyList(iw_env$cfg, c(list(input_file = NULL, depth_col = "depth", make_plots = FALSE,
                                          tz = cfg$tz, day_start_hour = cfg$day_start_hour,
                                          min_coverage = cfg$min_coverage), cfg$iw))
  message("--- Internal-wave correction ---")
  iw <- iw_env$correct_internal_waves(iw_env$prepare_input(hr[, .(time, DO, T, depth)], iw_cfg), iw_cfg)
  hr[iw$data, on = .(depth, time), `:=`(DO_clean = i.DO_clean, day = i.day)]
  hr[is.na(day), day := day_of(time, cfg)]

  # 2) light
  pa <- par_at_depths(inp$par, depths, grid, cfg)
  hr[pa$par, on = .(time, depth), PAR := i.PAR]

  # 3) physical terms
  ph <- physical_terms(hr, met, cfg)
  hr <- merge(hr, ph$terms, by = c("time", "depth"))
  hr[met, on = "time", z_th := i.z_th]
  setorder(hr, depth, time)

  # 4) observed rate and nep = dDO/dt - gas - diffusion (= NEP + residual)
  hr[, th := as.numeric(time) / 3600]
  hr[, dDOdt := centered_diff(th, DO_clean), by = .(depth, day)]
  hr[, nep_obs := dDOdt - gas - diffusion]

  # 5) usable days (Text S2) and Eq. S7 fits
  idl <- iw$daily[, .(depth, day, coverage, r2_IW = r2, corrected, near_surface, iw_decision)]
  idl[, usable := coverage >= cfg$min_coverage &
                  !(corrected & !is.null(cfg$r2_IW_max) & is.finite(r2_IW) & r2_IW > if (is.null(cfg$r2_IW_max)) Inf else cfg$r2_IW_max) &
                  !(!corrected & !near_surface)]
  hr[idl, on = .(depth, day), usable := i.usable]
  hr[is.na(usable), usable := FALSE]
  d0 <- min(hr$day)
  hr[, grp := as.integer(floor(as.numeric(day - d0) / cfg$fit_days))]

  fits <- hr[, {
    f <- fit_eq7(nep_obs[usable], PAR[usable], T[usable], cfg$theta, cfg$par_min)
    if (is.null(f)) list(Pmax = NA_real_, alpha = NA_real_, R20 = NA_real_, r2_fit = NA_real_, aphotic = NA)
    else list(Pmax = f$par[["Pmax"]], alpha = f$par[["alpha"]], R20 = f$par[["R20"]], r2_fit = f$r2, aphotic = f$aphotic)
  }, by = .(depth, grp)]
  hr[fits, on = .(depth, grp), `:=`(Pmax = i.Pmax, alpha = i.alpha, R20 = i.R20)]
  hr[, gpp_mod := fifelse(Pmax > 0, Pmax * tanh(alpha * PAR / Pmax), 0)]
  hr[, r_mod := R20 * cfg$theta^(T - 20)]

  # 6) daily budget (rates x 24 -> per day)
  daily <- hr[, .(n = sum(is.finite(nep_obs)), usable = usable[1], z_th = mean(z_th),
                  in_ml = mean(in_ml) > 0.5, h = mean(h, na.rm = TRUE),
                  dDOdt = mean(dDOdt, na.rm = TRUE) * 24, gas = mean(gas, na.rm = TRUE) * 24,
                  diffusion = mean(diffusion, na.rm = TRUE) * 24,
                  nep_obs = mean(nep_obs, na.rm = TRUE) * 24,
                  GPP = mean(gpp_mod, na.rm = TRUE) * 24, R = mean(r_mod, na.rm = TRUE) * 24),
              by = .(depth, day, grp)]
  daily[, NEP := GPP - R]
  daily[, flag := fifelse(!usable, "fallback:iw", fifelse(!is.finite(GPP) | !is.finite(R), "fallback:fit_failed",
                  fifelse(GPP < 0 | R < 0, "fallback:negative_rate", "ok")))]
  daily[flag != "ok", `:=`(NEP = nep_obs, GPP = pmax(nep_obs, 0), R = pmax(-nep_obs, 0))]   # Text S2
  daily[, residual := nep_obs - NEP]                     # unknown physics (+ model error)
  daily[idl, on = .(depth, day), `:=`(r2_IW = i.r2_IW, iw_decision = i.iw_decision)]
  daily[fits, on = .(depth, grp), `:=`(r2_fit = i.r2_fit, Pmax = i.Pmax, alpha = i.alpha, R20 = i.R20)]
  daily[pa$kd_day, on = "day", `:=`(kd = i.kd, z_eu = i.z_eu)]
  daily[, productive := in_ml | depth <= z_eu]
  daily[, layer := fifelse(in_ml, "mixed layer", "below mixed layer")]
  daily[, closure := dDOdt - (NEP + gas + diffusion + residual)]  # ~0 by construction (check)

  # hourly residual (fitted days only)
  hr[daily, on = .(depth, day), flag := i.flag]
  hr[, residual := fifelse(flag == "ok", nep_obs - (gpp_mod - r_mod), NA_real_)]

  # 7) productive layer, areal (mmol m-2 d-1): rate x layer thickness
  areal <- daily[productive == TRUE, .(
    n_sensors = .N, depths = paste(depth, collapse = "+"), thickness = sum(h),
    GPP = sum(GPP * h), R = sum(R * h), NEP = sum(NEP * h), dDOdt = sum(dDOdt * h),
    gas = sum(gas * h), diffusion = sum(diffusion * h), residual = sum(residual * h),
    all_ok = all(flag == "ok")), by = day]

  res <- list(hourly = hr, daily = daily, areal = areal, fits = fits, iw = iw, kz = ph$kz,
       gas = ph$Fgas, kd = pa$kd_day, z_bottom = ph$z_bottom)
  res$par_raw <- inp$par                                  ## NEW: kept for the spectral method
  if (isTRUE(cfg$run_spectral) && !is.null(sp_env)) res$spec <- run_spectral(res, cfg, sp_env)
  res
}

## ------------------------------------------- NEW: SPECTRAL METHOD --------
lowpass_gauss <- function(x, step_days, hw_days) {      # Gaussian low-pass of DO_budget.py
  k <- seq(-floor(hw_days / step_days), floor(hw_days / step_days))
  w <- exp(-((k * step_days)^2) / (0.5 * hw_days)^2)
  ok <- is.finite(x); pad <- length(k) %/% 2
  num <- stats::filter(c(rep(0, pad), ifelse(ok, x, 0), rep(0, pad)), w, sides = 2)
  den <- stats::filter(c(rep(0, pad), as.numeric(ok), rep(0, pad)), w, sides = 2)
  r <- as.numeric(num / den)[(pad + 1):(pad + length(x))]
  r[!is.finite(r)] <- NA; r
}

# night length (h) per day from the shallowest PAR sensor
night_length <- function(par, cfg) {
  p <- par[depth == min(depth)]
  st <- median(diff(sort(unique(as.numeric(p$time))))) / 3600
  p[, .(t_night = 24 - sum(PAR > cfg$par_min) * st), by = .(day = day_of(time, cfg))]
}

run_spectral <- function(out, cfg, sp_env) {
  hr <- out$hourly
  hr[, t_h := as.POSIXct(floor(as.numeric(time) / 3600) * 3600, origin = "1970-01-01", tz = cfg$tz)]
  h1 <- hr[, .(DO = mean(DO, na.rm = TRUE), T = mean(T, na.rm = TRUE), I0 = mean(PAR, na.rm = TRUE),
               gas_term = -mean(gas, na.rm = TRUE),        # spectral code: + = O2 lost to the air
               diffusion = mean(diffusion, na.rm = TRUE), in_ml = mean(in_ml, na.rm = TRUE) > 0.5),
           by = .(depth, time = t_h)]
  hr[, t_h := NULL]
  for (v in c("DO", "T", "I0", "gas_term", "diffusion")) h1[is.nan(get(v)), (v) := NA_real_]
  grid <- CJ(depth = unique(h1$depth), time = seq(min(h1$time), max(h1$time), by = 3600))
  h1 <- h1[grid, on = .(depth, time)]
  h1[is.na(in_ml), in_ml := FALSE]; h1[is.na(gas_term), gas_term := 0]
  if (uniqueN(h1$time) < 257) {
    message("Spectral method skipped: record shorter than one 257-h window.")
    return(NULL)
  }
  tn <- night_length(out$par_raw, cfg)
  sp <- modifyList(sp_env$spectral_defaults, cfg$spectral)
  iwd <- out$iw$daily[, .(ns = all(near_surface), frac_corr = mean(corrected)), by = depth]

  res <- rbindlist(lapply(sort(unique(h1$depth)), function(z) {
    d <- h1[depth == z]
    hw <- cfg$spec_smooth_hw_days
    d[, `:=`(sDO = lowpass_gauss(DO, 1 / 24, hw), sT = lowpass_gauss(T, 1 / 24, hw))]
    sgas <- lowpass_gauss(d$gas_term, 1 / 24, hw); sdif <- lowpass_gauss(d$diffusion, 1 / 24, hw)
    d[, NEP_budget := centered_diff(as.numeric(time) / 86400, sDO) + 24 * sgas - 24 * sdif]
    d[, day := day_of(time, cfg)][tn, on = "day", t_night := i.t_night]
    d[!is.finite(t_night), t_night := median(tn$t_night, na.rm = TRUE)]
    spz <- sp
    w <- iwd[depth == z]
    if (cfg$spectral_iw_rule == "follow_td" && nrow(w) && w$ns && w$frac_corr < 0.5) spz$remove_iw <- FALSE
    r <- copy(sp_env$spectral_metabolism(d, spz, z))    # copy: avoids a data.table selfref warning
    r[, `:=`(depth = z, iw_removed = spz$remove_iw)]
  }))
  res[, day := day_of(time, cfg)]
  daily <- res[, .(GPP = mean(GPP), R = mean(R), NEP = mean(NEP), A24 = mean(A24),
                   frac_no_diel = mean(var_do24 == 0), phase_h = median(phase_h),
                   iw_removed = iw_removed[1]), by = .(depth, day)]
  daily[out$daily, on = .(depth, day), `:=`(nep_obs = i.nep_obs, h = i.h, productive = i.productive)]
  daily[, departure := nep_obs - NEP]                    # daily budget - fortnight NEP
  areal <- daily[productive == TRUE, .(GPP = sum(GPP * h), R = sum(R * h), NEP = sum(NEP * h),
                                       departure = sum(departure * h)), by = day]
  list(hourly = res, daily = daily, areal = areal, settings = sp)
}

convert_units <- function(out, cfg) {
  if (cfg$out_units != "mg") return(out)
  f <- 32 / 1000                                          # mmol m-3 -> mg L-1 ; mmol m-2 -> g m-2
  rate_cols <- c("dDOdt", "gas", "diffusion", "nep_obs", "GPP", "R", "NEP", "residual", "closure")
  out$daily[, (intersect(rate_cols, names(out$daily))) := lapply(.SD, `*`, f), .SDcols = intersect(rate_cols, names(out$daily))]
  out$areal[, (intersect(rate_cols, names(out$areal))) := lapply(.SD, `*`, f), .SDcols = intersect(rate_cols, names(out$areal))]
  if (!is.null(out$spec)) for (tb in c("daily", "areal")) {          ## NEW
    cc <- intersect(c("GPP", "R", "NEP", "A24", "nep_obs", "departure"), names(out$spec[[tb]]))
    out$spec[[tb]][, (cc) := lapply(.SD, `*`, f), .SDcols = cc]
  }
  out
}

## --------------------------------------------------------------- PLOTS ----
plot_budget <- function(out, cfg, truth = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  library(ggplot2)
  u  <- if (cfg$out_units == "mg") "mg~L^-1~d^-1" else "mmol~m^-3~d^-1"
  ua <- if (cfg$out_units == "mg") "g~m^-2~d^-1" else "mmol~m^-2~d^-1"
  d <- melt(out$daily, id.vars = c("depth", "day"),
            measure.vars = c("dDOdt", "NEP", "gas", "diffusion", "residual"))
  p1 <- ggplot(d, aes(day, value, colour = variable)) + geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_line() + geom_point(size = 0.8) + facet_wrap(~depth, ncol = 1, scales = "free_y", labeller = label_both) +
    labs(x = NULL, y = parse(text = u), colour = "budget term",
         title = "Daily DO budget: dDO/dt = NEP + gas + diffusion + residual") + theme_bw()
  m <- melt(out$daily, id.vars = c("depth", "day", "flag"), measure.vars = c("GPP", "R", "NEP"))
  p2 <- ggplot(m, aes(day, value, colour = variable)) + geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_line() + geom_point(aes(shape = flag != "ok"), size = 1.5) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c("fitted", "fallback"), name = NULL) +
    facet_wrap(~depth, ncol = 1, scales = "free_y", labeller = label_both) +
    labs(x = NULL, y = parse(text = u), colour = NULL) + theme_bw()
  if (!is.null(truth)) {
    tl <- melt(truth, id.vars = c("depth", "day"), measure.vars = c("GPP", "R", "NEP"))
    p2 <- p2 + geom_line(data = tl, aes(day, value, colour = variable), linetype = "dashed", inherit.aes = FALSE) +
      labs(caption = "dashed = truth")
  }
  a <- melt(out$areal, id.vars = "day", measure.vars = c("GPP", "R", "NEP", "gas", "diffusion", "residual"))
  p3 <- ggplot(a, aes(day, value, colour = variable)) + geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_line() + geom_point(size = 0.8) +
    labs(x = NULL, y = parse(text = ua), colour = NULL, title = "Productive layer (areal)") + theme_bw()
  ## NEW: time domain vs spectral
  p4 <- NULL
  if (!is.null(out$spec)) {
    cm <- rbind(melt(out$daily[, .(depth, day, GPP, R, NEP)], id.vars = c("depth", "day"))[, method := "time domain"],
                melt(out$spec$daily[, .(depth, day, GPP, R, NEP)], id.vars = c("depth", "day"))[, method := "spectral"])
    if (!is.null(truth))
      cm <- rbind(cm, melt(truth[, .(depth, day, GPP, R, NEP)], id.vars = c("depth", "day"))[, method := "truth"])
    p4 <- ggplot(cm, aes(day, value, colour = method)) + geom_hline(yintercept = 0, linewidth = 0.2) +
      geom_line() + facet_grid(depth ~ variable, scales = "free_y", labeller = label_both) +
      scale_colour_manual(values = c(`time domain` = "#1b9e77", spectral = "#d95f02", truth = "black")) +
      labs(x = NULL, y = parse(text = u), colour = NULL) + theme_bw() + theme(legend.position = "bottom")
  }
  list(budget = p1, rates = p2, areal = p3, methods = p4)
}

## --------------------------------------------------- SYNTHETIC TEST ------
# Three-box lake (mixed layer 0-z_th, metalimnion z_th-14 m, hypolimnion 14-25 m)
# with sensors at 1, 9, 19 m, known GPP/R, gas exchange, Hondzo-Stefan
# diffusion, internal waves at 9 and 19 m, and UNKNOWN physics:
#   - lateral inflow of low-DO water into the mixed layer (day 4, 12 h)
#   - entrainment when the thermocline deepens from 5 to 8 m (days 7-8)
#   - intrusion of oxygenated water at 9 m (day 10)
make_synthetic <- function(cfg, days = 12, step_min = 30, seed = 7, event_offset = 0) {
  ## NEW: event_offset (days) moves the events later, for longer records
  set.seed(seed)
  time <- seq(as.POSIXct("2022-04-12", tz = "UTC"), by = step_min * 60, length.out = days * 1440 / step_min)
  nt <- length(time); th <- (seq_len(nt) - 1) * step_min / 60; dday <- th / 24
  dth <- step_min / 60
  red <- function(tau, sd) { e <- rnorm(nt); x <- numeric(nt)
    phi <- exp(-dth / tau); for (i in 2:nt) x[i] <- phi * x[i - 1] + sqrt(1 - phi^2) * e[i]; sd * x }
  hod <- as.numeric(format(time, "%H")) + as.numeric(format(time, "%M")) / 60
  cloud <- rep(runif(days, 0.4, 1), each = 1440 / step_min)
  PAR0 <- pmax(0, 1800 * cloud * sin(pi * (hod - 6.5) / 13)); PAR0[hod < 6.5 | hod > 19.5] <- 0
  kd <- 0.3; zs <- c(1, 9, 19)
  PARz <- sapply(zs, function(z) PAR0 * exp(-kd * z))
  eo <- event_offset
  z_th <- ifelse(dday < 6 + eo, 5, ifelse(dday < 7 + eo, 5 + 3 * (dday - 6 - eo), 8))
  wind <- pmax(0.3, 3 + 1.5 * sin(2 * pi * (hod - 14) / 24) + red(12, 1) + ifelse(dday >= 6 + eo & dday < 7 + eo, 5, 0))
  T_box <- cbind(17 + 0.08 * dday + 0.5 * sin(2 * pi * (hod - 10) / 24), 12.5 + 0.02 * dday, 10)
  bio <- list(Pmax = c(1.5, 0.5, 0.1), alpha = c(0.004, 0.006, 0.002), R20 = c(0.9, 0.25, 0.15))
  inflow_ml <- ifelse(dday >= 3.5 + eo & dday < 4 + eo, -0.4, 0)     # mmol m-3 h-1
  intru_m   <- ifelse(dday >= 9 + eo & dday < 10 + eo, 0.2, 0)

  C <- matrix(NA, nt, 3); C[1, ] <- c(300, 330, 250)
  gpp <- r <- unk <- gasT <- difT <- matrix(0, nt, 3)
  sub <- 5; hs <- dth / sub
  for (i in seq_len(nt)) {
    g <- sapply(1:3, function(b) bio$Pmax[b] * tanh(bio$alpha[b] * PARz[i, b] / bio$Pmax[b]))
    rr <- bio$R20 * cfg$theta^(T_box[i, ] - 20)
    gpp[i, ] <- g; r[i, ] <- rr
    if (i == nt) break
    x <- C[i, ]; zt0 <- z_th[i]; zt1 <- z_th[i + 1]
    acc_unk <- acc_gas <- acc_dif <- numeric(3)
    for (s in 1:sub) {
      zt <- zt0 + (zt1 - zt0) * (s - 1) / sub
      N2a <- 9.81 / 1000 * (rho_w(T_box[i, 2]) - rho_w(T_box[i, 1])) / 8
      N2b <- 9.81 / 1000 * (rho_w(T_box[i, 3]) - rho_w(T_box[i, 2])) / 10
      Fa <- kz_hs(N2a, cfg) * (x[1] - x[2]) / 8; Fb <- kz_hs(N2b, cfg) * (x[2] - x[3]) / 10
      Fg <- k_o2(wind[i], T_box[i, 1], cfg) * (x[1] - o2_sat(T_box[i, 1], cfg$elev))
      gas <- c(-Fg / zt, 0, 0); dif <- c(-Fa / zt, (Fa - Fb) / (14 - zt), Fb / 11)
      unk_r <- c(inflow_ml[i], intru_m[i], 0)
      x <- x + (g - rr + gas + dif + unk_r) * hs
      dz <- (zt1 - zt0) / sub                               # entrainment into the mixed layer
      if (dz > 0) { e <- (x[2] - x[1]) * dz / (zt + dz); x[1] <- x[1] + e; unk_r[1] <- unk_r[1] + e / hs }
      acc_unk <- acc_unk + unk_r / sub; acc_gas <- acc_gas + gas / sub; acc_dif <- acc_dif + dif / sub
    }
    C[i + 1, ] <- x; unk[i, ] <- acc_unk; gasT[i, ] <- acc_gas; difT[i, ] <- acc_dif
  }
  zeta <- 0.4 * sin(2 * pi * th / 7.3) + 0.25 * sin(2 * pi * th / 11.1 + 1) + red(4, 0.2)
  obs <- data.table(time = time,
    DO_1 = C[, 1] + rnorm(nt, 0, 0.3), DO_9 = C[, 2] - 1.5 * zeta + rnorm(nt, 0, 0.3),
    DO_19 = C[, 3] - 1.0 * 0.6 * zeta + rnorm(nt, 0, 0.3),
    T_1 = T_box[, 1] + rnorm(nt, 0, 0.01), T_9 = T_box[, 2] - 0.8 * zeta + rnorm(nt, 0, 0.01),
    T_19 = T_box[, 3] - 0.2 * 0.6 * zeta + rnorm(nt, 0, 0.01),
    PAR_1 = PARz[, 1], PAR_9 = PARz[, 2], wind = wind, z_th = z_th)
  obs[, (c("DO_1", "DO_9", "DO_19")) := lapply(.SD, function(v) v * 32 / 1000), .SDcols = c("DO_1", "DO_9", "DO_19")]
  tr <- rbindlist(lapply(1:3, function(b) data.table(depth = zs[b], day = as.Date(time, tz = "UTC"),
        GPP = gpp[, b], R = r[, b], unknown = unk[, b], gas = gasT[, b], diffusion = difT[, b])))
  truth <- tr[, lapply(.SD, function(v) mean(v) * 24), by = .(depth, day)]
  truth[, NEP := GPP - R]
  list(obs = obs, truth = truth)
}

## ------------------------------------------------------------------ RUN ----
## FIX: load a companion script without running its RUN block. The option is
## restored on exit even if sourcing fails, so it cannot stay TRUE in the session
## (which silently disabled the RUN blocks of all scripts).
load_as_library <- function(path) {
  old <- options(iw.library_mode = TRUE)
  on.exit(options(old))
  env <- new.env()
  sys.source(path, envir = env)
  env
}

if (!isTRUE(getOption("iw.library_mode", FALSE))) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) setwd(dirname(p))
  }
  if (!file.exists(cfg$iw_script)) stop("IW script not found: ", cfg$iw_script, call. = FALSE)
  iw_env <- load_as_library(cfg$iw_script)             ## FIX: option always restored
  sp_env <- NULL
  if (isTRUE(cfg$run_spectral)) {
    if (!file.exists(cfg$spectral_script)) stop("Spectral script not found: ", cfg$spectral_script, call. = FALSE)
    sp_env <- load_as_library(cfg$spectral_script)
  }

  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  truth <- NULL
  if (is.null(cfg$input_file)) {
    message("No input_file set: running the synthetic test (known truth).")
    cfg$z_bottom <- 25; cfg$do_units <- "mg/L"
    ## NEW: 40 days so that the spectral method has room (events on days 18-25)
    ev_off <- 15
    syn <- make_synthetic(cfg, days = 40, event_offset = ev_off); df <- syn$obs; truth <- syn$truth
  } else {
    if (!file.exists(cfg$input_file)) stop("Input file not found: ", cfg$input_file, call. = FALSE)
    df <- fread(cfg$input_file)
  }

  out <- run_budget(df, cfg, iw_env, sp_env)

  if (!is.null(truth)) {
    cmp_one <- function(o, lab) {
      m <- merge(o$daily, truth, by = c("depth", "day"), suffixes = c("", "_true"))
      m[, .(fit = lab, GPP = mean(GPP), GPP_true = mean(GPP_true), R = mean(R), R_true = mean(R_true),
            gas = mean(gas), gas_true = mean(gas_true), diff = mean(diffusion), diff_true = mean(diffusion_true),
            r_resid_vs_unknown = if (sd(residual) > 0 && sd(unknown) > 0) cor(residual, unknown) else NA_real_,
            rmse_resid = sqrt(mean((residual - unknown)^2)), fallback_days = sum(flag != "ok")), by = depth]
    }
    out1 <- if (cfg$fit_days == 1) out else run_budget(df, modifyList(cfg, list(fit_days = 1, run_spectral = FALSE)), iw_env)
    out3 <- if (cfg$fit_days == 3) out else run_budget(df, modifyList(cfg, list(fit_days = 3, run_spectral = FALSE)), iw_env)
    cat("\n==== Time domain: estimates vs truth (means over days, mmol m-3 d-1) ====\n")
    print(rbind(cmp_one(out1, "fit_days = 1"), cmp_one(out3, "fit_days = 3")), digits = 2)

    ## NEW: time domain vs spectral variants, interior days only
    variants <- list()
    variants[["TD fit_days=1"]] <- out1$daily[, .(depth, day, GPP, R, NEP, resid = residual)]
    variants[["TD fit_days=3"]] <- out3$daily[, .(depth, day, GPP, R, NEP, resid = residual)]
    if (!is.null(sp_env)) {
      spv <- list(`SP paper (trapz, step i everywhere)` = list(spectral = list(band_integration = "trapz"), spectral_iw_rule = "all"),
                  `SP trapz, follow_td`                  = list(spectral = list(band_integration = "trapz"), spectral_iw_rule = "follow_td"),
                  `SP rect, follow_td`                   = list(spectral = list(band_integration = "rect"),  spectral_iw_rule = "follow_td"))
      for (nm in names(spv)) {
        sv <- run_spectral(out1, modifyList(cfg, spv[[nm]]), sp_env)
        variants[[nm]] <- sv$daily[, .(depth, day, GPP, R, NEP, resid = departure)]
      }
    }
    allv <- rbindlist(variants, idcol = "method")
    allv <- merge(allv, truth, by = c("depth", "day"), suffixes = c("", "_true"))
    dmin <- min(truth$day) + 7; dmax <- max(truth$day) - 7
    allv <- allv[day >= dmin & day <= dmax]
    allv[, week := as.integer(as.numeric(day - dmin) %/% 7)]
    wk <- allv[, lapply(.SD, mean), by = .(method, depth, week),
               .SDcols = c("GPP", "GPP_true", "R", "R_true", "NEP", "NEP_true")]
    summ <- allv[, .(GPP_bias = mean(GPP - GPP_true), R_bias = mean(R - R_true), NEP_bias = mean(NEP - NEP_true),
                     NEP_rmse_daily = sqrt(mean((NEP - NEP_true)^2)),
                     r_resid_unknown = if (sd(resid) > 0 && sd(unknown) > 0) cor(resid, unknown) else NA_real_),
                 by = .(method, depth)]
    summ <- merge(summ, wk[, .(GPP_rmse_week = sqrt(mean((GPP - GPP_true)^2)),
                               R_rmse_week = sqrt(mean((R - R_true)^2))), by = .(method, depth)],
                  by = c("method", "depth"))
    tm <- truth[day >= dmin & day <= dmax, .(GPP_true = mean(GPP), R_true = mean(R), NEP_true = mean(NEP)), by = depth]
    cat("\n==== Time domain vs spectral, days", format(dmin), "to", format(dmax), "(mmol m-3 d-1) ====\n")
    cat("Truth means:\n"); print(tm, digits = 3)
    setorder(summ, depth, method); print(summ, digits = 2)
    fwrite(allv, file.path(cfg$out_dir, "synthetic_method_comparison_daily.csv"))
    fwrite(summ, file.path(cfg$out_dir, "synthetic_method_comparison_summary.csv"))

    evd <- as.Date("2022-04-12") + ev_off + c(3, 6, 9)
    cat("\nEvent days (inflow at 1 m, entrainment at 1 m, intrusion at 9 m): residual/departure vs true unknown\n")
    print(dcast(allv[(depth == 1 & day %in% evd[1:2]) | (depth == 9 & day == evd[3])],
                depth + day + unknown ~ method, value.var = "resid"), digits = 2)
  }

  out <- convert_units(out, cfg)
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  fmt <- function(d) { d <- copy(d); if ("time" %in% names(d)) d[, time := format(time, "%Y-%m-%dT%H:%M:%SZ", tz = cfg$tz)]; d }
  fwrite(fmt(out$hourly[, .(time, depth, day, DO, DO_clean, T, PAR, z_th, in_ml, h, dDOdt, gas, diffusion,
                            nep_obs, gpp_mod, r_mod, residual)]), file.path(cfg$out_dir, "budget_hourly_mmol.csv"))
  fwrite(out$daily, file.path(cfg$out_dir, "budget_daily.csv"))
  fwrite(out$areal, file.path(cfg$out_dir, "budget_productive_layer_areal.csv"))
  if (!is.null(out$kz)) fwrite(fmt(out$kz), file.path(cfg$out_dir, "Kz_interfaces_m2s.csv"))
  fwrite(fmt(out$gas), file.path(cfg$out_dir, "gas_flux_mmol_m2_h.csv"))
  fwrite(out$iw$daily, file.path(cfg$out_dir, "iw_correction_daily.csv"))
  if (!is.null(out$spec)) {                                ## NEW
    fwrite(fmt(out$spec$hourly), file.path(cfg$out_dir, "spectral_hourly.csv"))
    fwrite(out$spec$daily, file.path(cfg$out_dir, "spectral_daily.csv"))
    fwrite(out$spec$areal, file.path(cfg$out_dir, "spectral_productive_layer_areal.csv"))
  }

  cat(sprintf("\nClosure check (max |dDO/dt - sum of terms|): %.2e\n", max(abs(out$daily$closure), na.rm = TRUE)))
  print(out$daily[, .(days = .N, fitted = sum(flag == "ok"), GPP = mean(GPP), R = mean(R), NEP = mean(NEP),
                      gas = mean(gas), diffusion = mean(diffusion), residual = mean(residual)), by = depth],
        digits = 3)

  if (cfg$make_plots) {
    pl <- plot_budget(out, cfg, if (!is.null(truth) && cfg$out_units == "mmol") truth)
    if (!is.null(pl)) {
      nd <- uniqueN(out$daily$depth)
      ggplot2::ggsave(file.path(cfg$out_dir, "budget_terms_daily.png"), pl$budget, width = 9, height = 1.5 + 2.5 * nd)
      ggplot2::ggsave(file.path(cfg$out_dir, "metabolism_daily.png"), pl$rates, width = 9, height = 1.5 + 2.5 * nd)
      ggplot2::ggsave(file.path(cfg$out_dir, "budget_productive_layer.png"), pl$areal, width = 9, height = 4)
      if (!is.null(pl$methods))
        ggplot2::ggsave(file.path(cfg$out_dir, "metabolism_time_domain_vs_spectral.png"), pl$methods,
                        width = 10, height = 1.5 + 2.5 * nd)
    }
  }
  message("Outputs in ", normalizePath(cfg$out_dir))
}
