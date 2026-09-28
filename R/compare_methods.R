###############################################################################
# compare_methods.R
# Time-domain vs spectral (frequency-domain) diel methods for NEP, GPP, R,
# following Fernández Castro et al. (2021), on synthetic data with known truth.
#
# Needs in cfg$code_dir: iw_correction_time_domain.R, metabolism_functions.R,
#                        spectral_iw_correction.R, synthetic_data.R
# Replace the synthetic block by your own data (see "YOUR DATA" below).
###############################################################################

## ---------------------------------------------------------------- CONFIG ----
cfg <- list(
  code_dir       = ".",
  out_dir        = ".",
  cases          = c("metalimnion", "metalimnion_strongIW", "surface"),
  days           = 80,
  lat            = 46.5,    elev = 372,
  kd             = 0.3,             # m-1 (light attenuation, for Iz in Eq. S7)
  theta          = 1.07,            # Eq. S7 temperature coefficient
  tz             = "UTC",
  day_start_hour = 0,
  min_coverage   = 0.8,
  r2_IW_max      = 0.75,            # Text S2: day rejected if R2(DO'~T') > 0.75
  zmix_default   = NA_real_,
  smooth_hw_days = 30,              # <.> low-pass, as in DO_budget.py
  trim_days      = 6,               # ignore record edges in the metrics
  example_time   = "2019-07-02 00:00:00",
  coh_min        = 0.75,            # for the "coh-threshold" spectral variant (~95% null level, 3 segments)
  seed           = 42
)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
src <- function(f) {
  p <- file.path(cfg$code_dir, f)
  if (!file.exists(p)) stop("Cannot find ", p, " (set cfg$code_dir).", call. = FALSE)
  p
}
source(src("/R/metabolism_functions.R"))
source(src("/R/spectral_iw_correction.R"))
source(src("/R/synthetic_data.R"))
options(iw.library_mode = TRUE)   # load the IW scripts without running their RUN block
td_env <- new.env(); sys.source(src("/R/iw_correction_time_domain.R"), envir = td_env)

## ----------------------------------------------------- helper: diel var ----
# Diel-band (18-28 h) variance of any series' anomaly, same windows as the
# spectral method -> compares IW removal of both approaches on equal footing.
diel_band_var_series <- function(x, time, centres, sp) {
  px <- x - gauss_smooth(x, time, cfg$smooth_hw_days); px[!is.finite(px)] <- 0
  nt <- length(x); half <- sp$seglen %/% 2
  vapply(centres, function(j) {
    jj <- (j - half):(j + half)
    if (min(jj) < 1) jj <- jj + (1 - min(jj)); if (max(jj) > nt) jj <- jj - (max(jj) - nt)
    s <- welch_csd(px[jj], fs = 24, nperseg = sp$nperseg)
    band <- s$f >= 24 / sp$period_max_h & s$f <= 24 / sp$period_min_h
    band_var(s$f, s$S, band, sp$band_integration)
  }, numeric(1))
}

metrics <- function(est, obs) {
  ok <- is.finite(est) & is.finite(obs)
  data.table(n = sum(ok), bias = mean(est[ok] - obs[ok]),
             rmse = sqrt(mean((est[ok] - obs[ok])^2)),
             r = if (sum(ok) > 2) cor(est[ok], obs[ok]) else NA_real_)
}

## ------------------------------------------------------------ run cases ----
all_daily <- list(); all_var <- list(); all_ts <- list(); all_spec <- list()

for (case in cfg$cases) {
  message("\n=== Case: ", case, " ===")
  ## ---- SYNTHETIC DATA (replace this block with YOUR DATA) -------------
  syn <- make_synthetic_lake(case, days = cfg$days, lat = cfg$lat, elev = cfg$elev,
                             kd = cfg$kd, seed = cfg$seed)
  fwrite(syn[, time := format(time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")],
         file.path(cfg$out_dir, sprintf("synthetic_%s.csv", case)))
  syn[, time := as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")]
  truth <- truth_daily(syn)
  z_sensor <- syn$z_sensor[1]
  # YOUR DATA: a data.table with time (POSIXct), DO (mmol m-3), T (degC),
  #   I0 (W m-2), and for mixed-layer sensors zmix (m) + k_gas (m h-1) or wind
  #   (m s-1); optional 'transport' (mmol m-3 h-1). Then: dat <- regularize(x)
  dat <- regularize(syn[, .(time, DO, T, I0, wind, zmix)], 60, cfg$tz)
  ## -------------------------------------------------------------------

  ccfg <- c(cfg, list(z_sensor = z_sensor))
  dat <- add_physics_terms(dat, ccfg)
  dat[, t_night := 24 - day_length_h(time, cfg$lat)]
  bud <- nep_budget_smooth(dat, cfg$smooth_hw_days)
  dat[, `:=`(sDO = bud$sDO, NEP_budget = bud$NEP,
             sT = gauss_smooth(T, time, cfg$smooth_hw_days))]

  # ---- time-domain method
  td <- td_metabolism(dat, ccfg, td_env)

  # ---- spectral method (original integration, and rectangle integration)
  ex_t <- as.POSIXct(cfg$example_time, tz = cfg$tz)
  sp_tr <- spectral_defaults
  sp_re <- modifyList(spectral_defaults, list(band_integration = "rect"))
  spec_tr <- spectral_metabolism(dat, sp_tr, z_sensor, keep_spectra_at = ex_t)
  spec_re <- spectral_metabolism(dat, sp_re, z_sensor)
  sp_co <- modifyList(sp_re, list(coh_min = cfg$coh_min))
  spec_co <- spectral_metabolism(dat, sp_co, z_sensor)
  all_spec[[case]] <- attr(spec_tr, "spectra")[[1]]
  day_of <- function(x) as.Date(as.POSIXct(as.numeric(x$time) - cfg$day_start_hour * 3600,
                                           origin = "1970-01-01", tz = cfg$tz), tz = cfg$tz)
  agg <- function(s, lab) s[, .(GPP = mean(GPP), R = mean(R), NEP = mean(NEP)),
                            by = .(day = day_of(s))][, method := lab]

  daily <- rbindlist(list(
    td$daily[, .(day, GPP, R, NEP, method = "time domain")],
    agg(spec_tr, "spectral (trapz, original)"),
    agg(spec_re, "spectral (rect)"),
    agg(spec_co, "spectral (rect, coh-threshold)")), use.names = TRUE)
  daily <- merge(daily, truth, by = "day")
  daily[, case := case]
  all_daily[[case]] <- daily

  # ---- diel-band variance: raw, spectral S^, time-domain DO_clean, truth
  noon <- which(format(dat$time, "%H") == "12")
  sp_cmp <- sp_re
  x <- merge(dat[, .(time)], td$hourly[, .(time, DO_clean)], by = "time", all.x = TRUE)
  var_dt <- data.table(
    time      = dat$time[noon],
    raw       = diel_band_var_series(dat$DO, dat$time, noon, sp_cmp),
    time_dom  = diel_band_var_series(x$DO_clean, x$time, noon, sp_cmp),
    truth     = diel_band_var_series(syn$DO_bio, syn$time, noon, sp_cmp))
  var_dt <- merge(var_dt, spec_re[, .(time, spectral_hat = var_hat,
                                      spectral_final = var_do24)], by = "time")
  var_dt <- merge(var_dt, spec_co[, .(time, spectral_hat_coh = var_hat)], by = "time")
  var_dt[, case := case]
  all_var[[case]] <- var_dt

  all_ts[[case]] <- merge(syn[, .(time, DO_obs = DO, DO_bio)],
                          td$hourly[, .(time, DO_clean)], by = "time")[, case := case]

  message(sprintf("TD days flagged: %s", paste(capture.output(print(td$daily[, .N, by = flag])), collapse = "\n")))
}

daily <- rbindlist(all_daily)
vars  <- rbindlist(all_var)

## ------------------------------------------------------------- metrics ----
d0 <- min(daily$day) + cfg$trim_days; d1 <- max(daily$day) - cfg$trim_days
dd <- daily[day >= d0 & day <= d1]
long <- melt(dd, id.vars = c("case", "method", "day"),
             measure.vars = list(est = c("GPP", "R", "NEP"),
                                 obs = c("GPP_true", "R_true", "NEP_true")),
             variable.name = "var")
long[, var := c("GPP", "R", "NEP")[var]]
m_daily <- long[, metrics(est, obs), by = .(case, method, var)]
long[, week := as.integer(as.numeric(day - d0) %/% 7)]
m_week <- long[, .(est = mean(est, na.rm = TRUE), obs = mean(obs)), by = .(case, method, var, week)][
  , metrics(est, obs), by = .(case, method, var)]

cat("\n==== Daily estimates vs truth (mmol m-3 d-1), interior days ====\n")
print(m_daily[order(case, var, method)], digits = 3)
cat("\n==== Weekly means vs truth ====\n")
print(m_week[order(case, var, method)], digits = 3)

v_sum <- vars[time >= as.POSIXct(d0) & time <= as.POSIXct(d1),
              lapply(.SD, median, na.rm = TRUE), by = case,
              .SDcols = c("raw", "spectral_hat", "spectral_hat_coh", "spectral_final", "time_dom", "truth")]
cat("\n==== Median diel-band (18-28 h) DO variance, (mmol m-3)^2 ====\n")
print(v_sum, digits = 3)

fwrite(daily, file.path(cfg$out_dir, "compare_daily_estimates.csv"))
vars[, time := format(time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")]
fwrite(rbind(m_daily[, scale := "daily"], m_week[, scale := "weekly"]),
       file.path(cfg$out_dir, "compare_metrics.csv"))
fwrite(vars, file.path(cfg$out_dir, "compare_diel_variance.csv"))

## --------------------------------------------------------------- plots ----
cols <- c("truth" = "black", "time domain" = "#1b9e77",
          "spectral (trapz, original)" = "#d95f02", "spectral (rect)" = "#7570b3",
          "spectral (rect, coh-threshold)" = "#e7298a")

p_rates <- {
  tr <- unique(daily[, .(case, day, GPP = GPP_true, R = R_true, NEP = NEP_true)])[, method := "truth"]
  pl <- melt(rbind(daily[, .(case, day, GPP, R, NEP, method)], tr),
             id.vars = c("case", "day", "method"))
  ggplot(pl, aes(day, value, colour = method)) +
    geom_line(linewidth = 0.4) + geom_hline(yintercept = 0, linewidth = 0.2) +
    facet_grid(variable ~ case, scales = "free_y") +
    scale_colour_manual(values = cols) +
    labs(x = NULL, y = expression(mmol~m^-3~d^-1), colour = NULL) +
    theme_bw() + theme(legend.position = "bottom")
}
ggsave(file.path(cfg$out_dir, "compare_rates.png"), p_rates, width = 10, height = 7, dpi = 150)

p_ts <- {
  ts <- rbindlist(all_ts)[time >= as.POSIXct("2019-06-25", tz = "UTC") &
                            time <  as.POSIXct("2019-07-03", tz = "UTC")]
  pl <- melt(ts, id.vars = c("case", "time"))
  pl[, variable := factor(variable, c("DO_obs", "DO_clean", "DO_bio"),
                          c("observed", "time-domain DO_clean", "true biological DO"))]
  ggplot(pl, aes(time, value, colour = variable)) + geom_line(linewidth = 0.4) +
    facet_wrap(~case, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("grey60", "#1b9e77", "black")) +
    labs(x = NULL, y = expression(DO~(mmol~m^-3)), colour = NULL) +
    theme_bw() + theme(legend.position = "bottom")
}
ggsave(file.path(cfg$out_dir, "compare_timeseries.png"), p_ts, width = 10, height = 6, dpi = 150)

p_spec <- wrap_plots(lapply(names(all_spec), function(cs) {
  w <- all_spec[[cs]]
  pl <- data.table(f = w$f, S_DO = w$S_raw, `S^_DO` = w$S_hat, S_n = w$S_n)[f > 0]
  pl <- melt(pl, id.vars = "f")
  ggplot(pl[value > 0], aes(f, value, colour = variable, linewidth = variable)) +
    annotate("rect", xmin = 24 / 28, xmax = 24 / 18, ymin = 0, ymax = Inf, alpha = 0.15) +
    geom_line() + scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("grey60", "black", "grey40")) +
    scale_linewidth_manual(values = c(0.4, 0.7, 1.4)) +
    labs(title = sprintf("%s, window centred %s", cs, cfg$example_time),
         subtitle = sprintf("var S^ = %.1f, var Sn = %.1f, var DO24 = %.1f, phase = %.1f h",
                            w$var_hat, w$var_noise, w$var_do24, w$phase_h),
         x = "f (cpd)", y = expression(S[DO]~((mmol~m^-3)^2/cpd)), colour = NULL, linewidth = NULL) +
    theme_bw()
}), ncol = 1)
ggsave(file.path(cfg$out_dir, "compare_example_spectra.png"), p_spec, width = 7, height = 7, dpi = 150)

p_var <- {
  pl <- melt(vars, id.vars = c("case", "time"))
  ggplot(pl, aes(time, value, colour = variable)) + geom_line() +
    facet_wrap(~case, ncol = 1, scales = "free_y") + scale_y_log10() +
    labs(x = NULL, y = expression(diel~band~variance~((mmol~m^-3)^2)), colour = NULL) +
    theme_bw() + theme(legend.position = "bottom")
}
ggsave(file.path(cfg$out_dir, "compare_diel_variance.png"), p_var, width = 10, height = 6, dpi = 150)
message("\nDone. Outputs written to ", normalizePath(cfg$out_dir))
