###############################################################################
# spectral_iw_correction.R
# Frequency-domain (spectral) internal-wave correction and diel spectral
# method of Fernández Castro et al. (2021, WRR), Supplement Text S1.
# Port of the "[1] Simple Amplitude approach" block of DO_budget.py.
#
# For each centre time (sliding window of 257 hourly values, ~10.7 d):
#   DO'  = DO - <DO>,  T' = T - <T>      (<.> = 30-d Gaussian low-pass)
#   (mixed layer) DO'' = DO' + detrend(cumsum(detrend(F_gas/zmix)))   Eq. S6
#   Step i  : S^_DO = S_DO - |H|^2 S_T,  H = S_DO,T / S_T            (IW removal)
#   Step ii : baseline S_n = a f^b fitted on 12-48 h excl. 18-28 h;
#             sigma2_DO24 = int_{18-28h} S^ - int_{18-28h} S_n
#   Step iii: sigma2_DO24 = 0 if int S^ < 1.5 int S_n
#   Step iv : sigma2_DO24 = 0 if DO-light phase < -2 h
#   A_DO24 = sqrt(2 sigma2), r_night = 2 A / t_night, R = 24 r_night
#   GPP = R + NEP (NEP from the smoothed budget; see metabolism_functions.R)
#
# Welch/CSD reproduce scipy.signal.welch/csd defaults (periodic Hann,
# 50% overlap, constant detrend per segment, 'density' scaling, one-sided).
###############################################################################

## ------------------------------------------------------------- spectra ----
hann_periodic <- function(n) 0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / n)

# One-sided Welch cross-spectral density, conj(X)*Y convention (as scipy).
# y = NULL returns the (real) power spectral density of x.
welch_csd <- function(x, y = NULL, fs = 24, nperseg = 128,
                      noverlap = nperseg %/% 2) {
  auto <- is.null(y); if (auto) y <- x
  N <- length(x)
  if (length(y) != N) stop("welch_csd: x and y differ in length", call. = FALSE)
  if (N < nperseg) stop("welch_csd: series shorter than nperseg", call. = FALSE)
  w <- hann_periodic(nperseg)
  starts <- seq(1, N - nperseg + 1, by = nperseg - noverlap)
  P <- complex(nperseg)
  for (s in starts) {
    idx <- s:(s + nperseg - 1)
    X <- fft((x[idx] - mean(x[idx])) * w)
    Y <- fft((y[idx] - mean(y[idx])) * w)
    P <- P + Conj(X) * Y
  }
  P <- P / length(starts) / (fs * sum(w^2))
  nf <- nperseg %/% 2 + 1
  P <- P[1:nf]
  dbl <- if (nperseg %% 2 == 0) 2:(nf - 1) else 2:nf
  P[dbl] <- 2 * P[dbl]
  list(f = (0:(nf - 1)) * fs / nperseg, S = if (auto) Re(P) else P)
}

trapz <- function(x, y) if (length(x) < 2) 0 else sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

detrend_lin <- function(y) {
  i <- seq_along(y)
  y - stats::predict(stats::lm(y ~ i))
}

# Band variance: "trapz" = original code (integrates only between the first
# and last bin inside the band, ~66% of a pure sinusoid's variance with
# nperseg=128 h); "rect" = sum(S)*df over the band bins (~95%).
band_var <- function(f, S, sel, method) {
  if (method == "trapz") trapz(f[sel], S[sel]) else sum(S[sel]) * (f[2] - f[1])
}

## ---------------------------------------------- one window, steps i - iv ----
# xx : DO' in the window (gaps already set to 0)          [mmol m-3]
# tt : T'  in the window (gaps set to 0)                  [degC]
# irr: incident shortwave in the window (gaps set to 0)   [W m-2]
# gas: F_gas/zmix in the window (mmol m-3 d-1) or NULL if below the mixed layer
spectral_window <- function(xx, tt, irr, gas = NULL, sp) {
  fs <- 24
  nz <- sum(xx == 0)
  zero_scale <- length(xx) / (length(xx) - nz)

  xx1 <- xx
  if (!is.null(gas)) {                              # Eq. S6 (air-lake flux)
    g <- detrend_lin(gas)
    xx1 <- xx + detrend_lin(cumsum(g / fs))
  }

  S0  <- welch_csd(xx1, fs = fs, nperseg = sp$nperseg)
  f   <- S0$f
  ST  <- welch_csd(tt,  fs = fs, nperseg = sp$nperseg)$S
  CS  <- welch_csd(xx1, tt, fs = fs, nperseg = sp$nperseg)$S
  if (sp$remove_iw) {
    H2   <- Mod(CS / ST)^2                          # |S_DO,T / S_T|^2
    ## NEW: optional coherence threshold. With 3 Welch segments, the
    ## coherence of two INDEPENDENT series is ~0.33 on average, so step i
    ## removes ~1/3 of any DO variance. coh_min > 0 only removes variance at
    ## frequencies where squared coherence exceeds coh_min.
    coh2 <- Mod(CS)^2 / (S0$S * ST)
    H2[!is.finite(coh2) | coh2 < sp$coh_min] <- 0
    Shat <- S0$S - H2 * ST                          # Step i
  } else Shat <- S0$S
  Shat[!is.finite(Shat)] <- S0$S[!is.finite(Shat)]

  f_lo <- fs / sp$period_max_h; f_hi <- fs / sp$period_min_h   # 24/28, 24/18 cpd
  band <- f >= f_lo & f <= f_hi

  # Step ii: power-law baseline on 12-48 h excluding the diel band
  fit_sel <- f >= fs / sp$baseline_max_h & f < fs / sp$baseline_min_h &
             !(f > f_lo & f <= f_hi) & Shat > 0
  if (sum(fit_sel) >= 2) {
    pp <- stats::coef(stats::lm(log10(Shat[fit_sel]) ~ log10(f[fit_sel])))
    Sn <- 10^(pp[1] + pp[2] * log10(f))
  } else Sn <- rep(0, length(f))

  var_raw   <- zero_scale * band_var(f, S0$S, band, sp$band_integration)
  var_hat   <- zero_scale * band_var(f, Shat, band, sp$band_integration)
  var_noise <- zero_scale * band_var(f, Sn,   band, sp$band_integration)
  var_do24  <- max(var_hat - var_noise, 0)
  flag <- TRUE
  if (var_hat < sp$snr_min * var_noise) { var_do24 <- 0; flag <- FALSE }   # Step iii

  # Step iv: DO-light phase (positive = DO lags light), from DO' without gas term
  CL <- welch_csd(xx, irr, fs = fs, nperseg = sp$nperseg)$S
  ib <- which(band)
  ip <- if (sp$phase_bin == "original") ib[2] else ib[which.min(abs(f[ib] - 1))]
  phase_h <- Arg(CL[ip]) * 180 / pi * 12 / 180      # degrees -> hours (24 h cycle)
  if (is.finite(phase_h) && phase_h < sp$phase_min_h) var_do24 <- 0
  if (nz > sp$max_zero_frac_flag * length(xx)) flag <- FALSE

  list(var_raw = var_raw, var_hat = var_hat, var_noise = var_noise,
       var_do24 = var_do24, phase_h = phase_h, amp_flag = flag,
       f = f, S_raw = S0$S, S_hat = Shat, S_n = Sn)
}

## ---------------------------------------------------- full time series ----
# dt: REGULAR HOURLY data.table with time, DO, T, I0, gas_term, NEP_budget,
#     sDO, sT (low-passed), t_night (h). Returns one row per centre time.
spectral_metabolism <- function(dt, sp, z_sensor, keep_spectra_at = NULL) {
  stop_if_missing(dt, c("time", "DO", "T", "I0", "gas_term", "in_ml",
                        "NEP_budget", "sDO", "sT", "t_night"), "spectral_metabolism")
  step_h <- median(diff(as.numeric(dt$time))) / 3600
  if (abs(step_h - 1) > 1e-6)
    stop("spectral_metabolism expects hourly data (got ", step_h, " h). Use regularize().", call. = FALSE)

  nt <- nrow(dt); half <- sp$seglen %/% 2
  if (nt < sp$seglen + 1) stop("Record shorter than one spectral window (", sp$seglen + 1, " h).", call. = FALSE)
  pDO <- dt$DO - dt$sDO; pDO[!is.finite(pDO)] <- 0
  pT  <- dt$T  - dt$sT;  pT[!is.finite(pT)]  <- 0
  irr <- dt$I0;          irr[!is.finite(irr)] <- 0

  centres <- seq(1, nt, by = sp$step_h)
  out <- vector("list", length(centres)); spectra <- list()
  for (k in seq_along(centres)) {
    j  <- centres[k]
    jj <- (j - half):(j + half)
    if (min(jj) < 1)  jj <- jj + (1 - min(jj))      # shift window at the edges,
    if (max(jj) > nt) jj <- jj - (max(jj) - nt)     # as in the original code
    xx <- pDO[jj]
    if (sum(xx == 0) >= sp$max_zero_frac_skip * length(xx) || !is.finite(dt$DO[j])) next
    in_ml <- isTRUE(mean(dt$in_ml[jj], na.rm = TRUE) > 0.5)
    gas <- NULL
    if (in_ml && sp$gas_correction) {
      gas <- dt$gas_term[jj] * 24                  # mmol m-3 d-1
      gas[!is.finite(gas)] <- 0
    }
    w <- spectral_window(xx, pT[jj], irr[jj], gas, sp)
    if (!is.null(keep_spectra_at) && any(abs(as.numeric(dt$time[j]) - as.numeric(keep_spectra_at)) < 1800))
      spectra[[format(dt$time[j])]] <- w
    out[[k]] <- data.table(time = dt$time[j], var_raw = w$var_raw, var_hat = w$var_hat,
                           var_noise = w$var_noise, var_do24 = w$var_do24,
                           phase_h = w$phase_h, amp_flag = w$amp_flag)
  }
  res <- rbindlist(out)
  res <- merge(res, dt[, .(time, t_night, NEP = NEP_budget)], by = "time")
  res[, A24 := sqrt(2 * var_do24)]                   # zero-to-crest amplitude
  res[, R := 24 / t_night * 2 * A24]                  # Eqs. S1-S2, mmol m-3 d-1
  res[NEP < 0, R := R + abs(NEP)]                     # as in DO_budget.py
  res[, GPP := R + NEP]
  attr(res, "spectra") <- spectra
  res[]
}

# Default spectral settings (= DO_budget.py)
spectral_defaults <- list(
  nperseg = 128, seglen = 256, step_h = 1,
  period_min_h = 18, period_max_h = 28,        # diel band
  baseline_min_h = 12, baseline_max_h = 48,    # baseline fit range
  snr_min = 1.5, phase_min_h = -2,
  max_zero_frac_skip = 0.5, max_zero_frac_flag = 0.1,
  remove_iw = TRUE, gas_correction = TRUE,
  coh_min = 0,                                 # NEW: 0 = original; e.g. 0.7 = only significant coherence
  band_integration = "trapz",                  # "trapz" (original) | "rect"
  phase_bin = "original"                       # "original" (2nd band bin) | "nearest" (to 1 cpd)
)
