###############################################################################
# Frequency-domain (spectral) internal-wave correction of dissolved oxygen
# Time-series version of step (i) of Fernández Castro et al. (2021, WRR,
# Supplement Text S1): S^_DO = S_DO - |H|^2 S_T,  H = S_DO,T / S_T
#
# The paper applies step (i) to spectra only. Here the same transfer function
# is applied to the T' time series, so that a DO_clean series is returned:
#   DO'    = DO - <DO>,   T' = T - <T>        (<.> = 30-d Gaussian low-pass)
#   H(f)   = <conj(T'_f) DO'_f> / <|T'_f|^2>  (Welch, 257-h window, 128-h segments)
#   DO'_IM = IFFT[ H(f) * FFT(T') ]           (only the central block is kept)
#   DO_clean = DO - DO'_IM
# With one H per window (block_h ~ window_h) the spectrum of DO_clean matches
# the paper's S^_DO (diel-band variance ratio 1.03-1.09 on synthetic data).
# With the default 24-h blocks H adapts more locally and removes somewhat more
# (ratio ~0.7-0.8), which was closer to the truth in the synthetic tests.
#
# Differences with the time-domain correction (iw_correction_time_domain.R):
#   + H is complex and frequency-dependent: phase lags between DO and T and
#     different DO/T ratios for different wave periods are allowed
#   - linear in T' (no T'^2 term), needs >= ~11 d of data around each block
#
# Input : csv with time, DO, T (+ optional depth column for stacked sensors)
# Output: csv with DO_IM and DO_clean (+ a per-block diagnostics csv)
###############################################################################

## ---------------------------------------------------------------- CONFIG ----
cfg <- list(
  input_file    = NULL,              # path to csv; NULL -> run synthetic test
  output_file   = "DO_clean_spectral.csv",
  diag_file     = "DO_clean_spectral_blocks.csv",
  time_col      = "time",
  do_col        = "DO",
  temp_col      = "T",
  depth_col     = NULL,              # e.g. "depth" if several depths are stacked
  time_format   = NULL,              # e.g. "%Y-%m-%d %H:%M:%S"; NULL = auto
  tz            = "UTC",
  lowpass_days  = 30,                # half-width of the Gaussian low-pass (as DO_budget.py)
  window_h      = 256,               # analysis window (h); 257 hourly points in the paper
  nperseg_h     = 128,               # Welch sub-segment (h), 50% overlap, Hann
  block_h       = 24,                # DO_clean is built block by block (centre of each window)
  coh_min       = 0,                 # 0 = paper; e.g. 0.75 = only remove significantly coherent variance
  protect_diel  = FALSE,             # TRUE: do not remove T-coherent variance in the diel band
  diel_band_h   = c(18, 28),
  max_gap_frac  = 0.5,               # skip block if > this fraction of the window is missing
  on_fail       = "NA",              # "NA" or "raw": DO_clean for skipped blocks
  diel_coh_warn = 0.5,               # warn if mean coherence^2 in the diel band exceeds this
  make_plots    = TRUE
)

suppressPackageStartupMessages(library(data.table))

## ------------------------------------------------------------- FUNCTIONS ----

prepare_input <- function(df, cfg) {
  dt <- as.data.table(df)
  needed <- c(cfg$time_col, cfg$do_col, cfg$temp_col, cfg$depth_col)
  miss <- setdiff(needed, names(dt))
  if (length(miss))
    stop("Missing column(s) in input: ", paste(miss, collapse = ", "),
         "\n  Available: ", paste(names(dt), collapse = ", "), call. = FALSE)
  dt <- dt[, ..needed]
  setnames(dt, c(cfg$time_col, cfg$do_col, cfg$temp_col), c("time", "DO", "T"))
  if (!is.null(cfg$depth_col)) setnames(dt, cfg$depth_col, "depth") else dt[, depth := NA_real_]

  if (is.numeric(dt$time)) {
    dt[, time := as.POSIXct(time, origin = "1970-01-01", tz = cfg$tz)]
  } else if (!inherits(dt$time, "POSIXct")) {
    raw <- dt$time
    dt[, time := if (is.null(cfg$time_format)) as.POSIXct(raw, tz = cfg$tz,
                   tryFormats = c("%Y-%m-%dT%H:%M:%OSZ", "%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M", "%Y-%m-%d"))
                 else as.POSIXct(raw, format = cfg$time_format, tz = cfg$tz)]
    bad <- which(is.na(dt$time) & !is.na(raw))
    if (length(bad)) stop(length(bad), " time value(s) could not be parsed (e.g. '", raw[bad[1]],
                          "'). Set cfg$time_format.", call. = FALSE)
  }
  for (v in c("DO", "T"))
    if (!is.numeric(dt[[v]])) stop("Column '", v, "' is not numeric.", call. = FALSE)
  setorder(dt, depth, time)
  dt
}

# Regular grid per depth (median sampling step); duplicates within a bin are averaged
regular_grid <- function(d) {
  step <- median(diff(as.numeric(d$time)))
  if (!is.finite(step) || step <= 0) stop("Cannot determine the sampling interval.", call. = FALSE)
  if (step > 3600) stop("Sampling interval > 1 h: the spectral correction needs (sub)hourly data.", call. = FALSE)
  t0 <- min(as.numeric(d$time))
  d[, k := round((as.numeric(time) - t0) / step)]
  if (anyDuplicated(d$k)) warning("Several samples per time bin: averaged.", call. = FALSE)
  g <- d[, .(DO = mean(DO, na.rm = TRUE), T = mean(T, na.rm = TRUE), n_obs = .N), by = k]
  full <- data.table(k = 0:max(g$k))
  g <- g[full, on = "k"]
  g[, time := as.POSIXct(t0 + k * step, origin = "1970-01-01", tz = attr(d$time, "tzone"))]
  g[is.nan(DO), DO := NA]; g[is.nan(T), T := NA]
  list(grid = g[], step_s = step)
}

# Gaussian low-pass as in DO_budget.py: w = exp(-dt^2/(0.5 hw)^2), |dt| <= hw,
# normalised convolution so that gaps are handled
lowpass <- function(x, step_days, hw_days) {
  k <- seq(-floor(hw_days / step_days), floor(hw_days / step_days))
  w <- exp(-((k * step_days)^2) / (0.5 * hw_days)^2)
  ok <- is.finite(x); pad <- length(k) %/% 2
  xp <- c(rep(0, pad), ifelse(ok, x, 0), rep(0, pad))
  op <- c(rep(0, pad), as.numeric(ok), rep(0, pad))
  num <- stats::filter(xp, w, sides = 2); den <- stats::filter(op, w, sides = 2)
  r <- as.numeric(num / den)[(pad + 1):(pad + length(x))]
  r[!is.finite(r)] <- NA
  r
}

hann_periodic <- function(n) 0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / n)

# One-sided Welch cross-spectrum, conj(X)*Y (= scipy.signal.csd defaults)
welch_csd <- function(x, y = NULL, fs, nperseg, noverlap = nperseg %/% 2) {
  auto <- is.null(y); if (auto) y <- x
  w <- hann_periodic(nperseg)
  starts <- seq(1, length(x) - nperseg + 1, by = nperseg - noverlap)
  P <- complex(nperseg)
  for (s in starts) {
    i <- s:(s + nperseg - 1)
    P <- P + Conj(fft((x[i] - mean(x[i])) * w)) * fft((y[i] - mean(y[i])) * w)
  }
  P <- P / length(starts) / (fs * sum(w^2))
  nf <- nperseg %/% 2 + 1; P <- P[1:nf]
  dbl <- if (nperseg %% 2 == 0) 2:(nf - 1) else 2:nf
  P[dbl] <- 2 * P[dbl]
  list(f = (0:(nf - 1)) * fs / nperseg, S = if (auto) Re(P) else P)
}

# Apply a transfer function H (given at Welch frequencies fw) to x by FFT,
# zero-padded to 2N to avoid circular wrap-around. Returns the filtered series.
apply_transfer <- function(x, fw, H, fs) {
  N <- length(x); M <- 2 * N
  xp <- c(x - mean(x), rep(0, N))
  fg <- (0:(M - 1)) * fs / M
  neg <- fg > fs / 2
  fpos <- ifelse(neg, fs - fg, fg)
  Hr <- approx(fw, Re(H), fpos, rule = 2)$y
  Hi <- approx(fw, Im(H), fpos, rule = 2)$y
  Hf <- complex(real = Hr, imaginary = ifelse(neg, -Hi, Hi))   # Hermitian -> real output
  Re(fft(fft(xp) * Hf, inverse = TRUE) / M)[1:N]
}

correct_one_depth <- function(d, cfg) {
  rg <- regular_grid(copy(d)); g <- rg$grid; step_s <- rg$step_s
  sph <- 3600 / step_s                              # samples per hour
  fs  <- 86400 / step_s                             # samples per day -> f in cpd
  nw  <- round(cfg$window_h * sph) + 1
  nps <- round(cfg$nperseg_h * sph)
  nb  <- max(1, round(cfg$block_h * sph))
  nt  <- nrow(g)
  if (nt < nw) stop("Record (", round(nt / sph / 24, 1), " d) shorter than one analysis window (",
                    round(cfg$window_h / 24, 1), " d).", call. = FALSE)

  step_days <- step_s / 86400
  g[, DO_lp := lowpass(DO, step_days, cfg$lowpass_days)]
  g[, T_lp  := lowpass(T,  step_days, cfg$lowpass_days)]
  pDO <- g$DO - g$DO_lp; pT <- g$T - g$T_lp
  miss <- !is.finite(pDO) | !is.finite(pT)
  pDO[miss] <- 0; pT[miss] <- 0                    # gaps -> 0, as in DO_budget.py

  DO_IM <- rep(NA_real_, nt)
  f_lo <- 24 / cfg$diel_band_h[2]; f_hi <- 24 / cfg$diel_band_h[1]
  starts <- seq(1, nt, by = nb)
  diag <- vector("list", length(starts))
  for (b in seq_along(starts)) {
    blk <- starts[b]:min(starts[b] + nb - 1, nt)
    ctr <- round(mean(range(blk)))
    jj  <- (ctr - (nw - 1) %/% 2):(ctr - (nw - 1) %/% 2 + nw - 1)
    if (min(jj) < 1)  jj <- jj + (1 - min(jj))
    if (max(jj) > nt) jj <- jj - (max(jj) - nt)
    gap <- mean(miss[jj])
    status <- if (gap > cfg$max_gap_frac) "too_many_gaps" else "ok"
    coh_iw <- coh_diel <- frac_rm <- NA_real_
    if (status == "ok") {
      x <- pDO[jj]; y <- pT[jj]
      SD <- welch_csd(x, fs = fs, nperseg = nps)
      ST <- welch_csd(y, fs = fs, nperseg = nps)$S
      C  <- welch_csd(y, x, fs = fs, nperseg = nps)$S          # conj(T) * DO
      H  <- C / ST
      coh2 <- Mod(C)^2 / (ST * SD$S)
      H[!is.finite(H) | !is.finite(coh2) | coh2 < cfg$coh_min] <- 0
      H[1] <- 0                                                 # keep the mean
      diel <- SD$f >= f_lo & SD$f <= f_hi
      if (cfg$protect_diel) H[diel] <- 0
      im <- apply_transfer(y, SD$f, H, fs)
      pos <- match(blk, jj)
      DO_IM[blk] <- im[pos]
      coh_diel <- mean(coh2[diel], na.rm = TRUE)
      coh_iw   <- mean(coh2[SD$f > f_hi & SD$f <= fs / 2], na.rm = TRUE)
      frac_rm  <- var(im) / var(x)
    }
    diag[[b]] <- data.table(block_start = g$time[blk[1]], block_end = g$time[max(blk)],
                            gap_frac = gap, coh2_diel = coh_diel, coh2_subdiel = coh_iw,
                            frac_var_removed = frac_rm, status = status)
  }
  g[, DO_IM := DO_IM]
  g[miss & !is.finite(DO_IM), DO_IM := NA]
  g[!is.finite(T), DO_IM := NA]                     # no T -> no correction possible
  g[, DO_clean := DO - DO_IM]
  if (cfg$on_fail == "raw") g[is.na(DO_IM) & is.finite(DO), DO_clean := DO]
  list(data = g[n_obs > 0 & !is.na(n_obs), .(time, DO, T, DO_lp, T_lp, DO_IM, DO_clean)],
       blocks = rbindlist(diag))
}

correct_internal_waves_spectral <- function(dt, cfg) {
  if (cfg$window_h < 2 * cfg$nperseg_h)
    warning("window_h < 2 * nperseg_h: fewer than 3 Welch segments, H will be very noisy.", call. = FALSE)
  out <- lapply(split(dt, by = "depth", keep.by = TRUE, sorted = TRUE), function(d) {
    r <- correct_one_depth(d[, .(time, DO, T)], cfg)
    r$data[, depth := d$depth[1]]; r$blocks[, depth := d$depth[1]]
    r
  })
  data   <- rbindlist(lapply(out, `[[`, "data"))
  blocks <- rbindlist(lapply(out, `[[`, "blocks"))
  setcolorder(data, "depth"); setcolorder(blocks, "depth")
  if (all(is.na(data$depth))) { data[, depth := NULL]; blocks[, depth := NULL] }

  cat(sprintf("Blocks: %d | ok: %d | median fraction of DO' variance removed: %.2f\n",
              nrow(blocks), sum(blocks$status == "ok"),
              median(blocks$frac_var_removed, na.rm = TRUE)))
  if (!cfg$protect_diel && any(blocks$coh2_diel > cfg$diel_coh_warn, na.rm = TRUE))
    message(sum(blocks$coh2_diel > cfg$diel_coh_warn, na.rm = TRUE),
            " block(s) with diel-band coherence^2 > ", cfg$diel_coh_warn,
            ": T' and DO' covary at ~24 h. Fine if it is internal waves at diel periods;",
            " if it is diel heating (mixed layer), the correction also removes metabolic",
            " signal -> consider protect_diel = TRUE.")
  list(data = data[], blocks = blocks[])
}

plot_correction <- function(out) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) { message("ggplot2 not installed: no plots."); return(NULL) }
  d <- melt(out$data, id.vars = intersect(c("depth", "time"), names(out$data)),
            measure.vars = c("DO", "DO_clean"), variable.name = "series")
  p <- ggplot2::ggplot(d, ggplot2::aes(time, value, colour = series)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = c(DO = "grey60", DO_clean = "steelblue4")) +
    ggplot2::labs(x = NULL, y = "DO", colour = NULL) + ggplot2::theme_bw()
  if ("depth" %in% names(d)) p <- p + ggplot2::facet_wrap(~depth, ncol = 1, scales = "free_y")
  p
}

## ---------------------------------------------------- SYNTHETIC TEST CASE ----
# Known biological diel cycle + internal waves whose DO signature LAGS T by
# 45 min at the 8-h mode (e.g. sensor response / horizontal offset), which a
# pointwise DO'~T' regression cannot capture but a complex H(f) can.
make_synthetic <- function(days = 40, step_min = 10, seed = 3) {
  set.seed(seed)
  time <- seq(as.POSIXct("2024-07-01", tz = "UTC"), by = step_min * 60, length.out = days * 1440 / step_min)
  th <- as.numeric(time - time[1]) / 3600
  DO_bio <- 300 + 0.05 * th / 24 + 4 * sin(2 * pi * (th - 9) / 24)
  e <- rnorm(length(th)); phi <- exp(-(step_min / 60) / 10)
  red <- as.numeric(stats::filter(e * sqrt(1 - phi^2), phi, method = "recursive"))
  z_now <- 1.2 * red + 1.0 * sin(2 * pi * th / 8) + 0.6 * sin(2 * pi * th / 13.5 + 1)
  z_lag <- 1.2 * red + 1.0 * sin(2 * pi * (th - 0.75) / 8) + 0.6 * sin(2 * pi * th / 13.5 + 1)
  T  <- 14 - 1.0 * z_now + rnorm(length(th), 0, 0.01)
  DO <- DO_bio - 6 * z_lag + rnorm(length(th), 0, 0.5)
  list(df = data.frame(time = time, DO = DO, T = T), DO_bio = DO_bio)
}

## ------------------------------------------------------------------ RUN ----
## Runs with Rscript and with RStudio's "Source" button; skipped when another
## script sets options(iw.library_mode = TRUE) before sourcing this file.
if (!isTRUE(getOption("iw.library_mode", FALSE))) {
  if (is.null(cfg$input_file)) {
    message("No input_file set: running synthetic test.")
    syn <- make_synthetic(); df <- syn$df
  } else {
    if (!file.exists(cfg$input_file)) stop("Input file not found: ", cfg$input_file, call. = FALSE)
    df <- fread(cfg$input_file)
  }
  dt  <- prepare_input(df, cfg)
  out <- correct_internal_waves_spectral(dt, cfg)

  if (is.null(cfg$input_file)) {
    hp <- function(x) x - lowpass(x, 10 / 1440, cfg$lowpass_days)
    rmse <- function(a) sqrt(mean((hp(a) - hp(syn$DO_bio))^2, na.rm = TRUE))
    cat(sprintf("RMSE vs true biological anomaly  raw: %.2f  |  DO_clean: %.2f\n",
                rmse(out$data$DO), rmse(out$data$DO_clean)))
  }
  fwrite(out$data[, time := format(time, "%Y-%m-%dT%H:%M:%SZ", tz = cfg$tz)], cfg$output_file)
  fwrite(out$blocks[, `:=`(block_start = format(block_start, "%Y-%m-%dT%H:%M:%SZ", tz = cfg$tz),
                           block_end   = format(block_end,   "%Y-%m-%dT%H:%M:%SZ", tz = cfg$tz))],
         cfg$diag_file)
  if (cfg$make_plots) {
    out$data[, time := as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%SZ", tz = cfg$tz)]
    p <- plot_correction(out)
    if (!is.null(p)) ggplot2::ggsave("DO_clean_spectral.png", p, width = 10, height = 4)
  }
}
