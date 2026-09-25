###############################################################################
# Time-domain internal-wave (IW) correction of dissolved oxygen
# Following Fernández Castro et al. (2021, WRR, doi:10.1029/2020WR029283)
#
#   For each 24 h segment (and depth):
#     T'  = T  - <T>          (anomaly w.r.t. daily mean)
#     DO' = DO - <DO>
#     DO' ~ a0 + a1*T' + a2*T'^2        (Eq. 4, 2nd-order polynomial)
#     DO'_IM  = fitted part explained by vertical dislocations
#     DO_clean = DO' - DO'_IM + <DO>
#
# Input : one data frame with time, DO, T (single depth; optional depth col)
# Output: DO_clean (+ daily fit diagnostics)
###############################################################################

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(repo_root)

## ---------------------------------------------------------------- CONFIG ----
cfg <- list(
  input_file     = NULL,            # path to csv; NULL -> run synthetic test
  output_file    = "DO_clean.csv",
  daily_file     = "DO_clean_daily_fits.csv",
  time_col       = "time",
  do_col         = "DO",
  temp_col       = "T",
  depth_col      = NULL,            # e.g. "depth" if several depths stacked
  time_format    = NULL,            # e.g. "%Y-%m-%d %H:%M:%S"; NULL = auto
  tz             = "UTC",           # use LOCAL SOLAR-ish time if day_start != 0
  day_start_hour = 0,               # 24 h segments start at this hour (0 = midnight)
  min_coverage   = 0.8,             # min fraction of expected samples per day
  min_T_sd       = 1e-3,            # degC; below this, no IW signal to regress on
  poly_degree    = 2,               # Eq. 4 uses 2
  diel_flag_thr  = 0.5,             # flag day if >50% of var(T') is a 24 h harmonic
  on_fail        = "NA",            # "NA" or "raw": DO_clean for days not fitted
  make_plots     = TRUE
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
  setnames(dt, c(cfg$time_col, cfg$do_col, cfg$temp_col),
               c("time", "DO", "T"))
  if (!is.null(cfg$depth_col)) setnames(dt, cfg$depth_col, "depth") else dt[, depth := NA_real_]

  if (!inherits(dt$time, "POSIXct")) {
    raw <- dt$time
    dt[, time := if (is.null(cfg$time_format)) as.POSIXct(raw, tz = cfg$tz)
                 else as.POSIXct(raw, format = cfg$time_format, tz = cfg$tz)]
    n_bad <- sum(is.na(dt$time) & !is.na(raw))
    if (n_bad) stop(n_bad, " time value(s) could not be parsed (e.g. '",
                    raw[which(is.na(dt$time) & !is.na(raw))[1]],
                    "'). Set cfg$time_format.", call. = FALSE)
  } else {
    attr(dt$time, "tzone") <- cfg$tz
  }
  for (v in c("DO", "T")) {
    if (!is.numeric(dt[[v]]))
      stop("Column '", v, "' is not numeric (class: ", class(dt[[v]])[1], ").", call. = FALSE)
  }

  setorder(dt, depth, time)
  dups <- dt[, .N, by = .(depth, time)][N > 1]
  if (nrow(dups))
    stop(nrow(dups), " duplicated timestamp(s), first at ", format(dups$time[1]),
         ". Aggregate or remove duplicates first.", call. = FALSE)
  dt
}

# Fraction of var(x) explained by a 24 h harmonic (sin/cos of time of day).
diel_fraction <- function(x, time) {
  ok <- is.finite(x)
  if (sum(ok) < 5 || var(x[ok]) == 0) return(NA_real_)
  w  <- 2 * pi * (as.numeric(time[ok]) %% 86400) / 86400
  summary(lm(x[ok] ~ sin(w) + cos(w)))$r.squared
}

correct_internal_waves <- function(dt, cfg) {
  dt <- copy(dt)

  # median sampling interval (s) per depth -> expected samples per day
  dt[, dt_s := median(diff(as.numeric(time)), na.rm = TRUE), by = depth]
  if (any(!is.finite(dt$dt_s) | dt$dt_s <= 0))
    stop("Could not determine sampling interval (need >= 2 timestamps per depth).", call. = FALSE)
  if (any(dt$dt_s > 3600))
    warning("Median sampling interval > 1 h: the IW correction needs high-frequency data.", call. = FALSE)

  # 24 h segments starting at day_start_hour
  shift_s <- cfg$day_start_hour * 3600
  dt[, day := as.Date(as.POSIXct(as.numeric(time) - shift_s,
                                 origin = "1970-01-01", tz = cfg$tz), tz = cfg$tz)]

  deg <- cfg$poly_degree
  fit_day <- function(DO, T, time, dt_s) {
    ok       <- is.finite(DO) & is.finite(T)
    n_ok     <- sum(ok)
    coverage <- n_ok / (86400 / dt_s[1])
    DO_mean  <- mean(DO[ok]); T_mean <- mean(T[ok])
    Tp  <- T  - T_mean
    DOp <- DO - DO_mean
    T_sd <- if (n_ok > 1) sd(Tp[ok]) else NA_real_
    out_na <- list(DO_mean = DO_mean, T_prime = Tp, DO_prime = DOp,
                   DO_IM = NA_real_, fitted = FALSE)

    status <- if (coverage < cfg$min_coverage) "low_coverage"
              else if (!is.finite(T_sd) || T_sd < cfg$min_T_sd) "no_T_variability"
              else if (n_ok <= deg + 2) "too_few_points"
              else "ok"

    diag <- list(n = n_ok, coverage = coverage, T_mean = T_mean, DO_mean = DO_mean,
                 T_sd = T_sd, DO_sd = if (n_ok > 1) sd(DOp[ok]) else NA_real_,
                 diel_frac_T = diel_fraction(Tp, time),
                 a0 = NA_real_, a1 = NA_real_, a2 = NA_real_,
                 r2 = NA_real_, p_value = NA_real_, status = status)
    if (status != "ok") return(list(ts = out_na, diag = diag))

    # Eq. 4: DO' = a0 + a1 T' + a2 T'^2 (raw polynomial; T' is already centred)
    X   <- sapply(seq_len(deg), function(k) Tp[ok]^k)
    fit <- lm(DOp[ok] ~ X)
    cf  <- coef(fit)
    if (anyNA(cf)) {
      diag$status <- "singular_fit"
      return(list(ts = out_na, diag = diag))
    }
    Xall <- sapply(seq_len(deg), function(k) Tp^k)
    DO_IM <- as.vector(cbind(1, Xall) %*% cf)   # NA where T is missing

    s <- summary(fit)
    fs <- s$fstatistic
    diag$a0 <- cf[1]; diag$a1 <- cf[2]; if (deg >= 2) diag$a2 <- cf[3]
    diag$r2 <- s$r.squared
    diag$p_value <- unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE))
    list(ts = list(DO_mean = DO_mean, T_prime = Tp, DO_prime = DOp,
                   DO_IM = DO_IM, fitted = TRUE),
         diag = diag)
  }

  res <- dt[, {
    r <- fit_day(DO, T, time, dt_s)
    c(list(time = time, DO = DO, T = T), r$ts,
      lapply(r$diag, function(x) rep(x, .N)))
  }, by = .(depth, day)]

  # DO_clean = DO' - DO'_IM + <DO>
  res[, DO_clean := DO_prime - DO_IM + DO_mean]
  res[, time := as.POSIXct(as.numeric(time), origin = "1970-01-01", tz = cfg$tz)]  # clean POSIXct after by-group
  if (cfg$on_fail == "raw") res[fitted == FALSE, DO_clean := DO]

  diag_cols <- c("n", "coverage", "T_mean", "DO_mean", "T_sd", "DO_sd",
                 "diel_frac_T", "a0", "a1", "a2", "r2", "p_value", "status")
  daily <- unique(res[, c("depth", "day", diag_cols), with = FALSE])
  daily[, diel_flag := !is.na(diel_frac_T) & diel_frac_T > cfg$diel_flag_thr]

  ts <- res[, .(depth, time, day, DO, T, T_prime, DO_prime, DO_IM, DO_clean)]
  if (all(is.na(ts$depth))) { ts[, depth := NULL]; daily[, depth := NULL] }

  # console summary
  cat(sprintf("Days processed: %d | fitted: %d | median R2: %.2f | diel-flagged: %d\n",
              nrow(daily), sum(daily$status == "ok"),
              median(daily$r2, na.rm = TRUE), sum(daily$diel_flag)))
  if (any(daily$status != "ok")) print(daily[status != "ok", .N, by = status])
  if (any(daily$diel_flag))
    message("Some days have T' dominated by a 24 h cycle (diel_flag = TRUE): ",
            "regressing DO' on T' there may also remove the metabolic signal.")

  list(data = ts[], daily = daily[])
}

plot_correction <- function(out) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed: skipping plots."); return(invisible(NULL))
  }
  library(ggplot2)
  d <- melt(out$data, id.vars = intersect(c("depth", "time"), names(out$data)),
            measure.vars = c("DO", "DO_clean"), variable.name = "series")
  p1 <- ggplot(d, aes(time, value, colour = series)) +
    geom_line(linewidth = 0.3) +
    scale_colour_manual(values = c(DO = "grey60", DO_clean = "firebrick")) +
    labs(x = NULL, y = "DO", colour = NULL) + theme_bw()
  if ("depth" %in% names(d)) p1 <- p1 + facet_wrap(~depth, ncol = 1, scales = "free_y")
  p2 <- ggplot(out$data, aes(T_prime, DO_prime)) +
    geom_point(size = 0.4, alpha = 0.4) +
    geom_line(aes(y = DO_IM), colour = "firebrick") +
    facet_wrap(~day, scales = "free") +
    labs(x = "T' (°C)", y = "DO'") + theme_bw()
  list(timeseries = p1, fits = p2)
}

## ---------------------------------------------------- SYNTHETIC TEST CASE ----
# Known biological diel signal + IW-induced dislocations; checks recovery.
make_synthetic <- function(days = 10, step_min = 10, seed = 1) {
  set.seed(seed)
  time <- seq(as.POSIXct("2024-07-01", tz = "UTC"),
              by = step_min * 60, length.out = days * 1440 / step_min)
  tt   <- as.numeric(time - time[1]) / 3600                  # hours
  DO_bio <- 8 + 0.4 * sin(2 * pi * (tt - 6) / 24)             # "true" metabolism
  zeta <- 1.5 * sin(2 * pi * tt / 7.3) + 0.8 * sin(2 * pi * tt / 11.1 + 1) +
          as.numeric(stats::filter(rnorm(length(tt), 0, 0.3), rep(1/6, 6), circular = TRUE))
  dTdz <- -1.2; dDOdz <- -0.9                                 # gradients at sensor
  T  <- 15 + dTdz * zeta
  DO <- DO_bio + dDOdz * zeta + 0.1 * zeta^2 + rnorm(length(tt), 0, 0.02)
  list(df = data.frame(time = time, DO = DO, T = T), DO_bio = DO_bio)
}

## ------------------------------------------------------------------ RUN ----
if (sys.nframe() == 0L) {                      # only when run as a script
  if (is.null(cfg$input_file)) {
    message("No input_file set: running synthetic test.")
    syn <- make_synthetic()
    df  <- syn$df
  } else {
    if (!file.exists(cfg$input_file)) stop("Input file not found: ", cfg$input_file, call. = FALSE)
    df <- fread(cfg$input_file)
  }

  dt  <- prepare_input(df, cfg)
  out <- correct_internal_waves(dt, cfg)

  if (is.null(cfg$input_file)) {
    rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
    cat(sprintf("RMSE vs true bio signal  raw: %.3f  |  clean: %.3f\n",
                rmse(out$data$DO, syn$DO_bio), rmse(out$data$DO_clean, syn$DO_bio)))
  }

  setwd("results/")

  fwrite(out$data,  cfg$output_fil)
  fwrite(out$daily, cfg$daily_file)

  if (cfg$make_plots) {
    pl <- plot_correction(out)
    if (!is.null(pl)) {
      ggplot2::ggsave("DO_clean_timeseries.png", pl$timeseries, width = 10, height = 4)
      ggplot2::ggsave("DO_clean_daily_fits.png", pl$fits, width = 10, height = 7)
    }
  }
}
