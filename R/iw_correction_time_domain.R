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
# Input : one data frame with time, DO, T (+ depth col for stacked sensors)
# Output: DO_clean (+ daily fit diagnostics)
#
# NEW: near-surface sensors (depth <= surface_depth_max) follow surface_rule:
#   "always": corrected like the others (original behaviour)
#   "skip"  : never corrected (DO_clean = DO)
#   "hybrid": corrected on a given day only if the reference sensor below
#             (default: shallowest sensor deeper than surface_depth_max)
#             shows internal waves (DO'~T' R2 >= iw_r2_min and T' not
#             dominated by a 24 h cycle) AND the near-surface T' moves with
#             the reference T' (correlation of sub-diel T' >= link_cor_min).
#             Otherwise DO_clean = DO for that day.
###############################################################################

## ---------------------------------------------------------------- CONFIG ----
cfg <- list(
  input_file     = "C:/Projects/myGit/o2-internalwaves/data/Gergal_example_Camille_formatted.csv",            # path to csv; NULL -> run synthetic test
  output_file    = "DO_clean.csv",
  daily_file     = "DO_clean_daily_fits.csv",
  time_col       = "time",
  do_col         = "DO",
  temp_col       = "T",
  depth_col      = "depth",            # e.g. "depth" if several depths stacked
  time_format    = NULL,            # e.g. "%Y-%m-%d %H:%M:%S"; NULL = auto
  tz             = "UTC",           # use LOCAL SOLAR-ish time if day_start != 0
  day_start_hour = 0,               # 24 h segments start at this hour (0 = midnight)
  min_coverage   = 0.8,             # min fraction of expected samples per day
  min_T_sd       = 1e-3,            # degC; below this, no IW signal to regress on
  poly_degree    = 2,               # Eq. 4 uses 2
  diel_flag_thr  = 0.5,             # flag day if >50% of var(T') is a 24 h harmonic
  on_fail        = "NA",            # "NA" or "raw": DO_clean for days not fitted
  make_plots     = TRUE,
  out_dir        = "results",       ## NEW: outputs written here (no setwd)
  ## NEW: near-surface sensors -------------------------------------------
  surface_rule      = "hybrid",     # "always" | "skip" | "hybrid"
  surface_depth_max = 3,            # m; sensors at depth <= this are "near-surface"
  ref_depth         = NULL,         # reference sensor for "hybrid"; NULL = shallowest deeper sensor
  iw_r2_min         = 0.3,          # hybrid: IW present at reference if daily R2(DO'~T') >= this
  link_cor_min      = 0.5,          # hybrid: min correlation of sub-diel T' (surface vs reference)
  link_remove_diel  = TRUE          # hybrid: remove 24 h harmonic from T' before correlating
)

suppressPackageStartupMessages(library(data.table))

## FIX: the RUN block used `sys.nframe() == 0`, which is FALSE when the file is
## run with RStudio's "Source" button. Now it always runs, unless another
## script sets options(iw.library_mode = TRUE) before sourcing this file.
run_main <- !isTRUE(getOption("iw.library_mode", FALSE))

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

## NEW: correlation of sub-diel T' between two sensors for one day ---------
remove_diel <- function(x, time) {
  ok <- is.finite(x)
  if (sum(ok) < 5) return(x)
  w <- 2 * pi * (as.numeric(time) %% 86400) / 86400
  x[ok] <- residuals(lm(x[ok] ~ sin(w[ok]) + cos(w[ok])))
  x
}

link_cor_day <- function(a, b, bucket_s, remove_24h) {
  # a, b: data.tables with time and T_prime for one day at two depths
  agg <- function(d) d[, .(Tp = mean(T_prime, na.rm = TRUE)),
                       by = .(tb = round(as.numeric(time) / bucket_s))]
  m <- merge(agg(a), agg(b), by = "tb", suffixes = c("_s", "_r"))
  m <- m[is.finite(Tp_s) & is.finite(Tp_r)]
  if (nrow(m) < 10) return(NA_real_)
  if (remove_24h) {
    tt <- m$tb * bucket_s
    m[, `:=`(Tp_s = remove_diel(Tp_s, tt), Tp_r = remove_diel(Tp_r, tt))]
  }
  if (sd(m$Tp_s) == 0 || sd(m$Tp_r) == 0) return(NA_real_)
  cor(m$Tp_s, m$Tp_r)
}

## NEW: decide, per depth and day, whether the IW correction is applied -----
apply_surface_rule <- function(res, daily, cfg) {
  daily[, `:=`(near_surface = FALSE, ref_depth = NA_real_, ref_r2 = NA_real_,
               ref_diel_frac = NA_real_, link_cor = NA_real_,
               corrected = status == "ok",
               iw_decision = fifelse(status == "ok", "corrected", paste0("not_fitted:", status)))]
  if (all(is.na(daily$depth)) || cfg$surface_rule == "always") return(daily)
  if (!cfg$surface_rule %in% c("skip", "hybrid"))
    stop("cfg$surface_rule must be 'always', 'skip' or 'hybrid'.", call. = FALSE)

  depths <- sort(unique(daily$depth))
  surf <- depths[depths <= cfg$surface_depth_max]
  if (!length(surf)) return(daily)
  daily[depth %in% surf, near_surface := TRUE]

  rule <- cfg$surface_rule
  ref <- NA_real_
  if (rule == "hybrid") {
    deeper <- setdiff(depths, surf)
    ref <- if (!is.null(cfg$ref_depth)) cfg$ref_depth else if (length(deeper)) min(deeper) else NA_real_
    if (!is.finite(ref) || !ref %in% depths) {
      warning("hybrid rule: no valid reference sensor (ref_depth = ", ref,
              "); near-surface sensors are left uncorrected.", call. = FALSE)
      rule <- "skip"
    }
  }
  if (rule == "skip") {
    daily[near_surface == TRUE, `:=`(corrected = FALSE, iw_decision = "skipped:near_surface")]
    return(daily)
  }

  bucket_s <- max(res[, median(diff(as.numeric(time))), by = depth]$V1)
  refd <- daily[depth == ref, .(day, ref_status = status, r2, diel_frac_T)]
  for (s in surf) {
    for (d in daily[depth == s, day]) {
      rr <- refd[day == d]
      lc <- link_cor_day(res[depth == s & day == d, .(time, T_prime)],
                         res[depth == ref & day == d, .(time, T_prime)],
                         bucket_s, cfg$link_remove_diel)
      st <- daily[depth == s & day == d, status]
      why <- if (st != "ok") paste0("not_fitted:", st)
      else if (!nrow(rr) || rr$ref_status != "ok") "no_IW_evidence:ref_not_fitted"
      else if (!is.finite(rr$r2) || rr$r2 < cfg$iw_r2_min) "no_IW_evidence:ref_low_R2"
      else if (is.finite(rr$diel_frac_T) && rr$diel_frac_T > cfg$diel_flag_thr) "no_IW_evidence:ref_T_diel"
      else if (!is.finite(lc) || lc < cfg$link_cor_min) "no_IW_evidence:surface_not_linked"
      else "corrected"
      daily[depth == s & day == d, `:=`(
        ref_depth = ref, ref_r2 = if (nrow(rr)) rr$r2 else NA_real_,
        ref_diel_frac = if (nrow(rr)) rr$diel_frac_T else NA_real_,
        link_cor = lc, corrected = why == "corrected", iw_decision = why)]
    }
  }
  daily
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

  ## NEW: near-surface rule -> days not corrected keep the raw DO
  daily <- apply_surface_rule(res, daily, cfg)
  res <- res[, !duplicated(names(res)), with = FALSE]   ## FIX: DO_mean appeared twice
  res[daily, on = .(depth, day), `:=`(corrected = i.corrected, near_surface = i.near_surface)]
  res[near_surface == TRUE & corrected == FALSE, DO_clean := DO]
  setorder(res, depth, time)

  ts <- res[, .(depth, time, day, DO, T, T_prime, DO_prime, DO_IM, DO_clean, corrected)]
  if (all(is.na(ts$depth))) { ts[, depth := NULL]; daily[, depth := NULL] }

  ## NEW: console summary per depth
  by_d <- if ("depth" %in% names(daily)) "depth" else NULL
  print(daily[, .(days = .N, fitted = sum(status == "ok"), corrected = sum(corrected),
                  median_R2 = round(median(r2, na.rm = TRUE), 2),
                  diel_flagged = sum(diel_flag)), by = by_d])
  if ("iw_decision" %in% names(daily) && any(daily$iw_decision != "corrected"))
    print(daily[iw_decision != "corrected", .N, by = c(by_d, "iw_decision")])
  if (any(daily$diel_flag & daily$corrected))
    message("Some corrected days have T' dominated by a 24 h cycle (diel_flag = TRUE): ",
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
  ## FIX: facet the daily fits by depth x day (depths were mixed in one panel)
  fd <- copy(out$data)
  fd[, applied := factor(corrected, c(TRUE, FALSE), c("applied", "not applied"))]
  p2 <- ggplot(fd, aes(T_prime, DO_prime)) +
    geom_point(size = 0.4, alpha = 0.4) +
    geom_line(aes(y = DO_IM, linetype = applied), colour = "firebrick") +
    scale_linetype_manual(values = c(applied = "solid", `not applied` = "dashed"), drop = FALSE) +
    labs(x = "T' (°C)", y = "DO'", linetype = "IW correction") + theme_bw() +
    theme(legend.position = "bottom")
  p2 <- if ("depth" %in% names(fd)) p2 + facet_grid(depth ~ day, scales = "free")
  else p2 + facet_wrap(~day, scales = "free")
  list(timeseries = p1, fits = p2)
}

## ---------------------------------------------------- SYNTHETIC TEST CASE ----
## NEW: three sensors (1, 9, 19 m), 30-min data. IW at 9 and 19 m every day.
## 1 m: days 1-4 in a mixed layer (diel heating, no IW); days 5-8 a shallow
## thermocline brings coherent IW displacements up to 1 m.
make_synthetic <- function(days = 8, step_min = 30, seed = 1) {
  set.seed(seed)
  time <- seq(as.POSIXct("2022-04-12", tz = "UTC"),
              by = step_min * 60, length.out = days * 1440 / step_min)
  tt  <- as.numeric(time - time[1]) / 3600                    # hours
  n   <- length(tt)
  zeta <- 1.5 * sin(2 * pi * tt / 7.3) + 0.8 * sin(2 * pi * tt / 11.1 + 1) +
    as.numeric(stats::filter(rnorm(n, 0, 0.3), rep(1/6, 6), circular = TRUE))
  strat <- tt >= 4 * 24                                       # shallow thermocline from day 5
  heat  <- 0.6 * sin(2 * pi * (tt - 10) / 24)                 # diel heating at 1 m
  bio   <- function(A) A * sin(2 * pi * (tt - 12) / 24)       # metabolic diel cycle
  mk <- function(z, DO0, T0, A, zfac, dTdz, dDOdz, extraT = 0) {
    DO_bio <- DO0 + bio(A)
    data.table(time = time, depth = z,
               T  = T0 + extraT + dTdz * zfac * zeta + rnorm(n, 0, 0.01),
               DO = DO_bio + dDOdz * zfac * zeta + 0.05 * (zfac * zeta)^2 + rnorm(n, 0, 0.02),
               DO_bio = DO_bio)
  }
  d <- rbind(mk(1,  11.5, 19, 0.5, ifelse(strat, 0.6, 0), -1.5, -0.6, heat),
             mk(9,   7.5, 14, 0.15, 1.0, -1.0, -0.9),
             mk(19,  3.5, 12.3, 0.05, 0.7, -0.3, -0.4))
  list(df = d[, .(time, DO, T, depth)], truth = d[, .(time, depth, DO_bio)])
}

## ------------------------------------------------------------------ RUN ----
if (run_main) {                                ## FIX: see run_main above
  ## FIX: moved here from the top of the file, and only used inside RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(repo_root)
  }
  if (is.null(cfg$input_file)) {
    message("No input_file set: running synthetic test (3 depths).")
    syn <- make_synthetic()
    df  <- syn$df
    cfg$depth_col <- "depth"
  } else {
    if (!file.exists(cfg$input_file)) stop("Input file not found: ", cfg$input_file, call. = FALSE)
    df <- fread(cfg$input_file)
  }

  dt  <- prepare_input(df, cfg)
  out <- correct_internal_waves(dt, cfg)

  if (is.null(cfg$input_file)) {                ## NEW: RMSE per depth and per surface rule
    rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
    for (rule in c("always", "skip", "hybrid")) {
      o <- out
      if (rule != cfg$surface_rule)
        invisible(capture.output(o <- suppressMessages(
          correct_internal_waves(dt, modifyList(cfg, list(surface_rule = rule))))))
      m <- merge(o$data, syn$truth, by = c("depth", "time"))
      cat(sprintf("rule %-6s | RMSE vs true bio DO by depth: %s\n", rule,
                  paste(m[, sprintf("%g m raw %.3f clean %.3f", depth[1], rmse(DO, DO_bio), rmse(DO_clean, DO_bio)),
                          by = depth]$V1, collapse = " | ")))
    }
  }

  ## FIX: no setwd("results/") (it moved one level deeper at each run) and
  ## cfg$output_fil -> cfg$output_file
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(out$data,  file.path(cfg$out_dir, cfg$output_file))
  fwrite(out$daily, file.path(cfg$out_dir, cfg$daily_file))

  if (cfg$make_plots) {
    pl <- plot_correction(out)
    if (!is.null(pl)) {
      n_dep <- if ("depth" %in% names(out$data)) uniqueN(out$data$depth) else 1
      n_day <- uniqueN(out$data$day)
      ggplot2::ggsave(file.path(cfg$out_dir, "DO_clean_timeseries.png"), pl$timeseries,
                      width = 10, height = 1 + 2.5 * n_dep)
      ggplot2::ggsave(file.path(cfg$out_dir, "DO_clean_daily_fits.png"), pl$fits,
                      width = max(8, 1.8 * n_day), height = 1.5 + 2 * n_dep)
    }
  }
}
