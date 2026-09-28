# o2-internalwaves

Lake metabolism (GPP, R, NEP) from high-frequency dissolved oxygen (DO), temperature and light
records, **after removing the signature of internal waves** from the oxygen signal, and closure of
the DO budget at each sensor depth.

The methods follow Fernández Castro et al. (2021), *Primary and net ecosystem production in a large
lake diagnosed from high-resolution oxygen measurements*, Water Resources Research,
[doi:10.1029/2020WR029283](https://doi.org/10.1029/2020WR029283), and its Supporting Information
(Text S1: spectral method, Text S2: time-domain method). A few extensions and corrections are
described in [Known issues and deviations from the paper](#known-issues-and-deviations-from-the-paper).

---

## Contents

1. [Why internal waves matter](#why-internal-waves-matter)
2. [What is in the repository](#what-is-in-the-repository)
3. [Requirements](#requirements)
4. [Quick start](#quick-start)
5. [Preparing input for the DO budget / metabolism script](#preparing-input-for-the-do-budget--metabolism-script)
6. [Running only the time-domain internal-wave correction](#running-only-the-time-domain-internal-wave-correction)
7. [Methods in brief](#methods-in-brief)
8. [Outputs](#outputs)
9. [Known issues and deviations from the paper](#known-issues-and-deviations-from-the-paper)
10. [References](#references)

---

## Why internal waves matter

The diel oxygen method estimates metabolism from the day–night oscillation of DO: oxygen rises
during the day (photosynthesis > respiration) and falls at night (respiration only). In stratified
lakes, internal waves move isotherms, and the oxygen gradients with them, up and down past a fixed
sensor. A sensor in the metalimnion can therefore see DO changes of several mg L⁻¹ within hours
that have nothing to do with biology. If these are not removed, GPP and R are badly biased.

Because internal waves displace temperature and oxygen together, the temperature record can be
used to identify and remove the physical part of the DO signal, either in the time domain
(regression of DO on T) or in the frequency domain (co-spectrum of DO and T).

---

## What is in the repository

| File | What it does | Run on its own? |
|---|---|---|
| `do_budget_metabolism.R` | **Main pipeline.** Reads a wide table (DO, T, PAR at several depths + wind + thermocline depth), corrects DO for internal waves (time domain), computes gas exchange and vertical diffusion, estimates GPP/R/NEP with the time-domain method and (optionally) the spectral method, and closes the DO budget at each sensor. The budget residual is interpreted as unknown physics (entrainment, lateral advection, ...). | Yes |
| `iw_correction_time_domain.R` | Time-domain internal-wave correction (Text S2): daily regression DO′ ~ T′ + T′², for one or several stacked sensors, with a rule for near-surface sensors. Output: `DO_clean`. | Yes |
| `iw_correction_frequency_domain.R` | Frequency-domain internal-wave correction that returns a `DO_clean` **time series** (an extension: the paper only cleans spectra). | Yes |
| `spectral_iw_correction.R` | Functions of the spectral scale-separation method (Text S1: Welch spectra, IW removal, baseline, phase test, R/GPP). Used by `do_budget_metabolism.R` and `compare_methods.R`. | No (library) |
| `spectral_diel_method.py` | Clean Python 3 rewrite of the spectral part of the original `DO_budget.py`. Gives the same numbers as the R version (checked to ~1e-12). | Yes (command line) |
| `compare_methods.R` | Single-sensor comparison of time-domain vs spectral methods on synthetic data with known truth. | Yes |
| `metabolism_functions.R`, `synthetic_data.R` | Helpers and synthetic data generator for `compare_methods.R`. | No (library) |

Dependencies between scripts:

```
do_budget_metabolism.R ──► iw_correction_time_domain.R
                       └─► spectral_iw_correction.R        (if run_spectral = TRUE)

compare_methods.R ──► metabolism_functions.R, spectral_iw_correction.R,
                      synthetic_data.R, iw_correction_time_domain.R
```

**Keep all scripts in the same folder.** Scripts that load other scripts do so without running
them (they set the R option `iw.library_mode` temporarily; see
[Troubleshooting](#troubleshooting)).

---

## Requirements

**R** (≥ 4.1) with:

- `data.table` (required)
- `ggplot2` (figures; the scripts run without it and skip the plots)
- `patchwork` (only for `compare_methods.R`)
- `rstudioapi` (optional: when run from RStudio, the scripts set the working directory
  automatically)

```r
install.packages(c("data.table", "ggplot2", "patchwork"))
```

**Python** (only for `spectral_diel_method.py`): Python ≥ 3.9, `numpy` ≥ 2.0, `scipy`, `pandas`.

---

## Quick start

Every runnable script has a **CONFIG block at the top**. Setting `input_file = NULL` runs a
**synthetic test with a known answer**. This is the easiest way to check that everything works
before using real data.

```r
# 1. DO budget and metabolism, synthetic test (time domain + spectral, compared with the truth)
source("do_budget_metabolism.R")

# 2. Only the internal-wave correction, synthetic test with 3 sensors (1, 9, 19 m)
source("iw_correction_time_domain.R")
```

For real data, edit `input_file` (and the other settings described below) in the CONFIG block and
run the script again.

---

## Preparing input for the DO budget / metabolism script

`do_budget_metabolism.R` reads **one CSV file in wide format: one row per timestamp, one column
per variable and sensor depth**.

### Column overview

| Column | Content | Units | Required |
|---|---|---|---|
| `time` | Timestamp | see [Time](#time) | yes |
| `DO_<z>` | Dissolved oxygen at depth `<z>` m | mg L⁻¹ (default) or mmol m⁻³ | yes, ≥ 1 depth (≥ 2 for diffusion) |
| `T_<z>` | Water temperature at depth `<z>` m | °C | yes, one for **every** DO depth |
| `PAR_<z>` | Photosynthetically active radiation at depth `<z>` m | µmol photons m⁻² s⁻¹ | yes, ≥ 1 depth |
| `wind` | Wind speed | m s⁻¹ | yes |
| `z_th` | Thermocline / mixed-layer depth | m, positive downwards | yes |

Example for three sensors at 1, 9 and 19 m, with PAR measured at 1 and 9 m:

```
time,DO_1,DO_9,DO_19,T_1,T_9,T_19,PAR_1,PAR_9,wind,z_th
2022-04-12 00:00:00,9.61,10.32,7.80,17.02,12.61,10.00,0,0,2.8,5.0
2022-04-12 00:30:00,9.60,10.35,7.80,17.00,12.55,10.01,0,0,2.6,5.0
2022-04-12 01:00:00,9.58,10.30,7.79,16.98,12.70,10.00,0,0,3.1,5.0
...
```

### Column names and sensor depths

- The depth is read from the column name: `DO_9` → 9 m, `T_9.5` → 9.5 m. Use a **dot** as the
  decimal separator and **no unit** in the name (`DO_9m` or `DO_9,5` will not be recognised).
- The depth in the name must be the **actual sensor depth below the surface** (m), because it is
  used for layer thicknesses, gradients, light extrapolation and the near-surface rule.
- Each DO depth needs a temperature column **at the same depth** (`DO_9` ↔ `T_9`; `T_9.0` also
  matches). If a DO logger has no own thermistor, use the closest thermistor and name the column
  with the DO depth.
- Extra thermistors without a DO sensor (e.g. `T_5`) are ignored by the budget. They are useful
  for computing `z_th` beforehand (see below).
- The prefixes can be changed in the CONFIG block (`do_prefix`, `t_prefix`, `par_prefix`), as can
  the names of the `time`, `wind` and `z_th` columns (`time_col`, `wind_col`, `zth_col`).

### Time

- Accepted formats: `2022-04-12 00:30:00`, `2022-04-12 00:30`, `2022-04-12T00:30:00Z`,
  `12/04/2022 00:30`, or Unix time in seconds. For any other format, set `time_format`
  (e.g. `"%d.%m.%Y %H:%M"`).
- **Time zone:** set `tz`. Days are split at `day_start_hour` (default midnight) **in this time
  zone**. Use a time zone in which midnight falls during the night (local standard time or UTC for
  European lakes), and make sure DO, T and PAR share the same clock. A shift between PAR and DO
  directly biases the light–productivity fit and the spectral phase test.
- **Sampling interval:** constant, between a few minutes and 1 hour (e.g. 10 or 30 min). Missing
  timestamps are allowed: the script builds a regular grid from the median interval and fills
  missing rows with NA. Data coarser than hourly are not supported.
- One row per timestamp, no duplicates.

### Dissolved oxygen

- Concentration, not % saturation. Convert % saturation to concentration first, using the
  saturation concentration at the in-situ temperature and pressure.
- Set `do_units = "mg/L"` (default) or `"mmol/m3"`. Internally everything is converted to
  mmol m⁻³. Rates are reported in mmol m⁻³ d⁻¹ (`out_units = "mmol"`) or mg L⁻¹ d⁻¹
  (`out_units = "mg"`).
- Remove obvious sensor artefacts (cleaning spikes, drift after biofouling) before running. The
  internal-wave correction does not remove them.

### PAR (light)

- PAR in **µmol m⁻² s⁻¹**, at one or more depths. `PAR_0` can be used for a surface sensor.
- **Night-time values must be 0, not NA.** NA hours are dropped from the metabolism fit.
- DO depths without a PAR sensor get PAR extrapolated from the nearest PAR sensor with
  PAR(z) = PAR(z_ref) · exp(−kd (z − z_ref)), where kd is fitted each day on all PAR sensors
  (log-linear). With **a single PAR sensor**, kd cannot be fitted and `kd_default` (m⁻¹) is used:
  set it to a realistic value for your lake (e.g. from a PAR profile or Secchi depth,
  kd ≈ 1.7 / Secchi).
- If only **shortwave radiation** (W m⁻²) is available, convert it approximately with
  PAR ≈ 2.1 × SW (µmol m⁻² s⁻¹), i.e. ~45 % of shortwave is PAR, at ~4.6 µmol J⁻¹. Name it
  `PAR_0`.
- Night length for the spectral method is taken from the shallowest PAR sensor (hours with
  PAR > `par_min`), so no latitude is needed.

### Wind

- Wind speed in m s⁻¹ at the height given by `wind_height` (m). It is converted to 10 m with a
  neutral log profile.
- The gas-transfer velocity is computed with Cole & Caraco (1998, default) or Vachon & Prairie
  (2013, `k600_model = "vachon"`, which also uses `lake_area_km2`).
- Gas exchange is applied only to sensors in the mixed layer (depth < `z_th`), using the DO and T
  of the **shallowest** sensor.

### Thermocline / mixed-layer depth `z_th`

- One value per row (m, positive downwards). It may change in time: deepening of the mixed layer
  is what produces entrainment. Short gaps are filled by linear interpolation.
- It decides, at every time step, which sensors are in the mixed layer (gas exchange, well-mixed
  box) and which are below (diffusion between layers).
- If you do not have it, compute it from your thermistors, for example as the depth where
  temperature is 1 °C below the value of the top sensor (criterion used in the paper, Fig. S3),
  or with `rLakeAnalyzer::thermo.depth()`:

```r
library(data.table)
w <- fread("my_wide_file.csv")
Tcols <- grep("^T_[0-9.]+$", names(w), value = TRUE)
z <- as.numeric(sub("^T_", "", Tcols)); o <- order(z); z <- z[o]; Tcols <- Tcols[o]
Tm <- as.matrix(w[, ..Tcols])
w[, z_th := apply(Tm, 1, function(Tz) {
  ok <- is.finite(Tz); if (sum(ok) < 2) return(NA_real_)
  zz <- z[ok]; tt <- Tz[ok]; target <- tt[1] - 1
  i <- which(tt <= target)[1]
  if (is.na(i)) return(max(zz))       # no thermocline above the deepest sensor
  approx(tt[(i - 1):i], zz[(i - 1):i], xout = target)$y
})]
fwrite(w, "my_wide_file.csv")
```

  With only a few thermistors this estimate is coarse. Use a full thermistor chain if one is
  available.

### Site settings in the CONFIG block

| Setting | Meaning |
|---|---|
| `elev` | Lake elevation (m a.s.l.), for O₂ saturation |
| `lake_area_km2` | Lake surface area (km²), for the Hondzo & Stefan diffusivity (and Vachon k600) |
| `z_bottom` | Bottom of the deepest sensor's layer (m), e.g. the maximum depth or the depth where that sensor stops being representative |
| `wind_height`, `k600_model` | See [Wind](#wind) |
| `kz_const` | Constant Kz (m² s⁻¹) to use instead of Hondzo & Stefan, if you have an independent estimate (e.g. from a heat budget) |

### Method settings you will most likely touch

| Setting | Default | Meaning |
|---|---|---|
| `fit_days` | 3 | Number of consecutive days sharing the metabolic parameters (Pmax, α, R20). `1` reproduces the paper, but then the daily residual is ~0 by construction (see [Methods](#methods-in-brief)). |
| `r2_IW_max` | 0.75 | Days where internal waves explain more than this fraction of DO variance are not fitted (Text S2). `NULL` switches this off. |
| `iw` | hybrid, 3 m | Settings passed to the internal-wave correction (see next section). |
| `run_spectral` | TRUE | Also run the spectral method. Needs ≥ ~11 days of data. |
| `spectral` | trapz | `list(band_integration = "rect")` gives a less biased diel amplitude than the original code. |
| `spectral_iw_rule` | follow_td | Do not apply the spectral internal-wave removal at near-surface depths that the time-domain rule left uncorrected. |
| `spec_smooth_hw_days` | 30 | Low-pass window for the spectral NEP. Reduce it (e.g. 10) for records shorter than 2–3 months. |

### Record length

- Time-domain method: works on any number of whole days.
- Spectral method: needs at least one 257-h window (~11 days) and is skipped with a message
  otherwise. Its NEP uses a 30-day low-pass, so it is only meaningful for records much longer than
  that, or with a smaller `spec_smooth_hw_days`.

### Checklist before running

- [ ] One row per timestamp, constant interval ≤ 1 h, one clock for all variables
- [ ] `DO_<z>` and `T_<z>` for every sensor, depths in the names in metres
- [ ] DO as concentration, `do_units` set
- [ ] PAR in µmol m⁻² s⁻¹, zeros at night, `kd_default` set if only one PAR sensor
- [ ] `wind` in m s⁻¹ and `wind_height` set
- [ ] `z_th` in metres, no long gaps
- [ ] `elev`, `lake_area_km2`, `z_bottom` set
- [ ] Missing values as empty cells or `NA` (not −9999)

---

## Running only the time-domain internal-wave correction

`iw_correction_time_domain.R` can be used alone, for example to clean DO before your own
metabolism model.

### Input: one CSV in long format

One row per timestamp **and** sensor:

| Column | Content |
|---|---|
| `time` | Timestamp (same formats as above) |
| `DO` | Dissolved oxygen, any unit. `DO_clean` is returned in the same unit. |
| `T` | Water temperature (°C) at the same sensor |
| `depth` | Sensor depth (m). Only needed if several sensors are stacked in the file. |

```
time,DO,T,depth
2022-04-12 00:00:00,9.61,17.02,1
2022-04-12 00:30:00,9.60,17.00,1
...
2022-04-12 00:00:00,10.32,12.61,9
2022-04-12 00:30:00,10.35,12.55,9
...
```

Column names are set in the CONFIG block (`time_col`, `do_col`, `temp_col`, `depth_col`).
For a single sensor, set `depth_col = NULL`. Sensors may have different sampling intervals, which
are handled per depth.

To convert a wide file (as used by the budget script) to this long format:

```r
library(data.table)
w  <- fread("my_wide_file.csv")
do <- melt(w, id.vars = "time", measure.vars = patterns("^DO_"), variable.name = "s", value.name = "DO")
tt <- melt(w, id.vars = "time", measure.vars = patterns("^T_"),  variable.name = "s", value.name = "T")
do[, depth := as.numeric(sub("^DO_", "", s))]
tt[, depth := as.numeric(sub("^T_",  "", s))]
long <- merge(do[, .(time, depth, DO)], tt[, .(time, depth, T)], by = c("time", "depth"))
fwrite(long, "my_long_file.csv")
```

### Settings

| Setting | Default | Meaning |
|---|---|---|
| `input_file` | — | Path to the long CSV. `NULL` runs the synthetic 3-sensor test. |
| `tz`, `day_start_hour` | UTC, 0 | Time zone and start hour of the 24-h segments |
| `min_coverage` | 0.8 | Minimum fraction of a day's samples needed to fit that day |
| `on_fail` | `"NA"` | `DO_clean` on days that cannot be fitted: `NA` or the raw DO (`"raw"`) |
| `surface_rule` | `"hybrid"` | Treatment of near-surface sensors: `"always"` (correct like the others), `"skip"` (never correct), `"hybrid"` (correct only on days with evidence that internal waves reach the sensor) |
| `surface_depth_max` | 3 | Sensors at or above this depth (m) count as near-surface |
| `ref_depth` | `NULL` | Reference sensor for the hybrid rule. Default: the shallowest sensor deeper than `surface_depth_max`. It should sit in the metalimnion. |
| `iw_r2_min` | 0.3 | Hybrid: internal waves are considered present at the reference sensor if R²(DO′ ~ T′) ≥ this |
| `link_cor_min` | 0.5 | Hybrid: minimum correlation between the near-surface and reference temperature anomalies (24-h cycle removed) |
| `out_dir` | `results` | Output folder |

With `surface_rule = "hybrid"`, a near-surface sensor is corrected on a given day only if (1) the
reference sensor shows internal waves and (2) the near-surface temperature moves together with the
reference temperature. Otherwise `DO_clean = DO` for that day. This avoids the main failure of the
regression near the surface: diel heating and photosynthesis both follow the sun, so the
regression would remove the metabolic signal as if it were an internal wave.

### Output

- `DO_clean.csv`: `depth, time, day, DO, T, T_prime, DO_prime, DO_IM, DO_clean, corrected`
- `DO_clean_daily_fits.csv`: one row per depth and day. It contains the regression coefficients
  (`a0, a1, a2`), `r2`, `status`, `diel_frac_T` and `diel_flag` (share of T′ variance explained by
  a 24-h cycle; high values warn that the correction may also remove biology), and, for the
  near-surface rule, `near_surface, ref_depth, ref_r2, link_cor, corrected, iw_decision` (the
  reason why a day was or was not corrected).
- `DO_clean_timeseries.png`, `DO_clean_daily_fits.png`

When run from RStudio, the working directory is set to the **parent folder** of the script (the
repository root), so `results/` is created there.

---

## Methods in brief

### Time-domain internal-wave correction (Text S2)

For each 24-h segment and sensor, anomalies with respect to the daily means are computed
(T′ = T − ⟨T⟩, DO′ = DO − ⟨DO⟩). The part of DO′ explained by vertical displacements is modelled
as DO′_IM = a₀ + a₁T′ + a₂T′², and removed: DO_clean = DO′ − DO′_IM + ⟨DO⟩.

### Spectral scale separation (Text S1)

The DO record is split into its periods (Welch spectra on ~11-day windows). DO variance coherent
with temperature is removed as internal waves, a background baseline is subtracted around 24 h,
and the remaining 24-h peak gives the amplitude A of the diel oxygen cycle. Night-time respiration
is the night-time drop (2A) divided by the night length. NEP comes from the slow (30-day) DO
budget, and GPP = R + NEP. It gives robust ~10-day averages, not day-to-day values.

### DO budget at each sensor (`do_budget_metabolism.R`)

```
dDO/dt  =  NEP  +  gas  +  diffusion  +  residual
```

- **dDO/dt** from the internal-wave-corrected DO, computed within each day.
- **NEP = GPP − R**, from the fit of hourly rates to
  nep = Pmax · tanh(α · PAR / Pmax) − R20 · θ^(T − 20) (Eq. S7).
- **gas** = −F/z_th for sensors in the mixed layer, with F = k_O₂ (DO − DO_sat).
- **diffusion**: Kz ∂DO/∂z across the thermocline and between sensors below it. Kz follows
  Hondzo & Stefan (1993), Kz = a (N²)^−0.43. N² and the DO gradients use 24-h smoothed profiles.
  Mixed-layer sensors share one well-mixed box from the surface to z_th, and deeper sensors
  represent the layer between the midpoints to their neighbours.
- **residual**: everything else (entrainment, lateral advection, sediment fluxes, model error).

**Important.** Respiration is almost constant within a day. If the metabolic parameters are
fitted day by day (`fit_days = 1`, as in the paper), any physical flux that lasts the whole day is
absorbed into R or GPP, and the daily residual is ~0 by construction. With `fit_days > 1` the
parameters are shared over several days, and day-to-day physical events appear in the residual
instead. The residual then also contains genuine day-to-day changes in biology, so interpret it
with care.

### What the synthetic tests showed (summary)

- The time-domain method gave the least biased GPP and R in the productive layer.
- The spectral method, as in the original code, underestimates R and GPP, for three reasons:
  trapezoid integration, spurious coherence in the internal-wave removal, and the sinusoid
  assumption.
- The spectral method's NEP is smoother, and it correctly returns GPP ≈ 0 in the dark.
- Near the surface, both internal-wave corrections can remove the biological signal when diel
  heating makes temperature and oxygen co-vary. Hence the hybrid rule and `follow_td`.
- With `fit_days = 3`, the budget residual tracked imposed entrainment, inflow and intrusion
  events (r ≈ 0.7–0.8). With `fit_days = 1` it did not (r ≈ 0).

---

## Outputs

`do_budget_metabolism.R` writes to `out_dir` (default `results_budget/`):

| File | Content |
|---|---|
| `budget_hourly_mmol.csv` | Per time step and depth: `DO`, `DO_clean`, `T`, `PAR`, `z_th`, `in_ml`, layer thickness `h`, `dDOdt`, `gas`, `diffusion`, `nep_obs`, modelled `gpp_mod` and `r_mod`, `residual` (mmol m⁻³ h⁻¹) |
| `budget_daily.csv` | Per day and depth: `GPP`, `R`, `NEP`, `dDOdt`, `gas`, `diffusion`, `residual`, `flag` (`ok` or fallback reason), fit parameters, `r2_IW`, `iw_decision`, `kd`, `z_eu`, `productive`, `layer`, `closure` (≈ 0, a consistency check) |
| `budget_productive_layer_areal.csv` | Budget terms integrated over the productive layer (mixed layer + sensors above the 1 % light depth), per m² |
| `Kz_interfaces_m2s.csv` | Kz at each interface between sensors (m² s⁻¹) |
| `gas_flux_mmol_m2_h.csv` | Air–lake flux `F_gas` and `k_O2` |
| `iw_correction_daily.csv` | Daily diagnostics of the internal-wave correction |
| `spectral_daily.csv`, `spectral_hourly.csv`, `spectral_productive_layer_areal.csv` | Spectral-method results (if `run_spectral = TRUE`) |
| `*.png` | Budget terms, metabolism per depth, productive layer, time domain vs spectral |

**Sign conventions.** In the budget files every term is a contribution to dDO/dt: positive means
it adds oxygen (e.g. `gas > 0` means the lake takes up O₂ from the air). In
`gas_flux_mmol_m2_h.csv`, `F_gas > 0` means **outgassing**. R is reported as a positive number.

**Fallback days (`flag` ≠ `ok`)**, following Text S2: the internal-wave signal is too dominant,
the fit failed, or a rate came out negative. On these days, GPP = NEP and R = 0 if NEP > 0, or
R = −NEP and GPP = 0 if NEP < 0, with NEP from the daily mean of the observed rates.

---

## Known issues and deviations from the paper

Found while porting the original Python code (`DO_budget.py`) and testing on synthetic data:

1. **Band integration (spectral).** The original code integrates the diel peak with the trapezoid
   rule over three frequency bins, which captures ~66 % of the variance of a pure sinusoid
   (amplitude −19 %). `band_integration = "rect"` captures ~95 %. `"trapz"` is kept as the default
   for reproducibility.
2. **Spurious coherence (spectral).** With three Welch segments, two unrelated series have a mean
   squared coherence of ~0.33. The internal-wave removal therefore also removes ~1/3 of
   *unrelated* DO variance. `coh_min` restricts the removal to significantly coherent frequencies
   (not in the paper).
3. **Phase frequency (spectral).** The DO–light phase is evaluated at 1.125 cycles per day
   (21.3 h), not at exactly 24 h. `phase_bin = "nearest"` uses the bin closest to 24 h.
4. **GPP when NEP < 0 (spectral).** The code sets R = R_spectral + |NEP|, i.e. GPP = R_spectral.
   This is kept as in the code, whereas the Supplement text states GPP = NEP + R.
5. **Site-specific terms** of the original code (heat-budget Kz, Rhône advection, entrainment,
   Soloviev piston velocity) are replaced by generic ones (Hondzo & Stefan Kz, wind-based k,
   residual term).
6. **Near-surface sensors.** Both corrections can confuse diel heating with internal waves. See
   the hybrid rule and `spectral_iw_rule`.
7. **Main-paper equations.** Eqs. 5, 7, 8 and 11 of the main paper were not checked line by line.
   The sign conventions used are the ones given under [Outputs](#outputs).

---

## Troubleshooting

- **A script runs but does nothing.** The R option `iw.library_mode` is still `TRUE` in your
  session, which disables the RUN blocks. Run `options(iw.library_mode = NULL)` or restart R.
  Current versions always reset this option, even when loading fails.
- **"Missing column(s)" or "No DO_<depth> columns found".** Check the column names against
  [Column names and sensor depths](#column-names-and-sensor-depths).
- **"time value(s) not parsed".** Set `time_format` in the CONFIG block.
- **Most days flagged `fallback:iw`.** Internal waves explain > 75 % of the DO variance at that
  sensor. Check the `r2` column of `iw_correction_daily.csv`, and relax or disable `r2_IW_max` if
  you accept the risk.
- **Spectral method skipped.** The record is shorter than ~11 days.

---

## References

- Fernández Castro, B., Chmiel, H. E., Minaudo, C., Krishna, S., Perolo, P., Rasconi, S., &
  Wüest, A. (2021). Primary and net ecosystem production in a large lake diagnosed from
  high-resolution oxygen measurements. *Water Resources Research*, 57, e2020WR029283.
- Cole, J. J., & Caraco, N. F. (1998). Atmospheric exchange of carbon dioxide in a low-wind
  oligotrophic lake measured by the addition of SF6. *Limnology and Oceanography*, 43, 647–656.
- Garcia, H. E., & Gordon, L. I. (1992). Oxygen solubility in seawater: Better fitting equations.
  *Limnology and Oceanography*, 37, 1307–1312.
- Hanson, P. C., et al. (2008). Evaluation of metabolism models for free-water dissolved oxygen
  methods in lakes. *Limnology and Oceanography: Methods*, 6, 454–465.
- Hondzo, M., & Stefan, H. G. (1993). Lake water temperature simulation model. *Journal of
  Hydraulic Engineering*, 119, 1251–1273.
- Vachon, D., & Prairie, Y. T. (2013). The ecosystem size and shape dependence of gas transfer
  velocity versus wind speed relationships in lakes. *Canadian Journal of Fisheries and Aquatic
  Sciences*, 70, 1757–1764.
- Wanninkhof, R. (1992). Relationship between wind speed and gas exchange over the ocean.
  *Journal of Geophysical Research*, 97, 7373–7382.
- Welch, P. D. (1967). The use of fast Fourier transform for the estimation of power spectra.
  *IEEE Transactions on Audio and Electroacoustics*, 15, 70–73.
