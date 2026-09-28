#!/usr/bin/env python3
"""
Spectral (frequency-domain) diel method for GPP, R and NEP from a single
high-frequency DO/T record, with removal of internal-wave (IW) signals.

Clean Python 3 rewrite of the "[1] Simple Amplitude approach" of DO_budget.py
(B. Fernández Castro, used in Fernández Castro et al. 2021, WRR,
doi:10.1029/2020WR029283; Supplement Text S1).

Pipeline (per sensor, hourly data):
  <X>  = 30-day Gaussian low-pass of X (as in DO_budget.py)
  NEP  = d<DO>/dt + <F_gas/zmix> - <transport>                [mmol m-3 d-1]
  For each centre time, on a 257-h window of DO' = DO-<DO>, T' = T-<T>:
    (mixed layer) DO'' = DO' + detrend(cumsum(detrend(F_gas/zmix)))   (Eq. S6)
    Step i  : S^ = S_DO - |S_DO,T / S_T|^2 S_T          (IW removal, co-spectrum)
    Step ii : S_n = a f^b fitted on 12-48 h excluding 18-28 h;
              var_DO24 = int_18-28h S^ - int_18-28h S_n
    Step iii: var_DO24 = 0 if int S^ < 1.5 int S_n
    Step iv : var_DO24 = 0 if DO-light phase < -2 h
  A24 = sqrt(2 var_DO24); R = 24 * 2 A24 / t_night; if NEP < 0: R += |NEP|;
  GPP = R + NEP

Conventions: DO in mmol m-3, T in degC, I0 in W m-2, k_gas in m h-1,
F_gas = k_gas (DO - DO_sat) in mmol m-2 h-1, POSITIVE = O2 leaving the lake.

What was dropped from DO_budget.py (site-specific to Lake Geneva/LéXPLORE):
heat-budget Kz, Rhône advection, mixed-layer entrainment, Soloviev piston
velocity, plotting. Provide their net effect as an optional 'transport'
column (mmol m-3 h-1) and k_gas directly if you have a better model.

Usage:
  python spectral_diel_method.py --input data.csv --z-sensor 10 --lat 46.5 \
         --elev 372 --out spectral_results.csv
Input CSV columns: time, DO, T, I0 [, zmix, k_gas | wind, transport]
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy import signal


# ------------------------------------------------------------------ config --
@dataclass
class SpectralConfig:
    nperseg: int = 128              # Welch sub-segment length (h)
    seglen: int = 256               # analysis window = seglen + 1 hourly points
    step_h: int = 1                 # evaluate every step_h hours
    period_min_h: float = 18.0      # diel band
    period_max_h: float = 28.0
    baseline_min_h: float = 12.0    # baseline fit range
    baseline_max_h: float = 48.0
    snr_min: float = 1.5            # step iii
    phase_min_h: float = -2.0       # step iv
    max_zero_frac_skip: float = 0.5
    max_zero_frac_flag: float = 0.1
    remove_iw: bool = True          # step i on/off
    gas_correction: bool = True     # Eq. S6 for mixed-layer sensors
    coh_min: float = 0.0            # NEW: 0 = original; >0 keeps only coherent bins in step i
    band_integration: str = "trapz" # "trapz" (original) | "rect"
    phase_bin: str = "original"     # "original" (2nd bin in band) | "nearest" (to 1 cpd)
    smooth_hw_days: float = 30.0    # half-width of the Gaussian low-pass


# ------------------------------------------------------- physics helpers --
def o2_sat(T, elev=0.0):
    """O2 saturation (Garcia & Gordon 1992, Benson & Krause), S=0, mmol m-3."""
    T = np.asarray(T, float)
    Ts = np.log((298.15 - T) / (273.15 + T))
    lnC = (5.80871 + 3.20291 * Ts + 4.17887 * Ts**2 + 5.10006 * Ts**3
           - 9.86643e-2 * Ts**4 + 3.80369 * Ts**5)
    rho = 1000 * (1 - (T + 288.9414) / (508929.2 * (T + 68.12963)) * (T - 3.9863)**2)
    P = (1 - 2.25577e-5 * elev) ** 5.25588
    return np.exp(lnC) * rho / 1000 * P


def schmidt_o2(T):
    return 1800.6 - 120.10 * T + 3.7818 * T**2 - 0.047608 * T**3


def k_gas_cole(u10, T, n=-0.5):
    """Cole & Caraco (1998) k600 scaled to O2; returns m h-1."""
    k600 = 2.07 + 0.215 * np.asarray(u10, float) ** 1.7
    return k600 * (schmidt_o2(T) / 600.0) ** n / 100.0


def day_length_h(time: pd.DatetimeIndex, lat: float):
    doy = time.dayofyear.values
    decl = np.deg2rad(23.44) * np.sin(2 * np.pi * (284 + doy) / 365)
    phi = np.deg2rad(lat)
    cosw = (np.sin(np.deg2rad(-0.833)) - np.sin(phi) * np.sin(decl)) / (np.cos(phi) * np.cos(decl))
    return 2 * np.rad2deg(np.arccos(np.clip(cosw, -1, 1))) / 15


def gauss_smooth(x, t_days, hw_days=30.0):
    """Gaussian low-pass of DO_budget.py: w = exp(-dt^2/(0.5 hw)^2), |dt| <= hw."""
    x = np.asarray(x, float)
    out = np.full(x.size, np.nan)
    ok = np.isfinite(x)
    for i in range(x.size):
        d = t_days - t_days[i]
        j = ok & (np.abs(d) <= hw_days)
        if j.any():
            w = np.exp(-d[j] ** 2 / (0.5 * hw_days) ** 2)
            out[i] = np.sum(w * x[j]) / np.sum(w)
    return out


def _epoch_days(t: pd.Series):
    # days since 1970-01-01, same time axis as the R version (identical window edges)
    return (t - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().values / 86400


def centered_diff(x, y):
    out = np.full(y.size, np.nan)
    ok = np.where(np.isfinite(x) & np.isfinite(y))[0]
    if ok.size < 2:
        return out
    x0, y0 = x[ok], y[ok]
    d = np.empty(ok.size)
    d[0] = (y0[1] - y0[0]) / (x0[1] - x0[0])
    d[-1] = (y0[-1] - y0[-2]) / (x0[-1] - x0[-2])
    d[1:-1] = (y0[2:] - y0[:-2]) / (x0[2:] - x0[:-2])
    out[ok] = d
    return out


# ------------------------------------------------------------ data prep --
def regularize(df: pd.DataFrame, step_min=60) -> pd.DataFrame:
    """Average onto a regular grid (hourly by default); gaps become NaN."""
    if "time" not in df:
        raise KeyError(f"'time' column missing. Available: {list(df.columns)}")
    d = df.copy()
    t = d["time"]
    d["time"] = (pd.to_datetime(t, unit="s", utc=True) if pd.api.types.is_numeric_dtype(t)
                 else pd.to_datetime(t, utc=True))            # epoch seconds or ISO strings
    d = d.set_index("time").select_dtypes("number")
    return d.resample(f"{step_min}min").mean().reset_index()


def add_physics_terms(df, z_sensor, elev=0.0, zmix_default=np.nan):
    d = df.copy()
    for c in ("time", "DO", "T", "I0"):
        if c not in d:
            raise KeyError(f"Column '{c}' missing. Available: {list(d.columns)}")
    if "transport" not in d:
        d["transport"] = 0.0
    if "zmix" not in d:
        d["zmix"] = zmix_default
    if "k_gas" not in d:
        d["k_gas"] = k_gas_cole(d["wind"], d["T"]) if "wind" in d else np.nan
    d["in_ml"] = z_sensor < d["zmix"]
    if d["in_ml"].any() and d["k_gas"].isna().all():
        raise ValueError("Sensor in the mixed layer but neither 'k_gas' nor 'wind' given.")
    d["F_gas"] = np.where(d["in_ml"], d["k_gas"] * (d["DO"] - o2_sat(d["T"], elev)), 0.0)
    d["gas_term"] = np.where(d["in_ml"] & (d["zmix"] > 0), d["F_gas"] / d["zmix"], 0.0)
    return d


def nep_budget_smooth(d, hw_days=30.0):
    """Fortnight-scale NEP (mmol m-3 d-1) from the low-passed DO budget."""
    td = _epoch_days(d["time"])
    sDO = gauss_smooth(d["DO"].values, td, hw_days)
    sgas = gauss_smooth(d["gas_term"].values, td, hw_days) * 24
    strn = gauss_smooth(d["transport"].values, td, hw_days) * 24
    return sDO, centered_diff(td, sDO) + sgas - strn


# ------------------------------------------------------- spectral core --
def _band_var(f, S, sel, method):
    if method == "trapz":
        return np.trapezoid(S[sel], f[sel]) if sel.sum() > 1 else 0.0
    return S[sel].sum() * (f[1] - f[0])


def spectral_window(xx, tt, irr, gas, cfg: SpectralConfig, fs=24.0):
    """Steps i-iv on one window. gas = F_gas/zmix (mmol m-3 d-1) or None."""
    nzeros = np.sum(xx == 0)
    zero_scale = xx.size / (xx.size - nzeros)

    xx1 = xx + signal.detrend(np.cumsum(signal.detrend(gas) / fs)) if gas is not None else xx

    f, S0 = signal.welch(xx1, fs, nperseg=cfg.nperseg)
    _, ST = signal.welch(tt, fs, nperseg=cfg.nperseg)
    _, CS = signal.csd(xx1, tt, fs, nperseg=cfg.nperseg)
    if cfg.remove_iw:
        with np.errstate(divide="ignore", invalid="ignore"):
            H2 = np.abs(CS / ST) ** 2
            coh2 = np.abs(CS) ** 2 / (S0 * ST)
        H2[~np.isfinite(coh2) | (coh2 < cfg.coh_min)] = 0.0
        Shat = S0 - H2 * ST                                   # step i
        Shat[~np.isfinite(Shat)] = S0[~np.isfinite(Shat)]
    else:
        Shat = S0.copy()

    f_lo, f_hi = fs / cfg.period_max_h, fs / cfg.period_min_h
    band = (f >= f_lo) & (f <= f_hi)

    fit = ((f >= fs / cfg.baseline_max_h) & (f < fs / cfg.baseline_min_h)
           & ~((f > f_lo) & (f <= f_hi)) & (Shat > 0))
    if fit.sum() >= 2:                                        # step ii
        pp = np.polyfit(np.log10(f[fit]), np.log10(Shat[fit]), 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            Sn = 10 ** np.polyval(pp, np.log10(f))
        Sn[~np.isfinite(Sn)] = 0.0
    else:
        Sn = np.zeros_like(f)

    var_raw = zero_scale * _band_var(f, S0, band, cfg.band_integration)
    var_hat = zero_scale * _band_var(f, Shat, band, cfg.band_integration)
    var_noise = zero_scale * _band_var(f, Sn, band, cfg.band_integration)
    var_do24 = max(var_hat - var_noise, 0.0)
    flag = True
    if var_hat < cfg.snr_min * var_noise:                     # step iii
        var_do24, flag = 0.0, False

    _, CL = signal.csd(xx, irr, fs, nperseg=cfg.nperseg)      # step iv
    ib = np.where(band)[0]
    ip = ib[1] if cfg.phase_bin == "original" else ib[np.argmin(np.abs(f[ib] - 1))]
    phase_h = np.angle(CL[ip], deg=True) * 12 / 180
    if phase_h < cfg.phase_min_h:
        var_do24 = 0.0
    if nzeros > cfg.max_zero_frac_flag * xx.size:
        flag = False
    return dict(var_raw=var_raw, var_hat=var_hat, var_noise=var_noise,
                var_do24=var_do24, phase_h=phase_h, amp_flag=flag,
                f=f, S_raw=S0, S_hat=Shat, S_n=Sn)


def spectral_metabolism(d: pd.DataFrame, cfg: SpectralConfig, lat: float) -> pd.DataFrame:
    """d: hourly frame from add_physics_terms(). Returns one row per centre time."""
    step = d["time"].diff().dt.total_seconds().median() / 3600
    if abs(step - 1) > 1e-6:
        raise ValueError(f"Hourly data expected (got {step} h); use regularize().")
    nt, half = len(d), cfg.seglen // 2
    if nt < cfg.seglen + 1:
        raise ValueError(f"Record shorter than one window ({cfg.seglen + 1} h).")

    td = _epoch_days(d["time"])
    sDO, NEP = nep_budget_smooth(d, cfg.smooth_hw_days)
    sT = gauss_smooth(d["T"].values, td, cfg.smooth_hw_days)
    pDO = np.nan_to_num(d["DO"].values - sDO, nan=0.0)
    pT = np.nan_to_num(d["T"].values - sT, nan=0.0)
    irr = np.nan_to_num(d["I0"].values, nan=0.0)
    gas_all = np.nan_to_num(d["gas_term"].values * 24, nan=0.0)
    in_ml = d["in_ml"].values.astype(float)
    t_night = 24 - day_length_h(pd.DatetimeIndex(d["time"]), lat)

    rows = []
    for j in range(0, nt, cfg.step_h):
        jj = np.arange(j - half, j + half + 1)
        if jj.min() < 0:
            jj += -jj.min()
        elif jj.max() >= nt:
            jj -= jj.max() - nt + 1
        xx = pDO[jj]
        if np.sum(xx == 0) >= cfg.max_zero_frac_skip * xx.size or not np.isfinite(d["DO"].iloc[j]):
            continue
        gas = gas_all[jj] if (cfg.gas_correction and np.nanmean(in_ml[jj]) > 0.5) else None
        w = spectral_window(xx, pT[jj], irr[jj], gas, cfg)
        rows.append(dict(time=d["time"].iloc[j], t_night=t_night[j], NEP=NEP[j],
                         **{k: w[k] for k in ("var_raw", "var_hat", "var_noise",
                                              "var_do24", "phase_h", "amp_flag")}))
    res = pd.DataFrame(rows)
    res["A24"] = np.sqrt(2 * res["var_do24"])
    res["R"] = 24 / res["t_night"] * 2 * res["A24"]
    res.loc[res["NEP"] < 0, "R"] += res.loc[res["NEP"] < 0, "NEP"].abs()
    res["GPP"] = res["R"] + res["NEP"]
    return res


# ------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", default="spectral_results.csv")
    ap.add_argument("--z-sensor", type=float, required=True, help="sensor depth (m)")
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--elev", type=float, default=0.0)
    ap.add_argument("--zmix-default", type=float, default=np.nan)
    ap.add_argument("--band-integration", choices=["trapz", "rect"], default="trapz")
    ap.add_argument("--coh-min", type=float, default=0.0)
    ap.add_argument("--no-iw-removal", action="store_true")
    ap.add_argument("--step-h", type=int, default=1)
    a = ap.parse_args()

    cfg = SpectralConfig(band_integration=a.band_integration, coh_min=a.coh_min,
                         remove_iw=not a.no_iw_removal, step_h=a.step_h)
    raw = pd.read_csv(a.input)
    keep = [c for c in ("time", "DO", "T", "I0", "zmix", "k_gas", "wind", "transport") if c in raw]
    d = add_physics_terms(regularize(raw[keep]), a.z_sensor, a.elev, a.zmix_default)
    res = spectral_metabolism(d, cfg, a.lat)
    res.to_csv(a.out, index=False)
    daily = res.assign(day=res["time"].dt.date).groupby("day")[["GPP", "R", "NEP"]].mean()
    print(f"Config: {asdict(cfg)}")
    print(f"{len(res)} windows | amp_flag False: {(~res['amp_flag']).sum()} | "
          f"var_DO24 = 0: {(res['var_do24'] == 0).sum()}")
    print("Mean daily rates (mmol m-3 d-1):")
    print(daily.mean().round(2).to_string())


if __name__ == "__main__":
    main()
