"""
Estimate T_sid (Moon sidereal period) directly from Tycho's Moon declinations,
without any orbit fit.

Idea: dec(t) is nearly sinusoidal with period T_sid.  With ~70 points spread
over 275 days (~10 cycles), a Lomb-Scargle periodogram finds the period to
high precision.

Inputs:  moon_dec_observations.dat  (iso_utc source tag dec_deg jd_tt)
"""
from __future__ import annotations

import numpy as np
from astropy.timeseries import LombScargle

from fit_moon_orbit import load_observations


def main() -> None:
    jds, decs, times, _ = load_observations("moon_dec_observations.dat")
    t_days = jds - jds[0]
    y = np.asarray(decs, dtype=float)

    print(f"{len(t_days)} obs; baseline {t_days[-1]:.2f} days")

    # Search T_sid in [20, 35] days -> frequencies [1/35, 1/20] cycles/day
    f_min, f_max = 1.0 / 35.0, 1.0 / 20.0
    freqs = np.linspace(f_min, f_max, 200_000)
    power = LombScargle(t_days, y).power(freqs)

    f_peak = freqs[np.argmax(power)]
    T_peak = 1.0 / f_peak

    # Parabolic refinement around peak
    k = int(np.argmax(power))
    if 0 < k < len(power) - 1:
        y0, y1, y2 = power[k - 1], power[k], power[k + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom != 0.0:
            dk = 0.5 * (y0 - y2) / denom
            f_peak = freqs[k] + dk * (freqs[1] - freqs[0])
            T_peak = 1.0 / f_peak

    print(f"Lomb-Scargle peak:  T_sid = {T_peak:.5f} days   "
          f"(f = {f_peak:.7f} cycles/day)")
    print(f"Reference modern   T_sid = 27.32166 days")

    # Bootstrap uncertainty
    rng = np.random.default_rng(0)
    N = 500
    Ts = np.empty(N)
    for b in range(N):
        idx = rng.integers(0, len(t_days), len(t_days))
        ts = t_days[idx]
        ys = y[idx]
        p = LombScargle(ts, ys).power(freqs)
        k = int(np.argmax(p))
        f_b = freqs[k]
        if 0 < k < len(p) - 1:
            y0, y1, y2 = p[k - 1], p[k], p[k + 1]
            d = y0 - 2.0 * y1 + y2
            if d != 0.0:
                dk = 0.5 * (y0 - y2) / d
                f_b = freqs[k] + dk * (freqs[1] - freqs[0])
        Ts[b] = 1.0 / f_b
    print(f"Bootstrap sigma(T_sid) = {Ts.std():.5f} days  "
          f"({Ts.std() * 86400:.2f} s)")


if __name__ == "__main__":
    main()
