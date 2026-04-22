"""
Extended Kalman Filter fit of the Moon's Keplerian elements to Tycho
Brahe's 1586 declination observations, using only observed quantities.

Design
------
State x (6):  [a, e, i, Omega, omega, M0]    defined at EPOCH = 1586-06-15 TT.
These are constants of motion in the 2-body approximation, so the
dynamical propagation is the identity map (F = I).  Small process noise
Q is added each step to allow the filter to absorb unmodelled effects
(perturbations, measurement bias) rather than locking up.

Measurement h(x):  topocentric apparent declination at Uraniborg,
produced by ``predicted_dec_array`` from ``fit_moon_orbit``.  No ephemeris,
no Moon mass, no prior distance are used.  The mean motion n = 2*pi/T_sid
with T_sid = 27.321661 d is fixed from naked-eye star-transit timings.

Each observation is processed sequentially in time:

    x-  = x+_(k-1)                              (propagate: F = I)
    P-  = P+_(k-1) + Q
    H   = dh/dx   (finite differences)
    y   = z_k - h(x-)                           (innovation, arcmin)
    S   = H P- H' + R
    K   = P- H' / S
    x+  = x- + K y
    P+  = (I - K H) P-

We pass through all observations, then take the final state as the new
initial state and iterate.  If the filter is well-conditioned the
outer iteration converges; if not, the sequence of (x, RMS) tells us so.

The result is an initial guess for a Keplerian orbit; the user can
inspect the trace of a, e, i, RMS across outer iterations to see whether
the estimator is converging or drifting.
"""
import warnings
import numpy as np
from astropy.time import Time
from astropy.utils.iers import conf as iers_conf

from fit_moon_orbit import (
    load_observations,
    predicted_dec_array,
    N_MEAN,
    EPOCH,
    URANIBORG,
)

warnings.filterwarnings('ignore')
iers_conf.auto_download = False

R_EARTH_KM = 6378.137


def h_scalar(x, t_single: Time) -> float:
    """Predicted topocentric apparent dec (deg) at a single time."""
    return float(predicted_dec_array(tuple(x), t_single)[0])


def jacobian(x, t_single: Time, steps=None) -> np.ndarray:
    """Finite-difference dh/dx as a (1,6) row vector of arcmin/unit."""
    if steps is None:
        # Scaled step per element
        steps = np.array([
            2_000.0,          # a, km
            1e-4,             # e
            np.radians(0.02), # i
            np.radians(0.02), # Omega
            np.radians(0.02), # omega
            np.radians(0.02), # M0
        ])
    h0 = h_scalar(x, t_single)
    H = np.zeros(6)
    for k in range(6):
        xp = x.copy(); xp[k] += steps[k]
        xm = x.copy(); xm[k] -= steps[k]
        hp = h_scalar(xp, t_single)
        hm = h_scalar(xm, t_single)
        # dh/dx in deg/unit  -> convert to arcmin/unit
        H[k] = (hp - hm) / (2.0*steps[k]) * 60.0
    return H


def project_state(x):
    """Project state to a physically meaningful region.  Instead of hard
    clamping ``e`` at 0 (which traps the filter), a negative ``e`` is
    reflected: negate it and shift ``omega`` by pi (geometrically the
    same orbit, viewed from the opposite apsis).  This lets the filter
    pass smoothly through e=0 in either direction."""
    x = x.copy()
    # a > 0; physical floor is R_earth (Moon can't be inside Earth).
    # Upper bound is loose -- any realistic lunar distance is < 10 R_lunar_SOI.
    x[0] = np.clip(x[0], 3*R_EARTH_KM, 10_000_000.0)
    # e reflection through 0
    if x[1] < 0.0:
        x[1] = -x[1]
        x[4] = x[4] + np.pi        # omega -> omega + pi
        x[5] = x[5] + np.pi        # M0    -> M0 + pi (apsis flipped)
    x[1] = min(x[1], 0.5)
    # i in [0, 60 deg] (Moon won't go higher)
    x[2] = np.clip(x[2], 0.0, np.radians(60.0))
    # angles mod 2*pi
    for k in (3, 4, 5):
        x[k] = x[k] % (2.0 * np.pi)
    return x


def ekf_pass(x0, P0, Q, R_meas, times: Time, decs_deg: np.ndarray,
             verbose=False):
    """Run one forward pass of the EKF.  Returns (x_final, P_final,
    residuals_before, residuals_after).  Residuals are in arcmin.
    `decs_deg` are measurements in deg, `R_meas` is variance in arcmin^2."""
    x = x0.copy()
    P = P0.copy()
    n = len(decs_deg)
    res_before = np.zeros(n)
    res_after  = np.zeros(n)
    for k in range(n):
        t_k = times[k:k+1]
        z = decs_deg[k]
        # Predict
        P = P + Q
        # Innovation (arcmin)
        h0 = h_scalar(x, t_k)
        y_before = (z - h0) * 60.0
        res_before[k] = y_before
        # Measurement Jacobian (arcmin/unit)
        H = jacobian(x, t_k)
        S = float(H @ P @ H.T + R_meas)
        K = (P @ H.T) / S                      # (6,)
        # Clip innovation to avoid huge excursions on first look-at of a
        # new epoch (Mahalanobis gating at 5-sigma)
        y_clip = float(np.clip(y_before, -5.0*np.sqrt(S), 5.0*np.sqrt(S)))
        x = x + K * y_clip
        x = project_state(x)
        P = (np.eye(6) - np.outer(K, H)) @ P
        # Symmetrize + enforce PSD by floor on diag
        P = 0.5 * (P + P.T)
        # Posterior residual
        h1 = h_scalar(x, t_k)
        res_after[k] = (z - h1) * 60.0
        if verbose:
            print(f'    k={k:3d}  t={t_k[0].iso[:19]}  '
                  f'innov={y_before:+7.2f}\'  '
                  f'post={res_after[k]:+6.2f}\'  '
                  f'a={x[0]:7.0f}  e={x[1]:.4f}')
    return x, P, res_before, res_after


def run(n_outer=6, verbose_inner=False):
    jds, decs, times, sources = load_observations()
    # Sort by time
    order = np.argsort(jds)
    jds = jds[order]; decs = decs[order]
    times = Time(jds, format='jd', scale='tt')
    sources = [sources[i] for i in order]
    n_obs = len(decs)
    print(f'{n_obs} observations from {times[0].iso} to {times[-1].iso}')

    # --- Initial state (loose but bounded; no ephemeris) ------------------
    x = np.array([
        0.0,               # a, km    (in ballpark)
        0.00,                    # e
        np.radians(0.0),         # i
        np.radians(0.0),         # Omega
        np.radians(0.0),         # omega
        np.radians(0.0),         # M0
    ])

    # Initial covariance: wide but not absurd.  These widths set how much
    # the filter is willing to move based on the first few observations.
    P0 = np.diag([
        (100_000.0)**2,            # a: +/- 100 000 km
        (0.10)**2,                 # e: +/- 0.10
        np.radians(15.0)**2,       # i: +/- 15 deg
        np.radians(60.0)**2,       # Omega: +/- 60 deg
        np.radians(60.0)**2,       # omega: +/- 60 deg
        np.radians(180.0)**2,      # M0: uniform-ish
    ])

    # Process noise per step.  Tiny on invariants; measurable on M0 to
    # absorb mean-anomaly drift from unmodelled perturbations.
    Q = np.diag([
        (5.0)**2,                  # a: +/- 5 km/step
        (1e-5)**2,                 # e
        np.radians(1e-4)**2,       # i
        np.radians(1e-4)**2,       # Omega
        np.radians(1e-4)**2,       # omega
        np.radians(2e-3)**2,       # M0: the one angle that may drift
    ])

    # Measurement noise (arcmin^2).  Tycho naked-eye + 2-body modelling
    # error (evection up to 78', variation up to 40') => a few arcmin
    # typical post-fit.  Starting with a larger R damps gain so early
    # updates can't drive state out of physical range; we'll decrease
    # across outer iterations (simulated annealing).
    R_meas_start = (20.0) ** 2
    R_meas_end   = (1.5)  ** 2

    trace = []
    P = P0.copy()
    print(f'initial : a={x[0]:8.1f} km  e={x[1]:6.4f}  '
          f'i={np.degrees(x[2]):5.2f}d  Omega={np.degrees(x[3]):5.1f}d  '
          f'omega={np.degrees(x[4]):5.1f}d  M0={np.degrees(x[5]):5.1f}d')
    for it in range(n_outer):
        # Anneal measurement noise from loose to tight
        frac = it / max(1, n_outer - 1)
        R_meas = R_meas_start * (R_meas_end / R_meas_start) ** frac
        x_new, P_new, res_b, res_a = ekf_pass(
            x, P, Q, R_meas, times, decs, verbose=verbose_inner)
        rms_b = float(np.sqrt(np.mean(res_b**2)))
        rms_a = float(np.sqrt(np.mean(res_a**2)))
        par = np.degrees(np.arcsin(R_EARTH_KM / x_new[0])) * 60.0
        sig_a = float(np.sqrt(max(P_new[0, 0], 0.0)))
        print(f'outer {it+1:2d}: R={np.sqrt(R_meas):4.1f}\'  '
              f'a={x_new[0]:8.1f}+/-{sig_a:6.0f}km  '
              f'e={x_new[1]:6.4f}  i={np.degrees(x_new[2]):5.2f}d  '
              f'pi={par:5.2f}\'  '
              f'RMS_pre={rms_b:5.2f}\'  RMS_post={rms_a:5.2f}\'')
        trace.append((x_new.copy(), rms_b, rms_a))
        x = x_new
        # Re-inflate P a bit for next outer pass so filter can still
        # move, but less than P0 -- this plays the role of a smoother.
        P = P_new + 0.1 * P0
        P = 0.5 * (P + P.T)

    x_final, rms_b, rms_a = trace[-1]
    a, e, i, Om, om, M0 = x_final
    n_rad_per_s = N_MEAN / 86400.0
    mu = n_rad_per_s**2 * a**3
    par = np.degrees(np.arcsin(R_EARTH_KM / a)) * 60.0

    print()
    print('===============  FINAL ESTIMATE  ===============')
    print(f'  a            = {a:8.1f} km    ({a/R_EARTH_KM:5.2f} R_earth)')
    print(f'  e            = {e:8.5f}')
    print(f'  i            = {np.degrees(i):7.3f} deg')
    print(f'  Omega        = {np.degrees(Om) % 360:7.3f} deg')
    print(f'  omega        = {np.degrees(om) % 360:7.3f} deg')
    print(f'  M0  @ epoch  = {np.degrees(M0) % 360:7.3f} deg')
    print(f'  horiz pi     = {par:6.2f} arcmin  (= arcsin(R_earth / a))')
    print(f'  mu = n^2 a^3 = {mu:8.0f} km^3/s^2')
    print(f'  RMS (pred)   = {rms_b:5.2f} arcmin')
    print(f'  RMS (post)   = {rms_a:5.2f} arcmin')
    print()
    print('Reference (not used in fit):')
    print('  a_modern_mean    = 384400 km   (60.27 R_earth)')
    print('  horiz pi modern  =  57.04 arcmin')
    print('  GM_earth+GM_moon =  403503 km^3/s^2')
    return trace


if __name__ == '__main__':
    import sys
    n_outer = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    run(n_outer=n_outer)
