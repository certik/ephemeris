"""
Keplerian fit of the Moon's orbit to Tycho Brahe's 1586 declination
observations at Uraniborg, using ONLY observed quantities.

What is used as input
---------------------
* ``moon_dec_observations.dat``: pair-averaged Moon CENTER declinations
  (no ephemeris, no cusp assumption -- see ``extract_moon_dec.py``).
* The sidereal period T_sid = 27.321661 days.  This has been known with
  high accuracy since antiquity from repeated observations of the Moon
  passing given reference stars -- it is not an ephemeris.
* The observer's latitude and longitude (Uraniborg).
* Standard Earth-orientation models (precession, nutation, sidereal time)
  that relate topocentric alt/az to equator-of-date RA/dec.  These depend
  only on Earth rotation, not on any solar-system ephemeris.

What is NOT used
----------------
* No JPL kernel (DE441 etc.).
* No prior knowledge of the Moon's distance, mass, or position.

State vector (6 free parameters at epoch 1586-06-15 TT)
-------------------------------------------------------
    a       semi-major axis, km   (constrained by diurnal parallax)
    e       eccentricity
    i       inclination to equator-of-date at epoch, rad
    Omega   RA of ascending node, rad
    omega   argument of periapsis, rad
    M0      mean anomaly at epoch, rad

The mean motion n = 2*pi / T_sid is fixed.  This replaces the usual
(mu = n^2 a^3) link -- neither mu nor a is needed to *propagate* the
angular part of a Keplerian orbit once n is known; a enters the fit only
through the geocentric-to-topocentric parallax.
"""
import warnings

import numpy as np
import astropy.units as u
from astropy.coordinates import (
    SkyCoord, EarthLocation, GCRS, AltAz, PrecessedGeocentric,
    CartesianRepresentation,
)
from astropy.time import Time
from astropy.utils.iers import conf as iers_conf
from scipy.optimize import least_squares, differential_evolution

warnings.filterwarnings('ignore')
iers_conf.auto_download = False

T_SID_DAYS = 27.321661
N_MEAN = 2.0 * np.pi / T_SID_DAYS  # rad / day

URANIBORG = EarthLocation(lat=(55 + 54/60)*u.deg,
                          lon=(12 + 42/60)*u.deg,
                          height=40*u.m)
EPOCH = Time('1586-06-15 00:00:00', scale='tt')


# --- Kepler propagator ----------------------------------------------------

def solve_kepler(M, e, tol=1e-12, maxiter=80):
    M = np.atleast_1d(M).astype(float)
    E = np.where(e < 0.8, M, np.pi * np.ones_like(M))
    for _ in range(maxiter):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def kepler_xyz(elements, times_jdtt: np.ndarray) -> np.ndarray:
    """Return geocentric Moon position (km) in inertial equatorial frame
    (equator-of-date at EPOCH) for an array of TT Julian dates."""
    a, e, inc, Om, om, M0 = elements
    dt_days = np.asarray(times_jdtt) - EPOCH.tt.jd
    M = M0 + N_MEAN * dt_days
    E = solve_kepler(M, e)
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2),
                      np.sqrt(1-e)*np.cos(E/2))
    r = a*(1 - e*np.cos(E))
    xp = r*np.cos(nu); yp = r*np.sin(nu)
    cO, sO = np.cos(Om), np.sin(Om)
    ci, si = np.cos(inc), np.sin(inc)
    co, so = np.cos(om), np.sin(om)
    R11 = cO*co - sO*so*ci;  R12 = -cO*so - sO*co*ci
    R21 = sO*co + cO*so*ci;  R22 = -sO*so + cO*co*ci
    R31 = so*si;             R32 = co*si
    return np.array([R11*xp + R12*yp,
                     R21*xp + R22*yp,
                     R31*xp + R32*yp])   # shape (3, N)


# --- Observation model -----------------------------------------------------

def _bennett_deg(alt_deg):
    a = np.atleast_1d(alt_deg)
    out = np.where(a > -1.0,
                   (1.0/np.tan(np.radians(a + 7.31/(a + 4.4))))/60.0,
                   0.0)
    return out


def predicted_dec_array(elements, times: Time) -> np.ndarray:
    """Batched: return predicted apparent topocentric declination in deg."""
    xyz = kepler_xyz(elements, times.tt.jd)      # (3, N) km, equator@EPOCH
    # Treat this as GCRS-like inertial; the frame drift from EPOCH to each
    # observation (precession) is tiny (arcsec/year) over a few months, and
    # is the same precession astropy will apply to reach equator-of-date.
    # By using GCRS(obstime=t) astropy handles the frame transform internally.
    cart = CartesianRepresentation(xyz[0]*u.km, xyz[1]*u.km, xyz[2]*u.km)
    sc = SkyCoord(cart, frame=GCRS(obstime=times))
    aa = sc.transform_to(AltAz(obstime=times, location=URANIBORG))
    R = _bennett_deg(aa.alt.deg)
    app = SkyCoord(alt=(aa.alt.deg + R)*u.deg, az=aa.az.deg*u.deg,
                   frame=AltAz(obstime=times, location=URANIBORG))
    return app.transform_to(PrecessedGeocentric(
        equinox=times, obstime=times)).dec.deg


# --- Observations ---------------------------------------------------------

def load_observations(path='moon_dec_observations.dat'):
    jds, decs, sources = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            sources.append(parts[2])
            decs.append(float(parts[4]))
            jds.append(float(parts[5]))
    return (np.array(jds), np.array(decs),
            Time(np.array(jds), format='jd', scale='tt'), sources)


# --- Residuals ------------------------------------------------------------

def residuals(x, times, decs):
    try:
        pred = predicted_dec_array(tuple(x), times)
    except Exception:
        return np.full_like(decs, 1e4)
    return (pred - decs) * 60.0  # arcmin


# --- Global bootstrap -----------------------------------------------------

def global_bootstrap(times, decs, seed=0):
    """Find a reasonable starting point using observational bounds only."""
    # Loose bounds:
    #   a: 200_000 - 800_000 km (covers everything from sub-geostationary to
    #      well beyond any naked-eye hypothesis)
    #   e: 0 - 0.25
    #   i: 5 - 35 deg
    #   Om, om, M0: 0 - 2*pi
    bounds = [(200_000.0, 800_000.0),
              (0.0, 0.25),
              (np.radians(5.0), np.radians(35.0)),
              (0.0, 2*np.pi),
              (0.0, 2*np.pi),
              (0.0, 2*np.pi)]

    def cost(x):
        r = residuals(x, times, decs)
        return np.sqrt(np.mean(r*r))

    print('Global search (differential_evolution) -- this takes a minute ...',
          flush=True)
    result = differential_evolution(
        cost, bounds,
        seed=seed, tol=1e-4, popsize=25, mutation=(0.5, 1.3),
        recombination=0.7, maxiter=200, polish=False, workers=1)
    print(f'  DE best cost: {result.fun:.2f}\' at {result.x}')
    return result.x


# --- Main -----------------------------------------------------------------

def main():
    jds, decs, times, sources = load_observations()
    print(f'Loaded {len(decs)} Moon CENTER declination observations.')
    print(f'Sidereal period fixed: T = {T_SID_DAYS} days '
          f'(n = 2*pi/T = {N_MEAN:.6f} rad/day).\n')

    x0 = global_bootstrap(times, decs)
    labels = ['a(km)', 'e', 'i(deg)', 'Om(deg)', 'om(deg)', 'M0(deg)']

    def pretty(x):
        return [x[0], x[1], np.degrees(x[2]), np.degrees(x[3]),
                np.degrees(x[4]), np.degrees(x[5]) % 360]

    def show(x, tag):
        print(tag)
        for lab, v in zip(labels, pretty(x)):
            print(f'   {lab:>8s} = {v:12.6f}')

    # ---- Fit A: all 6 parameters free (a, e, i, Om, om, M0) ----
    scale = np.array([1000.0, 0.01, 0.01, 0.01, 0.01, 0.01])

    def scaled(y, x_anchor):
        return residuals(x_anchor + y*scale, times, decs)

    print('\n' + '='*70)
    print('FIT A: all 6 Keplerian elements free')
    print('='*70)
    result_A = least_squares(lambda y: scaled(y, x0), x0=np.zeros(6),
                             method='lm', xtol=1e-10, ftol=1e-10,
                             max_nfev=400)
    x_A = x0 + result_A.x * scale
    show(x_A, '\nFitted elements:')
    rA = result_A.fun
    print(f'\nRMS = {np.sqrt(np.mean(rA**2)):.2f}\',  '
          f'max |residual| = {np.max(np.abs(rA)):.2f}\'')
    muA = (N_MEAN/86400.0)**2 * x_A[0]**3
    print(f'Derived GM (Earth+Moon) = n^2 a^3 = {muA:.1f} km^3/s^2  '
          f'(modern: 403503)')

    # ---- Fit B: a fixed at Tycho's measured value (60 Earth radii),
    #              which he himself derived from his own lunar-parallax
    #              observations (not from any ephemeris) ----
    a_tycho = 60.0 * 6378.137   # 60 Earth radii -- Tycho's own value
    print('\n' + '='*70)
    print(f'FIT B: a = {a_tycho:.0f} km fixed (Tycho\'s 60 Earth radii),'
          ' 5 params free')
    print('='*70)
    x0B = np.concatenate(([a_tycho], x_A[1:]))
    scale5 = scale[1:]

    def scaledB(y):
        x = np.concatenate(([a_tycho], x0B[1:] + y*scale5))
        return residuals(x, times, decs)

    result_B = least_squares(scaledB, x0=np.zeros(5),
                             method='lm', xtol=1e-10, ftol=1e-10,
                             max_nfev=400)
    x_B = np.concatenate(([a_tycho], x0B[1:] + result_B.x * scale5))
    show(x_B, '\nFitted elements (a frozen):')
    rB = result_B.fun
    print(f'\nRMS = {np.sqrt(np.mean(rB**2)):.2f}\',  '
          f'max |residual| = {np.max(np.abs(rB)):.2f}\'')
    muB = (N_MEAN/86400.0)**2 * x_B[0]**3
    print(f'Derived GM (Earth+Moon) = n^2 a^3 = {muB:.1f} km^3/s^2  '
          f'(modern: 403503)')

    # Per-script residuals for the physically anchored fit
    print('\nPer-script residual statistics (Fit B):')
    src_arr = np.array(sources)
    for src in sorted(set(sources), key=lambda s: int(s.replace('tycho',''))):
        mask = src_arr == src
        chunk = rB[mask]
        print(f'   {src:<8s} n={mask.sum():3d}  '
              f'mean={chunk.mean():+5.2f}\'  '
              f'rms={np.sqrt((chunk**2).mean()):5.2f}\'')

    print('\n' + '='*70)
    print('INTERPRETATION')
    print('='*70)
    print("""
Fit A lets the semi-major axis a absorb unmodelled lunar perturbations
(evection ~78', variation ~40', parallactic inequality ~2') by
adjusting the diurnal-parallax amplitude -- it minimizes RMS at the
cost of an unphysical distance.

Fit B holds a at Tycho's own parallax-derived value.  The residual
~11' floor then reflects the genuine limit of a single-epoch two-body
Keplerian model against 10 months of real-Moon motion, not
measurement error.  Going lower requires adding evection / variation
terms (the "second lunar anomaly" discovered by Ptolemy + Tycho).
""")


if __name__ == '__main__':
    main()
