"""
Least-squares fit of a Keplerian orbit (Moon around Earth) to Tycho Brahe's
1586 declination observations at Uraniborg.

State vector (6 free parameters):
    a       semi-major axis, km
    e       eccentricity
    i       inclination to J2000 equator, rad
    Omega   RA of ascending node, rad
    omega   argument of periapsis, rad
    M0      mean anomaly at epoch (J2000 TT), rad

mu = GM_Earth + GM_Moon held fixed at 403503.2 km^3/s^2 (standard).

Predicted declination for each observation:
    Kepler-propagate the state to the observation epoch (J2000 equatorial
    Cartesian), then transform through astropy to apparent (equator-of-date)
    declination at Uraniborg, with refraction.  This automatically includes
    topocentric parallax, precession, and nutation.
"""
import warnings
from datetime import datetime

import numpy as np
import astropy.units as u
from astropy.coordinates import (
    SkyCoord, EarthLocation, GCRS, AltAz, PrecessedGeocentric,
    CartesianRepresentation,
)
from astropy.time import Time, TimeDelta
from astropy.utils.iers import conf as iers_conf
from scipy.optimize import least_squares

warnings.filterwarnings('ignore')
iers_conf.auto_download = False

MU_EARTH_MOON = 403503.2358         # km^3/s^2 (GM_Earth + GM_Moon)
URANIBORG = EarthLocation(lat=(55 + 54/60)*u.deg,
                          lon=(12 + 42/60)*u.deg,
                          height=40*u.m)
EPOCH = Time('1586-06-15 00:00:00', scale='tt')  # mid-1586, in DE441 part-1 range


# --- Kepler propagator ----------------------------------------------------

def solve_kepler(M, e, tol=1e-12, maxiter=60):
    M = np.atleast_1d(M).astype(float)
    E = np.where(e < 0.8, M, np.pi)
    for _ in range(maxiter):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def kepler_xyz_j2000(elements, t: Time) -> np.ndarray:
    """Return geocentric Moon position (km) in J2000 equatorial Cartesian."""
    a, e, inc, Om, om, M0 = elements
    dt = (t.tt.jd - EPOCH.tt.jd) * 86400.0  # seconds
    n = np.sqrt(MU_EARTH_MOON / a**3)
    M = M0 + n * dt
    E = solve_kepler(M, e)[0] if np.ndim(M) else solve_kepler(M, e)[0]
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
                     R31*xp + R32*yp])


# --- Observation model -----------------------------------------------------

def _bennett_deg(alt):
    if alt <= -1.0:
        return 0.0
    return (1.0/np.tan(np.radians(alt + 7.31/(alt + 4.4))))/60.0


def predicted_dec_deg(elements, t: Time) -> float:
    r = kepler_xyz_j2000(elements, t)
    cart = CartesianRepresentation(r[0]*u.km, r[1]*u.km, r[2]*u.km)
    sc = SkyCoord(cart, frame=GCRS(obstime=t))
    aa = sc.transform_to(AltAz(obstime=t, location=URANIBORG))
    R = _bennett_deg(float(aa.alt.deg))
    app = SkyCoord(alt=(aa.alt.deg + R)*u.deg, az=aa.az.deg*u.deg,
                   frame=AltAz(obstime=t, location=URANIBORG))
    return float(app.transform_to(
        PrecessedGeocentric(equinox=t, obstime=t)).dec.deg)


# --- Load observations -----------------------------------------------------

def load_observations(path='moon_dec_observations.dat'):
    times, decs, sources = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Format: "YYYY-MM-DD HH:MM:SS.SSS  source  horn  raw  center  jd"
            date, tod = line.split()[0], line.split()[1]
            parts = line.split()
            iso = f'{date} {tod}'
            src = parts[2]
            dec_center = float(parts[5])
            jd_tt = float(parts[6])
            times.append(jd_tt)
            decs.append(dec_center)
            sources.append(src)
    return np.array(times), np.array(decs), sources


# --- Residuals & fit -------------------------------------------------------

def residuals(x, times_jdtt, decs):
    elems = tuple(x)
    out = np.empty_like(decs)
    for k, (jd, d_obs) in enumerate(zip(times_jdtt, decs)):
        t = Time(jd, format='jd', scale='tt')
        out[k] = predicted_dec_deg(elems, t) - d_obs
    return out * 60.0  # arcmin


def initial_guess_from_de441():
    """Osculating Kepler elements derived from DE441 Moon state at J2000."""
    from jplephem.spk import SPK
    kern = SPK.open('/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp')
    jd = EPOCH.tt.jd
    # Earth-Moon barycenter -> geocentric Moon (segments 3->301 give EMB-Moon;
    # 3->399 gives EMB-Earth geocenter).  r(Moon wrt Earth) = r(3,301)-r(3,399)
    pos, vel = kern[3, 301].compute_and_differentiate(jd)
    pos_e, vel_e = kern[3, 399].compute_and_differentiate(jd)
    r = pos - pos_e     # km
    v = (vel - vel_e) / 86400.0  # km/s (compute_and_differentiate returns /day)
    mu = MU_EARTH_MOON
    h = np.cross(r, v)
    r_mag = np.linalg.norm(r); v_mag = np.linalg.norm(v); h_mag = np.linalg.norm(h)
    e_vec = np.cross(v, h)/mu - r/r_mag
    e = np.linalg.norm(e_vec)
    energy = v_mag**2/2 - mu/r_mag
    a = -mu / (2*energy)
    inc = np.arccos(h[2]/h_mag)
    n_vec = np.array([-h[1], h[0], 0.0])
    n_mag = np.linalg.norm(n_vec)
    Om = np.arccos(n_vec[0]/n_mag) if n_mag > 0 else 0.0
    if n_vec[1] < 0: Om = 2*np.pi - Om
    om = np.arccos(np.dot(n_vec, e_vec) / (n_mag*e))
    if e_vec[2] < 0: om = 2*np.pi - om
    nu = np.arccos(np.dot(e_vec, r) / (e*r_mag))
    if np.dot(r, v) < 0: nu = 2*np.pi - nu
    E = 2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2),
                     np.sqrt(1+e)*np.cos(nu/2))
    M0 = E - e*np.sin(E)
    return np.array([a, e, inc, Om, om, M0 % (2*np.pi)])


# --- Main ------------------------------------------------------------------

def main():
    times_jdtt, decs, sources = load_observations()
    x0 = initial_guess_from_de441()
    labels = ['a(km)', 'e', 'i(deg)', 'Om(deg)', 'om(deg)', 'M0(deg)']

    def pretty(x):
        return [x[0], x[1], np.degrees(x[2]), np.degrees(x[3]),
                np.degrees(x[4]), np.degrees(x[5]) % 360]

    def show(x, tag):
        vals = pretty(x)
        print(f'{tag}:')
        for lab, v in zip(labels, vals):
            print(f'   {lab:>8s} = {v:12.6f}')

    show(x0, 'Initial guess (DE441 osculating at 1586-06-15 TT)')
    r0 = residuals(x0, times_jdtt, decs)
    print(f'\ninitial-guess RMS dec residual = {np.sqrt(np.mean(r0**2)):.2f}\'')
    print(f'initial-guess max |residual|   = {np.max(np.abs(r0)):.2f}\'')

    print('\nFitting ...', flush=True)
    # Scale parameters so LS behaves well
    scale = np.array([1000.0, 0.01, 0.01, 0.01, 0.01, 0.01])

    def scaled_residuals(y):
        return residuals(x0 + y*scale, times_jdtt, decs)

    result = least_squares(scaled_residuals, x0=np.zeros(6),
                           method='lm', xtol=1e-10, ftol=1e-10, max_nfev=400)
    x_fit = x0 + result.x * scale

    show(x_fit, '\nFitted elements')
    rf = result.fun
    print(f'\nfitted RMS dec residual = {np.sqrt(np.mean(rf**2)):.2f}\'')
    print(f'fitted max |residual|   = {np.max(np.abs(rf)):.2f}\'')

    # Breakdown by source (tycho file)
    print('\nPer-script residual statistics (fitted):')
    src_arr = np.array(sources)
    for src in sorted(set(sources), key=lambda s: int(s.replace('tycho',''))):
        mask = src_arr == src
        chunk = rf[mask]
        print(f'   {src:<8s} n={mask.sum():3d}  '
              f'mean={chunk.mean():+5.2f}\'  rms={np.sqrt((chunk**2).mean()):5.2f}\'')


if __name__ == '__main__':
    main()
