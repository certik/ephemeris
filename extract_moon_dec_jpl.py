"""
Generate moon_dec_observations_jpl.dat: identical to
moon_dec_observations.dat but with the declination column replaced by the
exact JPL DE441 apparent topocentric declination at Uraniborg at Tycho's
timestamps.

The goal is to isolate whether the residual ~arcmin scatter the EKF sees
on the real Tycho data is due to Tycho's instrumental noise or to the
physics not captured by our 2-body Kepler model.  Feeding perfect JPL
declinations into ekf_moon_orbit.py and watching what it converges to
answers that question.

Usage:
    uv run extract_moon_dec_jpl.py
    uv run ekf_moon_orbit.py      # uses moon_dec_observations.dat
    # then edit load path or call the helper below:
    uv run ekf_moon_orbit.py moon_dec_observations_jpl.dat
"""
import warnings
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import (
    solar_system_ephemeris, EarthLocation, SkyCoord, get_body,
    PrecessedGeocentric, AltAz,
)
from astropy.utils.iers import conf as iers_conf
from astropy.utils.exceptions import AstropyWarning
import erfa

warnings.filterwarnings('ignore', category=erfa.ErfaWarning)
warnings.filterwarnings('ignore', category=AstropyWarning)
iers_conf.auto_download = False

DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)

URANIBORG_LAT_DEG = 55 + 54/60
URANIBORG_LON_DEG = 12 + 42/60
loc = EarthLocation(lat=URANIBORG_LAT_DEG*u.deg,
                    lon=URANIBORG_LON_DEG*u.deg,
                    height=40*u.m)


def _bennett_refraction_deg(alt_deg: float) -> float:
    if alt_deg <= -1.0:
        return 0.0
    return (1.0/np.tan(np.radians(alt_deg + 7.31/(alt_deg + 4.4))))/60.0


def apparent_moon_dec_deg(t: Time) -> float:
    """Topocentric of-date Moon declination with Bennett refraction.
    Matches tycho*.py's apparent_of_date()."""
    mo = get_body('moon', t, location=loc)
    aa = mo.transform_to(AltAz(obstime=t, location=loc))
    R = _bennett_refraction_deg(float(aa.alt.deg))
    app = SkyCoord(alt=(aa.alt.deg + R)*u.deg, az=aa.az.deg*u.deg,
                   frame=AltAz(obstime=t, location=loc))
    return float(app.transform_to(
        PrecessedGeocentric(equinox=t, obstime=t)).dec.deg)


def main():
    src = 'moon_dec_observations.dat'
    dst = 'moon_dec_observations_jpl.dat'
    rows = []
    header_lines = []
    with open(src) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.rstrip('\n'))
                continue
            if not line.strip():
                continue
            parts = line.split()
            iso = ' '.join(parts[:2])
            source = parts[2]
            tag = parts[3]
            jd_tt = float(parts[5])
            rows.append((iso, source, tag, jd_tt))

    print(f'Computing JPL apparent dec for {len(rows)} observations ...')
    out_rows = []
    for i, (iso, source, tag, jd_tt) in enumerate(rows):
        t = Time(jd_tt, format='jd', scale='tt')
        dec = apparent_moon_dec_deg(t)
        out_rows.append((iso, source, tag, dec, jd_tt))
        if (i+1) % 10 == 0 or i == len(rows)-1:
            print(f'  {i+1:3d}/{len(rows)}  {iso}  {source}  dec={dec:+9.5f}')

    with open(dst, 'w') as f:
        f.write('# Moon CENTER apparent topocentric declination observations\n')
        f.write('# IDENTICAL timestamps to moon_dec_observations.dat, but the\n')
        f.write('# declination column is replaced by the exact JPL DE441\n')
        f.write('# topocentric apparent declination at Uraniborg (Bennett\n')
        f.write('# refraction applied, matching tycho*.py apparent_of_date()).\n')
        f.write(f'# {len(out_rows)} observations.\n')
        f.write('# Columns: iso_utc  source  tag    dec_center_deg  jd_tt\n')
        for iso, source, tag, dec, jd in out_rows:
            f.write(f'{iso:<26s}  {source:<8s}  {tag:<6s}  '
                    f'{dec:+11.6f}  {jd:.8f}\n')
    print(f'wrote {dst}')


if __name__ == '__main__':
    main()
