"""
Modern topocentric ephemeris reproducing Tycho Brahe's
21 January 1586 (Julian) Moon observations at Uraniborg (Hven, Denmark).

Source: Tycho Brahe, *Loca Lunae*, 1586, p. 30.

Observation details
-------------------
  Date:                 21 Jan 1586 Julian  =  31 Jan 1586 Gregorian
  Clock times (H.M):    astronomical day counted from **noon** of Jan 21 Julian
                        H. 9      --> civil 21:00 Uraniborg local mean time
  Location:             Uraniborg, lat 55d54'N, lon 12d42'E

Quantities measured by Tycho:
  * Declination of Moon center
  * "Distantia aequatorea" Moon-Jupiter  = difference in RIGHT ASCENSION
    (29 deg 4', after applying Moon semi-diameter 16')
  * "Distantia aequatorea" Moon-Aldebaran = difference in RIGHT ASCENSION
    (25 deg 24.5')
  * Right ascensions of Moon, Jupiter, Aldebaran at H. 9 completa

Bugs fixed vs the original tycho.py:
  1. Calendar:  Jan 21 Julian = Jan 31 Gregorian (10-day shift).
  2. Clock:     hours are astronomical-day (counted from local apparent noon
                of Jan 21 Julian), not civil morning of Jan 21.
  3. Longitude: Uraniborg lon = 12d42' E -> UT = local - 50.8 min.
  4. Aldebaran: J2000 RA is 4h 35m 55s, not 4h 00m.
  5. "Distantia aequatorea" is Delta-RA, not great-circle separation.
  6. Comparisons are in **equinox of date (1586)**, not ICRS/J2000, because
     Tycho's RAs are naturally in the equinox of his own epoch.
  7. Ephemeris: long-span DE441 kernel covering years 13200 BC - 1969 AD
     (the default de440s.bsp / de441s.bsp only span 1849-2150).
"""

import warnings
import erfa
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import (
    solar_system_ephemeris, EarthLocation, SkyCoord, get_body,
    PrecessedGeocentric,
)
from astropy.coordinates.baseframe import NonRotationTransformationWarning
from astropy.utils.iers import conf as iers_conf
from astropy.utils.exceptions import AstropyWarning
import pandas as pd

# Silence harmless astropy/ERFA warnings for this 1586-CE calculation:
#   * ERFA 'dubious year' (date is far outside its calibrated window)
#   * IERS polar-motion fallback to 50-yr mean (sub-arcsec effect)
#   * NonRotationTransformationWarning (ICRS->GCRS for separation)
warnings.filterwarnings('ignore', category=erfa.ErfaWarning)
warnings.filterwarnings('ignore', category=NonRotationTransformationWarning)
warnings.filterwarnings('ignore', category=AstropyWarning,
                        message='.*polar motions.*')
iers_conf.auto_download = False

# ----------------------------------------------------------------------------
# Observation setup
# ----------------------------------------------------------------------------
URANIBORG_LAT_DEG = 55 + 54/60
URANIBORG_LON_DEG = 12 + 42/60
loc = EarthLocation(lat=URANIBORG_LAT_DEG * u.deg,
                    lon=URANIBORG_LON_DEG * u.deg,
                    height=40 * u.m)
LOCAL_OFFSET = TimeDelta(URANIBORG_LON_DEG / 15.0 * 3600.0, format='sec')

# Astronomical day for 21 Jan 1586 Julian starts at local noon on
# 31 Jan 1586 Gregorian.
NOON_LOCAL_JAN21_JULIAN = Time('1586-01-31 12:00:00', scale='ut1') - LOCAL_OFFSET


def astro_hour(h: int, m: float) -> Time:
    """'H. h M. m' (astronomical reckoning) -> UT1."""
    return NOON_LOCAL_JAN21_JULIAN + TimeDelta((h * 60 + m) * 60.0, format='sec')


t_decl      = astro_hour(9, 0)       # "H. 9 completa"
t_moon_jup  = astro_hour(8, 55.5)    # "H. 8 M. 55 1/2"
t_moon_alde = astro_hour(9, 2.5)     # "H. 9 M. 2 1/2"

# ----------------------------------------------------------------------------
# Ephemeris: long-span DE441
# ----------------------------------------------------------------------------
DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)

# ICRS / J2000 apparent topocentric positions
moon_d  = get_body('moon',    t_decl,     location=loc)
moon_j  = get_body('moon',    t_moon_jup, location=loc)
jupiter = get_body('jupiter', t_moon_jup, location=loc)
moon_a  = get_body('moon',    t_moon_alde, location=loc)
# Aldebaran, J2000 ICRS (Hipparcos)
aldebaran = SkyCoord(ra='04h35m55.24s', dec='+16d30m33.5s', frame='icrs')

# Precess to equinox of date for direct comparison with Tycho's numbers
def of_date(sc, t):
    return sc.transform_to(PrecessedGeocentric(equinox=t, obstime=t))

moon_d_d  = of_date(moon_d,    t_decl)
moon_j_d  = of_date(moon_j,    t_moon_jup)
jup_d     = of_date(jupiter,   t_moon_jup)
moon_a_d  = of_date(moon_a,    t_moon_alde)
ald_d     = of_date(aldebaran, t_moon_alde)

# ----------------------------------------------------------------------------
# Observables (NB: "Distantia aequatorea" = RA difference, not separation)
# ----------------------------------------------------------------------------
def wrap180(deg):
    return (deg + 180) % 360 - 180

dec_moon       = moon_d_d.dec.deg
sep_moon_jup   = moon_j.separation(jupiter).deg
sep_moon_alde  = moon_a.separation(aldebaran).deg

# Tycho's measured values
tycho = {
    'Declination of Moon center':         18 + 37/60,    # +18d 37'
    'Moon - Jupiter   (angular dist.)':   29 + 4/60,     # 29d 4'
    'Moon - Aldebaran (angular dist.)':   25 + 24.5/60,  # 25d 24.5'
}
mine = {
    'Declination of Moon center':         dec_moon,
    'Moon - Jupiter   (angular dist.)':   sep_moon_jup,
    'Moon - Aldebaran (angular dist.)':   sep_moon_alde,
}

rows = []
for k in tycho:
    err_arcmin = (tycho[k] - mine[k]) * 60
    rows.append((k + ' (deg)', round(tycho[k], 4), round(mine[k], 4), round(err_arcmin, 2)))
print(pd.DataFrame(rows, columns=['Quantity', 'Tycho', 'DE441 (of date)',
                                  'Error Tycho-Mine (arcmin)']).to_string(index=False))

# Also the three RAs Tycho lists
print()
print('Right ascensions at H. 9 completa (equinox of date):')
def fmt_dms(deg):
    d = int(deg); m = (deg - d) * 60
    return f"{d:3d} deg {m:5.1f} min"
m_RA_at_9 = of_date(get_body('moon', t_decl, location=loc), t_decl).ra.deg
print(f'  Moon       Tycho  88 deg 29.0 min   DE441 {fmt_dms(m_RA_at_9)}   (diff {( 88 + 29/60 - m_RA_at_9)*60:+6.2f}\')')
print(f'  Jupiter    Tycho  59 deg 25.0 min   DE441 {fmt_dms(jup_d.ra.deg)}   (diff {( 59 + 25/60 - jup_d.ra.deg)*60:+6.2f}\')')
print(f'  Aldebaran  Tycho  63 deg  5.0 min   DE441 {fmt_dms(ald_d.ra.deg)}   (diff {( 63 +  5/60 - ald_d.ra.deg)*60:+6.2f}\')')
