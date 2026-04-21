"""
Modern topocentric ephemeris for Tycho Brahe's 21 January 1586 (Julian)
Moon observation at Uraniborg.

Fixes applied vs. the original tycho.py:
  1. Date: 21 Jan 1586 Julian = 31 Jan 1586 Gregorian (10-day calendar shift).
  2. Clock times are local mean solar time at Uraniborg (lon = 12°42' E);
     convert to UT by subtracting lon/15 hours (~50.8 min).
  3. Aldebaran J2000 RA is 4h35m55s, not 4h00m.

Ephemeris: astropy 'builtin' (ERFA analytical). Note this is extrapolated
outside 1900-2100, so planetary positions (e.g. Jupiter) at 1586 are only
accurate to ~degrees. Moon and Sun remain good to ~arcmin. For full accuracy
in the 16th century, load the long DE441 kernel (de441_part-1.bsp).
"""

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import (
    solar_system_ephemeris, EarthLocation, SkyCoord, get_body,
)
import pandas as pd

# Uraniborg observatory (Hven, Denmark)
URANIBORG_LAT_DEG = 55 + 54/60        # 55° 54' N
URANIBORG_LON_DEG = 12 + 42/60        # 12° 42' E
loc = EarthLocation(
    lat=URANIBORG_LAT_DEG * u.deg,
    lon=URANIBORG_LON_DEG * u.deg,
    height=40 * u.m,
)

# Local mean solar time offset from UT (east positive)
LOCAL_OFFSET = TimeDelta(URANIBORG_LON_DEG / 15.0 * 3600.0, format='sec')


def local_to_ut(date_greg: str, h: int, m: float):
    """Convert Uraniborg local mean solar time to UT1."""
    s = int(round((m - int(m)) * 60))
    t_local = Time(f"{date_greg} {h:02d}:{int(m):02d}:{s:02d}", scale='ut1')
    return t_local - LOCAL_OFFSET


# 21 Jan 1586 Julian calendar = 31 Jan 1586 Gregorian
DATE = "1586-01-31"
t_decl     = local_to_ut(DATE, 9, 0)      # ~H. 9 completa
t_moon_jup = local_to_ut(DATE, 8, 55.5)   # H. 8 M. 55½
t_moon_alde = local_to_ut(DATE, 9, 2.5)   # H. 9 M. 2½

# Ephemeris
solar_system_ephemeris.set('builtin')

moon_decl = get_body('moon',    t_decl,     location=loc)
moon_jup  = get_body('moon',    t_moon_jup, location=loc)
jupiter   = get_body('jupiter', t_moon_jup, location=loc)
moon_alde = get_body('moon',    t_moon_alde, location=loc)

# Aldebaran (α Tauri), J2000 ICRS (Hipparcos)
aldebaran = SkyCoord(ra='04h35m55.24s', dec='+16d30m33.5s', frame='icrs')

dec_modern_deg         = moon_decl.dec.deg
dist_moon_jup_deg      = moon_jup.separation(jupiter).deg
dist_moon_alde_deg     = moon_alde.separation(aldebaran).deg

tycho_values = {
    'Declination (Moon center)':    18 + 37/60,        # +18° 37'
    'Angular dist. Moon-Jupiter':   29 + 4/60,         # 29° 04'
    'Angular dist. Moon-Aldebaran': 25 + 24.5/60,      # 25° 24.5'
}
modern_values = {
    'Declination (Moon center)':    dec_modern_deg,
    'Angular dist. Moon-Jupiter':   dist_moon_jup_deg,
    'Angular dist. Moon-Aldebaran': dist_moon_alde_deg,
}

data = {
    'Quantity': [
        'Declination of Moon center (°)',
        'Angular distance Moon-Jupiter (°)',
        'Angular distance Moon-Aldebaran (°)',
    ],
    'Tycho Brahe':        [tycho_values[k]  for k in tycho_values],
    'Modern (builtin)':   [round(modern_values[k], 6) for k in modern_values],
    'Error (Tycho - Modern) (arcmin)': [
        round((tycho_values[k] - modern_values[k]) * 60, 1) for k in tycho_values
    ],
}
print(pd.DataFrame(data).to_string(index=False))
print("\nNote: 'builtin' ephemeris extrapolates outside 1900-2100;")
print("Jupiter at 1586 is only accurate to ~1 degree. Use DE441 long kernel")
print("for full ~arcmin-level accuracy on all three quantities.")
