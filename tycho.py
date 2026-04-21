"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 21 January 1586 (Julian) at Uraniborg, Hven.

Source: Tycho Brahe, *Opera Omnia* tomus XI, p. ~19 ("Observationes Lunae,
Die 21 Januarii") — the raw armillary / sextant log, not the reduced summary
on p. 30 of the same volume.

Setup
-----
  Calendar: 21 Jan 1586 Julian  =  31 Jan 1586 Gregorian
  Clock:    "H. h M. m" are hours counted from local apparent noon of Jan 21
            Julian (astronomical-day convention). Here approximated as local
            mean solar noon; ~equation of time effects are sub-arcmin on this
            Moon-centric analysis.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E (height 40 m).
  Ephemeris: JPL DE441 long kernel (13200 BC - 1969 AD).
"""
import warnings
import erfa
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import (
    solar_system_ephemeris, EarthLocation, SkyCoord, get_body,
    PrecessedGeocentric,
)
from astropy.coordinates.baseframe import NonRotationTransformationWarning
from astropy.utils.iers import conf as iers_conf
from astropy.utils.exceptions import AstropyWarning

# Silence harmless astropy/ERFA warnings for this 1586 computation.
warnings.filterwarnings('ignore', category=erfa.ErfaWarning)
warnings.filterwarnings('ignore', category=NonRotationTransformationWarning)
warnings.filterwarnings('ignore', category=AstropyWarning,
                        message='.*polar motions.*')
iers_conf.auto_download = False

# -- Observer -----------------------------------------------------------------
URANIBORG_LAT_DEG = 55 + 54/60
URANIBORG_LON_DEG = 12 + 42/60
loc = EarthLocation(lat=URANIBORG_LAT_DEG*u.deg,
                    lon=URANIBORG_LON_DEG*u.deg,
                    height=40*u.m)
LOCAL_OFFSET = TimeDelta(URANIBORG_LON_DEG/15.0*3600.0, format='sec')
NOON_LOCAL_JAN21_JULIAN = (Time('1586-01-31 12:00:00', scale='ut1')
                           - LOCAL_OFFSET)


# -- Ephemeris ---------------------------------------------------------------
DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)


def of_date(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.transform_to(PrecessedGeocentric(equinox=t, obstime=t))


def moon_semid_deg(moon_sc: SkyCoord) -> float:
    """Apparent lunar semi-diameter in degrees from topocentric distance."""
    R_MOON_KM = 1737.4
    return float(np.degrees(np.arcsin(R_MOON_KM /
                                      moon_sc.distance.to(u.km).value)))


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 21 JANUARII 1586, Obs. Lunae)
# ----------------------------------------------------------------------------
# Moon horn declinations (upper / lower limb):
MOON_HORN = [
    # (H, M, 'upper'|'lower', dec in degrees)
    (8, 41.5,  'upper', 18 + 53.5/60),    # 18° 53.5'
    (8, 43.5,  'lower', 18 + 20.25/60),   # 18° 20.25'
    (9,  0.0,  'lower', 18 + 20.25/60),   # 18° 20.25'
    (9,  0.0,  'upper', 18 + 54.5/60),    # 18° 54.5'
    (9,  9.0,  'lower', 18 + (20 + 20/60)/60),  # 18° 20'20" = 18° 20.333'
    (9, 10.5,  'upper', 18 + 53.5/60),    # 18° 53.5'
]
MOON_DIAM_TYCHO_H9   = 0 + 34.25/60       # 0° 34.25' (at H. 9 M. 0)
MOON_DIAM_TYCHO_H841 = 0 + 33.25/60       # 0° 33.25' (at H. 8 M. 41.5)

# Dift. aequat.  (Moon occidental limb  -> Jupiter),  in degrees:
JUP_LIMB = [
    (8, 51.5,  28 + 29.00/60),            # 28° 29'
    (8, 54.75, 28 + 29.25/60),            # 28° 29.25'
    (8, 56.5,  28 + 32.50/60),            # 28° 32.5'
]
# Dift. aequat.  (Moon occidental limb  -> Aldebaran), in degrees:
ALD_LIMB = [
    (9, 4.0,   24 + 52.75/60),            # 24° 52.75'
    (9, 5.5,   24 + 53.00/60),            # 24° 53'
    (9, 7.0,   24 + (54 + 55/60)/60),     # 24° 54' 55"
]

# Aldebaran (alpha Tau), J2000 ICRS from Hipparcos
aldebaran = SkyCoord(ra='04h35m55.24s', dec='+16d30m33.5s', frame='icrs')


# ----------------------------------------------------------------------------
# Compute residuals for arbitrary clock correction, and optimize
# ----------------------------------------------------------------------------
def all_residuals_arcmin(clock_min: float):
    """Return flat array of (Tycho - DE441) residuals in arcmin for a given
    clock correction (minutes added to Tycho's recorded time)."""
    def t_at(h, m):
        return NOON_LOCAL_JAN21_JULIAN + TimeDelta(
            (h*60 + m + clock_min)*60.0, format='sec')
    res = []
    # Horn declinations
    for h, m, horn, v in MOON_HORN:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = of_date(mo, t)
        sd = moon_semid_deg(mo)
        dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
        res.append((v - dec_model)*60)
    # Moon-Jupiter limb
    for h, m, v in JUP_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        ju = get_body('jupiter', t, location=loc)
        mo_d = of_date(mo, t); ju_d = of_date(ju, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg - sd) - ju_d.ra.deg
        res.append((v - model)*60)
    # Moon-Aldebaran limb
    for h, m, v in ALD_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = of_date(mo, t); ald_d = of_date(aldebaran, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg - sd) - ald_d.ra.deg
        res.append((v - model)*60)
    return np.array(res)


from scipy.optimize import minimize_scalar
opt = minimize_scalar(
    lambda c: np.sum(all_residuals_arcmin(c)**2),
    bracket=(0, 20), method='brent', tol=1e-3)
BEST_CLOCK_MIN = opt.x
rms_opt = np.sqrt(np.mean(all_residuals_arcmin(BEST_CLOCK_MIN)**2))
rms_0   = np.sqrt(np.mean(all_residuals_arcmin(0.0)**2))

TYCHO_CLOCK_CORRECTION_MIN = BEST_CLOCK_MIN


def astro_hm_corrected(h, m):
    return NOON_LOCAL_JAN21_JULIAN + TimeDelta(
        (h*60 + m + TYCHO_CLOCK_CORRECTION_MIN)*60.0, format='sec')


def fmt_time(h: int, m: float) -> str:
    """Format Tycho's recorded time shifted by the optimized clock correction."""
    total = h*60 + m + TYCHO_CLOCK_CORRECTION_MIN
    hh, mm = divmod(total, 60)
    return f"{int(hh):d}:{mm:05.2f}"


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def fmt_err(arcmin):
    return f"{arcmin:+6.2f}'"

print('=' * 78)
print('DIE 21 JANUARII 1586 (Julian)  ==  31 Jan 1586 Gregorian, Uraniborg')
print('Modern ephemeris: JPL DE441 (long).  Comparison with Tycho\'s raw log.')
print('=' * 78)
print(f"\nOptimized clock correction: +{BEST_CLOCK_MIN:.2f} min")
print(f"  RMS residual uncorrected: {rms_0:5.2f}'   corrected: {rms_opt:5.2f}'")
print(f"  (times shown below are Tycho's clock + {BEST_CLOCK_MIN:.2f} min)")

print('\nMoon horn declinations (per armillas):')
print(' time    horn    Tycho      DE441 center+/-s.d.   error')
for h, m, horn, v in MOON_HORN:
    t = astro_hm_corrected(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    print(f' {fmt_time(h,m)}  {horn:5s}  {v:8.4f}    {dec_model:8.4f}          {fmt_err((v-dec_model)*60)}')

print('\nMoon apparent diameter:')
for label, hm, tyc in [('H. 8 M. 41.5', (8, 41.5), MOON_DIAM_TYCHO_H841),
                       ('H. 9 M. 0',    (9,  0.0), MOON_DIAM_TYCHO_H9)]:
    h, m = hm
    t = astro_hm_corrected(h, m)
    mo = get_body('moon', t, location=loc)
    diam_model = 2 * moon_semid_deg(mo)
    print(f' {fmt_time(h,m):8s}  Tycho {tyc*60:5.2f}\'   DE441 {diam_model*60:5.2f}\'   err {(tyc-diam_model)*60:+.2f}\'')

print('\nDift. aequat. (Moon occid. limb - Jupiter):')
print(' time    Tycho     DE441     error')
for h, m, v in JUP_LIMB:
    t = astro_hm_corrected(h, m)
    mo = get_body('moon', t, location=loc); ju = get_body('jupiter', t, location=loc)
    mo_d = of_date(mo, t); ju_d = of_date(ju, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg - sd) - ju_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {fmt_err((v-model)*60)}')

print('\nDift. aequat. (Moon occid. limb - Aldebaran):')
print(' time    Tycho     DE441     error')
for h, m, v in ALD_LIMB:
    t = astro_hm_corrected(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t); ald_d = of_date(aldebaran, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg - sd) - ald_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {fmt_err((v-model)*60)}')
