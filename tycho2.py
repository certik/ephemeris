"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 22 January 1586 (Julian) at Uraniborg, Hven.

Source: Tycho Brahe, *Opera Omnia* ("Observationes Lunae, Die 22 Januarii").
This is the raw armillary / sextant / quadrant log (not any reduced summary).

Setup
-----
  Calendar: 22 Jan 1586 Julian  =  1 Feb 1586 Gregorian
  Clock:    "H. h M. m" are hours counted from local (mean solar) noon of
            Jan 22 Julian, i.e. astronomical-day reckoning.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E (height 40 m).
  Ephemeris: JPL DE441 long kernel (13200 BC - 1969 AD).

Tycho notes: "Horologij circa horam 8 correcti hucusque nullus error"
  -- he re-synchronised the clock near H. 8; so expect the residual clock
  drift on this night to be small compared to the ~10.7 min seen on Jan 21.
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
from jplephem.spk import SPK
from scipy.optimize import brentq

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

# Local mean solar noon of 22 Jan 1586 Julian = 1 Feb 1586 Gregorian.
NOON_LOCAL = Time('1586-02-01 12:00:00', scale='ut1') - LOCAL_OFFSET

DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)
# astropy's get_body('sun', ...) misbehaves on the DE441 long kernel at 1586;
# we read the Sun directly via jplephem for the equation-of-time computation.
_SPK = SPK.open(DE441_LONG)


def sun_geocentric_icrs(t: Time) -> SkyCoord:
    jd = t.tt.jd
    pos = (_SPK[0, 10].compute(jd)
           - _SPK[0, 3].compute(jd)
           - _SPK[3, 399].compute(jd))
    r = float(np.linalg.norm(pos))
    ra = float(np.degrees(np.arctan2(pos[1], pos[0])) % 360.0)
    dec = float(np.degrees(np.arcsin(pos[2]/r)))
    return SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')


def equation_of_time_shift_min(date_str: str) -> float:
    """Minutes by which local apparent noon lags local mean noon at Uraniborg
    on the given Gregorian date (i.e. the offset to ADD to Tycho's
    apparent-noon-based hour count to get local mean solar time)."""
    def ha_deg(sec):
        t = (Time(date_str + ' 00:00:00', scale='ut1')
             + TimeDelta((12 - URANIBORG_LON_DEG/15.0)*3600 + sec, format='sec'))
        sun = sun_geocentric_icrs(t).transform_to(
            PrecessedGeocentric(equinox=t, obstime=t))
        lst = t.sidereal_time('apparent', longitude=loc.lon).deg
        return ((lst - sun.ra.deg + 540.0) % 360.0) - 180.0
    sec = brentq(ha_deg, -3600.0, 3600.0, xtol=1e-3)
    return sec / 60.0


def of_date(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.transform_to(PrecessedGeocentric(equinox=t, obstime=t))


def moon_semid_deg(moon_sc: SkyCoord) -> float:
    R_MOON_KM = 1737.4
    return float(np.degrees(np.arcsin(R_MOON_KM /
                                      moon_sc.distance.to(u.km).value)))


# Aldebaran (alpha Tau), J2000 ICRS from Hipparcos
aldebaran = SkyCoord(ra='04h35m55.24s', dec='+16d30m33.5s', frame='icrs')


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 22 IANUARIJ 1586)
# ----------------------------------------------------------------------------
# Moon horn declinations (per armillas):
MOON_HORN = [
    # (H, M, 'upper'|'lower', dec in degrees)
    (9, 55.5,  'lower', 17 +  9.0/60),    # 17 deg  9'
    (9, 57.0,  'upper', 17 + 45.0/60),    # 17 deg 45'
    (10, 10.5, 'lower', 17 +  8.5/60),    # 17 deg  8.5'
    (10, 13.5, 'upper', 17 + 44.0/60),    # 17 deg 44'
    (10, 31.5, 'lower', 17 +  5.0/60),    # 17 deg  5' (cornu)
    (10, 31.5, 'upper', 17 + 42.5/60),    # 17 deg 42.5'
]

# Moon apparent diameters (degrees), recorded alongside each dec pair:
MOON_DIAM = [
    (9,  57.0, 36.0/60),    # 0 deg 36'
    (10, 13.5, 35.5/60),    # 0 deg 35.5'
    (10, 31.5, 37.5/60),    # 0 deg 37.5'
]

# Moon-limb -> Aldebaran, "Dift. aequat." (= RA(moon occid.limb) - RA(Aldebaran)):
# Times given as H. M S (M here fractional minutes).
# From the two tables on the page (oculus Taurij = Aldebaran):
ALD_LIMB = [
    # (H, M, dRA_deg)
    (10,  0 + 40/60,  40 + 26.0/60),   # 10h 00m 40s
    (10,  2 + 10/60,  40 + 28.0/60),   # 10h 02m 10s
    (10,  3 + 40/60,  40 + 30.0/60),   # 10h 03m 40s
    (10,  6.0,        40 + 32.0/60),   # 10h 06m
    (10,  7.5,        40 + 33.0/60),   # 10h 07m 30s
    (10, 17.0,        40 + 35.5/60),   # 10h 17m
    (10, 18.5,        40 + 37.5/60),   # 10h 18m 30s
    (10, 20.0,        40 + 37.5/60),   # 10h 20m
    (10, 21.5,        40 + 38.5/60),   # 10h 21m 30s
]

# Moon-limb -> Jupiter, "Dift. aequat." (= RA(moon occid.limb) - RA(Jupiter)):
JUP_LIMB = [
    (10, 25.5, 44 + 20.5/60),    # 44 deg 20.5'
    (10, 27.0, 44 + 21.5/60),    # 44 deg 21.5'
    (10, 29.0, 44 + 23.0/60),    # 44 deg 23'
]

# Clock check: Aldebaran meridian transit recorded by Tycho.
ALDEBARAN_TRANSIT_HM = (7, 11.5)


# ----------------------------------------------------------------------------
# Residuals with EoT-only shift (NO fitted clock parameter)
# ----------------------------------------------------------------------------
# Tycho's H.0 convention = local APPARENT noon (sundial). His clocks run at
# mean rate (he re-synchronised them on stellar meridian transits throughout
# the night). The equation of time contributes a fixed offset between the
# apparent-noon-based hour count and the instant on a mean-solar clock.
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-02-01')
CLOCK_SHIFT_MIN = EOT_SHIFT_MIN   # no clock fit; only EoT is applied


def t_at(h, m):
    return NOON_LOCAL + TimeDelta((h*60 + m + CLOCK_SHIFT_MIN)*60.0, format='sec')


def all_residuals_arcmin():
    """Return flat array of (Tycho - DE441) residuals in arcmin."""
    res = []
    for h, m, horn, v in MOON_HORN:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = of_date(mo, t)
        sd = moon_semid_deg(mo)
        dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
        res.append((v - dec_model)*60)
    for h, m, v in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        diam_model = 2*moon_semid_deg(mo)
        res.append((v - diam_model)*60)
    for h, m, v in ALD_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = of_date(mo, t); ald_d = of_date(aldebaran, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg - sd) - ald_d.ra.deg
        res.append((v - model)*60)
    for h, m, v in JUP_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        ju = get_body('jupiter', t, location=loc)
        mo_d = of_date(mo, t); ju_d = of_date(ju, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg - sd) - ju_d.ra.deg
        res.append((v - model)*60)
    return np.array(res)


_r = all_residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


# ----------------------------------------------------------------------------
# Aldebaran meridian transit prediction (clock check)
# ----------------------------------------------------------------------------
def aldebaran_hour_angle_deg(clock_min_from_noon):
    """Aldebaran's hour angle at the time NOON_LOCAL + (offset_min) minutes."""
    t = NOON_LOCAL + TimeDelta(clock_min_from_noon*60, format='sec')
    lst = t.sidereal_time('apparent', longitude=loc.lon)
    ald_d = of_date(aldebaran, t)
    ha = (lst.deg - ald_d.ra.deg + 540.0) % 360.0 - 180.0
    return ha


transit_min = brentq(aldebaran_hour_angle_deg, 6*60, 9*60, xtol=1e-4)
transit_h = int(transit_min // 60)
transit_m = transit_min - transit_h*60
tycho_transit_min = ALDEBARAN_TRANSIT_HM[0]*60 + ALDEBARAN_TRANSIT_HM[1]


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 22 IANUARIJ 1586 (Julian)  ==  1 Feb 1586 Gregorian, Uraniborg')
print('Modern ephemeris: JPL DE441 (long).  Comparison with Tycho raw log.')
print('=' * 78)
print()
print('Aldebaran meridian transit (clock check):')
print(f'  Tycho recorded: H.{ALDEBARAN_TRANSIT_HM[0]} M.{ALDEBARAN_TRANSIT_HM[1]:.1f}'
      f'  (= {tycho_transit_min:.2f} min past local noon)')
print(f'  DE441 predicts: H.{transit_h} M.{transit_m:.2f}'
      f'  (= {transit_min:.2f} min past local noon)')
print(f'  Clock shift implied (Tycho - DE441): '
      f'{tycho_transit_min - transit_min:+.2f} min')
print()
print('Single fixed shift applied: equation of time (no clock fit).')
print(f'  EoT shift (apparent noon -> mean noon) = {EOT_SHIFT_MIN:+.2f} min')
print(f'  Aldebaran transit offset implied (Tycho - DE441): '
      f'{tycho_transit_min - transit_min:+.2f} min')
print(f'  Residual clock error (after EoT applied): '
      f'{(tycho_transit_min - transit_min) + EOT_SHIFT_MIN:+.2f} min')
print(f'  All-observation RMS residual: {rms:.2f}\'')

print('\nMoon horn declinations (per armillas):')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_HORN:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    print(f' {fmt_time(h,m)}  {horn:5s}  {v:8.4f}  {dec_model:8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    diam_model = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  Tycho {v*60:5.2f}\'   '
          f'DE441 {diam_model*60:5.2f}\'   err {(v-diam_model)*60:+.2f}\'')

print('\nDift. aequat. (Moon occid. limb - Aldebaran):')
print(' time     Tycho     DE441     error')
for h, m, v in ALD_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t); ald_d = of_date(aldebaran, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg - sd) - ald_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDift. aequat. (Moon occid. limb - Jupiter):')
print(' time     Tycho     DE441     error')
for h, m, v in JUP_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    ju = get_body('jupiter', t, location=loc)
    mo_d = of_date(mo, t); ju_d = of_date(ju, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg - sd) - ju_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

# Per-channel summary
r = _r
dec_r   = r[0:6]
diam_r  = r[6:9]
ald_r   = r[9:18]
jup_r   = r[18:21]
print('\nPer-channel residual summary (EoT shift applied, NO clock fit):')
for name, x in [('6 horn declinations', dec_r),
                ('3 Moon diameters   ', diam_r),
                ('9 Moon-Aldebaran dRA', ald_r),
                ('3 Moon-Jupiter  dRA ', jup_r)]:
    print(f'  {name}  mean = {x.mean():+.2f}\'   '
          f'scatter(std) = {x.std(ddof=1):.2f}\'')
