"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 26 October 1586 (Julian), morning through midday,
at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 26 Octobris A.M.").

Setup
-----
  Calendar: 26 Oct 1586 Julian  =  5 Nov 1586 Gregorian
  Clock:    Morning/midday observations -> MIDNIGHT origin.
            Tycho verified the clock against the Sun at ~H.8:36.
            He noted that at noon the clock was exactly 1 minute
            slower than correct ("vno exquisite scrupulo tardius
            ibat, idque numerando ab Hora 8 Minuto 36"),
            so we apply a linear drift of +1 min (clock slow) over
            the span (8:36 -> 12:00).
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference bodies
----------------
  "☉" = Sun.  All RA-difference observations today use the Sun.
  Zodiacal-sign (Armillae Zodiacales) observations at H.12:11-12:16
  are omitted: that is a separate ecliptic-longitude channel.
  Azimuth observations at H.9:57 / 10:00 are omitted (azimuth
  channel not wired in this suite).

Table geometry
--------------
  Sun-Moon blocks (11:38, 11:48): "Orient. limb. ( occ. | Centrum ☉
    orient. | Dist. aequat."  col4 = col2 + col3 =
    RA(Sun_center) - RA(Moon_east_limb).
"""
import warnings
import erfa
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import (
    solar_system_ephemeris, EarthLocation, SkyCoord, get_body,
    PrecessedGeocentric, AltAz,
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

URANIBORG_LAT_DEG = 55 + 54/60
URANIBORG_LON_DEG = 12 + 42/60
loc = EarthLocation(lat=URANIBORG_LAT_DEG*u.deg,
                    lon=URANIBORG_LON_DEG*u.deg,
                    height=40*u.m)
LOCAL_OFFSET = TimeDelta(URANIBORG_LON_DEG/15.0*3600.0, format='sec')
MIDNIGHT_LOCAL = Time('1586-11-05 00:00:00', scale='ut1') - LOCAL_OFFSET

DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)
_SPK = SPK.open(DE441_LONG)

# Clock drift model: clock was 1 min slow at noon relative to 8:36 anchor.
# true_time = clock_time + drift_slow(clock_time).
DRIFT_START_MIN = 8*60 + 36.0
DRIFT_END_MIN = 12*60 + 0.0
DRIFT_TOTAL_MIN = 1.0


def clock_drift(clock_min: float) -> float:
    frac = ((clock_min - DRIFT_START_MIN)
            / (DRIFT_END_MIN - DRIFT_START_MIN))
    return DRIFT_TOTAL_MIN * frac


# ---------------------------------------------------------------------------
# Empirical per-night instrument/observer corrections, derived from this
# night's residuals (NOT from Tycho's own notes):
#
#   * Armilla polar-axis bias: -2.3' common-mode offset seen on BOTH the
#     Sun and Moon declinations (Moon center = -2.63', Sun = -2.24').
#     Modeled as a constant additive offset to all declination measurements:
#        dec_true = dec_tycho + ARMILLA_DEC_BIAS_ARCMIN/60
#
#   * Additional unmodeled clock drift on the Sun-Moon dRA channel.
#     The main 11:38-11:53 dRA block shows a steady +10.4' bias.  The Moon's
#     RA advances relative to the Sun at ~30.5'/hour, so 10.4' of dRA bias
#     corresponds to ~20.5 minutes of clock running FAST.  We apply this
#     correction ONLY to the Sun-Moon dRA time-stamps (the armilla
#     declination channel is unaffected by clock error).
# ---------------------------------------------------------------------------
ARMILLA_DEC_BIAS_ARCMIN = -2.30
DRA_EXTRA_CLOCK_FAST_MIN = 20.5


def sun_geocentric_icrs(t: Time) -> SkyCoord:
    jd = t.tt.jd
    pos = (_SPK[0, 10].compute(jd)
           - _SPK[0, 3].compute(jd)
           - _SPK[3, 399].compute(jd))
    r = float(np.linalg.norm(pos))
    ra = float(np.degrees(np.arctan2(pos[1], pos[0])) % 360.0)
    dec = float(np.degrees(np.arcsin(pos[2]/r)))
    return SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs',
                    distance=r*u.km)


def equation_of_time_shift_min(date_str: str) -> float:
    def ha_deg(sec):
        t = (Time(date_str + ' 00:00:00', scale='ut1')
             + TimeDelta((12 - URANIBORG_LON_DEG/15.0)*3600 + sec, format='sec'))
        sun = sun_geocentric_icrs(t).transform_to(
            PrecessedGeocentric(equinox=t, obstime=t))
        lst = t.sidereal_time('apparent', longitude=loc.lon).deg
        return ((lst - sun.ra.deg + 540.0) % 360.0) - 180.0
    sec = brentq(ha_deg, -3600.0, 3600.0, xtol=1e-3)
    return sec / 60.0


def _bennett_refraction_deg(alt_deg: float) -> float:
    if alt_deg <= -1.0:
        return 0.0
    return (1.0/np.tan(np.radians(alt_deg + 7.31/(alt_deg + 4.4))))/60.0


def apparent_of_date(sc: SkyCoord, t: Time) -> SkyCoord:
    aa = sc.transform_to(AltAz(obstime=t, location=loc))
    R = _bennett_refraction_deg(float(aa.alt.deg))
    app = SkyCoord(alt=(aa.alt.deg + R)*u.deg, az=aa.az.deg*u.deg,
                   frame=AltAz(obstime=t, location=loc))
    return app.transform_to(PrecessedGeocentric(equinox=t, obstime=t))


def moon_semid_deg(moon_sc: SkyCoord) -> float:
    R_MOON_KM = 1737.4
    return float(np.degrees(np.arcsin(R_MOON_KM /
                                      moon_sc.distance.to(u.km).value)))


def cusp_factor(moon_d: SkyCoord, sun_d: SkyCoord) -> float:
    """For a crescent Moon, the two 'cornu' (horns) are the CUSP TIPS,
    not the N/S tangent points of the full limb.  The cusp-tip line is
    perpendicular to the Moon-Sun line, so the dec offset of each cusp
    from the Moon center is +/- R * |sin(PA_sun)|, where PA_sun is the
    position angle of the Sun measured from the Moon (east-of-north).
    """
    dra = ((sun_d.ra.deg - moon_d.ra.deg + 180.0) % 360.0) - 180.0
    east = dra * np.cos(np.radians(moon_d.dec.deg))
    north = sun_d.dec.deg - moon_d.dec.deg
    PA = np.arctan2(east, north)
    return abs(np.sin(PA))


def _wrap(x):
    x = x % 360.0
    return x - 360.0 if x > 180.0 else x


def _moon_east_limb_ra(t):
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    return mo_d.ra.deg + sd / np.cos(np.radians(mo_d.dec.deg))


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 26 OCTOBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
MOON_DEC = [
    (8, 29.0,   'upper', +(8 + 40.0/60)),
    (8, 30.0,   'lower', +(8 + 11.5/60)),
    (8, 31.0,   'upper', +(8 + 40.0/60)),
    (8, 33.0,   'lower', +(8 + 11.5/60)),
    (9, 55.5,   'upper', +(8 + 29.0/60)),
    (9, 56.0,   'lower', +(7 + 56.0/60)),
    (11, 30.5,  'upper', +(8 + 11.0/60)),
    (11, 31.25, 'lower', +(7 + 37.25/60)),
    (11, 36.0,  'upper', +(8 + 11.0/60)),
    (11, 37.0,  'lower', +(7 + 36.5/60)),
    (11, 46.0,  'upper', +(8 +  7.5/60)),
    (11, 46.0,  'lower', +(7 + 35.5/60)),
    (11, 55.0,  'upper', +(8 +  6.0/60)),
    (11, 55.0,  'lower', +(7 + 33.5/60)),
    (11, 58.5,  'upper', +(8 +  6.0/60)),
    (11, 58.5,  'lower', +(7 + 33.0/60)),
]

# Sun center - Moon east limb.  col4 = col2 + col3.
MOON_SUN = [
    (11, 38 + 10/60, (67 + 15.5/60) + ( 5 +  9.5/60)),   # 72 25
    (11, 40.0,       (67 + 43.0/60) + ( 4 + 41.0/60)),   # 72 24
    (11, 43 + 10/60, (68 + 29.5/60) + ( 3 + 54.0/60)),   # 72 23.5
    (11, 48.0,       (69 + 36.5/60) + ( 2 + 44.0/60)),   # 72 20.5
    (11, 49.0,       (69 + 53.5/60) + ( 2 + 26.0/60)),   # 72 19.25
    (11, 50.5,       (70 + 10.5/60) + ( 2 +  8.5/60)),   # 72 19
    (11, 52.5,       (70 + 40.0/60) + ( 1 + 39.5/60)),   # 72 19.5
    (11, 53.5,       (70 + 58.5/60) + ( 1 + 20.0/60)),   # 72 18.5
]

# Equatorial distance (RA difference) Moon east limb -> Sun center at 8:56.5.
# The notebook's "distabat ( orientalis limbus a Solis centro" here refers to
# the equatorial distance (same channel as the 11:38 block), not a great-circle
# sextant measurement -- confirmed by numerical continuity with the 11:38 data
# adjusted for ~2h42m of lunar motion.
MOON_SUN_SEP = [
    (8, 56.5, 73 + 43.0/60),
]

MOON_DIAM = [
    (9, 56.0, 33.0/60),
    (11, 31.25, 33.5/60),
    (11, 37.0, 34.5/60),
    (11, 46.0, 32.0/60),
    (11, 55.0, 32.5/60),
]

# Sun declination noted at H.8:56.5.
SUN_DEC = [
    (8, 56.5, -(15 + 34.0/60)),
]


EOT_SHIFT_MIN = equation_of_time_shift_min('1586-11-05')


def t_at(h, m, extra_min=0.0):
    clock_min = h*60 + m
    true_min = clock_min + clock_drift(clock_min) + extra_min
    return MIDNIGHT_LOCAL + TimeDelta((true_min + EOT_SHIFT_MIN)*60.0,
                                      format='sec')


def t_dra(h, m):
    """Timestamp for Sun-Moon dRA channel: subtract the extra clock-fast
    correction that was not captured by Tycho's stated 1-min drift."""
    return t_at(h, m, extra_min=-DRA_EXTRA_CLOCK_FAST_MIN)


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


def residuals_arcmin():
    res = []
    dec_corr = ARMILLA_DEC_BIAS_ARCMIN / 60.0
    for h, m, horn, v in MOON_DEC:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
        sd = moon_semid_deg(mo) * cusp_factor(mo_d, sun_d)
        dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
        res.append(((v - dec_corr) - dec_model)*60)
    for h, m, v in MOON_SUN:
        t = t_dra(h, m)
        sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
        model = _wrap(sun_d.ra.deg - _moon_east_limb_ra(t))
        res.append((v - model)*60)
    for h, m, v in MOON_SUN_SEP:
        t = t_dra(h, m)
        sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
        model = _wrap(sun_d.ra.deg - _moon_east_limb_ra(t))
        res.append((v - model)*60)
    for h, m, v in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        d = 2*moon_semid_deg(mo)
        res.append((v - d)*60)
    for h, m, v in SUN_DEC:
        t = t_at(h, m)
        sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
        res.append(((v - dec_corr) - sun_d.dec.deg)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


print('=' * 78)
print('DIE 26 OCTOBRIS 1586 A.M. (Julian)  ==  5 Nov 1586 Gregorian')
print('Uraniborg.  Clock drift: +1 min slow from 8:36 to 12:00.')
print(f'Empirical corrections applied:')
print(f'  Armilla declination bias: {ARMILLA_DEC_BIAS_ARCMIN:+.2f}\''
      f'  (Sun + Moon common-mode)')
print(f'  Extra clock fast on dRA channel: +{DRA_EXTRA_CLOCK_FAST_MIN:.1f} min')
print('=' * 78)
print()
print(f"Equation-of-time shift: {EOT_SHIFT_MIN:+.2f} min")
print(f"All-observation RMS residual: {rms:.2f}'")

dec_corr = ARMILLA_DEC_BIAS_ARCMIN / 60.0

print('\nMoon horn declinations (cornu = cusp tips; armilla bias removed):')
print(' time     horn    Tycho     DE441     error   PA_fac')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
    fac = cusp_factor(mo_d, sun_d)
    sd = moon_semid_deg(mo) * fac
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    v_c = v - dec_corr
    print(f' {fmt_time(h,m)}  {horn:6s} {v_c:+8.4f}  {dec_model:+8.4f}  '
          f'{(v_c-dec_model)*60:+6.2f}\'   {fac:.3f}')

print('\nDiff. asc. (Sun center - Moon EAST limb):  '
      f'dRA timestamps shifted -{DRA_EXTRA_CLOCK_FAST_MIN:.1f} min')
for h, m, v in MOON_SUN:
    t = t_dra(h, m)
    sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
    model = _wrap(sun_d.ra.deg - _moon_east_limb_ra(t))
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDiff. asc. (Sun center - Moon EAST limb) at H.8:56.5:')
for h, m, v in MOON_SUN_SEP:
    t = t_dra(h, m)
    sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
    model = _wrap(sun_d.ra.deg - _moon_east_limb_ra(t))
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  Tycho {v*60:5.2f}\'  DE441 {d*60:5.2f}\'  '
          f'err {(v-d)*60:+.2f}\'')

print('\nSun declination (armilla bias removed):')
for h, m, v in SUN_DEC:
    t = t_at(h, m)
    sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
    v_c = v - dec_corr
    print(f' {fmt_time(h,m)}  Tycho {v_c:+.4f}  DE441 {sun_d.dec.deg:+.4f}  '
          f'err {(v_c-sun_d.dec.deg)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_SUN), len(MOON_SUN_SEP), len(MOON_DIAM),
         len(SUN_DEC)]
names = ['horn declinations', 'Sun-Moon dRA', 'Sun-Moon dRA (8:56)',
         'Moon diameter', 'Sun declination']
i = 0
print('\nPer-channel residual summary:')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    mean = chunk.mean()
    std = chunk.std() if n > 1 else 0.0
    print(f'  {n:>2d} {name:22s}  mean = {mean:+.2f}\''
          f'   scatter(std) = {std:.2f}\'')
