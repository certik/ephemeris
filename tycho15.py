"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 24 October 1586 (Julian), morning, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 24 Octobris A.M.").

Setup
-----
  Calendar: 24 Oct 1586 Julian  =  3 Nov 1586 Gregorian
  Clock:    Morning observations (A.M.) -> MIDNIGHT origin.
            Tycho notes the clock was verified via Cor Leonis at ~5 A.M.
            In the block at 6:46 he measures the Moon's east limb 21d13'
            west of the meridian as a secondary verification, consistent
            with the Cor Leonis calibration.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference bodies
----------------
  "Cor ☉" = Cor Leonis = Regulus (alpha Leo).
  "η"     = Jupiter. Tycho alternates the symbol with "Iouis"
            ("Iouis & (" sunt repetitae) in lines introducing
            the η tables.  Jupiter was near Gemini at this epoch
            with dec ~+22.5 deg, consistent with Tycho's
            "fuit Declin. η 22 36 B".
  (Sextant great-circle distances -- "per Sext." -- are skipped here
   as they use a different measurement channel than the equatorial
   differences this suite analyses.)

Table geometry
--------------
  Regulus blocks (5:27 / 6:12): "Cor ☉ orient. | Orient. limb. ( occid.",
    col4 = col2 + col3 = RA(Regulus) - RA(Moon_east_limb).
  Jupiter blocks (6:49, 7:20): "η occident. | Or. limbus ( occid.",
    col4 = col2 - col3 = RA(Moon_east_limb) - RA(Jupiter)
    (both objects west of meridian; Moon further east -> smaller HA).
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
MIDNIGHT_LOCAL = Time('1586-11-03 00:00:00', scale='ut1') - LOCAL_OFFSET

DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)
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


regulus = SkyCoord(
    ra='10h08m22.31s', dec='+11d58m01.9s',
    pm_ra_cosdec=-249.4*u.mas/u.yr, pm_dec=4.91*u.mas/u.yr,
    distance=(1.0/0.04113)*u.pc, obstime='J2000.0', frame='icrs',
)


def star_at(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.apply_space_motion(new_obstime=t)


def _wrap(x):
    x = x % 360.0
    return x - 360.0 if x > 180.0 else x


def _moon_east_limb_ra(t):
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    return mo_d.ra.deg + sd / np.cos(np.radians(mo_d.dec.deg))


def jupiter_apparent(t):
    j = get_body('jupiter', t, location=loc)
    return apparent_of_date(j, t)


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 24 OCTOBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
MOON_DEC = [
    (5, 20.0,       'upper', +(15 + (27 + 1/3)/60)),
    (5, 20.0,       'lower', +(14 + 55.0/60)),
    (5, 44 + 1/3,   'upper', +(15 + 24.0/60)),
    (5, 45.5,       'lower', +(14 + 51.0/60)),
    (6,  6.5,       'upper', +(15 + 21.0/60)),
    (6,  8.25,      'lower', +(14 + 49.5/60)),
    (6, 17.0,       'upper', +(15 + 22.0/60)),
    (6, 19.5,       'lower', +(14 + 56.25/60)),       # paired with 6:17 upper
    (6, 25.0,       'upper', +(15 + 20.0/60)),
    (6, 26.5,       'lower', +(14 + 48.0/60)),
    (6, 28.0,       'upper', +(15 + (19 + 1/4)/60)),
    (6, 30.0,       'lower', +(14 + 48.0/60)),
    (6, 31.5,       'upper', +(15 + (18 + 2/3)/60)),
    (6, 31.5,       'lower', +(14 + (47 + 2/3)/60)),
    (6, 53.0,       'upper', +(15 + 16.5/60)),
    (6, 55.0,       'lower', +(14 + 46.0/60)),
    (6, 56.0,       'upper', +(15 + 16.5/60)),
    (6, 57.5,       'lower', +(14 + 46.5/60)),
    (6, 59.0,       'upper', +(15 + 16.5/60)),
    (7,  0.0,       'lower', +(14 + 45.5/60)),
    (7, 32.5,       'upper', +(15 + 12.0/60)),
    (7, 33.0,       'lower', +(14 + 43.5/60)),
    (7, 36.0,       'upper', +(15 + 12.5/60)),
    (7, 37.5,       'lower', +(14 + 42.0/60)),
    (7, 44.5,       'upper', +(15 + 13.0/60)),
    (7, 45.5,       'lower', +(14 + 41.0/60)),
    (7, 47.5,       'upper', +(15 + 12.0/60)),
    (7, 48.5,       'lower', +(14 + 41.0/60)),
]

# Block 1 (Regulus east, Moon east limb west): col4 = col2 + col3
MOON_REG = [
    (5, 27.0,     (26 + 46.0/60)     + ( 1 + 56.5/60)),
    (5, 29.0,     (26 + 19.5/60)     + ( 2 + 22.0/60)),
    (5, 32.5,     (25 + 30.5/60)     + ( 3 +  9.25/60)),
    (5, 34.0,     (25 +  2.0/60)     + ( 3 + 36.5/60)),
    (5, 35.5,     (24 + 42.5/60)     + ( 3 + 56.0/60)),
    (5, 37.0,     (24 + 19.25/60)    + ( 4 + 19.0/60)),
    (5, 39.0,     (23 + 50.0/60)     + ( 4 + 47.25/60)),
    (5, 40.75,    (23 + 21.0/60)     + ( 5 + 16.0/60)),
    # Block 2 at 6:12 - same geometry; header swapped cols but same sum
    (6, 12.5,     (12 + 37.0/60)     + (15 + 44.5/60)),
    (6, 13 + 1/3, (12 + 58.0/60)     + (15 + 22.25/60)),
    (6, 15.0,     (13 + 14.5/60)     + (15 +  6.25/60)),
]

# Jupiter ("η") - Moon_east_limb dRA (Jupiter west of Moon).
# col4 = col2 - col3 = RA(Moon_east_limb) - RA(Jupiter)
MOON_JUP = [
    (6, 49.5,     (33 + 42.0/60)    - (21 + 55.5/60)),
    (6, 50 + 2/3, (33 + 59.5/60)    - (22 + 16.0/60)),
    (6, 51 + 1/3, (34 + 19.0/60)    - (22 + 35.0/60)),
    (7, 20.0,     (41 + 17.0/60)    - (29 + 18.5/60)),
    (7, 21.5,     (41 + 42.5/60)    - (29 + 44.5/60)),
    (7, 23.0,     (42 +  3.0/60)    - (30 +  0.0/60)),     # "dubia"
    (7, 26.5,     (42 + 59.25/60)   - (30 + 59.5/60)),     # "bona" = 12 00.5
]

MOON_DIAM = [
    (5, 20.0, 32.5/60),
    (5, 45.5, 33.0/60),
    (6, 55.0, 30.5/60),
    (7, 31.5, 31.5/60),
    (7, 37.5, 30.5/60),
    (7, 45.5, 32.0/60),
    (7, 48.5, 31.0/60),
]

# Jupiter declinations recorded by Tycho.
JUP_DEC = [
    (7,  2.0,  22 + 36.0/60),
    (7,  3.5,  22 + 36.0/60),
]


EOT_SHIFT_MIN = equation_of_time_shift_min('1586-11-03')


def t_at(h, m):
    return MIDNIGHT_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0,
                                      format='sec')


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


def residuals_arcmin():
    res = []
    for h, m, horn, v in MOON_DEC:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sd = moon_semid_deg(mo)
        dec_model = mo_d.dec.deg + (sd if horn == 'upper'
                                    else -sd if horn == 'lower' else 0.0)
        res.append((v - dec_model)*60)
    for h, m, v in MOON_REG:
        t = t_at(h, m)
        reg_d = apparent_of_date(star_at(regulus, t), t)
        model = _wrap(reg_d.ra.deg - _moon_east_limb_ra(t))
        res.append((v - model)*60)
    for h, m, v in MOON_JUP:
        t = t_at(h, m)
        jd = jupiter_apparent(t)
        model = _wrap(_moon_east_limb_ra(t) - jd.ra.deg)
        res.append((v - model)*60)
    for h, m, v in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        d = 2*moon_semid_deg(mo)
        res.append((v - d)*60)
    for h, m, v in JUP_DEC:
        t = t_at(h, m)
        jd = jupiter_apparent(t)
        res.append((v - jd.dec.deg)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


print('=' * 78)
print('DIE 24 OCTOBRIS 1586 A.M. (Julian)  ==  3 Nov 1586 Gregorian')
print('Uraniborg.  Moon near max epicycle remoteness.')
print('=' * 78)
print()
print(f"Equation-of-time shift: {EOT_SHIFT_MIN:+.2f} min")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn declinations:')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    print(f' {fmt_time(h,m)}  {horn:6s} {v:+8.4f}  {dec_model:+8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nDiff. asc. (Regulus - Moon EAST limb):')
for h, m, v in MOON_REG:
    t = t_at(h, m)
    reg_d = apparent_of_date(star_at(regulus, t), t)
    model = _wrap(reg_d.ra.deg - _moon_east_limb_ra(t))
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDiff. asc. (Moon EAST limb - Jupiter = "η"):')
for h, m, v in MOON_JUP:
    t = t_at(h, m)
    jd = jupiter_apparent(t)
    model = _wrap(_moon_east_limb_ra(t) - jd.ra.deg)
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  Tycho {v*60:5.2f}\'  DE441 {d*60:5.2f}\'  '
          f'err {(v-d)*60:+.2f}\'')

print('\nJupiter ("η") declination:')
for h, m, v in JUP_DEC:
    t = t_at(h, m)
    jd = jupiter_apparent(t)
    print(f' {fmt_time(h,m)}  Tycho {v:+.4f}  DE441 {jd.dec.deg:+.4f}  '
          f'err {(v-jd.dec.deg)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_REG), len(MOON_JUP), len(MOON_DIAM),
         len(JUP_DEC)]
names = ['horn declinations', 'Regulus-Moon dRA', 'Moon-Jupiter dRA',
         'Moon diameter', 'Jupiter dec']
i = 0
print('\nPer-channel residual summary:')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
          f'   scatter(std) = {chunk.std():.2f}\'')
