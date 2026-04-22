"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 25 October 1586 (Julian), morning, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 25 Octobris").

Setup
-----
  Calendar: 25 Oct 1586 Julian  =  4 Nov 1586 Gregorian
  Clock:    Morning observations -> MIDNIGHT origin.
            Clock was rectified at H.4:30 A.M.  Tycho checked the
            clock against the Sun at H.10:19:30 and found true time
            = H.10:16 (clock 5.5 min fast per his note; 3.5 min fast
            as computed from the noted Sun distance).  This script
            applies a linear drift of 5.5 min over (10:19.5 - 4:30)
            as Tycho explicitly stated.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference bodies
----------------
  "Cor ☉" = Cor Leonis = Regulus (alpha Leo).
  "η"     = Jupiter (confirmed by DE441 match to <0.01' for its
            declination in tycho15 on the previous night, and by
            Tycho's own alternation with "Iouis").

Table geometry
--------------
  Regulus block (6:18): "Cor ☉ orient. | Lunae orient. limbus occid."
    col4 = col2 + col3 = RA(Regulus) - RA(Moon_east_limb).
  Jupiter blocks (7:09, 7:26): "η occid. | ( or. limbus occ."
    col4 = col2 - col3 = RA(Moon_east_limb) - RA(Jupiter).
  (Azimuth/altitude via Quadrans Volubilis and the clock-check Sun
   observations at H.10:19-10:21 are omitted: azimuth is not in this
   suite's channels, and the Sun check is what defines the drift model.)
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
MIDNIGHT_LOCAL = Time('1586-11-04 00:00:00', scale='ut1') - LOCAL_OFFSET

DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)
_SPK = SPK.open(DE441_LONG)

# Clock drift model (Tycho's own words):
#   rectified at H.4:30, fast by 5.5 min at H.10:19:30 -> linear.
DRIFT_START_MIN = 4*60 + 30.0
DRIFT_END_MIN = 10*60 + 19.5
DRIFT_TOTAL_MIN = 5.5


def clock_drift(clock_min_from_midnight: float) -> float:
    frac = ((clock_min_from_midnight - DRIFT_START_MIN)
            / (DRIFT_END_MIN - DRIFT_START_MIN))
    return DRIFT_TOTAL_MIN * frac


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
# Tycho's raw notebook data (DIE 25 OCTOBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
MOON_DEC = [
    (6, 15.0,       'upper', +(12 + 36.5/60)),
    (6, 15.0,       'lower', +(12 +  4.0/60)),
    (7, 13.5,       'upper', +(12 + 27.0/60)),
    (7, 15.0,       'lower', +(11 + 58.0/60)),
    (7, 16.0,       'upper', +(12 + 27.0/60)),
    (7, 17.0,       'lower', +(11 + 58.5/60)),
    (7, 35.0,       'upper', +(12 + 24.5/60)),
    (7, 37.0,       'lower', +(11 + 55.5/60)),
    (7, 38.0,       'upper', +(12 + 24.5/60)),
    (7, 40.0,       'lower', +(11 + 55.25/60)),
]

# Regulus - Moon east limb.  Using Tycho's recorded Dist. aequat. (col4)
# directly, since row 3 has an internal inconsistency between (col2+col3)
# and col4 in the printed notebook.
MOON_REG = [
    (6, 18.5,       14 + (42 + 1/3)/60),
    (6, 19 + 35/60, 14 + 42.75/60),
    (6, 22 +  5/60, 14 + (40 + 1/3)/60),
]

# Jupiter (η) - Moon east limb.  col4 = col2 - col3 = RA(Moon) - RA(Jup).
MOON_JUP = [
    (7,  9 + 45/60, (39 + 46.0/60)   - (14 + 16.0/60)),
    (7, 11 + 28/60, (40 + 11.5/60)   - (14 + 40.5/60)),
    (7, 26 + 52/60, (44 +  2.0/60)   - (18 + 24.5/60)),
    (7, 29.0,       (44 + 34.0/60)   - (18 + 54.0/60)),
]

# Moon apparent diameter
MOON_DIAM = [
    (7, 17.0, 28.5/60),
    (7, 35.0, 29.0/60),
    (7, 40.0, 29.0/60),
]

# Transit altitude of east limb (meridian) at 6:12.
MOON_TRANSIT_ALT = [
    (6, 12.0, 'upper-horn', 46 + (42 + 1/3)/60),
    (6, 12.0, 'lower-limb', 46 + 11.5/60),
]


EOT_SHIFT_MIN = equation_of_time_shift_min('1586-11-04')


def t_at(h, m):
    clock_min = h*60 + m
    true_min = clock_min - clock_drift(clock_min)
    return MIDNIGHT_LOCAL + TimeDelta((true_min + EOT_SHIFT_MIN)*60.0,
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
        dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
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
    for h, m, kind, v in MOON_TRANSIT_ALT:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        aa = mo_d.transform_to(AltAz(obstime=t, location=loc))
        sd = moon_semid_deg(mo)
        alt_model = aa.alt.deg + (sd if kind == 'upper-horn' else -sd)
        res.append((v - alt_model)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


print('=' * 78)
print('DIE 25 OCTOBRIS 1586 A.M. (Julian)  ==  4 Nov 1586 Gregorian')
print('Uraniborg.  Clock drift: linear 5.5 min from 4:30 to 10:19.5.')
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

print('\nMeridian transit altitude (Moon):')
for h, m, kind, v in MOON_TRANSIT_ALT:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    aa = mo_d.transform_to(AltAz(obstime=t, location=loc))
    sd = moon_semid_deg(mo)
    alt_model = aa.alt.deg + (sd if kind == 'upper-horn' else -sd)
    print(f' {fmt_time(h,m)}  {kind:11s} Tycho {v:+.4f}  DE441 '
          f'{alt_model:+.4f}  err {(v-alt_model)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_REG), len(MOON_JUP), len(MOON_DIAM),
         len(MOON_TRANSIT_ALT)]
names = ['horn declinations', 'Regulus-Moon dRA', 'Moon-Jupiter dRA',
         'Moon diameter', 'Transit altitude']
i = 0
print('\nPer-channel residual summary:')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    mean = chunk.mean() if n else 0.0
    std = chunk.std() if n > 1 else 0.0
    print(f'  {n:>2d} {name:22s}  mean = {mean:+.2f}\''
          f'   scatter(std) = {std:.2f}\'')
