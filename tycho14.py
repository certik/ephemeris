"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 23 October 1586 (Julian), morning, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 23 Octobris A.M.").

Setup
-----
  Calendar: 23 Oct 1586 Julian  =  2 Nov 1586 Gregorian
  Clock:    Morning observations (A.M.) -> MIDNIGHT origin.
            Tycho notes that the clock was corrected via Cor Leonis
            and was running well thereafter.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference bodies
----------------
  "Cor ☉" = Cor Leonis = Regulus (alpha Leo).
  "☉"     = Sun (for the equatorial distances at H.9:14+).

Table geometry
--------------
  Block 1 (Regulus): header "Cor ☉ orient. | Orient. limbus ( occid.";
    Regulus is EAST of Moon, Moon east limb WEST; col4 = col2 + col3 =
    RA(Regulus) - RA(Moon_east_limb).
  Sun block:  "Dist. aequat. centri ☉ & orient. limbi (" --
    RA(Sun) - RA(Moon_east_limb).
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
MIDNIGHT_LOCAL = Time('1586-11-02 00:00:00', scale='ut1') - LOCAL_OFFSET

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


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 23 OCTOBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
MOON_DEC = [
    (6,  4.0,       'upper', +(17 + (12 + 2/3)/60)),
    (6,  4.0,       'lower', +(16 + 42.5/60)),
    (6, 40.0,       'upper', +(17 + 10.0/60)),
    (6, 40.0,       'lower', +(16 + 38.0/60)),
    (9,  9.0,       'upper', +(16 + 57.25/60)),
    (9, 10.0,       'lower', +(16 + (26 + 1/3)/60)),
    (9, 24.0,       'upper', +(16 + 55.0/60)),
    (9, 26.0,       'lower', +(16 + 25.0/60)),
    (9, 51 + 1/3,   'upper', +(16 + (51 + 1/3)/60)),
]

# Regulus - Moon_east_limb dRA.  col4 = col2 + col3.
MOON_REG = [
    (6, 15 + 1/3,  (26 + (16 + 5/8)/60) + (15 + 41.25/60)),    # 41 58.5
    (6, 17 + 1/3,  (26 + (45 + 2/3)/60) + (15 + 11.0/60)),     # 41 56.67
    (6, 19.0,      (27 + (14 + 1/6)/60) + (14 + 41.0/60)),     # 41 55.17
    (6, 20.5,      (27 + 35.0/60)       + (14 + 20.0/60)),     # 41 55
    (6, 25.0,      (28 + (37 + 2/3)/60) + (13 + 14.5/60)),     # 41 52.17
    (6, 30 + 1/3,  (29 + (58 + 2/3)/60) + (11 + 52.25/60)),    # 41 50.33
    (6, 34.5,      (30 + 48.0/60)       + (11 +  0.25/60)),    # 41 48.25
]

# Sun center - Moon east limb dRA ("Dist. aequat. centri ☉ & orient. limbi (")
MOON_SUN = [
    (9, 14.0,     111 + 18.0/60),
    (9, 16.0,     111 + 16.5/60),
    (9, 19.0,     111 + 15.0/60),
    (9, 21.0,     111 + 13.5/60),
    (9, 42.5,     111 +  3.0/60),  # "I"
    (9, 45.0,     111 +  2.0/60),  # "II"
]

MOON_DIAM = [
    (9, 10.0, (30 + 55.0/3600)/60, 'armilla (30\'55")'),
]


EOT_SHIFT_MIN = equation_of_time_shift_min('1586-11-02')


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
    for h, m, v in MOON_SUN:
        t = t_at(h, m)
        sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
        model = _wrap(sun_d.ra.deg - _moon_east_limb_ra(t))
        res.append((v - model)*60)
    for h, m, v, _ in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        d = 2*moon_semid_deg(mo)
        res.append((v - d)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


print('=' * 78)
print('DIE 23 OCTOBRIS 1586 A.M. (Julian)  ==  2 Nov 1586 Gregorian')
print('Uraniborg.  Morning/daytime Moon observations.')
print('=' * 78)
print()
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn/limb declinations:')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper'
                                else -sd if horn == 'lower' else 0.0)
    print(f' {fmt_time(h,m)}  {horn:6s} {v:+8.4f}  {dec_model:+8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nDiff. asc. (Regulus - Moon EAST limb = "Cor ☉"):')
for h, m, v in MOON_REG:
    t = t_at(h, m)
    reg_d = apparent_of_date(star_at(regulus, t), t)
    model = _wrap(reg_d.ra.deg - _moon_east_limb_ra(t))
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDiff. asc. (Sun center - Moon EAST limb):')
for h, m, v in MOON_SUN:
    t = t_at(h, m)
    sun_d = apparent_of_date(sun_geocentric_icrs(t), t)
    model = _wrap(sun_d.ra.deg - _moon_east_limb_ra(t))
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:18s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_REG), len(MOON_SUN), len(MOON_DIAM)]
names = ['horn declinations', 'Regulus-Moon dRA', 'Sun-Moon dRA',
         'Moon diameter']
i = 0
print('\nPer-channel residual summary:')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    if n > 1:
        print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
              f'   scatter(std) = {chunk.std():.2f}\'')
    else:
        print(f'  {n:>2d} {name:22s}  value = {chunk[0]:+.2f}\'')
