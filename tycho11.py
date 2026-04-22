"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 16 October 1586 (Julian), evening, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 16 Octobris").

Setup
-----
  Calendar: 16 Oct 1586 Julian  =  26 Oct 1586 Gregorian
  Clock:    Evening observation with NOON origin (post-meridiem).
            Before observing, Tycho calibrated the clock against
            Aldebaran's meridian distance:
              "Ergo M.1 Sec.24 deficiebat Horologium, tantundem autem
               in continenti promotus est eius index, vt sequentia
               tempora observationum vera essent."
            i.e. he advanced the clock by 1m24s, so the recorded
            H.9:13+ etc. should be on the correct local apparent
            solar time.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference star
--------------
  Both blocks use Markab (alpha Pegasi): "Prima alae Pegasi" and
  "Marcab" are Tycho's two names for the same star.

Limb convention
---------------
  "Occid. limb. ( or." / "Occid. limb. ( orient." -- the Moon's WEST
  limb (the limb on the eastern side of the disc, which rises first).
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
NOON_LOCAL = Time('1586-10-26 12:00:00', scale='ut1') - LOCAL_OFFSET

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


# Markab = alpha Pegasi.  Hipparcos: mu_a* = 61.10, mu_d = -42.56 mas/yr,
# parallax 23.36 mas.
markab = SkyCoord(
    ra='23h04m45.65s', dec='+15d12m18.96s',
    pm_ra_cosdec=61.10*u.mas/u.yr,
    pm_dec=-42.56*u.mas/u.yr,
    distance=(1.0/0.02336)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)


def star_at(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.apply_space_motion(new_obstime=t)


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 16 OCTOBRIS 1586 P.M.)
# All decs NORTH (B. = boream).
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    (9, 29.5,    'upper', +(8 + 47.5/60)),
    (9, 31.0,    'lower', +(8 + 16.5/60)),
    (9, 48.75,   'upper', +(8 + 50.25/60)),
    (9, 50 + 1/6,'lower', +(8 + 18.25/60)),
]

# Moon-Markab right-ascension difference.  Both notebook blocks are
# Moon vs Markab; in both blocks col4 = col2 + col3 =
# RA(Moon_west_limb) - RA(Markab).  "Prima alae Pegasi" and "Marcab"
# are two names for the same star.
MOON_MARKAB = [
    (9, 13 + 5/8,  ( 8 + 30.0/60) + (34 + 10.0/60)),   # 42 40
    (9, 17 + 1/3,  ( 9 + 27.0/60) + (33 + 16.0/60)),   # 42 43
    (9, 19 + 1/4,  (10 +  1.25/60) + (32 + 42.5/60)),  # 42 43.75
    (9, 23 + 1/3,  (10 + 56.0/60) + (31 + 50.0/60)),   # 42 46
    (9, 26.5,      (11 + (41 + 1/3)/60) + (31 + 5.0/60)),  # 42 46 1/3
    (9, 33.5,      (13 + 29.0/60) + (29 + 19.0/60)),   # 42 48
    (9, 36.5,      (14 + 12.0/60) + (28 + 39.0/60)),   # 42 51
    (9, 38 + 1/3,  (14 + 43.25/60) + (28 + 8.5/60)),   # 42 51.75
    (9, 39 + 5/8,  (15 +  4.0/60) + (27 + 49.75/60)),  # 42 53.75
    (9, 45.0,      (16 + 19.0/60) + (26 + 36.0/60)),   # 42 55  ("mediocriter bona")
]

MOON_DIAM = [
    (9, 31.0,     30.625/60,  'armilla'),
    (9, 50 + 1/6, 32.0/60,    'armilla'),
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-26')


def t_at(h, m):
    return NOON_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0, format='sec')


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
        if horn == 'upper':
            dec_model = mo_d.dec.deg + sd
        elif horn == 'lower':
            dec_model = mo_d.dec.deg - sd
        else:
            dec_model = mo_d.dec.deg
        res.append((v - dec_model)*60)
    for h, m, v in MOON_MARKAB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        mk_d = apparent_of_date(star_at(markab, t), t)
        sd = moon_semid_deg(mo)
        model = ((mo_d.ra.deg - sd) - mk_d.ra.deg) % 360.0
        res.append((v - model)*60)
    for h, m, v, _ in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        d = 2*moon_semid_deg(mo)
        res.append((v - d)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 16 OCTOBRIS 1586 P.M. (Julian)  ==  26 Oct 1586 Gregorian')
print('Uraniborg.  Evening Moon observations (post-meridiem).')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local NOON (apparent).")
print("Clock was corrected by +1m24s via Aldebaran meridian distance")
print("before the observation sequence started.")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982) and stellar")
print("proper motion for Markab (J2000 -> 1586).")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn/limb declinations:')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    if horn == 'upper':
        dec_model = mo_d.dec.deg + sd
    elif horn == 'lower':
        dec_model = mo_d.dec.deg - sd
    else:
        dec_model = mo_d.dec.deg
    print(f' {fmt_time(h,m)}  {horn:6s} {v:+8.4f}  {dec_model:+8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nDiff. asc. (Moon WEST limb - Markab = "Prima alae Pegasi"/"Marcab"):')
for h, m, v in MOON_MARKAB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    mk_d = apparent_of_date(star_at(markab, t), t)
    sd = moon_semid_deg(mo)
    model = ((mo_d.ra.deg - sd) - mk_d.ra.deg) % 360.0
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:14s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_MARKAB), len(MOON_DIAM)]
names = ['horn declinations', 'Moon-Markab dRA', 'Moon diameter']
i = 0
print('\nPer-channel residual summary (refraction + EoT applied):')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    if n > 1:
        print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
              f'   scatter(std) = {chunk.std():.2f}\'')
    else:
        print(f'  {n:>2d} {name:22s}  value = {chunk[0]:+.2f}\'')
