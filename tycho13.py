"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 22 October 1586 (Julian), morning, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 22 Octobris A.M.").

Setup
-----
  Calendar: 22 Oct 1586 Julian  =  1 Nov 1586 Gregorian
  Clock:    Morning observation with MIDNIGHT origin (A.M.).
            Recorded times H.3:32 -- H.5:48 are pre-dawn.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference stars
---------------
  "Oculus ♀" = Oculus Tauri = Aldebaran (alpha Tau).
  "Cor ☉"    = Cor Leonis   = Regulus  (alpha Leo).

Table geometry
--------------
  Moon sits between Aldebaran (RA ~69 deg) and Regulus (RA ~152 deg)
  at Moon RA ~ 96 deg this morning, so:
  Block 1 (Aldebaran west of Moon): header "( occ. | Oculus ♀ occ.";
    col4 = col2 - col3 = RA(Moon_east_limb) - RA(Aldebaran).
  Block 2 (Regulus east of Moon):   header "( occ. | Cor ☉ orient.";
    col4 = col2 + col3 = RA(Regulus) - RA(Moon_east_limb).

Limb convention
---------------
  "Orient. limb." throughout = Moon's EAST limb (limbus orientalis).
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
MIDNIGHT_LOCAL = Time('1586-11-01 00:00:00', scale='ut1') - LOCAL_OFFSET

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


# Aldebaran (alpha Tau).
aldebaran = SkyCoord(
    ra='04h35m55.24s', dec='+16d30m33.5s',
    pm_ra_cosdec=63.45*u.mas/u.yr,
    pm_dec=-188.94*u.mas/u.yr,
    distance=(1.0/0.04894)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)

# Regulus (alpha Leo).  Hipparcos: mu_a* = -249.4, mu_d = 4.91 mas/yr,
# parallax 41.13 mas.
regulus = SkyCoord(
    ra='10h08m22.31s', dec='+11d58m01.9s',
    pm_ra_cosdec=-249.4*u.mas/u.yr,
    pm_dec=4.91*u.mas/u.yr,
    distance=(1.0/0.04113)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)


def star_at(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.apply_space_motion(new_obstime=t)


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 22 OCTOBRIS 1586 A.M.)
# Decs NORTH (B. = boream).
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    (5,  8.0, 'upper', +(18 +  8.0/60)),
    (5,  8.0, 'lower', +(17 + 37.0/60)),  # "infer. limbi" alongside upper
    (5, 48.0, 'upper', +(18 + 11.0/60)),
    (5, 48.0, 'lower', +(17 + 38.0/60)),
]

# Moon-Aldebaran RA difference ("Oculus ♀").  Moon is EAST of Aldebaran.
# col4 = col2 - col3 = RA(Moon_east_limb) - RA(Aldebaran).
MOON_ALD = [
    (5, 20.0, (55 + (16 + 1/3)/60) - (27 + 39.5/60)),      # 27 37.25
    (5, 22.0, (55 + 49.0/60)       - (28 + 11.5/60)),      # 27 37.5
    (5, 25.0, (56 + 27.5/60)       - (28 + (48 + 5/8)/60)),# 27 38.375
]

# Moon-Regulus RA difference ("Cor ☉" = Cor Leonis).  Moon is WEST of
# Regulus.  col4 = col2 + col3 = RA(Regulus) - RA(Moon_east_limb).
MOON_REG = [
    (5, 42.0, (32 + 46.0/60) + (22 + 54.0/60)),    # 55 40
    (5, 46.0, (33 + 34.25/60) + (22 +  3.5/60)),   # 55 37.375
]

# Moon apparent diameter (from horn pairs measured close together in time).
MOON_DIAM = [
    (5,  8.0, (18 + 8.0/60) - (17 + 37.0/60),  'upper-lower (5:08)'),
    (5, 48.0, (18 + 11.0/60) - (17 + 38.0/60), 'upper-lower (5:48)'),
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-11-01')


def t_at(h, m):
    return MIDNIGHT_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0, format='sec')


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


def _moon_east_limb_ra(t):
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    return mo_d.ra.deg + sd / np.cos(np.radians(mo_d.dec.deg))


def _wrap(x):
    x = x % 360.0
    return x - 360.0 if x > 180.0 else x


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
    for h, m, v in MOON_ALD:
        t = t_at(h, m)
        ald_d = apparent_of_date(star_at(aldebaran, t), t)
        model = _wrap(_moon_east_limb_ra(t) - ald_d.ra.deg)
        res.append((v - model)*60)
    for h, m, v in MOON_REG:
        t = t_at(h, m)
        reg_d = apparent_of_date(star_at(regulus, t), t)
        model = _wrap(reg_d.ra.deg - _moon_east_limb_ra(t))
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
print('DIE 22 OCTOBRIS 1586 A.M. (Julian)  ==  1 Nov 1586 Gregorian')
print('Uraniborg.  Pre-dawn Moon observations.')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local MIDNIGHT (apparent).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982) and stellar")
print("proper motion for Aldebaran and Regulus (J2000 -> 1586).")
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

print('\nDiff. asc. (Moon EAST limb - Aldebaran = "Oculus ♀"):')
for h, m, v in MOON_ALD:
    t = t_at(h, m)
    ald_d = apparent_of_date(star_at(aldebaran, t), t)
    model = _wrap(_moon_east_limb_ra(t) - ald_d.ra.deg)
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDiff. asc. (Regulus - Moon EAST limb = "Cor ☉" = Cor Leonis):')
for h, m, v in MOON_REG:
    t = t_at(h, m)
    reg_d = apparent_of_date(star_at(regulus, t), t)
    model = _wrap(reg_d.ra.deg - _moon_east_limb_ra(t))
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter (upper-lower horn pair):')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:22s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_ALD), len(MOON_REG), len(MOON_DIAM)]
names = ['horn declinations', 'Moon-Aldebaran dRA', 'Regulus-Moon dRA',
         'Moon diameter']
i = 0
print('\nPer-channel residual summary (refraction + EoT applied):')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    if n > 1:
        print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
              f'   scatter(std) = {chunk.std():.2f}\'')
    else:
        print(f'  {n:>2d} {name:22s}  value = {chunk[0]:+.2f}\'')
