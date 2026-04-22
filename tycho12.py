"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 18 October 1586 (Julian), late evening, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 18 Octobris").

Setup
-----
  Calendar: 18 Oct 1586 Julian  =  28 Oct 1586 Gregorian
  Clock:    Evening observation with NOON origin (post-meridiem).
            Recorded times H.11:24 -- H.12:00 are late evening /
            approaching local midnight.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference stars
---------------
  "Oculus ♀" = Oculus Tauri = Aldebaran (alpha Tau).
  "Lucida ♀" = Lucida Arietis = Hamal (alpha Ari), as used by
               Tycho in the Sep 23-24 observations.  Identified
               from Tycho's col4 = 23.4 deg, which matches the
               RA(Moon) - RA(Hamal) separation to ~0.2 deg at
               this date (Moon at RA ~55 deg, Hamal at 31.79 deg).

Table geometry
--------------
  On this date the Moon sits BETWEEN Hamal and Aldebaran in RA
  (Moon RA ~ 55 deg, Hamal ~ 32 deg, Aldebaran ~ 69 deg), so:
  Block 1 (Aldebaran is EAST of Moon): col4 = col3 - col2 =
    RA(Aldebaran) - RA(Moon_east_limb).
  Block 2 (Hamal is WEST of Moon):     col4 = col2 + col3 =
    RA(Moon_east_limb) - RA(Hamal).

Limb convention
---------------
  "Orient. limb." throughout = the Moon's EAST limb (limbus
  orientalis, the leading limb in diurnal motion -- larger RA).
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
NOON_LOCAL = Time('1586-10-28 12:00:00', scale='ut1') - LOCAL_OFFSET

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


# Aldebaran (alpha Tau).  Hipparcos: mu_a* = 63.45, mu_d = -188.94 mas/yr,
# parallax 48.94 mas.
aldebaran = SkyCoord(
    ra='04h35m55.24s', dec='+16d30m33.5s',
    pm_ra_cosdec=63.45*u.mas/u.yr,
    pm_dec=-188.94*u.mas/u.yr,
    distance=(1.0/0.04894)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)

# Hamal (alpha Ari, Lucida Arietis).  Hipparcos: mu_a* = 188.55,
# mu_d = -148.08 mas/yr, parallax 14.66 mas.
hamal = SkyCoord(
    ra='02h07m10.41s', dec='+23d27m44.7s',
    pm_ra_cosdec=188.55*u.mas/u.yr,
    pm_dec=-148.08*u.mas/u.yr,
    distance=(1.0/0.01466)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)


def star_at(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.apply_space_motion(new_obstime=t)


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 18 OCTOBRIS 1586 P.M.)
# All decs NORTH (B. = boream).
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    (11, 24.0, 'upper', +(14 + 47.0/60)),
    (11, 26.0, 'lower', +(14 + 13.5/60)),
    (11, 49.0, 'upper', +(14 + 48.5/60)),
    (11, 51.0, 'lower', +(14 + 16.0/60)),
]

# Moon-Aldebaran RA difference ("Oculus ♀").  On this date the Moon is
# WEST of Aldebaran, so col4 = col3 - col2 = RA(Aldebaran) - RA(Moon_east_limb).
MOON_ALD = [
    (11, 33.0, (36 + 57.5/60) - (23 + 11.0/60)),   # 13 46.5
    (11, 34.0, (36 + 13.0/60) - (22 + 28.0/60)),   # 13 45
    (11, 41.0, (34 + 54.0/60) - (21 + 10.0/60)),   # 13 44
    (11, 43.5, (34 + 17.0/60) - (20 + 34.5/60)),   # 13 42.5
    (11, 45.0, (33 + 53.5/60) - (20 + 11.5/60)),   # 13 42
]

# Moon-Hamal RA difference ("Lucida ♀" = Lucida Arietis).  Moon is EAST
# of Hamal here, so col4 = col2 + col3 = RA(Moon_east_limb) - RA(Hamal).
MOON_HAM = [
    (11, 55.0,  (17 + 51.0/60) + ( 5 + 33.0/60)),  # 23 24
    (11, 57.5,  (17 +  6.5/60) + ( 6 + 15.5/60)),  # 23 22  (clouds; "non satis certa")
    (12,  0.0,  (16 + 40.5/60) + ( 6 + 44.5/60)),  # 23 25
]

MOON_DIAM = [
    (11, 26.0, 33.5/60,  'armilla'),
    (11, 51.0, 32.5/60,  'armilla'),
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-28')


def t_at(h, m):
    return NOON_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0, format='sec')


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


def _moon_east_limb_ra(t):
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    return mo_d.ra.deg + sd / np.cos(np.radians(mo_d.dec.deg))


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
        model = (ald_d.ra.deg - _moon_east_limb_ra(t)) % 360.0
        if model > 180.0:
            model -= 360.0
        res.append((v - model)*60)
    for h, m, v in MOON_HAM:
        t = t_at(h, m)
        ham_d = apparent_of_date(star_at(hamal, t), t)
        model = (_moon_east_limb_ra(t) - ham_d.ra.deg) % 360.0
        if model > 180.0:
            model -= 360.0
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
print('DIE 18 OCTOBRIS 1586 P.M. (Julian)  ==  28 Oct 1586 Gregorian')
print('Uraniborg.  Late-evening Moon observations (post-meridiem).')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local NOON (apparent).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982) and stellar")
print("proper motion for Aldebaran and Hamal (J2000 -> 1586).")
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

print('\nDiff. asc. (Aldebaran - Moon EAST limb = "Oculus ♀"):')
for h, m, v in MOON_ALD:
    t = t_at(h, m)
    ald_d = apparent_of_date(star_at(aldebaran, t), t)
    model = (ald_d.ra.deg - _moon_east_limb_ra(t)) % 360.0
    if model > 180.0:
        model -= 360.0
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDiff. asc. (Moon EAST limb - Hamal = "Lucida ♀" = Lucida Arietis):')
for h, m, v in MOON_HAM:
    t = t_at(h, m)
    ham_d = apparent_of_date(star_at(hamal, t), t)
    model = (_moon_east_limb_ra(t) - ham_d.ra.deg) % 360.0
    if model > 180.0:
        model -= 360.0
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:14s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_ALD), len(MOON_HAM), len(MOON_DIAM)]
names = ['horn declinations', 'Aldebaran-Moon dRA', 'Moon-Hamal dRA',
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
