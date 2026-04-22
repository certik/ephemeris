"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 13 October 1586 (Julian), evening, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 13 Octobris, Ad Vesperas").

Setup
-----
  Calendar: 13 Oct 1586 Julian  =  23 Oct 1586 Gregorian
  Clock:    Evening observation with NOON origin (post-meridiem).
            Clock zeroed at H.6:42 PM against Aldebaran's distance from
            the meridian.  Tycho himself flags a +3.5 min clock drift
            over the 5.5-hour session:
              "Horologium ab eo tempore quo coepimus observare hucusque
               M.3.5 celerius iusto movebatur"
            We apply only the standard equation-of-time correction; the
            clock drift is left visible in the residuals (roughly linear
            in time since 6:42 PM).
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference star
--------------
  "Lucida Vulturis" (Luc. Vult.) = Altair (alpha Aql), as on Oct 10.

Limb convention
---------------
  Tycho explicitly states: "In omnibus distantiis accipiebatur occidentalis
  limbus (" -- the Moon's WEST limb was used for all distance measurements.
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
NOON_LOCAL = Time('1586-10-23 12:00:00', scale='ut1') - LOCAL_OFFSET

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


# Reference star (J2000 ICRS)
altair = SkyCoord(ra='19h50m46.99s', dec='+08d52m05.9s', frame='icrs')


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 13 OCTOBRIS 1586 P.M.)
# All decs SOUTH (Mer. = meridionalis).
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower'|'center', dec_deg)
    (8, 30.0,   'upper',  -(2 + 21.0/60)),      # Per Armillas subterr.
    (8, 50.0,   'upper',  -(2 + 20.0/60)),      # per Armillas maior. subterr.
    (8, 51 + 1/6, 'lower', -(2 + 49.0/60)),
    (8, 52.5,   'lower',  -(2 + 51.0/60)),
    (8, 53.5,   'lower',  -(2 + 51.0/60)),
    (9, 26.0,   'upper',  -(2 + 11.25/60)),     # Moon on meridian
    (9, 26.0,   'lower',  -(2 + 44.0/60)),
]

# Moon-Altair right-ascension difference.  "Occid. limb. ( ad ort." --
# west limb of Moon, ad ortum (the limb that rose first).  Tycho's
# explicit N.B. confirms west limb used throughout.
# col4 = col2 + col3 = RA_Moon_west_limb - RA_Altair
MOON_ALT_LIMB = [
    (8, 35.5,     56 + 46.0/60),
    (8, 39.0,     56 + 45.75/60),
    (8, 42.0,     56 + (47 + 2/3)/60),
    (8, 44 + 1/6, 56 + 47.0/60),
    (8, 47.5,     56 + 48.0/60),
]

# Moon apparent diameter (explicit measurements)
MOON_DIAM = [
    (9, 26.0, (32 + 1/3)/60,  'armilla (merid.)'),
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-23')


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
    for h, m, v in MOON_ALT_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        alt_d = apparent_of_date(altair, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg - sd) - alt_d.ra.deg
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
print('DIE 13 OCTOBRIS 1586 P.M. (Julian)  ==  23 Oct 1586 Gregorian')
print('Uraniborg.  Evening Moon observations (post-meridiem).')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local NOON (apparent).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982).")
print(f"All-observation RMS residual: {rms:.2f}'")
print("Tycho flags a +3.5 min clock drift over the 5.5-h session")
print("(clock ran fast; Aldebaran check at H.12:07.5).")

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
    print(f' {fmt_time(h,m)}  {horn:6s} {v:8.4f}  {dec_model:8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nDiff. asc. (Moon WEST limb - Altair = "Luc. Vult."):')
for h, m, v in MOON_ALT_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    alt_d = apparent_of_date(altair, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg - sd) - alt_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:18s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_ALT_LIMB), len(MOON_DIAM)]
names = ['horn declinations', 'Moon-Altair dRA', 'Moon diameter']
i = 0
print('\nPer-channel residual summary (refraction + EoT applied):')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    if n > 1:
        print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
              f'   scatter(std) = {chunk.std():.2f}\'')
    else:
        print(f'  {n:>2d} {name:22s}  value = {chunk[0]:+.2f}\'')
