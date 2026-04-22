"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 14 October 1586 (Julian), evening, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 14 Octobris, Ad Vesperas").

Setup
-----
  Calendar: 14 Oct 1586 Julian  =  24 Oct 1586 Gregorian
  Clock:    Evening observation with NOON origin (post-meridiem).
            Tycho verified at the end of the session (by Aldebaran's
            distance from the meridian) that his clock was running
            +3.5 min fast over the interval since mean noon:
              "Ergo tunc visum est horologium iusto celerius promotum
               fuisse M.3.5 habita ratione temporis a Meridie..."
            We apply only the standard equation-of-time correction;
            the clock gain is left visible in the residuals.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference stars
---------------
  "Prima alae Pegasi" = Markab (alpha Pegasi).
  "Lucida Vulturis"   = Altair (alpha Aquilae).
  Both are propagated from J2000 back to 1586 with proper motion.

Limb convention
---------------
  Both "Occid. limb. ( ad or." (Markab block) and "Occid. limb. (
  orient." (Altair block) = the Moon's WEST limb (the limb rising
  first, on the eastern side of the lunar disc).
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
NOON_LOCAL = Time('1586-10-24 12:00:00', scale='ut1') - LOCAL_OFFSET

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


# Reference stars (J2000 ICRS) with proper motion.
# Markab (alpha Peg): Hipparcos mu_a* = 61.1 mas/yr, mu_d = -42.56 mas/yr,
# parallax 23.36 mas (~43 pc).
markab = SkyCoord(
    ra='23h04m45.65s', dec='+15d12m18.96s',
    pm_ra_cosdec=61.10*u.mas/u.yr,
    pm_dec=-42.56*u.mas/u.yr,
    distance=(1.0/0.02336)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)

# Altair (alpha Aql): mu_a* = 536.23, mu_d = 385.29 mas/yr, parallax 194.95 mas.
altair = SkyCoord(
    ra='19h50m46.99s', dec='+08d52m05.9s',
    pm_ra_cosdec=536.23*u.mas/u.yr,
    pm_dec=385.29*u.mas/u.yr,
    distance=(1.0/0.19495)*u.pc,
    obstime='J2000.0',
    frame='icrs',
)


def star_at(sc: SkyCoord, t: Time) -> SkyCoord:
    """Apply proper motion to a star from its J2000 epoch to time t."""
    return sc.apply_space_motion(new_obstime=t)


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 14 OCTOBRIS 1586 P.M.)
# All decs NORTH (Bor. = boream).
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    (7, 28 + 1/6, 'upper', +(1 + 15.0/60)),
    (7, 30 + 1/3, 'lower', +(    43.5/60)),
    (7, 46 + 2/3, 'upper', +(1 + 17.5/60)),
    (7, 48 + 1/3, 'lower', +(    46.0/60)),
    (8,  4 + 1/6, 'upper', +(1 + 23.0/60)),
    (8,  6.0,     'lower', +(    50.0/60)),
    (8,  7.0,     'upper', +(1 + 21.0/60)),
]

# Moon-Markab right-ascension difference.  The table header is
# "Prima alae Pegasi or. | Occid. limb. ( ad or. | Dist. aequat.";
# with col2 = star's armillary reading "ad orientem" and
# col3 = Moon west-limb reading "ad orientem", the equatorial
# distance is col3 - col2 = RA(Moon_west_limb) - RA(Markab).
MOON_MARKAB = [
    (7, 35 + 2/3,  (37 + 57.5/60) - (18 + 43.25/60)),  # 19 14.25
    (7, 37 + 1/3,  (37 + 34.5/60) - (18 + 18.0/60)),   # 19 16.5
    (7, 38 + 1/6,  (37 + 10.5/60) - (17 + 53.5/60)),   # 19 17.0
    (7, 40.0,      (36 + 48.5/60) - (17 + 31.0/60)),   # 19 17.5
    (7, 43 + 1/6,  (36 +  6.0/60) - (16 + 48.5/60)),   # 19 17.5
]

# Moon-Altair right-ascension difference.  Header "Luc. Vult. occid. |
# Occid. limb. ( orient. | Dist. aequat."; col2 = Altair "ad occasum",
# col3 = Moon west limb "ad orientem", so col4 = col2 + col3 =
# RA(Moon_west_limb) - RA(Altair).
MOON_ALTAIR = [
    (7, 54.5,     (34 + 24.0/60) + (33 + 22.0/60)),    # 67 46
    (7, 58.5,     (35 + 23.0/60) + (32 + (24 + 1/3)/60)),  # 67 47 2/3
    (8,  0.5,     (35 + 55.0/60) + (31 + 54.0/60)),    # 67 49
]

# Moon apparent diameter
MOON_DIAM = [
    (7, 30 + 1/3, 31.5/60,  'armilla subterr.'),
    (7, 48 + 1/3, 31.5/60,  'Q. min.'),
    (8,  7.0,     31.0/60,  'armilla subterr.'),
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-24')


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
    for h, m, v in MOON_ALTAIR:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        al_d = apparent_of_date(star_at(altair, t), t)
        sd = moon_semid_deg(mo)
        model = ((mo_d.ra.deg - sd) - al_d.ra.deg) % 360.0
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
print('DIE 14 OCTOBRIS 1586 P.M. (Julian)  ==  24 Oct 1586 Gregorian')
print('Uraniborg.  Evening Moon observations (post-meridiem).')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local NOON (apparent).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982) and stellar")
print("proper motion for Markab and Altair (J2000 -> 1586).")
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

print('\nDiff. asc. (Moon WEST limb - Markab = "Prima alae Pegasi"):')
for h, m, v in MOON_MARKAB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    mk_d = apparent_of_date(star_at(markab, t), t)
    sd = moon_semid_deg(mo)
    model = ((mo_d.ra.deg - sd) - mk_d.ra.deg) % 360.0
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDiff. asc. (Moon WEST limb - Altair = "Luc. Vult."):')
for h, m, v in MOON_ALTAIR:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    al_d = apparent_of_date(star_at(altair, t), t)
    sd = moon_semid_deg(mo)
    model = ((mo_d.ra.deg - sd) - al_d.ra.deg) % 360.0
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:18s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

sizes = [len(MOON_DEC), len(MOON_MARKAB), len(MOON_ALTAIR), len(MOON_DIAM)]
names = ['horn declinations', 'Moon-Markab dRA', 'Moon-Altair dRA',
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
