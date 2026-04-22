"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 26 September 1586 (Julian), morning + daytime, at
Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 26 Septembris A.M.").

Setup
-----
  Calendar: 26 Sep 1586 Julian  =  6 Oct 1586 Gregorian
  Clock:    Morning/daytime observation with civil-midnight origin.
            Tycho explicitly notes "Tempora assignata sunt correcta"
            -- times were calibrated.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Only Sun observations are used (no reference stars).  Tycho notes:
  "orientalem limbum ( observatum vbique esse a centro Solis"
  -- Moon's EAST limb was always observed from the Sun's center.

He also cautions: "habenda etiam est ratio parallaxeos lunaris, eo
quod luna a 90 Gradu plurimum removeatur versus occasum" -- lunar
parallax must be accounted for because the Moon is far west of the
meridian.  Our model uses full topocentric coordinates so this is
handled automatically.
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
MIDNIGHT_LOCAL = Time('1586-10-06 00:00:00', scale='ut1') - LOCAL_OFFSET

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


def sun_topo(t: Time) -> SkyCoord:
    return sun_geocentric_icrs(t)


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 26 SEPTEMBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    (8, 32.0,     'lower', 16 +  7.0/60),      # non satis exacta
    (8, 58.0,     'lower', 16 +  4.5/60),      # melior priori
    (9,  0.0,     'upper', 16 + (30 + 1/3)/60),
    (9,  2.0,     'upper', 16 + (31 + 1/3)/60),  # melior
    (9,  3.0,     'upper', 16 + (31 + 1/3)/60),  # bona
    (9,  4.0,     'lower', 16 +  3.5/60),
    (9,  5.0,     'lower', 16 +  3.0/60),       # bona
    (9, 17.5,     'upper', 16 + 29.0/60),
    (9, 19.5,     'upper', 16 + 29.25/60),
    (9, 20.5,     'lower', 16 +  2.0/60),
    (9, 22.0,     'lower', 16 +  1.75/60),
    (9, 50 + 1/3, 'upper', 16 + 26.5/60),
    (9, 52.0,     'lower', 15 + 59.25/60),
    (10, 24.5,    'upper', 16 + 22.5/60),
    (10, 27 + 5/8,'lower', 15 + (54 + 1/3)/60),
    (10, 29.0,    'lower', 15 + (54 + 2/3)/60),
]

# Sun - Moon east limb equatorial distance.  Tycho explicitly states the
# east limb ("orientalem limbum ( ... a centro Solis") was used throughout.
# col3 = col1 + col2 = RA_Sun - RA_Moon_east_limb (Moon east of mer., Sun
# east too, i.e., HA_Moon > HA_Sun).
#
# Block 1: morning, both bodies east of meridian.
SUN_MOON_1 = [
    (9,  6.75,   82 +  3.0/60),
    (9,  9.0,    82 +  3.5/60),
    (9, 10 + 2/3, 82 +  2.0/60),
    (9, 12 + 2/3, 82 +  0.75/60),
    (9, 14 + 1/6, 82 +  0.0/60),
]

# Block 2: Moon still east / just past transit.  Starred "Hae tres erant
# bonae" applies to the first three of the preceding block (9:10-14).
SUN_MOON_2 = [
    (9, 23 + 1/3, 81 + 55.0/60),
    (9, 26.0,     81 + 53.5/60),
    (9, 36.0,     81 + 53.5/60),
    (9, 38.0,     81 + 52.0/60),
    (9, 42.0,     81 + 50.5/60),
    (9, 45 + 1/3, 81 + 49.0/60),   # bonae
]

# Block 3: header "( occ. / Sun orient." -- Moon past meridian to the west.
SUN_MOON_3 = [
    (10,  0.0,    81 + 43.5/60),
    (10,  5.0,    81 + 41.0/60),
    (10,  7.0,    81 + 39.5/60),
    (10,  8 + 2/3, 81 + 38.0/60),
    # (10, 15 + 2/3, 81 + 37.5/60),  # marked "non bona" -- skip
    (10, 20.5,    81 + 34.0/60),
    (10, 21.5,    81 + 33.5/60),
]

SUN_MOON_NONBONA = [
    (10, 15 + 2/3, 81 + 37.5/60),  # shown separately
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-06')


def t_at(h, m):
    return MIDNIGHT_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0, format='sec')


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
    for block in (SUN_MOON_1, SUN_MOON_2, SUN_MOON_3):
        for h, m, v in block:
            t = t_at(h, m)
            mo = get_body('moon', t, location=loc)
            mo_d = apparent_of_date(mo, t)
            sn = sun_topo(t)
            sn_d = apparent_of_date(sn, t)
            sd = moon_semid_deg(mo)
            model = sn_d.ra.deg - (mo_d.ra.deg + sd)
            res.append((v - model)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 26 SEPTEMBRIS 1586 A.M. (Julian)  ==  6 Oct 1586 Gregorian')
print('Uraniborg.  Morning armillary + DAYTIME Sun/Moon observations.')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local MIDNIGHT (morning/day obs.).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982).")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn/limb declinations:')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    print(f' {fmt_time(h,m)}  {horn:5s}  {v:8.4f}  {dec_model:8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')


def print_sun_moon_block(label, block):
    print(f'\n{label}:')
    for h, m, v in block:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sn = sun_topo(t); sn_d = apparent_of_date(sn, t)
        sd = moon_semid_deg(mo)
        model = sn_d.ra.deg - (mo_d.ra.deg + sd)
        print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')


print_sun_moon_block('Sun - Moon east limb dRA (block 1, both east of mer.)',
                     SUN_MOON_1)
print_sun_moon_block('Sun - Moon east limb dRA (block 2)', SUN_MOON_2)
print_sun_moon_block('Sun - Moon east limb dRA (block 3, Moon past meridian)',
                     SUN_MOON_3)
print_sun_moon_block('Sun - Moon east limb dRA (flagged "non bona")',
                     SUN_MOON_NONBONA)

sizes = [len(MOON_DEC), len(SUN_MOON_1), len(SUN_MOON_2), len(SUN_MOON_3)]
names = ['horn declinations', 'Sun-Moon dRA (b1)', 'Sun-Moon dRA (b2)',
         'Sun-Moon dRA (b3)']
i = 0
print('\nPer-channel residual summary (refraction + EoT applied):')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    if n > 1:
        print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
              f'   scatter(std) = {chunk.std():.2f}\'')
    else:
        print(f'  {n:>2d} {name:22s}  value = {chunk[0]:+.2f}\'')
