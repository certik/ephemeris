"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 24 September 1586 (Julian), morning + daytime, at
Uraniborg.  This is an exceptional entry -- simultaneous Sun and Moon
observations in broad daylight.

Source: Tycho Brahe, *Opera Omnia* ("Die 24 Septembris A.M.").

Setup
-----
  Calendar: 24 Sep 1586 Julian  =  4 Oct 1586 Gregorian
  Clock:    Morning/daytime observation with civil-midnight origin.
            Tycho explicitly notes "tempora annotata sunt verificata ad
            Solem" -- times were verified against Sun transit, so the
            clock was calibrated that very day.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E.
  Ephemeris: JPL DE441 long kernel.

Reference stars & Sun
---------------------
  "Luc. ♀"  = Lucida Arietis = Hamal (alpha Ari)     ["II" in OCR]
  "Oculus ♀" = Oculus Tauri  = Aldebaran (alpha Tau)
  Sun observations are given with the Moon east limb (block 1) and west
  limb (block 2) separately.

Tycho's own note on this entry:
  "Ambo luminaria tam alta supra Horizontem, ut nulla sensibilis
  refractio insinuari potuerit" -- he considered both bodies high enough
  that refraction was negligible.  For the morning block Moon altitude
  reached ~52 deg at transit (refraction ~0.8'), but for the DAYTIME
  Sun-Moon block the Moon had sunk to altitude ~20 deg (refraction ~2.6')
  -- so in fact refraction was not negligible there.  Our model applies
  it; see whether Tycho's intuition held up.
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
MIDNIGHT_LOCAL = Time('1586-10-04 00:00:00', scale='ut1') - LOCAL_OFFSET

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


def sun_semid_deg(sun_sc: SkyCoord) -> float:
    R_SUN_KM = 695700.0
    return float(np.degrees(np.arcsin(R_SUN_KM /
                                      sun_sc.distance.to(u.km).value)))


# Reference stars (J2000 ICRS)
hamal = SkyCoord(ra='02h07m10.41s', dec='+23d27m44.7s', frame='icrs')
aldebaran = SkyCoord(ra='04h35m55.24s', dec='+16d30m33.5s', frame='icrs')


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 24 SEPTEMBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    # Morning armillary:
    (4,  8.0,  'upper', 18 +  8.0/60),
    (4,  8.0,  'lower', 17 + 37.0/60),
    (4, 32.0,  'upper', 18 +  8.0/60),
    (4, 32.0,  'lower', 17 + 37.25/60),
    # H.4:46 meridian transit altitudes (convert alt -> dec at meridian):
    (4, 46.0,  'upper', 52 + 15.0/60 - (90 - URANIBORG_LAT_DEG)),
    (4, 46.0,  'lower', 51 + 45.0/60 - (90 - URANIBORG_LAT_DEG)),
    # Daytime armillary (first burst, "super. cornu"):
    (9, 35.0,  'upper', 18 +  3.5/60),
    (9, 36.0,  'upper', 18 +  5.0/60),
    (9, 38.0,  'upper', 18 +  3.5/60),
    # Daytime armillary (second burst, "super. limbi"):
    (9, 55.0,  'upper', 18 +  6.0/60),
    (9, 56.5,  'upper', 18 +  6.0/60),
    (9, 57.75, 'upper', 18 +  5.625/60),
    # Daytime armillary (third burst):
    (10, 14.0, 'upper', 18 +  5.0/60),
    (10, 15.5, 'upper', 18 +  4.5/60),
    (10, 17.0, 'upper', 18 +  5.0/60),
]

MOON_DIAM = [
    (4,  8.0, 31.0/60,   'armillas'),
    (4, 32.0, 30.75/60,  'armillas'),
    (4, 46.0, 30.0/60,   'meridian'),
]

# Moon east limb -> Hamal (Lucida Arietis).  Moon EAST of meridian, star
# WEST of meridian; col3 = col1 + col2 = RA_Moon_east_limb - RA_Hamal.
HAM_LIMB = [
    (4, 15.0,     54 + 13.0/60),
    (4, 18.0,     54 + 14.25/60),
    (4, 20 + 1/3, 54 + 14.75/60),
]

# Moon east limb -> Aldebaran (Oculus Tauri).
ALD_LIMB = [
    (4, 23.0, 17 + 11.25/60),
    (4, 25.0, 17 + 12.375/60),   # 17 12 3/8
    (4, 27.0, 17 + 12.75/60),
    (4, 29.0, 17 + 14.0/60),
]

# Sun - Moon equatorial distance during DAYTIME.  Tycho's own note at the
# end of this entry says explicitly: "Vbique in distantia aequatoria
# accepta est orientalis limbus ( Soli proximus" -- "EVERYWHERE in the
# equatorial distance we took the Moon's EASTERN limb (nearest the Sun)".
# So both blocks use the east limb.  The "( occid." header of block 2
# refers to the MOON (in the western sky past meridian), not the limb.
SUN_MOON_1 = [
    (9, 47.0, 107 +  8.5/60),
    (9, 50.0, 107 +  7.0/60),
    (9, 52.0, 107 +  6.0/60),
]

SUN_MOON_2 = [
    (10,  0.0, 107 +  3.0/60),
    (10,  2.0, 107 +  1.0/60),
    (10,  4.0, 107 +  0.5/60),
    (10,  6.0, 107 +  0.5/60),
    (10,  7.5, 106 + 59.0/60),
]


# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-04')


def t_at(h, m):
    return MIDNIGHT_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0, format='sec')


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


def sun_of_date(t: Time) -> SkyCoord:
    """Geocentric Sun (of date), used to build a topocentric SkyCoord."""
    return sun_geocentric_icrs(t)


def sun_topo(t: Time) -> SkyCoord:
    """For Sun, parallax between geocentric and topocentric is ~9 arcsec
    (observer offset 6400 km vs 1 AU baseline) -- negligible at our ~0.5'
    precision.  Just return the geocentric ICRS Sun."""
    return sun_geocentric_icrs(t)


def residuals_arcmin():
    res = []
    for h, m, horn, v in MOON_DEC:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sd = moon_semid_deg(mo)
        dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
        res.append((v - dec_model)*60)
    for h, m, v, _src in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        res.append((v - 2*moon_semid_deg(mo))*60)
    for h, m, v in HAM_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        ham_d = apparent_of_date(hamal, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg + sd) - ham_d.ra.deg
        res.append((v - model)*60)
    for h, m, v in ALD_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        ald_d = apparent_of_date(aldebaran, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg + sd) - ald_d.ra.deg
        res.append((v - model)*60)
    for h, m, v in SUN_MOON_1:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sn = sun_topo(t)
        sn_d = apparent_of_date(sn, t)
        sd = moon_semid_deg(mo)
        # RA_Sun - RA_Moon_east_limb
        model = sn_d.ra.deg - (mo_d.ra.deg + sd)
        res.append((v - model)*60)
    for h, m, v in SUN_MOON_2:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sn = sun_topo(t)
        sn_d = apparent_of_date(sn, t)
        sd = moon_semid_deg(mo)
        # Block 2 also east limb per Tycho's own note at end of entry
        model = sn_d.ra.deg - (mo_d.ra.deg + sd)
        res.append((v - model)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 24 SEPTEMBRIS 1586 A.M. (Julian)  ==  4 Oct 1586 Gregorian')
print('Uraniborg.  Morning armillary + DAYTIME Sun/Moon observations.')
print('Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local MIDNIGHT (morning/day obs.).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982).")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn declinations (armillary + meridian):')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    print(f' {fmt_time(h,m)}  {horn:5s}  {v:8.4f}  {dec_model:8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    d = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:10s}  Tycho {v*60:5.2f}\'  '
          f'DE441 {d*60:5.2f}\'  err {(v-d)*60:+.2f}\'')

print('\nDift. aequat. (Moon east limb - Hamal = "Lucida Arietis"):')
for h, m, v in HAM_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    ham_d = apparent_of_date(hamal, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg + sd) - ham_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDift. aequat. (Moon east limb - Aldebaran = "Oculus Tauri"):')
for h, m, v in ALD_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    ald_d = apparent_of_date(aldebaran, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg + sd) - ald_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nSun - Moon east limb equatorial distance (DAYTIME, block 1):')
for h, m, v in SUN_MOON_1:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sn = sun_topo(t); sn_d = apparent_of_date(sn, t)
    sd = moon_semid_deg(mo)
    model = sn_d.ra.deg - (mo_d.ra.deg + sd)
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nSun - Moon east limb equatorial distance (DAYTIME, block 2):')
for h, m, v in SUN_MOON_2:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    sn = sun_topo(t); sn_d = apparent_of_date(sn, t)
    sd = moon_semid_deg(mo)
    model = sn_d.ra.deg - (mo_d.ra.deg + sd)
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

sizes = [len(MOON_DEC), len(MOON_DIAM), len(HAM_LIMB), len(ALD_LIMB),
         len(SUN_MOON_1), len(SUN_MOON_2)]
names = ['horn declinations', 'Moon diameters', 'Moon-Hamal   dRA',
         'Moon-Aldebaran dRA', 'Sun-Moon dRA (b1)', 'Sun-Moon dRA (b2)']
i = 0
print('\nPer-channel residual summary (refraction + EoT applied):')
for n, name in zip(sizes, names):
    chunk = _r[i:i+n]; i += n
    if n > 1:
        print(f'  {n:>2d} {name:22s}  mean = {chunk.mean():+.2f}\''
              f'   scatter(std) = {chunk.std():.2f}\'')
    else:
        print(f'  {n:>2d} {name:22s}  value = {chunk[0]:+.2f}\'')
