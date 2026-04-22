"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 23 January 1586 (Julian) at Uraniborg, Hven.

Source: Tycho Brahe, *Opera Omnia* ("Observationes Lunae, Die 23 Januarii").

Setup
-----
  Calendar: 23 Jan 1586 Julian  =  2 Feb 1586 Gregorian
  Clock:    Tycho's "H. h M. m" counts from local APPARENT noon (sundial
            noon). His clocks were stellar-calibrated to run at mean rate,
            but the day-origin follows Renaissance astronomical convention.
            The fixed equation-of-time shift (apparent noon - mean noon)
            is the ONLY time correction applied here (no clock fit).
  Location: Uraniborg, lat 55d54'N, lon 12d42'E (height 40 m).
  Ephemeris: JPL DE441 long kernel (13200 BC - 1969 AD).

Reference stars
---------------
  "Oculus Taurij" (Oculus ♀) = Aldebaran (alpha Tau).
  "Cor ☉" in the transcription is almost certainly "Cor ♌" = Cor Leonis =
    Regulus (alpha Leo); the ☉/♌ glyph confusion is an OCR/transcription
    artefact. Verified numerically: Regulus-Moon dRA matches Tycho's
    27.058 deg to ~2 arcmin.
"""
import warnings
import erfa
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import (
    solar_system_ephemeris, EarthLocation, SkyCoord, get_body,
    PrecessedGeocentric,
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

# -- Observer -----------------------------------------------------------------
URANIBORG_LAT_DEG = 55 + 54/60
URANIBORG_LON_DEG = 12 + 42/60
loc = EarthLocation(lat=URANIBORG_LAT_DEG*u.deg,
                    lon=URANIBORG_LON_DEG*u.deg,
                    height=40*u.m)
LOCAL_OFFSET = TimeDelta(URANIBORG_LON_DEG/15.0*3600.0, format='sec')
NOON_LOCAL = Time('1586-02-02 12:00:00', scale='ut1') - LOCAL_OFFSET

DE441_LONG = "/Users/ondrej/repos/python-skyfield/examples/de441_part-1.bsp"
solar_system_ephemeris.set(DE441_LONG)
# astropy's get_body('sun', ...) misbehaves on DE441 long kernel at 1586;
# we read the Sun directly via jplephem for the equation-of-time computation.
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
    """Minutes by which local apparent noon lags local mean noon at Uraniborg
    on the given Gregorian date (= amount to ADD to Tycho's apparent-noon-
    based hour count to obtain local mean solar time)."""
    def ha_deg(sec):
        t = (Time(date_str + ' 00:00:00', scale='ut1')
             + TimeDelta((12 - URANIBORG_LON_DEG/15.0)*3600 + sec, format='sec'))
        sun = sun_geocentric_icrs(t).transform_to(
            PrecessedGeocentric(equinox=t, obstime=t))
        lst = t.sidereal_time('apparent', longitude=loc.lon).deg
        return ((lst - sun.ra.deg + 540.0) % 360.0) - 180.0
    sec = brentq(ha_deg, -3600.0, 3600.0, xtol=1e-3)
    return sec / 60.0


def of_date(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.transform_to(PrecessedGeocentric(equinox=t, obstime=t))


def moon_semid_deg(moon_sc: SkyCoord) -> float:
    R_MOON_KM = 1737.4
    return float(np.degrees(np.arcsin(R_MOON_KM /
                                      moon_sc.distance.to(u.km).value)))


# Reference stars (J2000 ICRS, Hipparcos)
aldebaran = SkyCoord(ra='04h35m55.24s', dec='+16d30m33.5s', frame='icrs')
regulus   = SkyCoord(ra='10h08m22.31s', dec='+11d58m01.9s', frame='icrs')


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 23 IANUARIJ 1586)
# ----------------------------------------------------------------------------
# Moon horn declinations (per armillas). "vno"/"altero pinnacidio" = two
# independent sightings (same target, different pinnule) -- both kept as
# separate measurements.
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)
    (11, 38.0,    'lower', 14 + 30.5/60),   # 14 deg 30.5'  (vno pinn.)
    (11, 39.5,    'lower', 14 + 31.5/60),   # 14 deg 31.5'  (alt. pinn.)
    (12,  6.0,    'upper', 15 +  5.625/60), # 15 deg  5 5/8'
    (12,  6.5,    'upper', 15 +  6.5/60),   # 15 deg  6.5'  (alt. pinn.)
    (12, 24+5/8,  'lower', 14 + 25.25/60),  # 14 deg 25 1/4'
    (12, 25+5/8,  'lower', 14 + 26.0/60),   # 14 deg 26'    (alt. pinn.)
    (12, 30.0,    'upper', 15 +  1.0/60),   # 15 deg  1'
    (12, 31+1/8,  'upper', 15 +  0.5/60),   # 15 deg  0.5'  (alt. pinn.)
]

# Moon apparent diameter: "Ergo diameter C Min. 35." = 35', approx at H.12:30
MOON_DIAM = [
    (12, 30.0, 35.0/60),
]

# Moon occidental limb -> Regulus ("Cor [Leonis]"), Dift. aequat.:
REG_LIMB = [
    # (H, M, dRA_deg)   Times given as H M with fractional minutes.
    (11, 23+1/8,   27 +  3.5/60),     # 27 deg  3.5'
    (11, 26+5/8,   27 +  1.625/60),   # 27 deg  1 5/8'
    (11, 29+1/3,   27 +  0.5/60),     # 27 deg  0.5'    (footnote 2)
    (11, 32+1/3,   26 + 57 + 2/3*1/60/60 * 60 ),  # placeholder, will fix below
]
# Fix last entry: 26 deg 57 2/3'
REG_LIMB[-1] = (11, 32+1/3, 26 + (57 + 2/3)/60)

# Later pair of Regulus observations:
REG_LIMB += [
    (12, 13.0,     26 + 39.5/60),     # 26 deg 39.5'    (footnote 3)
    (12, 16+2/3,   26 + (40 + 1/3)/60),  # 26 deg 40 1/3'
]

# Moon occidental limb -> Aldebaran (Oculus ♀), Dift. aequat.:
ALD_LIMB = [
    (11, 52+1/6,   56 + 39.75/60),    # 56 deg 39 3/4'
    (11, 55+2/3,   56 + 40.5/60),     # 56 deg 40.5'
    (11, 59+1/6,   56 + 40.5/60),     # 56 deg 40.5'
    (12,  2+1/3,   56 + 41.5/60),     # 56 deg 41.5'
]

# ----------------------------------------------------------------------------
# Time construction with EoT-only shift (NO fitted clock parameter)
# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-02-02')


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
        mo_d = of_date(mo, t)
        sd = moon_semid_deg(mo)
        dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
        res.append((v - dec_model)*60)
    for h, m, v in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        res.append((v - 2*moon_semid_deg(mo))*60)
    for h, m, v in REG_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = of_date(mo, t); reg_d = of_date(regulus, t)
        sd = moon_semid_deg(mo)
        model = reg_d.ra.deg - (mo_d.ra.deg - sd)
        res.append((v - model)*60)
    for h, m, v in ALD_LIMB:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = of_date(mo, t); ald_d = of_date(aldebaran, t)
        sd = moon_semid_deg(mo)
        model = (mo_d.ra.deg - sd) - ald_d.ra.deg
        res.append((v - model)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 23 IANUARIJ 1586 (Julian)  ==  2 Feb 1586 Gregorian, Uraniborg')
print('Modern ephemeris: JPL DE441 (long).  Comparison with Tycho raw log.')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local APPARENT noon.")
print(f"Equation-of-time shift to mean solar: +{EOT_SHIFT_MIN:.2f} min")
print("No further clock adjustment is applied.")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn declinations (per armillas):')
print(' time     horn    Tycho     DE441     error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t)
    sd = moon_semid_deg(mo)
    dec_model = mo_d.dec.deg + (sd if horn == 'upper' else -sd)
    print(f' {fmt_time(h,m)}  {horn:5s}  {v:8.4f}  {dec_model:8.4f}  '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
for h, m, v in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    diam_model = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  Tycho {v*60:5.2f}\'   '
          f'DE441 {diam_model*60:5.2f}\'   err {(v-diam_model)*60:+.2f}\'')

print('\nDift. aequat. (Moon occid. limb - Regulus = "Cor Leonis"):')
print(' time     Tycho     DE441     error')
for h, m, v in REG_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t); reg_d = of_date(regulus, t)
    sd = moon_semid_deg(mo)
    # Regulus has HIGHER RA than Moon -> dRA = RA_Reg - RA_Moon_limb
    model = reg_d.ra.deg - (mo_d.ra.deg - sd)
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nDift. aequat. (Moon occid. limb - Aldebaran = "Oculus Taurij"):')
print(' time     Tycho     DE441     error')
for h, m, v in ALD_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t); ald_d = of_date(aldebaran, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg - sd) - ald_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

# Per-channel summary
n_dec = len(MOON_DEC); n_dia = len(MOON_DIAM)
n_reg = len(REG_LIMB); n_ald = len(ALD_LIMB)
i = 0
dec_r  = _r[i:i+n_dec]; i += n_dec
diam_r = _r[i:i+n_dia]; i += n_dia
reg_r  = _r[i:i+n_reg]; i += n_reg
ald_r  = _r[i:i+n_ald]
print('\nPer-channel residual summary (EoT shift applied, NO clock fit):')
for name, x in [(f'{n_dec} horn declinations ', dec_r),
                (f'{n_dia} Moon diameter      ', diam_r),
                (f'{n_reg} Moon-Regulus dRA   ', reg_r),
                (f'{n_ald} Moon-Aldebaran dRA ', ald_r)]:
    if len(x) >= 2:
        print(f'  {name}  mean = {x.mean():+.2f}\'   '
              f'scatter(std) = {x.std(ddof=1):.2f}\'')
    else:
        print(f'  {name}  value = {x[0]:+.2f}\'')
