"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 3 February 1586 (Julian) at Uraniborg, Hven.

Source: Tycho Brahe, *Opera Omnia*
        ("De Luna pro ipsius parallaxi habenda. Die 3 Februarii").

Setup
-----
  Calendar: 3 Feb 1586 Julian  =  13 Feb 1586 Gregorian
  Clock:    This is a MORNING observation ("Sol ortus erat" -- Sun had
            risen), so Tycho's hour count runs from local MIDNIGHT (civil
            origin), not from the usual astronomical noon.  Verified: his
            H.7:56 meridian transit matches the Moon's computed transit
            near 08:00 local apparent time at Uraniborg on Feb 13 Greg.
            Only the equation-of-time shift is applied (no fitted clock).
  Location: Uraniborg, lat 55d54'N, lon 12d42'E (height 40 m).
  Ephemeris: JPL DE441 long kernel.

This entry is unusual: no star-distance observations -- the whole page is
dedicated to cross-checking the Moon's declination at meridian transit
(iuxta H.7:56) using THREE different instruments simultaneously, as a
parallax-investigation prelude.  The Moon is near 26 deg Aquarius, so the
declination is SOUTHERN (negative in the modern sign convention).

  * Armillae magnae  -> horn declinations at several times,
                        plus derived Moon diameter
  * Sextants         -> horn altitudes at meridian -> center dec
  * Mural quadrant   -> horn altitudes at meridian -> center dec

The sextant/mural measurements are reported by Tycho already reduced to
center declination; we compare those to our topocentric of-date Moon dec
directly (horn == 'center').

Note on the transcription: the phrase "in 26 videlicet Aquarii" ("in 26
deg Aquarius") is almost certainly an OCR/typographic misreading of
"26 Sagittarii" (glyphs aquarius vs sagittarii are very similar).  The
Moon's actual ecliptic longitude on this date is ~266 deg = 26 deg
Sagittarius; in Aquarius it would have been ~55 deg east of its true
position, incompatible with meridian transit at H.7:56 local apparent.
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
# Morning observation -> Tycho's H.0 is local MIDNIGHT (civil), not noon.
MIDNIGHT_LOCAL = Time('1586-02-13 00:00:00', scale='ut1') - LOCAL_OFFSET

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


def of_date(sc: SkyCoord, t: Time) -> SkyCoord:
    return sc.transform_to(PrecessedGeocentric(equinox=t, obstime=t))


def moon_semid_deg(moon_sc: SkyCoord) -> float:
    R_MOON_KM = 1737.4
    return float(np.degrees(np.arcsin(R_MOON_KM /
                                      moon_sc.distance.to(u.km).value)))


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 3 FEBRUARIJ 1586, all southern declinations)
# ----------------------------------------------------------------------------
# Convention: dec values stored as SIGNED degrees (negative = southern).
MOON_DEC = [
    # --- Large armillary sphere (horn readings) ---
    (7, 34.0, 'upper',  -(19 + 49.0/60)),
    (7, 36.0, 'lower',  -(20 + 18.0/60)),
    (7, 43.0, 'upper',  -(19 + 49.0/60)),
    (7, 45.0, 'lower',  -(20 + 17.0/60)),
    (7, 56.0, 'upper',  -(19 + 48.0/60)),   # "quasi in Merid."
    (7, 56.0, 'lower',  -(20 + 15.0/60)),
    (8,  2.0, 'upper',  -(19 + 48.5/60)),
    (8,  2.0, 'lower',  -(20 + 15.0/60)),
    # --- Sextants at meridian, center dec (derived by Tycho from alt) ---
    (7, 56.0, 'center', -(19 + 59.5/60)),
    # --- Mural quadrant at meridian, center dec ---
    (7, 56.0, 'center', -(19 + 58.5/60)),   # 19 deg 58' 30"
]

# Moon apparent diameters recorded through the evening
MOON_DIAM = [
    # (H, M, diameter_deg, source)
    (7, 36.0, 29.0/60, 'armillas'),
    (7, 45.0, 28.0/60, 'armillas'),
    (7, 56.0, 27.0/60, 'armillas'),   # "differ. M.27"
    (8,  2.0, 26.5/60, 'armillas'),
    (7, 56.0, 26.0/60, 'sextantes'),
    (7, 56.0, 30.0/60, 'muralem'),
]

# ----------------------------------------------------------------------------
# Time construction with EoT-only shift (NO fitted clock parameter)
# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-02-13')


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
        mo_d = of_date(mo, t)
        sd = moon_semid_deg(mo)
        if horn == 'upper':
            dec_model = mo_d.dec.deg + sd
        elif horn == 'lower':
            dec_model = mo_d.dec.deg - sd
        else:   # 'center'
            dec_model = mo_d.dec.deg
        res.append((v - dec_model)*60)
    for h, m, v, _src in MOON_DIAM:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        res.append((v - 2*moon_semid_deg(mo))*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 3 FEBRUARIJ 1586 (Julian)  ==  13 Feb 1586 Gregorian, Uraniborg')
print('Modern ephemeris: JPL DE441 (long).  Comparison with Tycho raw log.')
print('(Parallax cross-check: three instruments at meridian transit.)')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local MIDNIGHT (morning obs.).")
print(f"Equation-of-time shift to mean solar: +{EOT_SHIFT_MIN:.2f} min")
print("No further clock adjustment is applied.")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon declinations (southern, degrees):')
print(' time     horn     Tycho       DE441       error')
for h, m, horn, v in MOON_DEC:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = of_date(mo, t)
    sd = moon_semid_deg(mo)
    if horn == 'upper':
        dec_model = mo_d.dec.deg + sd
    elif horn == 'lower':
        dec_model = mo_d.dec.deg - sd
    else:
        dec_model = mo_d.dec.deg
    print(f' {fmt_time(h,m)}  {horn:6s}  {v:9.4f}   {dec_model:9.4f}   '
          f'{(v-dec_model)*60:+6.2f}\'')

print('\nMoon apparent diameter:')
print(' time     source      Tycho    DE441    error')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    diam_model = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:10s}  {v*60:5.2f}\'  {diam_model*60:5.2f}\'  '
          f'{(v-diam_model)*60:+6.2f}\'')

# Per-channel / per-instrument summary
n_dec = len(MOON_DEC); n_dia = len(MOON_DIAM)
dec_r  = _r[:n_dec]
diam_r = _r[n_dec:n_dec+n_dia]

# Break declinations by horn type
arm_r  = np.array([r for (row, r) in zip(MOON_DEC, dec_r) if row[2] != 'center'])
sext_r = dec_r[-2]     # sextants (first 'center' entry)
mur_r  = dec_r[-1]     # mural (second 'center' entry)

# Break diameters by source
diam_arm  = np.array([r for (row, r) in zip(MOON_DIAM, diam_r) if row[3] == 'armillas'])
diam_sext = np.array([r for (row, r) in zip(MOON_DIAM, diam_r) if row[3] == 'sextantes'])
diam_mur  = np.array([r for (row, r) in zip(MOON_DIAM, diam_r) if row[3] == 'muralem'])

print('\nPer-channel residual summary (EoT shift applied, NO clock fit):')
print(f'  {len(arm_r):>2d} armillae horn decs   mean = {arm_r.mean():+.2f}\''
      f'   scatter(std) = {arm_r.std():.2f}\'')
print(f'   1 sextantes  center dec   value = {sext_r:+.2f}\'')
print(f'   1 muralem    center dec   value = {mur_r:+.2f}\'')
print(f'  {len(diam_arm):>2d} diameter (armillae)  mean = {diam_arm.mean():+.2f}\''
      f'   scatter(std) = {diam_arm.std():.2f}\'')
print(f'   1 diameter (sextantes)   value = {diam_sext[0]:+.2f}\'')
print(f'   1 diameter (muralem)     value = {diam_mur[0]:+.2f}\'')
