"""
Modern topocentric ephemeris reproducing Tycho Brahe's raw notebook
observations for 23 September 1586 (Julian), morning, at Uraniborg.

Source: Tycho Brahe, *Opera Omnia* ("Die 23 Septembris A.M.").

Setup
-----
  Calendar: 23 Sep 1586 Julian  =  3 Oct 1586 Gregorian
  Clock:    Morning observation.  Tycho's "H.h M.m A.M." = local civil time
            counted from local MIDNIGHT (confirmed by the note
            "H.2 58' A.M. vel post M.N." -- 'before noon or after
            midnight').  Only the equation-of-time shift is applied.
  Location: Uraniborg, lat 55d54'N, lon 12d42'E (height 40 m).
  Ephemeris: JPL DE441 long kernel.

Reference stars (verified numerically)
--------------------------------------
  "Lucida ♀" in the transcription is "Lucida ♈" = Hamal (alpha Ari).
  Verified: his Asc.R. 26d02.5' vs Hamal's of-date RA 26d03' -- agreement
  to 0.5'.  The ♀ ⇄ ♈ glyph confusion is a standard OCR artefact.

  "Inferius caput ♊" ('lower head of the Twins') = Pollux (beta Gem).
  Verified: his sextant angular distance to Moon 41d09.5' matches Pollux
  to 9' (vs 97' for Castor) -- clear identification.
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

# -- Observer -----------------------------------------------------------------
URANIBORG_LAT_DEG = 55 + 54/60
URANIBORG_LON_DEG = 12 + 42/60
loc = EarthLocation(lat=URANIBORG_LAT_DEG*u.deg,
                    lon=URANIBORG_LON_DEG*u.deg,
                    height=40*u.m)
LOCAL_OFFSET = TimeDelta(URANIBORG_LON_DEG/15.0*3600.0, format='sec')
# Morning, civil-midnight origin:
MIDNIGHT_LOCAL = Time('1586-10-03 00:00:00', scale='ut1') - LOCAL_OFFSET

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
    """Bennett (1982) mean refraction at 1010 mbar, 10 C; alt in deg."""
    if alt_deg <= -1.0:
        return 0.0
    return (1.0/np.tan(np.radians(alt_deg + 7.31/(alt_deg + 4.4))))/60.0


def apparent_of_date(sc: SkyCoord, t: Time) -> SkyCoord:
    """Topocentric of-date RA/Dec with atmospheric refraction applied."""
    aa = sc.transform_to(AltAz(obstime=t, location=loc))
    R = _bennett_refraction_deg(float(aa.alt.deg))
    app = SkyCoord(alt=(aa.alt.deg + R)*u.deg, az=aa.az.deg*u.deg,
                   frame=AltAz(obstime=t, location=loc))
    return app.transform_to(PrecessedGeocentric(equinox=t, obstime=t))


def moon_semid_deg(moon_sc: SkyCoord) -> float:
    R_MOON_KM = 1737.4
    return float(np.degrees(np.arcsin(R_MOON_KM /
                                      moon_sc.distance.to(u.km).value)))


# Reference stars (J2000 ICRS)
hamal  = SkyCoord(ra='02h07m10.41s', dec='+23d27m44.7s', frame='icrs')
pollux = SkyCoord(ra='07h45m18.95s', dec='+28d01m34.3s', frame='icrs')


# ----------------------------------------------------------------------------
# Tycho's raw notebook data (DIE 23 SEPTEMBRIS 1586 A.M.)
# ----------------------------------------------------------------------------
# All decs are NORTHERN (positive in modern convention).
MOON_DEC = [
    # (H, M, 'upper'|'lower', dec_deg)  -- from large armillary spheres
    (2, 58.0, 'upper', 17 + 18.0/60),      # primo pinnacidio
    (2, 58.0, 'upper', 17 + 17.625/60),    # altero pinnacidio
    (3,  7.0, 'upper', 17 + 18.5/60),
    (3,  7.0, 'lower', 16 + 45.0/60),
    (3, 27.5, 'lower', 16 + 45.75/60),
    (3, 27.5, 'upper', 17 + 19.0/60),
    (3, 50.0, 'upper', 17 + 20.5/60),
    (3, 50.0, 'lower', 16 + 47.5/60),
    # H.3:53 meridian transit of Moon EASTERN LIMB: horn altitudes -> dec.
    # At meridian, dec = alt - (90 - phi), with phi=Uraniborg lat.
    # These are altitude-based (refraction 1:1 in dec); our 'apparent_of_date'
    # already accounts for refraction via alt/az round-trip, so they are
    # consistent with the armillary decs if we treat as regular horn decs.
    (3, 53.0, 'upper', 51 + 25.5/60 - (90 - URANIBORG_LAT_DEG)),  # sextant
    (3, 53.0, 'lower', 50 + 51.5/60 - (90 - URANIBORG_LAT_DEG)),  # Q.Volub.
]

# Moon apparent diameter measurements through the morning.
MOON_DIAM = [
    (3,  7.0, 33.5/60, 'armillas'),
    (3, 27.5, 33.25/60, 'armillas'),
    (3, 50.0, 33.0/60, 'armillas'),
    (3, 53.0, 34.0/60, 'meridian'),    # from meridian altitudes
]

# Moon EASTERN limb -> Hamal (Lucida Arietis), "Dist. aequat." = dRA.
# Moon east of meridian, Hamal west of meridian -> col1 + col2 = RA_Moon - RA_Hamal.
HAM_LIMB = [
    (3, 11 + 1/3,  40 + 39.5/60),     # 40 deg 39.5'
    (3, 16.5,      40 + 40.0/60),     # 40 deg 40'
    (3, 18.0,      40 + 41.25/60),    # 40 deg 41 1/4'
    (3, 20 + 5/8,  40 + 42.5/60),     # 40 deg 42.5'
    (3, 22 + 5/8,  40 + 43.0/60),     # 40 deg 43'
    (3, 25.5,      40 + 44.75/60),    # 40 deg 44 3/4'
]

# Moon EASTERN limb -> Pollux (inferius caput Geminorum).
# "Dist." via sextant -> TRUE ANGULAR DISTANCE on sky, not dRA.
POL_LIMB_SEP = [
    (3, 34 + 5/8,  41 + 9.75/60),     # 41 deg 9 3/4'
    (3, 36 + 5/8,  41 + 9.5/60),      # 41 deg 9.5'
    (3, 39.0,      41 + 9.5/60),      # 41 deg 9.5'
]


# ----------------------------------------------------------------------------
# Time construction
# ----------------------------------------------------------------------------
EOT_SHIFT_MIN = equation_of_time_shift_min('1586-10-03')


def t_at(h, m):
    return MIDNIGHT_LOCAL + TimeDelta((h*60 + m + EOT_SHIFT_MIN)*60.0, format='sec')


def fmt_time(h, m):
    total = h*60 + m
    hh, mm = divmod(total, 60)
    return f"{int(hh):2d}:{mm:05.2f}"


def moon_limb_on_sky(mo_apparent: SkyCoord, limb_dir: SkyCoord, sd_deg: float):
    """Return a SkyCoord at the Moon's apparent limb in the direction of
    'limb_dir' (another SkyCoord).  Great-circle displacement by sd_deg
    from the Moon's apparent center toward limb_dir."""
    # Unit vectors on celestial sphere
    def uvec(sc):
        ra = np.radians(sc.ra.deg); de = np.radians(sc.dec.deg)
        return np.array([np.cos(de)*np.cos(ra), np.cos(de)*np.sin(ra),
                         np.sin(de)])
    m = uvec(mo_apparent); d = uvec(limb_dir)
    # Component of d perpendicular to m, normalized
    tang = d - m*np.dot(m, d)
    tang /= np.linalg.norm(tang)
    ang = np.radians(sd_deg)
    v = m*np.cos(ang) + tang*np.sin(ang)
    ra = np.degrees(np.arctan2(v[1], v[0])) % 360.0
    dec = np.degrees(np.arcsin(v[2]))
    return SkyCoord(ra=ra*u.deg, dec=dec*u.deg,
                    frame=mo_apparent.frame.replicate_without_data())


def residuals_arcmin():
    res = []
    for h, m, horn, v in MOON_DEC:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        sd = moon_semid_deg(mo)
        if horn == 'upper':
            dec_model = mo_d.dec.deg + sd
        else:
            dec_model = mo_d.dec.deg - sd
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
        # Moon east limb -> add semidiameter to RA
        model = (mo_d.ra.deg + sd) - ham_d.ra.deg
        res.append((v - model)*60)
    for h, m, v in POL_LIMB_SEP:
        t = t_at(h, m)
        mo = get_body('moon', t, location=loc)
        mo_d = apparent_of_date(mo, t)
        pol_d = apparent_of_date(pollux, t)
        sd = moon_semid_deg(mo)
        limb = moon_limb_on_sky(mo_d, pol_d, sd)
        # Great-circle distance from Moon east limb to Pollux
        model = limb.separation(pol_d).deg
        res.append((v - model)*60)
    return np.array(res)


_r = residuals_arcmin()
rms = np.sqrt(np.mean(_r**2))


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print('=' * 78)
print('DIE 23 SEPTEMBRIS 1586 A.M. (Julian)  ==  3 Oct 1586 Gregorian (morning)')
print('Uraniborg.  Modern ephemeris: JPL DE441 (long).')
print('=' * 78)
print()
print("Tycho's clock (H.0) convention: local MIDNIGHT (morning obs.).")
print(f"Equation-of-time shift to mean solar: {EOT_SHIFT_MIN:+.2f} min")
print("Model includes atmospheric refraction (Bennett 1982).")
print(f"All-observation RMS residual: {rms:.2f}'")

print('\nMoon horn declinations (armillary + meridian altitudes):')
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
print(' time     source      Tycho    DE441    error')
for h, m, v, src in MOON_DIAM:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    diam_model = 2*moon_semid_deg(mo)
    print(f' {fmt_time(h,m)}  {src:10s}  {v*60:5.2f}\'  {diam_model*60:5.2f}\'  '
          f'{(v-diam_model)*60:+6.2f}\'')

print('\nDift. aequat. (Moon east limb - Hamal = "Lucida Arietis"):')
print(' time     Tycho     DE441     error')
for h, m, v in HAM_LIMB:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    ham_d = apparent_of_date(hamal, t)
    sd = moon_semid_deg(mo)
    model = (mo_d.ra.deg + sd) - ham_d.ra.deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

print('\nAngular distance via sextant (Moon east limb - Pollux = "inf. caput Gem."):')
print(' time     Tycho     DE441     error')
for h, m, v in POL_LIMB_SEP:
    t = t_at(h, m)
    mo = get_body('moon', t, location=loc)
    mo_d = apparent_of_date(mo, t)
    pol_d = apparent_of_date(pollux, t)
    sd = moon_semid_deg(mo)
    limb = moon_limb_on_sky(mo_d, pol_d, sd)
    model = limb.separation(pol_d).deg
    print(f' {fmt_time(h,m)}  {v:8.4f}  {model:8.4f}  {(v-model)*60:+6.2f}\'')

# Per-channel summary
n_dec = len(MOON_DEC); n_dia = len(MOON_DIAM)
n_ham = len(HAM_LIMB); n_pol = len(POL_LIMB_SEP)
dec_r  = _r[:n_dec]
diam_r = _r[n_dec:n_dec+n_dia]
ham_r  = _r[n_dec+n_dia:n_dec+n_dia+n_ham]
pol_r  = _r[n_dec+n_dia+n_ham:n_dec+n_dia+n_ham+n_pol]

print('\nPer-channel residual summary (refraction + EoT applied):')
print(f'  {n_dec:>2d} horn declinations   mean = {dec_r.mean():+.2f}\''
      f'   scatter(std) = {dec_r.std():.2f}\'')
print(f'  {n_dia:>2d} Moon diameters     mean = {diam_r.mean():+.2f}\''
      f'   scatter(std) = {diam_r.std():.2f}\'')
print(f'  {n_ham:>2d} Moon-Hamal   dRA    mean = {ham_r.mean():+.2f}\''
      f'   scatter(std) = {ham_r.std():.2f}\'')
print(f'  {n_pol:>2d} Moon-Pollux  sep    mean = {pol_r.mean():+.2f}\''
      f'   scatter(std) = {pol_r.std():.2f}\'')
