# Proving that historical transit circle data is geocentric

This document presents three independent lines of evidence that the USNO
and Greenwich transit circle observations we use are **geocentric** (parallax
removed), not topocentric (as observed from the telescope).

---

## 1. The obs_type flag in trnstswash6a.txt

Source: <https://ssd.jpl.nasa.gov/dat/planets/trnstswash6a.txt>
Format: O'Handley (1968), JPL Technical Report 32-1296
        <https://ssd.jpl.nasa.gov/dat/planets/optical_format.txt>

### Planet code 010 is the Sun, not the Moon

The file uses NAIF target IDs. In NAIF numbering, **10 = Sun** (not Moon,
which is 301). The planet codes present in the file are:

| Code | NAIF body | Count |
|------|-----------|------:|
| 001  | Mercury   | 1,873 |
| 002  | Venus     | 3,158 |
| 004  | Mars      |   654 |
| 005  | Jupiter   |   852 |
| 006  | Saturn    |   839 |
| 007  | Uranus    |   872 |
| 008  | Neptune   |   790 |
| 010  | **Sun**   | 6,587 |

There is no Moon data in this file.

### All observations are geocentric by flag

Every record has `obs_type = 1` (transit, geocentric). There are **no**
topocentric observations (`obs_type = 4`) in this file.

---

## 2. Byte-level comparison: trnstswash6a.txt vs USNO catalog PDF

The Sun positions in `trnstswash6a.txt` were compared against positions
extracted by OCR from the definitive USNO catalog
(1982PUSNO\_\_23\_\_165H.pdf, stored in `sun_transit_observations.csv`).

**Example — 1963 March 8:**

Raw line from `trnstswash6a.txt`:
```
P010243809722168000007866   191002313427080      0-045815190            1963USNO
```

Decoded: JD 2438097.22168, RA = 23h 13m 42.708s = 348.42795°, Dec = −4° 58′ 15.19″ = −4.970886°

From the PDF (`sun_transit_observations.csv`):
JD 2438097.22167, RA = 23h 13m 42.712s = 348.427967°, Dec = −4° 58′ 15.20″ = −4.970889°

**Difference: ΔRA = +0.06″, ΔDec = −0.01″** — pure fixed-width rounding.

Over 743 matched observations (1963–1971):

| Statistic    |  ΔRA (″) | ΔDec (″) |
|--------------|----------|-----------|
| Mean         |   +0.003 |    +0.004 |
| Std dev      |    0.018 |     0.006 |
| Max absolute |    0.060 |     0.012 |

The two datasets are identical to within parsing precision. Both contain
geocentric positions with parallax already removed, as confirmed by the
USNO catalog text (page 183): *"all positions given in this catalog are
geocentric and refer to the center of the object observed."*

---

## 3. Computational proof: Greenwich 1948 Moon observation

The most direct proof: we compute the Moon's position both from the Earth's
center (geocentric) and from the observatory surface (topocentric) and see
which one matches the published observation.

### Method

A transit circle measures the instant an object crosses the local meridian.
The **right ascension is the time measurement** — it equals the local
sidereal time at the moment of transit. So the published RA tells us the
exact time of observation. We:

1. Use the observed RA to find the precise UT of transit (by bisection)
2. At that UT, compute the Moon's **declination** two ways:
   - **Geocentric**: from the Earth's center
   - **Topocentric**: from the Greenwich Observatory surface (51.48° N)
3. Compare each with the published declination

### Test case

From *Greenwich Observations 1948*, page A12 (Jan 18, 1948):

- **Observed RA** = 1h 13m 27.56s
- **Observed Dec** = +4° 52′ 52.65″

Computation using DE440s ephemeris with IAU 2000A precession/nutation
(program: `app/greenwich_test.f90`):

```
=================================================================
 Greenwich Transit Circle - Moon - 1948 Jan 18
 Observed: RA = 1h 13m 27.56s, Dec = +4d 52' 52.65"
=================================================================

  TOPOCENTRIC (from Greenwich surface):
    Transit UT =    17.41978831 h
    Dec = +04d 12'  7.44"
    Dec residual (obs - comp):    2445.21"

  GEOCENTRIC (from Earth center):
    Transit UT =    17.42036521 h
    Dec = +04d 52' 53.96"
    Dec residual (obs - comp):      -1.31"

 =========================================
  Topocentric Dec residual:    2445.21"
  Geocentric  Dec residual:      -1.31"
  Topo - Geo (parallax):      -2446.52"
 =========================================
```

### Result

| Model       | Dec residual (obs − computed) |
|-------------|-------------------------------|
| Geocentric  | **−1.31″**                    |
| Topocentric | +2445″ (= +40.8′)            |

The observation matches the **geocentric** computation to **1.3 arcseconds**,
consistent with the ~1″ precision of transit circle data.

The topocentric position differs by **40.8 arcminutes** — exactly the
expected lunar diurnal parallax at Greenwich's latitude:

  Δδ ≈ HP × sin(φ − δ) ≈ 57′ × sin(51.5° − 5°) ≈ 57′ × 0.72 ≈ 41′

### Conclusion

The Greenwich 1948 transit circle data is **geocentric beyond any doubt**.
The parallax has been removed using the observatory's known position and
the Moon's distance from the contemporary ephemeris.

This is consistent with standard practice for all major transit circle
catalogs of the 19th and 20th centuries: parallax, refraction, aberration,
and orbital motion corrections were applied during data reduction, and only
geocentric positions were published.
