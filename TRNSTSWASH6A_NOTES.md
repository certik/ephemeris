# Notes on trnstswash6a.txt

Source: <https://ssd.jpl.nasa.gov/dat/planets/trnstswash6a.txt>
Format: O'Handley (1968), JPL Technical Report 32-1296
        <https://ssd.jpl.nasa.gov/dat/planets/optical_format.txt>

## Planet code 010 is the Sun, not the Moon

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

## All observations are geocentric

Every record has `obs_type = 1` (transit, geocentric). There are **no**
topocentric observations (`obs_type = 4`) in this file.

## Comparison with the USNO catalog PDF

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
geocentric positions with parallax already removed, as confirmed by:

1. The `obs_type = 1` flag (transit, geocentric)
2. The USNO catalog text (page 183): *"all positions given in this catalog
   are geocentric and refer to the center of the object observed"*
