"""
Convert Sun transit circle observations into an observations file for
orbit-determination programs.

Reads:
  - sun_transit_observations.csv  (Sun RA/Dec from USNO Pub. Vol. 23, 1963–1971)

Note: trnstswash6a.txt planet code 010 was initially assumed to be the Moon
but is actually the Sun (RA rate ~0.91°/day, Dec within ±23.44°).  That file
contains no Moon transit data — only planets + Sun.

Sorts by JD and writes transit_observations.dat.
"""

import csv
import sys


def load_sun(path: str):
    """Load Sun transit observations, return list of (jd, ra_deg, dec_deg, dist_au)."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip incomplete rows
            if not row['ra_h'] or not row['ra_s'] or not row['dec_deg'] or not row['dec_sec']:
                continue
            jd = float(row['jd_ephem']) + 2400000.0
            ra = (float(row['ra_h']) * 15.0
                  + float(row['ra_m']) * (15.0 / 60.0)
                  + float(row['ra_s']) * (15.0 / 3600.0))
            sign = -1.0 if row['dec_sign'] == '-' else 1.0
            dec = sign * (float(row['dec_deg'])
                          + float(row['dec_min']) / 60.0
                          + float(row['dec_sec']) / 3600.0)
            dist = float(row['distance_au'])
            rows.append((jd, ra, dec, dist))
    return rows


def main():
    sun_csv = 'sun_transit_observations.csv'
    outfile = 'transit_observations.dat'

    sun = load_sun(sun_csv)
    if not sun:
        print("ERROR: no Sun observations loaded", file=sys.stderr)
        sys.exit(1)

    sun.sort(key=lambda r: r[0])
    jd_min = sun[0][0]
    jd_max = sun[-1][0]
    print(f"Sun: {len(sun)} observations, JD {jd_min:.1f} – {jd_max:.1f}")

    with open(outfile, 'w') as f:
        f.write('# Transit circle observations — USNO 6-inch, Washington D.C.\n')
        f.write(f'# Sun observations:  {len(sun)}\n')
        f.write(f'# JD range: {jd_min:.5f} – {jd_max:.5f}\n')
        f.write('#\n')
        f.write('# Columns:\n')
        f.write('#  1  JD_UTC          Julian date (UTC)\n')
        f.write('#  2  body            S=Sun\n')
        f.write('#  3  RA   (deg)      right ascension (apparent, of date)\n')
        f.write('#  4  Dec  (deg)      declination (apparent, of date)\n')
        f.write('#  5  dist (AU)       geocentric distance\n')
        for jd, ra, dec, dist in sun:
            f.write(f' {jd:17.5f}  S  {ra:12.6f}  {dec:+12.6f}  {dist:14.7E}\n')

    print(f"✅ Saved → {outfile}")


if __name__ == '__main__':
    main()
