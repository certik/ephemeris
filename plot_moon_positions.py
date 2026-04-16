"""
Plot Moon RA/Dec positions from USNO transit circle observations.

Reads moon_transit_observations.csv (extracted from USNO Pub. Vol. 23) and
generates a figure with:
  - Top panel:    RA vs time
  - Bottom panel: Dec vs time
  - Saved to:     moon_positions.png
"""

import csv
import os
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np


def jd_to_datetime(jd_ephem: float) -> datetime:
    """Convert JD ephemeris (JD - 2400000 approx) to Python datetime."""
    jd = jd_ephem + 2400000.0
    j2000_epoch = datetime(2000, 1, 1, 12, 0, 0)
    return j2000_epoch + timedelta(days=jd - 2451545.0)


def load_moon_data(csv_path: str):
    dates, ra_deg, dec_deg = [], [], []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row['ra_h'] or not row['ra_s'] or not row['dec_deg'] or not row['dec_sec']:
                continue
            jd = float(row['jd_ephem'])
            dates.append(jd_to_datetime(jd))

            ra = float(row['ra_h']) * 15.0 \
                + float(row['ra_m']) * (15.0 / 60.0) \
                + float(row['ra_s']) * (15.0 / 3600.0)
            ra_deg.append(ra)

            sign = -1.0 if row['dec_sign'] == '-' else 1.0
            dec = float(row['dec_deg']) \
                + float(row['dec_min']) / 60.0 \
                + float(row['dec_sec']) / 3600.0
            dec_deg.append(sign * dec)
    return np.array(dates), np.array(ra_deg), np.array(dec_deg)


def main():
    csv_path = 'moon_transit_observations.csv'
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"{csv_path} not found")

    dates, ra, dec = load_moon_data(csv_path)
    print(f"Loaded {len(dates):,} Moon observations ({dates[0]:%Y} – {dates[-1]:%Y})")

    fig, (ax_ra, ax_dec) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    fig.suptitle(
        f'Moon Transit Observations — USNO 6-inch, Washington D.C.\n'
        f'{len(dates):,} observations, {dates[0]:%Y}–{dates[-1]:%Y}',
        fontsize=14
    )

    dot_kw = dict(s=3, alpha=0.5, color='#2060c0', rasterized=True)

    ax_ra.scatter(dates, ra, **dot_kw)
    ax_ra.set_ylabel('Right Ascension (°)')
    ax_ra.set_ylim(0, 360)
    ax_ra.set_yticks(range(0, 361, 60))
    ax_ra.grid(True, alpha=0.3)

    ax_dec.scatter(dates, dec, **dot_kw)
    ax_dec.set_ylabel('Declination (°)')
    ax_dec.set_ylim(-32, 32)
    ax_dec.axhline(28.5, color='red', ls='--', lw=0.8, alpha=0.5, label='±28.5° (max)')
    ax_dec.axhline(-28.5, color='red', ls='--', lw=0.8, alpha=0.5)
    ax_dec.set_xlabel('Year')
    ax_dec.legend(loc='upper right', fontsize=9)
    ax_dec.grid(True, alpha=0.3)

    ax_dec.xaxis.set_major_locator(mdates.YearLocator())
    ax_dec.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.autofmt_xdate(rotation=30)

    fig.tight_layout()
    out = 'moon_positions.png'
    fig.savefig(out, dpi=150)
    print(f"✅ Saved → {out}")
    plt.close(fig)


if __name__ == '__main__':
    main()
