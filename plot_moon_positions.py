"""
Plot Moon RA/Dec positions from parsed USNO 6-inch transit observations.

Reads the CSV produced by parse_usno_transit.py and generates a figure with:
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


def jd_to_datetime(jd: float) -> datetime:
    j2000_epoch = datetime(2000, 1, 1, 12, 0, 0)
    return j2000_epoch + timedelta(days=jd - 2451545.0)


def load_moon_data(csv_path: str):
    dates, ra, dec = [], [], []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['planet_code'] != '010':
                continue
            dates.append(jd_to_datetime(float(row['jd'])))
            ra.append(float(row['ra_deg']))
            dec.append(float(row['dec_deg']))
    return np.array(dates), np.array(ra), np.array(dec)


def main():
    csv_path = 'trnstswash6a_parsed.csv'
    if not os.path.exists(csv_path):
        print(f"CSV not found, generating: {csv_path}")
        from parse_usno_transit import parse_file_to_csv
        parse_file_to_csv('trnstswash6a.txt', csv_path)

    dates, ra, dec = load_moon_data(csv_path)
    print(f"Loaded {len(dates):,} Moon observations ({dates[0]:%Y} – {dates[-1]:%Y})")

    fig, (ax_ra, ax_dec) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    fig.suptitle(
        f'Moon Transit Observations — USNO 6-inch, Washington D.C.\n'
        f'{len(dates):,} observations, {dates[0]:%Y}–{dates[-1]:%Y}',
        fontsize=14
    )

    dot_kw = dict(s=1.5, alpha=0.4, color='#2060c0', rasterized=True)

    ax_ra.scatter(dates, ra, **dot_kw)
    ax_ra.set_ylabel('Right Ascension (°)')
    ax_ra.set_ylim(0, 360)
    ax_ra.set_yticks(range(0, 361, 60))
    ax_ra.grid(True, alpha=0.3)

    ax_dec.scatter(dates, dec, **dot_kw)
    ax_dec.set_ylabel('Declination (°)')
    ax_dec.set_xlabel('Year')
    ax_dec.grid(True, alpha=0.3)

    ax_dec.xaxis.set_major_locator(mdates.YearLocator(5))
    ax_dec.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.autofmt_xdate(rotation=30)

    fig.tight_layout()
    out = 'moon_positions.png'
    fig.savefig(out, dpi=150)
    print(f"✅ Saved → {out}")
    plt.close(fig)


if __name__ == '__main__':
    main()
