"""
Data from https://ssd.jpl.nasa.gov/planets/obs_data.html

USNO 6″, 1913-76
https://ssd.jpl.nasa.gov/dat/planets/trnstswash6a.txt

Fixed-width 80-column format per O'Handley (1968),
JPL Technical Report 32-1296.
https://ssd.jpl.nasa.gov/dat/planets/optical_format.txt
"""

import csv
import os
from datetime import datetime, timedelta

# USNO Washington, D.C. coordinates (observatory code 7866)
USNO_LAT_DEG = 38 + 55/60 + 17/3600.0   # ≈ 38.9214° N
USNO_LON_DEG = -(77 + 4/60 + 1/3600.0)  # ≈ -77.0669° W


def jd_to_date(jd: float) -> str:
    """Convert Julian Date to ISO calendar date string."""
    j2000_epoch = datetime(2000, 1, 1, 12, 0, 0)
    delta = timedelta(days=jd - 2451545.0)
    dt = j2000_epoch + delta
    return dt.strftime('%Y-%m-%d %H:%M:%S')


def _field_int(s: str) -> int:
    """Parse a fixed-width numeric sub-field, treating spaces as zeros."""
    return int(s.replace(' ', '0') or '0')


def parse_ra(raw: str) -> float | None:
    """Parse RA from hhmmss.sss (9-char fixed-width) to decimal degrees."""
    if len(raw) != 9:
        return None
    try:
        hh = _field_int(raw[0:2])
        mm = _field_int(raw[2:4])
        ss = _field_int(raw[4:6]) + _field_int(raw[6:9]) / 1000.0
        return (hh + mm / 60.0 + ss / 3600.0) * 15.0
    except ValueError:
        return None


def parse_dec(raw: str) -> float | None:
    """Parse Dec from ±ddmmss.ss (9-char fixed-width) to decimal degrees."""
    if len(raw) != 9:
        return None
    try:
        sign = -1 if raw[0] == '-' else 1
        dd = _field_int(raw[1:3])
        mm = _field_int(raw[3:5])
        ss = _field_int(raw[5:7]) + _field_int(raw[7:9]) / 100.0
        return sign * (dd + mm / 60.0 + ss / 3600.0)
    except ValueError:
        return None


def parse_usno_transit_line(line: str):
    """Parse one 80-column fixed-width transit observation record."""
    if len(line.rstrip()) < 80 or line.startswith('#'):
        return None

    # Fixed-width columns (1-based per format spec → Python 0-based slicing)
    planet_code = line[1:4].strip()       # cols  2-4
    jd5_raw     = line[4:16].strip()      # cols  5-16  (JD × 10^5)
    observatory = line[21:25].strip()     # cols 22-25
    obs_type    = line[28:29].strip()     # col  29
    ra_raw      = line[33:42]             # cols 34-42  (hhmmss.sss)
    ra_equinox  = line[42:43].strip()     # col  43
    ra_residual = line[43:49].strip()     # cols 44-49
    dec_raw     = line[50:59]             # cols 51-59  (±ddmmss.ss)
    dec_equinox = line[59:60].strip()     # col  60
    dec_residual = line[60:65].strip()    # cols 61-65
    year        = line[72:76].strip()     # cols 73-76
    source      = line[76:80].strip()     # cols 77-80

    try:
        jd = int(jd5_raw) / 100000.0
    except ValueError:
        return None

    ra_deg = parse_ra(ra_raw)
    dec_deg = parse_dec(dec_raw)
    if ra_deg is None or dec_deg is None:
        return None

    return {
        'planet_code': planet_code,
        'jd': round(jd, 5),
        'date_utc': jd_to_date(jd),
        'observatory': observatory,
        'obs_type': obs_type,
        'ra_deg': round(ra_deg, 6),
        'dec_deg': round(dec_deg, 6),
        'ra_equinox': ra_equinox,
        'dec_equinox': dec_equinox,
        'ra_residual': ra_residual,
        'dec_residual': dec_residual,
        'year': year,
        'source': source,
        'latitude_deg': round(USNO_LAT_DEG, 6),
        'longitude_deg': round(USNO_LON_DEG, 6),
    }


def parse_file_to_csv(input_file: str, output_csv: str = None):
    if output_csv is None:
        base = os.path.splitext(input_file)[0]
        output_csv = f"{base}_parsed.csv"

    records = []
    skipped = 0
    with open(input_file, 'r', encoding='ascii') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            rec = parse_usno_transit_line(line)
            if rec:
                records.append(rec)
            else:
                skipped += 1

    if records:
        fieldnames = records[0].keys()
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        print(f"✅ Parsed {len(records):,} observations → {output_csv}")
        if skipped:
            print(f"⚠️  Skipped {skipped} unparseable lines")
    else:
        print("No valid records found.")

    return records


# ====================== EXAMPLE USAGE ======================
if __name__ == "__main__":
    filename = "trnstswash6a.txt"
    records = parse_file_to_csv(filename)

    for rec in records[:3]:
        print(rec)
        print(f"→ Planet {rec['planet_code']}, {rec['date_utc']}, "
              f"RA {rec['ra_deg']:.4f}°, Dec {rec['dec_deg']:.4f}°")
        print("-" * 80)
