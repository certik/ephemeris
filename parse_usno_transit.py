import re
from datetime import timedelta
import csv
import os

# USNO Washington, D.C. coordinates (fixed for all these observations)
USNO_LAT_DEG = 38 + 55/60 + 17/3600.0   # ≈ 38.9213889° N
USNO_LON_DEG = -(77 + 4/60 + 1/3600.0)  # ≈ -77.0669444° W

def parse_usno_transit_line(line: str):
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    # Split on any whitespace (the files use variable spacing)
    parts = re.split(r'\s+', line)
    if len(parts) < 7:
        return None

    obs_id = parts[0]
    category = parts[1]          # usually "1"
    body_code = parts[2]         # e.g. "003" — maps to Sun or Moon inside each file
    time_raw = parts[3]          # e.g. "058084000" → milliseconds past reference (≈0h)
    quality = parts[4]
    zd_raw = parts[5]            # zenith distance / circle reading in 0.001 arcsec
    label = parts[6]             # e.g. "1913USNO"

    # Convert time
    time_ms = int(time_raw)
    time_seconds = time_ms / 1000.0
    time_hms = str(timedelta(seconds=time_seconds))

    # Convert zenith distance
    zd_arcsec = int(zd_raw) / 1000.0
    zd_deg = zd_arcsec / 3600.0

    # Approximate declination (upper culmination, northern hemisphere)
    # Real JPL reduction also corrects for refraction, semi-diameter, etc.
    approx_dec_deg = USNO_LAT_DEG - zd_deg

    return {
        'obs_id': obs_id,
        'year_label': label,
        'body_code': body_code,
        'time_seconds_day': round(time_seconds, 3),
        'time_hms': time_hms,
        'zd_arcsec': round(zd_arcsec, 3),
        'zd_deg': round(zd_deg, 6),
        'approx_dec_deg': round(approx_dec_deg, 6),
        'quality': quality,
        'latitude_deg': round(USNO_LAT_DEG, 6),
        'longitude_deg': round(USNO_LON_DEG, 6)
    }


def parse_file_to_csv(input_file: str, output_csv: str = None):
    if output_csv is None:
        base = os.path.splitext(input_file)[0]
        output_csv = f"{base}_parsed.csv"

    records = []
    with open(input_file, 'r', encoding='ascii') as f:
        for line in f:
            rec = parse_usno_transit_line(line)
            if rec:
                records.append(rec)

    # Write clean CSV
    if records:
        fieldnames = records[0].keys()
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        print(f"✅ Parsed {len(records):,} observations → {output_csv}")
    else:
        print("No valid records found.")

    return records


# ====================== EXAMPLE USAGE ======================
if __name__ == "__main__":
    filename = "trnstswash6a.txt"          # ← change to your downloaded file
    records = parse_file_to_csv(filename)

    # Quick preview of first few lines
    for rec in records[:3]:
        print(rec)
        print(f"→ Body code {rec['body_code']}, Time {rec['time_hms']}, "
              f"ZD {rec['zd_deg']:.4f}°, Approx Dec {rec['approx_dec_deg']:.4f}°")
        print("-" * 80)
