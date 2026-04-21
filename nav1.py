#!/usr/bin/env python3
"""
Celestial Navigation Fix Calculator: Latitude from Polaris altitude + Longitude from Jupiter altitude
using Astropy for daily planetary ephemerides (builtin low-precision ephemeris, sufficient for navigation).
Accepts Polaris altitude (≈ latitude in Northern Hemisphere) and Jupiter observed altitude at given UTC.
Outputs latitude and the two mathematically possible longitudes (western/eastern).
"""

import argparse
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body, solar_system_ephemeris
import astropy.units as u
from astropy.coordinates import Angle

def parse_dms(dms_str):
    """Parse DMS string like '35d53m', '35d 53m', '35 53', or decimal degrees to float degrees."""
    # Normalize input for astropy Angle parser
    dms_str = dms_str.replace("deg", "d").replace("'", "m").replace('"', "").replace(" ", "")
    if "d" not in dms_str and "m" not in dms_str:
        return float(dms_str)  # already decimal
    return Angle(dms_str).deg

def main():
    parser = argparse.ArgumentParser(
        description="Compute exact observer latitude and longitude from Polaris & Jupiter altitudes "
                    "using daily ephemerides (Astropy).",
        epilog="Example: python script.py --date 2026-04-21 --time 05:43:00 --alt_polaris 35d53m --alt_jupiter 20d29m"
    )
    parser.add_argument('--date', required=True, help='UTC date in YYYY-MM-DD format')
    parser.add_argument('--time', required=True, help='UTC time in HH:MM:SS format (e.g. 05:43:00)')
    parser.add_argument('--alt_polaris', required=True, help='Polaris observed altitude (e.g. "35d53m" or 35.8833)')
    parser.add_argument('--alt_jupiter', required=True, help='Jupiter observed altitude (e.g. "20d29m" or 20.4833)')
    args = parser.parse_args()

    # Parse inputs
    utc_str = f"{args.date} {args.time}"
    t = Time(utc_str, scale='utc')
    lat = parse_dms(args.alt_polaris)          # Latitude ≈ Polaris altitude (Northern Hemisphere)
    alt_j = parse_dms(args.alt_jupiter)

    print(f"Input UTC: {utc_str}")
    print(f"Latitude (from Polaris alt): {lat:.4f}° N")
    print(f"Jupiter observed altitude: {alt_j:.4f}°")

    # Get Jupiter geocentric position from daily ephemerides
    with solar_system_ephemeris.set('builtin'):
        jup = get_body('jupiter', t)
        ra = jup.ra
        dec = jup.dec.deg

    # Greenwich Mean Sidereal Time (GMST) and GHA of Jupiter
    gmst = t.sidereal_time('mean', 'greenwich')
    gha = (gmst - ra).wrap_at(360 * u.deg).deg

    print(f"Jupiter (ephemeris): Dec = {dec:.4f}°, GHA = {gha:.4f}°")

    # Spherical trigonometry: solve for Local Hour Angle (LHA)
    # cos(LHA) = [sin(alt) - sin(lat)·sin(dec)] / [cos(lat)·cos(dec)]
    sa = np.sin(np.deg2rad(alt_j))
    sl = np.sin(np.deg2rad(lat))
    sd = np.sin(np.deg2rad(dec))
    cl = np.cos(np.deg2rad(lat))
    cd = np.cos(np.deg2rad(dec))

    cos_lha = (sa - sl * sd) / (cl * cd)
    cos_lha = np.clip(cos_lha, -1.0, 1.0)  # numerical safety

    lha1 = np.rad2deg(np.arccos(cos_lha))
    lha2 = 360.0 - lha1

    # Longitude = GHA - LHA (mod 360); convert to E/W
    lon1 = (gha - lha1) % 360
    lon2 = (gha - lha2) % 360

    def format_lon(lon_deg):
        if lon_deg > 180:
            return f"{360 - lon_deg:.4f}° E"
        else:
            return f"{lon_deg:.4f}° W"

    print("\n=== Computed Positions ===")
    print(f"1. Longitude (LHA = {lha1:.2f}°): {format_lon(lon1)}")
    print(f"2. Longitude (LHA = {lha2:.2f}°): {format_lon(lon2)}")
    print("\nNote: The western longitude is usually the practical nighttime solution; the eastern is the alternate daytime fix.")
    print("      Polaris altitude is used directly as latitude (standard approximation; tables add <1° correction if needed).")

if __name__ == "__main__":
    main()
