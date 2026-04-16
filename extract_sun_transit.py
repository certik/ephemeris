#!/usr/bin/env python3
"""
Extract Sun transit circle observations (1963-1971) from scanned USNO Publication Vol. 23.

Uses macOS Vision framework (via PyObjC) to OCR pages 207-220 of the PDF
'1982PUSNO__23__165H.pdf'. The approach:
  1. Extract page images from PDF at 400 DPI
  2. For each page, OCR the JD column strip to find row positions
  3. For each row, OCR a horizontal strip to get all column values
  4. Clean OCR artifacts, parse RA/Dec, and write structured CSV

Output: sun_transit_observations.csv with columns:
  greenwich_date, jd_ephem, observer, clamp,
  ra_h, ra_m, ra_s, ra_o_minus_c,
  dec_sign, dec_deg, dec_min, dec_sec, dec_o_minus_c,
  distance_au

Requires: macOS with Vision framework, pymupdf, pillow, pyobjc-framework-vision
"""

import csv
import os
import re
import sys
import tempfile

import fitz  # PyMuPDF
from PIL import Image

# Only import Vision/Quartz at module level (macOS-only)
import Quartz
import Vision
from Foundation import NSURL

# --- Configuration ---
PDF_FILE = "1982PUSNO__23__165H.pdf"
FIRST_PAGE = 207  # PDF page number (document numbering)
LAST_PAGE = 220
PDF_PAGE_OFFSET = 165  # PDF page 0 = document page 165
DPI = 400
OUTPUT_CSV = "sun_transit_observations.csv"

# Column x-boundaries (normalized 0-1 of full page width)
COL_BOUNDS = {
    "date": (0.04, 0.175),
    "jd": (0.175, 0.270),
    "obsr": (0.270, 0.310),
    "cl": (0.310, 0.342),
    "ra": (0.342, 0.440),
    "ra_oc": (0.440, 0.498),
    "dec": (0.498, 0.600),
    "dec_oc": (0.600, 0.660),
    "dist": (0.660, 0.730),
}

# Data region (normalized y from top of page)
Y_START = 0.12
Y_END = 0.88


# --- OCR helpers ---

def ocr_image_file(path):
    """OCR an image file using macOS Vision. Returns [(x_center, y_center, text)]."""
    url = NSURL.fileURLWithPath_(path)
    source = Quartz.CGImageSourceCreateWithURL(url, None)
    cg_image = Quartz.CGImageSourceCreateImageAtIndex(source, 0, None)

    request = Vision.VNRecognizeTextRequest.alloc().init()
    request.setRecognitionLevel_(Vision.VNRequestTextRecognitionLevelAccurate)
    request.setUsesLanguageCorrection_(False)

    handler = Vision.VNImageRequestHandler.alloc().initWithCGImage_options_(cg_image, {})
    handler.performRequests_error_([request], None)

    items = []
    for obs in request.results():
        text = obs.topCandidates_(1)[0].string()
        bbox = obs.boundingBox()
        x_center = bbox.origin.x + bbox.size.width / 2
        y_center = 1.0 - bbox.origin.y - bbox.size.height / 2
        items.append((x_center, y_center, text))
    return items


def get_jd_positions(img):
    """OCR the JD column strip to find row y-positions and JD values."""
    W, H = img.size
    x1, x2 = int(0.18 * W), int(0.27 * W)
    y1, y2 = int(Y_START * H), int(Y_END * H)
    crop = img.crop((x1, y1, x2, y2))

    tmp = tempfile.mktemp(suffix=".png")
    crop.save(tmp)
    items = ocr_image_file(tmp)
    os.unlink(tmp)

    crop_h = y2 - y1
    positions = []
    for _, y_local, text in items:
        # Skip header text
        if text in ("24",) or "Eph" in text or "Jul" in text:
            continue
        # Convert y from crop-local to full-image normalized
        y_pixel = y1 + y_local * crop_h
        y_norm = y_pixel / H
        positions.append((y_norm, text))

    positions.sort(key=lambda p: p[0])
    return positions


def ocr_row_strip(img, y_norm, strip_half_height=18):
    """OCR a horizontal strip at the given y-position. Returns {col: text}."""
    W, H = img.size
    y_pixel = int(y_norm * H)
    strip_y1 = max(0, y_pixel - strip_half_height)
    strip_y2 = min(H, y_pixel + strip_half_height)

    x_left = int(0.04 * W)
    x_right = int(0.73 * W)
    strip = img.crop((x_left, strip_y1, x_right, strip_y2))

    tmp = tempfile.mktemp(suffix=".png")
    strip.save(tmp)
    items = ocr_image_file(tmp)
    os.unlink(tmp)

    # Map x from crop-local back to full-page coordinates
    x_offset = 0.04
    x_scale = 0.73 - 0.04

    result = {}
    for x_local, _, text in items:
        x_page = x_offset + x_local * x_scale
        for col, (lo, hi) in COL_BOUNDS.items():
            if lo <= x_page < hi:
                if col in result:
                    result[col].append((x_page, text))
                else:
                    result[col] = [(x_page, text)]
                break

    # Merge text within each column (sorted by x)
    merged = {}
    for col, entries in result.items():
        entries.sort(key=lambda e: e[0])
        merged[col] = " ".join(t for _, t in entries)
    return merged


# --- Cleaning helpers ---

def clean_artifacts(s):
    """Remove OCR artifacts and normalize text."""
    if not s:
        return s
    s = s.replace("\u2022", " ").replace("'", " ").replace('"', " ")
    s = s.replace(",", " ").replace(";", " ")
    s = s.replace("$", "+").replace("*", "+").replace("#", "")
    # Cyrillic lookalikes
    s = s.replace("\u0417", "3")  # З -> 3
    s = s.replace("\u041e", "0")  # О -> 0
    s = s.replace("\u0435", "e")  # е -> e
    s = re.sub(r"\s+", " ", s).strip()
    return s


def clean_ra(ra_str):
    """Clean and fix RA string."""
    ra = clean_artifacts(ra_str)
    if not ra:
        return ra
    ra = re.sub(r"[^0-9.\s]", "", ra).strip()
    ra = re.sub(r"\s+", " ", ra)
    parts = ra.split()
    if len(parts) == 2:
        # Missing hour — prepend 0 (Sun near RA=0h)
        try:
            mm = int(parts[0])
            if 0 <= mm <= 59:
                ra = "0 " + ra
        except ValueError:
            pass
    return ra


def clean_signed(s):
    """Clean a signed decimal value like +0.035 or -1.31."""
    s = clean_artifacts(s)
    if not s:
        return s
    s = re.sub(r"[^0-9.\s+\-]", "", s).strip()
    s = re.sub(r"([+-])\s+", r"\1", s)
    return s


def clean_distance(s):
    """Clean distance value."""
    s = clean_artifacts(s)
    if not s:
        return s
    return re.sub(r"[^0-9.]", "", s)


# --- Main pipeline ---

def extract_pages(pdf_path):
    """Extract page images from PDF at high DPI."""
    doc = fitz.open(pdf_path)
    pages = {}
    for page_num in range(FIRST_PAGE, LAST_PAGE + 1):
        idx = page_num - PDF_PAGE_OFFSET
        page = doc[idx]
        mat = fitz.Matrix(DPI / 72, DPI / 72)
        pix = page.get_pixmap(matrix=mat)
        tmp = tempfile.mktemp(suffix=".png")
        pix.save(tmp)
        pages[page_num] = tmp
    doc.close()
    return pages


def extract_all_rows(page_images):
    """OCR all pages and return raw row dicts."""
    all_rows = []
    for page_num in sorted(page_images.keys()):
        path = page_images[page_num]
        print(f"  Processing page {page_num}...", file=sys.stderr)
        img = Image.open(path)
        jd_positions = get_jd_positions(img)

        rows = []
        for y_norm, jd_text in jd_positions:
            row_data = ocr_row_strip(img, y_norm)
            row_data["jd"] = jd_text  # use reliable column-OCR JD
            rows.append(row_data)

        print(f"    {len(rows)} rows", file=sys.stderr)
        all_rows.extend(rows)
    return all_rows


def parse_ra(ra_str):
    """Parse RA string into (h, m, s) or None."""
    parts = ra_str.split() if ra_str else []
    if len(parts) != 3:
        return None
    try:
        h, m, s = int(parts[0]), int(parts[1]), float(parts[2])
        if 0 <= h <= 23 and 0 <= m <= 59 and 0 <= s < 60:
            return (str(h), str(m), parts[2])
    except ValueError:
        pass
    return None


def parse_dec(dec_str):
    """Parse Dec string into (sign, deg, min, sec) or partial."""
    dec = clean_artifacts(dec_str)
    if not dec:
        return None
    dec = re.sub(r"[^0-9.\s+\-]", "", dec).strip()

    sign = "+"
    if dec.startswith("-"):
        sign = "-"
        dec = dec[1:].strip()
    elif dec.startswith("+"):
        dec = dec[1:].strip()

    parts = dec.split()
    if len(parts) == 3:
        try:
            d, m, s = int(parts[0]), int(parts[1]), float(parts[2])
            if 0 <= d <= 90 and 0 <= m <= 59 and 0 <= s < 60:
                return (sign, str(d), str(m), parts[2])
        except ValueError:
            pass
    elif len(parts) == 2:
        try:
            d, m = int(parts[0]), int(parts[1])
            if 0 <= d <= 90 and 0 <= m <= 59:
                return (sign, str(d), str(m), "")
        except ValueError:
            pass
    return None


def clean_jd(jd_str):
    """Clean and validate JD value."""
    jd = clean_artifacts(jd_str)
    if not jd:
        return None
    # Remove trailing non-numeric chars
    jd = re.sub(r"[^0-9.]", "", jd)
    if not jd:
        return None
    # Fix missing dot: "3910420881" -> "39104.20881"
    if "." not in jd and len(jd) == 10:
        jd = jd[:5] + "." + jd[5:]
    # Validate range (should be ~38000-42000 for 1963-1971)
    try:
        v = float(jd)
        if 38000 < v < 42000:
            return jd
    except ValueError:
        pass
    return None


def clean_and_write(raw_rows, output_path):
    """Clean OCR output and write structured CSV."""
    fieldnames = [
        "greenwich_date", "jd_ephem", "observer", "clamp",
        "ra_h", "ra_m", "ra_s", "ra_o_minus_c",
        "dec_sign", "dec_deg", "dec_min", "dec_sec", "dec_o_minus_c",
        "distance_au",
    ]

    current_year = None
    current_month = None
    n_good_ra = 0
    n_good_dec = 0
    clean_rows = []

    for row in raw_rows:
        # Validate and clean JD — skip rows with invalid JD
        jd = clean_jd(row.get("jd", ""))
        if jd is None:
            continue

        # Date tracking
        date_str = row.get("date", "")
        year_match = re.search(r"(\d{4})", date_str)
        if year_match:
            current_year = year_match.group(1)
        month_match = re.search(
            r"(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)", date_str, re.IGNORECASE
        )
        if month_match:
            current_month = month_match.group(1).capitalize()
        day_match = re.search(r"(\d+\.\d+)", date_str)
        day_str = day_match.group(1) if day_match else ""

        full_date = ""
        if current_year and current_month and day_str:
            full_date = f"{current_year} {current_month}. {day_str}"

        # Parse RA
        ra = clean_ra(row.get("ra", ""))
        ra_parsed = parse_ra(ra)
        ra_h = ra_m = ra_s = ""
        if ra_parsed:
            ra_h, ra_m, ra_s = ra_parsed
            n_good_ra += 1

        # Parse Dec
        dec_parsed = parse_dec(row.get("dec", ""))
        dec_sign = dec_deg = dec_min = dec_sec = ""
        if dec_parsed:
            dec_sign, dec_deg, dec_min = dec_parsed[0], dec_parsed[1], dec_parsed[2]
            dec_sec = dec_parsed[3] if len(dec_parsed) > 3 else ""
            n_good_dec += 1

        clean_rows.append({
            "greenwich_date": full_date,
            "jd_ephem": jd,
            "observer": clean_artifacts(row.get("obsr", "")),
            "clamp": clean_artifacts(row.get("cl", "")),
            "ra_h": ra_h,
            "ra_m": ra_m,
            "ra_s": ra_s,
            "ra_o_minus_c": clean_signed(row.get("ra_oc", "")),
            "dec_sign": dec_sign,
            "dec_deg": dec_deg,
            "dec_min": dec_min,
            "dec_sec": dec_sec,
            "dec_o_minus_c": clean_signed(row.get("dec_oc", "")),
            "distance_au": clean_distance(row.get("dist", "")),
        })

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in clean_rows:
            writer.writerow(r)

    return clean_rows, n_good_ra, n_good_dec


def print_summary(clean_rows, n_good_ra, n_good_dec):
    """Print extraction summary."""
    n = len(clean_rows)
    n_jd = sum(1 for r in clean_rows if r["jd_ephem"])
    n_date = sum(1 for r in clean_rows if r["greenwich_date"])
    n_dist = sum(1 for r in clean_rows if r["distance_au"])

    print("\n" + "=" * 70)
    print("SUN TRANSIT CIRCLE OBSERVATIONS — EXTRACTION SUMMARY")
    print("=" * 70)
    print(f"  Source:  USNO Publication Vol. 23 (1982), pages {FIRST_PAGE}-{LAST_PAGE}")
    print(f"  Period:  1963 Mar – 1971 Dec")
    print(f"  Output:  {OUTPUT_CSV}")
    print(f"\n  Total observations:  {n:,}")
    print(f"  Complete JD:         {n_jd:,} ({100*n_jd/n:.1f}%)")
    print(f"  Complete dates:      {n_date:,} ({100*n_date/n:.1f}%)")
    print(f"  Parsed RA (h m s):   {n_good_ra:,} ({100*n_good_ra/n:.1f}%)")
    print(f"  Parsed Dec (° ′ ″): {n_good_dec:,} ({100*n_good_dec/n:.1f}%)")
    print(f"  Distance (AU):       {n_dist:,} ({100*n_dist/n:.1f}%)")

    # JD range
    jds = []
    for r in clean_rows:
        try:
            jds.append(float(r["jd_ephem"]))
        except ValueError:
            pass
    if jds:
        jds_valid = sorted(jds)
        print(f"\n  JD range: {jds_valid[0]:.5f} – {jds_valid[-1]:.5f}")
        print(f"  ({jds_valid[-1] - jds_valid[0]:.0f} days span)")

    print("\n" + "-" * 70)
    print("  CSV COLUMNS")
    print("-" * 70)
    print("  greenwich_date    Greenwich date (YYYY Mon. DD.dd)")
    print("  jd_ephem          Julian Ephemeris Date (24xxxxx.xxxxx)")
    print("  observer          Observer code (2 chars)")
    print("  clamp             Clamp position (E or W)")
    print("  ra_h, ra_m, ra_s  Right Ascension (hours, minutes, seconds)")
    print("  ra_o_minus_c      RA observed minus computed (arcsec)")
    print("  dec_sign/deg/min/sec  Declination (+/-, degrees, minutes, seconds)")
    print("  dec_o_minus_c     Dec observed minus computed (arcsec)")
    print("  distance_au       Earth-Sun distance (AU)")
    print("=" * 70)


def main():
    if not os.path.exists(PDF_FILE):
        print(f"Error: {PDF_FILE} not found", file=sys.stderr)
        sys.exit(1)

    print("Extracting Sun transit circle observations from USNO Pub. Vol. 23...",
          file=sys.stderr)
    print(f"  Pages {FIRST_PAGE}-{LAST_PAGE}, {DPI} DPI\n", file=sys.stderr)

    # Step 1: Extract page images
    print("Step 1: Extracting page images from PDF...", file=sys.stderr)
    page_images = extract_pages(PDF_FILE)

    # Step 2: OCR all pages
    print("\nStep 2: OCR extraction (column-strip + row-strip)...", file=sys.stderr)
    raw_rows = extract_all_rows(page_images)

    # Cleanup temp images
    for path in page_images.values():
        os.unlink(path)

    # Step 3: Clean and write CSV
    print(f"\nStep 3: Cleaning and writing {OUTPUT_CSV}...", file=sys.stderr)
    clean_rows, n_good_ra, n_good_dec = clean_and_write(raw_rows, OUTPUT_CSV)

    print_summary(clean_rows, n_good_ra, n_good_dec)


if __name__ == "__main__":
    main()
