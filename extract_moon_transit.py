#!/usr/bin/env python3
"""
Extract Moon transit circle observations (1963-1971) from scanned USNO Pub. Vol. 23.

Uses macOS Vision framework (via PyObjC) to OCR pages 221-230 of the PDF
'1982PUSNO__23__165H.pdf'.  Same column-strip + row-strip approach as
extract_sun_transit.py, but adapted for the Moon table which has:
  - Limb number (1/2) column
  - RA limb correction + limb side (N/S)
  - Dec limb correction
  - No distance column

Output: moon_transit_observations.csv
"""

import csv
import os
import re
import sys
import tempfile

import fitz  # PyMuPDF
from PIL import Image

import Quartz
import Vision
from Foundation import NSURL

# --- Configuration ---
PDF_FILE = "1982PUSNO__23__165H.pdf"
FIRST_PAGE = 221
LAST_PAGE = 230
PDF_PAGE_OFFSET = 165
DPI = 400
OUTPUT_CSV = "moon_transit_observations.csv"

# Column x-boundaries (normalized 0-1 of full page width)
# Moon table is wider than Sun table with extra Limb columns
COL_BOUNDS = {
    "lun":      (0.02,  0.060),
    "date":     (0.060, 0.215),
    "jd":       (0.230, 0.310),
    "obsr":     (0.310, 0.348),
    "cl":       (0.348, 0.378),
    "limb_num": (0.378, 0.410),
    "ra":       (0.410, 0.485),
    "ra_oc":    (0.485, 0.560),
    "ra_limb_corr": (0.560, 0.605),
    "limb_side": (0.605, 0.640),
    "dec":      (0.640, 0.725),
    "dec_oc":   (0.725, 0.795),
    "dec_limb_corr": (0.795, 0.890),
}

Y_START = 0.12
Y_END = 0.88

# Row grouping tolerance (normalized y distance to consider same row)
ROW_TOL = 0.005


# --- OCR helpers ---

def ocr_image_file(path):
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


def ocr_full_page(img):
    """OCR the full data region of a page and return items in page coordinates."""
    W, H = img.size
    y1, y2 = int(Y_START * H), int(Y_END * H)
    data_region = img.crop((0, y1, W, y2))

    tmp = tempfile.mktemp(suffix=".png")
    data_region.save(tmp)
    items = ocr_image_file(tmp)
    os.unlink(tmp)

    crop_h = y2 - y1
    page_items = []
    for x, y_local, text in items:
        y_page = (y1 + y_local * crop_h) / H
        page_items.append((x, y_page, text))
    return page_items


def group_into_rows(items):
    """Group OCR text fragments into rows by y-position proximity.
    Uses mean-based clustering to prevent chaining across adjacent rows."""
    if not items:
        return []
    items_sorted = sorted(items, key=lambda i: i[1])

    rows = []
    current_row = [items_sorted[0]]
    current_y_sum = items_sorted[0][1]
    for item in items_sorted[1:]:
        current_y_mean = current_y_sum / len(current_row)
        if abs(item[1] - current_y_mean) < ROW_TOL:
            current_row.append(item)
            current_y_sum += item[1]
        else:
            rows.append(current_row)
            current_row = [item]
            current_y_sum = item[1]
    rows.append(current_row)
    return rows


def row_to_columns(row_items):
    """Assign text fragments in a row to columns by x-position."""
    result = {}
    for x, _, text in row_items:
        for col, (lo, hi) in COL_BOUNDS.items():
            if lo <= x < hi:
                if col in result:
                    result[col].append((x, text))
                else:
                    result[col] = [(x, text)]
                break

    merged = {}
    for col, entries in result.items():
        entries.sort(key=lambda e: e[0])
        merged[col] = " ".join(t for _, t in entries)
    return merged


# --- Cleaning helpers ---

def clean_artifacts(s):
    if not s:
        return s
    s = s.replace("\u2022", " ").replace("'", " ").replace('"', " ")
    s = s.replace(",", " ").replace(";", " ")
    s = s.replace("$", "+").replace("*", "+").replace("#", "")
    s = s.replace("\u0417", "3")
    s = s.replace("\u041e", "0")
    s = s.replace("\u0435", "e")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def clean_ra(ra_str):
    ra = clean_artifacts(ra_str)
    if not ra:
        return ra
    ra = re.sub(r"[^0-9.\s]", "", ra).strip()
    ra = re.sub(r"\s+", " ", ra)
    return ra


def clean_signed(s):
    s = clean_artifacts(s)
    if not s:
        return s
    s = re.sub(r"[^0-9.\s+\-]", "", s).strip()
    s = re.sub(r"([+-])\s+", r"\1", s)
    return s


def clean_jd(jd_str):
    jd = clean_artifacts(jd_str)
    if not jd:
        return None
    jd = re.sub(r"[^0-9.]", "", jd)
    if not jd:
        return None
    if "." not in jd and len(jd) == 10:
        jd = jd[:5] + "." + jd[5:]
    try:
        v = float(jd)
        if 38000 < v < 42000:
            return jd
    except ValueError:
        pass
    return None


# --- Main pipeline ---

def extract_pages(pdf_path):
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
    all_rows = []
    for page_num in sorted(page_images.keys()):
        path = page_images[page_num]
        print(f"  Processing page {page_num}...", file=sys.stderr)
        img = Image.open(path)
        W, H = img.size

        # Step 1: Full-page OCR to find row y-positions from JD values
        page_items = ocr_full_page(img)
        jd_rows = {}
        for x, y, text in page_items:
            # JD values are ~11-digit numbers with decimal at x~0.265
            if 0.23 < x < 0.31 and re.search(r'\d{4,5}\.\d{3,}', text):
                jd_rows[y] = text

        # Merge nearby y-positions (within ROW_TOL)
        y_positions = sorted(jd_rows.keys())
        merged_positions = []
        for y in y_positions:
            if merged_positions and abs(y - merged_positions[-1][0]) < ROW_TOL:
                continue
            merged_positions.append((y, jd_rows[y]))

        # Step 2: Row-strip OCR for each detected position
        rows = []
        for y_norm, jd_text in merged_positions:
            row_data = ocr_row_strip(img, y_norm)
            row_data["jd"] = jd_text  # use full-page JD (more reliable)
            rows.append(row_data)

        print(f"    {len(rows)} rows", file=sys.stderr)
        all_rows.extend(rows)
    return all_rows


def ocr_row_strip(img, y_norm, strip_half_height=22):
    """OCR a horizontal strip at the given y-position. Returns {col: text}."""
    W, H = img.size
    y_pixel = int(y_norm * H)
    strip_y1 = max(0, y_pixel - strip_half_height)
    strip_y2 = min(H, y_pixel + strip_half_height)

    x_left = int(0.02 * W)
    x_right = int(0.89 * W)
    strip = img.crop((x_left, strip_y1, x_right, strip_y2))

    tmp = tempfile.mktemp(suffix=".png")
    strip.save(tmp)
    items = ocr_image_file(tmp)
    os.unlink(tmp)

    x_offset = 0.02
    x_scale = 0.89 - 0.02

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

    merged = {}
    for col, entries in result.items():
        entries.sort(key=lambda e: e[0])
        merged[col] = " ".join(t for _, t in entries)
    return merged


def parse_ra(ra_str):
    parts = ra_str.split() if ra_str else []
    if len(parts) == 3:
        try:
            h, m, s = int(parts[0]), int(parts[1]), float(parts[2])
            if 0 <= h <= 23 and 0 <= m <= 59 and 0 <= s < 60:
                return (str(h), str(m), parts[2])
        except ValueError:
            pass
    elif len(parts) == 2:
        # Missing hours digit — store minutes and seconds, leave hours empty
        try:
            m, s = int(parts[0]), float(parts[1])
            if 0 <= m <= 59 and 0 <= s < 60:
                return ("", str(m), parts[1])
        except ValueError:
            pass
    return None


def parse_dec(dec_str):
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


def clean_limb_side(s):
    """Extract N or S from limb side text."""
    s = clean_artifacts(s)
    if not s:
        return ""
    if "N" in s.upper():
        return "N"
    if "S" in s.upper():
        return "S"
    return s


def clean_limb_num(s):
    """Extract limb number (1 or 2)."""
    s = clean_artifacts(s)
    if not s:
        return ""
    s = re.sub(r"[^12]", "", s)
    return s[:1] if s else ""


def clean_and_write(raw_rows, output_path):
    fieldnames = [
        "greenwich_date", "jd_ephem", "observer", "clamp", "limb",
        "ra_h", "ra_m", "ra_s", "ra_o_minus_c", "ra_limb_corr",
        "limb_side",
        "dec_sign", "dec_deg", "dec_min", "dec_sec", "dec_o_minus_c",
        "dec_limb_corr",
    ]

    current_year = None
    current_month = None
    n_good_ra = 0
    n_good_dec = 0
    clean_rows = []

    for row in raw_rows:
        jd = clean_jd(row.get("jd", ""))
        if jd is None:
            continue

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

        ra = clean_ra(row.get("ra", ""))
        ra_parsed = parse_ra(ra)
        ra_h = ra_m = ra_s = ""
        if ra_parsed:
            ra_h, ra_m, ra_s = ra_parsed
            n_good_ra += 1

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
            "limb": clean_limb_num(row.get("limb_num", "")),
            "ra_h": ra_h,
            "ra_m": ra_m,
            "ra_s": ra_s,
            "ra_o_minus_c": clean_signed(row.get("ra_oc", "")),
            "ra_limb_corr": clean_signed(row.get("ra_limb_corr", "")),
            "limb_side": clean_limb_side(row.get("limb_side", "")),
            "dec_sign": dec_sign,
            "dec_deg": dec_deg,
            "dec_min": dec_min,
            "dec_sec": dec_sec,
            "dec_o_minus_c": clean_signed(row.get("dec_oc", "")),
            "dec_limb_corr": clean_signed(row.get("dec_limb_corr", "")),
        })

    # --- Infer missing RA hours from neighboring rows ---
    # Moon moves ~0.55h RA per day. For rows with minutes+seconds but no hours,
    # interpolate from nearest complete neighbors.
    for i, r in enumerate(clean_rows):
        if r["ra_m"] and r["ra_s"] and not r["ra_h"]:
            jd_i = float(r["jd_ephem"])
            # Find nearest complete RA rows before and after
            best_h = None
            best_dist = 999
            for j in range(max(0, i - 10), min(len(clean_rows), i + 10)):
                if clean_rows[j]["ra_h"] and clean_rows[j]["ra_m"] and clean_rows[j]["ra_s"]:
                    jd_j = float(clean_rows[j]["jd_ephem"])
                    dist = abs(jd_i - jd_j)
                    if dist < best_dist:
                        best_dist = dist
                        # Estimate: Moon moves ~0.55h/day
                        ref_h = int(clean_rows[j]["ra_h"])
                        ref_m = int(clean_rows[j]["ra_m"])
                        ref_ra_h = ref_h + ref_m / 60.0
                        est_ra_h = ref_ra_h + (jd_i - jd_j) * 0.55
                        est_ra_h = est_ra_h % 24.0
                        # The missing hour should be close to est
                        m_i = int(r["ra_m"])
                        # Try each possible hour and pick closest
                        best_hour = round(est_ra_h) % 24
                        # Check if minute value makes sense
                        for candidate in [best_hour, (best_hour - 1) % 24, (best_hour + 1) % 24]:
                            cand_total = candidate + m_i / 60.0
                            if abs(((cand_total - est_ra_h + 12) % 24) - 12) < 1.5:
                                best_h = str(candidate)
                                best_dist = dist
                                break
            if best_h is not None:
                r["ra_h"] = best_h

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in clean_rows:
            writer.writerow(r)

    return clean_rows, n_good_ra, n_good_dec


def print_summary(clean_rows, n_good_ra, n_good_dec):
    n = len(clean_rows)
    if n == 0:
        print("No rows extracted!")
        return
    n_jd = sum(1 for r in clean_rows if r["jd_ephem"])
    n_date = sum(1 for r in clean_rows if r["greenwich_date"])

    print("\n" + "=" * 70)
    print("MOON TRANSIT CIRCLE OBSERVATIONS — EXTRACTION SUMMARY")
    print("=" * 70)
    print(f"  Source:  USNO Publication Vol. 23 (1982), pages {FIRST_PAGE}-{LAST_PAGE}")
    print(f"  Output:  {OUTPUT_CSV}")
    print(f"\n  Total observations:  {n:,}")
    print(f"  Complete JD:         {n_jd:,} ({100*n_jd/n:.1f}%)")
    print(f"  Complete dates:      {n_date:,} ({100*n_date/n:.1f}%)")
    print(f"  Parsed RA (h m s):   {n_good_ra:,} ({100*n_good_ra/n:.1f}%)")
    print(f"  Parsed Dec (° ′ ″): {n_good_dec:,} ({100*n_good_dec/n:.1f}%)")

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
    print("=" * 70)


def main():
    if not os.path.exists(PDF_FILE):
        print(f"Error: {PDF_FILE} not found", file=sys.stderr)
        sys.exit(1)

    print("Extracting Moon transit circle observations from USNO Pub. Vol. 23...",
          file=sys.stderr)
    print(f"  Pages {FIRST_PAGE}-{LAST_PAGE}, {DPI} DPI\n", file=sys.stderr)

    print("Step 1: Extracting page images from PDF...", file=sys.stderr)
    page_images = extract_pages(PDF_FILE)

    print("\nStep 2: OCR extraction (column-strip + row-strip)...", file=sys.stderr)
    raw_rows = extract_all_rows(page_images)

    for path in page_images.values():
        os.unlink(path)

    print(f"\nStep 3: Cleaning and writing {OUTPUT_CSV}...", file=sys.stderr)
    clean_rows, n_good_ra, n_good_dec = clean_and_write(raw_rows, OUTPUT_CSV)

    print_summary(clean_rows, n_good_ra, n_good_dec)


if __name__ == "__main__":
    main()
