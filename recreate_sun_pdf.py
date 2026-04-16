#!/usr/bin/env python3
"""
Recreate the Sun transit circle observation tables as a PDF matching
the original format from USNO Publication Vol. 23 (1982), pages 207-220.

Reads sun_transit_observations.csv and produces sun_transit_recreated.pdf
with the same visual layout as the scanned original for side-by-side
comparison and manual correction.
"""

import csv
import sys
from fpdf import FPDF

INPUT_CSV = "sun_transit_observations.csv"
OUTPUT_PDF = "sun_transit_recreated.pdf"

# Page layout (US Letter, landscape-like density)
PAGE_W = 215.9  # mm (8.5")
PAGE_H = 279.4  # mm (11")
MARGIN_LEFT = 12.0   # mm
MARGIN_RIGHT = 8.0   # mm
MARGIN_TOP = 10.0    # mm

FONT_SIZE = 7.5       # pt — matches original density
LINE_HEIGHT = 2.8     # mm — tight like original (~75 rows/page)
GROUP_GAP = 1.4       # mm extra between 5-row groups
HEADER_GAP = 3.0      # mm after column headers

ROWS_PER_GROUP = 5
GROUPS_PER_PAGE = 15   # 15 groups × 5 rows = 75 rows per page

# Month display names matching the original publication
MONTH_DISPLAY = {
    "Jan": "Jan.", "Feb": "Feb.", "Mar": "Mar.", "Apr": "Apr.",
    "May": "May",  "Jun": "June", "Jul": "July", "Aug": "Aug.",
    "Sep": "Sept.", "Oct": "Oct.", "Nov": "Nov.", "Dec": "Dec.",
}

# Column x-positions (mm from left margin) — tuned to match original
# These are left-edge positions for left-aligned fields,
# or right-edge positions for right-aligned fields
COL = {
    "year":     0.0,     # left-aligned
    "month":    9.0,     # left-aligned
    "day":      26.0,    # right-aligned (to this x)
    "jd":       56.0,    # right-aligned
    "obsr":     60.0,    # left-aligned
    "cl":       69.0,    # left-aligned
    "ra_h":     80.0,    # right-aligned
    "ra_m":     87.0,    # right-aligned
    "ra_s":     100.0,   # right-aligned
    "ra_oc":    116.0,   # right-aligned
    "dec_sign": 122.0,   # centered
    "dec_deg":  130.0,   # right-aligned
    "dec_min":  137.0,   # right-aligned
    "dec_sec":  150.0,   # right-aligned
    "dec_oc":   163.0,   # right-aligned
    "dist":     186.0,   # right-aligned
}


def sanitize(text):
    """Strip non-ASCII characters from text (OCR artifacts)."""
    if not text:
        return text
    return text.encode("ascii", errors="ignore").decode("ascii")


class SunTransitPDF(FPDF):
    def __init__(self):
        super().__init__(orientation="P", unit="mm", format="Letter")
        self.page_number_val = 207  # starting page number to match original
        self.set_auto_page_break(auto=False)

    def new_data_page(self):
        """Start a new page with headers."""
        self.add_page()
        x0 = MARGIN_LEFT

        # Top line: publication title + page number
        self.set_font("Courier", "B", 8)
        self.set_xy(x0, MARGIN_TOP)
        self.cell(0, 3, "OBSERVATIONS OF THE SUN, MOON, AND PLANETS", align="C")
        # Page number at far right
        self.set_xy(PAGE_W - MARGIN_RIGHT - 15, MARGIN_TOP)
        self.cell(15, 3, str(self.page_number_val), align="R")
        self.page_number_val += 1

        # "SUN" centered
        y = MARGIN_TOP + 7
        self.set_font("Courier", "B", 10)
        self.set_xy(x0, y)
        self.cell(0, 4, "SUN", align="C")

        # Column headers — line 1
        y += 8
        self.set_font("Courier", "", FONT_SIZE)
        self._text_at(x0 + COL["year"], y, "Greenwich Date")
        self._text_at(x0 + 30, y, "Julian")
        self._rtext_at(x0 + COL["ra_s"], y, "Right Ascension")
        self._rtext_at(x0 + COL["dec_sec"], y, "Declination")

        # Column headers — line 2
        y += LINE_HEIGHT + 0.5
        self._text_at(x0 + 28, y, "Ephemeris Date")
        self._text_at(x0 + COL["obsr"], y, "Obsr.")
        self._text_at(x0 + COL["cl"], y, "Cl.")
        self._rtext_at(x0 + COL["ra_s"], y, "Observed")
        self._rtext_at(x0 + COL["ra_oc"], y, "(O-C)")
        self._rtext_at(x0 + COL["dec_sec"], y, "Observed")
        self._rtext_at(x0 + COL["dec_oc"], y, "(O-C)")
        self._rtext_at(x0 + COL["dist"], y, "Distance")

        # Sub-headers — line 3 (units)
        y += LINE_HEIGHT + 0.5
        self.set_font("Courier", "", FONT_SIZE - 1)
        self._rtext_at(x0 + 38, y, "24")
        self._text_at(x0 + COL["ra_h"] - 5, y, "h")
        self._text_at(x0 + COL["ra_m"] - 3, y, "m")
        self._text_at(x0 + COL["ra_s"] - 9, y, "s")
        self._rtext_at(x0 + COL["ra_oc"], y, "''")
        self._text_at(x0 + COL["dec_sign"] + 3, y, "o")
        self._text_at(x0 + COL["dec_min"] - 3, y, "'")
        self._text_at(x0 + COL["dec_sec"] - 9, y, "''")
        self._rtext_at(x0 + COL["dec_oc"], y, "''")
        self.set_font("Courier", "", FONT_SIZE)

        y += HEADER_GAP
        return y

    def _text_at(self, x, y, text):
        """Draw left-aligned text at (x, y)."""
        self.set_xy(x, y)
        self.cell(0, LINE_HEIGHT, text)

    def _rtext_at(self, x, y, text):
        """Draw right-aligned text ending at x."""
        w = self.get_string_width(text)
        self.set_xy(x - w, y)
        self.cell(w, LINE_HEIGHT, text)

    def _ctext_at(self, x, y, text):
        """Draw centered text at x."""
        w = self.get_string_width(text)
        self.set_xy(x - w / 2, y)
        self.cell(w, LINE_HEIGHT, text)


def format_month(month_abbr):
    """Convert 3-letter month to display format matching original."""
    return MONTH_DISPLAY.get(month_abbr, month_abbr + ".")


def parse_date(date_str):
    """Parse 'YYYY Mon. DD.dd' into (year, month, day)."""
    if not date_str or len(date_str) < 5:
        return None, None, None
    parts = date_str.split()
    if len(parts) < 3:
        return None, None, None
    year = parts[0]
    month = parts[1].rstrip(".")  # Remove period for lookup
    day = parts[2]
    return year, month, day


def main():
    # Read CSV
    with open(INPUT_CSV) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    print(f"Read {len(rows)} observations from {INPUT_CSV}")

    pdf = SunTransitPDF()
    x0 = MARGIN_LEFT

    row_idx = 0
    group_count = 0
    prev_year = None
    prev_month = None

    y = pdf.new_data_page()

    for row in rows:
        # Group gap every ROWS_PER_GROUP rows
        if row_idx > 0 and row_idx % ROWS_PER_GROUP == 0:
            group_count += 1
            # New page after GROUPS_PER_PAGE groups
            if group_count >= GROUPS_PER_PAGE:
                y = pdf.new_data_page()
                group_count = 0
                row_idx = 0
                prev_year = None
            else:
                y += GROUP_GAP
                prev_year = None  # show year at start of each group

        # Parse date
        year, month, day = parse_date(row["greenwich_date"])

        # --- Greenwich Date ---
        if year and year != prev_year:
            pdf._text_at(x0 + COL["year"], y, year)
            prev_year = year

        if month:
            month_str = format_month(month)
            pdf._text_at(x0 + COL["month"], y, month_str)
            prev_month = month

        if day:
            pdf._rtext_at(x0 + COL["day"], y, day)

        # --- Julian Ephemeris Date ---
        jd = row["jd_ephem"]
        if jd:
            pdf._rtext_at(x0 + COL["jd"], y, jd)

        # --- Observer ---
        obsr = sanitize(row["observer"])
        if obsr:
            pdf._text_at(x0 + COL["obsr"], y, obsr)

        # --- Clamp ---
        cl = sanitize(row["clamp"])
        if cl:
            pdf._text_at(x0 + COL["cl"], y, cl)

        # --- Right Ascension ---
        ra_h = row["ra_h"]
        ra_m = row["ra_m"]
        ra_s = row["ra_s"]
        if ra_h:
            pdf._rtext_at(x0 + COL["ra_h"], y, ra_h)
        if ra_m:
            pdf._rtext_at(x0 + COL["ra_m"], y, f"{int(ra_m):02d}")
        if ra_s:
            # Pad seconds to XX.XXX format
            try:
                s_val = float(ra_s)
                s_str = f"{s_val:06.3f}"
            except ValueError:
                s_str = ra_s
            pdf._rtext_at(x0 + COL["ra_s"], y, s_str)

        # --- RA (O-C) ---
        ra_oc = row["ra_o_minus_c"]
        if ra_oc:
            pdf._rtext_at(x0 + COL["ra_oc"], y, ra_oc)

        # --- Declination ---
        dec_sign = row["dec_sign"]
        dec_deg = row["dec_deg"]
        dec_min = row["dec_min"]
        dec_sec = row["dec_sec"]

        if dec_sign:
            sign_char = "-" if dec_sign == "-" else "+"
            pdf._ctext_at(x0 + COL["dec_sign"], y, sign_char)

        if dec_deg:
            pdf._rtext_at(x0 + COL["dec_deg"], y, dec_deg)
        if dec_min:
            pdf._rtext_at(x0 + COL["dec_min"], y, f"{int(dec_min):02d}")
        if dec_sec:
            try:
                sec_val = float(dec_sec)
                sec_str = f"{sec_val:05.2f}"
            except ValueError:
                sec_str = dec_sec
            pdf._rtext_at(x0 + COL["dec_sec"], y, sec_str)

        # --- Dec (O-C) ---
        dec_oc = row["dec_o_minus_c"]
        if dec_oc:
            pdf._rtext_at(x0 + COL["dec_oc"], y, dec_oc)

        # --- Distance ---
        dist = row["distance_au"]
        if dist:
            pdf._rtext_at(x0 + COL["dist"], y, dist)

        y += LINE_HEIGHT
        row_idx += 1

    # Track final group
    group_count += 1

    pdf.output(OUTPUT_PDF)
    n_pages = pdf.page_number_val - 207
    print(f"Wrote {OUTPUT_PDF} ({n_pages} pages, {len(rows)} observations)")


if __name__ == "__main__":
    main()
