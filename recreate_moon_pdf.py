#!/usr/bin/env python3
"""
Recreate the Moon transit circle observation tables as a PDF matching
the original format from USNO Publication Vol. 23 (1982), pages 221-230.

Reads moon_transit_observations.csv and produces moon_transit_recreated.pdf
for side-by-side comparison with the original scanned pages.
"""

import csv
import sys
from fpdf import FPDF

INPUT_CSV = "moon_transit_observations.csv"
OUTPUT_PDF = "moon_transit_recreated.pdf"

PAGE_W = 215.9
PAGE_H = 279.4
MARGIN_LEFT = 8.0
MARGIN_RIGHT = 6.0
MARGIN_TOP = 10.0

FONT_SIZE = 6.5
LINE_HEIGHT = 2.6
GROUP_GAP = 1.3
HEADER_GAP = 3.0

ROWS_PER_GROUP = 5
GROUPS_PER_PAGE = 15

MONTH_DISPLAY = {
    "Jan": "Jan.", "Feb": "Feb.", "Mar": "Mar.", "Apr": "Apr.",
    "May": "May",  "Jun": "June", "Jul": "July", "Aug": "Aug.",
    "Sep": "Sept.", "Oct": "Oct.", "Nov": "Nov.", "Dec": "Dec.",
}

# Column positions (mm from left margin) — Moon table is wider than Sun
COL = {
    "year":     0.0,
    "month":    9.0,
    "day":      24.0,
    "jd":       51.0,
    "obsr":     54.0,
    "cl":       61.0,
    "limb":     66.0,
    "ra_h":     74.0,
    "ra_m":     80.0,
    "ra_s":     91.0,
    "ra_oc":    103.0,
    "ra_limb":  115.0,
    "limb_side": 121.0,
    "dec_sign": 125.0,
    "dec_deg":  131.0,
    "dec_min":  137.0,
    "dec_sec":  148.0,
    "dec_oc":   159.0,
    "dec_limb": 172.0,
}


def sanitize(text):
    if not text:
        return text
    return text.encode("ascii", errors="ignore").decode("ascii")


class MoonTransitPDF(FPDF):
    def __init__(self):
        super().__init__(orientation="P", unit="mm", format="Letter")
        self.page_number_val = 221
        self.set_auto_page_break(auto=False)

    def new_data_page(self):
        self.add_page()
        x0 = MARGIN_LEFT

        self.set_font("Courier", "B", 7)
        self.set_xy(x0, MARGIN_TOP)
        self.cell(0, 3, "OBSERVATIONS OF THE SUN, MOON, AND PLANETS", align="C")
        self.set_xy(PAGE_W - MARGIN_RIGHT - 15, MARGIN_TOP)
        self.cell(15, 3, str(self.page_number_val), align="R")
        self.page_number_val += 1

        y = MARGIN_TOP + 7
        self.set_font("Courier", "B", 9)
        self.set_xy(x0, y)
        self.cell(0, 4, "MOON", align="C")

        y += 8
        self.set_font("Courier", "", FONT_SIZE)
        self._text_at(x0 + COL["year"], y, "Greenwich Date")
        self._text_at(x0 + 28, y, "Julian")
        self._rtext_at(x0 + COL["ra_s"], y, "Right Ascension")
        self._rtext_at(x0 + COL["dec_sec"], y, "Declination")

        y += LINE_HEIGHT + 0.5
        self._text_at(x0 + 26, y, "Ephemeris Date")
        self._text_at(x0 + COL["obsr"], y, "Obsr.")
        self._text_at(x0 + COL["cl"], y, "Cl.")
        self._text_at(x0 + COL["limb"], y, "Limb")
        self._rtext_at(x0 + COL["ra_s"], y, "Observed")
        self._rtext_at(x0 + COL["ra_oc"], y, "(O-C)")
        self._text_at(x0 + COL["ra_limb"] - 5, y, "Limb Corr.")
        self._rtext_at(x0 + COL["dec_sec"], y, "Observed")
        self._rtext_at(x0 + COL["dec_oc"], y, "(O-C)")
        self._text_at(x0 + COL["dec_limb"] - 3, y, "Limb Corr.")

        y += LINE_HEIGHT + 0.5
        self.set_font("Courier", "", FONT_SIZE - 1)
        self._rtext_at(x0 + 36, y, "24")
        self._text_at(x0 + COL["ra_h"] - 4, y, "h")
        self._text_at(x0 + COL["ra_m"] - 2, y, "m")
        self._text_at(x0 + COL["ra_s"] - 8, y, "s")
        self._rtext_at(x0 + COL["ra_oc"], y, "s")
        self._text_at(x0 + COL["dec_sign"] + 2, y, "o")
        self._text_at(x0 + COL["dec_min"] - 2, y, "'")
        self._text_at(x0 + COL["dec_sec"] - 8, y, "\"")
        self._rtext_at(x0 + COL["dec_oc"], y, "\"")
        self._rtext_at(x0 + COL["dec_limb"] + 8, y, "\"")
        self.set_font("Courier", "", FONT_SIZE)

        y += HEADER_GAP
        return y

    def _text_at(self, x, y, text):
        self.set_xy(x, y)
        self.cell(0, LINE_HEIGHT, text)

    def _rtext_at(self, x, y, text):
        w = self.get_string_width(text)
        self.set_xy(x - w, y)
        self.cell(w, LINE_HEIGHT, text)

    def _ctext_at(self, x, y, text):
        w = self.get_string_width(text)
        self.set_xy(x - w / 2, y)
        self.cell(w, LINE_HEIGHT, text)


def format_month(month_abbr):
    return MONTH_DISPLAY.get(month_abbr, month_abbr + ".")


def parse_date(date_str):
    if not date_str or len(date_str) < 5:
        return None, None, None
    parts = date_str.split()
    if len(parts) < 3:
        return None, None, None
    return parts[0], parts[1].rstrip("."), parts[2]


def main():
    with open(INPUT_CSV) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    print(f"Read {len(rows)} observations from {INPUT_CSV}")

    pdf = MoonTransitPDF()
    x0 = MARGIN_LEFT

    row_idx = 0
    group_count = 0
    prev_year = None

    y = pdf.new_data_page()

    for row in rows:
        if row_idx > 0 and row_idx % ROWS_PER_GROUP == 0:
            group_count += 1
            if group_count >= GROUPS_PER_PAGE:
                y = pdf.new_data_page()
                group_count = 0
                row_idx = 0
                prev_year = None
            else:
                y += GROUP_GAP
                prev_year = None

        year, month, day = parse_date(row["greenwich_date"])

        if year and year != prev_year:
            pdf._text_at(x0 + COL["year"], y, year)
            prev_year = year

        if month:
            pdf._text_at(x0 + COL["month"], y, format_month(month))

        if day:
            pdf._rtext_at(x0 + COL["day"], y, day)

        jd = row["jd_ephem"]
        if jd:
            pdf._rtext_at(x0 + COL["jd"], y, jd)

        obsr = sanitize(row["observer"])
        if obsr:
            pdf._text_at(x0 + COL["obsr"], y, obsr)

        cl = sanitize(row["clamp"])
        if cl:
            pdf._text_at(x0 + COL["cl"], y, cl)

        limb = row["limb"]
        if limb:
            pdf._ctext_at(x0 + COL["limb"] + 2, y, limb)

        # RA
        ra_h, ra_m, ra_s = row["ra_h"], row["ra_m"], row["ra_s"]
        if ra_h:
            pdf._rtext_at(x0 + COL["ra_h"], y, ra_h)
        if ra_m:
            pdf._rtext_at(x0 + COL["ra_m"], y, f"{int(ra_m):02d}")
        if ra_s:
            try:
                pdf._rtext_at(x0 + COL["ra_s"], y, f"{float(ra_s):06.3f}")
            except ValueError:
                pdf._rtext_at(x0 + COL["ra_s"], y, ra_s)

        ra_oc = row["ra_o_minus_c"]
        if ra_oc:
            pdf._rtext_at(x0 + COL["ra_oc"], y, ra_oc)

        ra_limb = row["ra_limb_corr"]
        if ra_limb:
            pdf._rtext_at(x0 + COL["ra_limb"], y, ra_limb)

        limb_side = row["limb_side"]
        if limb_side:
            pdf._text_at(x0 + COL["limb_side"], y, limb_side)

        # Dec
        dec_sign = row["dec_sign"]
        dec_deg = row["dec_deg"]
        dec_min = row["dec_min"]
        dec_sec = row["dec_sec"]

        if dec_sign:
            pdf._ctext_at(x0 + COL["dec_sign"], y, "-" if dec_sign == "-" else "+")
        if dec_deg:
            pdf._rtext_at(x0 + COL["dec_deg"], y, dec_deg)
        if dec_min:
            pdf._rtext_at(x0 + COL["dec_min"], y, f"{int(dec_min):02d}")
        if dec_sec:
            try:
                pdf._rtext_at(x0 + COL["dec_sec"], y, f"{float(dec_sec):05.2f}")
            except ValueError:
                pdf._rtext_at(x0 + COL["dec_sec"], y, dec_sec)

        dec_oc = row["dec_o_minus_c"]
        if dec_oc:
            pdf._rtext_at(x0 + COL["dec_oc"], y, dec_oc)

        dec_limb = row["dec_limb_corr"]
        if dec_limb:
            pdf._rtext_at(x0 + COL["dec_limb"] + 8, y, dec_limb)

        y += LINE_HEIGHT
        row_idx += 1

    pdf.output(OUTPUT_PDF)
    n_pages = pdf.page_number_val - 221
    print(f"Wrote {OUTPUT_PDF} ({n_pages} pages, {len(rows)} observations)")


if __name__ == "__main__":
    main()
