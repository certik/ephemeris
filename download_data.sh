#!/usr/bin/env bash
set -euo pipefail
URL=https://github.com/certik/ephemeris/releases/download/v0.1.0
curl -fSLO --progress-bar "$URL/de441s.bsp"
curl -fSLO --progress-bar "$URL/latest_eop2.long"
curl -fSLO --progress-bar "$URL/nutation.dat"
shasum -a 256 -c <<EOF
e51e348e42dfc68d58b0aafb656836515e890cc1410931473d0cc6ebae0f140e  de441s.bsp
a0271731924132875456e6e714360f786bf0a261e73fe39c8e9d463934ad5881  latest_eop2.long
f7d8f4d701ef7d0887511c7014b0b5db2fbea1d558d0b6a6e29d57469a6c2d2f  nutation.dat
EOF
