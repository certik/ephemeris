"""
Collect Moon-declination observations from the tycho*.py suite using ONLY
observed quantities to derive the Moon's center declination.

Method
------
For every pair of (upper cornu/limb, lower cornu/limb) observations taken at
the SAME moment (or within a few minutes, during which the Moon's declination
barely moves), the line joining the two tangent points passes through the
Moon's geometric center (the two tips are diametrically opposite points on
the Moon's limb).  Therefore

    dec_center = (dec_upper + dec_lower) / 2

is a purely geometric identity -- it does NOT depend on any knowledge of the
Moon's distance, phase, or cusp tilt.

Lone 'upper' or 'lower' readings without a pairable partner are DROPPED.
'center' rows (where Tycho himself derived the center via some other means)
are passed through unchanged.

The only ephemeris-derived quantity we retain is the equation-of-time shift
used inside each script's ``t_at()`` to convert Tycho's local apparent solar
time to UT.  That correction uses only the SUN's position (measured by Tycho
himself in the notebook via transit observations), not the Moon's, so it is
independent of the Moon orbit we are about to fit.

Output:  moon_dec_observations.dat
"""
import contextlib
import importlib
import io
import sys
import warnings

warnings.filterwarnings('ignore')

MODULES = [f'tycho{i}' for i in range(3, 18)]
PAIR_TOLERANCE_MIN = 30.0  # generous: linear dec-rate errors cancel in the mean


def pair_rows(rows):
    """Return (list_of_center_rows, list_of_dropped_rows).

    Each element of list_of_center_rows is (h, m, dec_center, tag) where
    ``tag`` is 'pair' when derived from an (upper, lower) pair (averaged
    time-stamp) or 'center' for rows Tycho already labelled 'center'.
    """
    used = [False] * len(rows)
    paired = []
    for i, (h1, m1, horn1, v1) in enumerate(rows):
        if used[i]:
            continue
        if horn1 == 'center':
            paired.append((h1, m1, v1, 'center'))
            used[i] = True
            continue
        if horn1 not in ('upper', 'lower'):
            continue
        want = 'lower' if horn1 == 'upper' else 'upper'
        t1_min = h1*60 + m1
        # Find closest-in-time row with the complementary horn.
        best_j = -1
        best_dt = PAIR_TOLERANCE_MIN + 1e-9
        for j, (h2, m2, horn2, _) in enumerate(rows):
            if used[j] or j == i or horn2 != want:
                continue
            dt = abs((h2*60 + m2) - t1_min)
            if dt < best_dt:
                best_dt = dt
                best_j = j
        if best_j < 0:
            continue
        h2, m2, _, v2 = rows[best_j]
        used[i] = True
        used[best_j] = True
        h_mid_min = (t1_min + h2*60 + m2) / 2.0
        hh = int(h_mid_min // 60)
        mm = h_mid_min - hh*60
        center = (v1 + v2) / 2.0
        paired.append((hh, mm, center, 'pair'))
    dropped = [rows[k] for k, u in enumerate(used) if not u]
    return paired, dropped


def main():
    all_obs = []
    totals = {'pairs': 0, 'center': 0, 'dropped': 0}
    for name in MODULES:
        print(f'loading {name} ...', file=sys.stderr)
        with contextlib.redirect_stdout(io.StringIO()):
            mod = importlib.import_module(name)
        if not hasattr(mod, 'MOON_DEC'):
            continue
        armilla_bias_deg = getattr(
            mod, 'ARMILLA_DEC_BIAS_ARCMIN', 0.0) / 60.0
        paired, dropped = pair_rows(mod.MOON_DEC)
        for h, m, dec_center, tag in paired:
            t = mod.t_at(h, m)
            v_corrected = dec_center - armilla_bias_deg
            all_obs.append({
                'utc': t.utc.iso,
                'jd_tt': float(t.tt.jd),
                'source': name,
                'tag': tag,
                'dec_center_deg': v_corrected,
            })
            totals['pairs' if tag == 'pair' else 'center'] += 1
        totals['dropped'] += len(dropped)
        if dropped:
            print(f'  [{name}] dropped {len(dropped)} unpaired row(s): '
                  f'{dropped}', file=sys.stderr)
    all_obs.sort(key=lambda o: o['jd_tt'])

    out = 'moon_dec_observations.dat'
    with open(out, 'w') as f:
        f.write('# Moon CENTER apparent topocentric declination observations\n')
        f.write('# from Tycho Brahe, Uraniborg 1586 (tycho*.py suite).\n')
        f.write('# Derivation: dec_center = (upper + lower)/2 from paired\n')
        f.write('# horn/limb observations at the same time; NO ephemeris\n')
        f.write('# knowledge of the Moon is used to compute center dec.\n')
        f.write(f'# {len(all_obs)} observations '
                f'({totals["pairs"]} pairs, {totals["center"]} direct centers, '
                f'{totals["dropped"]} unpaired rows dropped).\n')
        f.write('# Columns: iso_utc  source  tag    dec_center_deg  jd_tt\n')
        for o in all_obs:
            f.write(f'{o["utc"]:<26s}  {o["source"]:<8s}  {o["tag"]:<6s}  '
                    f'{o["dec_center_deg"]:+11.6f}  {o["jd_tt"]:.8f}\n')
    print(f'wrote {len(all_obs)} observations to {out}', file=sys.stderr)
    print(f'  {totals["pairs"]} from (upper, lower) pairs', file=sys.stderr)
    print(f'  {totals["center"]} Tycho-derived centers', file=sys.stderr)
    print(f'  {totals["dropped"]} unpaired rows dropped', file=sys.stderr)


if __name__ == '__main__':
    main()
