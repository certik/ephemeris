"""
Collect every Moon-declination measurement from the tycho*.py suite
into a single file for orbit fitting.

For each observation we compute:
   - UTC timestamp (Tycho's clock corrected for all documented drifts
     and for the equation of time, as already encoded in each script's
     t_at() helper);
   - Moon CENTER apparent topocentric declination in degrees
     (the raw limb/horn/cornu reading is corrected by +/- sd, where sd
     is the model semidiameter multiplied by the cusp-tip factor
     sin(PA_sun)).

Output:  moon_dec_observations.dat  (whitespace-separated, one obs/line).
Columns:
   iso_utc           source     horn     dec_raw_deg  dec_center_deg  jd_tt
"""
import contextlib
import importlib
import io
import sys
import warnings

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, get_body, AltAz, PrecessedGeocentric
from astropy.time import Time

warnings.filterwarnings('ignore')

MODULES = [f'tycho{i}' for i in range(3, 18)]  # tycho3..tycho17


def moon_geo_sd_cusp(mod, t: Time):
    """Compute apparent Moon center (dec, sd_deg, cusp_factor) at time t
    using the module's loc and helpers."""
    mo = get_body('moon', t, location=mod.loc)
    mo_d = mod.apparent_of_date(mo, t)
    sd = mod.moon_semid_deg(mo)
    # cusp factor requires Sun position
    sun_d = mod.apparent_of_date(mod.sun_geocentric_icrs(t), t)
    if hasattr(mod, 'cusp_factor'):
        fac = mod.cusp_factor(mo_d, sun_d)
    else:
        dra = ((sun_d.ra.deg - mo_d.ra.deg + 180.0) % 360.0) - 180.0
        east = dra * np.cos(np.radians(mo_d.dec.deg))
        north = sun_d.dec.deg - mo_d.dec.deg
        fac = abs(np.sin(np.arctan2(east, north)))
    return mo_d.dec.deg, sd, fac


def main():
    obs = []
    for name in MODULES:
        print(f'loading {name} ...', file=sys.stderr)
        with contextlib.redirect_stdout(io.StringIO()):
            mod = importlib.import_module(name)
        if not hasattr(mod, 'MOON_DEC'):
            print(f'  no MOON_DEC in {name}, skipping', file=sys.stderr)
            continue
        armilla_bias_deg = getattr(
            mod, 'ARMILLA_DEC_BIAS_ARCMIN', 0.0) / 60.0
        for row in mod.MOON_DEC:
            h, m, horn, v = row
            t = mod.t_at(h, m)
            dec_app, sd, fac = moon_geo_sd_cusp(mod, t)
            sd_eff = sd * fac  # cornu cusp tips
            # Reduce raw reading to Moon-center declination.
            # 'upper' / 'lower' are horn (cornu) OR limb -- we treat both as
            # the tangent extremum of the illuminated disk, which for
            # crescent/gibbous phases is the cusp line.  For near-full
            # phases cusp_factor ~ 1 so sd_eff ~ sd anyway.
            if horn == 'upper':
                v_center = v - sd_eff
            elif horn == 'lower':
                v_center = v + sd_eff
            else:  # 'center' (derived by Tycho himself)
                v_center = v
            v_center -= armilla_bias_deg
            obs.append({
                'utc': t.utc.iso,
                'jd_tt': float(t.tt.jd),
                'source': name,
                'horn': horn,
                'dec_raw_deg': v,
                'dec_center_deg': v_center,
            })
    obs.sort(key=lambda o: o['jd_tt'])
    out_path = 'moon_dec_observations.dat'
    with open(out_path, 'w') as f:
        f.write('# Moon CENTER apparent topocentric declination measurements\n')
        f.write('# from Tycho Brahe, Uraniborg 1586 (tycho*.py suite).\n')
        f.write(f'# {len(obs)} observations.\n')
        f.write('# Columns: iso_utc  source  horn    dec_raw_deg  '
                'dec_center_deg  jd_tt\n')
        for o in obs:
            f.write(f'{o["utc"]:<26s}  {o["source"]:<8s}  {o["horn"]:<6s}  '
                    f'{o["dec_raw_deg"]:+11.6f}  '
                    f'{o["dec_center_deg"]:+11.6f}  '
                    f'{o["jd_tt"]:.8f}\n')
    print(f'wrote {len(obs)} observations to {out_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
