"""R1 on all Q = 50 grid runs common to a variant and its 1-layer counterpart only (no intersection with JH/J1; the
same paired statistics as tools/report.py), per scenario and R_th, optionally per setting.

    python3 tools/r1_pair.py X [--by-setting] [--drops 901,902,...]
"""
import os, sys, argparse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'harness')); sys.path.insert(0, os.path.join(ROOT, 'tools'))
import variants as V, report as RP


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('variant'); ap.add_argument('--by-setting', action='store_true')
    ap.add_argument('--drops', default=None); ap.add_argument('--budget', type=int, default=50)
    a = ap.parse_args()
    X = a.variant; X1 = V.VARIANTS[X]['one_layer']
    dx, d1 = RP.load(X, a.budget), RP.load(X1, a.budget)
    drops = None if a.drops is None else {int(x) for x in a.drops.split(',')}
    for sc in RP.SCEN:
        if sc not in dx or sc not in d1: continue
        keys = [k for k in set(dx[sc]) & set(d1[sc]) if drops is None or k[1] in drops]
        groups = [('all settings', keys)] + ([(f'setting {st}', [k for k in keys if k[0] == st]) for st in sorted({k[0] for k in keys})]
                                             if a.by_setting else [])
        print(f'\n===== scenario {sc}, Q = {a.budget}: {len(keys)} runs common to {X} and {X1}, drops {sorted({k[1] for k in keys})}')
        for title, ks in groups:
            print(f'-- {sc}, {title}')
            print(f"{'R_th':>6} {'n':>4} {X:>9} {X1:>9} {X + '-' + X1:>17}   QoS met {X}/{X1}")
            for R in sorted({k[2] for k in ks}) + ['all']:
                kk = [k for k in ks if R == 'all' or k[2] == R]
                m, e = RP.stats([dx[sc][k]['rate'] - d1[sc][k]['rate'] for k in kk])
                print(f"{(f'{R:g}' if R != 'all' else 'all'):>6} {len(kk):>4} {np.mean([dx[sc][k]['rate'] for k in kk]):9.2f} "
                      f"{np.mean([d1[sc][k]['rate'] for k in kk]):9.2f} {m:+8.2f}±{e:5.2f}   "
                      f"{sum(dx[sc][k]['qos_met'] for k in kk)}/{sum(d1[sc][k]['qos_met'] for k in kk)}")


if __name__ == '__main__':
    main()
