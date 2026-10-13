"""Compact Q = 5 screen summary of every variant tried (DESIGN.md, "Other methods tried"):
per scenario and R_th, the mean rate of X, X - X1 (its 1-layer counterpart) and X - J (paired, mean +- SE),
on the runs common to X, X1, JH and J1.

    python3 tools/summary_tried.py X1,X2,... [--budget 5]
"""
import os, sys, json, argparse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'harness')); sys.path.insert(0, os.path.join(ROOT, 'tools'))
import variants as V, report as RP


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('variants'); ap.add_argument('--budget', type=int, default=5)
    a = ap.parse_args()
    for X in a.variants.split(','):
        X1 = V.VARIANTS[X].get('one_layer')
        d = {n: RP.load(n, a.budget) for n in (X, X1, 'JH', 'J1')}
        for sc in RP.SCEN:
            if any(sc not in d[n] for n in d): continue
            keys = set(d[X][sc]) & set(d[X1][sc]) & set(d['JH'][sc]) & set(d['J1'][sc])
            if not keys: continue
            row = []
            for R in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5):
                ks = [k for k in keys if k[2] == R]
                x = [d[X][sc][k]['rate'] for k in ks]
                m1 = RP.stats([d[X][sc][k]['rate'] - d[X1][sc][k]['rate'] for k in ks])
                mj = RP.stats([d[X][sc][k]['rate'] - RP.derive_j(d['JH'][sc][k], d['J1'][sc][k])['rate'] for k in ks])
                row.append(f'{np.mean(x):6.1f} {m1[0]:+5.2f}±{m1[1]:4.2f} {mj[0]:+6.2f}')
            print(f'{X:>4} {sc} n={len(keys):3d} | ' + ' | '.join(row))


if __name__ == '__main__':
    main()
