"""Behaviour of a discrete design with many BSs scattered over a large area (design only; no AO, no search).

    python3 tools/large_network.py MODULE:FUNCTION [--B 9] [--side 300] [--K 40] [--N 4] [--Smax 10] [--seed 1]

BSs on a jittered square grid over a side x side area (height 15 m), K UEs uniform (height 1.5 m), the channel model
of the thesis (hrs.make_network with explicit BS/UE positions). Prints the private-layer choice, how many BSs carry each
core user's private stream, and per BS the streams it carries and its distance to the UEs it serves.
"""
import os, sys, argparse, importlib
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))
import hrs as H_, run as R, discrete as DC


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('design'); ap.add_argument('--B', type=int, default=9); ap.add_argument('--side', type=float, default=300.0)
    ap.add_argument('--K', type=int, default=40); ap.add_argument('--N', type=int, default=4); ap.add_argument('--Smax', type=int, default=10)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    g = int(np.ceil(np.sqrt(a.B))); step = a.side / g
    cells = [(i, j) for i in range(g) for j in range(g)][:a.B]
    BS = [((i + 0.5 + 0.3 * rng.uniform(-1, 1)) * step, (j + 0.5 + 0.3 * rng.uniform(-1, 1)) * step, 15.0) for i, j in cells]
    UE = [(x, y, 1.5) for x, y in rng.uniform(0, a.side, (a.K, 2))]
    p = H_.params(K=a.K, L=3, N=a.N, M=100, Smax=a.Smax, Pmax=1.0, Dmax=10.0, Rth=1.0, bank=R.BANK, B=a.B)
    bank = dict(R.BANK, bsSetups={'even': dict(R.BANK['bsSetups']['even'], UE={'K%d' % a.K: {'random': UE}})})
    net = H_.make_network(p, bank, drop=1, kind='random', BS=BS, UE=UE)
    mod, fn = a.design.split(':')
    d = getattr(importlib.import_module(mod), fn)(net, p)
    DC.check(d, p)
    bm = np.ones((a.B, d.nS), bool) if d.bsmask is None else np.asarray(d.bsmask, bool)
    dist = np.linalg.norm(np.asarray(net['UE'])[:, None, :2] - np.asarray(net['BS'])[None, :, :2], axis=2)
    print(f'B = {a.B} BSs over {a.side:.0f} m x {a.side:.0f} m, K = {a.K}, N = {a.N}, S_max = {a.Smax}: '
          f'{d.nS} distinct streams ({len(d.private_served)} private, {len(d.served_groups)} group-common), rule {d.rule}')
    lay = getattr(d, 'layer', None) or (dict(n_private=d.dof[0], s=d.dof[1]) if hasattr(d, 'dof') else None)
    if lay: print(f"private layer: {lay['n_private']} private streams, each carried by s = {lay['s']} of the {a.B} BSs")
    ncar = [int(bm[:, d.col[1 + d.G + k]].sum()) for k in d.core]
    rank = [int(sorted(range(a.B), key=lambda b: dist[k, b]).index(b)) for k in d.core
            for b in np.flatnonzero(bm[:, d.col[1 + d.G + k]])]
    print(f'core users: {len(d.core)}; BSs per core private stream: {np.bincount(ncar).tolist()} (count by number of BSs); '
          f'distance rank of a carrying BS among the BSs of its UE: {np.bincount(rank, minlength=a.B).tolist()}')
    for b in range(a.B):
        js = [d.served[i] for i in np.flatnonzero(bm[b])]
        ues = sorted({k for j in js if j > 0 for k in d.decoders[j]})
        print(f'  BS {b} at ({BS[b][0]:5.0f},{BS[b][1]:5.0f}): {len(js)} streams '
              f'({sum(1 <= j <= d.G for j in js)} group-common, {sum(j > d.G for j in js)} private), '
              f'serves {len(ues)} UEs at median distance {np.median(dist[ues, b]) if ues else float("nan"):.0f} m '
              f'(all UEs: {np.median(dist[:, b]):.0f} m)')


if __name__ == '__main__':
    main()
