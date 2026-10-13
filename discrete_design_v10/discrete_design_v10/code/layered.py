"""Discrete design with a stream set per BS (design round v10, BRIEF.md): the layered design.

    design(net, p)          the design (harness variant SP)
    design_1layer(net, p)   its 1-layer RSMA counterpart: the same private layer, no group layer (variant SP1)

1. Order. The users in the order of semi-orthogonal user selection (SUS) on the stacked direct channels, as in the
   thesis (hrs.sus): c_1, ..., c_|C|, |C| = min(K, BN).
2. Private layer. The private streams of the first n users of the order, the stream of c_i carried by the s BSs with
   the largest direct-channel energy to c_i among the BSs with a free fronthaul slot (each BS keeps one of its S_max
   slots for the super-common stream). (n, s) maximize the zero-forcing sum rate of this layer,
       R_ZF(n, s) = sum_i log2(1 + (P/n) / [(H_i H_i^H)^-1]_ii),   P = B P_max,
   H_i the direct channels of c_1..c_n on the antennas of the BSs that carry c_i's stream (R_ZF = -inf if those BSs
   have fewer than n antennas: the stream cannot be zero-forced), over n = 1..|C| and s = 1..B (ties: the larger n,
   then the larger s).
3. Groups. Every UE without a private stream joins the private-layer user with the most correlated direct channel;
   every private-layer user with at least one follower heads a group (disjoint, sizes 2..K-1, at most K/2 groups).
4. Group layer. Every BS fills the slots left by the private layer with group-common streams, in decreasing order
   of its direct-channel energy to the weakest member of the group.
5. Every BS carries the super-common stream.
Deterministic, from the direct channels and the settings (B, N, K, S_max, P_max) only.
"""
import numpy as np
import hrs as H_


def bs_energy(Hd, B, N):
    """E[k, b] = ||h_{k,b}||^2: direct-channel energy from BS b to UE k (noise-normalized channels)."""
    return np.stack([np.sum(np.abs(Hd[:, b * N:(b + 1) * N]) ** 2, axis=1) for b in range(B)], axis=1)


def zf_rate(Hd, users, carriers, N, P):
    """ZF sum rate of the private streams of `users`, the stream of user k precoded on the antennas of the BSs
    carriers[k], with the power P shared equally; -inf if a stream cannot be zero-forced."""
    n = len(users); r = 0.0
    for i, k in enumerate(users):
        ant = np.concatenate([np.arange(b * N, (b + 1) * N) for b in carriers[k]])
        if len(ant) < n:
            return -np.inf
        Hk = Hd[users][:, ant]
        try:
            g = 1.0 / max(np.real(np.linalg.inv(Hk @ Hk.conj().T)[i, i]), 1e-300)
        except np.linalg.LinAlgError:
            return -np.inf
        r += np.log2(1 + P / n * g)
    return r


def private_layer(Hd, p, order, E):
    """Steps 1-2: (users of the private layer, their carrying BSs, s)."""
    B, N, S = p['B'], p['N'], p['Smax']
    P = B * p['Pmax']
    best = None
    for s in range(1, B + 1):
        load = np.zeros(B, int); car = {}
        for n, k in enumerate(order, 1):
            bs = [b for b in np.argsort(-E[k], kind='stable') if load[b] < S - 1][:s]
            if len(bs) < s:
                break
            car[k] = bs; load[bs] += 1
            r = zf_rate(Hd, order[:n], car, N, P)
            if best is None or r > best[0] + 1e-9 or (abs(r - best[0]) <= 1e-9 and (n, s) > (len(best[1]), best[3])):
                best = (r, order[:n], dict(car), s)
    return best[1], best[2], best[3]


def design(net, p, one_layer=False):
    K, B, N, S = p['K'], p['B'], p['N'], p['Smax']
    Hd = net['Hd']
    E = bs_energy(Hd, B, N)
    order = H_.sus(Hd, p['BN'])                                   # step 1
    core, car, s = private_layer(Hd, p, order, E)                 # step 2
    groups = []
    if not one_layer:                                             # step 3
        energy = np.sum(np.abs(Hd) ** 2, axis=1)
        Hn = Hd / np.sqrt(energy)[:, None]
        corr = np.abs(Hn @ Hn.conj().T) ** 2
        follow = {c: [] for c in core}
        for k in range(K):
            if k not in core:
                follow[core[int(np.argmax(corr[k, core]))]].append(k)
        groups = [[c] + follow[c] for c in core if follow[c]]
    G = len(groups)
    sets = [[0] for _ in range(B)]                                # step 5
    for k in core:
        for b in car[k]:
            sets[b].append(1 + G + k)
    weakest = {1 + g: E[groups[g]].min(0) for g in range(G)}      # step 4
    for b in range(B):
        sets[b] += sorted(weakest, key=lambda j: -weakest[j][b])[:S - len(sets[b])]
    served = sorted(set().union(*sets))
    d = H_.Design(groups, served, K)
    d.core = list(core)
    d.rule = 'layered' + ('_1layer' if one_layer else '')
    others = [k for k in range(K) if k not in core]
    d.priority = [0] + [1 + G + k for k in core] + sorted(weakest, key=lambda j: (-len(groups[j - 1]), j)) \
        + [1 + G + k for k in sorted(others, key=lambda k: (-E[k].sum(), k))]
    d.layer = dict(n_private=len(core), s=s)
    if any(set(x) != set(served) for x in sets):
        d.bsmask = np.array([[j in sets[b] for j in d.served] for b in range(B)])
        d.bs_sets = [sorted(x) for x in sets]
    return d


def design_1layer(net, p):
    return design(net, p, one_layer=True)
