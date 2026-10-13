"""The other discrete designs tried in design round v10 (BRIEF.md; DESIGN.md, "Other methods tried"), kept so that
every run in results/ can be reproduced. The final design is in code/layered.py.

  tiered        TA   every BS fills its slots tier by tier (core privates, group-common, other privates), within a
                     tier by its direct-channel energy to the stream's weakest decoder
  cluster       TB   every BS serves clusters (core private followed by the core user's group-common stream)
  coreonly      TC   TA without private streams of non-core UEs
  dof           DF   private layer by the spatial DoF under the fronthaul (dof_split), groups on the free slots
                DG   (dof_small) DF with DoF ties broken towards fewer BSs per stream
                DFF  (dof_fill) DF, then private streams of the other UEs in the slots still free
  dof_zf        ZG   DF on the users of greedy ZF user selection (hrs.zf_greedy)
  zfs           ZS   greedy ZF selection with s BSs per stream, s by the ZF sum rate
  sus_zf        SZ   the private layer of the final design with the thesis' groups (sus_zf_regroup = the final
                     design, verified identical to layered.design on every configuration of the harness)
Each function takes one_layer=True for its 1-layer counterpart (no group-common streams).
"""
import numpy as np
import hrs as H_


def bs_energy(net, p):
    """E[k, b] = ||h_{k,b}||^2, the direct-channel energy from BS b to UE k (noise-normalized)."""
    N, B = p['N'], p['B']
    Hd = net['Hd']
    return np.stack([np.sum(np.abs(Hd[:, b * N:(b + 1) * N]) ** 2, axis=1) for b in range(B)], axis=1)


def _finish(groups, tiers, aff, p, core, pr, rule):
    """Every BS carries the super-common stream and fills its remaining slots tier by tier, within a tier in
    decreasing order of its affinity to the stream; the served streams are those carried by some BS."""
    K, B, S = p['K'], p['B'], p['Smax']
    sets = []
    for b in range(B):
        lst = []
        for t in tiers:
            lst += sorted(t, key=lambda j: -aff(j, b))      # stable: ties keep the tier order
        sets.append(lst[:S])
    served = sorted(set().union(*sets))
    d = H_.Design(groups, served, K)
    d.core = [k for k in core if (1 + d.G + k) in d.col]
    d.priority = pr; d.rule = rule
    if any(set(s) != set(served) for s in sets):
        d.bsmask = np.array([[j in sets[b] for j in d.served] for b in range(B)])
        d.bs_sets = [sorted(s) for s in sets]
    return d


def tiered(net, p, one_layer=False):
    K = p['K']
    order, groups, energy, _ = H_.joint_structure(net['Hd'], p, 'sus')
    if one_layer: groups = []
    G = len(groups)
    pr = H_.priority_list(order, groups, energy, K, 'private_first')
    E = bs_energy(net, p)
    dec = H_.Design(groups, pr, K).decoders
    aff = lambda j, b: np.min(E[dec[j], b])
    core_priv = [1 + G + k for k in order]
    grp = [j for j in pr if 1 <= j <= G]
    other = [j for j in pr if j > G and j not in core_priv]
    return _finish(groups, [[0], core_priv, grp, other], aff, p, order, pr, 'tiered')


def tiered_1layer(net, p):
    return tiered(net, p, one_layer=True)


def _parts(net, p, one_layer=False):
    K = p['K']
    order, groups, energy, _ = H_.joint_structure(net['Hd'], p, 'sus')
    if one_layer: groups = []
    G = len(groups)
    pr = H_.priority_list(order, groups, energy, K, 'private_first')
    E = bs_energy(net, p)
    dec = H_.Design(groups, pr, K).decoders
    aff = lambda j, b: np.min(E[dec[j], b])
    core_priv = [1 + G + k for k in order]
    grp = [j for j in pr if 1 <= j <= G]
    other = [j for j in pr if j > G and j not in core_priv]
    head = {groups[g][0]: 1 + g for g in range(G)}
    return dict(order=order, groups=groups, G=G, pr=pr, aff=aff, core_priv=core_priv, grp=grp, other=other, head=head)


def _build(lists, groups, p, core, pr, rule):
    K, B, S = p['K'], p['B'], p['Smax']
    sets = [l[:S] for l in lists]
    served = sorted(set().union(*sets))
    d = H_.Design(groups, served, K)
    d.core = [k for k in core if (1 + d.G + k) in d.col]
    d.priority = pr; d.rule = rule
    if any(set(s) != set(served) for s in sets):
        d.bsmask = np.array([[j in sets[b] for j in d.served] for b in range(B)])
        d.bs_sets = [sorted(s) for s in sets]
    return d


def cluster(net, p, one_layer=False):
    """Per BS: super-common; then the clusters (core user and its group) in decreasing affinity of the BS to the core
    user, each with the core user's private stream followed by its group-common stream; then non-core privates."""
    P = _parts(net, p, one_layer); G = P['G']
    lists = []
    for b in range(p['B']):
        lst = [0]
        for j in sorted(P['core_priv'], key=lambda j: -P['aff'](j, b)):
            lst.append(j); k = j - 1 - G
            if k in P['head']: lst.append(P['head'][k])
        lst += sorted(P['other'], key=lambda j: -P['aff'](j, b))
        lists.append(lst)
    return _build(lists, P['groups'], p, P['order'], P['pr'], 'cluster')


def cluster_1layer(net, p):
    return cluster(net, p, one_layer=True)


def coreonly(net, p, one_layer=False):
    """Per BS: super-common, core privates, group-common streams (each tier in decreasing affinity); no private
    streams of non-core UEs."""
    P = _parts(net, p, one_layer)
    lists = []
    for b in range(p['B']):
        lists.append([0] + sorted(P['core_priv'], key=lambda j: -P['aff'](j, b)) + sorted(P['grp'], key=lambda j: -P['aff'](j, b)))
    return _build(lists, P['groups'], p, P['order'], P['pr'], 'coreonly')


def coreonly_1layer(net, p):
    return coreonly(net, p, one_layer=True)


def dof_split(K, B, N, S, tie='large'):
    """(n_p, s): the number of private streams and the number of BSs carrying each, maximizing the spatial DoF
    n_p = min(K, sN, floor(B(S-1)/s)) (a stream on s BSs can be zero-forced at n_p - 1 users iff sN >= n_p; the
    private streams use n_p s of the B(S-1) slots left by the super-common stream); ties: the larger s."""
    sgn = 1 if tie == 'large' else -1
    n_p, s = max(((min(K, s * N, B * (S - 1) // s), sgn * s) for s in range(1, B + 1)))
    return n_p, sgn * s


def dof(net, p, one_layer=False, tie='large', fill=False, core_rule='sus'):
    K, B, N, S = p['K'], p['B'], p['N'], p['Smax']
    Hd = net['Hd']
    if core_rule == 'sus':
        n_p, s = dof_split(K, B, N, S, tie)
        core = H_.sus(Hd, n_p)
    else:                                            # greedy ZF selection, then the DoF split of the selected users
        sel = H_.zf_greedy(Hd, B * p['Pmax'])
        n_p, s = dof_split(len(sel), B, N, S, tie)
        core = sel[:n_p]
    energy = np.sum(np.abs(Hd) ** 2, axis=1)
    Hn = Hd / np.sqrt(energy)[:, None]
    corr = np.abs(Hn @ Hn.conj().T) ** 2
    follow = {c: [] for c in core}
    for k in range(K):
        if k not in core:
            follow[core[int(np.argmax(corr[k, core]))]].append(k)
    groups = [] if one_layer else [[c] + follow[c] for c in core if follow[c]]
    G = len(groups)
    E = bs_energy(net, p)
    sets = [[0] for _ in range(B)]
    for k in core:                                   # private streams: the s strongest BSs with a free slot
        for b in [b for b in np.argsort(-E[k], kind='stable') if len(sets[b]) < S][:s]:
            sets[b].append(1 + G + k)
    gaff = {1 + g: E[groups[g]].min(0) for g in range(G)}
    for b in range(B):                               # group-common streams: each BS fills its free slots
        for j in sorted(gaff, key=lambda j: -gaff[j][b])[:S - len(sets[b])]:
            sets[b].append(j)
    others = sorted([k for k in range(K) if k not in core], key=lambda k: (-energy[k], k))
    if fill:                                         # then private streams of the other UEs
        for b in range(B):
            for k in sorted(others, key=lambda k: -E[k, b])[:S - len(sets[b])]:
                sets[b].append(1 + G + k)
    pr = [0] + [1 + G + k for k in core] + sorted(gaff, key=lambda j: (-len(groups[j - 1]), j)) + [1 + G + k for k in others]
    return _build(sets, groups, p, core, pr, 'dof')


def dof_1layer(net, p):
    return dof(net, p, one_layer=True)


def dof_small(net, p):
    return dof(net, p, tie='small')


def dof_small_1layer(net, p):
    return dof(net, p, one_layer=True, tie='small')


def dof_fill(net, p):
    return dof(net, p, fill=True)


def dof_fill_1layer(net, p):
    return dof(net, p, one_layer=True, fill=True)


def dof_zf(net, p):
    return dof(net, p, core_rule='zf')


def dof_zf_1layer(net, p):
    return dof(net, p, one_layer=True, core_rule='zf')


def zf_layer(net, p):
    """The private layer: greedy zero-forcing user selection [Dimic-Sidiropoulos] in which the private stream of every
    selected UE is carried by its s strongest BSs that have a free fronthaul slot, for s = 1..B; the s whose selection
    has the largest ZF sum rate is kept (ties: the larger s). ZF sum rate of a selection C with carriers A_k:
    sum_k log2(1 + (P/|C|) / [(H_k H_k^H)^-1]_kk), H_k the channels of C on the antennas of A_k (P = B P_max, as in
    hrs.zf_greedy; a stream whose BSs have fewer antennas than |C| cannot be zero-forced, so the selection stops there).
    Returns (core users in selection order, carriers {k: [BSs]}, s)."""
    K, B, N, S = p['K'], p['B'], p['N'], p['Smax']
    Hd = net['Hd']; P = B * p['Pmax']
    E = bs_energy(net, p)
    order = [np.argsort(-E[k], kind='stable') for k in range(K)]

    def rate(sel, car):
        n = len(sel); r = 0.0
        for i, k in enumerate(sel):
            ant = np.concatenate([np.arange(b * N, (b + 1) * N) for b in car[k]])
            if len(ant) < n: return -np.inf
            Hk = Hd[sel][:, ant]
            try:
                g = 1.0 / max(np.real(np.linalg.inv(Hk @ Hk.conj().T)[i, i]), 1e-300)
            except np.linalg.LinAlgError:
                return -np.inf
            r += np.log2(1 + P / n * g)
        return r

    best = None
    for s in range(1, B + 1):
        sel, car, load, cur = [], {}, np.zeros(B, int), 0.0
        while True:
            cand = None
            for k in range(K):
                if k in sel: continue
                bs = [b for b in order[k] if load[b] < S - 1][:s]
                if len(bs) < s: continue
                r = rate(sel + [k], {**car, k: bs})
                if r > cur + 1e-9 and (cand is None or r > cand[0]): cand = (r, k, bs)
            if cand is None: break
            cur, k, bs = cand; sel.append(k); car[k] = bs; load[bs] += 1
        if best is None or cur >= best[0] - 1e-9:
            best = (cur, sel, car, s)
    return best[1], best[2], best[3]


def zfs(net, p, one_layer=False):
    """Discrete design ZS: private layer by zf_layer; every other UE joins the core user with the most correlated
    direct channel (groups); each BS carries the super-common stream, the private streams assigned to it, and fills
    its free slots with group-common streams in decreasing order of its direct-channel energy to the weakest member."""
    K, B, S = p['K'], p['B'], p['Smax']
    Hd = net['Hd']
    core, car, s = zf_layer(net, p)
    energy = np.sum(np.abs(Hd) ** 2, axis=1)
    Hn = Hd / np.sqrt(energy)[:, None]
    corr = np.abs(Hn @ Hn.conj().T) ** 2
    follow = {c: [] for c in core}
    for k in range(K):
        if k not in core:
            follow[core[int(np.argmax(corr[k, core]))]].append(k)
    groups = [] if one_layer else [[c] + follow[c] for c in core if follow[c]]
    G = len(groups)
    E = bs_energy(net, p)
    sets = [[0] for _ in range(B)]
    for k in core:
        for b in car[k]: sets[b].append(1 + G + k)
    gaff = {1 + g: E[groups[g]].min(0) for g in range(G)}
    for b in range(B):
        for j in sorted(gaff, key=lambda j: -gaff[j][b])[:S - len(sets[b])]:
            sets[b].append(j)
    others = sorted([k for k in range(K) if k not in core], key=lambda k: (-energy[k], k))
    pr = [0] + [1 + G + k for k in core] + sorted(gaff, key=lambda j: (-len(groups[j - 1]), j)) + [1 + G + k for k in others]
    d = _build(sets, groups, p, core, pr, 'zfs')
    d.dof = (len(core), s)
    return d


def zfs_1layer(net, p):
    return zfs(net, p, one_layer=True)


def zf_rate(Hd, sel, car, N, P):
    """ZF sum rate of the private streams of the UEs sel, the stream of UE k precoded on the antennas of the BSs car[k]:
    sum_k log2(1 + (P/|sel|) / [(H_k H_k^H)^-1]_kk), H_k the channels of sel on those antennas; -inf if some stream
    cannot be zero-forced (its BSs have fewer antennas than |sel|)."""
    n = len(sel); r = 0.0
    for i, k in enumerate(sel):
        ant = np.concatenate([np.arange(b * N, (b + 1) * N) for b in car[k]])
        if len(ant) < n: return -np.inf
        Hk = Hd[sel][:, ant]
        try:
            g = 1.0 / max(np.real(np.linalg.inv(Hk @ Hk.conj().T)[i, i]), 1e-300)
        except np.linalg.LinAlgError:
            return -np.inf
        r += np.log2(1 + P / n * g)
    return r


def sus_zf(net, p, one_layer=False, regroup=False):
    """Discrete design SZ.
    Core users and groups as in the thesis: SUS order C on the stacked direct channels (|C| = min(K, BN)), every other
    UE joins the core user with the most correlated direct channel, every core user with a follower heads a group.
    Private layer: the private streams of the first n users of C, each carried by its s strongest BSs with a free slot
    (S_max - 1 slots per BS besides the super-common stream); (n, s) maximize the ZF sum rate (zf_rate, P = B P_max)
    over n = 1..|C| and s = 1..B (ties: the larger n, then the larger s).
    Group layer: every BS fills its free slots with group-common streams in decreasing order of its direct-channel
    energy to the weakest member. Super-common stream at every BS."""
    K, B, N, S = p['K'], p['B'], p['N'], p['Smax']
    Hd = net['Hd']; P = B * p['Pmax']
    order, groups, energy, _ = H_.joint_structure(Hd, p, 'sus')
    if one_layer: groups = []
    G = len(groups)
    E = bs_energy(net, p)
    best = None
    for s in range(1, B + 1):
        load = np.zeros(B, int); car = {}
        for n, k in enumerate(order, 1):
            bs = [b for b in np.argsort(-E[k], kind='stable') if load[b] < S - 1][:s]
            if len(bs) < s: break
            car[k] = bs; load[bs] += 1
            r = zf_rate(Hd, order[:n], car, N, P)
            if best is None or r > best[0] + 1e-9 or (abs(r - best[0]) <= 1e-9 and (n, s) > (len(best[1]), best[3])):
                best = (r, order[:n], dict(car), s)
    _, core, car, s = best
    if regroup and not one_layer:                    # groups around the private-layer users
        Hn = Hd / np.sqrt(energy)[:, None]
        corr = np.abs(Hn @ Hn.conj().T) ** 2
        follow = {c: [] for c in core}
        for k in range(K):
            if k not in core:
                follow[core[int(np.argmax(corr[k, core]))]].append(k)
        groups = [[c] + follow[c] for c in core if follow[c]]
        G = len(groups)
    sets = [[0] for _ in range(B)]
    for k in core:
        for b in car[k]: sets[b].append(1 + G + k)
    gaff = {1 + g: E[groups[g]].min(0) for g in range(G)}
    for b in range(B):
        for j in sorted(gaff, key=lambda j: -gaff[j][b])[:S - len(sets[b])]:
            sets[b].append(j)
    pr = H_.priority_list(order, groups, energy, K, 'private_first')
    d = _build(sets, groups, p, core, pr, 'sus_zf')
    d.dof = (len(core), s)
    return d


def sus_zf_1layer(net, p):
    return sus_zf(net, p, one_layer=True)


def sus_zf_regroup(net, p):
    return sus_zf(net, p, regroup=True)
