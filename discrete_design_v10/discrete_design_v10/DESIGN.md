# Design return (session `cloud_medium`): the layered discrete design

**Variant `SP`** (`code/layered.py: design`), its 1-layer RSMA counterpart **`SP1`** (`code/layered.py:
design_1layer`). All runs: `results/`, made with the harness of package v10 on this machine (self-test passed, every
trace identical).

## Summary

A private layer, a group layer and the super-common stream, each with its own fronthaul rule:

* **Private layer.** The private streams of the first *n* users in the SUS order of the thesis. Each stream is carried
  by its user's *s* strongest BSs that still have a free slot. (*n*, *s*) maximize the zero-forcing (ZF) sum rate of
  the layer under the fronthaul limit, with one criterion that contains the spatial DoF (a stream on *s* BSs can be
  zero-forced only if *sN* ≥ *n*), the fronthaul (*n·s* ≤ *B*(*S*max − 1)), the array gain of cooperation, and the
  diminishing return of one more private stream.
* **Group layer.** Every UE without a private stream joins the private-layer user with the most correlated direct
  channel (the thesis' nearest-core rule). Each BS fills the slots the private layer leaves with the group-common
  streams it reaches best (its direct-channel energy to the weakest member).
* **Super-common stream** at every BS.

No test of alternatives by the AO, no tuned constant. The design costs 2–3 ms per drop and setting.

At Q = 50 against the current framework J (reference only): SP is better in E3 at every R_th (+6 to +11). In T2 it
is better at R_th ≤ 1 (+1) and worse at R_th ≥ 1.5 (−1.5 to −4.6): there it serves fewer private streams in
settings a and f (§4). Its own group layer beats its 1-layer counterpart by more than the current design's does:
J − J1 ≤ 1.3 everywhere, SP − SP1 up to +12.5 (T2) and +7.9 (E3) at R_th = 2.5.

Requirement status (details in "Results"):

| | T2 | E3 |
|---|---|---|
| R1 (Q = 50, confirm tier): SP − SP1 at R_th = 0 / 0.5 / 1 / 1.5 / 2 / 2.5 | −0.34 / +1.53 / **+1.48 / +2.09 / +4.78 / +12.50** | −1.04 / −0.79 / **+0.22 / +2.24 / +4.64 / +7.86** |
| R1 (Q = 50, 6 drops: confirm tier + drops 904–906) | +0.37 / +1.15 / **+0.67 / +2.24 / +5.38 / +14.77** | −0.42 / −0.20 / **+1.73 / +2.89 / +5.31 / +9.44** |
| R1 (Q = 5, screen tier), same R_th | 0.00 / +0.15 / **+1.27 / +2.47 / +5.43 / +18.76** | 0.00 / +0.01 / **+0.63 / +2.13 / +5.00 / +7.55** |
| For reference: SP − J (current framework), Q = 50 confirm tier | +0.90 / +1.28 / +1.03 / −1.48 / −4.59 / −3.42 | +6.36 / +7.84 / +6.38 / +6.91 / +6.49 / +11.39 |
| R2 (bench, Q = 50): SP vs sparse BF / vs stream clustering on SP's streams | **+5.66 ± 2.58** / **+23.33 ± 2.49** | — |
| R3 (bench): shares super-common / group-common | **16.3 % / 9.6 %** | — |
| R4 (smax, Q = 50) | rises with S_max, then flat (two steps of −0.3 ± 0.5) | — |
| BS-specific group-common or private streams | 18 of 144 screen runs (group-common) | 36 of 144 (private and group-common) |

## 1. The design

### 1.1 Notation

$\mathbf h_{d,k}\in\mathbb C^{1\times BN}$ the stacked direct channel of UE $k$ (noise-normalized), $\mathbf f_{b,k}$ its
block of BS $b$, $E_{b,k}=\|\mathbf f_{b,k}\|^2$. Each BS has $S^{\max}$ fronthaul slots, one of which carries the
super-common stream, so $S^{\max}-1$ remain for the private and group-common streams it carries. $M_{b,j}\in\{0,1\}$:
BS $b$ carries stream $j$.

### 1.2 Private layer

**Order.** The semi-orthogonal user selection of the thesis (eq. core_selection), on the stacked direct channels,
gives the order $c_1,c_2,\dots,c_{|\mathcal C|}$, $|\mathcal C|=\min\{K,BN\}$.

**Association for given $(n,s)$.** The private streams of $c_1,\dots,c_n$ are associated in this order: stream $c_i$
is carried by the $s$ BSs with the largest $E_{b,c_i}$ among those that still carry fewer than $S^{\max}-1$ private
streams; $\mathcal B_i(s)$ denotes this set. (If fewer than $s$ BSs have a free slot, $n$ is not admissible.)

**Choice of $(n,s)$.** With $\mathbf H_i\in\mathbb C^{n\times|\mathcal B_i|N}$ the direct channels of
$c_1,\dots,c_n$ on the antennas of the BSs in $\mathcal B_i(s)$, the zero-forcing sum rate of the layer, with the
network power $P=BP^{\max}$ shared equally, is

$$
R_{\mathrm{ZF}}(n,s)=\sum_{i=1}^{n}\log_2\!\Bigl(1+\frac{P}{n}\,\frac{1}{[(\mathbf H_i\mathbf H_i^H)^{-1}]_{ii}}\Bigr),
\qquad R_{\mathrm{ZF}}(n,s)=-\infty\ \text{if } sN<n ,
$$

where $1/[(\mathbf H_i\mathbf H_i^H)^{-1}]_{ii}$ is the gain of the ZF beam of $c_i$ restricted to the antennas of
its BSs (the norm of the part of its channel there that is orthogonal to the other $n-1$ users' channels there). The
design takes

$$
(n^\star,s^\star)=\arg\max_{1\le n\le|\mathcal C|,\ 1\le s\le B} R_{\mathrm{ZF}}(n,s)
$$

(ties: the larger $n$, then the larger $s$), and serves the private streams of $c_1,\dots,c_{n^\star}$ with the
association $\mathcal B_i(s^\star)$.

### 1.3 Group layer

Every UE $k$ without a private stream joins the private-layer user with the most correlated direct channel
(eq. nearest_core of the thesis, over $\{c_1,\dots,c_{n^\star}\}$); every private-layer user with at least one
follower heads the group $\mathcal U_a$. Every BS $b$ fills the slots left by the private layer with group-common
streams in decreasing order of $\min_{k\in\mathcal U_g}E_{b,k}$ (its channel to the weakest member, the member that
limits the group-common rate).

### 1.4 Super-common stream

Carried by every BS.

### 1.5 Algorithm

```
Input: direct channels h_{d,k}, B, N, K, S_max, P_max
1  c_1..c_|C| <- SUS order on the stacked direct channels, |C| = min(K, BN)
2  for s = 1..B:
3      load_b <- 0 for every BS b
4      for n = 1..|C|:
5          B_n(s) <- the s BSs with largest E_{b,c_n} among {b : load_b < S_max - 1}; stop if fewer than s
6          load_b <- load_b + 1 for b in B_n(s)
7          R(n, s) <- R_ZF of c_1..c_n with carriers B_1(s)..B_n(s)       (-inf if sN < n)
8  (n*, s*) <- argmax R(n, s)  (ties: larger n, then larger s)
9  private layer: stream of c_i on B_i(s*), i = 1..n*
10 groups: every other UE joins argmax_{a in c_1..c_n*} |h_k h_a^H|^2 / (|h_k|^2 |h_a|^2); heads with followers
11 every BS b: super-common stream; its private streams; then group-common streams in decreasing
   min_{k in U_g} E_{b,k} until S_max streams
```

Complexity: $B\cdot|\mathcal C|$ evaluations of $R_{\mathrm{ZF}}$, each $n$ inversions of $n\times n$ matrices, i.e.
$O(B|\mathcal C|^5)$ in the worst case; 2–3 ms per drop and setting in the test grid (`cpu_design` in the records),
about 1 s with 25 BSs and 60 UEs.

### 1.6 Why this design

* **What the private layer is for.** The private streams carry the spatially multiplexed rate. How many can be
  multiplexed, and from how many BSs, is decided by the fronthaul. A stream carried by $s$ BSs is precoded on $sN$
  antennas, so it can be zero-forced at the other $n-1$ private users only if $sN\ge n$. And the private layer
  occupies $n s$ of the $B(S^{\max}-1)$ slots left by the super-common stream. Both limits are inside
  $R_{\mathrm{ZF}}$: the first makes it $-\infty$, and the second makes $(n,s)$ inadmissible because the association
  runs out of slots. So does the array gain that extra cooperation buys (more antennas per stream, larger ZF gain),
  and the effect that stops greedy ZF user selection: a user whose channel is nearly in the span of the others adds
  little and costs every other stream ZF gain. The ZF sum rate at the direct channels is the standard model-based
  measure of a private layer (it is the criterion of greedy ZF user selection, `hrs.zf_greedy`, with the user order of
  SUS). It is evaluated in closed form, never by the AO, and the AO then optimizes the precoders freely within the
  chosen association.
* **Association.** Each private stream goes to its user's strongest BSs. In a large network these are the
  neighbouring BSs, so a core UE is associated with only some BSs; *s* is chosen by the rate criterion, not fixed.
  When every BS–UE link is strong and the fronthaul allows it, full cooperation maximizes $R_{\mathrm{ZF}}$, and the
  design associates every core UE with every BS by itself. This happens in every T2 run ($s^\star = B = 2$), and in
  E3 settings a, e and f and 15 of 16 drops of c ($s^\star = 3$). In E3 b and d, where zero-forcing more users with fewer BSs each wins,
  $s^\star = 2$. With 9–25 BSs scattered over 300–1000 m, $s^\star$ is 5–7 and each stream is carried by exactly
  its nearest $s^\star$ BSs (`results/large_network_*.txt`, `tools/large_network.py`). Different BSs then carry
  different private streams: some UEs are core UEs of one group of BSs, others of another, and a BS's private streams
  are those of the core UEs it is among the strongest BSs for.
* **What the group layer is for.** UEs without a private stream must get their QoS rate from common streams. Without
  groups, only the super-common stream reaches them, and its rate is limited by the weakest of all *K* UEs. A
  group-common stream reaches only its members, who are spatially correlated with their head and served by it, and
  several groups multiplex their common messages. The group layer gets exactly the fronthaul the private layer does
  not use. Each BS picks the groups it reaches best, so the BSs carry different group-common streams when the free
  slots do not admit every group everywhere (T2 settings c and e; E3 b, c, d and e in some drops).
* **Why the group layer beats the 1-layer counterpart (R1).** The private layer, its association and the
  super-common stream do not depend on the groups. The group layer only takes slots that the 1-layer design leaves
  empty. Every stream set of SP1 is therefore a subset of SP's with the same carriers, and a group-common stream with
  zero power is the same as not serving it. So the optimum of SP is at least that of SP1 at every R_th. The group
  layer pays off when the QoS load of the UEs without private streams is large (high R_th, many such UEs), and costs
  nothing when it is not (R_th → 0: the AO can turn the group streams down, so a tie).
* **Fronthaul (R4).** A larger S_max relaxes only the slot constraint of the association. More (n, s) become
  admissible and the streams get their strongest BSs more often, so the private layer grows, with more streams or
  more BSs per stream, until it reaches its ZF-best form. From then on every additional slot goes to the group layer,
  and once every group is served the design no longer changes (the rate is flat in S_max). The stream sets are not
  nested in S_max, though: when the private layer grows, the groups around it are formed anew (in 13 % of the S_max
  steps of T2 and 22 % of E3 over K ∈ {8, 12, 16}, N ∈ {3, 5}, 16 drops). So R4 is a measured property (§2.3), not a
  structural one.

### 1.7 No empirical tuning

The only inputs are the channels and the settings ($B$, $N$, $K$, $S^{\max}$, $P^{\max}$). The rules are: the SUS
order of the thesis; the ZF sum rate with equal power (a model expression, not a fitted one; $P=BP^{\max}$ is the
network power); "strongest BS" by direct-channel energy; and the thesis' nearest-core grouping. Ties are broken
deterministically. Nothing is chosen with the test drops, the R_th grid or the scenarios: there is no threshold, no
weight and no count, and no AO evaluation of alternatives. The same code runs for any B, any BS positions and any
setting.

### 1.8 The 1-layer RSMA counterpart (SP1)

The same order, the same private layer $(n^\star,s^\star)$ and association, the super-common stream at every BS, and
no group layer. The UEs without private streams are served by the super-common stream only.

### 1.9 Evaluation, search and refinement

Unchanged: `discrete.run` with the default options (AO start `core`, search R trimmed, Q = 50, D_max = 10 m). No
switch is used by SP. (`algo.init_core_mask`, a start for BS-specific associations, was added as an option for the
tried variant DFm. It is off by default and changed nothing measurable, so it is not part of the design.)

## 2. Results

All rates in bit/s/Hz, zero when the QoS shortfall exceeds 1e-2; differences paired (mean ± SE). Report outputs
are pasted as printed by `tools/report.py`; the files are in `results/` (`report_*.txt`, `*report*.json`).

### 2.1 R1: SP against its 1-layer counterpart SP1

**Confirm tier (Q = 50, both scenarios, 6 settings × 6 R_th × 3 drops)**, the primary result:

| R_th | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 |
|---|---|---|---|---|---|---|
| T2: SP − SP1 | −0.34 ± 0.65 | +1.53 ± 1.06 | +1.48 ± 0.90 | +2.09 ± 1.63 | +4.78 ± 2.52 | +12.50 ± 4.64 |
| E3: SP − SP1 | −1.04 ± 0.67 | −0.79 ± 0.59 | +0.22 ± 0.86 | +2.24 ± 1.11 | +4.64 ± 2.18 | +7.86 ± 3.37 |
| T2: SP − J (reference) | +0.90 ± 1.36 | +1.28 ± 1.73 | +1.03 ± 1.53 | −1.48 ± 1.89 | −4.59 ± 1.48 | −3.42 ± 1.38 |
| E3: SP − J (reference) | +6.36 ± 2.02 | +7.84 ± 2.18 | +6.38 ± 2.45 | +6.91 ± 2.77 | +6.49 ± 3.02 | +11.39 ± 5.39 |
| T2: J − J1, for comparison | +0.92 | +0.55 | +0.59 | +0.02 | +0.32 | +1.25 |
| E3: J − J1, for comparison | +0.41 | +0.70 | +0.60 | +0.79 | −0.15 | +0.27 |

At every R_th from 1 to 2.5 the mean difference is positive in both scenarios. It is a tie at R_th = 0 and 0.5
(within one SE, except E3 at 0, −1.04 ± 0.67). In settings where the two designs are identical (T2 b, c, d; E3 c,
e: no slot left for a group, §4) the runs are identical; the margins come from T2 a, e, f and E3 b, d, f. At Q = 50
the paired differences are noisy, because SP and SP1 have different evaluators, so the deployment search follows
different paths: single paired runs differ by up to ±10 at low R_th. E3 at R_th = 1 (+0.22 ± 0.86) is therefore
not significant on 3 drops; see the 6-drop extension below and the screen tier.

**Extension to 6 drops (Q = 50, drops 901–906, 36 runs per R_th and scenario).** To resolve the search noise, the
confirm tier of SP and SP1 was repeated on drops 904–906 (`--tier confirm --drops 6`; the Q = 50 result files hold
all 6 drops). Paired SP − SP1 on all 6 drops, computed by `tools/r1_pair.py` (the statistics of `tools/report.py`,
without intersecting with the JH/J1 runs, which exist for drops 901–903 only):

| R_th | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 |
|---|---|---|---|---|---|---|
| T2: SP − SP1 (6 drops) | +0.37 ± 0.53 | +1.15 ± 0.73 | +0.67 ± 0.57 | +2.24 ± 1.00 | +5.38 ± 1.85 | +14.77 ± 4.60 |
| E3: SP − SP1 (6 drops) | −0.42 ± 0.59 | −0.20 ± 0.50 | +1.73 ± 0.87 | +2.89 ± 1.04 | +5.31 ± 1.87 | +9.44 ± 2.73 |
| T2, drops 904–906 only | +1.08 ± 0.82 | +0.77 ± 1.02 | −0.15 ± 0.67 | +2.38 ± 1.20 | +5.98 ± 2.77 | +17.05 ± 8.06 |
| E3, drops 904–906 only | +0.21 ± 0.97 | +0.39 ± 0.81 | +3.24 ± 1.44 | +3.53 ± 1.79 | +5.97 ± 3.09 | +11.02 ± 4.37 |

On 6 drops SP beats SP1 at every R_th from 1 to 2.5 in both scenarios: by 2 SE or more at R_th ≥ 1.5 in both, and at
R_th = 1 in E3. **The weakest point is T2 at R_th = 1** (+0.67 ± 0.57; drops 904–906 alone −0.15 ± 0.67). There, the
group layer exists in 3 of the 6 settings (a, e, f) and is worth about 1 bit/s/Hz (+1.27 ± 0.61 at Q = 5, where the
search noise is small), while one Q = 50 run of SP and one of SP1 differ by several bit/s/Hz through the search
alone. R_th = 0 and 0.5 are ties in both scenarios. Per-setting tables: `results/report_SP_R1_Q50_6drops.txt` and
`results/report_SP_R1_Q50_drops904-906.txt`.

**Screen tier (Q = 5, 4 drops)**: near the initial deployment the search noise is small, and the pairing shows the
effect of the group layer alone:

```
-- T2, all settings
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0   24    108.40    108.40    104.77    104.76    104.41    -0.00± 0.00       +3.63± 1.73      24/24/24/24/24
    0.5   24    105.50    105.35    102.52    102.51    102.16    +0.15± 0.10       +2.98± 1.73      24/24/24/24/24
      1   24     97.34     96.08     95.64     95.63     95.34    +1.27± 0.61       +1.70± 1.39      24/24/24/24/24
    1.5   24     84.93     82.47     87.23     87.23     86.94    +2.47± 1.32       -2.30± 1.07      24/24/24/24/24
      2   24     65.95     60.52     70.13     70.40     69.41    +5.43± 2.31       -4.17± 1.53      20/20/20/20/20
    2.5   24     52.39     33.63     55.88     55.88     54.42   +18.76± 6.83       -3.49± 1.73      16/12/16/16/16
    all  144     85.75     81.07     86.03     86.07     85.45    +4.68± 1.32       -0.27± 0.67      132/128/132/132/132
-- E3, all settings
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0   24    134.73    134.73    128.91    128.82    128.91    +0.00± 0.00       +5.83± 1.75      24/24/24/24/24
    0.5   24    129.54    129.54    122.90    122.83    122.90    +0.01± 0.01       +6.64± 1.89      24/24/24/24/24
      1   24    120.25    119.62    113.25    113.24    113.25    +0.63± 0.29       +7.01± 2.11      24/24/24/24/24
    1.5   24    108.21    106.07    101.61    101.60    101.71    +2.13± 1.06       +6.60± 2.18      24/24/24/24/24
      2   24     90.32     85.32     83.68     83.68     83.66    +5.00± 2.08       +6.65± 2.43      20/20/20/20/20
    2.5   24     80.95     73.39     70.64     70.64     70.63    +7.55± 3.01      +10.30± 4.08      20/20/16/16/16
    all  144    110.67    108.11    103.50    103.47    103.51    +2.55± 0.67       +7.17± 1.02      136/136/132/132/132
```

(`results/report_SP_Q5.txt` has the tables per setting.) The J columns are the current framework (JH or J1 by the
stream-set test), derived from the JH and J1 runs.

<details><summary>Confirm tier, tools/report.py SP --by-setting (Q = 50, drops 901–903, with the JH/J1 reference runs), as printed</summary>

```
===== scenario T2, Q = 50: 108 common runs (SP: 216, SP1: 216, JH: 108, J1: 108)
-- T2, all settings
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0   18    125.58    125.92    124.68    123.91    123.76    -0.34± 0.65       +0.90± 1.36      18/18/18/18/18
    0.5   18    122.68    121.15    121.41    120.44    120.85    +1.53± 1.06       +1.28± 1.73      18/18/18/18/18
      1   18    114.57    113.09    113.54    113.02    112.95    +1.48± 0.90       +1.03± 1.53      18/18/18/18/18
    1.5   18    102.13    100.04    103.61    103.47    103.59    +2.09± 1.63       -1.48± 1.89      18/18/18/18/18
      2   18     88.95     84.17     93.54     94.06     93.22    +4.78± 2.52       -4.59± 1.48      18/18/18/18/18
    2.5   18     73.44     60.94     76.86     77.34     75.61   +12.50± 4.64       -3.42± 1.38      15/15/15/15/15
    all  108    104.56    100.89    105.61    105.37    105.00    +3.67± 1.02       -1.05± 0.67      105/105/105/105/105
-- T2, setting a
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    104.72    105.55    110.24    106.84    106.33    -0.83± 1.10       -5.52± 1.52      3/3/3/3/3
    0.5    3    103.38    104.00    111.23    107.15    110.08    -0.63± 0.84       -7.85± 1.42      3/3/3/3/3
      1    3    104.80    103.39    111.42    107.38    109.07    +1.41± 0.36       -6.62± 1.36      3/3/3/3/3
    1.5    3     99.66    100.21    108.84    105.70    108.44    -0.55± 1.24       -9.18± 3.33      3/3/3/3/3
      2    3     98.96     96.12    110.08    108.48    108.73    +2.84± 1.76      -11.12± 3.69      3/3/3/3/3
    2.5    3     97.27     91.51    106.47    109.33    104.77    +5.76± 0.42       -9.20± 4.38      3/3/3/3/3
    all   18    101.47    100.13    109.71    107.48    107.90    +1.33± 0.68       -8.25± 1.09      18/18/18/18/18
-- T2, setting b
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3     84.75     84.75     84.75     84.75     84.75    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    0.5    3     78.06     78.06     78.06     78.06     78.06    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      1    3     71.01     71.01     71.01     71.01     71.01    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    1.5    3     62.39     62.39     62.39     62.39     62.39    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      2    3     53.41     53.41     53.41     53.41     53.41    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    2.5    3     42.77     42.77     42.77     42.77     42.77    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    all   18     65.40     65.40     65.40     65.40     65.40    +0.00± 0.00       +0.00± 0.00      18/18/18/18/18
-- T2, setting c
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    117.51    117.51    117.51    117.51    117.51    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    0.5    3    111.22    111.22    111.22    111.22    111.22    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      1    3     95.29     95.29     95.29     95.29     95.29    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    1.5    3     74.50     74.50     74.50     74.50     74.50    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      2    3     49.93     49.93     49.93     49.93     49.93    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    2.5    3      0.00      0.00      0.00      0.00      0.00    +0.00± 0.00       +0.00± 0.00      0/0/0/0/0
    all   18     74.74     74.74     74.74     74.74     74.74    +0.00± 0.00       +0.00± 0.00      15/15/15/15/15
-- T2, setting d
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    129.06    129.06    129.06    129.06    129.06    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    0.5    3    125.63    125.63    125.63    125.63    125.63    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      1    3    117.95    117.95    117.95    117.95    117.95    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    1.5    3    108.74    108.74    108.74    108.74    108.74    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      2    3     94.69     94.69     94.69     94.69     94.69    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    2.5    3     81.65     81.65     81.65     81.65     81.65    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    all   18    109.62    109.62    109.62    109.62    109.62    +0.00± 0.00       +0.00± 0.00      18/18/18/18/18
-- T2, setting e
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    161.40    161.30    158.62    158.62    158.62    +0.10± 4.14       +2.78± 5.52      3/3/3/3/3
    0.5    3    163.57    153.40    155.20    155.20    155.20   +10.17± 3.28       +8.38± 5.45      3/3/3/3/3
      1    3    152.73    148.47    143.46    143.46    143.46    +4.26± 5.57       +9.26± 2.34      3/3/3/3/3
    1.5    3    134.58    131.26    137.38    137.38    137.38    +3.31± 7.56       -2.80± 6.98      3/3/3/3/3
      2    3    117.69    110.11    127.31    127.31    127.31    +7.58±10.14       -9.61± 3.77      3/3/3/3/3
    2.5    3    104.38     80.88    110.35    110.35    110.35   +23.50±12.11       -5.97± 3.59      3/3/3/3/3
    all   18    139.06    130.90    138.72    138.72    138.72    +8.15± 3.25       +0.34± 2.39      18/18/18/18/18
-- T2, setting f
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    156.05    157.35    147.89    146.67    146.28    -1.30± 1.48       +8.16± 3.38      3/3/3/3/3
    0.5    3    154.24    154.60    147.09    145.38    144.94    -0.36± 1.09       +7.16± 5.83      3/3/3/3/3
      1    3    145.64    142.44    142.09    143.02    140.89    +3.21± 1.41       +3.55± 6.65      3/3/3/3/3
    1.5    3    132.92    123.13    129.85    132.11    130.11    +9.80± 6.00       +3.07± 8.81      3/3/3/3/3
      2    3    119.01    100.77    125.83    130.55    125.27   +18.24± 9.38       -6.83± 4.00      3/3/3/3/3
    2.5    3    114.56     68.85    119.92    119.92    114.14   +45.72± 9.14       -5.36± 5.05      3/3/3/3/3
    all   18    137.07    124.52    135.45    136.27    133.60   +12.55± 4.44       +1.63± 2.46      18/18/18/18/18
-- SP: mean served streams 7.89, group-common served 1.44, private served 5.44; BS-specific in 12/108 runs; per run, carried by a subset of the BSs: group-common 0.22, private 0.00
-- J: stream-set test made in 36/108 runs, 1-layer kept in 12
-- CPU per run (s): SP 39.0, SP1 19.8, J 39.7, JH 43.2, J1 29.1; SP discrete design 0.002

===== scenario E3, Q = 50: 108 common runs (SP: 216, SP1: 216, JH: 108, J1: 108)
-- E3, all settings
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0   18    145.86    146.91    139.50    139.60    139.09    -1.04± 0.67       +6.36± 2.02      18/18/18/18/18
    0.5   18    140.86    141.65    133.02    133.37    132.32    -0.79± 0.59       +7.84± 2.18      18/18/18/18/18
      1   18    131.37    131.15    124.99    126.26    124.38    +0.22± 0.86       +6.38± 2.45      18/18/18/18/18
    1.5   18    119.69    117.46    112.79    112.47    112.00    +2.24± 1.11       +6.91± 2.77      18/18/18/18/18
      2   18    101.94     97.30     95.46     95.46     95.61    +4.64± 2.18       +6.49± 3.02      16/16/16/16/16
    2.5   18     89.82     81.96     78.43     79.12     78.16    +7.86± 3.37      +11.39± 5.39      15/15/12/12/12
    all  108    121.59    119.40    114.03    114.38    113.59    +2.19± 0.77       +7.56± 1.28      103/103/100/100/100
-- E3, setting a
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    142.06    143.07    137.13    137.13    137.13    -1.01± 0.77       +4.94± 4.68      3/3/3/3/3
    0.5    3    139.15    142.02    135.09    135.09    135.09    -2.88± 2.55       +4.06± 2.11      3/3/3/3/3
      1    3    140.52    141.24    143.92    143.92    143.92    -0.72± 0.47       -3.40± 4.18      3/3/3/3/3
    1.5    3    138.64    139.16    138.21    138.21    138.21    -0.52± 0.43       +0.43± 3.75      3/3/3/3/3
      2    3    138.43    139.23    143.74    143.74    143.74    -0.80± 1.20       -5.31± 5.36      3/3/3/3/3
    2.5    3    137.21    137.19    138.79    138.79    138.79    +0.02± 1.23       -1.59± 4.81      3/3/3/3/3
    all   18    139.33    140.32    139.48    139.48    139.48    -0.99± 0.50       -0.15± 1.72      18/18/18/18/18
-- E3, setting b
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3     87.30     86.64     76.67     76.67     76.67    +0.66± 0.31      +10.63± 2.43      3/3/3/3/3
    0.5    3     85.82     84.64     68.77     68.77     68.77    +1.17± 0.40      +17.04± 1.43      3/3/3/3/3
      1    3     76.73     76.93     60.77     60.77     60.77    -0.20± 0.97      +15.96± 0.20      3/3/3/3/3
    1.5    3     67.30     66.19     51.52     51.52     51.52    +1.10± 1.89      +15.78± 3.42      3/3/3/3/3
      2    3     58.79     54.88     42.80     42.80     42.80    +3.91± 4.03      +15.98± 4.58      3/3/3/3/3
    2.5    3     52.63     41.27      0.00      0.00      0.00   +11.36± 5.76      +52.63± 5.01      3/3/0/0/0
    all   18     71.43     68.43     50.09     50.09     50.09    +3.00± 1.41      +21.34± 3.61      18/18/15/15/15
-- E3, setting c
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    122.35    122.35    122.35    122.35    122.35    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    0.5    3    105.73    105.73    105.73    105.73    105.73    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      1    3     85.46     85.46     85.46     85.46     85.46    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    1.5    3     64.42     64.42     64.42     64.42     64.42    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      2    3     12.73     12.73     12.73     12.73     12.73    +0.00± 0.00       +0.00± 0.00      1/1/1/1/1
    2.5    3      0.00      0.00      0.00      0.00      0.00    +0.00± 0.00       +0.00± 0.00      0/0/0/0/0
    all   18     65.12     65.12     65.12     65.12     65.12    +0.00± 0.00       +0.00± 0.00      13/13/13/13/13
-- E3, setting d
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    127.90    128.81    124.94    124.94    124.94    -0.91± 0.91       +2.96± 2.96      3/3/3/3/3
    0.5    3    120.82    121.16    116.46    116.46    116.46    -0.34± 0.34       +4.36± 4.36      3/3/3/3/3
      1    3    115.41    115.73    107.64    107.64    107.64    -0.32± 0.32       +7.77± 7.77      3/3/3/3/3
    1.5    3    106.17    105.19     95.37     95.37     95.37    +0.99± 0.99      +10.80±10.80      3/3/3/3/3
      2    3     94.71     93.46     81.89     81.89     81.89    +1.26± 1.26      +12.83±12.83      3/3/3/3/3
    2.5    3     81.36     78.96     67.18     67.18     67.18    +2.40± 2.40      +14.18±14.18      3/3/3/3/3
    all   18    107.73    107.22     98.91     98.91     98.91    +0.51± 0.51       +8.82± 3.49      18/18/18/18/18
-- E3, setting e
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    195.79    195.79    195.79    195.79    195.79    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    0.5    3    191.44    191.44    191.44    191.44    191.44    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      1    3    173.39    173.39    173.39    173.39    173.39    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    1.5    3    150.58    150.58    150.58    150.58    150.58    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
      2    3    124.31    124.31    124.31    124.31    124.31    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    2.5    3     98.11     98.11     98.11     98.11     98.11    +0.00± 0.00       +0.00± 0.00      3/3/3/3/3
    all   18    155.60    155.60    155.60    155.60    155.60    +0.00± 0.00       +0.00± 0.00      18/18/18/18/18
-- E3, setting f
   R_th    n        SP       SP1         J        JH        J1            SP-SP1              SP-J   QoS met SP/SP1/J/JH/J1
      0    3    199.77    204.78    180.11    180.70    177.67    -5.01± 3.35      +19.66± 5.08      3/3/3/3/3
    0.5    3    202.19    204.88    180.64    182.69    176.42    -2.68± 2.02      +21.56± 2.28      3/3/3/3/3
      1    3    196.72    194.13    178.75    186.37    175.12    +2.59± 5.72      +17.97± 4.95      3/3/3/3/3
    1.5    3    191.06    179.19    176.63    174.75    171.90   +11.87± 1.35      +14.44±10.22      3/3/3/3/3
      2    3    182.69    159.19    167.27    167.27    168.19   +23.50± 1.96      +15.42± 5.72      3/3/3/3/3
    2.5    3    169.62    136.23    166.51    170.67    164.90   +33.39± 9.96       +3.11±11.35      3/3/3/3/3
    all   18    190.34    179.74    174.98    177.07    172.37   +10.61± 3.81      +15.36± 2.90      18/18/18/18/18
-- SP: mean served streams 8.89, group-common served 1.00, private served 6.89; BS-specific in 24/108 runs; per run, carried by a subset of the BSs: group-common 0.22, private 1.06
-- J: stream-set test made in 18/108 runs, 1-layer kept in 9
-- CPU per run (s): SP 67.9, SP1 39.8, J 59.5, JH 53.5, J1 59.8; SP discrete design 0.003

written /home/user/Thesis-Test-medium/discrete_design_v10/discrete_design_v10/results/SP__report.json
```
</details>

### 2.2 R2 and R3: the benchmark setting (T2, K = 16, N = 3, S_max = 12, R_th = 2, Q = 50)

```
===== benchmark setting (T2, K = 16, N = 3, S_max = 12, R_th = 2, Q = 50): 15 common runs (SP: 15, sparse_bf: 15, SP@SC: 15, J: 15, J@SC: 15)
sparse_bf, J (the current framework) and J@SC (stream clustering on its streams): stored runs of Fig. 5.10
-- by M
      M    n        SP sparse_bf     SP@SC         J      J@SC      SP-sparse_bf          SP-SP@SC   QoS met SP/sparse_bf/SP@SC/J/J@SC
    100    5     71.99     63.08     43.60     73.51     44.18    +8.92± 4.20      +28.39± 3.17      5/5/5/5/5
    300    5     78.76     75.09     57.98     81.95     54.38    +3.67± 5.50      +20.78± 5.33      5/5/5/5/5
    500    5     85.43     81.03     64.61     88.72     64.51    +4.40± 4.25      +20.82± 4.16      5/5/5/5/5
    all   15     78.73     73.07     55.40     81.39     54.36    +5.66± 2.58      +23.33± 2.49      15/15/15/15/15
-- SP: mean served streams 11.00, group-common served 4.20, private served 5.80; BS-specific in 0/15 runs; per run, carried by a subset of the BSs: group-common 0.00, private 0.00
-- SP: shares of the sum rate over its 15 runs meeting the QoS: super-common 16.3%, group-common 9.6%, private 74.0%

written /home/user/Thesis-Test-medium/discrete_design_v10/discrete_design_v10/results/SP__report_bench.json
```

* **R2.** SP 78.73 against sparse beamforming 73.07: **+5.66 ± 2.58**, above it at every M (+8.92 / +3.67 /
  +4.40). Against stream clustering on SP's own streams (SP@SC): 55.40, **+23.33 ± 2.49**. For reference, the current
  framework J: 81.39. In 12 of 15 drops SP serves exactly J's 6 private streams; in drop 802 (all M) its ZF criterion
  serves 5, which costs about 13 at M = 100 and 300 there.
* **R3.** Shares of the sum rate over the 15 runs: **super-common 16.3 %, group-common 9.6 %**, private 74.0 %
  (Fig. 5.7 of the thesis at R_th = 2: 16.7 % and 10.8 %).

### 2.3 R4: sum rate versus S_max (T2, N = 3, M = 300, Q = 50)

```
===== sum rate versus S_max (T2, N = 3, M = 300, Q = 50): 126 runs of SP
-- K = 8: mean sum rate over the drops (QoS met / runs)
  R_th          S=4          S=6          S=8         S=10         S=12         S=14         S=16
     0   83.82  2/2   101.66  2/2   100.65  2/2   100.65  2/2   100.65  2/2   100.65  2/2   100.65  2/2 
     1   79.28  2/2    97.24  2/2    99.39  2/2    99.39  2/2    99.39  2/2    99.39  2/2    99.39  2/2 
     2   69.89  2/2    90.46  2/2    91.69  2/2    91.69  2/2    91.69  2/2    91.69  2/2    91.69  2/2 
   all   77.66  6/6    96.45  6/6    97.24  6/6    97.24  6/6    97.24  6/6    97.24  6/6    97.24  6/6 
   change to the next S_max (paired over R_th and drops): 4->6 +18.79±0.82, 6->8 +0.79±0.85, 8->10 +0.00±0.00, 10->12 +0.00±0.00, 12->14 +0.00±0.00, 14->16 +0.00±0.00
-- K = 12: mean sum rate over the drops (QoS met / runs)
  R_th          S=4          S=6          S=8         S=10         S=12         S=14         S=16
     0   84.37  2/2   113.44  2/2   122.98  2/2   123.07  2/2   121.72  2/2   121.72  2/2   121.72  2/2 
     1   68.85  2/2   106.67  2/2   111.13  2/2   116.95  2/2   117.14  2/2   117.14  2/2   117.14  2/2 
     2   51.03  2/2    82.86  2/2    96.87  2/2   104.38  2/2   104.64  2/2   104.64  2/2   104.64  2/2 
   all   68.08  6/6   100.99  6/6   110.33  6/6   114.80  6/6   114.50  6/6   114.50  6/6   114.50  6/6 
   change to the next S_max (paired over R_th and drops): 4->6 +32.91±1.76, 6->8 +9.33±2.72, 8->10 +4.47±2.51, 10->12 -0.30±0.49, 12->14 +0.00±0.00, 14->16 +0.00±0.00
-- K = 16: mean sum rate over the drops (QoS met / runs)
  R_th          S=4          S=6          S=8         S=10         S=12         S=14         S=16
     0   84.23  2/2   116.64  2/2   122.32  2/2   123.59  2/2   122.63  2/2   122.63  2/2   122.63  2/2 
     1   63.37  2/2    94.37  2/2   112.06  2/2   110.32  2/2   115.43  2/2   115.43  2/2   115.43  2/2 
     2    0.00  0/2    50.59  2/2    66.83  2/2    76.74  2/2    91.99  2/2    91.99  2/2    91.99  2/2 
   all   49.20  4/6    87.20  6/6   100.40  6/6   103.55  6/6   110.02  6/6   110.02  6/6   110.02  6/6 
   change to the next S_max (paired over R_th and drops): 4->6 +38.00±3.99, 6->8 +13.20±3.19, 8->10 +3.15±2.54, 10->12 +6.47±4.07, 12->14 +0.00±0.00, 14->16 +0.00±0.00

written /home/user/Thesis-Test-medium/discrete_design_v10/discrete_design_v10/results/SP__report_smax.json
```

The rate rises with S_max for every K. Once the private layer has its ZF-best size and every group is served, the
design no longer changes, and the rows are flat and identical (K = 8 from S_max = 8, K = 12 and 16 from S_max = 12).
The two negative steps, −0.27 ± 0.47 at K = 8 (6 → 8) and −0.30 ± 0.49 at K = 12 (10 → 12), are well within one SE
(the search noise of a changed stream set).

### 2.4 Where the BSs carry different streams (design requirement)

| | carried by a subset of the BSs | runs (confirm tier) | runs (screen tier) |
|---|---|---|---|
| T2 | group-common streams (setting e in the confirm drops, c and e over drops 901–916: each BS takes the groups it reaches best when the free slots do not admit all) | 12 / 108 | 18 / 144 |
| E3 | private streams on 2 of the 3 BSs (s* = 2: settings b and d) and group-common streams (confirm drops: b 18, d 6 runs) | 24 / 108 | 36 / 144 |

Over the 16 drops 901–916 (design only): T2 c 3/16 and e 14/16 drops with BS-specific group-common streams; E3 b
15/16 and d 10/16 drops with BS-specific private and group-common streams, E3 c 1/16, e 3/16. Every core UE is
associated with all BSs in T2 (s* = 2), with all 3 in E3 a, c, e, f, and with 2 of 3 in E3 b (15/16 drops) and
d (10/16). Large networks: §1.6 and `results/large_network_*.txt`.

### 2.5 Cost

CPU per run (confirm tier, Q = 50, mean over the runs; the report): T2 SP 39.0 s, SP1 19.8 s, J 39.7 s (JH 43.2, J1
29.1); E3 SP 67.9 s, SP1 39.8 s, J 59.5 s (JH 53.5, J1 59.8). The discrete design itself takes 2–3 ms.

## 3. Other methods tried

All in `code/perbs.py` (variants in `harness/variants.py`), in the order they were tried. Each was screened with
`--tier screen` (Q = 5) together with its 1-layer counterpart (the same function with `one_layer=True`). The table is
`results/summary_tried_Q5.txt` (`tools/summary_tried.py`). Per R_th: the mean rate of X, then X − X1 ± SE, then X − J,
paired over the 144 screen runs of a scenario (DG: T2 settings b and d only, 48 runs; DFm and ZS: E3 only).

```
R_th:              0                        0.5                      1                        1.5                      2                        2.5
  SP T2 n=144 |  108.4 -0.00±0.00  +3.63 |  105.5 +0.15±0.10  +2.98 |   97.3 +1.27±0.61  +1.70 |   84.9 +2.47±1.32  -2.30 |   66.0 +5.43±2.31  -4.17 |   52.4 +18.76±6.83  -3.49
  SP E3 n=144 |  134.7 +0.00±0.00  +5.83 |  129.5 +0.01±0.01  +6.64 |  120.3 +0.63±0.29  +7.01 |  108.2 +2.13±1.06  +6.60 |   90.3 +5.00±2.08  +6.65 |   80.9 +7.55±3.01 +10.30
  TA T2 n=144 |   91.1 +0.35±0.30 -13.66 |   90.3 +0.35±0.30 -12.18 |   87.3 +0.29±0.29  -8.36 |   80.8 +0.29±0.22  -6.47 |   70.7 +0.99±0.59  +0.61 |   62.3 +1.46±1.03  +6.45
  TA E3 n=144 |  114.1 -0.08±0.08 -14.80 |  114.0 -0.07±0.06  -8.85 |  112.3 -0.01±0.01  -0.95 |  107.1 -0.11±0.10  +5.52 |  101.5 +0.01±0.01 +17.86 |   94.7 +0.01±0.01 +24.02
  TB T2 n=144 |   85.1 -5.68±2.55 -19.69 |   84.6 -5.39±2.34 -17.92 |   82.1 -4.87±2.74 -13.51 |   78.5 -1.97±1.79  -8.72 |   70.5 +0.80±2.99  +0.42 |   60.0 -0.90±2.11  +4.09
  TB E3 n=144 |  113.2 -1.03±0.95 -15.74 |  111.8 -2.36±0.86 -11.14 |  110.0 -2.34±0.79  -3.28 |  105.0 -2.20±0.83  +3.43 |  100.6 -0.96±0.52 +16.88 |   93.9 -0.72±0.75 +23.29
  DF T2 n=144 |  104.4 -0.24±0.24  -0.36 |  102.2 -0.24±0.24  -0.36 |   95.3 -0.05±0.05  -0.30 |   87.0 +0.08±0.08  -0.20 |   70.5 +1.04±0.61  +0.34 |   55.9 +1.57±1.02  +0.06
  DF E3 n=144 |  121.9 -0.09±0.08  -7.03 |  120.0 -0.62±0.77  -2.85 |  115.5 +0.79±0.60  +2.26 |  108.8 +2.75±1.22  +7.17 |  101.7 +6.91±2.78 +18.00 |   91.1 +7.45±3.13 +20.48
 DFF T2 n=144 |  104.8 +0.35±0.30  -0.01 |  102.5 +0.35±0.30  -0.01 |   95.6 +0.29±0.29  -0.01 |   87.2 +0.29±0.22  -0.01 |   70.4 +0.99±0.59  +0.28 |   55.9 +1.46±1.03  -0.00
 DFF E3 n=144 |  121.8 -0.08±0.08  -7.06 |  120.0 -1.01±0.82  -2.88 |  115.5 +1.96±1.38  +2.23 |  108.7 +1.19±1.28  +7.07 |  101.7 +1.04±1.34 +18.00 |   91.1 +1.14±2.30 +20.47
  ZG T2 n=144 |  106.5 -0.00±0.00  +1.74 |  103.6 -0.08±0.23  +1.12 |   95.4 +0.91±0.43  -0.21 |   82.5 +2.65±1.38  -4.77 |   64.1 +6.44±2.30  -6.02 |   52.0 +21.17±6.68  -3.91
  ZG E3 n=144 |  125.4 +0.03±0.03  -3.51 |  123.2 -1.07±0.77  +0.32 |  118.3 +1.42±0.47  +5.06 |  112.1 +5.02±1.20 +10.45 |  104.6 +13.25±3.54 +20.97 |   91.6 +12.84±3.54 +20.99
  DG T2 n= 48 |   65.2 +1.06±0.75 -28.23 |   68.1 +2.65±1.40 -20.21 |   64.4 +2.68±2.45 -15.39 |   59.3 +5.40±1.71 -10.22 |   54.6 +4.93±2.17  -4.00 |   42.6 +17.17±5.70 +12.09
 DFm E3 n=144 |  121.8 +0.64±0.94  -7.13 |  120.2 +0.83±0.78  -2.73 |  115.5 +1.04±0.43  +2.27 |  109.2 +3.06±1.13  +7.64 |  101.5 +6.69±2.67 +17.87 |   91.5 +8.28±3.22 +20.89
  ZS E3 n=144 |  135.6 -0.00±0.00  +6.67 |  130.9 +0.20±0.16  +7.98 |  121.7 +0.50±0.24  +8.50 |  110.3 +2.35±1.00  +8.68 |   92.7 +5.23±1.88  +9.00 |   82.3 +6.93±3.14 +11.63
```

1. **TA, tiered per-BS filling.** Every BS fills its slots with the current design's streams tier by tier (core
   privates, group-common, other privates), within a tier in decreasing order of its direct-channel energy to the
   stream's weakest decoder. If the slots suffice, every BS carries everything (joint); otherwise each BS keeps the
   streams it serves best. Clearly better than J at high R_th (more UEs get a private stream), clearly worse at low
   R_th in T2, where BS-specific core streams lose joint zero-forcing (as the thesis' 51–78 % remark predicts). But the
   group-common tier is almost never reached, so TA ≈ TA1: **R1 fails**. TC (TA without other UEs' privates) gave
   the same numbers (`results/TC__*__Q5.jsonl`; TC1 was stopped as uninformative).
2. **TB, per-BS clusters.** Every BS serves the clusters it is strongest to: a core user's private stream, then its
   group-common stream. Groups now displace private streams. This wins where the displaced privates cannot be zero-
   forced anyway (T2 setting b: +5 to +15 over TB1 at every R_th) and loses where they can (T2 c and e, −11 to −19 at
   low R_th): **R1 fails** at R_th ≤ 1.5. Lesson: a group-common stream must not take a slot from a private stream
   that zero-forcing can use.
3. **DF, DoF private layer.** $(n_p,s)$ maximize the zero-forcing DoF $\min(K,sN,\lfloor B(S^{\max}-1)/s\rfloor)$
   (ties: larger $s$); the private streams of the first $n_p$ SUS users go to their $s$ strongest BSs, the other UEs
   join the nearest of these core users, and the group layer fills the free slots. This is the first design with
   R1 margins (E3, R_th ≥ 1.5). Bench: 80.60 (vs sparse BF +7.54 ± 2.65, vs DF@SC +26.25 ± 1.76; shares 15.6 % /
   9.2 %). But in T2 it is the current design (s = B, n_p = min(2N, S_max − 1)): the group layer exists only in
   setting f, and **R1 fails in T2**. In E3 e it takes 10 streams on 2 BSs (critically loaded ZF) over 9 on 3, which
   costs about 16 against J. Q = 50 subset (`results/DF_Q50_subset.txt`): the same picture.
   * **DG** (DoF ties to fewer BSs per stream, T2 b and d only): it frees slots for groups, and so R1 holds, but
     single-BS private streams cost 20–28 against J at low R_th. Cooperation is worth too much here.
   * **DFm** (DF with a start that computes the RZF direction on the antennas of the carrying BSs only, the option
     `init='core_mask'` in `algo.py`): no measurable change (E3 109.5 vs 109.4). Dropped; SP uses the default start.
   * **DFF** (DF, then the other UEs' private streams in the slots still free, in the design and in its 1-layer):
     the 1-layer then refills the group slots with private streams, as J1 does. T2 becomes the current design and
     R1 ties everywhere: **fails**.
4. **ZG, greedy ZF user selection.** The core users are those of greedy ZF user selection (`hrs.zf_greedy`, stopped
   when one more user no longer raises the ZF sum rate; 4–8 instead of BN in many settings), then the DoF split, and
   groups around them. **R1 holds in both scenarios** (T2 +0.91 ± 0.43 at R_th 1, E3 +1.42 ± 0.47), but **R2 fails**.
   At the bench, 4–6 private streams and only 3–4 groups leave 3 slots per BS idle in several drops. On all 15 bench
   runs: 74.4 vs sparse BF 73.1 (+1.3 ± 3.2; `results/ZG__T2__bench.jsonl`).
5. **ZS, ZG with the cooperation size chosen by the ZF sum rate** (E3 only): fixes E3 e (9 streams on 3 BSs) and is
   the best E3 design (135.6 at R_th 0, +6.7 over J); R1 holds. Same bench problem as ZG (same T2 private layer).
6. **SZ / SP.** The private layer of ZS, but with the users in SUS order: the ZF-best prefix, which gives 6 private
   streams at the bench in 12 of 15 drops. SZ kept the thesis' groups (heads: all BN SUS users); in E3, K ≤ BN leaves
   no non-core UE and hence no groups, so it was not run. **SP** forms the groups around the private-layer users and
   is the final design.

## 4. Weak points

* **T2 settings b, c, d (S_max ≤ 2N − 1 with B = 2) and E3 c, e.** The private layer takes every slot (zero-forcing
  wants all S_max − 1 streams on all BSs), so there is no group layer: the design equals its 1-layer counterpart and,
  in these settings, the current joint design (identical runs). The R1 margins come from the other settings. This
  is structural: in these settings a group-common stream could only take a slot from a private stream that zero-
  forcing uses, which costs at low R_th (method TB above).
* **Fewer private streams than J where the fronthaul is generous** (T2 a: 4–6 instead of 6 + 2 non-core privates;
  the bench drop 802: 5 instead of 6). The ZF criterion at the direct channels with equal power is a conservative
  model of what the AO achieves. In T2 a the design is 7 below J at Q = 5, and at the bench 2.7 below J. At high R_th
  in T2 (1.5–2.5) it is 1.5–4.6 below J at Q = 50 (confirm tier), while it is 1 above J at R_th ≤ 1. Slots can stay unused when the groups are fewer than the free slots (T2 a, f;
  by design: no private stream beyond the ZF-best layer).
* **The criterion uses the direct channels only** (no DS-STAR cascade, which the deployment search optimizes later),
  equal power and the network power B P_max with no per-BS split: a model, not the AO.
* **Large networks.** The ZF criterion treats every BS–UE link as interfering, which is right for the channel model
  of the thesis: with α = 2.2, β0 = 10⁻³ and σ² = 10⁻¹² W, a BS 1 km away is still about 30 dB above the noise at a
  UE. With few slots per BS, the private layer then takes most of the fronthaul (16 BSs, S_max = 8: 15 private
  streams on 7 BSs each), and little is left for groups. A BS with many free slots also carries groups that are not
  near it (it ranks all groups by its channel to the weakest member). This is harmless to the problem (the AO can
  give them zero power), but it is not "local".
* **Interaction with the deployment search at Q = 50.** The rate score of search R averages over the served streams,
  so SP's group-common streams also steer the search, including at low R_th, where they do not matter for F. Together
  with the path noise of two different evaluators, this gives small negative means at R_th ≤ 0.5 at Q = 50 (E3 f:
  −5.0 ± 3.4 at R_th = 0; E3 a: about −1 at every R_th, where the group layer is one small group). At Q = 5 the same
  pairs are exact ties. The framework is fixed (ground rule 8), and the same holds for the current design against J1.
* **Not nested in S_max.** When the private layer grows with S_max, the groups are formed anew around it, so a larger
  S_max does not always carry a superset of streams (§1.6). The measured rate still rises (§2.3), with two steps of
  −0.3 ± 0.5.
* **R2 margin.** +5.7 ± 2.6 over sparse beamforming (J: +8.3), positive at every M (+8.9 / +3.7 / +4.4).
* **R3.** The group-common share at the bench is 9.6 % ("around 10 %"), the super-common share 16.3 %.
* **CPU.** No cost issue: the design takes 2–3 ms, and SP's runs cost about the same as J's (fewer streams than J in
  T2 a and f, more in E3 b and d).

## 5. Files

* `code/layered.py`: the design (`design`, `design_1layer`).
* `code/perbs.py`: the other methods tried.
* `code/algo.py`: one added option, `init_core_mask` (`INITS['core_mask']`), off by default and used by DFm only;
  JH and J1 unchanged (`harness/selftest.py` passed again with the final code).
* `harness/variants.py`: SP, SP1 and the variants tried.
* `tools/large_network.py`: the association of a design with B BSs scattered over a large area (no AO);
  `tools/summary_tried.py`: the summary table above; `tools/r1_pair.py`: R1 on the runs common to SP and SP1 only
  (the 6-drop extension).
* `results/`: every run (harness format, unedited), the `tools/report.py` outputs (`*report*.json`, `report_*.txt`),
  `large_network_*.txt`, `summary_tried_Q5.txt`, `DF_Q50_subset.txt`, `report_SP_R1_Q50_*.txt`. The Q = 50 grid files
  of SP and SP1 (`results/SP__T2.jsonl` etc.) hold drops 901–906: 901–903 are the confirm tier, 904–906 the extension.
* `package.diff`: unified diff (`diff -ruN`) of `code/`, `harness/` and `tools/` against the package v10.
* Notes. The Q = 5 runs of SP (`results/SP__*__Q5.jsonl`) were made with the same function under its exploratory name
  (`perbs.sus_zf_regroup`). `layered.design` returns the identical design on every configuration of the harness (333
  checked: grid, bench, smax), and three of these runs re-run with `layered.design` gave identical traces.
  `SP1__*__Q5` and `SZ1__*__Q5` are identical runs (288 of 288). ZG@SC (bench) and the Q = 50 runs of ZG were
  stopped when R2 ruled ZG out (partial files).
