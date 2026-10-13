"""Discrete-design variants for the harness: name -> entry.

Entry keys:
  design      'module:function' in code/ ; the function design(net, p) -> hrs.Design (see code/discrete.py)
  one_layer   the name of the variant that is this design's 1-layer RSMA counterpart (requirement R1 of BRIEF.md)
  opt         (optional) dict of options passed to discrete.run (keys of newprop.NEW); leave it out unless your design
              needs a changed evaluation, and then explain why in DESIGN.md
Use a new name instead of editing a variant that already has results.
KEY@SC (no entry needed; harness/grid.py and tools/report.py) is the stream-clustering benchmark (code/benchmarks.py)
on the streams of KEY's design.

The baselines (the current thesis):
  JH   current design (joint transmission, HRS order)           discrete.current
  J1   its 1-layer RSMA (group-common streams removed)          discrete.current_1layer
  J    the thesis framework = JH or J1 by the stream-set test at the initial design of the search; derived by
       tools/report.py from the runs of JH and J1 (not run separately)
Design round v10 (DESIGN.md): SP / SP1 the final design and its 1-layer counterpart (code/layered.py); the other
entries are the methods tried on the way (code/perbs.py).
"""

VARIANTS = {
    'JH': dict(design='discrete:current', one_layer='J1'),
    'J1': dict(design='discrete:current_1layer'),
    # ---- design round v10 -------------------------------------------------------------------------------------
    # the final design (code/layered.py; DESIGN.md) and its 1-layer RSMA counterpart
    'SP': dict(design='layered:design', one_layer='SP1'),
    'SP1': dict(design='layered:design_1layer'),
    # the other methods tried (code/perbs.py; DESIGN.md, "Other methods tried"), in the order they were tried
    'TA': dict(design='perbs:tiered', one_layer='TA1'),               # per-BS filling by affinity, tier by tier
    'TA1': dict(design='perbs:tiered_1layer'),
    'TB': dict(design='perbs:cluster', one_layer='TB1'),              # per-BS clusters: core private + its group
    'TB1': dict(design='perbs:cluster_1layer'),
    'TC': dict(design='perbs:coreonly', one_layer='TC1'),             # TA without non-core privates (TC1 not run)
    'TC1': dict(design='perbs:coreonly_1layer'),
    'DF': dict(design='perbs:dof', one_layer='DF1'),                  # DoF private layer, groups in the free slots
    'DF1': dict(design='perbs:dof_1layer'),
    'DG': dict(design='perbs:dof_small', one_layer='DG1'),            # DF, DoF ties to fewer BSs (T2 b, d only)
    'DG1': dict(design='perbs:dof_small_1layer'),
    'DFm': dict(design='perbs:dof', one_layer='DF1m', opt=dict(init='core_mask')),   # DF, masked AO start (E3)
    'DF1m': dict(design='perbs:dof_1layer', opt=dict(init='core_mask')),
    'DFF': dict(design='perbs:dof_fill', one_layer='DFF1'),           # DF, then non-core privates in free slots
    'DFF1': dict(design='perbs:dof_fill_1layer'),
    'ZG': dict(design='perbs:dof_zf', one_layer='ZG1'),               # greedy ZF users, DoF split, groups on them
    'ZG1': dict(design='perbs:dof_zf_1layer'),
    'ZS': dict(design='perbs:zfs', one_layer='ZS1'),                  # ZG with s by the ZF rate (E3 only)
    'ZS1': dict(design='perbs:zfs_1layer'),
    'SZ': dict(design='perbs:sus_zf', one_layer='SZ1'),               # SP private layer, thesis groups (not run)
    'SZ1': dict(design='perbs:sus_zf_1layer'),                        # = SP1 (identical designs)
}
