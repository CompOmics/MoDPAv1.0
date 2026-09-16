"""
Step 6, part 2: monotonicity of support across the |SDC| strata.

Reads:  work/s15_background.csv.
Writes: work/s24_monotonicity.csv.
Feeds:  report section 8, monotonicity table.

Split out of s15 so the trend statistics can be recomputed without re-running the 8-minute
permutation.

Why this is not scipy's default spearmanr p-value: that p comes from a t-approximation whose
denominator vanishes as |rho| approaches 1. At rho = +1.000 across ten strata it returned
p = 6.6e-64, and at rho = -1.000 it returned p = 0. Both are artefacts of the approximation, not
results. With ten points the smallest attainable two-sided p is 2/10! = 5.5e-7, so any figure below
that is impossible by construction.

This script uses an exact permutation test over all orderings when n! is small enough to enumerate
(n <= 8, that is 40,320 orderings), and a Monte Carlo permutation test with 200,000 shuffles
otherwise. Both are two-sided on |rho|. The Monte Carlo floor is 1/200,001, about 5.0e-6, and is
reported as such rather than as an exact value.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import pandas as pd, numpy as np, itertools, math

rng = np.random.default_rng(42)
R = pd.read_csv(WORK/"s15_background.csv")
order = {s:i for i,s in enumerate(sorted(R.stratum.unique()))}

def perm_p(x, y, n_mc=200_000):
    """Two-sided permutation p-value for Spearman rho, exact where feasible.

    Permuting one vector against the other is the correct null here: it asks how often a random
    pairing of the same stratum values reaches a rank correlation at least as extreme as observed.

    Spearman rho is Pearson correlation on ranks, and permuting y only reorders its ranks, so the
    means and standard deviations of both rank vectors are invariant across permutations. Only the
    cross term changes, which makes rho a linear function of dot(rank_x, permuted rank_y). The whole
    null distribution is therefore one matrix product, rather than 200,000 calls to spearmanr.

    Returns (rho, p, mode, floor).
    """
    n = len(x)
    # .copy() because pandas can hand back a read-only view, which the in-place centring below
    # cannot write to.
    rx = pd.Series(x).rank().to_numpy(float).copy()
    ry = pd.Series(y).rank().to_numpy(float).copy()
    rx -= rx.mean(); ry -= ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    rho = float(rx @ ry / denom)
    obs = abs(rx @ ry) - 1e-9
    if math.factorial(n) <= 40_320:
        P = np.array(list(itertools.permutations(ry)))
        mode, floor = f"exact ({math.factorial(n):,} orderings)", 1/math.factorial(n)
        return rho, float((np.abs(P @ rx) >= obs).mean()), mode, floor
    P = np.tile(ry, (n_mc, 1))
    P = rng.permuted(P, axis=1)
    hits = int((np.abs(P @ rx) >= obs).sum())
    return rho, (1+hits)/(1+n_mc), f"Monte Carlo ({n_mc:,} shuffles)", 1/(1+n_mc)

log("=== STEP 6b: MONOTONICITY, permutation test ===")
rows = []
for sg in ["pos","neg"]:
    sub = R[R.sign==sg].copy()
    sub["k"] = sub.stratum.map(order)
    sub = sub.sort_values("k")
    if len(sub) < 4: continue
    for col, lab in [("exact_T1_4_full_stratum","exact T1-4 rate, full stratum"),
                     ("obs_T1_4","observed T1-4 rate, sample"),
                     ("ratio_T1_4","ratio to matched null"),
                     ("obs_T1_2","observed T1-2 rate"),
                     ("obs_T1_5","observed T1-5 rate, type (a)")]:
        if col not in sub.columns or sub[col].isna().all(): continue
        if sub[col].nunique() < 2:
            log(f"  {sg}  {lab:<32} constant, no rank correlation defined"); continue
        rho, p, mode, floor = perm_p(sub.k.to_numpy(), sub[col].to_numpy())
        rows.append(dict(sign=sg, metric=col, label=lab, n_strata=len(sub),
                         rho=rho, p=p, mode=mode, p_floor=floor))
        at = " (at the floor)" if p <= floor*1.5 else ""
        log(f"  {sg}  {lab:<32} rho={rho:+.3f}  p={p:.3g}{at}  [{mode}]")
M = pd.DataFrame(rows)
M.to_csv(WORK/"s24_monotonicity.csv", index=False)
log(f"  -> work/s24_monotonicity.csv ({len(M)} tests)")
