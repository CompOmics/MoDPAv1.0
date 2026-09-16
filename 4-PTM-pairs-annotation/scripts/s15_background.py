"""
Step 6: matched permutation null and the monotonicity test.

Reads:  work/pairs_tiered.parquet, work/node_ann.parquet, the pair_t*.parquet tables,
        work/ev_tryptic.parquet.
Writes: work/s15_background.csv (one row per stratum and sign, with observed and null rates,
        the ratio, an empirical p-value, and a sensitivity value).
Feeds:  report section 8.

An enrichment claim needs a matched null, because well-studied proteins accumulate annotations and
detections together. The null here is node-level: each PTM event is binned on the three properties
that could produce support on their own, and each endpoint of an observed pair is replaced by a
random node from its own bin. That holds those properties fixed while destroying the association
between the specific pair and its score.

The observed sample is drawn with a REPEATABLE seed, so the whole script is deterministic. An
earlier version used an unseeded sample at MMAX=25,000, and the strata below 0.40 moved by more
than their own effect size between runs. exact_T1_4_full_stratum records the support rate over the
complete stratum with no sampling at all, and is the value to quote.

Matching bins are (modification type, detection-frequency decile, annotation-density decile),
coarsened to (modification, detection decile) and then to (modification) when a bin holds fewer
than 5 nodes, so that rare modifications still have somewhere to draw from. The full three-way bin
covers 7,173 of 7,644 events.

Performance notes, which are the reason the code looks the way it does:
  * pack() and packp() encode an unordered pair of node or protein indices as a single int64, so
    membership of the Tier 1, Tier 2, Tier 4 and co-peptide sets becomes one sorted-array lookup
    through isin_sorted rather than a Python set lookup per pair.
  * bflat, bstart and bsize flatten the per-node bins into one array with offsets, so an entire
    permutation is drawn with two vectorised index computations. The earlier list-comprehension
    version took roughly 25 s per permutation at this sample size, which is 7 hours for the full
    run; this version takes about 5 s per stratum for all 1,000 permutations.
  * tier_of() reproduces the SQL tier assignment from s13. The agreement check at the top of the
    run is there to catch any drift between the two implementations and reports 100.0000% on a
    200,000-pair sample.

duckdb applies USING SAMPLE before WHERE, so the observed sample is drawn from a subquery that has
already been filtered. Sampling the table directly and then filtering returns a few hundred rows
instead of 25,000.

P-values are empirical permutation p-values, not t-tests. p_T1_4 is one-tailed on the upper tail
(enrichment), p_T1_4_lower one-tailed on the lower tail (depletion), and p_T1_4_two_sided the
doubled smaller tail. The monotonicity Spearman tests at the end are scipy defaults, two-sided.

No pair is excluded from the observed set. The two final columns record the co-quantifiable rate in
the tested sample and what the observed rate would become if those pairs were dropped, so the size
of that effect is visible without a filter being imposed.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np, json, time

rng = np.random.default_rng(42)
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()

A = pd.read_parquet(WORK/"node_ann.parquet").reset_index(drop=True)
N = len(A); idx = {n:i for i,n in enumerate(A.name)}
acc_codes, acc_uniq = pd.factorize(A.acc)
ann  = A.ann_mod.to_numpy(); expv = A.ann_exp.to_numpy(); htpv = A.ann_htp.to_numpy()
mod  = A.unimod.to_numpy(); nann = A.n_ann_on_protein.to_numpy(); nruns = A.n_runs.to_numpy()
P = len(acc_uniq)

def pack(i, j):
    """Encode an unordered pair of node indices as one int64: min*N + max.

    N is the number of PTM events (7,644), so min*N + max is unique per unordered pair and fits
    comfortably in int64. Sorting the two indices first makes the encoding order-independent, which
    matters because the evidence tables and the association list do not agree on which endpoint is
    written first. Packing turns "is this pair in that set" into a single sorted-array lookup, which
    is what makes 25 million membership tests per stratum affordable.
    """
    lo = np.minimum(i, j); hi = np.maximum(i, j)
    return lo.astype(np.int64) * N + hi

def packp(i, j):
    """Same encoding for a pair of protein indices, used for the Tier 4 protein-pair set."""
    lo = np.minimum(i, j); hi = np.maximum(i, j)
    return lo.astype(np.int64) * P + hi

def load_pairs(df, a="nodeA", b="nodeB"):
    """Map a two-column frame of PTM event names to the integer node indices used everywhere here."""
    return (df[a].map(idx).to_numpy(), df[b].map(idx).to_numpy())

t1 = pd.read_parquet(WORK/"pair_t1a.parquet")[["nodeA","nodeB"]]
S_T1 = np.unique(pack(*load_pairs(t1)))
S_T2 = np.unique(pack(*load_pairs(pd.read_parquet(WORK/"pair_t2.parquet"))))
t4 = pd.read_parquet(WORK/"pair_t4.parquet")
acc2i = {a:i for i,a in enumerate(acc_uniq)}
S_T4 = np.unique(packp(t4.accA.map(acc2i).to_numpy(), t4.accB.map(acc2i).to_numpy()))
tr = pd.read_parquet(WORK/"ev_tryptic.parquet")
tr = tr[tr.same_peptide_std]
S_CP = np.unique(pack(*load_pairs(tr)))

def isin_sorted(v, s):
    """Vectorised membership of v in the sorted array s, faster here than np.isin.

    searchsorted gives the insertion point of each value; the value is present exactly when the
    element sitting at that index equals it. Insertion points past the end are clamped to 0 purely
    so the subsequent fancy-index is in range, and that row then fails the equality test anyway.
    """
    if len(s)==0: return np.zeros(len(v), bool)
    p = np.searchsorted(s, v); p[p>=len(s)] = 0
    return s[p] == v

def tier_of(i, j):
    """Vectorised copy of the tier assignment in s13, for arrays of node indices.

    The assignments run from weakest to strongest so that a later, stronger tier overwrites an
    earlier one. That reproduces the first-match-wins semantics of the SQL CASE expression in s13
    while staying a handful of whole-array operations. The two implementations are checked against
    one another on a 200,000-pair sample at the top of the run.
    """
    k = pack(i, j); kp = packp(acc_codes[i], acc_codes[j])
    same_prot = acc_codes[i] == acc_codes[j]
    t = np.full(len(i), 9, np.int8)                                      # 9 = unsupported
    t5 = ann[i] & ann[j] & (expv[i]|htpv[i]) & (expv[j]|htpv[j]);        t[t5] = 5
    t4m = (~same_prot) & isin_sorted(kp, S_T4);                           t[t4m] = 4
    t3 = same_prot & ann[i] & ann[j];                                     t[t3] = 3
    t2 = isin_sorted(k, S_T2);                                            t[t2] = 2
    t1 = isin_sorted(k, S_T1);                                            t[t1] = 1
    return t

def is_cop(i, j):
    """Whether each pair can fall on one tryptic peptide. Reported, never used to filter."""
    return isin_sorted(pack(i, j), S_CP)

log("=== STEP 6: MATCHED PERMUTATION BACKGROUND ===")
# validation of the vectorised tier function against the SQL assignment
chk = con.execute(f"SELECT * FROM (SELECT nodeA,nodeB,tier,same_peptide FROM '{pt}' WHERE is_sig) USING SAMPLE 200000 ROWS").df()
ci, cj = load_pairs(chk)
log(f"  vectorised tier agrees with the SQL assignment on {100*(tier_of(ci,cj)==chk.tier.to_numpy()).mean():.4f}% of 200k sampled pairs")
log(f"  co-peptide flag agreement: {100*(is_cop(ci,cj)==chk.same_peptide.to_numpy()).mean():.4f}%")

# ---- node bins: modification type x detection-frequency decile x annotation-density decile
#
# These are the three properties that could produce annotation support on their own, independently
# of the SDC. Modification type, because phosphorylation is far better annotated than the rest.
# Detection frequency, because a site seen in more runs has more chance of being studied. And
# annotation density of the protein, because well-studied proteins accumulate both annotations and
# detections. Replacing an endpoint by another node from the same bin holds all three fixed.
def dec(x):
    """Decile index 0-9 by rank, so bins are equally populated regardless of the distribution shape."""
    r = pd.Series(x).rank(method="average", pct=True).to_numpy()
    return np.clip((r*10).astype(int), 0, 9)
d_runs, d_ann = dec(nruns), dec(nann)
# Three nested bin definitions, from most to least specific. The coarser two are fallbacks for rare
# modifications: succinylation has three events in the matrix, so its full three-way bins would hold
# one or two nodes and the "random" draw would keep returning the node itself.
keys_full = np.array([f"{m}|{a}|{b}" for m,a,b in zip(mod, d_runs, d_ann)])
keys_mid  = np.array([f"{m}|{a}"     for m,a   in zip(mod, d_runs)])
keys_low  = np.array([f"{m}"         for m     in mod])
members = {}
for lvl,kk in (("full",keys_full),("mid",keys_mid),("low",keys_low)):
    dd = {}
    for i,k in enumerate(kk): dd.setdefault(k, []).append(i)
    members[lvl] = {k:np.array(v) for k,v in dd.items()}
# Per node, take the most specific bin that holds at least 5 candidates. The for/else fires only if
# even the modification-only bin is too small, in which case there is nothing coarser to fall back to.
bin_of = np.empty(N, object); lvl_used = np.empty(N, object)
for i in range(N):
    for lvl,kk in (("full",keys_full),("mid",keys_mid),("low",keys_low)):
        if len(members[lvl][kk[i]]) >= 5: bin_of[i]=members[lvl][kk[i]]; lvl_used[i]=lvl; break
    else: bin_of[i]=members["low"][keys_low[i]]; lvl_used[i]="low"
log(f"  matching-bin level used per node: {pd.Series(lvl_used).value_counts().to_dict()}")
log(f"  median bin size: {np.median([len(b) for b in bin_of]):.0f}; min: {min(len(b) for b in bin_of)}")
# Flatten the ragged per-node bins into one array plus offsets, so a whole permutation can be drawn
# with array arithmetic instead of a Python loop over pairs. bflat holds every bin's members end to
# end; bstart[i] is where node i's bin begins and bsize[i] how long it is, so
# bflat[bstart[i] + floor(u * bsize[i])] with u uniform in [0,1) is a uniform draw from that bin.
bsize = np.array([len(b) for b in bin_of], np.int64)
bstart = np.concatenate([[0], np.cumsum(bsize)[:-1]])
bflat = np.concatenate([np.asarray(b, np.int64) for b in bin_of])

# MMAX is the observed-sample size per cell. At the ~0.5% support rate of the low strata, 25,000
# pairs rests on about 125 supported pairs, giving roughly +-11% Poisson noise, which is the same
# size as the effect being tested there. 200,000 cuts that noise by a factor of about 2.8.
NPERM, MMAX = 1000, 200000
rows=[]
strata = con.execute(f"""SELECT stratum, sign, count(*) n,
   count(*) FILTER (same_peptide) n_cop FROM '{pt}' WHERE is_sig GROUP BY 1,2 ORDER BY 1,2""").df()
for _,st in strata.iterrows():
    if st.n < 20: continue
    # Exact T1-4 rate over the WHOLE stratum, no sampling. The permutation below runs on a sample,
    # so this is the number to quote; the sampled observed rate is only there to be comparable with
    # the null, which is computed on the same sample.
    exact = con.execute(f"""SELECT count(*) FILTER (tier<=4)::DOUBLE/count(*) FROM '{pt}'
        WHERE is_sig AND stratum='{st.stratum}' AND sign='{st['sign']}'""").fetchone()[0]
    obs = con.execute(f"""SELECT * FROM (SELECT nodeA,nodeB FROM '{pt}' WHERE is_sig
        AND stratum='{st.stratum}' AND sign='{st['sign']}')
        USING SAMPLE reservoir({MMAX} ROWS) REPEATABLE (42)""").df()
    oi, oj = load_pairs(obs); oi=oi.astype(np.int64); oj=oj.astype(np.int64); m = len(oi)
    ot = tier_of(oi, oj)
    o14 = float((ot<=4).mean()); o12 = float((ot<=2).mean()); o15 = float((ot<=5).mean())
    n14 = np.empty(NPERM); n12 = np.empty(NPERM); n15 = np.empty(NPERM)
    # One permutation = replace both endpoints of every observed pair by a random node from the
    # same bin, then recompute the tier. The offsets and sizes are looked up once outside the loop.
    si0, sz0 = bstart[oi], bsize[oi]; si1, sz1 = bstart[oj], bsize[oj]
    for p in range(NPERM):
        si = bflat[si0 + (rng.random(m)*sz0).astype(np.int64)]
        sj = bflat[si1 + (rng.random(m)*sz1).astype(np.int64)]
        # Drop self pairs only. Co-quantifiable resamples are kept, at their natural rate of about
        # 0.03%, so the null is not filtered in a way the observed set is not.
        keep = (si!=sj)
        if keep.sum()==0:
            n14[p]=n12[p]=n15[p]=np.nan; continue
        tt = tier_of(si[keep].astype(np.int64), sj[keep].astype(np.int64))
        n14[p]=(tt<=4).mean(); n12[p]=(tt<=2).mean(); n15[p]=(tt<=5).mean()
    def pv(o,nul):
        """One-tailed empirical permutation p-value, UPPER tail: is the observed rate higher
        than the matched null? Add-one correction on numerator and denominator, so the floor is
        1/(NPERM+1) and a p-value of exactly 0 cannot be reported. No distributional assumption
        is made; the null standard deviation is stored for display only."""
        nul=nul[~np.isnan(nul)]
        return (1+np.sum(nul>=o))/(1+len(nul))
    def pv_lo(o,nul):
        """Same, LOWER tail: is the observed rate lower than the matched null? Needed because the
        upper-tail test alone cannot distinguish 'no enrichment' from 'significant depletion',
        which is the situation in every negative-SDC stratum."""
        nul=nul[~np.isnan(nul)]
        return (1+np.sum(nul<=o))/(1+len(nul))
    def pv_two(o,nul):
        """Two-sided version, the doubled smaller tail capped at 1. Reported alongside the
        one-tailed values so the direction of each test is explicit."""
        return min(1.0, 2*min(pv(o,nul), pv_lo(o,nul)))
    ocop = is_cop(oi, oj)
    keep_ns = ~ocop
    o14_ns = float((ot[keep_ns]<=4).mean()) if keep_ns.sum() else float('nan')
    rows.append(dict(stratum=st.stratum, sign=st["sign"], n_stratum=int(st.n), n_copeptide=int(st.n_cop),
        m_tested=m, exact_T1_4_full_stratum=exact, frac_copeptide_in_sample=float(ocop.mean()),
        obs_T1_4_sensitivity_copeptide_dropped=o14_ns,
        obs_T1_4=o14, null_T1_4=np.nanmean(n14), null_sd_T1_4=np.nanstd(n14), p_T1_4=pv(o14,n14),
        p_T1_4_lower=pv_lo(o14,n14), p_T1_4_two_sided=pv_two(o14,n14),
        ratio_T1_4=o14/max(np.nanmean(n14),1e-12),
        obs_T1_2=o12, null_T1_2=np.nanmean(n12), p_T1_2=pv(o12,n12),
        p_T1_2_lower=pv_lo(o12,n12), p_T1_2_two_sided=pv_two(o12,n12),
        obs_T1_5=o15, null_T1_5=np.nanmean(n15), p_T1_5=pv(o15,n15)))
    log(f"  {st.stratum:<15}{st['sign']}  m={m:>7,}  exact={exact:.5f} obs={o14:.5f} null={np.nanmean(n14):.5f}+-{np.nanstd(n14):.5f} "
        f"ratio={o14/max(np.nanmean(n14),1e-12):>6.2f} p={pv(o14,n14):.4f} | T1-2 obs={o12:.6f} null={np.nanmean(n12):.6f} p={pv(o12,n12):.4f}"
        f" | p_lo={pv_lo(o14,n14):.4f} p_2s={pv_two(o14,n14):.4f}"
        f" | copep={100*ocop.mean():5.2f}% sens_T1-4={o14_ns:.5f}")
R = pd.DataFrame(rows); R.to_csv(WORK/"s15_background.csv", index=False)

import scipy.stats as sst
order = {s:i for i,s in enumerate(sorted(R.stratum.unique()))}
log("=== MONOTONICITY (Spearman of ratio vs stratum rank) ===")
for sg in ["pos","neg"]:
    sub = R[R.sign==sg].copy(); sub["k"]=sub.stratum.map(order)
    if len(sub)>=4:
        for col in ["obs_T1_4","ratio_T1_4","obs_T1_2","obs_T1_5"]:
            r = sst.spearmanr(sub.k, sub[col])
            log(f"  {sg}  {col:<12} rho={r.statistic:+.3f}  p={r.pvalue:.4g}  (n strata={len(sub)})")
