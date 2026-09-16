"""
Assorted facts quoted in the report that do not belong to any single step.

Reads:  work/pairs_tiered.parquet, work/node_ann.parquet.
Writes: work/s20_tier1_pairs.csv (all 85 Tier 1 significant pairs).
Feeds:  report sections 1, 3, 5 and 10.

Covers: the two pairs at SDC <= -0.6 that the signed cutoff excluded from the network; the
filtering cascade with pair and event counts at each stage; the full Tier 1 list, which shows that
61 of the 85 are co-quantifiable and is the direct evidence that a proximity filter would discard
documented crosstalk; per-modification annotation coverage; and the share of site-level annotation
that comes only from large-scale MS, which is the quantitative basis for the circularity argument.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np
con = duckdb.connect(); con.execute("PRAGMA threads=12;")
pt = (WORK/"pairs_tiered.parquet").as_posix()
log("=== MISCELLANEOUS FACTS FOR THE REPORT ===")
d = con.execute(f"""SELECT nodeA,nodeB,geneA,geneB,Score,qvalue,tier,in_network
   FROM '{pt}' WHERE Score<=-0.6""").df()
log(f"  pairs with SDC <= -0.6 (excluded by the signed >=0.6 cutoff): {len(d)}")
for r in d.itertuples(): log(f"    {r.nodeA} ({r.geneA}) x {r.nodeB} ({r.geneB})  SDC={r.Score:.4f} q={r.qvalue:.3g} tier={r.tier} in_network={r.in_network}")

log("  --- filtering cascade ---")
for lab,w in [("all pairs in the raw list","TRUE"),
              ("q < 0.05","is_sig"),
              ("q < 0.05 and not co-peptide","is_sig AND NOT same_peptide"),
              ("SDC >= 0.6 (published network)","in_network"),
              ("SDC >= 0.6 and not co-peptide","in_network AND NOT same_peptide")]:
    n = con.execute(f"SELECT count(*) FROM '{pt}' WHERE {w}").fetchone()[0]
    nn = con.execute(f"SELECT count(DISTINCT x) FROM (SELECT nodeA x FROM '{pt}' WHERE {w} UNION ALL SELECT nodeB FROM '{pt}' WHERE {w})").fetchone()[0]
    log(f"    {lab:<34} pairs={n:>12,}  PTM events={nn:>7,}")

log("  --- tier 1 (site-level explicit crosstalk) significant pairs, all 85 ---")
t1 = con.execute(f"""SELECT nodeA,nodeB,geneA,geneB,Score,stratum,in_network,same_peptide,position_gap
   FROM '{pt}' WHERE is_sig AND tier=1 ORDER BY abs(Score) DESC""").df()
t1.to_csv(WORK/"s20_tier1_pairs.csv", index=False)
for r in t1.head(20).itertuples():
    log(f"    {r.geneA:<10}{r.nodeA.split('|',1)[1]:<12} x {r.nodeB.split('|',1)[1]:<12} SDC={r.Score:+.3f} {r.stratum} net={r.in_network} copep={r.same_peptide} gap={r.position_gap}")
log(f"    ... {len(t1)} total; in network {int(t1.in_network.sum())}; co-peptide {int(t1.same_peptide.sum())}")

log("  --- per-modification annotation coverage (denominator for the discussion) ---")
A = pd.read_parquet(WORK/"node_ann.parquet")
for u,g in A.groupby("unimod"):
    log(f"    {UNIMOD[u][0]:<14} events={len(g):>5}  annotated_same_mod={100*g.ann_samemod.mean():>5.1f}%  "
        f"manual_experimental={100*g.ann_exp.mean():>5.1f}%  large_scale_MS={100*g.ann_htp.mean():>5.1f}%")
log(f"    {'ALL':<14} events={len(A):>5}  annotated_same_mod={100*A.ann_samemod.mean():>5.1f}%  "
    f"manual_experimental={100*A.ann_exp.mean():>5.1f}%  large_scale_MS={100*A.ann_htp.mean():>5.1f}%")

log("  --- circularity: share of site-level annotation that is large-scale MS ---")
both = A[A.ann_mod]
log(f"    events on an annotated modified residue: {len(both):,}")
log(f"      with manual experimental evidence (ECO:0000269):       {int(both.ann_exp.sum()):,} ({100*both.ann_exp.mean():.1f}%)")
log(f"      with large-scale MS evidence (ECO:0007744/0000244):    {int(both.ann_htp.sum()):,} ({100*both.ann_htp.mean():.1f}%)")
log(f"      large-scale MS only, no manual experimental evidence:  {int((both.ann_htp & ~both.ann_exp).sum()):,} ({100*(both.ann_htp & ~both.ann_exp).mean():.1f}%)")
t5 = con.execute(f"""SELECT count(*) n, count(*) FILTER ((NOT expA) AND (NOT expB)) both_htp_only
   FROM '{pt}' WHERE is_sig AND tier=5""").fetchone()
log(f"    tier-5 significant pairs: {t5[0]:,}; both sites supported only by large-scale MS: {t5[1]:,} ({100*t5[1]/t5[0]:.1f}%)")
