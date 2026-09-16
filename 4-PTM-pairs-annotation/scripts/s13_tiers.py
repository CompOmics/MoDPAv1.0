"""
Step 4, part 3: assign an evidence tier to every pair in the raw list.

Reads:  work/raw_pairs.parquet, work/node_ann.parquet, the four work/pair_t*.parquet tables,
        work/ev_tryptic.parquet, CLUSTERS, FILTERED.
Writes: work/pairs_tiered.parquet (all 29,211,546 pairs, about 1.4 GB) - the central table that
        every later step queries.
Feeds:  report sections 6, 8 and 9, and both delivered TSVs.

Tiers are assigned from annotation alone, before any score is consulted, so the classification
cannot be influenced by the SDC. The CASE expression is ordered from strongest to weakest evidence
and duckdb takes the first branch that matches, so each pair receives the best tier it qualifies
for:
    1  site-level explicit crosstalk   (in pair_t1a)
    2  same named enzyme on both sites (in pair_t2)
    3  same protein, both annotated
    4  different proteins, interaction or shared complex (in pair_t4)
    5  both sites annotated with an experimental code, no documented link
    9  unsupported

Tier 5 is type (a) evidence only and is reported separately throughout; it is never presented as
support for an association.

Both non-significant and significant pairs are tiered, so section 6.1 of the report can compare the
two and show that significance alone carries no annotation signal.

Cluster labels come from CLUSTERS, the weighted Leiden partition, joined on the PTM event
identifier. The withdrawn unweighted partition is not used. same_cluster is null-safe: a pair is
only "same cluster" when both endpoints actually carry a label.

The co-quantification flag from s10 is joined in as same_peptide but is never used as a filter here
or anywhere downstream.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, time

con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pq = (WORK/"raw_pairs.parquet").as_posix()

A   = pd.read_parquet(WORK/"node_ann.parquet")
T1A = pd.read_parquet(WORK/"pair_t1a.parquet")[["nodeA","nodeB"]]
T1B = pd.read_parquet(WORK/"pair_t1b.parquet")[["nodeA","nodeB"]]
T2  = pd.read_parquet(WORK/"pair_t2.parquet")
T4  = pd.read_parquet(WORK/"pair_t4.parquet")
TRY = pd.read_parquet(WORK/"ev_tryptic.parquet")[["nodeA","nodeB","position_gap","same_peptide_std","same_peptide_modblocked"]]
cl  = pl.read_csv(CLUSTERS).select(["name","LeidenCluster"]).to_pandas()
net = pl.read_csv(FILTERED, infer_schema_length=10000).select(["nodeA","nodeB"]).to_pandas()
net["in_network"] = True
for d in (T1A,T1B,T2,T4,TRY,net):
    pass
for n,d in [("A",A),("T1A",T1A),("T1B",T1B),("T2",T2),("T4",T4),("TRY",TRY),("CL",cl),("NET",net)]:
    con.register(n, d)

out = (WORK/"pairs_tiered.parquet").as_posix()
BAND = """CASE WHEN absScore>=0.70 THEN 'S10_0.70+' WHEN absScore>=0.60 THEN 'S09_0.60-0.70'
 WHEN absScore>=0.55 THEN 'S08_0.55-0.60' WHEN absScore>=0.50 THEN 'S07_0.50-0.55'
 WHEN absScore>=0.45 THEN 'S06_0.45-0.50' WHEN absScore>=0.40 THEN 'S05_0.40-0.45'
 WHEN absScore>=0.35 THEN 'S04_0.35-0.40' WHEN absScore>=0.30 THEN 'S03_0.30-0.35'
 WHEN absScore>=0.25 THEN 'S02_0.25-0.30' ELSE 'S01_<0.25' END"""

t0=time.time(); log("=== STEP 4c: TIER ASSIGNMENT over all pairs ===")
con.execute(f"""
COPY (
-- base: the raw list plus the three things derived from the score alone. is_sig is carried as a
-- column rather than applied as a filter, so the non-significant pairs stay available as the
-- comparison group for report section 6.1.
WITH base AS (
  SELECT nodeA, nodeB, Score, pval, qvalue, PCC, abs(Score) AS absScore,
         CASE WHEN Score>=0 THEN 'pos' ELSE 'neg' END AS sign,
         qvalue < 0.05 AS is_sig
  FROM '{pq}'
-- j: attach both endpoints' annotation (inner joins, since every node is in node_ann), then every
-- pair-level evidence table as a LEFT JOIN so a pair with no evidence survives with nulls.
), j AS (
  SELECT b.*, {BAND} AS stratum,
    a1.acc accA, a1.pos posA, a1.res resA, a1.unimod modA, a1.gene geneA, a1.n_runs nrunsA,
    a2.acc accB, a2.pos posB, a2.res resB, a2.unimod modB, a2.gene geneB, a2.n_runs nrunsB,
    a1.ann_mod annA, a2.ann_mod annB,
    a1.ann_samemod annsmA, a2.ann_samemod annsmB,
    a1.ann_exp expA, a2.ann_exp expB, a1.ann_htp htpA, a2.ann_htp htpB,
    a1.ann_alternate altA, a2.ann_alternate altB,
    a1.n_ann_on_protein nannA, a2.n_ann_on_protein nannB,
    (t1a.nodeA IS NOT NULL) AS is_t1a,
    (t1b.nodeA IS NOT NULL) AS ptm_comment_names_partner,
    t2.shared_enzymes,
    t4.pp_evidence,
    tr.position_gap, coalesce(tr.same_peptide_std,false) AS same_peptide,
    coalesce(tr.same_peptide_modblocked,false) AS same_peptide_modblocked,
    coalesce(n.in_network,false) AS in_network,
    cA.LeidenCluster AS clusterA, cB.LeidenCluster AS clusterB
  FROM base b
  JOIN A a1 ON a1.name=b.nodeA
  JOIN A a2 ON a2.name=b.nodeB
  LEFT JOIN T1A t1a ON t1a.nodeA=b.nodeA AND t1a.nodeB=b.nodeB
  LEFT JOIN T1B t1b ON t1b.nodeA=b.nodeA AND t1b.nodeB=b.nodeB
  LEFT JOIN T2  t2  ON t2.nodeA =b.nodeA AND t2.nodeB =b.nodeB
  LEFT JOIN TRY tr  ON tr.nodeA =b.nodeA AND tr.nodeB =b.nodeB
  LEFT JOIN NET n   ON n.nodeA  =b.nodeA AND n.nodeB  =b.nodeB
  LEFT JOIN CL  cA  ON cA.name  =b.nodeA
  LEFT JOIN CL  cB  ON cB.name  =b.nodeB
  LEFT JOIN T4  t4  ON t4.accA = least(a1.acc,a2.acc) AND t4.accB = greatest(a1.acc,a2.acc)
)
SELECT *,
  (accA=accB) AS same_protein, (accA=accB AND posA=posB) AS same_site, (modA=modB) AS same_mod,
  -- Branches run strongest evidence first and CASE takes the first match, so each pair gets the
  -- best tier it qualifies for. Nothing here reads Score or qvalue: tiers are assigned before any
  -- score is consulted, so the classification cannot be influenced by the quantity being tested.
  CASE
    WHEN is_t1a                                   THEN 1  -- both sites cited in one CC PTM block
    WHEN shared_enzymes IS NOT NULL               THEN 2  -- both sites "; by <same enzyme>"
    WHEN accA=accB AND annA AND annB              THEN 3  -- same protein, both annotated
    WHEN accA<>accB AND pp_evidence IS NOT NULL   THEN 4  -- interaction or shared complex
    WHEN annA AND annB AND (expA OR htpA) AND (expB OR htpB) THEN 5  -- type (a) only, no link
    ELSE 9 END AS tier,
  -- Null-safe: a pair counts as within-cluster only when both endpoints carry a label. Labels come
  -- from CLUSTERS, the weighted Leiden partition; the withdrawn unweighted file is not used.
  (clusterA IS NOT NULL AND clusterA=clusterB) AS same_cluster
FROM j
) TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000);
""")
log(f"  written in {time.time()-t0:.0f}s, {pathlib.Path(out).stat().st_size/1e9:.2f} GB")

log("  tier distribution, significant vs non-significant:")
d = con.execute(f"""SELECT tier, count(*) FILTER (is_sig) n_sig, count(*) FILTER (NOT is_sig) n_ns
                    FROM '{out}' GROUP BY 1 ORDER BY 1""").df()
tot_s = d.n_sig.sum(); tot_n = d.n_ns.sum()
LBL={1:"1 explicit crosstalk",2:"2 shared regulator",3:"3 same protein",4:"4 complex/interaction",5:"5 both annotated",9:"unsupported"}
log(f"    {'tier':<24}{'significant':>14}{'%':>8}{'non-signif':>13}{'%':>8}")
for _,r in d.iterrows():
    log(f"    {LBL[int(r.tier)]:<24}{int(r.n_sig):>14,}{100*r.n_sig/tot_s:>7.3f}%{int(r.n_ns):>13,}{100*r.n_ns/tot_n:>7.3f}%")
log(f"    {'TOTAL':<24}{int(tot_s):>14,}{'':>8}{int(tot_n):>13,}")
sup_s = d[d.tier<=5].n_sig.sum(); sup_n = d[d.tier<=5].n_ns.sum()
log(f"    supported (tier 1-5): significant {sup_s:,} ({100*sup_s/tot_s:.3f}%)  non-significant {sup_n:,} ({100*sup_n/tot_n:.3f}%)")
