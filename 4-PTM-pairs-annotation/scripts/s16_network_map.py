"""
Step 8: map the findings onto the published network and the Leiden clusters.

Reads:  work/pairs_tiered.parquet, CLUSTERS.
Writes: work/s16_retention_all.csv, work/s16_retention_nocop.csv, work/s16_cluster_capture.csv,
        work/s16_cluster_table.csv, work/s16_discarded.csv, work/s16_modtype_composition.csv.
Feeds:  report section 9.

Four separate questions, deliberately not pooled:
  8a retention     - for each tier, how much survives the 0.6 cutoff, in absolute counts on both
                     sides. The co-quantification view is printed as a sensitivity check only.
  8b cluster capture - among pairs that do survive, whether both events land in one cluster. A
                     supported pair split across two clusters is a different kind of miss from a
                     pair lost at the score cutoff.
  8c cluster level - which clusters hold supported pairs, tested against the global rate with a
                     binomial test. Cluster labels come from CLUSTERS, the weighted partition.
  8d and 8e        - what the cutoff discards, by stratum and by modification pair.
  8f               - what the cutoff retains that has no annotation support at any tier.

Superseded for the cluster-level part by s23, which re-derives the join from the cluster file alone
and adds Benjamini-Hochberg correction across the 49 testable clusters. Use s23 for section 9.3.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()
LBL={1:"T1 crosstalk",2:"T2 regulator",3:"T3 same protein",4:"T4 complex/PPI",5:"T5 both annot",9:"unsupported"}

log("=== STEP 8: MAPPING ONTO THE PUBLISHED NETWORK AND CLUSTERS ===")
log("--- 8a RETENTION at SDC >= 0.6 (headline: no filtering beyond q<0.05) ---")
for tag, extra in [("HEADLINE: all significant pairs",""), ("sensitivity only: co-peptide pairs dropped"," AND NOT same_peptide")]:
    d = con.execute(f"""SELECT tier,
          count(*) n, count(*) FILTER (in_network) n_in, count(*) FILTER (NOT in_network) n_out
        FROM '{pt}' WHERE is_sig{extra} GROUP BY 1 ORDER BY 1""").df()
    log(f"  [{tag}]")
    log(f"    {'tier':<18}{'significant':>14}{'in network':>12}{'below 0.6':>14}{'retained %':>12}")
    for _,r in d.iterrows():
        log(f"    {LBL[int(r.tier)]:<18}{int(r.n):>14,}{int(r.n_in):>12,}{int(r.n_out):>14,}{100*r.n_in/r.n:>11.4f}%")
    d.to_csv(WORK/f"s16_retention_{'nocop' if extra else 'all'}.csv", index=False)

log("--- 8b CLUSTER CAPTURE among pairs that survive 0.6 ---")
d = con.execute(f"""SELECT tier, count(*) n, count(*) FILTER (same_cluster) n_same,
     count(*) FILTER (same_peptide) n_cop
   FROM '{pt}' WHERE in_network GROUP BY 1 ORDER BY 1""").df()
log(f"    {'tier':<18}{'in network':>12}{'same cluster':>14}{'%':>9}{'co-peptide':>12}")
for _,r in d.iterrows():
    log(f"    {LBL[int(r.tier)]:<18}{int(r.n):>12,}{int(r.n_same):>14,}{100*r.n_same/r.n:>8.2f}%{int(r.n_cop):>12,}")
d.to_csv(WORK/"s16_cluster_capture.csv", index=False)

log("--- 8c CLUSTER-LEVEL: supported pairs per Leiden cluster ---")
cl = pl.read_csv(CLUSTERS).select(["name","LeidenCluster","UniModID","Gene","UniAcc"]).to_pandas()
size = cl.groupby("LeidenCluster").size().rename("n_nodes")
dom = (cl.groupby(["LeidenCluster","UniModID"]).size().rename("n").reset_index()
         .sort_values(["LeidenCluster","n"], ascending=[True,False])
         .groupby("LeidenCluster").first())
dom["dominant_mod"] = dom.UniModID.map(lambda u: UNIMOD.get(int(u),("?",))[0])
dom["dominant_frac"] = dom.n/size
w = con.execute(f"""SELECT clusterA AS cluster, tier, count(*) n FROM '{pt}'
   WHERE in_network AND same_cluster GROUP BY 1,2""").df()
piv = w.pivot_table(index="cluster", columns="tier", values="n", fill_value=0)
for t in [1,2,3,4,5,9]:
    if t not in piv.columns: piv[t]=0
piv["edges_within"] = piv[[1,2,3,4,5,9]].sum(axis=1)
piv["supported_14"] = piv[[1,2,3,4]].sum(axis=1)
piv["supported_15"] = piv[[1,2,3,4,5]].sum(axis=1)
piv = piv.join(size).join(dom[["dominant_mod","dominant_frac"]])
piv["frac_sup14"] = piv.supported_14/piv.edges_within
piv["frac_sup15"] = piv.supported_15/piv.edges_within
piv = piv.sort_values("edges_within", ascending=False)
piv.to_csv(WORK/"s16_cluster_table.csv")
g14 = piv.supported_14.sum()/piv.edges_within.sum()
log(f"  within-cluster network edges = {int(piv.edges_within.sum()):,}; tier1-4 supported = {int(piv.supported_14.sum()):,} ({100*g14:.2f}%)")
log(f"  clusters with >=20 within-cluster edges: {(piv.edges_within>=20).sum()}")
big = piv[piv.edges_within>=20].copy()
from scipy.stats import binomtest
big["p_binom"] = [binomtest(int(r.supported_14), int(r.edges_within), g14).pvalue for _,r in big.iterrows()]

big["enrich"] = big.frac_sup14/g14
big = big.sort_values("frac_sup14", ascending=False)
log(f"  top 12 clusters by tier1-4 supported fraction (>=20 within-cluster edges):")
log(f"    {'clu':>5}{'nodes':>7}{'edges':>8}{'sup14':>7}{'%sup':>8}{'x avg':>7}{'p':>10}  dominant PTM")
for c,r in big.head(12).iterrows():
    log(f"    {int(c):>5}{int(r.n_nodes):>7}{int(r.edges_within):>8}{int(r.supported_14):>7}{100*r.frac_sup14:>7.2f}%{r.enrich:>7.2f}{r.p_binom:>10.2g}  {r.dominant_mod} ({100*r.dominant_frac:.0f}%)")
log(f"  bottom 6:")
for c,r in big.tail(6).iterrows():
    log(f"    {int(c):>5}{int(r.n_nodes):>7}{int(r.edges_within):>8}{int(r.supported_14):>7}{100*r.frac_sup14:>7.2f}%{r.enrich:>7.2f}{r.p_binom:>10.2g}  {r.dominant_mod} ({100*r.dominant_frac:.0f}%)")

log("--- 8d WHAT THE CUTOFF DISCARDS: tier1-4 supported pairs below 0.6 (no filtering) ---")
d = con.execute(f"""SELECT stratum, sign, tier, count(*) n FROM '{pt}'
   WHERE is_sig AND NOT in_network AND tier<=4 GROUP BY 1,2,3 ORDER BY 1,2,3""").df()
d.to_csv(WORK/"s16_discarded.csv", index=False)
tot = con.execute(f"""SELECT count(*) FROM '{pt}' WHERE is_sig AND NOT in_network AND tier<=4""").fetchone()[0]
inn = con.execute(f"""SELECT count(*) FROM '{pt}' WHERE in_network AND tier<=4""").fetchone()[0]
log(f"  tier1-4 pairs (no filtering): below 0.6 = {tot:,}; at/above 0.6 = {inn:,}")
tc = con.execute(f"""SELECT count(*) FILTER (same_peptide) FROM '{pt}' WHERE is_sig AND NOT in_network AND tier<=4""").fetchone()[0]
log(f"    of the discarded tier1-4 pairs, {tc:,} are co-quantifiable (flag only, not filtered)")
p = d.pivot_table(index="stratum", columns="tier", values="n", aggfunc="sum", fill_value=0)
log(f"    {'stratum':<15}" + "".join(f"{'T'+str(t):>9}" for t in p.columns) + f"{'total':>10}")
for s,r in p.iterrows():
    log(f"    {s:<15}" + "".join(f"{int(v):>9,}" for v in r.values) + f"{int(r.sum()):>10,}")

log("--- 8e modification-type composition, network vs discarded supported pairs ---")
mm = con.execute(f"""SELECT
     CASE WHEN modA<=modB THEN modA ELSE modB END m1,
     CASE WHEN modA<=modB THEN modB ELSE modA END m2,
     count(*) FILTER (in_network) n_net,
     count(*) FILTER (is_sig AND NOT in_network AND tier<=4) n_disc_sup
   FROM '{pt}' GROUP BY 1,2 HAVING n_net>0 OR n_disc_sup>0 ORDER BY n_disc_sup DESC""").df()
mm["pair"] = [f"{UNIMOD[a][0]}-{UNIMOD[b][0]}" for a,b in zip(mm.m1,mm.m2)]
mm.to_csv(WORK/"s16_modtype_composition.csv", index=False)
log(f"    {'modification pair':<26}{'in network':>12}{'discarded T1-4':>16}")
for _,r in mm.head(14).iterrows():
    log(f"    {r['pair']:<26}{int(r.n_net):>12,}{int(r.n_disc_sup):>16,}")

log("--- 8f WHAT THE CUTOFF RETAINS THAT IS UNSUPPORTED ---")
u = con.execute(f"""SELECT count(*) n, count(*) FILTER (NOT same_peptide) n_nocop,
   count(*) FILTER (same_cluster) n_samecl FROM '{pt}' WHERE in_network AND tier=9""").fetchone()
tot_net = con.execute(f"SELECT count(*) FROM '{pt}' WHERE in_network").fetchone()[0]
log(f"  network edges with no annotation support at any tier: {u[0]:,} of {tot_net:,} ({100*u[0]/tot_net:.1f}%)")
log(f"    of those, not co-peptide: {u[1]:,}; within a single cluster: {u[2]:,}")
