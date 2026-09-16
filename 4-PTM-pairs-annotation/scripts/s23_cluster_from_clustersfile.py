"""
Step 8 re-derived: cluster analysis from the weighted partition alone.

Reads:  FILTERED, CLUSTERS, work/pairs_tiered.parquet.
Writes: work/s23_cluster_table.csv (all 556 clusters: size, within-cluster edges, per-tier counts,
        supported fraction, dominant PTM).
Feeds:  report sections 9.2 and 9.3. This supersedes the cluster-level part of s16.

Joins the network edges to the cluster file explicitly on the PTM event identifier, rather than
relying on the cluster columns already stored in pairs_tiered, and then checks the two agree. They
match on 100.0000% of the 14,591 edges with no unlabelled endpoint, which is what confirms that the
withdrawn unweighted partition never entered the analysis.

Two changes relative to s16: the per-cluster binomial test is corrected with Benjamini-Hochberg
across the 49 clusters that carry at least 20 within-cluster edges, and enriched and depleted
clusters are summarised by their dominant modification. The corrected result is that 10 clusters
are enriched and 13 depleted, and that the depleted set is not confined to methylation and
acetylation clusters as an earlier uncorrected reading of the top and bottom six suggested.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np
from scipy.stats import binomtest
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()
LBL={1:"T1 crosstalk",2:"T2 regulator",3:"T3 same protein",4:"T4 complex/PPI",5:"T5 both annot",9:"unsupported"}

log("=== STEP 8 RE-DERIVED: edges joined to Leiden-clusters.csv explicitly ===")
fd = pl.read_csv(FILTERED, infer_schema_length=10000).select(["nodeA","nodeB"]).to_pandas()
cl = pl.read_csv(CLUSTERS).select(["name","LeidenCluster","UniModID","Gene"]).to_pandas()
log(f"  edges from {FILTERED.name}: {len(fd):,}")
log(f"  clusters from {CLUSTERS.name}: {len(cl):,} nodes, {cl.LeidenCluster.nunique()} clusters")
lab = dict(zip(cl.name, cl.LeidenCluster))
fd["clusterA"]=fd.nodeA.map(lab); fd["clusterB"]=fd.nodeB.map(lab)
log(f"  edge endpoints with no cluster label: {int(fd.clusterA.isna().sum()+fd.clusterB.isna().sum())}")
fd["same_cluster"]=fd.clusterA==fd.clusterB
log(f"  within-cluster edges: {int(fd.same_cluster.sum()):,} / {len(fd):,} ({100*fd.same_cluster.mean():.2f}%)")

tiers = con.execute(f"SELECT nodeA,nodeB,tier,same_peptide,same_cluster AS sc_stored FROM '{pt}' WHERE in_network").df()
m = fd.merge(tiers, on=["nodeA","nodeB"], how="inner")
log(f"  joined to tier table: {len(m):,} rows; agreement with the stored same_cluster column: "
    f"{100*(m.same_cluster==m.sc_stored).mean():.4f}%")

log("  --- 8b CLUSTER CAPTURE (clusters file) ---")
log(f"    {'tier':<18}{'in network':>12}{'same cluster':>14}{'%':>9}{'co-peptide':>12}")
for t,g in m.groupby("tier"):
    log(f"    {LBL[int(t)]:<18}{len(g):>12,}{int(g.same_cluster.sum()):>14,}{100*g.same_cluster.mean():>8.2f}%{int(g.same_peptide.sum()):>12,}")

w = m[m.same_cluster].copy(); w["cluster"]=w.clusterA
piv = w.pivot_table(index="cluster", columns="tier", values="nodeA", aggfunc="count", fill_value=0)
for t in [1,2,3,4,5,9]:
    if t not in piv.columns: piv[t]=0
piv["edges_within"]=piv[[1,2,3,4,5,9]].sum(axis=1)
piv["sup14"]=piv[[1,2,3,4]].sum(axis=1)
piv["sup15"]=piv[[1,2,3,4,5]].sum(axis=1)
size = cl.groupby("LeidenCluster").size().rename("n_nodes")
dom = (cl.groupby(["LeidenCluster","UniModID"]).size().rename("n").reset_index()
         .sort_values(["LeidenCluster","n"],ascending=[True,False]).groupby("LeidenCluster").first())
dom["dominant_mod"]=dom.UniModID.map(lambda u: UNIMOD.get(int(u),("?",))[0])
dom["dominant_frac"]=dom.n/size
piv = piv.join(size).join(dom[["dominant_mod","dominant_frac"]])
piv["frac_sup14"]=piv.sup14/piv.edges_within
piv["frac_sup15"]=piv.sup15/piv.edges_within
g14 = piv.sup14.sum()/piv.edges_within.sum()
piv["enrich"]=piv.frac_sup14/g14
big = piv[piv.edges_within>=20].copy()
big["p_binom"]=[binomtest(int(r.sup14),int(r.edges_within),g14).pvalue for _,r in big.iterrows()]
from scipy.stats import false_discovery_control as fdr
big["q_binom"]=fdr(big.p_binom.values)
big=big.sort_values("frac_sup14",ascending=False)
piv.sort_values("edges_within",ascending=False).to_csv(WORK/"s23_cluster_table.csv")
log("  --- 8c CLUSTER-LEVEL ANNOTATION (clusters file) ---")
log(f"    within-cluster network edges={int(piv.edges_within.sum()):,}; tier1-4 supported={int(piv.sup14.sum()):,} ({100*g14:.2f}%)")
log(f"    clusters with >=1 within-cluster edge = {len(piv)}; with >=20 = {len(big)}")
log(f"    clusters containing at least one tier1-4 supported edge = {int((piv.sup14>0).sum())}")
log(f"    {'clu':>5}{'nodes':>7}{'edges':>7}{'sup14':>7}{'%sup':>8}{'x avg':>7}{'q(BH)':>10}  dominant PTM")
for c,r in big.head(10).iterrows():
    log(f"    {int(c):>5}{int(r.n_nodes):>7}{int(r.edges_within):>7}{int(r.sup14):>7}{100*r.frac_sup14:>7.2f}%{r.enrich:>7.2f}{r.q_binom:>10.2g}  {r.dominant_mod} ({100*r.dominant_frac:.0f}%)")
log("      ...")
for c,r in big.tail(8).iterrows():
    log(f"    {int(c):>5}{int(r.n_nodes):>7}{int(r.edges_within):>7}{int(r.sup14):>7}{100*r.frac_sup14:>7.2f}%{r.enrich:>7.2f}{r.q_binom:>10.2g}  {r.dominant_mod} ({100*r.dominant_frac:.0f}%)")
log(f"    significantly enriched at q<0.05: {int(((big.q_binom<0.05)&(big.enrich>1)).sum())}; "
    f"significantly depleted: {int(((big.q_binom<0.05)&(big.enrich<1)).sum())}")
log("  --- dominant PTM of enriched vs depleted clusters ---")
en = big[(big.q_binom<0.05)&(big.enrich>1)]; de = big[(big.q_binom<0.05)&(big.enrich<1)]
log(f"    enriched clusters dominant PTM: {en.dominant_mod.value_counts().to_dict()}")
log(f"    depleted clusters dominant PTM: {de.dominant_mod.value_counts().to_dict()}")
