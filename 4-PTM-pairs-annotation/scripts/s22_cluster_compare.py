"""
Side-by-side comparison of the two Leiden partitions.

Reads:  CLUSTERS, EDGES, work/pairs_tiered.parquet.
Writes: work/s22_cluster_table_clusters.csv, work/s22_cluster_table_edges.csv.

Was intended to quantify how much the Step 8 cluster results depend on which partition is used.

DOES NOT RUN: the EDGES file, which held the unweighted partition, has been withdrawn from the
folder. The comparison is therefore not part of the report, and s23 supersedes this script by
re-deriving the cluster analysis from the weighted partition alone.

If the unweighted file is restored, the right comparison is an adjusted Rand index and normalized
mutual information between the two partitions. Comparing raw cluster identifiers, as the first
inventory did, measures nothing, because the identifiers are arbitrary across independent runs.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np
from scipy.stats import binomtest
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()
LBL={1:"T1 crosstalk",2:"T2 regulator",3:"T3 same protein",4:"T4 complex/PPI",5:"T5 both annot",9:"unsupported"}

log("=== STEP 8 RECOMPUTED: clusters-file join, side by side with the edges-file labels ===")
cl = pl.read_csv(CLUSTERS).select(["name","LeidenCluster","UniModID"]).to_pandas()
ed = pl.read_csv(EDGES, infer_schema_length=10000).to_pandas()

# node -> label from each source
lab_cl = dict(zip(cl.name, cl.LeidenCluster))
lab_ed = {}
for a,ca in zip(ed.nodeA, ed.LeidenClusterA): lab_ed.setdefault(a, ca)
for b,cb in zip(ed.nodeB, ed.LeidenClusterB): lab_ed.setdefault(b, cb)
log(f"  clusters file: {len(lab_cl):,} nodes, {len(set(lab_cl.values()))} distinct labels")
log(f"  edges file   : {len(lab_ed):,} nodes, {len(set(lab_ed.values()))} distinct labels")
log(f"  node sets identical: {set(lab_cl)==set(lab_ed)}")

net = con.execute(f"""SELECT nodeA,nodeB,tier,same_peptide FROM '{pt}' WHERE in_network""").df()
net["cl_A"]=net.nodeA.map(lab_cl); net["cl_B"]=net.nodeB.map(lab_cl)
net["ed_A"]=net.nodeA.map(lab_ed); net["ed_B"]=net.nodeB.map(lab_ed)
net["same_cl"]=net.cl_A==net.cl_B
net["same_ed"]=net.ed_A==net.ed_B
log(f"  network edges: {len(net):,}")
log(f"  within-cluster edges: clusters file {int(net.same_cl.sum()):,}  |  edges file {int(net.same_ed.sum()):,}")
log(f"  edges where the two sources disagree on within/between: {int((net.same_cl!=net.same_ed).sum()):,}")

log("  --- 8b cluster capture by tier ---")
log(f"    {'tier':<18}{'in net':>8}{'same cl (clusters file)':>26}{'same cl (edges file)':>23}")
for t,g in net.groupby("tier"):
    log(f"    {LBL[int(t)]:<18}{len(g):>8,}{int(g.same_cl.sum()):>18,} ({100*g.same_cl.mean():5.1f}%){int(g.same_ed.sum()):>15,} ({100*g.same_ed.mean():5.1f}%)")

def cluster_table(labA, labB, tag):
    w = net[net[labA]==net[labB]].copy()
    w["cluster"]=w[labA]
    piv = w.pivot_table(index="cluster", columns="tier", values="nodeA", aggfunc="count", fill_value=0)
    for t in [1,2,3,4,5,9]:
        if t not in piv.columns: piv[t]=0
    piv["edges_within"]=piv[[1,2,3,4,5,9]].sum(axis=1)
    piv["sup14"]=piv[[1,2,3,4]].sum(axis=1)
    size = pd.Series({c:sum(1 for n,l in (lab_cl if tag=="clusters" else lab_ed).items() if l==c)
                      for c in piv.index})
    piv["n_nodes"]=size
    dom = {}
    src = lab_cl if tag=="clusters" else lab_ed
    um = dict(zip(cl.name, cl.UniModID))
    tmp = pd.DataFrame({"node":list(src), "c":[src[n] for n in src]})
    tmp["u"]=tmp.node.map(um)
    for c,g in tmp.groupby("c"):
        vc=g.u.value_counts(); dom[c]=(UNIMOD.get(int(vc.index[0]),("?",))[0], vc.iloc[0]/len(g))
    piv["dominant_mod"]=[dom.get(c,("?",0))[0] for c in piv.index]
    piv["dominant_frac"]=[dom.get(c,("?",0))[1] for c in piv.index]
    piv["frac_sup14"]=piv.sup14/piv.edges_within
    return piv

for tag,(a,b) in [("clusters",("cl_A","cl_B")), ("edges",("ed_A","ed_B"))]:
    piv = cluster_table(a,b,tag)
    g14 = piv.sup14.sum()/piv.edges_within.sum()
    big = piv[piv.edges_within>=20].copy()
    big["enrich"]=big.frac_sup14/g14
    big["p"]=[binomtest(int(r.sup14),int(r.edges_within),g14).pvalue for _,r in big.iterrows()]
    big=big.sort_values("frac_sup14",ascending=False)
    log(f"  --- 8c cluster-level, source = {tag} file ---")
    log(f"    within-cluster edges={int(piv.edges_within.sum()):,}  tier1-4 supported={int(piv.sup14.sum()):,} ({100*g14:.2f}%)"
        f"  clusters with >=20 edges={len(big)}")
    log(f"    {'clu':>5}{'nodes':>7}{'edges':>7}{'sup14':>7}{'%sup':>8}{'x avg':>7}{'p':>10}  dominant PTM")
    for c,r in big.head(6).iterrows():
        log(f"    {int(c):>5}{int(r.n_nodes):>7}{int(r.edges_within):>7}{int(r.sup14):>7}{100*r.frac_sup14:>7.2f}%{r.enrich:>7.2f}{r.p:>10.2g}  {r.dominant_mod} ({100*r.dominant_frac:.0f}%)")
    log("      ...")
    for c,r in big.tail(6).iterrows():
        log(f"    {int(c):>5}{int(r.n_nodes):>7}{int(r.edges_within):>7}{int(r.sup14):>7}{100*r.frac_sup14:>7.2f}%{r.enrich:>7.2f}{r.p:>10.2g}  {r.dominant_mod} ({100*r.dominant_frac:.0f}%)")
    piv.to_csv(WORK/f"s22_cluster_table_{tag}.csv")
