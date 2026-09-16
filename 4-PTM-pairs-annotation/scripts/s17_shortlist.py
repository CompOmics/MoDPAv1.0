"""
Step 7: build the shortlist of verification-ready supported associations.

Reads:  work/pairs_tiered.parquet, work/node_ann.parquet, the pair_t*.parquet tables,
        work/ev_htp_pmids.json.
Writes: work/tier14_annotated.parquet (all 156,024 Tier 1 to 4 significant pairs with evidence
        text), work/shortlist.parquet (the 165-pair shortlist).
Feeds:  report section 10 and the delivered shortlist TSV.

The shortlist is sampled up to 10 per (stratum, sign) cell and only then ranked within each cell,
by tier and then |SDC|. Ranking by |SDC| across the whole list would simply reproduce the published
network and would say nothing about what lies below the cutoff.

No filtering beyond q < 0.05 and tier <= 4. Co-quantifiable pairs are flagged, not removed.

evidence() resolves, for each pair, which table supplied its tier and returns the exact UniProt text
and the PubMed ids attached to it, so every row can be traced back to a line in the flat file. PMIDs
are then split into those UniProt attaches with a high-throughput ECO code and the rest, which is
what all_evidence_is_large_scale_MS reports: a pair supported only by large-scale MS surveys is not
independent of the data MoDPA was built from.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np, json

con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()
A   = pd.read_parquet(WORK/"node_ann.parquet").set_index("name")
T1A = pd.read_parquet(WORK/"pair_t1a.parquet")
T1B = pd.read_parquet(WORK/"pair_t1b.parquet")
T2  = pd.read_parquet(WORK/"pair_t2.parquet")
T4  = pd.read_parquet(WORK/"pair_t4.parquet")
htp = json.load(open(WORK/"ev_htp_pmids.json")); HTP=set(htp["counts"]); TITLES=htp["titles"]

log("=== STEP 7: SHORTLIST ===")
cand = con.execute(f"""SELECT nodeA,nodeB,Score,qvalue,stratum,sign,tier,in_network,same_cluster,
    clusterA,clusterB,geneA,geneB,modA,modB,same_protein,same_peptide,position_gap,
    annsmA,annsmB,expA,expB,htpA,htpB,ptm_comment_names_partner
  FROM '{pt}' WHERE is_sig AND tier<=4""").df()
log(f"  tier1-4 significant pairs (candidate pool): {len(cand):,}; non-co-peptide: {int((~cand.same_peptide).sum()):,}")

k1a = {(r.nodeA,r.nodeB):(r.ev_text,r.ev_pmids) for r in T1A.itertuples()}
k1b = {(r.nodeA,r.nodeB):(r.ev_text,r.ev_pmids) for r in T1B.itertuples()}
k2  = {(r.nodeA,r.nodeB):r.shared_enzymes for r in T2.itertuples()}
k4  = {(r.accA,r.accB):r.pp_evidence for r in T4.itertuples()}

def evidence(r):
    key=(r.nodeA,r.nodeB)
    if key in k1a: return "CC_PTM_same_protein", k1a[key][0], k1a[key][1]
    if key in k2:
        na, nb = A.loc[r.nodeA], A.loc[r.nodeB]
        return "shared_enzyme", f"shared enzyme(s): {k2[key]} || A: {na.ann_note} || B: {nb.ann_note}", \
               ";".join(sorted(set((na.ann_pmids+";"+nb.ann_pmids).strip(";").split(";"))-{""}))
    if r.tier==3:
        na, nb = A.loc[r.nodeA], A.loc[r.nodeB]
        return "same_protein_both_annotated", f"A: {na.ann_note} || B: {nb.ann_note}", \
               ";".join(sorted(set((na.ann_pmids+";"+nb.ann_pmids).strip(";").split(";"))-{""}))
    if key in k1b:
        na, nb = A.loc[r.nodeA], A.loc[r.nodeB]
        return "protein_level_crosstalk_statement", k1b[key][0], k1b[key][1]
    acc = tuple(sorted([r.nodeA.split("|")[0], r.nodeB.split("|")[0]]))
    na, nb = A.loc[r.nodeA], A.loc[r.nodeB]
    return "protein_interaction_or_complex", f"protein-pair evidence: {k4.get(acc,'')} || A: {na.ann_note} || B: {nb.ann_note}", \
           ";".join(sorted(set((na.ann_pmids+";"+nb.ann_pmids).strip(";").split(";"))-{""}))

ev = [evidence(r) for r in cand.itertuples()]
cand["evidence_class"] = [e[0] for e in ev]
cand["evidence_text"]  = [e[1] for e in ev]
cand["pmids"]          = [e[2] for e in ev]
cand["pmids_large_scale_MS"] = [";".join(p for p in s.split(";") if p in HTP) if s else "" for s in cand.pmids]
cand["pmids_other"]          = [";".join(p for p in s.split(";") if p and p not in HTP) if s else "" for s in cand.pmids]
cand["all_evidence_is_large_scale_MS"] = (cand.pmids!="") & (cand.pmids_other=="")
cand["no_pmid"] = cand.pmids==""
cand["weak_eco_only"] = (~cand.expA & ~cand.htpA) | (~cand.expB & ~cand.htpB)
cand["modA_name"]=cand.modA.map(lambda u:UNIMOD[u][0]); cand["modB_name"]=cand.modB.map(lambda u:UNIMOD[u][0])
cand["captured_by_network"]  = cand.in_network
cand["captured_by_cluster"]  = cand.in_network & cand.same_cluster
cand.to_parquet(WORK/"tier14_annotated.parquet", index=False)

pool = cand.copy()
log(f"  shortlist pool (no filtering beyond q<0.05 and tier<=4): {len(pool):,}; "
    f"co-quantifiable among them: {int(pool.same_peptide.sum()):,} (flagged, not removed)")
pool["rank_key"] = list(zip(pool.tier, -pool.Score.abs()))
parts=[]
for (s,sg), g in pool.groupby(["stratum","sign"]):
    g = g.sort_values(["tier","evidence_class"], ascending=[True,True])
    g = g.assign(absS=g.Score.abs()).sort_values(["tier","absS"], ascending=[True,False])
    parts.append(g.head(10))
short = pd.concat(parts).sort_values(["tier","stratum","sign"], ascending=[True,True,True])
log(f"  shortlist size = {len(short)} across {short.stratum.nunique()} strata and signs {sorted(short['sign'].unique())}")
log(f"  tier composition: {short.tier.value_counts().sort_index().to_dict()}")
log(f"  sign composition: {short['sign'].value_counts().to_dict()}")
log(f"  captured by published network: {int(short.captured_by_network.sum())} of {len(short)}")
log(f"  captured by network AND same cluster: {int(short.captured_by_cluster.sum())}")
log(f"  evidence backed only by large-scale MS studies: {int(short.all_evidence_is_large_scale_MS.sum())}")
log(f"  no PubMed identifier attached to the annotation: {int(short.no_pmid.sum())}")
short.to_parquet(WORK/"shortlist.parquet", index=False)

log("  tier-1 pairs in the shortlist (explicit crosstalk), all signs and strata:")
for r in short[short.tier==1].head(25).itertuples():
    log(f"    {r.nodeA} ({r.geneA}) x {r.nodeB} ({r.geneB})  SDC={r.Score:+.3f} {r.stratum} net={r.in_network} "
        f"class={r.evidence_class} pmids={r.pmids[:60]}")
