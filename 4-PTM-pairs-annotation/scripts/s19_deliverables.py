"""
Step 9: write deliverables 1 and 2.

Reads:  work/pairs_tiered.parquet, work/tier14_annotated.parquet.
Writes: 20260911-modpa-significant-pairs-annotated.tsv.gz (all 28,321,700 significant pairs,
        about 1.5 GB compressed, roughly 3.5 minutes to write),
        20260911-modpa-shortlist.tsv (165 rows).

Evidence text and PubMed ids are joined in only for Tier 1 to 4 rows, since the other 28 million
rows have nothing to attach. Everything else, including the control flags, the stratum, network
membership and both Leiden cluster labels, is present on every row.

For anything programmatic, use work/pairs_tiered.parquet instead of the TSV. The TSV exists because
it is the requested deliverable format, not because it is a convenient one at this size.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np, time
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()
T14 = pd.read_parquet(WORK/"tier14_annotated.parquet")[
   ["nodeA","nodeB","evidence_class","evidence_text","pmids","pmids_large_scale_MS","pmids_other",
    "all_evidence_is_large_scale_MS","no_pmid"]]
con.register("T14", T14)
MODNAME = "CASE " + " ".join(f"WHEN {u}={u} AND m={u} THEN '{n[0]}'" for u,n in UNIMOD.items()) + " END"
umcase = lambda c: "CASE " + " ".join(f"WHEN {c}={u} THEN '{v[0]}'" for u,v in UNIMOD.items()) + " ELSE '?' END"

out1 = BASE/f"{DATE}-modpa-significant-pairs-annotated.tsv.gz"
t0=time.time(); log("=== STEP 9: DELIVERABLE 1 (all significant pairs) ===")
con.execute(f"""
COPY (
 SELECT p.nodeA, p.nodeB, p.geneA, p.geneB,
        {umcase('p.modA')} AS modA_name, {umcase('p.modB')} AS modB_name,
        p.Score AS SDC, p.qvalue, p.absScore, p.sign, p.stratum, p.tier,
        CASE p.tier WHEN 1 THEN 'T1_explicit_crosstalk' WHEN 2 THEN 'T2_shared_regulator'
                    WHEN 3 THEN 'T3_same_protein_both_annotated' WHEN 4 THEN 'T4_complex_or_interaction'
                    WHEN 5 THEN 'T5_both_sites_annotated' ELSE 'unsupported' END AS tier_label,
        coalesce(t.evidence_class,'') AS evidence_class,
        coalesce(t.evidence_text,'')  AS evidence_text,
        coalesce(t.pmids,'')          AS pmids,
        coalesce(t.pmids_large_scale_MS,'') AS pmids_large_scale_MS,
        coalesce(t.pmids_other,'')          AS pmids_primary_or_other,
        coalesce(t.all_evidence_is_large_scale_MS,false) AS all_evidence_is_large_scale_MS,
        p.ptm_comment_names_partner,
        p.same_protein, p.same_site, p.same_mod, p.position_gap,
        p.same_peptide AS ctrl_shared_tryptic_peptide,
        p.same_peptide_modblocked AS ctrl_shared_peptide_modblocked,
        p.nrunsA, p.nrunsB, p.nannA AS n_annot_on_proteinA, p.nannB AS n_annot_on_proteinB,
        p.in_network AS in_published_network_SDC_ge_0_6,
        p.clusterA AS leiden_clusterA, p.clusterB AS leiden_clusterB, p.same_cluster
 FROM '{pt}' p LEFT JOIN T14 t ON t.nodeA=p.nodeA AND t.nodeB=p.nodeB
 WHERE p.is_sig
) TO '{out1.as_posix()}' (FORMAT CSV, DELIMITER '\t', HEADER, COMPRESSION GZIP);
""")
log(f"  {out1.name}: {out1.stat().st_size/1e6:.0f} MB in {time.time()-t0:.0f}s")

log("=== DELIVERABLE 2 (shortlist) ===")
sh = pd.read_parquet(WORK/"shortlist.parquet")
cols = ["nodeA","nodeB","geneA","geneB","modA_name","modB_name","Score","qvalue","stratum","sign","tier",
        "evidence_class","evidence_text","pmids","pmids_large_scale_MS","pmids_other",
        "all_evidence_is_large_scale_MS","no_pmid","weak_eco_only","ptm_comment_names_partner",
        "same_protein","position_gap","same_peptide","captured_by_network","captured_by_cluster",
        "clusterA","clusterB"]
sh = sh[cols].rename(columns={"Score":"SDC","same_peptide":"ctrl_shared_tryptic_peptide",
        "pmids_other":"pmids_primary_or_other","weak_eco_only":"weak_eco_only_on_one_site"})
sh = sh.sort_values(["tier","stratum","sign","SDC"], ascending=[True,True,True,False])
out2 = BASE/f"{DATE}-modpa-shortlist.tsv"
sh.to_csv(out2, sep="\t", index=False)
log(f"  {out2.name}: {len(sh)} rows")
log(f"  tier composition {sh.tier.value_counts().sort_index().to_dict()}; "
    f"captured by network {int(sh.captured_by_network.sum())}; "
    f"evidence only from large-scale MS {int(sh.all_evidence_is_large_scale_MS.sum())}")
