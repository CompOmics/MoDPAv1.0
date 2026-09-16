"""
Step 4, part 4 and Step 5, part 1: stratified tier counts and the two built-in controls.

Reads:  work/pairs_tiered.parquet.
Writes: work/s14_tier_stratum_sign.csv, work/s14_tier_stratum_pivot.csv,
        work/s14_same_residue_pairs.csv, work/s14_copeptide_by_stratum.csv.
Feeds:  report sections 6.2, 7.1 and 7.2.

Control 1, competing modifications at one residue. Lysine can carry acetyl, ubiquityl, methyl and
succinyl marks that cannot coexist at a given moment, so if MoDPA detects mutual exclusivity these
pairs should be enriched for negative SDC. Only 66 of 7,578 sites carry more than one modification,
giving 71 testable pairs, and they turn out to be depleted for negative SDC rather than enriched.
The formal comparison against matched reference classes is in s18.

Control 2, shared measurement. Counts how the co-quantification flag from s10 distributes across
the |SDC| strata. It is negligible below 0.55 and reaches 61% above 0.70, which is reported rather
than filtered.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()
LBL={1:"T1 crosstalk",2:"T2 regulator",3:"T3 same protein",4:"T4 complex/PPI",5:"T5 both annot",9:"unsupported"}

log("=== STEP 4d: TIER x STRATUM x SIGN (significant pairs) ===")
d = con.execute(f"""SELECT stratum, sign, tier, count(*) n FROM '{pt}' WHERE is_sig GROUP BY 1,2,3""").df()
d["tierlab"] = d.tier.map(LBL)
d.to_csv(WORK/"s14_tier_stratum_sign.csv", index=False)
tot = d.groupby(["stratum","sign"]).n.sum().rename("N")
piv = d.pivot_table(index=["stratum","sign"], columns="tier", values="n", fill_value=0).join(tot)
for t in [1,2,3,4,5,9]:
    if t not in piv.columns: piv[t]=0
piv["supported"] = piv[[1,2,3,4,5]].sum(axis=1)
piv["pct_sup"] = 100*piv.supported/piv.N
piv["pct_t14"] = 100*piv[[1,2,3,4]].sum(axis=1)/piv.N
log(f"  {'stratum':<15}{'sgn':<5}{'N':>12}{'T1':>6}{'T2':>5}{'T3':>7}{'T4':>8}{'T5':>11}{'%sup':>8}{'%T1-4':>8}")
for (s,sg),r in piv.iterrows():
    log(f"  {s:<15}{sg:<5}{int(r.N):>12,}{int(r[1]):>6}{int(r[2]):>5}{int(r[3]):>7,}{int(r[4]):>8,}{int(r[5]):>11,}{r.pct_sup:>7.2f}%{r.pct_t14:>7.3f}%")
piv.to_csv(WORK/"s14_tier_stratum_pivot.csv")

log("=== STEP 5a: competing modifications at the same residue ===")
q = con.execute(f"""
 SELECT count(*) n, count(*) FILTER (sign='neg') n_neg, avg(Score) mean_score, median(Score) med
 FROM '{pt}' WHERE is_sig AND same_site AND modA<>modB""").fetchone()
log(f"  same-residue, different-modification significant pairs: n={q[0]}, negative={q[1]} ({100*q[1]/q[0]:.1f}%), mean SDC={q[2]:.4f}")
allp = con.execute(f"SELECT count(*) n, count(*) FILTER (sign='neg') neg, avg(Score) m FROM '{pt}' WHERE is_sig").fetchone()
log(f"  all significant pairs:                                  n={allp[0]:,}, negative={allp[1]:,} ({100*allp[1]/allp[0]:.1f}%), mean SDC={allp[2]:.4f}")
sp = con.execute(f"""SELECT count(*) n, count(*) FILTER (sign='neg') neg, avg(Score) m
                     FROM '{pt}' WHERE is_sig AND same_protein AND NOT same_site""").fetchone()
log(f"  same protein, different residue:                        n={sp[0]:,}, negative={sp[1]:,} ({100*sp[1]/sp[0]:.1f}%), mean SDC={sp[2]:.4f}")
det = con.execute(f"""SELECT nodeA,nodeB,modA,modB,Score,qvalue,in_network,same_peptide
                      FROM '{pt}' WHERE is_sig AND same_site AND modA<>modB ORDER BY Score""").df()
det["modA_name"]=det.modA.map(lambda u:UNIMOD[u][0]); det["modB_name"]=det.modB.map(lambda u:UNIMOD[u][0])
det.to_csv(WORK/"s14_same_residue_pairs.csv", index=False)
log(f"  modification-pair composition: {det.groupby(['modA_name','modB_name']).size().to_dict()}")
log(f"  of these, in published network: {int(det.in_network.sum())}")

log("=== STEP 5b: shared-measurement artefact pairs (own tryptic rule) ===")
a = con.execute(f"""
 SELECT count(*) n_same_prot,
        count(*) FILTER (same_peptide) n_copeptide,
        count(*) FILTER (same_peptide AND same_mod) n_copeptide_samemod,
        count(*) FILTER (same_protein AND position_gap<=5) n_gap5
 FROM '{pt}' WHERE is_sig AND same_protein""").fetchone()
log(f"  significant same-protein pairs: {a[0]:,}")
log(f"    can share a tryptic peptide (<=2 MC): {a[1]:,}")
log(f"    of those, same modification type:     {a[2]:,}")
log(f"    within 5 residues (the provided rule):{a[3]:,}")
b = con.execute(f"""SELECT count(*) FILTER (same_peptide) cop, count(*) tot,
   avg(Score) FILTER (same_peptide) m_cop, avg(Score) FILTER (same_protein AND NOT same_peptide) m_nocop
   FROM '{pt}' WHERE is_sig""").fetchone()
log(f"  co-peptide pairs as fraction of all significant pairs: {b[0]:,}/{b[1]:,}")
log(f"  mean SDC: co-peptide {b[2]:.4f} vs same-protein non-co-peptide {b[3]:.4f}")
c = con.execute(f"""SELECT stratum, count(*) n, count(*) FILTER (same_peptide) cop
   FROM '{pt}' WHERE is_sig GROUP BY 1 ORDER BY 1""").df()
log("  co-peptide enrichment by stratum:")
for _,r in c.iterrows():
    log(f"    {r.stratum:<15} {int(r.n):>12,}  co-peptide {int(r.cop):>6,}  ({100*r.cop/r.n:.4f}%)")
c.to_csv(WORK/"s14_copeptide_by_stratum.csv", index=False)
