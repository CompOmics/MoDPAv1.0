"""
Step 5 continued: formal tests for the two built-in controls.

Reads:  work/pairs_tiered.parquet.
Writes: work/s18_tier_stratum_nocop.csv (the section 6.2 table with co-quantifiable pairs dropped,
        kept as a sensitivity check only), work/s18_sign_by_copeptide.csv.
Feeds:  report sections 6.2, 7.1 and 7.2.

The same-residue control is tested against four matched reference classes rather than against the
global negative rate, because same-residue pairs are by construction co-quantified and share a
protein. Fisher exact tests are run in both directions, so the report can state not only that the
predicted enrichment is absent but that the opposite is significant.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, pandas as pd, numpy as np
from scipy.stats import fisher_exact
con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pt = (WORK/"pairs_tiered.parquet").as_posix()

log("=== STEP 4d-bis: TIER x STRATUM x SIGN, CO-PEPTIDE PAIRS EXCLUDED (headline) ===")
d = con.execute(f"""SELECT stratum, sign, tier, count(*) n FROM '{pt}'
   WHERE is_sig AND NOT same_peptide GROUP BY 1,2,3""").df()
piv = d.pivot_table(index=["stratum","sign"], columns="tier", values="n", fill_value=0)
for t in [1,2,3,4,5,9]:
    if t not in piv.columns: piv[t]=0
piv["N"]=piv[[1,2,3,4,5,9]].sum(axis=1)
piv["pct_T1_4"]=100*piv[[1,2,3,4]].sum(axis=1)/piv.N
piv["pct_T1_5"]=100*piv[[1,2,3,4,5]].sum(axis=1)/piv.N
piv.to_csv(WORK/"s18_tier_stratum_nocop.csv")
log(f"  {'stratum':<15}{'sgn':<5}{'N':>12}{'T1':>5}{'T2':>5}{'T3':>7}{'T4':>8}{'T5':>11}{'%T1-4':>9}{'%T1-5':>8}")
for (s,sg),r in piv.iterrows():
    log(f"  {s:<15}{sg:<5}{int(r.N):>12,}{int(r[1]):>5}{int(r[2]):>5}{int(r[3]):>7,}{int(r[4]):>8,}{int(r[5]):>11,}{r.pct_T1_4:>8.3f}%{r.pct_T1_5:>7.2f}%")

log("=== STEP 5a-bis: competing modifications at the same residue, matched comparison ===")
q = lambda w: con.execute(f"SELECT count(*) n, count(*) FILTER (sign='neg') neg, avg(Score) m FROM '{pt}' WHERE is_sig AND {w}").fetchone()
rows = [
  ("same residue, different modification",        "same_site AND modA<>modB"),
  ("same protein, co-peptide, different mod",     "same_protein AND same_peptide AND NOT same_site AND modA<>modB"),
  ("same protein, co-peptide, same mod",          "same_protein AND same_peptide AND NOT same_site AND modA=modB"),
  ("same protein, not co-peptide, different mod", "same_protein AND NOT same_peptide AND modA<>modB"),
  ("different proteins, different mod",           "NOT same_protein AND modA<>modB"),
  ("all significant pairs",                       "TRUE"),
]
res={}
for lab,w in rows:
    n,neg,m = q(w); res[lab]=(n,neg)
    log(f"  {lab:<44} n={n:>12,}  negative={neg:>10,} ({100*neg/max(n,1):>5.1f}%)  mean SDC={m:+.4f}")
a = res["same residue, different modification"]
for ref in ["same protein, co-peptide, different mod","same protein, not co-peptide, different mod","different proteins, different mod"]:
    b = res[ref]
    odds,p = fisher_exact([[a[1], a[0]-a[1]], [b[1], b[0]-b[1]]], alternative="greater")
    log(f"  Fisher (same-residue MORE negative than '{ref}'): OR={odds:.3f} p={p:.4g}")
    odds2,p2 = fisher_exact([[a[1], a[0]-a[1]], [b[1], b[0]-b[1]]], alternative="less")
    log(f"  Fisher (same-residue LESS negative than that reference):      OR={odds2:.3f} p={p2:.4g}")

log("=== STEP 5b-bis: sign composition by co-peptide status ===")
d2 = con.execute(f"""SELECT same_peptide, same_protein, count(*) n, count(*) FILTER (sign='neg') neg,
   avg(Score) m FROM '{pt}' WHERE is_sig GROUP BY 1,2 ORDER BY 1,2""").df()
for _,r in d2.iterrows():
    log(f"  co-peptide={bool(r.same_peptide)} same_protein={bool(r.same_protein)}: n={int(r.n):>12,} negative={100*r.neg/r.n:>5.1f}% mean SDC={r.m:+.4f}")
d2.to_csv(WORK/"s18_sign_by_copeptide.csv", index=False)

log("=== network-level co-peptide contamination ===")
n = con.execute(f"""SELECT count(*) n, count(*) FILTER (same_peptide) cop, count(*) FILTER (same_protein) sp
   FROM '{pt}' WHERE in_network""").fetchone()
log(f"  published network edges: {n[0]:,}; same protein {n[2]:,}; can share a tryptic peptide {n[1]:,} ({100*n[1]/n[0]:.1f}%)")
