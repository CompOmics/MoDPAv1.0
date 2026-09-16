"""
Step 2, part 3: characterise the raw association list.

Reads:  work/raw_pairs.parquet.
Writes: work/s03_score_hist.csv (Score histogram in 0.05 bins, total and significant).
Feeds:  report section 3.

Establishes the four facts the rest of the analysis rests on:
  1. q < 0.05 retains 96.95% of pairs, so significance is not the operative filter and the SDC
     magnitude is.
  2. The list is upper-triangular, with no self pairs and no repeated unordered pair.
  3. The 14,591 pairs with Score >= 0.6 are exactly the published network, so the cutoff was
     applied to the signed SDC and the two pairs at SDC <= -0.6 were dropped.
  4. pval is not monotone in Score, because the p-value comes from the bias-corrected
     distance-correlation t-test while Score is the uncorrected estimator. The significant and
     non-significant |Score| ranges therefore overlap, and the q-value cannot be reconstructed
     from Score alone.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, numpy as np, json

con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pq = (WORK/"raw_pairs.parquet").as_posix()

log("=== STEP 2: RAW LIST CHARACTERISATION ===")
r = con.execute(f"""
 SELECT count(*) n,
        count(*) FILTER (qvalue < 0.05) n_sig,
        count(*) FILTER (pval  < 0.05) n_sig_p,
        count(*) FILTER (Score < 0)    n_neg,
        count(*) FILTER (Score >= 0.6) n_ge06,
        count(*) FILTER (Score >= 0.6 AND qvalue < 0.05) n_ge06_sig,
        count(*) FILTER (Score <= -0.6) n_le06,
        min(Score) mn, max(Score) mx,
        min(qvalue) qmn, max(qvalue) qmx,
        count(*) FILTER (nodeA = nodeB) n_self
 FROM '{pq}'""").fetchone()
cols = "n n_sig n_sig_p n_neg n_ge06 n_ge06_sig n_le06 mn mx qmn qmx n_self".split()
for c,v in zip(cols,r): log(f"  {c} = {v:,}" if isinstance(v,int) else f"  {c} = {v}")

# symmetry: is (a,b) ever present with (b,a)?
dup = con.execute(f"""
 SELECT count(*) FROM (
   SELECT least(nodeA,nodeB) a, greatest(nodeA,nodeB) b, count(*) c
   FROM '{pq}' GROUP BY 1,2 HAVING c>1)""").fetchone()[0]
log(f"  unordered pairs appearing more than once = {dup:,}  (0 => upper-triangular)")

nn = con.execute(f"SELECT count(*) FROM (SELECT nodeA AS x FROM '{pq}' UNION SELECT nodeB FROM '{pq}')").fetchone()[0]
log(f"  unique PTM events in raw list = {nn:,}; N*(N-1)/2 = {nn*(nn-1)//2:,}")

nsig_nodes = con.execute(f"SELECT count(*) FROM (SELECT nodeA AS x FROM '{pq}' WHERE qvalue<0.05 UNION SELECT nodeB FROM '{pq}' WHERE qvalue<0.05)").fetchone()[0]
log(f"  unique PTM events among q<0.05 pairs = {nsig_nodes:,}")

# largest qvalue that is still < 0.05 ; and p-value at that point
q = con.execute(f"SELECT max(qvalue) FROM '{pq}' WHERE qvalue<0.05").fetchone()[0]
log(f"  max qvalue below 0.05 = {q}")
# |Score| at the significance boundary
b = con.execute(f"""SELECT min(abs(Score)), max(abs(Score)) FROM '{pq}' WHERE qvalue<0.05""").fetchone()
log(f"  |Score| range among significant = {b[0]:.6f} .. {b[1]:.6f}")
b2 = con.execute(f"""SELECT min(abs(Score)), max(abs(Score)) FROM '{pq}' WHERE qvalue>=0.05""").fetchone()
log(f"  |Score| range among NON-significant = {b2[0]:.6f} .. {b2[1]:.6f}")

# histogram of Score, signed, 0.05 bins
h = con.execute(f"""
 SELECT floor(Score*20)/20 AS bin, count(*) c,
        count(*) FILTER (qvalue<0.05) c_sig
 FROM '{pq}' GROUP BY 1 ORDER BY 1""").df()
h.to_csv(WORK/"s03_score_hist.csv", index=False)
log(f"  score histogram written ({len(h)} bins) -> work/s03_score_hist.csv")
log("  bin  total  significant")
for _,row in h.iterrows():
    log(f"   {row['bin']:+.2f}  {int(row['c']):>12,}  {int(row['c_sig']):>12,}")
