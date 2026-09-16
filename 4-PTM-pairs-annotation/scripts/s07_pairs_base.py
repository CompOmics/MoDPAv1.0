"""
Step 2, part 6: materialise the significant pairs with strata and identifier-derived flags.

Reads:  work/raw_pairs.parquet.
Writes: work/sig_pairs_base.parquet (28,321,700 rows).
Feeds:  report section 3.

Splits each node identifier into accession, position, residue and Unimod inside duckdb, assigns the
|SDC| band, and records the flags that follow from the identifiers alone: same protein, same site,
same modification.

Superseded by s13, which repeats this join and adds the tier assignment. It is kept because it is
the cheapest way to get stratum counts without first building the evidence base.

The bands are fixed rather than quantile-based, for the reason given in s06.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, time

con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pq  = (WORK/"raw_pairs.parquet").as_posix()
out = (WORK/"sig_pairs_base.parquet").as_posix()

# |SDC| bands chosen from the observed shape (deciles are all inside 0.26-0.37 and
# would not resolve the published-network region at all)
BAND = """
 CASE WHEN abs(Score) >= 0.70 THEN 'S10_0.70+'
      WHEN abs(Score) >= 0.60 THEN 'S09_0.60-0.70'
      WHEN abs(Score) >= 0.55 THEN 'S08_0.55-0.60'
      WHEN abs(Score) >= 0.50 THEN 'S07_0.50-0.55'
      WHEN abs(Score) >= 0.45 THEN 'S06_0.45-0.50'
      WHEN abs(Score) >= 0.40 THEN 'S05_0.40-0.45'
      WHEN abs(Score) >= 0.35 THEN 'S04_0.35-0.40'
      WHEN abs(Score) >= 0.30 THEN 'S03_0.30-0.35'
      WHEN abs(Score) >= 0.25 THEN 'S02_0.25-0.30'
      ELSE 'S01_<0.25' END
"""
t0=time.time()
log("=== STEP 2f: significant pairs + strata + identifier-derived control flags ===")
con.execute(f"""
COPY (
  SELECT nodeA, nodeB, Score, pval, qvalue, PCC,
         abs(Score) AS absScore,
         CASE WHEN Score>=0 THEN 'pos' ELSE 'neg' END AS sign,
         {BAND} AS stratum,
         split_part(nodeA,'|',1) AS accA,
         CAST(split_part(nodeA,'|',2) AS INTEGER) AS posA,
         split_part(nodeA,'|',3) AS resA,
         CAST(split_part(nodeA,'|',4) AS INTEGER) AS modA,
         split_part(nodeB,'|',1) AS accB,
         CAST(split_part(nodeB,'|',2) AS INTEGER) AS posB,
         split_part(nodeB,'|',3) AS resB,
         CAST(split_part(nodeB,'|',4) AS INTEGER) AS modB
  FROM '{pq}' WHERE qvalue < 0.05
) TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000);
""")
log(f"  written {out} in {time.time()-t0:.0f}s, {pathlib.Path(out).stat().st_size/1e9:.2f} GB")

log("  stratum x sign counts (significant pairs):")
df = con.execute(f"SELECT stratum, sign, count(*) n FROM '{out}' GROUP BY 1,2 ORDER BY 1,2").df()
df.to_csv(WORK/"s07_stratum_sign.csv", index=False)
tot=0
for _,r in df.iterrows():
    log(f"    {r.stratum:<15} {r['sign']}  {int(r.n):>12,}"); tot+=int(r.n)
log(f"    TOTAL {tot:,}")

log("  identifier-derived control flags over significant pairs:")
r = con.execute(f"""
 SELECT count(*) n,
   count(*) FILTER (accA=accB) same_protein,
   count(*) FILTER (accA=accB AND posA=posB) same_site,
   count(*) FILTER (modA=modB) same_mod,
   count(*) FILTER (accA=accB AND abs(posA-posB)<=50) within50
 FROM '{out}'""").fetchone()
for k,v in zip("n same_protein same_site same_mod within50".split(), r):
    log(f"    {k:<14} {v:,}")
