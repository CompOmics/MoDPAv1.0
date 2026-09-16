"""
Step 2, part 5: |SDC| strata and the detection-frequency covariate.

Reads:  work/raw_pairs.parquet, NRUNS, work/ptm_events.parquet.
Writes: work/ptm_events_nruns.parquet (PTM events joined to n_runs and classification).
Feeds:  report sections 2 and 3.

Computes deciles of |SDC| among significant pairs and shows they all fall between 0.261 and 0.367.
That is why the analysis uses fixed |SDC| bands rather than deciles: deciles cannot resolve the
region above 0.6 where the published network lives.

Also derives the Unimod id to name mapping from the ptm_name column of N-runs-per-ptm.csv.gz, and
confirms that all 7,644 PTM events match that file on (accession, position, residue, Unimod), which
makes n_runs usable as the detection-frequency covariate in the s15 matched null.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, numpy as np, pandas as pd

con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pq = (WORK/"raw_pairs.parquet").as_posix()

log("=== STEP 2d: |SDC| quantiles among significant pairs ===")
qs = [0.1*i for i in range(1,10)]
sel = ", ".join(f"quantile_cont(abs(Score), {q}) q{int(q*100)}" for q in qs)
row = con.execute(f"SELECT {sel} FROM '{pq}' WHERE qvalue<0.05").fetchdf().iloc[0]
log("  deciles of |SDC| (significant): " + ", ".join(f"{v:.4f}" for v in row.values))
for sgn,w in [("positive","Score>0"),("negative","Score<0")]:
    r2 = con.execute(f"SELECT count(*) n, min(abs(Score)) mn, max(abs(Score)) mx, avg(abs(Score)) av FROM '{pq}' WHERE qvalue<0.05 AND {w}").fetchone()
    log(f"  {sgn}: n={r2[0]:,} |SDC| {r2[1]:.4f}..{r2[2]:.4f} mean {r2[3]:.4f}")

log("=== STEP 2e: N-runs file cross-check ===")
nr = pl.read_csv(NRUNS)
# unimod id + name from ptm_name  '[21]Phospho'
nr = nr.with_columns([
    pl.col("ptm_name").str.extract(r"^\[(\d+)\]", 1).cast(pl.Int64).alias("unimod"),
    pl.col("ptm_name").str.extract(r"^\[\d+\](.*)$", 1).alias("modname"),
])
mp = (nr.select(["unimod","modname"]).unique().drop_nulls()
        .sort("unimod").to_pandas())
log(f"  distinct unimod ids in N-runs file = {len(mp)}")
sub = mp[mp.unimod.isin(UNIMOD.keys())]
log("  names for the 11 analysed modifications (authoritative, from ptm_name):")
for _,r in sub.iterrows():
    log(f"    {int(r.unimod):>4}  {r.modname}")
log(f"  classification values: {nr['classification'].unique().to_list()}")

ev = pd.read_parquet(WORK/"ptm_events.parquet")
nrp = nr.select(["UniAcc","ptm_loc","ptm_res","unimod","n_runs","mean_psm_counts","median_psm_counts","classification"]).to_pandas()
m = ev.merge(nrp, left_on=["acc","pos","res","unimod"], right_on=["UniAcc","ptm_loc","ptm_res","unimod"], how="left")
log(f"  events matched to N-runs on (acc,pos,res,unimod): {m.n_runs.notna().sum():,} / {len(ev):,}")
if m.n_runs.isna().any():
    mm = m[m.n_runs.isna()]
    log(f"  UNMATCHED examples: {mm.name.head(5).tolist()}")
    m2 = ev.merge(nrp.drop(columns=['ptm_res']), left_on=["acc","pos","unimod"], right_on=["UniAcc","ptm_loc","unimod"], how="left")
    log(f"  matched ignoring residue: {m2.n_runs.notna().sum():,}")
log(f"  n_runs among matched: min={m.n_runs.min()} median={m.n_runs.median()} max={m.n_runs.max()}")
log(f"  classification of matched events: {m.classification.value_counts().to_dict()}")
m = m.drop(columns=["UniAcc","ptm_loc","ptm_res"])
m.to_parquet(WORK/"ptm_events_nruns.parquet", index=False)
log("  -> work/ptm_events_nruns.parquet")
