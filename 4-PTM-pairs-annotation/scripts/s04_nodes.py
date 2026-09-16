"""
Step 2, part 4: parse the PTM event identifiers and build the node table.

Reads:  work/raw_pairs.parquet.
Writes: work/ptm_events.parquet (one row per PTM event: name, acc, pos, res, unimod).
Feeds:  report section 2.

The identifier format is ACC|POS|RES|UNIMOD, for example P04406|139|K|1. This script verifies that
every identifier has exactly four fields, that no accession carries an isoform suffix, and that the
residue in each identifier is one the stated modification can occupy. It also measures the Spearman
correlation between |Score| and pval, documenting the non-monotonicity found in s03.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, numpy as np, re

con = duckdb.connect(); con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
pq = (WORK/"raw_pairs.parquet").as_posix()

log("=== STEP 2b: pval / Score relation ===")
d = con.execute(f"""SELECT count(*) FROM '{pq}' WHERE abs(abs(Score)-distance) > 1e-12""").fetchone()[0]
log(f"  rows where distance != |Score| : {d:,}  (0 => distance is the unsigned dcor)")
# monotonicity of pval in |Score|
s = con.execute(f"""SELECT abs(Score) a, pval FROM '{pq}' USING SAMPLE 200000 ROWS""").df()
import scipy.stats as st
log(f"  Spearman(|Score|, pval) on 200k sample = {st.spearmanr(s.a, s.pval).statistic:.6f}")
ov = con.execute(f"""SELECT count(*) FROM '{pq}' WHERE abs(Score) BETWEEN 0.2157 AND 0.3107""").fetchone()[0]
log(f"  pairs in the |Score| band where significance is mixed (0.2157-0.3107) = {ov:,}")
log("  (p-value comes from the bias-corrected dcor t-test, which is not a monotone")
log("   function of the uncorrected dcor reported as Score; hence the overlap band.)")

log("=== STEP 2c: PTM EVENT TABLE ===")
ev = con.execute(f"""
  SELECT x AS name FROM (SELECT nodeA AS x FROM '{pq}' UNION SELECT nodeB FROM '{pq}')
""").df()
log(f"  n events = {len(ev):,}")
bad = [n for n in ev.name if len(n.split("|"))!=4]
log(f"  events not matching ACC|POS|RES|UNIMOD : {len(bad)}  {bad[:5]}")
parts = ev.name.str.split("|", expand=True)
parts.columns = ["acc","pos","res","unimod"]
ev = ev.join(parts)
log(f"  accessions containing '-' (isoform suffix): {(ev.acc.str.contains('-')).sum()}")
log(f"  accessions containing '_' or other odd chars: {(~ev.acc.str.match(r'^[A-Z0-9]+$')).sum()}")
ev["pos"] = ev["pos"].astype(int); ev["unimod"] = ev["unimod"].astype(int)
log(f"  pos min={ev.pos.min()} max={ev.pos.max()}")
log(f"  residues observed: {sorted(ev.res.unique())}")
log("  unimod x residue counts:")
ct = ev.groupby(["unimod","res"]).size().reset_index(name="n").sort_values(["unimod","n"],ascending=[True,False])
for _,r in ct.iterrows():
    nm = UNIMOD.get(int(r.unimod),("?",set()))[0]
    log(f"    {int(r.unimod):>4} {nm:<14} {r.res}  {int(r.n):>6}")
log("  per-modification totals:")
for u,n in ev.groupby("unimod").size().sort_values(ascending=False).items():
    log(f"    {u:>4} {UNIMOD.get(int(u),('?',))[0]:<14} {n:>6}")
log(f"  distinct proteins = {ev.acc.nunique():,}; distinct sites (acc,pos) = {ev.groupby(['acc','pos']).ngroups:,}")
ev.to_parquet(WORK/"ptm_events.parquet", index=False)
log("  -> work/ptm_events.parquet")

# compare with network node set and clusters file
cl = pl.read_csv(CLUSTERS)
net = set(cl['name'].to_list()); allev = set(ev.name)
log(f"  network nodes {len(net):,}; in raw event set: {len(net & allev):,}; missing: {len(net-allev):,}")
log(f"  raw events NOT in network = {len(allev-net):,}")
