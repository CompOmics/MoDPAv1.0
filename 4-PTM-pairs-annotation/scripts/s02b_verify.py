"""
Step 2, part 2: confirm the Parquet conversion preserved every row.

Reads:  work/raw_pairs.parquet.
Writes: nothing; the row count goes to work/runlog.md.

Expected: 29,211,546 rows, which is 7,644 x 7,643 / 2, the complete upper triangle over the 7,644
PTM events.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, time
con = duckdb.connect(); con.execute("PRAGMA threads=12;")
pq = (WORK/"raw_pairs.parquet").as_posix()
n = con.execute(f"SELECT count(*) FROM '{pq}'").fetchone()[0]
log(f"parquet rows = {n:,}")
