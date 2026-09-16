"""
Step 2, part 1: convert the 1.3 GB gzipped association list to Parquet once.

Reads:  RAW (20260911-1042-keen_raman-signed-distances.csv.gz).
Writes: work/raw_pairs.parquet (about 0.96 GB, zstd).

Every later step queries the Parquet copy with duckdb instead of re-parsing the CSV. Decompressing
and parsing the CSV costs about 18 s per pass; the Parquet copy is scanned in a second or two and
supports column pruning, which matters because most queries touch two or three of the seven
columns. The conversion is skipped if the output already exists.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import duckdb, time, os

out = WORK/"raw_pairs.parquet"
if out.exists():
    log(f"parquet exists, skipping conversion: {out}")
else:
    t0=time.time()
    con = duckdb.connect()
    con.execute("PRAGMA threads=12; PRAGMA memory_limit='24GB';")
    log("converting raw csv.gz -> parquet (zstd) ...")
    con.execute(f"""
        COPY (
          SELECT nodeA, nodeB,
                 CAST(Score AS DOUBLE)    AS Score,
                 CAST(pval AS DOUBLE)     AS pval,
                 CAST(distance AS DOUBLE) AS distance,
                 CAST(PCC AS DOUBLE)      AS PCC,
                 CAST(qvalue AS DOUBLE)   AS qvalue
          FROM read_csv_auto('{RAW.as_posix()}', header=true)
        ) TO '{out.as_posix()}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000);
    """)
    log(f"done in {time.time()-t0:.0f}s; parquet size={out.stat().st_size/1e9:.2f} GB")
