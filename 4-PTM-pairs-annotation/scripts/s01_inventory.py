"""
Step 1, part 1: inventory every input file before any analysis code is written.

Reads:  all five input files.
Writes: nothing; everything goes to work/runlog.md through log().

Reports column names, dtypes, row counts, the cluster-size distribution, the node sets of the
network and of the cluster file, and the first few lines of the raw association list, so that the
PTM event identifier format can be read off the data rather than assumed.

REQUIRES the withdrawn EDGES file and will fail with FileNotFoundError until it is restored.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import numpy as np, gzip, itertools, json

log("=== STEP 1 INVENTORY ===")

# ---------- Leiden clusters ----------
cl = pl.read_csv(CLUSTERS)
log(f"CLUSTERS {CLUSTERS.name}: shape={cl.shape} cols={cl.columns}")
sizes = cl.group_by('LeidenCluster').len().sort('len',descending=True)
s = sizes['len'].to_numpy()
log(f"  n nodes={cl.height} unique names={cl['name'].n_unique()} n clusters={len(s)}")
log(f"  cluster sizes: top20={s[:20].tolist()}")
for thr in [1,2,3,5,10,11,12,20,50]:
    log(f"  clusters size>={thr}: {(s>=thr).sum()} covering {s[s>=thr].sum()} nodes")

# ---------- filtered network ----------
fd = pl.read_csv(FILTERED, infer_schema_length=10000)
log(f"FILTERED {FILTERED.name}: shape={fd.shape}")
log(f"  dtypes={dict(zip(fd.columns,[str(d) for d in fd.dtypes]))}")
log(f"  Score min={fd['Score'].min():.8f} max={fd['Score'].max():.8f}")
log(f"  n negative Score = {(fd['Score']<0).sum()}")
log(f"  qvalue min={fd['qvalue'].min()} max={fd['qvalue'].max()}")
nodes_fd = set(fd['nodeA'].to_list())|set(fd['nodeB'].to_list())
log(f"  unique nodes in filtered network = {len(nodes_fd)}")

ed = pl.read_csv(EDGES, infer_schema_length=10000)
log(f"EDGES {EDGES.name}: shape={ed.shape} extra cols={[c for c in ed.columns if c not in fd.columns]}")
nodes_ed = set(ed['nodeA'].to_list())|set(ed['nodeB'].to_list())
log(f"  unique nodes in edges file = {len(nodes_ed)}")
log(f"  edges==filtered on shared cols: {fd.select(['nodeA','nodeB','Score']).equals(ed.select(['nodeA','nodeB','Score']))}")

nodes_cl = set(cl['name'].to_list())
log(f"  nodes in clusters file = {len(nodes_cl)}")
log(f"  network nodes NOT in clusters file = {len(nodes_fd-nodes_cl)}")
log(f"  clusters-file nodes NOT in network = {len(nodes_cl-nodes_fd)}")

# cluster assignment consistency between edges file and clusters file
m = dict(zip(cl['name'].to_list(), cl['LeidenCluster'].to_list()))
mismatch = sum(1 for a,ca,b,cb in zip(ed['nodeA'],ed['LeidenClusterA'],ed['nodeB'],ed['LeidenClusterB'])
               if m.get(a)!=ca or m.get(b)!=cb)
log(f"  edges-file cluster labels disagreeing with clusters file: {mismatch}")

# ---------- N runs ----------
nr = pl.read_csv(NRUNS)
log(f"NRUNS {NRUNS.name}: shape={nr.shape} cols={nr.columns}")
log(f"  dtypes={dict(zip(nr.columns,[str(d) for d in nr.dtypes]))}")
log(f"  head:\n{nr.head(3)}")

# ---------- raw association list: header + a few lines, streamed ----------
with gzip.open(RAW,'rt') as fh:
    head = [next(fh) for _ in range(4)]
log(f"RAW {RAW.name} first lines:")
for h in head: log("   "+h.rstrip()[:300])
