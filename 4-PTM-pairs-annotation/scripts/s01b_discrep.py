"""
Step 1, part 2: follow-up checks on the discrepancies that s01 surfaced.

Reads:  CLUSTERS, EDGES, FILTERED.
Writes: nothing; output goes to work/runlog.md.

Checks whether the two Leiden partitions are a refinement of one another, whether the edges-file
labelling is self-consistent per node, whether Score equals distance, and decodes the control
flags shipped inside the filtered-network file.

Caveat on the "cluster labels disagreeing" count this script prints: cluster identifiers are
arbitrary across independent clustering runs, so comparing them directly measures nothing. The
meaningful comparison is an adjusted Rand index, which is not computed here. That line of the run
log should be ignored.

REQUIRES the withdrawn EDGES file and will fail with FileNotFoundError until it is restored.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import numpy as np
cl = pl.read_csv(CLUSTERS); ed = pl.read_csv(EDGES, infer_schema_length=10000)
fd = pl.read_csv(FILTERED, infer_schema_length=10000)

log("--- cluster label discrepancy ---")
ea = np.unique(ed['LeidenClusterA'].to_numpy()); eb=np.unique(ed['LeidenClusterB'].to_numpy())
log(f"edges LeidenClusterA range {ea.min()}..{ea.max()} n={len(ea)}; B range {eb.min()}..{eb.max()} n={len(eb)}")
allE = np.unique(np.concatenate([ed['LeidenClusterA'].to_numpy(), ed['LeidenClusterB'].to_numpy()]))
log(f"edges: distinct cluster ids used = {len(allE)} -> {allE.tolist()}")
log(f"clusters-file: distinct = {cl['LeidenCluster'].n_unique()}")

# is edges-file labelling per-node consistent (same node always same label)?
lab = {}
inconsistent = 0
for a,ca in zip(ed['nodeA'],ed['LeidenClusterA']):
    if a in lab and lab[a]!=ca: inconsistent+=1
    lab.setdefault(a,ca)
for b,cb in zip(ed['nodeB'],ed['LeidenClusterB']):
    if b in lab and lab[b]!=cb: inconsistent+=1
    lab.setdefault(b,cb)
log(f"edges-file per-node label inconsistencies: {inconsistent}  (n labelled nodes={len(lab)})")

# contingency: does edges labelling refine or coarsen clusters-file labelling?
m = dict(zip(cl['name'].to_list(), cl['LeidenCluster'].to_list()))
import collections
pairs = collections.Counter((lab[n], m[n]) for n in lab)
byE = collections.defaultdict(set); byC = collections.defaultdict(set)
for (e,c) in pairs: byE[e].add(c); byC[c].add(e)
log(f"edges-cluster -> how many clusters-file ids: max={max(len(v) for v in byE.values())}")
log(f"clusters-file id -> how many edges-cluster ids: max={max(len(v) for v in byC.values())}")
sz_e = collections.Counter(lab.values())
log(f"edges-file cluster sizes (sorted desc): {sorted(sz_e.values(),reverse=True)}")

log("--- qvalue / pval in filtered ---")
log(f"filtered pval unique (first 10): {fd['pval'].unique().head(10).to_list()}")
log(f"filtered qvalue unique: {fd['qvalue'].unique().to_list()}")
log("--- row-order difference filtered vs edges ---")
log(f"same nodeA order: {(fd['nodeA']==ed['nodeA']).all()}")
log(f"pair_key sets equal: {set(fd['pair_key'])==set(ed['pair_key'])}")
log("--- sign / score in filtered ---")
log(f"Score==distance for all rows: {(fd['Score']==fd['distance']).all()}")
log(f"PCC<0 rows in filtered: {(fd['PCC']<0).sum()}")
log("--- precomputed control flags present in provided files ---")
for c in ['same_protein','same_site','same_mod','shared_peptide','potential_artefact']:
    log(f"  {c}: True={fd[c].sum()} / {fd.height}")
log(f"  position_gap non-null: {fd['position_gap'].is_not_null().sum()}")
