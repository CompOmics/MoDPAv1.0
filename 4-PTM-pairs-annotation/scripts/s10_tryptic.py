"""
Step 5, part 2: the shared-measurement flag.

Reads:  work/ptm_events_nruns.parquet, work/sprot_human.jsonl.gz, FILTERED.
Writes: work/ev_tryptic.parquet (every same-protein PTM event pair, with the position gap and
        whether the two sites can fall on one tryptic peptide).
Feeds:  report section 7.2.

Two sites that can fall on one tryptic peptide are quantified from overlapping PSM sets, so their
correlation may reflect shared measurement rather than coordination. This is computed as a flag and
is never used as a filter anywhere in the analysis: a sequence rule cannot separate a
co-quantification artefact from a genuine phosphosite cluster, and 61 of the 85 pairs that UniProt
documents as explicit crosstalk fall on one peptide.

Search settings taken from the manuscript: trypsin, up to two missed cleavages, no cleavage before
proline. Two variants are computed, the standard rule and an upper bound in which a modified lysine
or arginine is treated as non-cleavable, since a side-chain modification blocks trypsin.

The script also decodes the flags shipped inside the filtered-network file. Their shared_peptide is
reproduced exactly by "same protein and position gap <= 5", a fixed proximity rule rather than a
tryptic calculation, and their potential_artefact is exactly "shared_peptide and same modification".
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import gzip, json, pandas as pd, numpy as np, itertools, collections
from bisect import bisect_left, bisect_right

log("=== STEP 5b: tryptic co-peptide control ===")
ev = pd.read_parquet(WORK/"ptm_events_nruns.parquet")
accs = set(ev.acc); seqs = {}
with gzip.open(WORK/"sprot_human.jsonl.gz","rt",encoding="utf-8") as fh:
    for ln in fh:
        r = json.loads(ln)
        if r["acc"][0] in accs: seqs[r["acc"][0]] = r["seq"]
log(f"  sequences loaded for {len(seqs):,}/{len(accs):,} MoDPA proteins")

def cuts_of(seq, block=()):
    """Sorted 1-based cleavage positions for trypsin, bracketed by 0 and len(seq).

    A cut is placed after every K or R that is not followed by proline. Positions listed in
    `block` are treated as non-cleavable, which models a modified lysine or arginine whose side
    chain trypsin cannot process. The sentinels at 0 and len(seq) let the peptide that starts at
    residue 1, and the one that ends at the C-terminus, be addressed the same way as any other.
    """
    n=len(seq); c=[0]
    for i,aa in enumerate(seq):
        if aa in "KR" and not (i+1<n and seq[i+1]=="P") and (i+1) not in block: c.append(i+1)
    if c[-1]!=n: c.append(n)
    return c

def shares(c, p1, p2, max_mc=2):
    """True if some tryptic peptide with <=max_mc missed cleavages covers both p1<=p2.

    Rather than enumerating peptides, this counts cleavage sites. A peptide running from cut index
    i to cut index j spans residues c[i]+1 .. c[j] and contains j-i-1 internal cleavage sites, that
    is j-i-1 missed cleavages. The shortest peptide containing both positions starts at the last
    cut before p1 and ends at the first cut at or after p2, so those two indices fix the minimum
    missed-cleavage count and the test is j-i-1 <= max_mc.

    Enumerating peptides instead is O(peptides) per pair and made the whole script take minutes on
    proteins carrying many PTM events; this is two binary searches.
    """
    i = bisect_right(c, p1-1) - 1      # last cut strictly before p1
    j = bisect_left(c, p2)             # first cut at or after p2
    return (j - i) <= max_mc + 1

rows=[]; by_prot=collections.defaultdict(list)
for t in ev.itertuples(): by_prot[t.acc].append((t.pos, t.name))
for acc,lst in by_prot.items():
    if len(lst)<2: continue
    s = seqs.get(acc)
    if not s: continue
    modpos = {p for p,_ in lst}
    c_std = cuts_of(s); c_blk = cuts_of(s, modpos)
    for (p1,n1),(p2,n2) in itertools.combinations(sorted(lst), 2):
        a,b = (n1,n2) if n1<n2 else (n2,n1)
        rows.append((a,b,acc,abs(p2-p1), shares(c_std,p1,p2), shares(c_blk,p1,p2)))
T = pd.DataFrame(rows, columns=["nodeA","nodeB","acc","position_gap","same_peptide_std","same_peptide_modblocked"])
T.to_parquet(WORK/"ev_tryptic.parquet", index=False)
log(f"  same-protein PTM-event pairs (all possible, not only significant): {len(T):,} over {T.acc.nunique():,} proteins")
log(f"  can share a tryptic peptide (<=2 missed cleavages, standard rule):  {int(T.same_peptide_std.sum()):,}")
log(f"  same, treating a modified K/R as non-cleavable (upper bound):       {int(T.same_peptide_modblocked.sum()):,}")
log(f"  position_gap percentiles 5/25/50/75/95: {np.percentile(T.position_gap,[5,25,50,75,95]).round(0).tolist()}")

fd = pl.read_csv(FILTERED, infer_schema_length=10000).to_pandas()
fd["k"]=fd.nodeA+"__"+fd.nodeB; T["k"]=T.nodeA+"__"+T.nodeB
mm = fd.merge(T[["k","same_peptide_std","same_peptide_modblocked","position_gap"]], on="k", how="left", suffixes=("","_mine"))
sp = mm[mm.same_protein]
log(f"  provided filtered file: same_protein={int(fd.same_protein.sum()):,} shared_peptide={int(fd.shared_peptide.sum()):,} potential_artefact={int(fd.potential_artefact.sum()):,}")
log(f"  my recomputation on those same-protein rows: std={int(sp.same_peptide_std.fillna(False).sum()):,} modblocked={int(sp.same_peptide_modblocked.fillna(False).sum()):,}")
log(f"  agreement with provided shared_peptide (standard rule):   {100*(sp.shared_peptide==sp.same_peptide_std.fillna(False)).mean():.2f}%")
log(f"  agreement with provided shared_peptide (mod-blocked rule): {100*(sp.shared_peptide==sp.same_peptide_modblocked.fillna(False)).mean():.2f}%")
log(f"  position_gap agreement: {100*(sp.position_gap==sp.position_gap_mine).mean():.2f}%")
log(f"  rows where provided shared_peptide=True but my standard rule=False: {int((sp.shared_peptide & ~sp.same_peptide_std.fillna(False)).sum()):,}")
log(f"  rows where provided shared_peptide=False but my standard rule=True: {int((~sp.shared_peptide & sp.same_peptide_std.fillna(False)).sum()):,}")
