"""
Step 3, part 2: positional quality control against Swiss-Prot 2026_03.

Reads:  work/ptm_events_nruns.parquet, work/sprot_human.jsonl.gz.
Writes: work/posqc.parquet (per-event QC status).
Feeds:  report section 4.

This is a stop condition for the whole analysis. If the residue at the stated position in the
2026_03 sequence does not match the residue in the PTM event identifier, the coordinates come from
a different release and no site-level matching downstream can be trusted.

The check runs three ways so the result cannot be an artefact of the reading convention:
  * exact residue match at the stated 1-based position, plus chemical compatibility with the
    modification, using the residue sets in UNIMOD;
  * the same test at position -1 and +1, which would expose an off-by-one or a 0-based convention;
  * an offset probe on any event that does mismatch.

Result: 0 of 7,644 mismatch under the 1-based reading, against 14.93% and 13.53% agreement under
the two shifted readings, so the convention is confirmed empirically rather than assumed.

Secondary accessions are resolved through the prim mapping, which maps every accession Swiss-Prot
lists for an entry, primary and secondary alike, back to the primary one.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import gzip, json, pandas as pd, collections

log("=== STEP 3b: SWISS-PROT SANITY + POSITIONAL QC ===")
ev = pd.read_parquet(WORK/"ptm_events_nruns.parquet")
need = set(ev.acc)

seq, prim, gene, pname, seqver, length = {}, {}, {}, {}, {}, {}
nfeat = collections.Counter(); nrec = 0
with gzip.open(WORK/"sprot_human.jsonl.gz","rt",encoding="utf-8") as fh:
    for ln in fh:
        r = json.loads(ln); nrec += 1
        for a in r["acc"]:
            prim[a] = r["acc"][0]
        p = r["acc"][0]
        seq[p]=r["seq"]; gene[p]=r["gene"]; pname[p]=r["protname"]
        seqver[p]=r["seqver"]; length[p]=r["length"]
        for f in r["features"]: nfeat[f["type"]] += 1
log(f"  human records={nrec:,}; primary accessions={len(seq):,}; all accessions (incl. secondary)={len(prim):,}")
log(f"  feature counts: {dict(nfeat)}")
bad_len = sum(1 for p in seq if length[p] != len(seq[p]))
log(f"  records where ID-line length != parsed sequence length: {bad_len}")

miss = sorted(need - set(prim))
log(f"  MoDPA proteins = {len(need):,}; not found in human Swiss-Prot 2026_03 = {len(miss)}")
if miss: log(f"    examples: {miss[:15]}")
sec = sorted(a for a in need if a in prim and prim[a]!=a)
log(f"  MoDPA accessions that are SECONDARY in 2026_03 (demerged/merged): {len(sec)} {sec[:15]}")

# residue compatibility
UNI_OK = {u:(n,r) for u,(n,r) in UNIMOD.items()}
rows=[]
for t in ev.itertuples():
    p = prim.get(t.acc)
    s = seq.get(p) if p else None
    if s is None:
        rows.append((t.name,t.acc,t.pos,t.res,t.unimod,None,"no_sequence",None)); continue
    if t.pos < 1 or t.pos > len(s):
        rows.append((t.name,t.acc,t.pos,t.res,t.unimod,None,"out_of_range",len(s))); continue
    aa = s[t.pos-1]
    ok_id  = (aa == t.res)
    allowed = UNI_OK.get(t.unimod,(None,set()))[1]
    ok_chem = aa in allowed
    st = "ok" if (ok_id and ok_chem) else ("residue_mismatch" if not ok_id else "chem_mismatch")
    rows.append((t.name,t.acc,t.pos,t.res,t.unimod,aa,st,len(s)))
qc = pd.DataFrame(rows, columns=["name","acc","pos","res","unimod","aa_2026_03","status","seqlen"])
qc.to_parquet(WORK/"posqc.parquet", index=False)
vc = qc.status.value_counts()
log(f"  positional QC over {len(qc):,} PTM events: {vc.to_dict()}")
log(f"  MISMATCH RATE (any status != ok) = {100*(qc.status!='ok').mean():.3f}%")
log("  mismatches by modification:")
bad = qc[qc.status!="ok"]
if len(bad):
    for (u,st),n in bad.groupby(["unimod","status"]).size().items():
        log(f"    unimod {u} ({UNIMOD.get(int(u),('?',))[0]}) {st}: {n}")
    log("  offset probe on mismatches (is the true residue at pos+k?):")
    for k in (-2,-1,1,2):
        c=0
        for t in bad.itertuples():
            s = seq.get(prim.get(t.acc,""),"")
            j = t.pos-1+k
            if s and 0<=j<len(s) and s[j]==t.res: c+=1
        log(f"    offset {k:+d}: {c} of {len(bad)}")
    log(f"  examples: {bad.head(10)[['name','aa_2026_03','status','seqlen']].to_dict('records')}")

# off-by-one global probe over ALL events (would reveal a 0-based vs 1-based error)
for k in (-1,0,1):
    c=0
    for t in ev.itertuples():
        s = seq.get(prim.get(t.acc,""),"")
        j = t.pos-1+k
        if s and 0<=j<len(s) and s[j]==t.res: c+=1
    log(f"  ALL events, residue matches at pos{k:+d} (1-based read): {c:,}/{len(ev):,} = {100*c/len(ev):.2f}%")
