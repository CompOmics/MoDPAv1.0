"""
Step 4, part 1: per-event annotation table.

Reads:  work/ptm_events_nruns.parquet, work/ev_sites.parquet, work/ev_ptmcc.parquet,
        work/ev_genes.json.
Writes: work/node_ann.parquet (one row per PTM event, with its UniProt annotation status).
Feeds:  report section 5, and the tier assignment in s13.

This is type (a) evidence only: whether each individual site is a known modified residue. It is
kept strictly separate from type (b) evidence, that the pair is linked, which is built in s12.

Columns worth explaining:
  * ann_mod       - a MOD_RES, CROSSLNK, CARBOHYD or LIPID feature exists at this exact position.
  * ann_samemod   - that annotation is the same chemistry as the MoDPA modification, tested with
                    the SAMEMOD regexes below. Swiss-Prot describes modifications in prose, so
                    matching is by keyword: acetylation appears as "N6-acetyllysine",
                    citrullination as "Citrulline", ubiquitination as a CROSSLNK "Glycyl lysine
                    isopeptide ... in ubiquitin", and so on.
  * ann_exp       - carries ECO:0000269, manual experimental evidence.
  * ann_htp       - carries ECO:0007744 or ECO:0000244, large-scale MS evidence. See s09.
  * ann_weak      - only ECO:0000250 (by similarity), ECO:0000255 (sequence model) or ECO:0000305
                    (curator inference), with no experimental code at all.
  * ann_alternate - the note carries the "alternate" qualifier, which is how UniProt records that a
                    residue can carry more than one modification.
  * enzymes       - the enzyme named in a "; by X" note, normalised by norm_enz. Composite names
                    such as "PKB/AKT1" and lists such as "by PKA, PKC and PKB/AKT1" are split, so
                    two sites written differently still match on a shared enzyme in s12.
  * n_ann_on_protein - annotation density of the protein, used as a matching covariate in the s15
                    null so that enrichment cannot be explained by well-studied proteins alone.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import pandas as pd, numpy as np, json, re, collections

log("=== STEP 4a: PER-EVENT ANNOTATION TABLE ===")
ev = pd.read_parquet(WORK/"ptm_events_nruns.parquet")
S  = pd.read_parquet(WORK/"ev_sites.parquet")
P  = pd.read_parquet(WORK/"ev_ptmcc.parquet")
genes = json.load(open(WORK/"ev_genes.json"))

SAMEMOD = {
 1:   re.compile(r"acetyl", re.I),
 7:   re.compile(r"citrullin|deimin", re.I),
 21:  re.compile(r"Phospho(serine|threonine|tyrosine)", re.I),
 23:  re.compile(r"dehydro(alanine|butyrine)", re.I),
 34:  re.compile(r"(N6-methyllysine|Omega-N-methylarginine|methylarginine|methyllysine)", re.I),
 36:  re.compile(r"dimethyl", re.I),
 37:  re.compile(r"trimethyl", re.I),
 53:  re.compile(r"hydroxynonenal", re.I),
 64:  re.compile(r"succinyl", re.I),
 299: re.compile(r"carboxy", re.I),
 535: re.compile(r"ubiquitin", re.I),
}
MODFEAT = {"MOD_RES","CROSSLNK","CARBOHYD","LIPID"}
AA1 = dict(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",His="H",Ile="I",
           Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V")

def norm_enz(s):
    """Split a MOD_RES '; by X' note into a set of normalised enzyme names.

    UniProt writes these as free prose, so the same enzyme appears in several shapes:
    "by CK2", "by PKB/AKT1", "by PKA, PKC and PKB/AKT1", "by RPS6KA1 and RPS6KB1". Splitting on
    "and", commas, slashes and semicolons turns each into a set of names, so two sites written
    differently still intersect in s12 if they share a kinase. Slashes are split rather than kept
    because "PKB/AKT1" is one enzyme under two names, and a site written "by AKT1" should match it.

    Case is normalised and a leading article stripped. "autocatalysis" is dropped: it names no
    external regulator, so two autocatalytic sites on different proteins share nothing. The length
    cap discards fragments of prose that survive the split without being enzyme names.
    """
    if not s: return set()
    out=set()
    for tok in re.split(r"\s+and\s+|,|/|;", s):
        t = tok.strip().strip(".").upper()
        t = re.sub(r"^(THE|A)\s+","",t)
        if t and t not in {"","AUTOCATALYSIS"} and len(t)<=20: out.add(t)
    return out

sidx = collections.defaultdict(list)
for r in S.itertuples(): sidx[(r.acc, r.pos)].append(r)

# CC PTM blocks: positions cited, per protein
pcc = collections.defaultdict(list)     # acc -> list of (block_id, set_of_positions, text, pmids)
for i,r in enumerate(P.itertuples()):
    pos = set()
    for tok in (r.positions.split(";") if r.positions else []):
        if not tok: continue
        aa3, n = tok.split("-"); pos.add(int(n))
    pcc[r.acc].append((i, pos, r.text, r.pmids))

rows=[]
for t in ev.itertuples():
    feats = sidx.get((t.acc, t.pos), [])
    modf  = [f for f in feats if f.ftype in MODFEAT]
    same  = [f for f in modf if SAMEMOD[t.unimod].search(f.note or "")]
    enz   = set()
    for f in modf: enz |= norm_enz(f.by)
    blocks = [b for b in pcc.get(t.acc, []) if t.pos in b[1]]
    rows.append(dict(
        name=t.name, acc=t.acc, pos=t.pos, res=t.res, unimod=t.unimod,
        gene=genes.get(t.acc), n_runs=t.n_runs,
        ann_any   = len(feats)>0,
        ann_mod   = len(modf)>0,
        ann_samemod = len(same)>0,
        ann_exp   = any(f.exp for f in modf),
        ann_htp   = any(f.htp for f in modf),
        ann_weak  = any(f.weak for f in modf) and not any(f.exp or f.htp for f in modf),
        ann_alternate = any(f.alternate for f in modf),
        enzymes   = ";".join(sorted(enz)),
        ann_note  = " | ".join(f"{f.ftype}:{f.note}" for f in modf)[:500],
        ann_pmids = ";".join(sorted({p for f in modf for p in (f.pmids.split(";") if f.pmids else []) if p})),
        ptmcc_blocks = ";".join(str(b[0]) for b in blocks),
        n_ann_on_protein = 0,
    ))
A = pd.DataFrame(rows)
nann = S[S.ftype.isin(MODFEAT)].groupby("acc").size()
A["n_ann_on_protein"] = A.acc.map(nann).fillna(0).astype(int)
A.to_parquet(WORK/"node_ann.parquet", index=False)

log(f"  events={len(A):,}")
for c in ["ann_any","ann_mod","ann_samemod","ann_exp","ann_htp","ann_weak","ann_alternate"]:
    log(f"    {c:<15} {int(A[c].sum()):>6,}  ({100*A[c].mean():5.1f}%)")
log(f"    has enzyme annotation      {int((A.enzymes!='').sum()):>6,}")
log(f"    cited in a CC PTM block    {int((A.ptmcc_blocks!='').sum()):>6,}")
log("  ann_samemod by modification:")
for u,g in A.groupby("unimod"):
    log(f"    {u:>4} {UNIMOD[u][0]:<13} n={len(g):>5}  ann_mod={int(g.ann_mod.sum()):>5}  ann_samemod={int(g.ann_samemod.sum()):>5}  ann_exp={int(g.ann_exp.sum()):>5}  ann_htp={int(g.ann_htp.sum()):>5}")
log(f"  distinct enzymes named: {len({e for s in A.enzymes if s for e in s.split(';')})}")
