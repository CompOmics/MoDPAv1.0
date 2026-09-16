"""
Step 3, part 3: build the offline evidence base from the parsed human entries.

Reads:  work/ptm_events_nruns.parquet, work/sprot_human.jsonl.gz.
Writes: work/ev_sites.parquet     - one row per annotated feature position on a MoDPA protein, with
                                    note text, ECO codes, PubMed ids, the "alternate" qualifier and
                                    the enzyme named in a "; by X" note.
        work/ev_ptmcc.parquet     - CC PTM comment blocks, with the residue positions they cite and
                                    their PubMed ids.
        work/ev_interact.json     - CC INTERACTION partners per protein.
        work/ev_complexes.json    - CORUM and ComplexPortal identifiers per protein.
        work/ev_subunit.json      - CC SUBUNIT free text per protein.
        work/ev_reactome.json     - Reactome cross-references per protein.
        work/ev_genes.json        - accession to gene symbol.
        work/ev_htp_pmids.json    - PubMed ids that Swiss-Prot attaches with a high-throughput ECO
                                    code, with titles and usage counts.
Feeds:  report sections 5 and 10.

Two ECO distinctions carry the circularity argument in the report and are recorded per feature:
  * ECO:0000269 is experimental evidence from a manual, usually dedicated, study.
  * ECO:0007744 and ECO:0000244 are combinatorial evidence, which in UniProt practice marks a site
    taken from a large-scale mass-spectrometry survey.
The second class is not independent of the data MoDPA was built from, so a pair supported only by
it cannot be treated as external validation.

The high-throughput PubMed set is defined from the data rather than by hand: a study counts as
large-scale MS if UniProt attaches it with one of those ECO codes anywhere in the human proteome.
That yields 2,390 studies, led by the well-known phosphoproteome and acetylome surveys.

A CC PTM block is the richest offline source of explicit crosstalk statements, so the positions
cited inside one are extracted with POSMENT, which matches the three-letter residue and position
form UniProt uses in prose, for example "Ser-123".
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import gzip, json, re, pandas as pd, collections, itertools

log("=== STEP 3c: EVIDENCE BASE TABLES ===")
ev = pd.read_parquet(WORK/"ptm_events_nruns.parquet")
need = set(ev.acc)

AA3 = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
POSMENT = re.compile(rf"\b({AA3})-(\d+)\b")
PMID_RE = re.compile(r"PubMed:(\d+)")
BY_RE   = re.compile(r";\s*by\s+([^;]+)")          # MOD_RES '; by CK2' / '; by PKA and PKC'
ECO_EXP = {"ECO:0000269"}                          # manual experimental
ECO_HTP = {"ECO:0007744","ECO:0000244"}            # combinatorial / large-scale MS
ECO_WEAK= {"ECO:0000250","ECO:0000255","ECO:0000305","ECO:0000312","ECO:0000303"}

# unimod -> regex over the Swiss-Prot feature note, for "same modification annotated"
SAMEMOD = {
 1:   re.compile(r"acetyl", re.I),
 7:   re.compile(r"citrullin|deimin", re.I),
 21:  re.compile(r"^Phospho(serine|threonine|tyrosine)", re.I),
 23:  re.compile(r"dehydro(alanine|butyrine)", re.I),
 34:  re.compile(r"\b(N6-methyllysine|Omega-N-methylarginine|methylarginine|methyllysine)\b", re.I),
 36:  re.compile(r"dimethyl", re.I),
 37:  re.compile(r"trimethyl", re.I),
 53:  re.compile(r"hydroxynonenal|HNE", re.I),
 64:  re.compile(r"succinyl", re.I),
 299: re.compile(r"carboxyglutamate|carboxy", re.I),
 535: re.compile(r"ubiquitin", re.I),
}

sites = []          # per (acc,pos) annotation
ptmcc = []          # CC PTM blocks
subunit, interact, complexes, reactome = {}, collections.defaultdict(set), collections.defaultdict(set), {}
seqs, genes, pnames = {}, {}, {}
ref_title = {}      # pmid -> title
htp_pmids = collections.Counter()
allec = collections.Counter()

with gzip.open(WORK/"sprot_human.jsonl.gz","rt",encoding="utf-8") as fh:
    for ln in fh:
        r = json.loads(ln); p = r["acc"][0]
        for ref in r["refs"]:
            t = re.sub(r'\s+'," ",ref["title"]).strip().strip('";')
            for pm in ref["pmids"]:
                if t: ref_title.setdefault(pm, t)
        if p not in need:
            for f in r["features"]:
                for e in f["eco"]:
                    if e in ECO_HTP: htp_pmids.update(PMID_RE.findall(f["evidence"]))
            continue
        seqs[p]=r["seq"]; genes[p]=r["gene"]; pnames[p]=r["protname"]
        for f in r["features"]:
            for e in f["eco"]:
                allec[e]+=1
                if e in ECO_HTP: htp_pmids.update(PMID_RE.findall(f["evidence"]))
            loc = f["loc"].split()[0] if f["loc"] else ""
            m = re.match(r"^(\d+)(?:\.\.(\d+))?$", loc)
            if not m: continue
            a = int(m.group(1)); b = int(m.group(2)) if m.group(2) else a
            eco = set(f["eco"])
            for pos in ({a,b} if f["type"]=="CROSSLNK" else range(a,b+1) if b-a<50 else {a}):
                sites.append(dict(acc=p, pos=pos, ftype=f["type"], note=f["note"],
                    eco=";".join(sorted(eco)), pmids=";".join(f["pmids"]),
                    exp=bool(eco & ECO_EXP), htp=bool(eco & ECO_HTP), weak=bool(eco & ECO_WEAK),
                    alternate=("alternate" in f["note"].lower()),
                    by=";".join(x.strip().rstrip(".") for x in BY_RE.findall(f["note"]))))
        for blk in r["cc"].get("PTM", []):
            ptmcc.append(dict(acc=p, text=blk,
                pmids=";".join(sorted(set(PMID_RE.findall(blk)))),
                positions=";".join(f"{aa}-{n}" for aa,n in POSMENT.findall(blk))))
        if r["cc"].get("SUBUNIT"): subunit[p] = " ".join(r["cc"]["SUBUNIT"])
        for line in r["cc"].get("INTERACTION", []):
            for m in re.finditer(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?[;:]", line):
                interact[p].add(m.group(1))
        for db in ("CORUM","ComplexPortal"):
            for line in r["dr"].get(db, []):
                complexes[p].add(line.split(";")[1].strip())
        reactome[p] = [l.split(";")[1].strip() for l in r["dr"].get("Reactome",[])]

S = pd.DataFrame(sites)
log(f"  annotated feature-positions on MoDPA proteins: {len(S):,} over {S.acc.nunique():,} proteins")
log(f"  by feature type: {S.ftype.value_counts().to_dict()}")
log(f"  ECO codes seen on those features: {dict(allec.most_common(12))}")
log(f"  'alternate' qualifier present on {int(S.alternate.sum()):,} feature-positions")
log(f"  features carrying an enzyme ('; by X'): {int((S.by!='').sum()):,}")
S.to_parquet(WORK/"ev_sites.parquet", index=False)

P = pd.DataFrame(ptmcc)
log(f"  CC PTM blocks on MoDPA proteins: {len(P):,} over {P.acc.nunique():,} proteins")
log(f"  CC PTM blocks citing >=2 explicit residue positions: {(P.positions.str.count(';')>=1).sum():,}")
P.to_parquet(WORK/"ev_ptmcc.parquet", index=False)

log(f"  proteins with SUBUNIT text: {len(subunit):,}; with CC INTERACTION partners: {len(interact):,}")
log(f"  proteins with CORUM/ComplexPortal xrefs: {sum(1 for v in complexes.values() if v):,}")
log(f"  proteins with Reactome xrefs: {sum(1 for v in reactome.values() if v):,}")
json.dump({k:sorted(v) for k,v in interact.items()},  open(WORK/"ev_interact.json","w"))
json.dump({k:sorted(v) for k,v in complexes.items()}, open(WORK/"ev_complexes.json","w"))
json.dump(subunit,  open(WORK/"ev_subunit.json","w"))
json.dump(reactome, open(WORK/"ev_reactome.json","w"))
json.dump({k:genes[k] for k in genes}, open(WORK/"ev_genes.json","w"))

# large-scale-MS PMIDs, defined data-driven: PMIDs that UniProt attaches with a
# combinatorial/high-throughput ECO code anywhere in the human proteome
htp = {pm:c for pm,c in htp_pmids.items()}
json.dump({"counts":htp, "titles":{pm:ref_title.get(pm,"") for pm in htp}},
          open(WORK/"ev_htp_pmids.json","w"))
log(f"  high-throughput-MS PMIDs (ECO:0007744/0000244 anywhere in human Swiss-Prot): {len(htp):,}")
top = sorted(htp.items(), key=lambda x:-x[1])[:8]
for pm,c in top: log(f"    PubMed:{pm} n={c:,}  {ref_title.get(pm,'')[:95]}")
