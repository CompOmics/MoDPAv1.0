"""
Step 4, part 2: pair-level evidence tables.

Reads:  work/node_ann.parquet, work/ev_ptmcc.parquet, work/ev_interact.json,
        work/ev_complexes.json, work/ev_subunit.json, work/ev_genes.json.
Writes: work/pair_t1a.parquet - Tier 1, site-level explicit crosstalk.
        work/pair_t1b.parquet - protein-level crosstalk statements, deliberately NOT Tier 1.
        work/pair_t2.parquet  - Tier 2, both sites attributed to the same enzyme.
        work/pair_t4.parquet  - Tier 4, protein pairs that interact or share a complex.
Feeds:  report section 5, and the tier assignment in s13.

This is type (b) evidence: that the two sites of a pair are linked. Tier 3 and Tier 5 need no table
here because they follow from the node-level columns in s11.

Tier 1a requires both positions of the pair to be cited inside one CC PTM block on the same
protein, which is genuine site-level evidence and yields 114 PTM event pairs.

Tier 1b is the weaker case in which a CC PTM block cites one site of the pair and names the gene
symbol of the partner protein, without naming the partner site. That is a statement about two
proteins, not two sites, so it fails the Tier 1 definition and is carried as a separate boolean
column (ptm_comment_names_partner) rather than folded into the tier ladder. An earlier version of
this analysis did fold it in, which inflated Tier 1 from 85 to 851 pairs and filled the shortlist
with pairs whose evidence did not actually concern the sites in question.

Tier 4 draws on three sources of decreasing reliability: CC INTERACTION lines, which name partner
accessions explicitly; shared CORUM or ComplexPortal identifiers; and gene symbols matched inside
CC SUBUNIT free text. The last can produce false positives, which is why Tier 4 is described in the
report as the weakest type (b) tier. Complexes with more than 200 MoDPA members are skipped, since
they would contribute a quadratic number of uninformative pairs.

STOP is a list of uppercase tokens that look like gene symbols in free text but are not, such as
amino-acid abbreviations and common words. Without it, matching gene symbols against prose produces
large numbers of spurious protein pairs.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import pandas as pd, numpy as np, json, re, itertools, collections

log("=== STEP 4b: PAIR-LEVEL EVIDENCE TABLES ===")
A  = pd.read_parquet(WORK/"node_ann.parquet")
P  = pd.read_parquet(WORK/"ev_ptmcc.parquet")
interact  = json.load(open(WORK/"ev_interact.json"))
complexes = json.load(open(WORK/"ev_complexes.json"))
subunit   = json.load(open(WORK/"ev_subunit.json"))
genes     = json.load(open(WORK/"ev_genes.json"))

accs = set(A.acc)
gene2acc = collections.defaultdict(set)
for a,g in genes.items():
    if g: gene2acc[g.upper()].add(a)
by_acc = collections.defaultdict(list)
for r in A.itertuples(): by_acc[r.acc].append(r)

STOP = set("""THE AND NOT FOR WITH FROM THIS THAT ARE MAY CAN ONE TWO SET MAX MIN PTM MS NMR ATP GTP DNA RNA
CTER NTER ECO PDB SDS TLC HPLC ALSO BUT WHEN ONLY HAS HAVE BEEN VIA ITS NON PRO ALA ARG ASN ASP CYS GLN GLU
GLY HIS ILE LEU LYS MET PHE SER THR TRP TYR VAL CELL TYPE FULL SITE SITES FORM FORMS ROLE HIGH LOW LIKE
INTO SUCH BOTH ALL ANY OUT OFF PER END ADP AMP CAMP""".split())

# ---------- Tier 1a: same protein, both positions cited in one CC PTM block ----------
# A block that cites two or more residue positions is making a statement about those positions
# together, which is what "explicit crosstalk" means here. Only MoDPA events whose position is
# actually named in the block qualify, so a block citing Ser-123 and Lys-456 links exactly those
# two sites and not every other annotated site on the protein. Blocks citing fewer than two
# positions cannot link anything and are skipped.
t1a = []
for r in P.itertuples():
    pos = {int(t.split("-")[1]) for t in (r.positions.split(";") if r.positions else []) if t}
    if len(pos) < 2: continue
    ev_here = [e for e in by_acc.get(r.acc, []) if e.pos in pos]
    # Sort by name before combining so the pair is emitted in the same orientation every time and
    # can be joined against the association list, which is itself stored upper-triangular.
    for a,b in itertools.combinations(sorted(ev_here, key=lambda x:x.name), 2):
        t1a.append((a.name, b.name, r.acc, r.text[:800], r.pmids))
T1A = pd.DataFrame(t1a, columns=["nodeA","nodeB","acc","ev_text","ev_pmids"]).drop_duplicates(["nodeA","nodeB"])
log(f"  Tier 1a (same protein, both sites cited in one CC PTM block): {len(T1A):,} PTM-event pairs")

# ---------- Tier 1b: cross-protein, block cites this site and names the partner gene ----------
t1b = []
for r in P.itertuples():
    pos = {int(t.split("-")[1]) for t in (r.positions.split(";") if r.positions else []) if t}
    if not pos: continue
    toks = {t for t in re.findall(r"\b[A-Z][A-Z0-9-]{2,}\b", r.text)} - STOP
    partners = {q for t in toks for q in gene2acc.get(t, ()) if q != r.acc}
    if not partners: continue
    src = [e for e in by_acc.get(r.acc, []) if e.pos in pos]
    for s in src:
        for q in partners:
            for e2 in by_acc.get(q, []):
                a,b = sorted([s.name, e2.name])
                t1b.append((a, b, r.acc, r.text[:800], r.pmids))
T1B = pd.DataFrame(t1b, columns=["nodeA","nodeB","src_acc","ev_text","ev_pmids"]).drop_duplicates(["nodeA","nodeB"])
log(f"  Tier 1b (cross-protein, CC PTM block cites the site and names the partner gene): {len(T1B):,}")

# ---------- Tier 2: shared named enzyme ----------
enzmap = collections.defaultdict(list)
for r in A.itertuples():
    for e in (r.enzymes.split(";") if r.enzymes else []):
        if e: enzmap[e].append(r.name)
t2 = {}
for e, names in enzmap.items():
    if len(names) < 2: continue
    for a,b in itertools.combinations(sorted(names), 2):
        t2.setdefault((a,b), set()).add(e)
T2 = pd.DataFrame([(a,b,";".join(sorted(v))) for (a,b),v in t2.items()],
                  columns=["nodeA","nodeB","shared_enzymes"])
log(f"  Tier 2 (both sites annotated '; by <same enzyme>'): {len(T2):,} pairs, "
    f"{len([e for e,n in enzmap.items() if len(n)>1])} enzymes with >=2 annotated MoDPA sites")

# ---------- Tier 4: protein-pair interaction / complex ----------
pp = collections.defaultdict(set)
for p, parts in interact.items():
    for q in parts:
        if q in accs and q != p:
            pp[tuple(sorted((p,q)))].add("CC_INTERACTION")
cx = collections.defaultdict(set)
for p, ids in complexes.items():
    for c in ids: cx[c].add(p)
for c, mem in cx.items():
    mem = sorted(m for m in mem if m in accs)
    if len(mem) < 2 or len(mem) > 200: continue
    for a,b in itertools.combinations(mem, 2): pp[(a,b)].add("COMPLEX")
gset = {g.upper():a for a,g in genes.items() if g}
for p, txt in subunit.items():
    if p not in accs: continue
    toks = {t for t in re.findall(r"\b[A-Z][A-Z0-9-]{2,}\b", txt)} - STOP
    for t in toks:
        for q in gene2acc.get(t, ()):
            if q in accs and q != p: pp[tuple(sorted((p,q)))].add("SUBUNIT_TEXT")
T4 = pd.DataFrame([(a,b,";".join(sorted(v))) for (a,b),v in pp.items()],
                  columns=["accA","accB","pp_evidence"])
log(f"  Tier 4 (protein pairs with interaction / complex / subunit-text evidence): {len(T4):,} protein pairs")
log(f"    by evidence source: " + str(collections.Counter(x for v in pp.values() for x in v)))

for df,nm in [(T1A,"t1a"),(T1B,"t1b"),(T2,"t2"),(T4,"t4")]:
    df.to_parquet(WORK/f"pair_{nm}.parquet", index=False)
log("  -> work/pair_t1a.parquet, pair_t1b.parquet, pair_t2.parquet, pair_t4.parquet")
