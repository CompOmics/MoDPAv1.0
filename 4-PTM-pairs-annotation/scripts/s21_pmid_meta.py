"""
Step 7, part 2: classify every PubMed id in the shortlist.

Reads:  work/sprot_human.jsonl.gz, work/ev_htp_pmids.json, work/shortlist.parquet.
Writes: work/pmid_meta.json (PubMed id to title and journal), updates work/shortlist.parquet,
        rewrites the delivered shortlist TSV.
Feeds:  report section 10.

Each PubMed id is labelled large_scale_MS, review or primary_or_other. Titles and journals come
from the RT and RL lines of the same Swiss-Prot entries, so the classification needs no network
access and stays traceable to the flat file. Reviews are detected by journal name and by title
phrasing; large-scale MS is detected by ECO code, as described in s09.

No web or PubMed tool was used anywhere in this analysis. Every identifier in the shortlist comes
verbatim from a UniProt evidence string, and the shortlist is therefore UniProt-derived and awaits
manual verification against the primary literature.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import gzip, json, re, pandas as pd

log("=== STEP 7b: PMID METADATA FROM THE SWISS-PROT REFERENCE LISTS ===")
meta = {}
with gzip.open(WORK/"sprot_human.jsonl.gz","rt",encoding="utf-8") as fh:
    for ln in fh:
        r = json.loads(ln)
        for ref in r["refs"]:
            t = re.sub(r"\s+"," ",ref["title"]).strip().strip('";')
            j = re.sub(r"\s+"," ",ref["loc"]).strip()
            for pm in ref["pmids"]:
                if pm not in meta or (not meta[pm][0] and t): meta[pm] = (t, j)
json.dump(meta, open(WORK/"pmid_meta.json","w"))
log(f"  PMIDs with a title from the human Swiss-Prot reference lists: {len(meta):,}")

htp = set(json.load(open(WORK/"ev_htp_pmids.json"))["counts"])
REVIEW_J = re.compile(r"Nat\. Rev\.|Trends |Curr\. Opin\.|Annu\. Rev\.|Physiol\. Rev\.|Pharmacol\. Rev\.|Chem\. Rev\.", re.I)
REVIEW_T = re.compile(r"\breview\b|\ban overview\b|\bperspectives?\b:", re.I)
def classify(pm):
    t,j = meta.get(pm, ("",""))
    if not t and not j: return "unknown_not_in_sprot_refs"
    if pm in htp: return "large_scale_MS"
    if REVIEW_J.search(j) or REVIEW_T.search(t): return "review"
    return "primary_or_other"

sh = pd.read_parquet(WORK/"shortlist.parquet")
def summarise(s):
    pms = [p for p in (s.split(";") if s else []) if p]
    if not pms: return "", "", ""
    cl = [classify(p) for p in pms]
    titles = " || ".join(f"PubMed:{p}: {meta.get(p,('',''))[0][:110]}" for p in pms[:4])
    return ";".join(f"{p}:{c}" for p,c in zip(pms,cl)), \
           ("yes" if all(c=="large_scale_MS" for c in cl) else "no"), titles
out = [summarise(s) for s in sh.pmids]
sh["pmid_classes"]          = [o[0] for o in out]
sh["only_large_scale_MS"]   = [o[1] for o in out]
sh["pmid_titles_first4"]    = [o[2] for o in out]
sh.to_parquet(WORK/"shortlist.parquet", index=False)
import collections
c = collections.Counter(x.split(":")[1] for s in sh.pmid_classes if s for x in s.split(";"))
log(f"  PMID classes across the shortlist: {dict(c)}")
log(f"  shortlist rows whose every PMID is a large-scale MS study: {(sh.only_large_scale_MS=='yes').sum()}")
log(f"  shortlist rows with no PMID at all: {int((sh.pmids=='').sum())}")
log(f"  shortlist rows citing at least one review: {sum('review' in s for s in sh.pmid_classes)}")

cols = ["nodeA","nodeB","geneA","geneB","modA_name","modB_name","Score","qvalue","stratum","sign","tier",
        "evidence_class","evidence_text","pmids","pmid_classes","pmid_titles_first4","only_large_scale_MS",
        "weak_eco_only","ptm_comment_names_partner","same_protein","position_gap","same_peptide",
        "captured_by_network","captured_by_cluster","clusterA","clusterB"]
o = sh[cols].rename(columns={"Score":"SDC","same_peptide":"ctrl_shared_tryptic_peptide",
                             "weak_eco_only":"weak_eco_only_on_one_site"})
o = o.sort_values(["tier","stratum","sign","SDC"], ascending=[True,True,True,False])
p = BASE/f"{DATE}-modpa-shortlist.tsv"; o.to_csv(p, sep="\t", index=False)
log(f"  rewrote {p.name} with PMID classification ({len(o)} rows)")
