"""
Step 3, part 1: stream Swiss-Prot and extract the human entries.

Reads:  SPROT (uniprot_2026_03_sprot.dat.gz, 668 MB gzipped, all species).
Writes: work/sprot_human.jsonl.gz (one JSON record per human entry, about 33 MB).
Feeds:  report section 5.

The flat file is never parsed into a single object. Lines are accumulated into a per-entry buffer,
the buffer is tested for OX NCBI_TaxID=9606, and only human entries are parsed and written out.
575,748 entries are scanned and 20,431 kept, in roughly 100 s.

Per entry the parser keeps: accessions (the first is primary, the rest secondary), gene names,
protein name, sequence and length, keywords, the reference list (PubMed ids, title, journal), the
CC comment blocks listed in CC_KEEP, the cross-references listed in DR_KEEP, and the sequence
features listed in FEAT_KEEP with their ECO evidence codes and PubMed ids.

Two flat-file details drive the parsing logic and are easy to get wrong:
  * FT records span several lines. A new feature begins when column 6 of an FT line is not a space;
    continuation lines carry qualifiers such as /note= and /evidence=, and a qualifier value can
    itself wrap onto the following line, which is what the final else branch handles by appending
    to whichever qualifier was opened last.
  * PubMed ids appear in two syntaxes. RX reference lines use "PubMed=12345678" while FT and CC
    evidence strings use "PubMed:12345678". Both regexes are needed; using only the colon form
    silently produces empty reference lists with no error.
"""
import sys; sys.path.insert(0,"scripts")
from _common import *
import gzip, json, re, time

OUT = WORK/"sprot_human.jsonl.gz"
FEAT_KEEP = {"MOD_RES","CROSSLNK","CARBOHYD","LIPID","BINDING","SITE","ACT_SITE","DISULFID"}
CC_KEEP   = {"PTM","SUBUNIT","INTERACTION","FUNCTION","SUBCELLULAR LOCATION","DOMAIN","MISCELLANEOUS"}
DR_KEEP   = {"Reactome","CORUM","STRING","IntAct","ComplexPortal","BioGRID","GlyConnect","iPTMnet","PhosphoSitePlus"}
PMID_RE   = re.compile(r"PubMed:(\d+)")
RX_PMID   = re.compile(r"PubMed=(\d+)")
ECO_RE    = re.compile(r"ECO:\d{7}")

def parse_entry(lines):
    rec = {"acc": [], "id": None, "gene": None, "genes": [], "protname": None,
           "seq": None, "seqver": None, "length": None, "taxid": None,
           "kw": [], "features": [], "cc": {}, "dr": {}, "refs": []}
    cc_key, cc_buf = None, []
    ft_cur = None
    seq_lines, in_seq = [], False
    ref_cur = None
    for ln in lines:
        tag, body = ln[:2], ln[5:].rstrip("\n")
        if in_seq:
            if ln.startswith("//"): break
            seq_lines.append(ln.strip().replace(" ", "")); continue
        if tag == "ID":
            rec["id"] = body.split()[0]
            m = re.search(r"(\d+)\s+AA", body)
            if m: rec["length"] = int(m.group(1))
        elif tag == "AC":
            rec["acc"] += [a.strip() for a in body.split(";") if a.strip()]
        elif tag == "DT":
            m = re.search(r"sequence version (\d+)", body)
            if m: rec["seqver"] = int(m.group(1))
        elif tag == "DE":
            if rec["protname"] is None and "Full=" in body:
                rec["protname"] = body.split("Full=",1)[1].split("{")[0].rstrip("; ").strip()
        elif tag == "GN":
            for m in re.finditer(r"Name=([^;{ ]+)", body):
                rec["genes"].append(m.group(1))
            for m in re.finditer(r"Synonyms=([^;]+)", body):
                rec["genes"] += [s.split("{")[0].strip() for s in m.group(1).split(",")]
        elif tag == "OX":
            m = re.search(r"NCBI_TaxID=(\d+)", body)
            if m: rec["taxid"] = int(m.group(1))
        elif tag == "KW":
            rec["kw"] += [k.strip().rstrip(".") for k in body.split(";") if k.strip()]
        elif tag == "RN":
            if ref_cur: rec["refs"].append(ref_cur)
            ref_cur = {"pmids": [], "title": "", "loc": "", "cmt": ""}
        elif tag == "RX" and ref_cur is not None:
            ref_cur["pmids"] += RX_PMID.findall(body) + PMID_RE.findall(body)
        elif tag == "RT" and ref_cur is not None:
            ref_cur["title"] += " " + body.strip().strip('";')
        elif tag == "RL" and ref_cur is not None:
            ref_cur["loc"] += " " + body
        elif tag == "RC" and ref_cur is not None:
            ref_cur["cmt"] += " " + body
        elif tag == "CC":
            if body.startswith("-!- "):
                if cc_key: rec["cc"].setdefault(cc_key, []).append(" ".join(cc_buf).strip())
                head = body[4:]
                cc_key = head.split(":",1)[0].strip()
                cc_buf = [head.split(":",1)[1].strip()] if ":" in head else [""]
                if cc_key not in CC_KEEP: cc_key, cc_buf = None, []
            elif cc_key is not None and not body.startswith("---") and not body.startswith("Copy"):
                cc_buf.append(body.strip())
        elif tag == "DR":
            db = body.split(";",1)[0].strip()
            if db in DR_KEEP:
                rec["dr"].setdefault(db, []).append(body.rstrip(".").strip())
        elif tag == "FT":
            if ln[5] != " ":                       # new feature line
                if ft_cur: rec["features"].append(ft_cur)
                parts = body.split(None, 1)
                ftype = parts[0]
                loc = parts[1].strip() if len(parts) > 1 else ""
                if ftype not in FEAT_KEEP:
                    ft_cur = None; continue
                ft_cur = {"type": ftype, "loc": loc, "note": "", "evidence": "", "ligand": ""}
            elif ft_cur is not None:
                # Continuation line of the feature opened above. It either starts a qualifier we
                # want, starts one we ignore, or continues the value of whichever qualifier was
                # opened last, because UniProt wraps long /note and /evidence values across lines.
                t = body.strip()
                if t.startswith("/note="):      ft_cur["note"]     = t[6:].strip('"')
                elif t.startswith("/evidence="):ft_cur["evidence"] = t[10:].strip('"')
                elif t.startswith("/ligand="):  ft_cur["ligand"]   = t[8:].strip('"')
                elif t.startswith("/"):         pass
                else:
                    # No leading slash, so this is a wrapped value. Append it to the first qualifier
                    # that already has content; only one can be open at a time, so the first hit is
                    # the right one. Without this, a wrapped evidence string loses its later PubMed
                    # ids and a wrapped note loses the text that identifies the modification.
                    for k in ("note","evidence","ligand"):
                        if ft_cur[k]: ft_cur[k] = ft_cur[k].rstrip('"') + " " + t.strip('"'); break
        elif tag == "SQ":
            in_seq = True
    if ft_cur: rec["features"].append(ft_cur)
    if cc_key: rec["cc"].setdefault(cc_key, []).append(" ".join(cc_buf).strip())
    if ref_cur: rec["refs"].append(ref_cur)
    rec["seq"] = "".join(seq_lines)
    for f in rec["features"]:
        f["pmids"] = sorted(set(PMID_RE.findall(f["evidence"])))
        f["eco"]   = sorted(set(ECO_RE.findall(f["evidence"])))
    rec["gene"] = rec["genes"][0] if rec["genes"] else None
    return rec

t0 = time.time(); n_tot = n_hum = 0
with gzip.open(SPROT, "rt", encoding="utf-8", errors="replace") as fh, \
     gzip.open(OUT, "wt", encoding="utf-8") as out:
    buf = []
    for ln in fh:
        if ln.startswith("//"):
            n_tot += 1
            if any(l.startswith("OX   NCBI_TaxID=9606") for l in buf):
                r = parse_entry(buf); n_hum += 1
                out.write(json.dumps(r, separators=(",", ":")) + "\n")
            buf = []
            if n_tot % 100000 == 0:
                log(f"  ...{n_tot:,} entries scanned, {n_hum:,} human, {time.time()-t0:.0f}s")
        else:
            buf.append(ln)
log(f"=== STEP 3: SWISS-PROT PARSE DONE: {n_tot:,} entries, {n_hum:,} human, {time.time()-t0:.0f}s -> {OUT.name} ({OUT.stat().st_size/1e6:.0f} MB)")
