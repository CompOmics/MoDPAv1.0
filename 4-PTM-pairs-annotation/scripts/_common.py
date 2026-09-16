"""
Shared configuration for the MoDPA literature-triage analysis.

Provides:
  * BASE / WORK      - project root and the work/ directory for intermediate tables.
  * DATE             - the YYYYMMDD prefix used on every delivered file.
  * log()            - prints to stdout and appends to work/runlog.md, so the run log ends up
                       a complete record of what was computed and in what order.
  * path constants   - RAW, CLUSTERS, EDGES, FILTERED, NRUNS, SPROT.
  * UNIMOD           - the eleven modifications in the MoDPA matrix, mapped to a display name
                       and the set of residues that modification can chemically occupy.

The residue sets in UNIMOD are used by the positional QC in s08 to test that each PTM event sits
on a chemically plausible residue. The names were cross-checked against the ptm_name column of
N-runs-per-ptm.csv.gz (format "[21]Phospho"), which is the authoritative mapping shipped with the
data, rather than against the manuscript, which is outdated. Note that Unimod 535 is the
Leu-Arg-Gly-Gly ubiquitin remnant, named LRGG in that column.

EDGES points at 20260911-True-data-edges-w-clusters.csv, which holds the unweighted Leiden
partition. That file has been withdrawn from the folder, so the scripts that read it (s01, s01b,
s22) no longer run. Every cluster assignment used in the report comes from CLUSTERS, the weighted
partition.
"""
import os, sys, datetime, pathlib
os.environ.setdefault("PYTHONIOENCODING","utf-8")
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import polars as pl
pl.Config.set_tbl_formatting("ASCII_FULL")
pl.Config.set_tbl_rows(60)
pl.Config.set_tbl_cols(40)
pl.Config.set_fmt_str_lengths(120)

BASE = pathlib.Path(__file__).resolve().parent.parent
WORK = BASE/"work"; WORK.mkdir(exist_ok=True)
DATE = "20260911"
RUNLOG = WORK/"runlog.md"

def log(msg):
    line = f"[{datetime.datetime.now():%H:%M:%S}] {msg}"
    print(line)
    with open(RUNLOG,"a",encoding="utf-8") as fh:
        fh.write(line+"\n")

RAW      = BASE/"20260911-1042-keen_raman-signed-distances.csv.gz"
CLUSTERS = BASE/"20260911-True-data-Leiden-clusters.csv"
EDGES    = BASE/"20260911-True-data-edges-w-clusters.csv"
FILTERED = BASE/"20260911-True-data-filtered-distances.csv"
NRUNS    = BASE/"N-runs-per-ptm.csv.gz"
SPROT    = BASE/"uniprot_2026_03_sprot.dat.gz"

# Unimod ID -> (name, allowed residues)
UNIMOD = {
    1:   ("Acetyl",        set("K")),
    7:   ("Deamidated",    set("R")),     # citrullination on R in this study
    21:  ("Phospho",       set("STY")),
    23:  ("Dehydrated",    set("STY")),
    34:  ("Methyl",        set("KR")),
    36:  ("Dimethyl",      set("KR")),
    37:  ("Trimethyl",     set("K")),
    53:  ("HNE",           set("CHK")),
    64:  ("Succinyl",      set("K")),
    299: ("Carboxy",       set("KDE")),
    535: ("GG_ubiquitin",  set("K")),
}
