# Analysis scripts

Literature triage of the unfiltered MoDPA PTM association list. Run from the project root, not from
this directory:

```
cd C:\Users\Enrico\Desktop\let-claude-cook
python scripts/s03_rawstats.py
```

Every script imports `scripts/_common.py`, which supplies the input paths, the `work/` output
directory, the `UNIMOD` table, and `log()`. `log()` prints and appends to `work/runlog.md`, so that
file is the complete record of the run in order.

Each script carries a header docstring saying what it reads, what it writes, and which section of
`20260911-modpa-literature-triage-report.md` it feeds.

## Order

Scripts are numbered in dependency order. Nothing runs in parallel, and re-running any script
overwrites its own outputs only.

| # | Script | Purpose | Approx. time |
|---|---|---|---|
| 01 | `s01_inventory.py` | inventory every input file | 30 s |
| 01b | `s01b_discrep.py` | follow up the discrepancies s01 found | 30 s |
| 02a | `s02a_to_parquet.py` | raw CSV to Parquet, once | 20 s |
| 02b | `s02b_verify.py` | row-count check on the Parquet copy | 10 s |
| 03 | `s03_rawstats.py` | characterise the raw list, score histogram | 10 s |
| 04 | `s04_nodes.py` | parse PTM event identifiers, build the node table | 10 s |
| 05 | `s05_parse_sprot.py` | stream Swiss-Prot, keep human entries | 100 s |
| 06 | `s06_strata_nruns.py` | |SDC| strata, join detection frequency | 15 s |
| 07 | `s07_pairs_base.py` | significant pairs with strata and identifier flags | 15 s |
| 08 | `s08_posqc.py` | positional QC against Swiss-Prot 2026_03 | 10 s |
| 09 | `s09_evbase.py` | build the offline evidence base | 10 s |
| 10 | `s10_tryptic.py` | shared-measurement flag (trypsin, 2 missed cleavages) | 10 s |
| 11 | `s11_nodeann.py` | per-event annotation, type (a) evidence | 5 s |
| 12 | `s12_pairev.py` | pair-level evidence tables, type (b) evidence | 5 s |
| 13 | `s13_tiers.py` | assign a tier to all 29.2 M pairs | 25 s |
| 14 | `s14_strat.py` | tier by stratum and sign, both built-in controls | 10 s |
| 15 | `s15_background.py` | matched permutation null, monotonicity test | 80 s |
| 16 | `s16_network_map.py` | map onto the 0.6 network and the clusters | 15 s |
| 17 | `s17_shortlist.py` | build the 165-pair shortlist | 80 s |
| 18 | `s18_controls2.py` | formal tests for the two controls | 10 s |
| 19 | `s19_deliverables.py` | write deliverables 1 and 2 | 220 s |
| 20 | `s20_misc.py` | assorted facts quoted in the report | 10 s |
| 21 | `s21_pmid_meta.py` | classify shortlist PubMed ids | 10 s |
| 22 | `s22_cluster_compare.py` | compare the two Leiden partitions | does not run, see below |
| 23 | `s23_cluster_from_clustersfile.py` | Step 8 cluster analysis, re-derived | 10 s |

To reproduce from scratch: 02a, then 03 through 21 in order, then 23. 01 and 01b are the initial
inventory and produce no files. 22 does not run.

## Scripts that no longer run

`20260911-True-data-edges-w-clusters.csv`, which held the unweighted Leiden partition, has been
withdrawn from the folder. `s01`, `s01b` and `s22` read it and now fail with `FileNotFoundError`.
They are kept because they document what was checked at the time, not because they are runnable.

No result in the report used the unweighted partition. `s23` re-derives the cluster analysis from
`20260911-True-data-Leiden-clusters.csv` alone and confirms it reproduces the stored `same_cluster`
column on 100.0000% of the 14,591 network edges.

## Principles the code follows

**No filtering beyond q < 0.05.** The co-quantification flag from `s10` is computed, carried and
reported, but never used to remove a pair. `s15` carries a sensitivity column and
`s18`/`s16` write the dropped-pair view to separate files, so the size of the effect is visible
without a filter being imposed anywhere.

**Tiers are assigned before any score is read.** The tier logic in `s13` touches no score column,
so the classification cannot be influenced by the quantity being tested.

**Type (a) and type (b) evidence never merge.** `s11` builds evidence that a site is a known
modified residue. `s12` builds evidence that a pair is linked. Tier 5 is type (a) only and is
reported in its own column throughout.

**Nothing large is loaded into memory whole.** The 1.3 GB association list is converted to Parquet
once and then queried with duckdb; Swiss-Prot is streamed entry by entry.

**Every citation is traceable.** No web or PubMed tool was used. Every PubMed id in the outputs
comes verbatim from a UniProt evidence string or `RX` line in `uniprot_2026_03_sprot.dat.gz`.

## Environment

Python 3.11.4, polars 1.44.2, pandas 3.0.3, duckdb 1.5.5, numpy 2.4.6, scipy 1.17.1,
matplotlib 3.10.8. Random seed 42 in `s15`.
