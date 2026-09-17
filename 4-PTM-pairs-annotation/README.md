# Step 4. PTM pair annotation and literature triage

Annotates the MoDPA association list against UniProt and the literature, and characterises the
edges that the score threshold keeps against the ones it discards.

## `MoDPA-UniProt-literature-search.ipynb`

Reproduces every number, table and stated result of the literature triage report
(`20260911-modpa-literature-triage-report.md`). It is written to be run top to bottom without
editing.

The notebook does not reimplement the analysis. It runs the scripts in `scripts/` in dependency
order, then loads the resulting tables from `work/` and renders each section of the report from
them. A step is skipped only when it previously exited 0 and all of its declared outputs are still
on disk, so a script that failed halfway is rerun rather than assumed complete. The full output of
every step is kept in `work/nb_logs/`, and `work/runlog.md` holds the same output interleaved in
run order.

Inputs are the study's own files plus Swiss-Prot release 2026_03. The release matters, because
tiers are assigned by site-level matching against it, so the notebook verifies the file by SHA-256
rather than by its name and downloads it from the UniProt archive if it is absent. The checksums of
all five inputs are recorded in the notebook. Set `VERIFY_INPUT_CHECKSUMS = True` to check the
local copies, which is worth doing once on a fresh copy of the data.

`scripts/README.md` documents the individual scripts, their order and their runtimes, and states
the principles the code follows: no filtering beyond q < 0.05, tiers assigned before any score is
read, type (a) and type (b) evidence never merged, nothing large loaded into memory whole, and
every PubMed identifier taken verbatim from a UniProt evidence string rather than from a web
search.

`scripts/s10_tryptic.py` holds the co-quantification test: it digests each protein sequence with
trypsin at up to two missed cleavages and asks whether one peptide can cover both sites, in a
standard variant and in an upper-bound variant that treats a modified lysine or arginine as
non-cleavable. This supersedes the `shared_peptide` and `potential_artefact` columns of
`parse_network`, which use a fixed `position_gap <= 5` proximity rule. The script also reports how
far the two rules disagree. The flag is carried and reported throughout stage 4 but is never used
to remove a pair, because 61 of the 85 pairs UniProt documents as explicit crosstalk fall on one
peptide.

## `Compare-retained-discarded-edges-nodes.ipynb`

Compares the edges and nodes that survive the `Score >= 0.6` threshold with those that do not,
across the model runs discovered under `MODEL_ROOT`. It imports `parse_network`, `build_graph`,
`discover_runs` and `jaccard` from `modpa_network_utilities.py` in the repository root, which it
adds to `sys.path`, so it must be started from within this directory.

The notebook is a copy of `../2-VAE-code/Sensitivity_analysis/MoDPA-sensisitivity-analysis.ipynb`
with the retained and discarded comparison as its focus. The sensitivity analysis README documents
the input layout, the configuration parameters, the filtering and the memory requirements, and all
of it applies here as well.

## Data

`input_and_expected_ouput.tar.gz` holds the inputs and the expected outputs of this step. It is not
tracked in git and is part of the Zenodo release. See the data availability section of the root
README.
