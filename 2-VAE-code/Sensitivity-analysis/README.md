# Latent-space sensitivity analysis of the PTM co-regulation network

`MoDPA-sensisitivity-analysis.ipynb` quantifies how the inferred PTM co-regulation
network changes with the latent-space dimensionality of the MoDPA model. Every model run
is compared against a single reference run at the dimensionality used in the manuscript
(128). Four replicate runs at that same dimensionality are carried through every
comparison, so the deviation caused by the latent dimensionality can be read against
ordinary run-to-run variation at a fixed dimensionality.

## Contents

| File | Description |
| --- | --- |
| `MoDPA-sensisitivity-analysis.ipynb` | The analysis notebook |
| `topology_compare.py` | Helper module with the topology comparison functions |
| `2026-09-15_sensitivity_output/` | Figures and tables produced by the run of 2026-09-15 |
| `Supplementary-Methods-Latent-Dimensionality-Sensitivity.docx` | Supplementary methods text |

## Dependencies outside this directory

The notebook imports `modpa_network_utilities.py` from the repository root. The root is
added to `sys.path` in the first cell, so the notebook must be started from within this
directory. Four functions are used: `discover_runs`, `parse_network`, `build_graph` and
`jaccard`.

`topology_compare.py` is imported inside the `topology_summary` function rather than at
the top of the notebook, and provides `compare_sets`, `threshold_proximity`,
`indirectness_report`, `common_neighbor_enrichment`, `partition_stability` and
`plot_cn_comparison`.

No FASTA file and no protein-level annotation are needed. The notebook works on node
identifiers only.

## Input data

The notebook expects one directory per model run under `MODEL_ROOT`, which defaults to
`../2-VAE-code/2026-09-10-HumanOnly`:

```
2026-09-10-HumanOnly/
├── 20260910-1839-s729727809-reverent_bouman/
│   ├── config.json                                  # contains "latent_dim"
│   └── ...-signed-distances.csv.gz
├── ...
├── 128d/                                            # replicate runs at the reference dimensionality
│   ├── <run>/
│   │   ├── config.json
│   │   └── ...-signed-distances.csv.gz
│   └── ...
└── 128d-REF/                                        # exactly one run: the reference
    └── <run>/
        ├── config.json
        └── ...-signed-distances.csv.gz
```

`discover_runs` scans three locations: `MODEL_ROOT` itself, where `EXCLUDED_DIMS`
applies; `MODEL_ROOT / REPLICATE_SUBDIR`; and `MODEL_ROOT / REFERENCE_SUBDIR`, which must
hold exactly one run. The reference and replicate directories are skipped during the scan
of `MODEL_ROOT`, so no run is counted twice.

A run directory is one that holds a `config.json` with a `latent_dim` field. If the
reference or replicate directory holds a `config.json` itself, it is treated as a single
run rather than as a parent of runs. A directory without a `config.json`, and a run
directory without a `*-signed-distances.csv.gz` file, are both skipped without an error.
A run directory holding more than one `*-signed-distances.csv.gz` file raises a
`ValueError`.

The signed-distance files have the columns `nodeA`, `nodeB`, `Score`, `pval`, `distance`,
`PCC` and `qvalue`. `Score` is the signed distance correlation (SDCor) between two PTMs,
and `distance` is its magnitude: `distance == abs(Score)` throughout, with the sign of
`Score` taken from `PCC`. Node identifiers use the format
`UniProtAccession|Position|Residue|UnimodID`, for example `P04637|15|S|21`.

Each compressed edge list is 1.2 GB to 1.3 GB, and the default `MODEL_ROOT` holds 14 runs
at the latent dimensionalities 32, 64, 96, 124, 128 (five runs: the reference and four
replicates), 132, 160, 192, 256 and 512. `REPLICATE_RUNS` therefore holds four labels and
`COMPARISON_RUNS` nine.

The model outputs are in `../MoDPA_models.tar.gz`, which is part of the Zenodo release of this
repository and is not tracked in git. See the data availability section of the root README.

## Installation

The notebook runs in the environment of the repository, created from the repository root:

```bash
conda env create -f env.yml
conda activate MoDPA-env
```

That environment covers everything this directory needs: `polars`, `pandas`, `matplotlib`,
`seaborn`, `scipy` and `jupyterlab` for the notebook, and `numpy`, `networkx`, `igraph`,
`leidenalg` and `scikit-learn` for `modpa_network_utilities.py` and `topology_compare.py`.

`leidenalg` and `igraph` ship wheels for the common platforms. If a source build is
triggered, a C compiler and the igraph C library headers are required; installing
through conda-forge avoids this.

## Running

1. Download the model outputs and set `MODEL_ROOT` in the configuration cell.
2. Confirm that `REFERENCE_SUBDIR` and `REPLICATE_SUBDIR` point at the intended
   directories.
3. Restart the kernel and run all cells.

`OUTPUT_DIR` defaults to `./{today}_sensitivity_output` and is created by the
configuration cell, so a run never writes into the output of an earlier run. All figures
and tables are written there.

The edge lists are parsed once with `min_score=0`, in `lazy_zero_filtered_edgelists`, and
`MIN_SCORE` is applied afterwards when the graphs are built, so that the same parsed
tables serve both the unfiltered and the filtered views.

## Run labels and the reference

The reference run is the single run under `MODEL_ROOT / REFERENCE_SUBDIR`, and
`discover_runs` raises an error if that directory does not hold exactly one run. Every
other run is labelled `{latent_dim:03d}_rep{n}`, with replicate numbers assigned in order
of directory name, so labels are stable across executions and across machines. The
reference is labelled `{latent_dim:03d}_ref` and sorts first within its dimensionality.
`discover_runs` returns the label to directory mapping as the `RUNS` manifest, indexed by
label, with the columns `latent_dim`, `replicate`, `is_reference`, `directory` and `path`.

Three run groups are derived from the manifest and used throughout:

| Variable | Contents |
| --- | --- |
| `REFERENCE_RUN` | the single reference label, `128_ref` |
| `REPLICATE_RUNS` | the four replicate runs at 128 dimensions |
| `COMPARISON_RUNS` | the nine runs at a dimensionality other than 128 |

`REPLICATE_RUNS` and `COMPARISON_RUNS` are selected by the literal label prefix `128_`.
Change those two lines as well if the reference dimensionality changes.

Most sections iterate over `REPLICATE_RUNS + COMPARISON_RUNS`, that is, 13 comparisons
against the reference. The replicates therefore appear in every figure and table as the
within-dimensionality baseline.

## Filtering

`parse_network` keeps an edge when `Score >= min_score` **and** `qvalue < 0.05`. The score
test uses the **signed** score, so the analysis covers positively correlated PTM pairs
only and discards strongly anti-correlated ones. The q-value test always applies, also in
the `min_score=0` tables, which are therefore q-value filtered but not score filtered.

Two filtering levels are used:

* `min_score=0`, in `lazy_zero_filtered_edgelists`, used for the score distributions and
  for the node-set comparison. At this level all 14 runs hold the same 7,644 nodes, and
  the node Jaccard index between any two runs is 1.000. The edge counts range from
  859,312 (32 dimensions) to 21,422,640 (512 dimensions), with 17,102,635 in the
  reference.
* `Score >= MIN_SCORE`, 0.6 by default, applied to the same tables when
  `replicate_graphs_filtered` is built, and used for every network comparison. The
  reference network then holds 6,491 nodes and 14,591 edges.

The score threshold is applied in two places, the cell that builds `REF_filtered` and the
cell that builds `replicate_graphs_filtered`. Both read `MIN_SCORE`, so the threshold is
changed in one place, in the configuration cell.

`parse_network` also annotates each edge with `same_protein`, `position_gap`, `same_site`,
`same_mod`, `shared_peptide` and `potential_artefact`. An edge is a potential artefact
when both PTMs carry the same modification and sit within `adjacency_window` (5) residues
of each other on the same protein, which is the configuration in which a shared peptide
can produce a correlation that is not biological. These edges are annotated but not
removed.

> **`shared_peptide` and `potential_artefact` are rough approximations and should be
> treated with caution.** The proximity rule behind them is not a statement about what was
> actually measured together. The tryptic test in
> `../../4-PTM-pairs-annotation/scripts/s10_tryptic.py` is the accurate one. See the root
> README. Both columns are scheduled for removal.

No section of this notebook reads either column, so the caveat does not affect any result
reported here.

## Analyses

| Section | Measure | Reference point |
| --- | --- | --- |
| Score distribution | `Score` of the unfiltered edges of each run, split by whether the edge is in the filtered reference network | replicate runs at 128 |
| Retained edge fraction | fraction of the unfiltered edges of a run that the reference network retains | replicate runs at 128 |
| Node sets | Jaccard index between node sets, unfiltered and at `Score >= MIN_SCORE` | all pairs of runs |
| Node degree distribution | degree of the filtered nodes of each run, split by presence in the reference network | replicate runs at 128 |
| Topology of differing edges | shortest-path distance and common neighbours of the reference edges missing from a run | random non-adjacent node pairs |
| Partition stability | ARI, NMI and AMI of the Leiden partitions on the common subgraph | replicate runs at 128 |
| Threshold proximity | median weight of the shared and of the run-specific edges in each network | the two networks against each other |
| Edge weight correlation | Pearson r between the `Score` values of the edges shared by a run and the reference | replicate runs at 128 |
| Node degree correlation | Spearman rho between the node degrees shared by a run and the reference | replicate runs at 128 |

The topology measures characterise the edges present in the reference network but absent
from a comparison network, always evaluated inside the network that lacks them, and
always against a null of random non-adjacent node pairs:

* `indirectness_report` gives the shortest-path distance between the endpoints of a
  differing edge. It covers triangle closure (distance 2) and square closure (distance 3)
  in one measure, so it reaches up to path length 3.
* `common_neighbor_enrichment` counts the shared neighbours of the endpoints, which is
  the number of triangles the edge would close. `prob_superiority` is the probability
  that a differing edge has more common neighbours than a random non-edge, with 0.5
  meaning no effect.
* `partition_stability` restricts both networks to their shared nodes, runs Leiden on
  each, and compares the two partitions. The columns are suffixed `_common_subgraph`
  because the measurement is defined on the common subgraph and not on the full networks.

Only the `G2_minus_G1` direction is computed, that is, the edges present in the reference
network and absent from the comparison network. The reverse direction is not reported,
and skipping it halves the runtime of that section.

## Outputs

Written to `OUTPUT_DIR`:

| File | Content |
| --- | --- |
| `retained-discarded-edges-scores.svg` | Score distribution per run, split by presence in the reference network |
| `retained_edge_fraction.csv` | Edge counts per run and the fraction retained by the reference |
| `retained-discarded-nodes-degree.svg` | Node degree distribution per run, split by presence in the reference network |
| `retained_node_fraction.csv` | Node counts per run and the fraction retained by the reference |
| `topology_comparison.csv` | All topology metrics, one row per run |
| `common_neighbour_enrichment_{run}_vs_128_ref.png` | Common neighbours of the differing edges against random non-edges, one figure per run |
| `differing_vs_random_triangle_fraction.svg` | Fraction of differing edges at distance 2 and distance 3, against the random baseline |
| `NMI_common_subgraph_across_runs.svg` | NMI of the Leiden partitions on the common subgraph |
| `edge_weight_correlation_{run}_vs_128_ref.png` | Scatter of the shared edge weights, one figure per run |
| `edge_weight_correlation_across_runs.svg` | Pearson r of the shared edge weights across runs |
| `node_degree_correlation_{run}_vs_128_ref.png` | Scatter of the shared node degrees, one figure per run |
| `node_degree_correlation_across_runs.svg` | Spearman rho of the shared node degrees across runs |
| `edge_weight_node_degree_correlation.csv` | Both correlations, their p-values and the shared counts, one row per run |

The two summary figures use fixed y-axis limits, `plt.ylim(0.7, 0.9)` for the Pearson r
and `plt.ylim(0.0, 0.7)` for the Spearman rho. Widen them if a run falls outside.

## Configuration

All parameters are in the second code cell:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `MIN_SCORE` | `0.6` | Score threshold of the filtered networks, read in both places that filter |
| `EXCLUDED_DIMS` | `{}` | Latent dimensionalities skipped when scanning `MODEL_ROOT` |
| `MODEL_ROOT` | `"../2-VAE-code/2026-09-10-HumanOnly"` | Directory holding one subdirectory per model run |
| `REFERENCE_SUBDIR` | `"128d-REF"` | Directory holding the single reference run |
| `REPLICATE_SUBDIR` | `"128d"` | Directory holding the replicate runs at the reference dimensionality |
| `LEIDEN_RESOLUTION` | `2` | Resolution parameter of the RBConfiguration objective, used by `partition_stability` |
| `LEIDEN_N_ITERATIONS` | `2` | Iterations of the Leiden optimiser |
| `LEIDEN_SEED` | `42` | Seed of the Leiden runs in the topology section |
| `N_TOPOLOGY_BASELINE` | `100000` | Random non-edges drawn as the null in the topology section |
| `TOPOLOGY_SEED` | `0` | Seed for drawing those non-edges |
| `FIG_DPI` | `300` | Resolution of the exported PNG figures |
| `OUTPUT_DIR` | `./{today}_sensitivity_output` | Output directory, created by the same cell |

No dimensionality is skipped by default. Note that the default `{}` is an empty dict
rather than an empty set; both are falsy, so both behave the same here, but use a set
literal or `set()` when excluding dimensionalities.

`LEIDEN_SEEDS`, `N_PERMUTATIONS` and `PERMUTATION_SEED` are still defined in that cell but
no current section reads them. They are left in place for the seed sweep and the
permutation baseline, which are not part of this version of the analysis.

## Memory and runtime

Memory is the binding constraint, not CPU time. The notebook holds all 14 unfiltered
graphs in `replicate_graphs` at the same time, which is about 165 million edges in total,
and then builds the 14 filtered graphs alongside them. The score distribution cell
concatenates the unfiltered edge tables of 13 runs into a single eager DataFrame, which is
of the same order. Run the notebook on a machine with a large amount of RAM, or delete
`replicate_graphs` and `tmp` once the sections that use them have run.

Reading the edge lists dominates the wall-clock time. Each of the 14 compressed files is
1.2 GB to 1.3 GB, and `parse_network` is called with `min_score=0`, so the full tables are
streamed. Run the notebook against fast local storage.

The topology section is comparatively cheap at the default `MIN_SCORE`. On the reference
network at `MIN_SCORE = 0.6`, which holds 6,491 nodes and 14,591 edges, drawing 100,000
random non-edges takes about 0.5 s, the endpoint distances about 0.5 s, the common
neighbours about 0.4 s, and `partition_stability` about 0.3 s, so roughly 2 s per
comparison and about 30 s for the 13 comparisons. The cost grows quickly as the threshold
is lowered and the networks become denser: at 512 dimensions the filtered network still
holds 398,405 edges. Lower `N_TOPOLOGY_BASELINE` for a quick pass; the baseline fractions
are stable well below 100,000 draws.

## Self test of `topology_compare.py`

`topology_compare.py` can be checked on its own without the model outputs:

```bash
python topology_compare.py
```

This runs an end-to-end test on two synthetic planted-partition graphs. Every step asserts
the property it is meant to hold, covering the edge-set partition returned by
`compare_sets`, the seeding and the attempt limit of `sample_non_edges`, the cutoff and
the skip counter of `endpoint_distances`, the distance classes of `summarise_distances`,
the tie correction in `cn_effect_size`, the direction handling of `indirectness_report`,
the `KeyError` raised by `threshold_proximity` on an unweighted graph, and an ARI of
exactly 1.0 when `partition_stability` is given the same graph twice. It prints a summary
and ends with `self test passed`; any failure raises an `AssertionError`. It takes a few
seconds and needs no model output.

Run it with a plain `python`. Under `python -O` the assertions are stripped and the test
checks nothing.

## Reproducibility

`../../env.yml` is the archived environment. A few of its entries are unpinned, `igraph`,
`leidenalg` and `scikit-learn` among them, so pin them before archiving:

```bash
conda env export --no-builds > env.yml
```

## License

Licensed under the Apache License, Version 2.0.

Unless required by applicable law or agreed to in writing, software distributed under the
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
either express or implied. See the License for the specific language governing permissions
and limitations under the License.
