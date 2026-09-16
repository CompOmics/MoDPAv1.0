# Step 5. Reactome pathway over-representation analysis

Tests the proteins carrying the PTM events of each network cluster for over-representation of
Reactome pathways.

## `20260901-reactome-enrichment.ipynb`

The pipeline is:

1. download the human Reactome pathway hierarchy once and cache it as
   `ReactomePathwaysHierarchyHuman-v<version>.csv`;
2. load the clustered PTM events and summarise their composition per cluster;
3. run an over-representation analysis (ORA) per cluster through the Reactome Analysis Service,
   caching the full result of every cluster in `per-cluster/`;
4. remove redundant pathways using the hierarchy, in both a keep-specific and a keep-general
   variant;
5. plot the results as `-log10(FDR)` heatmaps, optionally with the rows ordered by hierarchy
   distance.

Only the `CONFIG` cell needs editing to rerun the analysis on a new clustering.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `NODES_PATH` | a Leiden cluster CSV in the reference model folder | clustered PTM events, one row per event |
| `CLUSTER_COL` | `LeidenCluster` | cluster membership column |
| `REACTOME_VERSION` | `97` | Reactome release used for the hierarchy |
| `MIN_PROTEINS` | `20` | clusters with fewer distinct proteins are not tested |
| `FDR_CUTOFF` | `0.01` | significance threshold |
| `MIN_ENTITIES` | `5` | minimum number of pathway entities found, passed to Reactome |
| `OUT_DIR` | `<today>-keen_raman-Reactome` | output directory |

The ORA queries Reactome with `p_value="1"`, so every tested pathway is returned and cached, not
only the significant ones. Filtering to `FDR < FDR_CUTOFF` happens afterwards, before the
redundancy step. This matters for the sensitivity analysis below, which needs the full per-cluster
candidate set in order to correct correctly.

Outputs, written to `OUT_DIR`:

| File | Content |
| --- | --- |
| `nicely-formatted-nodes.csv` | the clustered PTM events with readable modification names |
| `cluster-sizes.csv` | events and distinct proteins per cluster |
| `cluster-composition-counts.csv`, `cluster-composition-percent.csv` | modification composition per cluster |
| `cluster-summary.csv` | per cluster, the number of pathways tested and kept by each redundancy rule |
| `enrichment-all-clusters.csv` | every tested pathway of every cluster, with both redundancy flags. Delivered as `SuppData5-enrichment-all-clusters.csv` |
| `per-cluster/cluster<N>.csv` | the PTM events submitted for cluster N |
| `per-cluster/cluster<N>-reac.csv` | every pathway Reactome tested for cluster N. Delete to force a rerun of that cluster |
| `ptms-per-cluster.png`, `.svg` | modification composition figure |
| `reactome-clustermap-*.png`, `.svg` | the `-log10(FDR)` heatmaps |

## `reactome-network-background-ora.py`

A supplementary sensitivity analysis. Reactome's Analysis Service always tests against the full
annotated human reference proteome; the v3 REST API offers no custom background parameter at all.
The manuscript reports the whole-proteome background, which is Reactome's standard behaviour. This
script quantifies how much that choice matters.

It submits the union of the network's proteins as a single ORA query to obtain, for every pathway,
how many network proteins are annotated to it. It then recomputes a one-tailed hypergeometric test
per pathway and cluster against that network-restricted background, and applies Benjamini-Hochberg
correction within each cluster over every pathway tested for that cluster, not only the ones that
were significant under the whole-proteome background.

```bash
python reactome-network-background-ora.py
```

`OUT_DIR` at the top of the file is the only thing to change when the pipeline is rerun on a new
clustering. The single Reactome query is cached in `OUT_DIR/network-background-pathway-sizes.json`,
and the script runs fully offline when that file is present, with no network access and no
`reactome2py` dependency. Delete the cache to force a refresh. On a cache miss the script prefers
`reactome2py` and falls back to a direct `urllib` call against the same endpoints.

Writes `OUT_DIR/network-background-reanalysis.csv`, which places the whole-proteome FDR, the
network-restricted raw p-value and the corrected network-restricted FDR side by side, for every
pathway tested.

## `2026-09-13-keen_raman-Reactome/`

The results of the run of 2026-09-13 on the reference network, including the cached Reactome query
that makes the sensitivity analysis reproducible offline.
