#!/usr/bin/env python
# coding: utf-8

# Utility functions to build and analyze the MoDPA network from PTM-PTM associations

import json
import os

import igraph as ig
import leidenalg
import networkx as nx
import numpy as np
import pandas as pd
import polars as pl


def build_graph(df: pl.DataFrame) -> nx.Graph:
    """Build a NetworkX graph from a filtered MoDPA edge table.

    Duplicate or reversed representations of the same undirected edge are
    collapsed explicitly. If their scores differ, the score with the largest
    absolute value is retained because Leiden clustering uses absolute SDC as
    its edge weight. Iterating over Polars rows avoids an unnecessary pandas /
    PyArrow conversion.
    """
    G = nx.Graph()
    for node_a, node_b, score in df.select(["nodeA", "nodeB", "Score"]).iter_rows():
        score = float(score)
        if G.has_edge(node_a, node_b):
            if abs(score) > abs(G[node_a][node_b]["Score"]):
                G[node_a][node_b]["Score"] = score
        else:
            G.add_edge(node_a, node_b, Score=score)
    return G


def run_leiden(
    G: nx.Graph,
    seed: int = 42,
    resolution: float = 2.0,
    n_iterations: int = 2,
) -> dict[str, int]:
    """Run Leiden clustering and return ``{node_name: cluster_id}``."""
    if G.number_of_nodes() == 0:
        return {}

    # Use absolute SDC as weights (Leiden requires non-negative weights)
    weights = [abs(G[u][v]["Score"]) for u, v in G.edges()]
    nx.set_edge_attributes(G, dict(zip(G.edges(), weights)), "abs_sdc")
    G_ig = ig.Graph.from_networkx(G)

    partition = leidenalg.RBConfigurationVertexPartition(
        G_ig,
        weights="abs_sdc",
        resolution_parameter=resolution,
    )
    optimiser = leidenalg.Optimiser()
    optimiser.set_rng_seed(seed)
    optimiser.optimise_partition(partition, n_iterations=n_iterations)

    return {
        G_ig.vs[i]["_nx_name"]: partition.membership[i]
        for i in range(G_ig.vcount())
    }


def build_clusters_df(leiden_membership: dict[str, int]) -> pd.DataFrame:
    """Convert Leiden membership into one annotated row per phosphosite node."""
    cluster_rows = []
    for node, cluster in leiden_membership.items():
        uniacc, pos, res, unimod_id = node.split("|")
        cluster_rows.append([node, cluster, uniacc, pos, res, unimod_id])
    return pd.DataFrame(
        cluster_rows,
        columns=["name", "LeidenCluster", "UniAcc", "POS", "RES", "UniModID"],
    )


def add_protein_annotations(clusters_df: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    """Add protein annotations to a phosphosite cluster table by UniProt accession."""
    if "UniAcc" not in clusters_df.columns:
        raise KeyError("clusters_df must contain a 'UniAcc' column")
    if "UniAcc" not in proteins.columns:
        raise KeyError("proteins must contain a 'UniAcc' column")

    protein_annotations = proteins.loc[proteins["UniAcc"].notna()].copy()
    duplicated_accessions = protein_annotations.loc[
        protein_annotations["UniAcc"].duplicated(keep=False), "UniAcc"
    ].unique()
    if len(duplicated_accessions):
        examples = ", ".join(map(str, duplicated_accessions[:5]))
        raise ValueError(
            "proteins must contain at most one annotation row per UniAcc; "
            f"duplicate accession(s) include: {examples}"
        )

    missing_accessions = clusters_df.loc[
        ~clusters_df["UniAcc"].isin(protein_annotations["UniAcc"]), "UniAcc"
    ].unique()
    if len(missing_accessions):
        examples = ", ".join(map(str, missing_accessions[:5]))
        raise ValueError(
            "Missing protein annotations for network accession(s); "
            f"examples: {examples}"
        )

    return clusters_df.merge(
        protein_annotations,
        on="UniAcc",
        how="left",
        validate="many_to_one",
    ).sort_values("LeidenCluster").reset_index(drop=True)


def parse_network(
    p,
    min_score: float,
    qvalue: float = 0.05,
    adjacency_window: int = 5,
) -> pl.LazyFrame:
    """Lazily parse, filter, and annotate a MoDPA edge table.

    CSV scanning, score/q-value filtering, node decomposition, and edge
    annotation are kept in a lazy query. A small validation query is
    materialised first so malformed node identifiers fail early and clearly.

    Node identifiers must have the form
    ``UniProt|Position|Residue|UnimodID``, with a numeric position and four
    non-empty components. Malformed identifiers raise an informative
    :class:`ValueError` before the annotated lazy query is returned.
    """
    # Filter as early as possible so rejected rows never reach the more
    # expensive string parsing and annotation expressions.
    edges = pl.scan_csv(p).filter(
        pl.col("Score") >= min_score,
        pl.col("qvalue") < qvalue,
    )

    # Fixed-size structs avoid eager Python-level parsing and remain compatible
    # with LazyFrame execution. The validation query below projects only the
    # node columns and stops after a few malformed examples are found.
    edges = (
        edges.with_columns(
            pl.col("nodeA")
            .cast(pl.String)
            .str.split_exact("|", 3)
            .struct.rename_fields(["protA", "posA", "resA", "modA"])
            .alias("_nodeA"),
            pl.col("nodeB")
            .cast(pl.String)
            .str.split_exact("|", 3)
            .struct.rename_fields(["protB", "posB", "resB", "modB"])
            .alias("_nodeB"),
        )
        .unnest(["_nodeA", "_nodeB"])
        .with_columns(
            pl.col("posA").cast(pl.Int64, strict=False),
            pl.col("posB").cast(pl.Int64, strict=False),
        )
    )

    valid_node_a = (
        (pl.col("nodeA").cast(pl.String).str.count_matches(r"\|") == 3)
        & pl.col("protA").is_not_null()
        & (pl.col("protA") != "")
        & pl.col("posA").is_not_null()
        & pl.col("resA").is_not_null()
        & (pl.col("resA") != "")
        & pl.col("modA").is_not_null()
        & (pl.col("modA") != "")
    ).fill_null(False)
    valid_node_b = (
        (pl.col("nodeB").cast(pl.String).str.count_matches(r"\|") == 3)
        & pl.col("protB").is_not_null()
        & (pl.col("protB") != "")
        & pl.col("posB").is_not_null()
        & pl.col("resB").is_not_null()
        & (pl.col("resB") != "")
        & pl.col("modB").is_not_null()
        & (pl.col("modB") != "")
    ).fill_null(False)

    malformed = (
        edges.filter(~(valid_node_a & valid_node_b))
        .select(["nodeA", "nodeB"])
        .limit(5)
        .collect(engine="streaming")
    )
    if malformed.height:
        examples = "; ".join(
            f"nodeA={node_a!r}, nodeB={node_b!r}"
            for node_a, node_b in malformed.iter_rows()
        )
        raise ValueError(
            "Malformed node identifier(s) after score/q-value filtering. "
            "Expected 'UniProt|Position|Residue|UnimodID' with exactly four "
            f"non-empty fields and an integer Position. Examples: {examples}"
        )

    same_protein = pl.col("protA") == pl.col("protB")
    edges = edges.with_columns(
        same_protein.alias("same_protein"),
        pl.when(same_protein).then(
            (pl.col("posB") - pl.col("posA")).abs()
        ).otherwise(None).alias("position_gap"),
        (same_protein & (pl.col("posA") == pl.col("posB"))).alias("same_site"),
        (pl.col("modA") == pl.col("modB")).alias("same_mod")
    ).with_columns(
        (
            pl.col("same_protein")
            & (pl.col("position_gap") <= adjacency_window)
        ).fill_null(False).alias("shared_peptide")
    ).with_columns(
        (pl.col("same_mod") & pl.col("shared_peptide")).alias("potential_artefact"),
        pl.concat_str(
            pl.when(
                pl.col("nodeA") < pl.col("nodeB")
            ).then(pl.col("nodeA")).otherwise(pl.col("nodeB")),
            pl.lit("__"),
            pl.when(
                pl.col("nodeA") < pl.col("nodeB")
            ).then(pl.col("nodeB")).otherwise(pl.col("nodeA"))
        ).alias("pair_key"),
        # pl.col("Score").rank(method="average").alias("Rank"),
    ).filter(
        ~pl.col("potential_artefact")
    ).sort(
        "Score"
    )

    return edges


def count_edges(network):
    return len(set(network['pair_key']))


def count_nodes(network):
    sources = list(network["nodeA"])
    targets = list(network["nodeB"])
    return len(set(sources+targets))


def jaccard(set_a: set, set_b: set) -> tuple[float, float]:
    """Return Jaccard similarity and the fraction of ``set_b`` that is shared.

    An undefined ratio is returned as ``numpy.nan``: Jaccard similarity is
    undefined when both sets are empty, and the shared fraction is undefined
    whenever ``set_b`` is empty.
    """
    shared = len(set_a & set_b)
    union_size = len(set_a | set_b)
    jaccard_index = shared / union_size if union_size else np.nan
    shared_fraction_b = shared / len(set_b) if set_b else np.nan
    return jaccard_index, shared_fraction_b


def edge_set(edges: pl.DataFrame) -> set[str]:
    return set(edges["pair_key"])


def node_set(edges: pl.DataFrame) -> set[str]:
    return set(edges["nodeA"]) | set(edges["nodeB"])


# --- List latent space test runs --- #
def _run_directories(parent, skip: set[str] = frozenset()):
    """Yield the model-run directories under `parent`.

    A run directory is one that contains a `config.json`. If `parent` itself
    contains one, `parent` is treated as a single run.
    """
    if not os.path.isdir(parent):
        return
    if os.path.isfile(os.path.join(parent, "config.json")):
        yield parent
        return
    for child_name in sorted(os.listdir(parent)):
        child = os.path.join(parent, child_name)
        if (
            child_name not in skip
            and os.path.isdir(child)
            and os.path.isfile(os.path.join(child, "config.json"))
        ):
            yield child


def discover_runs(
    model_root,
    reference_subdir: str,
    replicate_subdir: str,
    excluded_dims: set[int] = frozenset(),
) -> pd.DataFrame:
    """Locate the signed-distance file of every model run under `model_root`.

    Three locations are scanned: `model_root` itself, where `excluded_dims`
    applies; `os.path.join(model_root, replicate_subdir)`, holding replicate runs
    at one dimensionality; and `os.path.join(model_root, reference_subdir)`, which
    must hold exactly one run, the reference. Placing the reference in its own
    directory means the comparisons do not depend on the order in which the
    filesystem lists
    directories.

    Runs are labelled `{latent_dim:03d}_ref` for the reference and
    `{latent_dim:03d}_rep{n}` otherwise, with replicate numbers assigned in order
    of directory name so that the labels are stable across machines.

    Returns a manifest indexed by label, with columns `latent_dim`,
    `replicate`, `is_reference`, `directory` and `path`.
    """
    special = {reference_subdir, replicate_subdir}
    found = []
    reference_root = os.path.join(model_root, reference_subdir)
    replicate_root = os.path.join(model_root, replicate_subdir)
    for parent, apply_exclusion, is_reference, skip in [
        (model_root, True, False, special),
        (replicate_root, False, False, frozenset()),
        (reference_root, False, True, frozenset()),
    ]:
        for run_dir in _run_directories(parent, skip):
            matches = sorted(
                os.path.join(run_dir, filename)
                for filename in os.listdir(run_dir)
                if filename.endswith("-signed-distances.csv.gz")
                and os.path.isfile(os.path.join(run_dir, filename))
            )
            if not matches:
                continue
            if len(matches) > 1:
                raise ValueError(
                    f"{run_dir} contains {len(matches)} signed-distance files; "
                    "exactly one is expected"
                )
            config_path = os.path.join(run_dir, "config.json")
            with open(config_path, "r") as file:
                latent_dim = int(json.load(file)["latent_dim"])
            if apply_exclusion and latent_dim in excluded_dims:
                continue
            found.append(
                {
                    "latent_dim": latent_dim,
                    "is_reference": is_reference,
                    "directory": os.path.basename(run_dir),
                    "path": matches[0],
                }
            )

    if not found:
        raise FileNotFoundError(f"no model runs found under {model_root}")

    manifest = pd.DataFrame(found)
    n_reference = int(manifest["is_reference"].sum())
    if n_reference != 1:
        raise ValueError(
            f"expected exactly one run in {reference_root}, found "
            f"{n_reference}"
        )

    # The reference sorts first within its dimensionality; everything else by name.
    manifest = manifest.sort_values(
        ["latent_dim", "is_reference", "directory"], ascending=[True, False, True]
    )
    manifest["replicate"] = (
        manifest[~manifest["is_reference"]].groupby("latent_dim").cumcount() + 1
    )
    manifest["label"] = [
        f"{dim:03d}_ref" if is_ref else f"{dim:03d}_rep{rep:.0f}"
        for dim, is_ref, rep in zip(
            manifest["latent_dim"], manifest["is_reference"], manifest["replicate"]
        )
    ]
    return manifest.set_index("label")[
        ["latent_dim", "replicate", "is_reference", "directory", "path"]
    ]