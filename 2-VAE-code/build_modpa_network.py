#!/usr/bin/env python
# coding: utf-8
"""Filter a MoDPA association list, cluster it with Leiden, and write both tables.

This is the driver of the network step described in the root README. It takes
the signed-distance file written by ``calculate_sdcorr.py``, applies the score
and q-value filter, builds the network, partitions it with Leiden, annotates
the partition with protein-level information from a FASTA file, and writes two
CSV files into the output directory:

``<prefix>-filtered-distances.csv``
    one row per retained edge, with the annotation columns added by
    ``parse_network``

``<prefix>-Leiden-clusters.csv``
    one row per node, with its cluster, the protein annotations, and the name of
    its modification when ``--ptm-names`` is given

Both files are the input of steps 3, 4 and 5.

Clustering is weighted by default, using the absolute signed distance
correlation as the edge weight, which is the partition used throughout the
published analysis. ``--unweighted`` selects the alternative partition, in
which every retained edge counts equally, so the result depends on the topology
of the filtered network alone. The two partitions are not interchangeable and
the file name records which one was used.
"""

import argparse
import gzip
import os
import re
import sys
from datetime import date

import numpy as np
import pandas as pd

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
)

from modpa_network_utilities import (  # noqa: E402
    add_protein_annotations,
    build_clusters_df,
    build_graph,
    count_edges,
    count_nodes,
    parse_network,
    run_leiden,
    run_leiden_unweighted,
)

DEFAULT_FASTA = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "1-quant-pipeline-latest",
    "Human_2026_01_canonical.fasta.gz",
)


def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a MoDPA association list, cluster the resulting network "
            "with Leiden, and write the edge table and the cluster table as CSV."
        )
    )
    parser.add_argument(
        "associations", type=str, 
        help="Path to the association list written by calculate_sdcorr.py, that is <datetime>-<run name>-signed-distances.csv.gz."
    )
    parser.add_argument(
        "-o", "--output-dir", type=str, default=None,
        help="Directory for the two CSV files. Defaults to the directory of the association list, which is the model run folder."
    )
    parser.add_argument(
        "--prefix", type=str, default=date.today().isoformat(),
        help="Prefix of the two output file names. Defaults to <today>."
    )
    parser.add_argument(
        "--min-score", type=float, default=0.6,
        help="Minimum signed score of a retained edge. The test uses the signed score, so anti-correlated pairs are discarded. Default: 0.6."
    )
    parser.add_argument(
        "--qvalue", type=float, default=0.05,
        help="Maximum q-value of a retained edge. Default: 0.05.",
    )
    parser.add_argument("--adjacency-window", type=int, default=5, help="Flag PTMs closer than <window> amino acids.")
    parser.add_argument(
        "--resolution", type=float, default=2.0, help="Leiden resolution parameter. Default: 2.0.",
    )
    parser.add_argument("--n-iterations", type=int, default=2,help="Number of Leiden iterations. Default: 2.")
    parser.add_argument("--seed", type=int, default=42, help="Leiden random seed. Default: 42.")
    parser.add_argument(
        "--unweighted",
        action="store_true",
        help=(
            "Cluster the network without edge weights, so every retained edge "
            "counts equally. By default the absolute signed distance "
            "correlation is used as the edge weight."
        ),
    )
    parser.add_argument(
        "--fasta",
        type=str,
        default=DEFAULT_FASTA,
        help=(
            "FASTA file used to annotate the clustered nodes with their entry "
            "name and gene name. Default: the canonical human FASTA of step 1."
        ),
    )
    parser.add_argument(
        "--ptm-names",
        type=str,
        default=None,
        help=(
            "PTM list used to name the modification of every node, with the "
            "columns AA, unimod_id and ptm_name. It must be the list the MoDPA "
            "matrix was built with, because a node is named by its Unimod "
            "accession and its residue together; for the published run that is "
            "../1-quant-pipeline-latest/PTMs-of-interest-submission.csv. If "
            "omitted, the cluster table carries the Unimod accession alone and "
            "no PTM_name column."
        ),
    )
    parser.add_argument(
        "--min-proteins", type=int, default=20,
        help="Cluster size, in distinct proteins, reported as large enough for pathway over-representation in step 5. Default: 20."
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite the output files if they already exist.")
    return parser.parse_args()


def load_proteins(fasta_path: str) -> pd.DataFrame:
    """Read a UniProt FASTA file into one annotation row per accession.

    Headers are expected in UniProt format, `>db|Accession|EntryName
    Description`, with the gene name in the `GN=` field where it is present.
    """
    opener = gzip.open if fasta_path.endswith(".gz") else open
    proteins = []
    sequence_parts = []
    header = None

    def flush():
        if header is None:
            return
        _, uniacc, entry = header.split("|", 2)
        entry = entry.split(None, 1)[0]
        gene = re.search(r"(?:^| )GN=(\w+)", header)
        proteins.append(
            [
                uniacc,
                "".join(sequence_parts),
                entry,
                gene.group(1) if gene else np.nan,
                ">" + header,
            ]
        )

    with opener(fasta_path, "rt") as file:
        for line in file:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                sequence_parts = []
            elif header is not None:
                sequence_parts.append(line.strip())
    flush()

    if not proteins:
        raise ValueError(f"no FASTA records found in {fasta_path}")

    return pd.DataFrame(
        proteins, columns=["UniAcc", "Seq", "Entry", "Gene", "Header"]
    )


def load_ptm_names(ptm_names_path: str) -> dict[tuple[str, str], str]:
    """Read a PTM list into ``{(unimod_id, residue): ptm_name}``.

    The file is one of the `PTMs-of-interest*.csv` lists of step 1, with the
    columns `AA`, `unimod_id` and `ptm_name`. The key is the pair rather than
    the accession alone because `Generate_PTM_matrices.py` selects PTM events on
    residue and accession together, so the list also defines which residues each
    modification was allowed to occupy in the matrix.
    """
    ptm_names = pd.read_csv(ptm_names_path)
    missing = [
        column
        for column in ("AA", "unimod_id", "ptm_name")
        if column not in ptm_names.columns
    ]
    if missing:
        raise SystemExit(
            f"{ptm_names_path} lacks the column(s) {', '.join(missing)}; "
            "expected AA, unimod_id and ptm_name."
        )

    mapping = {}
    for residue, unimod_id, ptm_name in ptm_names[
        ["AA", "unimod_id", "ptm_name"]
    ].itertuples(index=False):
        key = (str(unimod_id).strip(), str(residue).strip())
        previous = mapping.get(key)
        if previous is not None and previous != ptm_name:
            raise SystemExit(
                f"{ptm_names_path} maps unimod_id {key[0]} on residue {key[1]} "
                f"to both {previous!r} and {ptm_name!r}."
            )
        mapping[key] = ptm_name
    return mapping


def add_ptm_names(
    clusters_df: pd.DataFrame,
    ptm_names: dict[tuple[str, str], str],
    ptm_names_path: str,
) -> pd.DataFrame:
    """Insert a `PTM_name` column, keyed on `UniModID` and `RES` together.

    An unmapped node is an error rather than a blank, because it means the PTM
    list is not the one the matrix was built with, which would mislabel part of
    the network without anything else failing.
    """
    keys = list(
        zip(
            clusters_df["UniModID"].astype(str).str.strip(),
            clusters_df["RES"].astype(str).str.strip(),
        )
    )
    unmapped = sorted({key for key in keys if key not in ptm_names})
    if unmapped:
        known_ids = {unimod_id for unimod_id, _ in ptm_names}
        unknown_ids = sorted(
            {unimod_id for unimod_id, _ in unmapped if unimod_id not in known_ids}
        )
        wrong_residue = [pair for pair in unmapped if pair[0] in known_ids]
        report = []
        if unknown_ids:
            report.append(
                "Unimod accession(s) absent from the list: " + ", ".join(unknown_ids)
            )
        if wrong_residue:
            report.append(
                "accession/residue combination(s) absent from the list: "
                + ", ".join(
                    f"{unimod_id} on {residue}" for unimod_id, residue in wrong_residue
                )
            )
        raise SystemExit(
            f"{ptm_names_path} does not name every PTM event in the network. "
            + "; ".join(report)
            + ". Pass the PTM list the MoDPA matrix was built with via --ptm-names."
        )

    annotated = clusters_df.copy()
    annotated.insert(
        annotated.columns.get_loc("UniModID") + 1,
        "PTM_name",
        [ptm_names[key] for key in keys],
    )
    return annotated


def main():
    args = parse_cli()

    output_dir = args.output_dir or os.path.dirname(
        os.path.abspath(args.associations)
    )
    os.makedirs(output_dir, exist_ok=True)
    prefix = args.prefix
    partition_name = "Leiden-unweighted-clusters" if args.unweighted else "Leiden-clusters"
    edges_path = os.path.join(output_dir, f"{prefix}-filtered-distances.csv")
    clusters_path = os.path.join(output_dir, f"{prefix}-{partition_name}.csv")

    existing = [p for p in (edges_path, clusters_path) if os.path.exists(p)]
    if existing and not args.overwrite:
        raise SystemExit(
            "Refusing to overwrite: "
            + ", ".join(existing)
            + ". Pass --overwrite, or change --prefix or --output-dir."
        )

    # Protein annotations are read before the network is built, so a missing or
    # unreadable FASTA file fails before the expensive steps rather than after.
    proteins = None
    proteins = load_proteins(args.fasta)
    print(f"Protein annotations: {len(proteins):,} records from {args.fasta}")
    ptm_names = None
    if args.ptm_names:
        ptm_names = load_ptm_names(args.ptm_names)
        print(
            f"PTM names: {len(set(ptm_names.values()))} modifications over "
            f"{len(ptm_names)} accession/residue combinations from "
            f"{args.ptm_names}"
        )

    edges = parse_network(
        args.associations,
        min_score=args.min_score,
        qvalue=args.qvalue,
        adjacency_window=args.adjacency_window,
    ).collect(engine="streaming")
    print(
        f"Filtered network (Score >= {args.min_score}, qvalue < {args.qvalue}): "
        f"{count_nodes(edges):,} nodes, {count_edges(edges):,} edges"
    )
    if not edges.height:
        raise SystemExit(
            "No edge passed the filter; nothing to cluster. Lower --min-score "
            "or raise --qvalue."
        )
    network = build_graph(edges)
    clustering = run_leiden_unweighted if args.unweighted else run_leiden
    membership = clustering(
        network,
        seed=args.seed,
        resolution=args.resolution,
        n_iterations=args.n_iterations,
    )
    clusters_df = build_clusters_df(membership)
    if ptm_names is not None:
        clusters_df = add_ptm_names(clusters_df, ptm_names, args.ptm_names)
    if proteins is not None:
        clusters_df = add_protein_annotations(clusters_df, proteins)

    # Both files are written only once every annotation step has succeeded, so a
    # PTM list or a FASTA file that does not match the network leaves no
    # half-written pair of outputs behind.
    edges.write_csv(edges_path)
    clusters_df.to_csv(clusters_path, index=False)

    weighting = "unweighted" if args.unweighted else "weighted by absolute SDCor"
    n_clusters = clusters_df["LeidenCluster"].nunique()
    proteins_per_cluster = clusters_df.groupby("LeidenCluster")["UniAcc"].nunique()
    n_large = int((proteins_per_cluster >= args.min_proteins).sum())
    print(
        f"Leiden ({weighting}, resolution {args.resolution}, "
        f"{args.n_iterations} iterations, seed {args.seed}): "
        f"{n_clusters} clusters, of which {n_large} hold at least "
        f"{args.min_proteins} distinct proteins"
    )
    if "PTM_name" in clusters_df.columns:
        print("\nPTM events per modification:")
        print(clusters_df["PTM_name"].value_counts().to_string())
    print("\nTen largest clusters, by number of PTM events:")
    print(clusters_df["LeidenCluster"].value_counts().head(10).to_string())
    print(f"\nWritten: {edges_path}")
    print(f"Written: {clusters_path}")


if __name__ == "__main__":
    main()
