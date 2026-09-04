"""
Core KEA3 pipeline used by the analysis and permutation notebooks.

The module only contains functionality needed by the notebook workflow:
cluster construction, HGNC/KEA3 symbol normalization, KEA3 querying,
deduplication of kinase aliases, optional rank plots, and aggregation of
per-cluster results.
"""

from __future__ import annotations

import json
import os
import time
from collections import defaultdict
from functools import lru_cache
from numbers import Real
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import requests

KEA3_URL = "https://maayanlab.cloud/kea3/api/enrich/"
REQUEST_SLEEP_SEC = 1.0
COMBINED_MEANRANK = "Combined_meanRank.csv"
COMBINED_TOPRANK = "Combined_topRank.csv"


@lru_cache(maxsize=None)
def _load_hgnc_maps(hgnc_complete_set_path: str):
    hgnc = pd.read_csv(hgnc_complete_set_path, sep="\t", low_memory=False)
    approved = set(hgnc["symbol"].dropna())
    eligible = set(
        hgnc.loc[hgnc["locus_group"] == "protein-coding gene", "symbol"].dropna()
    )

    claims = defaultdict(set)
    for _, row in hgnc.iterrows():
        official = row["symbol"]
        if official not in eligible:
            continue
        for col in ("alias_symbol", "prev_symbol"):
            value = row.get(col)
            if pd.isna(value):
                continue
            for alias in str(value).strip('"').split("|"):
                alias = alias.strip()
                if alias:
                    claims[alias].add(official)

    alias_map = {symbol: symbol for symbol in approved}
    ambiguous = set()
    for alias, owners in claims.items():
        if alias in approved:
            continue
        if len(owners) > 1:
            ambiguous.add(alias)
        else:
            alias_map[alias] = next(iter(owners))
    return alias_map, approved, ambiguous


@lru_cache(maxsize=1)
def load_kea3_kinase_alias_map() -> dict:
    """Map KEA3/Manning kinase nicknames to approved gene symbols."""
    import kinase_library

    table_path = os.path.join(
        os.path.dirname(kinase_library.__file__),
        "databases", "kinase_data", "kinome_information.tsv",
    )
    df = pd.read_csv(table_path, sep="\t")
    alias_map = {}
    for _, row in df.iterrows():
        gene = row["GENE_NAME"]
        if pd.isna(gene):
            continue
        alias_map[gene] = gene
        for col in ("KINASE", "MATRIX_NAME", "DISPLAY_NAME"):
            name = row.get(col)
            if pd.notna(name):
                alias_map[name] = gene
    return alias_map


def build_kea3_symbol_map(
    hgnc_complete_set_path: str,
    manual_overrides: dict | None = None,
    verbose: bool = True,
):
    """Build the symbol map for KEA3 kinase output.

    Kinase-library names take precedence over HGNC aliases because KEA3 uses
    Manning/KinBase-style kinase names that can collide with approved symbols.
    """
    hgnc_map, approved, ambiguous = _load_hgnc_maps(hgnc_complete_set_path)
    kinase_map = load_kea3_kinase_alias_map()
    symbol_map = {**hgnc_map, **kinase_map, **(manual_overrides or {})}

    rows = []
    for token in set(kinase_map) & set(hgnc_map):
        if kinase_map[token] != hgnc_map[token]:
            rows.append({
                "token": token,
                "kinase_library": kinase_map[token],
                "hgnc": hgnc_map[token],
                "resolved_to": symbol_map[token],
                "token_is_approved_symbol": token in approved,
            })
    audit = pd.DataFrame(
        rows,
        columns=["token", "kinase_library", "hgnc", "resolved_to", "token_is_approved_symbol"],
    ).sort_values("token").reset_index(drop=True)

    if verbose:
        print(
            f"[build_kea3_symbol_map] {len(symbol_map)} entries; "
            f"{len(ambiguous)} ambiguous HGNC aliases excluded; "
            f"{len(audit)} conflicts resolved in favour of kinase-library"
        )
    return symbol_map, audit


def normalize_query_genes(
    genes: Iterable[str],
    hgnc_complete_set_path: str,
    verbose: bool = True,
) -> list:
    """Normalize phosphoprotein symbols with HGNC only and remove duplicates."""
    alias_map, approved, _ = _load_hgnc_maps(hgnc_complete_set_path)
    normalized, unmapped = [], []
    for gene in genes:
        if not isinstance(gene, str) or not gene.strip():
            continue
        gene = gene.strip()
        if gene in approved:
            normalized.append(gene)
        elif gene in alias_map:
            normalized.append(alias_map[gene])
        else:
            normalized.append(gene)
            unmapped.append(gene)

    normalized = list(dict.fromkeys(normalized))
    if verbose and unmapped:
        print(
            f"[normalize_query_genes] {len(set(unmapped))} symbols absent from HGNC, "
            f"kept as-is: {sorted(set(unmapped))[:20]}"
        )
    return normalized


def query_kea3(genes: Iterable[str], query_name: str) -> dict:
    """Submit one gene set to the KEA3 API."""
    payload = {"gene_set": list(genes), "query_name": query_name}
    response = requests.post(KEA3_URL, data=json.dumps(payload), timeout=120)
    response.raise_for_status()
    time.sleep(REQUEST_SLEEP_SEC)
    return response.json()


def _integrated_df(results: dict, method: str) -> pd.DataFrame:
    key = f"Integrated--{method}"
    if key not in results:
        raise KeyError(f"'{key}' not found in KEA3 results. Available: {list(results)}")
    df = pd.DataFrame(results[key]).copy()
    df["Rank"] = pd.to_numeric(df["Rank"], errors="coerce")
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    return df.sort_values("Rank").reset_index(drop=True)


def _fmt(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def _merge_library_entries(entries: Iterable[str], score_mode: str):
    best = {}
    for entry in entries:
        for item in str(entry or "").split(";"):
            if not item:
                continue
            library, _, value = item.partition(",")
            try:
                value = float(value)
            except ValueError:
                continue
            library = library.strip()
            if library not in best or value < best[library]:
                best[library] = value

    if not best:
        return "", float("nan"), 0
    values = list(best.values())
    score = min(values) if score_mode == "min" else sum(values) / len(values)
    merged = ";".join(f"{library},{_fmt(value)}" for library, value in best.items())
    return merged, score, len(best)


def dedup_by_normalized_symbol(
    df: pd.DataFrame,
    symbol_map: dict,
    score_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge KEA3 rows that normalize to the same kinase symbol."""
    df = df.copy()
    df["symbol_norm"] = df["TF"].map(lambda gene: symbol_map.get(gene, gene))
    df = df.sort_values("Rank").reset_index(drop=True)

    grouped = df.groupby("symbol_norm")[["TF", "Rank"]].agg(list)
    collapse_log = grouped[grouped["TF"].map(len) > 1].reset_index()
    collapse_log = collapse_log.rename(columns={"TF": "source_names", "Rank": "source_ranks"})

    rows = []
    for symbol, group in df.groupby("symbol_norm", sort=False):
        row = group.iloc[0].to_dict()
        row["TF"] = ";".join(dict.fromkeys(group["TF"].astype(str)))
        row["symbol_norm"] = symbol
        row["Library"], row["Score"], row["n_libraries"] = _merge_library_entries(
            group["Library"], score_mode
        )
        genes = dict.fromkeys(
            gene.strip()
            for entry in group["Overlapping_Genes"].fillna("")
            for gene in str(entry).split(",")
            if gene.strip()
        )
        row["Overlapping_Genes"] = ",".join(genes)
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("Score").reset_index(drop=True)
    out["Rank"] = range(1, len(out) + 1)
    return out, collapse_log


def plot_bar(df: pd.DataFrame, title: str, outdir: str, n: int = 10) -> str:
    """Save a compact horizontal score plot for the top kinases."""
    os.makedirs(outdir, exist_ok=True)
    top = df.head(n)
    plt.figure(figsize=(6, 4))
    plt.barh(top["symbol_norm"][::-1], top["Score"][::-1])
    plt.xlabel("Score")
    plt.title(title)
    plt.tight_layout()
    outfile = os.path.join(outdir, f"{title}.png")
    plt.savefig(outfile, dpi=150)
    plt.close()
    return outfile


def _cluster_name(label, cluster_offset: int = 1) -> str:
    if isinstance(label, Real) and not pd.isna(label) and float(label).is_integer():
        label = int(label) + cluster_offset
    return f"Cluster_{label}"


def cluster_membership(
    data: str | pd.DataFrame,
    label_col: str = "LeidenCluster",
    gene_col: str = "Gene",
    minsize: int = 20,
    maxsize: int = 100,
    cluster_offset: int = 1,
    hgnc_complete_set_path: str | None = None,
) -> dict:
    """Return {cluster_name: unique submitted gene set} after the KEA3 size filter."""
    df = pd.read_csv(data) if isinstance(data, str) else data.copy()
    for col in (label_col, gene_col):
        if col not in df.columns:
            raise KeyError(f"'{col}' not found. Available columns: {list(df.columns)}")

    membership = {}
    for label, group in df[df[label_col].notna()].groupby(label_col):
        genes = {
            gene.strip()
            for gene in group[gene_col]
            if isinstance(gene, str) and gene.strip()
        }
        if hgnc_complete_set_path:
            genes = set(normalize_query_genes(genes, hgnc_complete_set_path, verbose=False))
        if minsize <= len(genes) <= maxsize:
            membership[_cluster_name(label, cluster_offset)] = genes
    return membership


def run_cluster(
    genes: Iterable[str],
    name: str,
    outdir: str,
    hgnc_complete_set_path: str = "hgnc_complete_set.txt",
    symbol_map: dict | None = None,
    n_top: int = 10,
    make_plots: bool = True,
) -> dict:
    """Query and save KEA3 results for one cluster."""
    os.makedirs(outdir, exist_ok=True)
    genes = sorted(normalize_query_genes(genes, hgnc_complete_set_path))
    if not genes:
        raise ValueError(f"{name} contains no valid gene symbols")

    if symbol_map is None:
        symbol_map, _ = build_kea3_symbol_map(hgnc_complete_set_path)

    with open(os.path.join(outdir, f"{name}_gene_list.txt"), "w") as handle:
        handle.write("\n".join(genes) + "\n")

    results = query_kea3(genes, query_name=name)
    meanrank = _integrated_df(results, "meanRank")
    toprank = _integrated_df(results, "topRank")

    meanrank, meanrank_log = dedup_by_normalized_symbol(meanrank, symbol_map, "mean")
    toprank, toprank_log = dedup_by_normalized_symbol(toprank, symbol_map, "min")
    for df in (meanrank, toprank):
        if "Query Name" not in df.columns:
            df["Query Name"] = name

    for label, log in (("MeanRank", meanrank_log), ("TopRank", toprank_log)):
        if len(log):
            log.to_csv(
                os.path.join(outdir, f"{name}_{label}_merges.tsv"),
                sep="\t",
                index=False,
            )

    meanrank.to_csv(os.path.join(outdir, f"{name}_MeanRank.tsv"), sep="\t", index=False)
    toprank.to_csv(os.path.join(outdir, f"{name}_TopRank.tsv"), sep="\t", index=False)
    if make_plots:
        plot_bar(meanrank, f"{name}_MeanRank", outdir, n_top)
        plot_bar(toprank, f"{name}_TopRank", outdir, n_top)
    return {"meanRank": meanrank, "topRank": toprank}


def run_clusters(
    membership: dict,
    outdir: str,
    hgnc_complete_set_path: str = "hgnc_complete_set.txt",
    n_top: int = 10,
    make_plots: bool = True,
) -> dict:
    """Run KEA3 for a prepared cluster membership dictionary and combine results."""
    if not membership:
        raise ValueError("No clusters passed the size filter")
    os.makedirs(outdir, exist_ok=True)

    symbol_map, audit = build_kea3_symbol_map(hgnc_complete_set_path, verbose=False)
    if len(audit):
        audit.to_csv(os.path.join(outdir, "symbol_conflicts.tsv"), sep="\t", index=False)

    results = {}
    for name, genes in membership.items():
        results[name] = run_cluster(
            genes,
            name=name,
            outdir=os.path.join(outdir, name),
            hgnc_complete_set_path=hgnc_complete_set_path,
            symbol_map=symbol_map,
            n_top=n_top,
            make_plots=make_plots,
        )

    meanrank = pd.concat([result["meanRank"] for result in results.values()], ignore_index=True)
    toprank = pd.concat([result["topRank"] for result in results.values()], ignore_index=True)
    meanrank.to_csv(os.path.join(outdir, COMBINED_MEANRANK), index=False)
    toprank.to_csv(os.path.join(outdir, COMBINED_TOPRANK), index=False)
    return {"clusters": results, "meanRank": meanrank, "topRank": toprank}


def run_kea_pipeline(
    path: str,
    outdir: str,
    clustering_column: str = "LeidenCluster",
    gene_column: str = "Gene",
    minsize: int = 20,
    maxsize: int = 100,
    cluster_offset: int = 1,
    hgnc_complete_set_path: str = "hgnc_complete_set.txt",
    n_top: int = 10,
    make_plots: bool = True,
) -> dict:
    """Build eligible clusters from a CSV file, run KEA3, and save combined tables."""
    membership = cluster_membership(
        path,
        label_col=clustering_column,
        gene_col=gene_column,
        minsize=minsize,
        maxsize=maxsize,
        cluster_offset=cluster_offset,
        hgnc_complete_set_path=hgnc_complete_set_path,
    )
    result = run_clusters(
        membership,
        outdir=outdir,
        hgnc_complete_set_path=hgnc_complete_set_path,
        n_top=n_top,
        make_plots=make_plots,
    )
    result["membership"] = membership
    return result
