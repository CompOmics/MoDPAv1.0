"""
Consensus scoring for combined KEA3 MeanRank results.

The score asks whether independent KEA3 libraries agree on the same kinase
within a cluster, rather than relying on KEA3's integrated MeanRank ordering.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

SUBSTRATE_LIBRARIES = {
    "The_Kinase_Library", "PhosDAll", "ChengKSIN", "PTMsigDB",
}
INTERACTION_LIBRARIES = {
    "prePPI", "STRING", "STRING.bind", "mentha", "HIPPIE", "BioGRID",
    "MINT", "ChengPPI",
}


def load_combined(data: str | pd.DataFrame) -> pd.DataFrame:
    """Load a combined KEA3 table from disk or copy an existing DataFrame."""
    df = pd.read_csv(data) if isinstance(data, str) else data.copy()
    required = {"Query Name", "Rank", "Library", "symbol_norm"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Combined KEA3 table is missing columns: {sorted(missing)}")
    df["Rank"] = pd.to_numeric(df["Rank"], errors="coerce")
    df["_cluster_order"] = pd.to_numeric(
        df["Query Name"].astype(str).str.extract(r"(\d+)")[0],
        errors="coerce",
    )
    return df.sort_values(["_cluster_order", "Query Name", "Rank"]).reset_index(drop=True)


def parse_library_ranks(df: pd.DataFrame) -> pd.DataFrame:
    """Explode Library into one row per cluster, kinase, library and rank."""
    records = []
    for cluster, symbol, entry in zip(
        df["Query Name"], df["symbol_norm"], df["Library"].fillna("")
    ):
        for item in str(entry).split(";"):
            if not item:
                continue
            library, _, rank = item.partition(",")
            try:
                records.append((cluster, symbol, library.strip(), float(rank)))
            except ValueError:
                continue
    return pd.DataFrame(
        records,
        columns=["Query Name", "symbol_norm", "library", "library_rank"],
    )


def consensus_table(
    meanrank: str | pd.DataFrame,
    top_frac: float = 0.10,
    min_libraries: int = 3,
) -> pd.DataFrame:
    """Score every cluster/kinase pair by cross-library agreement."""
    df = load_combined(meanrank)
    long = parse_library_ranks(df)
    if long.empty:
        raise ValueError("No per-library ranks could be parsed from the MeanRank table")

    library_sizes = long.groupby("library")["library_rank"].max()
    long["pct"] = long["library_rank"] / long["library"].map(library_sizes)
    long["is_top"] = long["pct"] <= top_frac
    long["kind"] = np.where(
        long["library"].isin(SUBSTRATE_LIBRARIES),
        "substrate",
        np.where(long["library"].isin(INTERACTION_LIBRARIES), "interaction", "other"),
    )

    out = (
        long.groupby(["Query Name", "symbol_norm"])
        .agg(
            n_libraries=("library", "size"),
            n_top=("is_top", "sum"),
            median_pct=("pct", "median"),
            min_pct=("pct", "min"),
        )
        .reset_index()
    )

    for kind in ("substrate", "interaction"):
        counts = (
            long[long["kind"] == kind]
            .groupby(["Query Name", "symbol_norm"])
            .agg(
                **{
                    f"n_{kind}": ("library", "size"),
                    f"n_top_{kind}": ("is_top", "sum"),
                }
            )
            .reset_index()
        )
        out = out.merge(counts, on=["Query Name", "symbol_norm"], how="left")

    count_cols = ["n_substrate", "n_top_substrate", "n_interaction", "n_top_interaction"]
    out[count_cols] = out[count_cols].fillna(0).astype(int)
    out["frac_top"] = out["n_top"] / out["n_libraries"]
    out["low_coverage"] = out["n_libraries"] < min_libraries
    out["evidence"] = np.select(
        [
            (out["n_top_substrate"] > 0) & (out["n_top_interaction"] > 0),
            out["n_top_substrate"] > 0,
            out["n_top_interaction"] > 0,
        ],
        ["both", "substrate only", "interaction only"],
        default="none",
    )

    ranks = df[["Query Name", "symbol_norm", "Rank", "_cluster_order"]].rename(
        columns={"Rank": "meanrank_rank"}
    )
    out = out.merge(ranks, on=["Query Name", "symbol_norm"], how="left")
    out = out.sort_values(
        ["_cluster_order", "Query Name", "n_top", "median_pct"],
        ascending=[True, True, False, True],
    )
    out["consensus_rank"] = out.groupby("Query Name").cumcount() + 1
    return out.drop(columns="_cluster_order").reset_index(drop=True)


def cross_cluster_frequency(
    consensus: pd.DataFrame,
    min_n_top: int = 3,
) -> pd.DataFrame:
    """Count in how many clusters each consensus kinase is called."""
    n_clusters = consensus["Query Name"].nunique()
    if n_clusters == 0:
        return pd.DataFrame(
            columns=["symbol_norm", "n_clusters_called", "mean_n_top", "mean_n_libraries", "specificity"]
        )

    called = consensus[consensus["n_top"] >= min_n_top]
    freq = (
        called.groupby("symbol_norm")
        .agg(
            n_clusters_called=("Query Name", "nunique"),
            mean_n_top=("n_top", "mean"),
            mean_n_libraries=("n_libraries", "mean"),
        )
        .reset_index()
    )
    if n_clusters == 1:
        freq["specificity"] = 1.0
    else:
        freq["specificity"] = (n_clusters - freq["n_clusters_called"]) / (n_clusters - 1)
    return freq.sort_values(
        ["n_clusters_called", "mean_n_top"],
        ascending=[False, False],
    ).reset_index(drop=True)


def top_calls_per_cluster(
    consensus: pd.DataFrame,
    n: int = 10,
    min_n_top: int = 3,
    exclude: set | None = None,
) -> pd.DataFrame:
    """Return the top consensus calls in each cluster."""
    selected = consensus[consensus["n_top"] >= min_n_top]
    if exclude:
        selected = selected[~selected["symbol_norm"].isin(exclude)]
    return selected.groupby("Query Name", sort=False).head(n).reset_index(drop=True)


def cluster_similarity(
    consensus: pd.DataFrame,
    n: int = 20,
    min_n_top: int = 3,
) -> pd.DataFrame:
    """Pairwise Jaccard overlap of each cluster's top consensus calls."""
    top = top_calls_per_cluster(consensus, n=n, min_n_top=min_n_top)
    names = list(dict.fromkeys(consensus["Query Name"]))
    sets = {
        name: set(top.loc[top["Query Name"] == name, "symbol_norm"])
        for name in names
    }
    matrix = pd.DataFrame(np.nan, index=names, columns=names, dtype=float)
    for name in names:
        matrix.loc[name, name] = 1.0 if sets[name] else np.nan
    for a, b in combinations(names, 2):
        union = sets[a] | sets[b]
        value = len(sets[a] & sets[b]) / len(union) if union else np.nan
        matrix.loc[a, b] = matrix.loc[b, a] = value
    return matrix
