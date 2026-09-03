"""
Label-shuffle null distribution for consensus-based KEA3 analysis.

Cluster labels are shuffled at phosphosite level. Every realization is passed
through the same KEA3 and consensus code as the observed clustering.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from itertools import combinations

import numpy as np
import pandas as pd

from kea3_consensus import consensus_table, cross_cluster_frequency, top_calls_per_cluster
from kea3_pipeline import COMBINED_MEANRANK, cluster_membership, run_clusters

STATISTIC_TAILS = {
    "evidence_concordance": "right",
    "mean_jaccard": "left",
    "background_fraction": "left",
    "frac_unique_kinases": "right",
    "mean_call_strength": "right",
    "self_membership_delta": "right",
}


def _two_way_residual(values: pd.Series, kinase: pd.Series, cluster: pd.Series) -> np.ndarray:
    """Remove kinase and cluster fixed effects by alternating demeaning."""
    residual = values.to_numpy(dtype=float) - float(values.mean())
    work = pd.DataFrame({"residual": residual, "kinase": kinase.to_numpy(), "cluster": cluster.to_numpy()})
    for _ in range(100):
        previous = work["residual"].to_numpy().copy()
        work["residual"] -= work.groupby("kinase")["residual"].transform("mean")
        work["residual"] -= work.groupby("cluster")["residual"].transform("mean")
        if np.max(np.abs(work["residual"].to_numpy() - previous)) < 1e-10:
            break
    return work["residual"].to_numpy()


def evidence_concordance(
    consensus: pd.DataFrame,
    min_substrate: int = 3,
    min_interaction: int = 6,
) -> float:
    """Correlation between substrate and interaction evidence after two-way residualisation."""
    data = consensus[
        (consensus["n_substrate"] >= min_substrate)
        & (consensus["n_interaction"] >= min_interaction)
    ].copy()
    if len(data) < 50 or data["Query Name"].nunique() < 3:
        return float("nan")

    data["substrate_rate"] = data["n_top_substrate"] / data["n_substrate"]
    data["interaction_rate"] = data["n_top_interaction"] / data["n_interaction"]
    substrate = _two_way_residual(data["substrate_rate"], data["symbol_norm"], data["Query Name"])
    interaction = _two_way_residual(data["interaction_rate"], data["symbol_norm"], data["Query Name"])
    if np.std(substrate) == 0 or np.std(interaction) == 0:
        return float("nan")
    return float(np.corrcoef(substrate, interaction)[0, 1])


def self_membership_delta(consensus: pd.DataFrame, membership: dict) -> float:
    """Median consensus-rank improvement for kinases in their own phosphoprotein cluster."""
    if not membership:
        return float("nan")

    ranks = consensus.pivot_table(
        index="symbol_norm",
        columns="Query Name",
        values="consensus_rank",
        aggfunc="first",
    )
    deltas = []
    for cluster, genes in membership.items():
        if cluster not in ranks.columns:
            continue
        for kinase in genes & set(ranks.index):
            own = ranks.loc[kinase, cluster]
            others = ranks.loc[kinase].drop(cluster).dropna()
            if pd.isna(own) or len(others) < 2:
                continue
            deltas.append(float(others.median() - own))
    return float(np.median(deltas)) if deltas else float("nan")


def realization_statistics(
    consensus: pd.DataFrame,
    top_n: int = 12,
    min_n_top: int = 4,
    background_frac: float = 0.5,
    membership: dict | None = None,
) -> dict:
    """Compute all summary statistics for one observed or shuffled partition."""
    cluster_names = list(dict.fromkeys(consensus["Query Name"]))
    n_clusters = len(cluster_names)
    top = top_calls_per_cluster(consensus, n=top_n, min_n_top=min_n_top)
    freq = cross_cluster_frequency(consensus, min_n_top=min_n_top)
    sets = {
        cluster: set(top.loc[top["Query Name"] == cluster, "symbol_norm"])
        for cluster in cluster_names
    }

    jaccard = []
    for a, b in combinations(sets.values(), 2):
        union = a | b
        if union:
            jaccard.append(len(a & b) / len(union))

    background = set(
        freq.loc[
            freq["n_clusters_called"] >= background_frac * n_clusters,
            "symbol_norm",
        ]
    )
    background_fractions = [
        len(genes & background) / len(genes)
        for genes in sets.values()
        if genes
    ]

    return {
        "n_clusters": n_clusters,
        "n_called_kinases": len(freq),
        "mean_jaccard": float(np.mean(jaccard)) if jaccard else np.nan,
        "frac_unique_kinases": float((freq["n_clusters_called"] == 1).mean()) if len(freq) else np.nan,
        "background_fraction": float(np.mean(background_fractions)) if background_fractions else np.nan,
        "mean_call_strength": float(top["n_top"].mean()) if len(top) else np.nan,
        "evidence_concordance": evidence_concordance(consensus),
        "self_membership_delta": self_membership_delta(consensus, membership or {}),
    }


def _statistics_for_partition(
    membership: dict,
    outdir: str,
    hgnc_complete_set_path: str,
    top_n: int,
    min_n_top: int,
    background_frac: float,
) -> dict:
    combined_path = os.path.join(outdir, COMBINED_MEANRANK)
    if not os.path.exists(combined_path):
        run_clusters(
            membership,
            outdir=outdir,
            hgnc_complete_set_path=hgnc_complete_set_path,
            make_plots=False,
        )
    consensus = consensus_table(combined_path)
    return realization_statistics(
        consensus,
        top_n=top_n,
        min_n_top=min_n_top,
        background_frac=background_frac,
        membership=membership,
    )


def _empirical_p(observed: float, null_values) -> float:
    """Two-tailed finite-sample empirical p-value."""
    values = np.asarray(null_values, dtype=float)
    values = values[~np.isnan(values)]
    if len(values) == 0 or np.isnan(observed):
        return float("nan")

    p_left = (1 + np.sum(values <= observed)) / (1 + len(values))
    p_right = (1 + np.sum(values >= observed)) / (1 + len(values))
    return float(min(1.0, 2 * min(p_left, p_right)))


def _cache_config(path: str, config: dict) -> None:
    config_path = os.path.join(path, "permutation_config.json")
    null_path = os.path.join(path, "null_statistics.csv")
    if os.path.exists(config_path):
        with open(config_path) as handle:
            cached = json.load(handle)
        if cached != config:
            raise ValueError(
                "This workdir contains permutation results generated with different settings. "
                "Use a new workdir or remove the existing permutation cache."
            )
    elif os.path.exists(null_path):
        raise ValueError(
            "This workdir contains a legacy permutation cache without configuration metadata. "
            "Use a new workdir or remove the old cache before running the cleaned pipeline."
        )
    else:
        with open(config_path, "w") as handle:
            json.dump(config, handle, indent=2, sort_keys=True)


def run_label_shuffle_null(
    clustered_phosphosites_path: str,
    n_permutations: int = 100,
    minsize: int = 40,
    maxsize: int = 100,
    label_col: str = "LeidenCluster",
    gene_col: str = "Gene",
    workdir: str = "kea3-comparison/perm_null",
    top_n: int = 12,
    min_n_top: int = 4,
    background_frac: float = 0.5,
    cluster_offset: int = 1,
    hgnc_complete_set_path: str = "hgnc_complete_set.txt",
    observed_meanrank_path: str | None = None,
    max_cluster_difference: int = 5,
    seed: int = 0,
    keep_kea3_output: bool = False,
    progress: bool = True,
):
    """Build a label-shuffle null and compare the observed clustering with it."""
    if minsize < 20:
        raise ValueError("KEA3 requires clusters with at least 20 unique genes; set minsize >= 20.")

    df = pd.read_csv(clustered_phosphosites_path)
    for col in (label_col, gene_col):
        if col not in df.columns:
            raise KeyError(f"'{col}' not found. Available columns: {list(df.columns)}")
    df = df[df[label_col].notna()].copy()
    os.makedirs(workdir, exist_ok=True)

    with open(clustered_phosphosites_path, "rb") as handle:
        input_hash = hashlib.sha256(handle.read()).hexdigest()

    config = {
        "cache_version": 3,
        "input_path": os.path.abspath(clustered_phosphosites_path),
        "input_sha256": input_hash,
        "minsize": minsize,
        "maxsize": maxsize,
        "label_col": label_col,
        "gene_col": gene_col,
        "top_n": top_n,
        "min_n_top": min_n_top,
        "background_frac": background_frac,
        "cluster_offset": cluster_offset,
        "hgnc_complete_set_path": os.path.abspath(hgnc_complete_set_path),
        "seed": seed,
        "observed_meanrank_path": os.path.abspath(observed_meanrank_path) if observed_meanrank_path else None,
        "max_cluster_difference": max_cluster_difference,
    }
    _cache_config(workdir, config)

    observed_membership = cluster_membership(
        df,
        label_col=label_col,
        gene_col=gene_col,
        minsize=minsize,
        maxsize=maxsize,
        cluster_offset=cluster_offset,
        hgnc_complete_set_path=hgnc_complete_set_path,
    )
    n_observed = len(observed_membership)
    if n_observed < 2:
        raise ValueError(f"Only {n_observed} observed clusters passed the size filter")

    if observed_meanrank_path:
        observed_consensus = consensus_table(observed_meanrank_path)
        observed = realization_statistics(
            observed_consensus,
            top_n=top_n,
            min_n_top=min_n_top,
            background_frac=background_frac,
            membership=observed_membership,
        )
        if observed["n_clusters"] != n_observed:
            raise ValueError(
                f"Observed MeanRank contains {observed['n_clusters']} clusters, "
                f"but {n_observed} clusters pass the current size filter"
            )
    else:
        observed = _statistics_for_partition(
            observed_membership,
            os.path.join(workdir, "observed"),
            hgnc_complete_set_path,
            top_n,
            min_n_top,
            background_frac,
        )

    observed["n_clusters_eligible"] = n_observed
    observed["n_clusters_analyzed"] = n_observed

    if progress:
        print(f"Observed: {n_observed} clusters eligible and analyzed")
        for key, value in observed.items():
            print(f"  {key:24s} {value}")

    labels = df[label_col].to_numpy()
    null_path = os.path.join(workdir, "null_statistics.csv")
    rows = []
    if os.path.exists(null_path):
        cached = pd.read_csv(null_path)
        done = {int(row.permutation): row._asdict() for row in cached.itertuples(index=False)}
        if progress:
            print(f"Resuming: {len(done)} permutations already computed")
    else:
        done = {}

    for i in range(n_permutations):
        if i in done:
            rows.append(done[i])
            continue

        rng = np.random.default_rng(seed + i)
        shuffled = df.copy()
        shuffled["_perm_label"] = rng.permutation(labels)
        membership = cluster_membership(
            shuffled,
            label_col="_perm_label",
            gene_col=gene_col,
            minsize=minsize,
            maxsize=maxsize,
            cluster_offset=cluster_offset,
            hgnc_complete_set_path=hgnc_complete_set_path,
        )
        n_clusters_eligible = len(membership)
        cluster_difference = n_observed - n_clusters_eligible
        if cluster_difference > max_cluster_difference:
            if progress:
                print(
                    f"  permutation {i + 1}/{n_permutations}: skipped "
                    f"({n_clusters_eligible} vs {n_observed} eligible clusters; "
                    f"difference {cluster_difference} > {max_cluster_difference})"
                )
            continue
        if len(membership) > n_observed:
            keep = rng.choice(list(membership), size=n_observed, replace=False)
            membership = {name: membership[name] for name in keep}

        n_clusters_analyzed = len(membership)
        perm_dir = os.path.join(workdir, f"perm_{i:04d}")
        stats = _statistics_for_partition(
            membership,
            perm_dir,
            hgnc_complete_set_path,
            top_n,
            min_n_top,
            background_frac,
        )
        stats["permutation"] = i
        stats["n_clusters_eligible"] = n_clusters_eligible
        stats["n_clusters_analyzed"] = n_clusters_analyzed
        rows.append(stats)

        pd.DataFrame(rows).sort_values("permutation").to_csv(
            null_path,
            index=False,
            float_format="%.17g",
        )
        if not keep_kea3_output:
            for entry in os.listdir(perm_dir):
                path = os.path.join(perm_dir, entry)
                if os.path.isdir(path):
                    shutil.rmtree(path)

        if progress:
            print(
                f"  permutation {i + 1}/{n_permutations}: "
                f"clusters={n_clusters_analyzed}/{n_clusters_eligible} analyzed/eligible, "
                f"concordance={stats['evidence_concordance']:.3f}, "
                f"jaccard={stats['mean_jaccard']:.3f}, "
                f"unique={stats['frac_unique_kinases']:.3f}"
            )

    null_df = pd.read_csv(null_path) if os.path.exists(null_path) else pd.DataFrame(rows)
    if "n_clusters_analyzed" not in null_df.columns and "n_clusters" in null_df.columns:
        null_df["n_clusters_analyzed"] = null_df["n_clusters"]
    if "n_clusters_eligible" not in null_df.columns:
        null_df["n_clusters_eligible"] = np.nan

    columns = [
        "permutation",
        "n_clusters_eligible",
        "n_clusters_analyzed",
        "n_clusters",
        "n_called_kinases",
        *STATISTIC_TAILS,
    ]
    null_df = null_df.reindex(columns=[col for col in columns if col in null_df.columns])
    if "permutation" in null_df.columns:
        null_df = null_df[null_df["permutation"] < n_permutations].copy()
        null_df["permutation"] = null_df["permutation"].astype(int)
        null_df = null_df.sort_values("permutation").reset_index(drop=True)
    for statistic in STATISTIC_TAILS:
        if statistic in null_df.columns:
            null_df[statistic] = null_df[statistic].astype(float)

    pvals = {
        statistic: _empirical_p(observed[statistic], null_df[statistic])
        for statistic in STATISTIC_TAILS
        if statistic in null_df.columns
    }

    if progress:
        print("\n=== Label-shuffle null ===")
        print(f"{'statistic':24s} {'observed':>10s} {'null mean':>10s} {'null sd':>9s} {'p (2-sided)':>12s}")
        for statistic in STATISTIC_TAILS:
            if statistic not in null_df.columns:
                continue
            values = null_df[statistic].dropna()
            if not len(values):
                continue
            print(
                f"{statistic:24s} {observed[statistic]:10.4f} "
                f"{values.mean():10.4f} {values.std():9.4f} {pvals[statistic]:12.4f}"
            )
        print(f"\nn_permutations used: {len(null_df)}")

    return null_df, observed, pvals


def plot_null(
    null_df: pd.DataFrame,
    observed: dict,
    pvals: dict,
    savepath: str | None = None,
):
    """Plot the null distributions with the observed statistic marked."""
    import matplotlib.pyplot as plt

    statistics = [
        statistic
        for statistic in STATISTIC_TAILS
        if statistic in null_df.columns and null_df[statistic].notna().any()
    ]
    if not statistics:
        raise ValueError(
            "No valid permutation statistics are available to plot. "
            "Check null_df and the permutation skip messages."
        )

    ncols = 3
    nrows = int(np.ceil(len(statistics) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, statistic in zip(axes, statistics):
        values = null_df[statistic].dropna()
        ax.hist(values, bins=25, alpha=0.75, label="Shuffled")
        ax.axvline(observed[statistic], lw=2, label=f"Observed = {observed[statistic]:.3f}")
        ax.set_title(
            f"{statistic}\ntwo-tailed p = "
            f"{pvals.get(statistic, float('nan')):.4f}",
            fontsize=10,
        )
        ax.set_xlabel(statistic)
        ax.set_ylabel("Permutations")
        ax.legend(fontsize=8)
    for ax in axes[len(statistics):]:
        ax.axis("off")

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=150)
    return fig
