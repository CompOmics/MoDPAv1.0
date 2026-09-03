#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_cli():
    parser = argparse.ArgumentParser(description="Calculate and plot MoDPA matrix sparsity.")
    parser.add_argument("input", help="Input pandas pickle file (.pkl or .pkl.gz).")
    parser.add_argument("--prefix", default="modpa-sparsity", help="Prefix for output figures.")
    parser.add_argument("--nbins", type=int, default=25, help="Number of bins for histograms.")
    return parser.parse_args()


def plot_coverage(modpa_matrix, prefix, nbins):
    data = pd.read_pickle(modpa_matrix)
    print(f"Matrix shape: {data.shape}")
    
    data_binary = data.notna()

    sparsity = 100 * (1 - data_binary.to_numpy().mean())
    print(f"Sparsity = {sparsity:.2f}%")

    observations_per_ptm = data_binary.sum(axis=1)
    sns.histplot(
        observations_per_ptm, 
        bins=np.logspace(0, np.log10(observations_per_ptm.max() + 1), nbins), 
        # stat="percent"
        )
    plt.xlabel("No. Observations per PTM event")
    plt.ylabel("No. PTM events")
    plt.savefig(f"{prefix}-per-run.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{prefix}-per-run.svg", dpi=300, bbox_inches="tight")
    plt.close()

    observations_per_run = data_binary.sum(axis=0)
    sns.histplot(
        observations_per_run, 
        bins=np.logspace(0, np.log10(observations_per_run.max() + 1), nbins), 
        # stat="percent",
        color='orange'
        )
    plt.xlabel("No. Observations per MS run")
    plt.ylabel("No. MS runs")
    plt.savefig(f"{prefix}-per-ptm.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{prefix}-per-ptm.svg", dpi=300, bbox_inches="tight")
    plt.close()


def main():
    args = parse_cli()
    plot_coverage(
        args.input,
        args.prefix,
        args.nbins
    )


if __name__ == "__main__":
    main()
