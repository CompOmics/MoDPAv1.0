#!/usr/bin/env python
# coding: utf-8
"""Inspect the raw KL divergence carried by each latent dimension of a saved VAE.

The script loads a model through ``VAE_bilayer.load_vae()``, evaluates
``VAE_bilayer.kl_per_dimension()`` on the supplied dataset, prints the full
per-dimension KL table and summary statistics, and writes the table to CSV.

The reported KL values are the *raw* KL divergences in nats. They are not
floored by ``free_bits``. This is intentional: the raw distribution is the
quantity to inspect when deciding whether a free-bits threshold is sensible.
"""

import argparse
import os

import numpy as np
import pandas as pd

from vae import VAE_bilayer


def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect raw KL divergence per latent dimension for a saved VAE. "
            "The model folder must contain config.json and vae.weights.h5."
        )
    )
    parser.add_argument(
        "model",
        type=str,
        help="Path to the saved model folder.",
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help=(
            "Path to the input matrix (.pkl or .pkl.gz). If omitted, the script "
            "uses original.pkl.gz from the parent of the model folder, matching "
            "the grid-search output layout."
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1024,
        help="Batch size used while encoding (default: 1024).",
    )
    parser.add_argument(
        "--active-threshold",
        type=float,
        default=0.01,
        help=(
            "Raw KL in nats above which a latent dimension is counted as active "
            "(default: 0.01)."
        ),
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "CSV output path. If omitted, writes kl_per_dimension.csv inside "
            "the model folder."
        ),
    )
    parser.add_argument(
        "--sort",
        choices=("dimension", "kl"),
        default="dimension",
        help="Print rows by latent dimension or descending KL (default: dimension).",
    )
    return parser.parse_args()


def resolve_data_path(model_folder: str, requested: str | None) -> str:
    if requested is not None:
        return requested

    model_folder = os.path.abspath(model_folder)
    parent_candidate = os.path.join(os.path.dirname(model_folder), "original.pkl.gz")
    local_candidate = os.path.join(model_folder, "original.pkl.gz")

    if os.path.isfile(parent_candidate):
        return parent_candidate
    if os.path.isfile(local_candidate):
        return local_candidate

    raise FileNotFoundError(
        "Could not find original.pkl.gz automatically. Expected it in the parent "
        f"folder ({parent_candidate}) or model folder ({local_candidate}). "
        "Pass the matrix explicitly with --data."
    )


def main() -> None:
    args = parse_cli()

    if args.batch_size <= 0:
        raise ValueError(f"--batch-size must be > 0; got {args.batch_size}.")
    if args.active_threshold < 0:
        raise ValueError(
            f"--active-threshold must be >= 0; got {args.active_threshold}."
        )

    model_folder = os.path.abspath(args.model)
    data_path = os.path.abspath(resolve_data_path(model_folder, args.data))
    output_path = (
        os.path.abspath(args.output)
        if args.output is not None
        else os.path.join(model_folder, "kl_per_dimension.csv")
    )

    print(f"Loading model from : {model_folder}")
    vae = VAE_bilayer.load_vae(model_folder)

    print(f"Loss type          : {vae.loss_type}")
    print(f"Latent dimensions  : {vae.dense_mu.units}")
    print(f"Saved free_bits    : {vae.free_bits:g} nats/dimension")

    if not vae.uses_kl:
        print(
            "WARNING: this model was trained without a KL term. The encoder still "
            "defines mu/logvar and KL can be calculated, but these values were not "
            "regularised during training and should not be interpreted as a "
            "posterior-collapse diagnostic in the same way as for a +KL model."
        )

    print(f"Loading data from  : {data_path}")
    data = pd.read_pickle(data_path)
    x = data.fillna(0).values.astype(np.float64)

    if x.ndim != 2:
        raise ValueError(f"Expected a 2D input matrix; got shape {x.shape}.")
    if x.shape[1] != vae.original_dim:
        raise ValueError(
            f"Input has {x.shape[1]} features, but the model expects "
            f"original_dim={vae.original_dim}."
        )

    print(f"Input shape        : {x.shape}")
    print("Computing raw KL per latent dimension...")
    kl = vae.kl_per_dimension(
        x,
        batch_size=args.batch_size,
        threshold=args.active_threshold,
    )

    values = np.asarray(kl["kl_per_dim"], dtype=float)
    table = pd.DataFrame(
        {
            "latent_dimension": np.arange(1, len(values) + 1),
            "kl_nats": values,
            "active": values > args.active_threshold,
        }
    )

    if vae.free_bits > 0:
        table["below_free_bits"] = values < vae.free_bits

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    table.to_csv(output_path, index=False)

    q = kl["quantiles"]
    print("\nKL summary (raw, unfloored):")
    print(f"  total KL                    {kl['kl_total']:.6g} nats")
    print(f"  mean KL / dimension         {values.mean():.6g} nats")
    print(f"  active threshold            {args.active_threshold:.6g} nats")
    print(
        f"  active dimensions           {kl['n_active']} / {kl['latent_dim']} "
        f"({kl['n_active'] / kl['latent_dim']:.1%})"
    )
    print(
        f"  informative dimensions      {kl['n_informative']} / {kl['latent_dim']} "
        "(mu variance > 1e-3)"
    )
    print(f"  p05                          {q[0.05]:.6g}")
    print(f"  p10                          {q[0.10]:.6g}")
    print(f"  p25                          {q[0.25]:.6g}")
    print(f"  median                       {q[0.50]:.6g}")
    print(f"  p75                          {q[0.75]:.6g}")
    print(f"  p95                          {q[0.95]:.6g}")

    if vae.free_bits > 0:
        n_below = int(np.sum(values < vae.free_bits))
        print("\nFree-bits comparison:")
        print(f"  saved free_bits              {vae.free_bits:.6g} nats/dimension")
        print(
            f"  dimensions below free_bits  {n_below} / {len(values)} "
            f"({n_below / len(values):.1%})"
        )
        print(
            "  NOTE: values above are raw KL. free_bits floors the KL penalty "
            "during optimization; it does not alter the raw KL reported here."
        )

    shown = table
    if args.sort == "kl":
        shown = table.sort_values("kl_nats", ascending=False)

    print("\nPer-dimension KL:")
    print(shown.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print(f"\nWritten: {output_path}")


if __name__ == "__main__":
    main()
