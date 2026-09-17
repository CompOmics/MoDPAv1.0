#!/usr/bin/env python
# coding: utf-8
import os, re, argparse
import numpy as np
import pandas as pd
from vae import VAE_bilayer

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Encode an input matrix with a trained VAE and save the latent space."
    )
    p.add_argument('model',  type=str, help="Path to the model folder (contains config.json + vae.weights.h5)")
    p.add_argument('data',   type=str, help="Path to input matrix (.pkl or .pkl.gz)")
    p.add_argument('--output', type=str, default=None,
                   help="Output path for the latent space pickle (default: Latent-space.pkl.gz inside the model folder)")
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_cli()

    # Default output path mirrors the convention used by New-VAE-train.py
    if args.output is None:
        args.output = os.path.join(args.model, "Latent-space.pkl.gz")

    print(f"Loading model from : {args.model}")
    vae = VAE_bilayer.load_vae(args.model)

    print(f"Loading data from  : {args.data}")
    data = pd.read_pickle(args.data)
    data.fillna(0, inplace=True)
    print(f"Input shape        : {data.shape}  |  min={data.values.min():.4f}  max={data.values.max():.4f}")

    # Strip bracket suffixes from the index (e.g. "PTM[mod]" → "PTMmod"),
    # matching the convention used in New-VAE-train.py.
    index_clean = [''.join(re.split(r'[][]', x)[:2]) for x in data.index]

    print("Encoding...")
    _batch = 256
    _vals  = data.values
    _mu_chunks = []
    for _i in range(0, len(_vals), _batch):
        _mu, _ = vae.encode(_vals[_i:_i+_batch])
        _mu_chunks.append(_mu.numpy())
    del _vals
    latent = pd.DataFrame(np.concatenate(_mu_chunks, axis=0), index=index_clean)

    print(f"Latent shape       : {latent.shape}")
    print(f"Saving to          : {args.output}")
    latent.to_pickle(args.output, compression='gzip')

    print("Reconstructing...")
    _rec_chunks = []
    for _i in range(0, len(latent), _batch):
        _rec = vae.decode(latent.values[_i:_i+_batch])
        _rec_chunks.append(_rec.numpy())
    reconstructed = pd.DataFrame(np.concatenate(_rec_chunks, axis=0), index=index_clean)

    rec_output = os.path.join(args.model, "Reconstruction.pkl.gz")
    print(f"Reconstr. shape    : {reconstructed.shape}")
    print(f"Saving to          : {rec_output}")
    reconstructed.to_pickle(rec_output, compression='gzip')
    print("Done.")
