#!/usr/bin/env python
# coding: utf-8
"""
Evaluate every trained VAE in a folder and rank them.

Two things distinguish this from a pure reconstruction sweep.

Distortion is measured per sample and restricted to observed entries. Indexing
a 2D array with a 2D boolean mask flattens it, so computing one cosine over all
observed entries at once gives a norm-weighted global statistic rather than the
mean per-sample cosine. The two disagree whenever sample norms are
heterogeneous, and can rank models in opposite order.

Rate is measured alongside it. A model that has collapsed to a handful of active
latent dimensions can reconstruct well while carrying a latent code that is
useless for downstream distances, so reconstruction alone applies a selection
pressure toward the least regularised model in the grid. n_active and kl_total
are reported for every model and must be read together with the distortion
columns.

Note that if the trainer used the full matrix for both training and validation,
every number here is in-sample. Reconstruction is close to monotonic in capacity
in-sample, so a latent_dim sweep ranked on distortion alone is largely
predetermined. Hold out a row subset in the trainer if the ranking is meant to
say something about generalisation.
"""
import argparse
import gc
import json
import os
import warnings

import numpy as np
import pandas as pd
import tensorflow as tf
from tqdm import tqdm

from vae import VAE_bilayer

# Config keys that should be numeric. from_dict(...).T yields object columns,
# which sort and group lexicographically ("100" before "16").
_NUMERIC_CFG = ("original_dim", "hidden_dim1", "hidden_dim2", "latent_dim",
                "rec_weight", "beta", "free_bits", "dropout_rate")

# Columns identifying one grid cell, for aggregating over replicates.
_GROUP_COLS = ("loss_type", "hidden_dim1", "hidden_dim2", "latent_dim",
               "rec_weight", "beta", "free_bits")


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('folder', type=str, help="Path to models folder")
    p.add_argument('--data', type=str, default=None,
                   help="Path to original matrix .pkl(.gz). If omitted, uses <folder>/original.pkl.gz")
    p.add_argument('--batch-size', dest='batch_size', type=int, default=512,
                   help="Batch size for encoding and decoding (default: 512)")
    p.add_argument('--active-threshold', dest='active_threshold', type=float, default=0.01,
                   help="Nats above which a latent dimension counts as active (default: 0.01)")
    return p.parse_args()


def is_model_folder(path):
    return (
        os.path.isdir(path)
        and os.path.isfile(os.path.join(path, "config.json"))
        and os.path.isfile(os.path.join(path, "vae.weights.h5"))
    )


def read_model_params(fld):
    with open(os.path.join(fld, "config.json"), "r") as f:
        return json.load(f)


def _reconstruction_metrics(vae, x_filled, mask, batch_size):
    """Per-sample distortion metrics, computed in batches.

    Uses mu rather than a sampled z, and training=False, so the result is
    deterministic and does not depend on dropout. This is therefore not the same
    quantity as model.evaluate(), which samples z.

    Per-sample cosine is accumulated from masked inner products and masked norms.
    Masking the reconstruction norm as well as the input norm matters: leaving it
    unmasked lets predictions at missing entries change the denominator.
    """
    n = x_filled.shape[0]
    acc = {k: 0.0 for k in ("num_m", "nx2_m", "nr2_m", "num_u", "nx2_u", "nr2_u")}
    rows = {k: [] for k in acc}
    sq_m = sq_u = ab_m = ab_u = 0.0
    n_obs = 0

    for start in range(0, n, batch_size):
        xb = x_filled[start:start + batch_size]
        mb = mask[start:start + batch_size]
        chunk = tf.convert_to_tensor(xb, dtype=tf.float64)
        mu, _ = vae.encode(chunk, training=False)
        rb = vae.decode(mu).numpy()

        mf = mb.astype(np.float64)
        rows["num_m"].append(np.sum(xb * rb * mf, axis=1))
        rows["nx2_m"].append(np.sum(xb * xb * mf, axis=1))
        rows["nr2_m"].append(np.sum(rb * rb * mf, axis=1))
        rows["num_u"].append(np.sum(xb * rb, axis=1))
        rows["nx2_u"].append(np.sum(xb * xb, axis=1))
        rows["nr2_u"].append(np.sum(rb * rb, axis=1))

        err = xb - rb
        sq_u += float(np.sum(err ** 2))
        ab_u += float(np.sum(np.abs(err)))
        sq_m += float(np.sum((err ** 2) * mf))
        ab_m += float(np.sum(np.abs(err) * mf))
        n_obs += int(mb.sum())

        del chunk, mu, rb

    cat = {k: np.concatenate(v) for k, v in rows.items()}

    def mean_cos(num, nx2, nr2):
        nx, nr = np.sqrt(nx2), np.sqrt(nr2)
        ok = (nx > 0) & (nr > 0)
        if not ok.any():
            return float('nan')
        return float(np.mean(num[ok] / (nx[ok] * nr[ok])))

    n_all = float(x_filled.size)
    return {
        # Primary distortion metric: mean per-sample cosine over observed entries.
        'cos_observed':   mean_cos(cat["num_m"], cat["nx2_m"], cat["nr2_m"]),
        'cos_unmasked':   mean_cos(cat["num_u"], cat["nx2_u"], cat["nr2_u"]),
        'mse_observed':   sq_m / n_obs if n_obs else float('nan'),
        'rmse_observed':  np.sqrt(sq_m / n_obs) if n_obs else float('nan'),
        'mae_observed':   ab_m / n_obs if n_obs else float('nan'),
        # Diagnostics only. These include entries imputed as 0 upstream, so they
        # partly measure agreement with the imputation rather than with data.
        'mse_unmasked':   sq_u / n_all,
        'rmse_unmasked':  np.sqrt(sq_u / n_all),
        'mae_unmasked':   ab_u / n_all,
        'frac_observed':  n_obs / n_all,
    }


def evaluate_model(model_folder, x_filled, mask, batch_size, active_threshold):
    """Return distortion and rate metrics for one saved model, or None if unloadable."""
    try:
        vae = VAE_bilayer.load_vae(model_folder)
    except (ValueError, KeyError, TypeError) as exc:
        warnings.warn(f"Skipping {os.path.basename(model_folder)}: {exc}", stacklevel=2)
        return None

    try:
        scores = _reconstruction_metrics(vae, x_filled, mask, batch_size)

        # --- Rate. Read together with the distortion columns above. ---
        kl = vae.kl_per_dimension(x_filled, batch_size=batch_size,
                                  threshold=active_threshold)
        scores.update({
            'kl_total':      kl['kl_total'],
            'n_active':      kl['n_active'],
            'n_informative': kl['n_informative'],
            # Fraction of the nominal latent_dim actually carrying information.
            'frac_active':   kl['n_active'] / kl['latent_dim'],
            'kl_median_dim': kl['quantiles'][0.50],
            'kl_p95_dim':    kl['quantiles'][0.95],
        })

        return scores
    finally:
        del vae
        tf.keras.backend.clear_session()
        gc.collect()


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    for col in _NUMERIC_CFG:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def aggregate_replicates(models: pd.DataFrame) -> pd.DataFrame:
    """Mean, sd and n per grid cell.

    Ranking individual folders selects the best replicate, which is upward
    biased and discards the run-to-run variance the replicates were run to
    measure. Selection should use the cell mean; the sd says whether the
    difference between two cells exceeds the noise.
    """
    group = [c for c in _GROUP_COLS if c in models.columns]
    if not group:
        return pd.DataFrame()
    metric_cols = ['cos_observed', 'mse_observed', 'mae_observed',
                   'kl_total', 'n_active', 'frac_active']
    metric_cols = [c for c in metric_cols if c in models.columns]
    agg = models.groupby(group, dropna=False)[metric_cols].agg(['mean', 'std', 'count'])
    agg.columns = ['_'.join(c) for c in agg.columns]
    rename = {f'{metric_cols[0]}_count': 'n_replicates'}
    return agg.rename(columns=rename).sort_values(f'{metric_cols[0]}_mean',
                                                  ascending=False)


def main():
    args = parse_cli()

    x_path = args.data or os.path.join(args.folder, 'original.pkl.gz')
    x_df = pd.read_pickle(x_path)
    mask = x_df.notna().values
    x_filled = x_df.fillna(0).values.astype(np.float64)
    print(f"Matrix: {x_filled.shape}, observed fraction {mask.mean():.3f}")

    folders = [e for e in os.scandir(args.folder) if is_model_folder(e.path)]
    if not folders:
        raise FileNotFoundError(
            f"No model folders with config.json + vae.weights.h5 found in: {args.folder}"
        )

    rows, skipped = [], 0
    for m in tqdm(sorted(folders, key=lambda e: e.name)):
        cfg = read_model_params(m.path)
        if cfg.get("original_dim") != x_filled.shape[1]:
            warnings.warn(
                f"Skipping {m.name}: config original_dim={cfg.get('original_dim')} "
                f"does not match the matrix ({x_filled.shape[1]} features).",
                stacklevel=2,
            )
            skipped += 1
            continue

        scores = evaluate_model(m.path, x_filled, mask,
                               args.batch_size, args.active_threshold)
        if scores is None:
            skipped += 1
            continue

        row = pd.DataFrame([{**cfg, **scores,
                             'model': m.name, 'model_path': m.path}])
        rows.append(row)

    if not rows:
        raise RuntimeError("No models could be evaluated.")

    models = coerce_numeric(pd.concat(rows, ignore_index=True)).set_index('model')
    models.sort_values('cos_observed', ascending=False, inplace=True)

    out = os.path.join(args.folder, 'models_summary.csv')
    models.to_csv(out)

    agg = aggregate_replicates(models.reset_index())
    if not agg.empty:
        agg_out = os.path.join(args.folder, 'models_summary_by_cell.csv')
        agg.to_csv(agg_out)

    # --- Report ---
    if skipped:
        print(f"\nSkipped {skipped} folder(s); see warnings above.")

    show = [c for c in ('loss_type', 'hidden_dim1', 'hidden_dim2', 'latent_dim', 
                        'n_active', 'cos_observed', 'mse_observed',
                        'kl_total', 'frac_active') if c in models.columns]
    print("\nBest 10 by mean per-sample cosine over observed entries:")
    print(models[show].head(10).to_string())

    # if not agg.empty:
    #     print("\nBy grid cell (mean over replicates):")
    #     print(agg.head(10).to_string())

    collapsed = models[models['frac_active'] < 0.5] if 'frac_active' in models else models.iloc[:0]
    if len(collapsed):
        print(f"\nWARNING: {len(collapsed)} model(s) use under half their nominal "
              "latent dimensions. High cosine with a low n_active means the code "
              "is narrower than latent_dim suggests.")

    print(f"\nWritten: {out}")
    if not agg.empty:
        print(f"Written: {agg_out}")


if __name__ == '__main__':
    main()
