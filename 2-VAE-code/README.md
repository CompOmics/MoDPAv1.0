# Step 2. VAE training and PTM association scoring

This folder embeds the MoDPA matrix into a latent space with a variational autoencoder, selects a
model, and scores every PTM pair by the signed distance correlation of their latent
representations.

A GPU is recommended for training. See `../env-wsl-gpu.yml` and `../setup-wsl-gpu.sh`.

## The model

`vae.py` defines `VAE_bilayer`, a two-layer VAE. The relevant hyperparameters are the two hidden
layer widths, the latent dimensionality, the loss type
(`mean_squared_error`, `cosine_similarity`, `MSE+KL` or `cos+KL`), `rec_weight`, `free_bits` and
`dropout_rate`.

The objective is `rec_weight * L_rec + KL`. There is no separate KL weight, because only the ratio
of the two terms moves the optimum, so `rec_weight` carries the whole reconstruction and
regularisation balance on its own. One consequence is that the total loss moves by orders of
magnitude across a `rec_weight` sweep, so any early stopping criterion with an absolute threshold
is not comparable between grid cells. `free_bits` sets the nats per latent dimension below which a
dimension is not penalised, which prevents posterior collapse once `rec_weight` is small enough for
collapse to be possible.

A saved model is a folder holding `config.json` and `vae.weights.h5`, which is what
`VAE_bilayer.load_vae()` reads.

## Steps

### 1. Grid search

```bash
python VAE_gridsearch_subprocess.py <MoDPA matrix .pkl.gz> -f <output folder>
```

Trains every hyperparameter combination. The grid is either the one written in `build_param_grid`
or a JSON file passed with `-p`. Either way the grid that was run is written to
`<output folder>/grid-search-params.json`.

Each combination is trained in its own subprocess. TensorFlow and CUDA do not return memory to the
operating system within a long-running process, so training the whole grid in one process
accumulates memory until the machine fails. A subprocess returns all of its memory when it exits.
The same file is both the driver and the worker; the driver never imports TensorFlow.

Progress is checkpointed to `<output folder>/.checkpoint.json` by parameter tuple, so rerunning the
same command resumes rather than restarting. `--restart-from <n>` reruns from a given iteration
onwards.

Useful options: `-e` early stopping patience, `-b` batch size, `-l` learning rate, `--free-bits`,
`--gpu`, and `--nondeterministic-ops` to trade bitwise reproducibility for speed.

Each trained model is written to `<output folder>/<timestamp>-s<seed>-<name>/`. The seed is derived
from the parameter tuple, so the folder name identifies the run rather than decorating it, and a
name collision is suffixed instead of merged. Alongside `config.json` and `vae.weights.h5` the
folder holds `provenance.json`, which records the seed, the parameter tuple, the git state, the
GPU, the determinism setting and a fingerprint of the weights. The training matrix is copied once
to `<output folder>/original.pkl.gz`.

### 2. Select a model

```bash
python VAE_validation.py <output folder> [--data <matrix .pkl.gz>]
```

Evaluates every model in the folder and ranks them. Distortion is measured per sample and only over
observed entries, because a single cosine over all observed entries at once is norm weighted and
can rank models in the opposite order. `n_active` and `kl_total` are reported next to it: a model
that has collapsed onto a few latent dimensions can reconstruct well while carrying a latent code
that is useless for distances, so reconstruction alone must not decide the ranking.

Writes `models_summary.csv` and, when replicates are present, `models_summary_by_cell.csv`.

If the trainer used the full matrix for training and validation, every number is in-sample.

### 3. Encode

```bash
python VAE_encode.py <model folder> <matrix .pkl.gz>
```

Writes `Latent-space.pkl.gz` and `Reconstruction.pkl.gz` into the model folder. The bracketed
modification name is stripped from the index here, so a row identified as
`P04637|15|S|[21]Phospho` in the matrix becomes the node `P04637|15|S|21`. That is the node format
used by every downstream step.

### 4. Score PTM pairs

```bash
python calculate_sdcorr.py <model folder>
```

Reads `Latent-space.pkl.gz` from the model folder and computes, for every pair of PTM events, the
distance correlation, its permutation-free t-test p-value, the Pearson correlation, and the signed
distance correlation (SDCor), which is the distance correlation carrying the sign of the Pearson
correlation. Pairs are processed in batches by a process pool, with partial results written to disk
and removed once combined. Benjamini-Hochberg correction is applied over all pairs.

Writes `<datetime>-<run name>-signed-distances.csv.gz` into the model folder, with the columns
`nodeA`, `nodeB`, `Score`, `pval`, `distance`, `PCC` and `qvalue`. `Score` is the SDCor and
`distance` is its magnitude.

This file is the MoDPA association list, and it is the input of steps 3, 4 and 5.

### 5. Build and cluster the network

```bash
python build_modpa_network.py <model folder>/<datetime>-<run name>-signed-distances.csv.gz
```

Filters the association list, builds the network, partitions it with Leiden, annotates the
partition with protein-level information from a FASTA file, and writes two CSV files into the
directory of the association list:

| File | Contents |
| --- | --- |
| `<prefix>-filtered-distances.csv` | one row per retained edge, with the annotation columns added by `parse_network` |
| `<prefix>-Leiden-clusters.csv` | one row per node, with its cluster and the protein annotations |

The prefix defaults to `<today>-True-data`, matching the published files. An edge is kept when
`Score >= --min-score` (0.6 by default) and `qvalue < --qvalue` (0.05 by default); edges flagged as
a potential artefact are annotated, not removed. Leiden runs on the RBConfiguration objective
with `--resolution` 2, `--n-iterations` 2 and `--seed` 42.

`--unweighted` clusters the network without edge weights, so every retained edge counts equally and
the partition depends on the topology alone rather than on the association scores. The default is
the weighted partition, using the absolute SDCor as the edge weight, which is the one used
throughout the published analysis. The unweighted partition is written to
`<prefix>-Leiden-unweighted-clusters.csv`, so the two never overwrite each other. Existing output is
not overwritten unless `--overwrite` is given.

Protein annotations are read from `--fasta`, which defaults to the canonical human FASTA of step 1.
The script fails if an accession in the network has no annotation, so the FASTA file must be the
one the matrix was built against.

Nodes carry the Unimod accession of their modification, not its name. `--ptm-names` adds a
`PTM_name` column to the cluster table from a PTM list of step 1, for example

```bash
python build_modpa_network.py <association list> --ptm-names ../1-quant-pipeline-latest/PTMs-of-interest-submission.csv
```

which is the list the published run was built with. The list has to be that same list, because
`Generate_PTM_matrices.py` selects PTM events on the accession and the residue together, so the
join is on both and an unnamed PTM event is an error rather than a blank. Without `--ptm-names` the
cluster table carries the accession alone, which is what every downstream step reads.

These two files are the network input of steps 3, 4 and 5.

## Diagnostics

| Script | Purpose |
| --- | --- |
| `VAE_inspect_kl_per_dimension.py` | raw KL in nats carried by each latent dimension, written to CSV. Not floored by `free_bits`, so it shows whether a `free_bits` threshold is sensible |
| `modpa_matrix_sparsity.py` | sparsity of a MoDPA matrix, with observation histograms per PTM event and per MS run |

## Model outputs

`MoDPA_models.tar.gz` (9.9 GB) holds the trained models, that is `config.json`, `vae.weights.h5`
and `provenance.json` per run. It is not tracked in git and is part of the Zenodo release.

**Weights only.** The latent spaces and the association lists were removed to bring the archive
within the Zenodo size limit, so steps 3 and 4 above have to be rerun on the unpacked models before
anything downstream can use them. Encoding is quick; scoring all PTM pairs is quadratic in the
number of PTM events and writes 1.2 GB to 1.3 GB per run, so regenerate only the runs that are
needed. See the data availability section of the root README.

## `Sensitivity-analysis/`

Quantifies how much the network changes with the latent dimensionality, against the run-to-run
variation at a fixed dimensionality. It has its own README.
