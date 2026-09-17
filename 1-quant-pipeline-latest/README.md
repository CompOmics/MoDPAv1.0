# Step 1. Preprocessing and PTM quantification

This folder turns peptidoform identifications and PSM counts into the PTM-by-MS-run matrix that
the VAE is trained on. Four scripts run in a fixed order. Each one writes a file that the next one
reads, so the order cannot be changed.

## Inputs

| Input | Description |
| --- | --- |
| Peptidoform IDs file | one row per identified peptidoform, with its modifications |
| Peptidoform counts file | PSM counts per peptidoform and per MS run |
| Target FASTA | proteins to keep. `Human_2026_01_canonical.fasta.gz` is included |
| Contaminant FASTA | proteins to remove. `MQcontaminants_2023_11_14.fasta.gz` is included |
| PTMs of interest | the modifications to extract, as a CSV (see below) |

The PTMs of interest file has three columns, one row per residue and modification pair:

| AA | unimod_id | ptm_name |
| --- | --- | --- |
| C | 4 | Carbamidomethyl |
| M | 35 | Oxidation |

Four lists are provided:

| File | Contents |
| --- | --- |
| `PTMs-of-interest-submission.csv` | the 20 residue and modification pairs used in the manuscript |
| `PTMs-of-interest-val.csv` | phosphorylation, dehydration and sulfation, for the validation runs |
| `PTMs-of-interest-phospho.csv` | phosphorylation only |
| `PTMs-of-interest.csv` | carbamidomethylation and oxidation, the default of every script |

## Steps

### 1. Map peptidoforms to proteins and count them

```bash
python Map_and_count_peptides.py <peptidoform_ids.csv.gz> <peptidoform_counts.csv.gz> <target.fasta.gz> -p <prefix>
```

Maps every peptide to the proteins that contain it, assigns a lead protein, groups the
identifications into peptidoforms, and joins the PSM counts. Peptides are kept when they are 7 to
30 residues long and their precursor m/z is between 500 and 6000.

Writes `<prefix>_Peptidoforms_counts_mapped.csv.gz`.

Use `-t` to cap the number of Polars threads. Polars scales its memory use with the thread count,
so lowering it is the first thing to try if the step runs out of memory.

### 2. Compute relative PTM counts

```bash
python Relative_counts.py <prefix>_Peptidoforms_counts_mapped.csv.gz
```

Expresses the PSM count of each modified site as a fraction of the counts observed at that position
in the same MS run. N-terminal modifications, semi-tryptic peptides and ragging products are
removed, and a peptidoform needs at least 2 PSMs to be counted.

The step partitions the data by MS run, writes one temporary file per run into a `counts-per-msrun`
directory next to the input, and concatenates them at the end. The temporary directory is removed
on success. The script refuses to start if it already exists, so delete it before a rerun.

Writes `<prefix>_PTMs_counts_relative.csv.gz`.

### 3. Extract one matrix per modification

```bash
python Generate_PTM_matrices.py <prefix>_PTMs_counts_relative.csv.gz \
    -t Human_2026_01_canonical.fasta.gz \
    -c MQcontaminants_2023_11_14.fasta.gz \
    -m PTMs-of-interest-submission.csv
```

Keeps the PTM events on target proteins, drops the ones on contaminants, and pivots the remaining
relative counts into one PTM-by-MS-run matrix per modification.

Writes `MoDPA_matrices_<prefix>/MoDPA_Rel_[<unimod_id>]<ptm_name>_<AA>.pkl.gz`, one file per row of
the PTMs of interest file. A modification with no usable data is skipped.

### 4. Combine the matrices

```bash
python Combine_PTM_matrices.py MoDPA_matrices_<prefix> -p <prefix> -m PTMs-of-interest-submission.csv -s 0.05 -t 5
```

Concatenates the per-modification matrices, drops PTM events whose relative counts vary less than
`-s` (standard deviation, default 0.05), and applies the `-t` observation threshold through the
`MoDPA` class in `MoDPA.py`: a column is kept when it holds at least `-t` non-zero values, then a
row is kept when it holds at least `-t` non-zero values among the surviving columns. Each column is
finally divided by its maximum.

Writes two files into the same folder:

| File | Content |
| --- | --- |
| `<prefix>-PTMs-thresh<t>r<t>c-std<s>.pkl.gz` | the MoDPA matrix, the input of step 2 |
| `<prefix>-PTMs-thresh<t>r<t>c-std<s>-analyzed-ptms.csv` | one row per retained PTM event |

Matrix rows are PTM events identified as `UniProtAccession|Position|Residue|[UnimodID]Name`, and
columns are MS runs.

## Running the four steps in one go

`quant_pipeline_executable.py` calls the four functions in order. Edit the paths at the top of the
file and run it:

```bash
python quant_pipeline_executable.py
```

`quant-pipeline-notebook.ipynb` is the same pipeline as a notebook.

## `MoDPA.py`

Holds the `MoDPA` class used by step 4 for filtering, naming and column normalisation. The class
also carries a column and row shuffled copy of the matrix and a pseudo-Jaccard correlation
implementation. Neither is used by the current pipeline, which measures association in the VAE
latent space instead (see `../2-VAE-code`).
