# Step 3. Pulsed SILAC validation

Pulsed SILAC gives a dataset in which the expected sign of an association is known in advance, so
the MoDPA network can be scored against ground truth rather than against expectation alone.

In pulsed SILAC, labelled amino acids are added to the growth medium for a limited period. Over
time the light forms decrease as the protein pool turns over, while the heavy forms increase.
Unlabelled arginine (R0) and lysine (K0) therefore follow the same trend as each other, the heavy
counterparts R10 and K8 follow the same trend as each other, and the two groups follow opposite
trends. Modifications such as methionine oxidation are not expected to track either group.

## Edge labelling rule

Each node is labelled from the Unimod identifier at the end of its node name:

| Node suffix | Label |
| --- | --- |
| `\|0` | light |
| `\|259` (Label:13C(6)15N(2), K8) | heavy |
| `\|267` (Label:13C(6)15N(4), R10) | heavy |
| anything else | unlabelled |

An edge is then counted as valid when

1. it is positive and joins two heavy or two light nodes, or
2. it is negative and joins one heavy node and one light node.

Every other edge between two labelled nodes is invalid. Edges touching an unlabelled node are
excluded from the denominator rather than counted as invalid.

## Randomised networks

The MoDPA network is compared against two null networks built from the same edge list.

```bash
python randomize_network.py <signed-distances.csv.gz> -o degree-preserved-random-network.csv.gz --swaps 10 --seed 42
```

Double-edge swaps that leave the degree of every node unchanged. `--swaps` sets the number of swap
attempts per edge. A swap is attempted only between two edges with four distinct endpoints and is
rejected when it would create an edge that already exists, so the number of swaps completed is
lower than the number attempted. Pass `--seed` for a reproducible result.

```bash
python generate_full_random_network.py <signed-distances.csv.gz> -o random-network.csv.gz --seed 42
```

Draws the same number of edges uniformly from all node pairs, with no constraint. The observed
scores are reassigned to the drawn edges. This enumerates every possible pair, so its memory use
grows with the square of the node count.

## Validation

```bash
python validate_edges.py \
    -r <signed-distances.csv.gz> \
    -n random-network.csv.gz \
    -d degree-preserved-random-network.csv.gz \
    -o validated-edges-plot.png
```

Labels the edges of the three networks, applies the rule above, and plots the proportion of valid
edges against the minimum absolute correlation, swept from 0.4 to 0.95 in steps of 0.05.

## Contents

| File | Description |
| --- | --- |
| `20260814-0939-vigilant_pike/` | the association list of the pulsed SILAC model run and the two randomised networks derived from it |
| `vigilant-pike-plot-reworked-model.png` | validation plot of the current model |
| `relaxed-carver-plot-1st-submission.png` | validation plot of the model of the first submission |

The three files in `20260814-0939-vigilant_pike/` are compressed data files and are therefore not
tracked in git. They are part of the Zenodo release. See the data availability section of the root
README.
