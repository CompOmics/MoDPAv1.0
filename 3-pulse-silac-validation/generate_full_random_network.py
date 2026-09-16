#!/usr/bin/env python
# coding: utf-8
import argparse
import numpy as np
import pandas as pd


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("input", type=str, help="Path to input file (nonrandom network)")
    p.add_argument("-o", "--output", dest='output', type=str, default='random-network.csv.gz',
                   help="Path to output file (default: random-network.csv.gz)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducible network randomization (default: None)")
    return p.parse_args()


def main():
    args = parse_cli()
    edges = pd.read_csv(args.input)
    all_nodes = pd.unique(pd.concat([edges.nodeA, edges.nodeB], ignore_index=True)).tolist()

    all_possible_edges = []
    for i in range(len(all_nodes)):
        for j in range(i+1,len(all_nodes)):
            all_possible_edges.append(f"{all_nodes[i]}__{all_nodes[j]}")

    rng = np.random.default_rng(args.seed)
    random_edges = rng.choice(all_possible_edges, len(edges), replace=False)

    fully_random = pd.DataFrame(zip(random_edges, edges.Score), columns=['edge','Score'])
    fully_random[['nodeA','nodeB']] = fully_random.edge.str.split("__", expand=True)
    fully_random['distance'] = fully_random.Score.abs()

    fully_random[['nodeA','nodeB','Score','distance']].to_csv(
        args.output,
        index=False,
        compression='gzip',
        encoding='utf8'
    )


if __name__ == "__main__":
    main()
