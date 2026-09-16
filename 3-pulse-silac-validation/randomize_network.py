#!/usr/bin/env python
# coding: utf-8
import argparse
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("input", type=str, help="Path to input file (nonrandom network)")
    p.add_argument("-o", "--output", dest='output', type=str, default='degree-preserved-random-network.csv.gz',
                   help="Path to output file (default: degree-preserved-random-network.csv.gz)")
    p.add_argument("-s", "--swaps", dest='swaps', type=int, default=10,
                   help="Number of swap attempts per edge (default: 10)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducible network randomization (default: None)")
    return p.parse_args()


def Degree_preserving_randomization(G, swaps_per_edge=10, seed=None):
    # Build edge list once; maintain a set for O(1) existence checks.
    edges = list(G.edges(data='weight'))
    n = len(edges)
    n_attempts = n * swaps_per_edge
    print(f"Edges: {n}, swap attempts: {n_attempts}")

    edge_set = set()
    for u, v, _ in edges:
        edge_set.add((u, v))
        edge_set.add((v, u))

    rng = np.random.default_rng(seed)
    i_all = rng.integers(0, n, n_attempts, dtype=np.int32)
    j_all = rng.integers(0, n, n_attempts, dtype=np.int32)

    swaps_done = 0
    for k in tqdm(range(n_attempts)):
        i, j = int(i_all[k]), int(j_all[k])
        if i == j:
            continue
        u1, v1, w1 = edges[i]
        u2, v2, w2 = edges[j]

        if len({u1, v1, u2, v2}) == 4:
            if (u1, u2) not in edge_set and (v1, v2) not in edge_set:
                edge_set.discard((u1, v1)); edge_set.discard((v1, u1))
                edge_set.discard((u2, v2)); edge_set.discard((v2, u2))
                edges[i] = (u1, u2, w1)
                edges[j] = (v1, v2, w2)
                edge_set.add((u1, u2)); edge_set.add((u2, u1))
                edge_set.add((v1, v2)); edge_set.add((v2, v1))
                swaps_done += 1

    print(f"Swaps completed: {swaps_done} / {n_attempts} attempts")

    G.clear()
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)


def main():
    args = parse_cli()
    edgelist = pd.read_csv(args.input, usecols=['nodeA','nodeB','Score'])
    G = nx.Graph()
    G.add_weighted_edges_from(zip(edgelist.nodeA, edgelist.nodeB, edgelist.Score))

    Degree_preserving_randomization(G, swaps_per_edge=args.swaps, seed=args.seed)

    edgelist2 = []
    for edge in G.edges():
        edgelist2.append([
            edge[0],
            edge[1],
            G.edges[edge]['weight'],
            abs(G.edges[edge]['weight'])
        ])
    edgelist2 = pd.DataFrame(edgelist2, columns=['nodeA','nodeB','Score','distance'])
    edgelist2.to_csv(args.output, index=False, encoding='utf8')


if __name__ == "__main__":
    main()
