#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import argparse


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("input", type=str, help="Path to input file (nonrandom network)")
    p.add_argument("--output", dest='output', type=str, default='random-network.csv.gz', 
                   help="Path to output file (default: random-network.csv.gz)")
    return p.parse_args()


args = parse_cli()
edges = pd.read_csv(args.input)
all_nodes = list(set(list(edges.nodeA) + list(edges.nodeB)))

all_possible_edges = []
for i in range(len(all_nodes)):
    for j in range(i+1,len(all_nodes)):
        all_possible_edges.append(f"{all_nodes[i]}__{all_nodes[j]}")
len(all_possible_edges)

random_edges = np.random.choice(all_possible_edges, len(edges), replace=False)

fully_random = pd.DataFrame(zip(random_edges, edges.Score), columns=['edge','Score'])
fully_random[['nodeA','nodeB']] = fully_random.edge.str.split("__", expand=True)
fully_random['distance'] = fully_random.Score.apply(abs)
fully_random

fully_random[['nodeA','nodeB','Score','distance']].to_csv(
    args.output,
    index=False,
    compression='gzip',
    encoding='utf8'
)
