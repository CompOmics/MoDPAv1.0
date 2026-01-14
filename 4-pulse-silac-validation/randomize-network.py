#!/usr/bin/env python
# coding: utf-8
import networkx as nx
import numpy as np
import pandas as pd
import argparse

def Degree_preserving_randomization(G):
    print(len(G.edges()))
    for _ in range(len(G.edges())):
        edgelist = list(G.edges())
        idx1, idx2 = np.random.choice(range(len(edgelist)), 2, replace=False)
        edge1, edge2 = edgelist[idx1], edgelist[idx2]
        weight1, weight2 = G.edges[edge1]['weight'], G.edges[edge2]['weight']
        
        if len(set([edge1[0], edge2[0], edge1[1], edge2[1]]))==4:
            # all 4 nodes must be different!
            newedge1, newedge2 = (edge1[0], edge2[0]), (edge1[1], edge2[1])
            if not (newedge1 in G.edges() or newedge2 in G.edges()):
                # skip if edges already there
                G.remove_edge(*edge1)
                G.remove_edge(*edge2)
                G.add_edge(*newedge1, weight=weight1)
                G.add_edge(*newedge2, weight=weight2)
    print(len(G.edges()))

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("input", type=str, help="Path to input file (nonrandom network)")
    p.add_argument("--output", dest='output', type=str, default='deegree-preserved-random-network.csv.gz', 
                   help="Path to output file (default: deegree-preserved-random-network.csv.gz)")
    return p.parse_args()



args = parse_cli()
edgelist = pd.read_csv(args.input, usecols=['nodeA','nodeB','Score'])
# edgelist = edgelist[edgelist.distance > 0.38].copy(deep=True)
G = nx.Graph()
G.add_weighted_edges_from(
    [_ for _ in zip(edgelist.nodeA, edgelist.nodeB, edgelist.Score)]
)
len(G.edges())

Degree_preserving_randomization(G)

edgelist2 = []
for edge in list(G.edges()):
    edgelist2.append([
        edge[0],
        edge[1],
        G.edges()[edge]['weight'],
        abs(G.edges()[edge]['weight'])
    ])
edgelist2 = pd.DataFrame(edgelist2, columns=['nodeA','nodeB','Score','distance'])
edgelist2.to_csv(args.output, index=False, encoding='utf8', compression='gzip')

