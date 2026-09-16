#!/usr/bin/env python
# coding: utf-8
import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def existing_file(path: str) -> str:
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"file not found: {path}")
    return path


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("-r", "--path_real", type=existing_file, required=True, help="Path to file with nonrandom network edges.")
    p.add_argument("-n", "--path_random", type=existing_file, required=True, help="Path to file with random network edges.")
    p.add_argument("-d", "--path_degree_preserving", type=existing_file, required=True, help="Path to file with degree-preserved random network edges.")
    p.add_argument("-o", "--path_output", help="Path to output file.", default="validated-edges-plot.png")
    return p.parse_args()


def get_label(node):
    if node.endswith("|0"):
        return "L"
    elif node.endswith("|259"):
        return "H"
    elif node.endswith("|267"):
        return "H"
    else:
        return "u" # unlabelled


def label_to_num(x):
    if x=='H':
        return 1
    elif x=='L':
        return -1
    else:
        return 0


def validate_edge(edge_row):
    a = label_to_num(edge_row.labelA)
    b = label_to_num(edge_row.labelB)
    c = a * b
    # c =  1 means they have the same label
    # c = -1 means they have different labels
    # c =  0 means either one is unlabelled
    if c==0:
        return 0
    elif c==1 and edge_row.Score > 0:
        return 1
    elif c==-1 and edge_row.Score < 0:
        return 1
    else:
        return 0


def count_validated(df, min_abs_corr, network):
    # A negative association between a H site and a L site is VALID
    # A positive association between 2 H sites OR 2 L sites is VALID
    # Unlabelled edges are excluded from the validation denominator
    df = df[(df.abs_corr >= min_abs_corr)&(df.network_type==network)]
    tmp = df.validated.dropna()
    if tmp.empty:
        return np.nan
    return tmp.mean()


def main():
    args = parse_cli()

    real = pd.read_csv(args.path_real)
    real['network_type'] = 'Signed distance correlation'

    random = pd.read_csv(args.path_degree_preserving)
    random['network_type'] = 'Degree-preserved random'

    fully_random = pd.read_csv(args.path_random)
    fully_random['network_type'] = 'Fully random'

    data = pd.concat([real, random, fully_random], ignore_index=True)
    data.rename(columns={'distance':'abs_corr'}, inplace=True)

    data['labelA'] = data.nodeA.apply(get_label)
    data['labelB'] = data.nodeB.apply(get_label)
    data['validated'] = data.apply(validate_edge, axis=1)
    label_pair = data[['labelA','labelB']].apply(lambda row: tuple(sorted(row)), axis=1)
    print(label_pair.value_counts(sort=False))

    validated_edges = []
    for i in np.arange(0.4, 1, .05):
        for j in ['Signed distance correlation', 'Degree-preserved random', 'Fully random']:
            validated_edges.append([i, j, count_validated(data, i, j)])

    validated_edges = pd.DataFrame(validated_edges, columns=['min_abs_corr','network','percent_validated'])
    sns.lineplot(data=validated_edges, x='min_abs_corr', y='percent_validated', hue='network')
    plt.xlabel('Min Abs Correlation')
    plt.ylabel('Proportion Validated Edges')
    plt.ylim(0,1)
    plt.savefig(args.path_output, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()
