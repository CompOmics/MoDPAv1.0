#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse, os


# In[2]:
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("path_real", type=str, help="Path to file with nonrandom network edges.")
    p.add_argument("path_random", type=str, help="Path to file with random network edges.")
    p.add_argument("path_degree_preserving", type=str, help="Path to file with degree-preserved random network edges.")
    return p.parse_args()


def get_label(node):
    if node.endswith("|0"):
        return "L"
    elif node.endswith("|259"):
        return "H"
    elif node.endswith("|267"):
        return "H"
    else:
        return np.nan

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
    # c=1 means they have the same label; c=-1 means they have different labels
    # c=0 means either one is unlabelled
    if c==1 and edge_row.Score > 0:
        return 1
    elif c==-1 and edge_row.Score < 0:
        return 1
    # elif c==0:
    #     return np.nan
    else:
        return 0

def count_validated(df, min_abs_corr, network):
    # A negative association between a H site and a L site is VALID
    # A positive association between 2 H sites OR 2 L sites is VALID
    # Anything else is NOT VALID
    df = df[(df.abs_corr >= min_abs_corr)&(df.network_type==network)]
    tmp = df.validated
    return np.sum(tmp)/len(tmp)


# In[3]:
args = parse_cli()

real   = pd.read_csv(args.path_real)
real['network_type'] = 'Signed distance correlation'

random = pd.read_csv(args.path_degree_preserving)
random['network_type'] = 'Degree-preserved random'

fully_random = pd.read_csv(args.path_random)
fully_random['network_type'] = 'Fully random'

# pearson = pd.read_csv('20251021-1101-relaxed_carver/20251021-1101-relaxed_carver-correlations-w-pval-v1.csv')
# pearson['network_type'] = 'Pearson correlation'
# pearson['distance'] = abs(pearson.Score)

# jaccard = pd.read_csv('20251021-1101-relaxed_carver/jaccard-similarities.csv.gz',
#                      usecols=['nodeA','nodeB','jaccard'])
# jaccard['Score'] = jaccard.jaccard * 2 - 1
# jaccard['distance'] = abs(jaccard.Score) 
# jaccard['network_type'] = 'Jaccard'


# In[8]:
data = pd.concat([
    real, 
    random, 
    # jaccard, 
    fully_random,
    # pearson
    ], ignore_index=True)
data.rename(columns={'distance':'abs_corr'}, inplace=True)

data['labelA'] = data.nodeA.apply(get_label)
data['labelB'] = data.nodeB.apply(get_label)
data['same_label'] = data['labelB']==data['labelA']
data['validated'] = data.apply(validate_edge, axis=1)
print(data[['labelA','labelB']].value_counts(sort=False))


# In[11]:
validated_edges = []
for i in np.arange(0.4, 1, .05):
    for j in [
        'Signed distance correlation',
        # 'Pearson correlation',
          'Degree-preserved random',
          'Fully random',
        # 'Jaccard'
             ]:
        validated_edges.append([i, j, count_validated(data, i, j)])

validated_edges = pd.DataFrame(validated_edges, columns=['min_abs_corr','network','percent_validated'])
sns.lineplot(data=validated_edges, x='min_abs_corr', y='percent_validated', hue='network')
plt.xlabel('Min Abs Correlation')
plt.ylabel('% Validated Edges')
plt.ylim(0,1)
plt.axvline(0.65, color='k', ls='--', lw=.7)
plt.savefig('20251021-1101-relaxed_carver/validated-edges-v3.png', dpi=300, bbox_inches='tight')