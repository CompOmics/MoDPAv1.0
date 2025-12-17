#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np

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
        return True
    elif c==-1 and edge_row.Score < 0:
        return True
    # elif c==0:
    #     return np.nan
    else:
        return False


# In[3]:
v2 = pd.read_csv("./20251021-1101-relaxed_carver/20251027-1609-relaxed_carver-signed-distances.csv.gz")

all_nodes = list(set(list(v2.nodeA) + list(v2.nodeB)))
len(all_nodes)

len(all_nodes) * (len(all_nodes)-1) /2


all_possible_edges = []
for i in range(len(all_nodes)):
    for j in range(i+1,len(all_nodes)):
        all_possible_edges.append(f"{all_nodes[i]}__{all_nodes[j]}")
len(all_possible_edges)

random_edges = np.random.choice(all_possible_edges, len(v2), replace=False)
len(random_edges)

fully_random = pd.DataFrame(zip(random_edges, v2.Score), columns=['edge','Score'])
fully_random[['nodeA','nodeB']] = fully_random.edge.str.split("__", expand=True)
fully_random['distance'] = fully_random.Score.apply(abs)
fully_random


# In[10]:
fully_random[['nodeA','nodeB','Score','distance']].to_csv(
    "./20251021-1101-relaxed_carver/20251027-1609-relaxed_carver-signed-distances-full-random.csv.gz",
    index=False,
    compression='gzip',
    encoding='utf8'
)
