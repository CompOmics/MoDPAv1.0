#!/usr/bin/env python
# coding: utf-8

# code to sort reactome pathways nicely into clusters
import os
import pandas as pd
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

data_path = "./without_sulfo/compassionate_buck-v2/Leiden-0.5/nodes.csv"
FLD = os.path.split(data_path)[0]

# In[2]: Links to download.reactome might need to be updated
pathways = pd.read_csv('https://download.reactome.org/95/ReactomePathways.txt', sep='\t', header=None)
pathways.columns = ['ID','name','species']
human_pathways = pathways[pathways.species=='Homo sapiens'].copy(deep=True)
id2name = human_pathways.set_index('ID').to_dict()['name']


# In[3]:
hierarchy = pd.read_csv('https://download.reactome.org/95/ReactomePathwaysRelation.txt', sep='\t', header=None)
hierarchy.columns = ['parent','child']
hierarchy_hs = hierarchy[hierarchy.child.isin(id2name)].copy(deep=True)
hierarchy_hs.to_csv('ReactomePathwaysHierarchyHuman.csv', index=False, encoding='utf-8')


# In[5]:
hierarchy_hs['parent_name'] = hierarchy_hs['parent'].map(id2name)
hierarchy_hs['child_name']  = hierarchy_hs['child'].map(id2name)
hierarchy_hs.sort_values('parent_name')


# In[6]:
import networkx as nx
G = nx.Graph()
G.add_edges_from(
    [_ for _ in zip(hierarchy_hs.parent, hierarchy_hs.child)]
)
print(len(G.edges()))


# In[7]:
def pathway_distance(G,a,b):
    try:
        return len(nx.shortest_path(G, a, b))
    except nx.NetworkXNoPath as E:
        return 100
    except:
        return 10000

# hierarchy_hs['distance'] = hierarchy_hs.apply(lambda row: pathway_distance(G, row.parent, row.child), axis=1)
# hierarchy_hs.sort_values('distance')


# In[8]:
all_pathways = list(set(pd.read_csv(os.path.join(FLD,'results_combo.csv')).stId))
all_distances = []
for a in all_pathways:
    for b in all_pathways:
        all_distances.append([a, b, pathway_distance(G,a,b)])
all_distances = pd.DataFrame(all_distances, columns=['pathwayA', 'pathwayB','distance'])
all_distances['pathwayA_name'] = '[' + all_distances.pathwayA + '] ' + all_distances.pathwayA.map(id2name)
all_distances['pathwayB_name'] = '[' + all_distances.pathwayB + '] ' + all_distances.pathwayB.map(id2name)
all_distances.tail()


# In[9]:
all_distances = all_distances.pivot(index='pathwayA_name', columns='pathwayB_name', values='distance')
all_distances.head()


# In[11]:
Z = linkage(all_distances, method='ward')
den = dendrogram(
    Z, 
    # truncate_mode='lastp'
    color_threshold=100
)
clusters = pd.DataFrame(([list(all_distances.index)[i] for i in den['leaves']], den['leaves_color_list'])).T
clusters.columns = ['name','cluster']
clusters.to_csv('ReactomePathsClusters.csv', index=False)
