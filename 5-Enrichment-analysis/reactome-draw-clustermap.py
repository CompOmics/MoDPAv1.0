#!/usr/bin/env python
# coding: utf-8
from reactome2py import content, analysis
import os, time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data_path = "./without_sulfo/compassionate_buck-v2/Leiden-0.5/nodes.csv"
FLD = os.path.split(data_path)[0]

combo = pd.read_csv(os.path.join(FLD,'results_combo.csv'))

results_heatmap = combo.pivot(
    index='name',
    columns='cluster',
    values='FDR'
)
results_heatmap = -1 * np.log10(results_heatmap)
results_heatmap.to_csv(os.path.join(FLD,'results_combo_heatmap.csv'))

reac_clusters = pd.read_csv('ReactomePathsClusters.csv', usecols=['name'])
reac_clusters = reac_clusters.merge(results_heatmap.reset_index(), on='name')
reac_clusters.set_index('name', inplace=True)

h = sns.clustermap(reac_clusters.fillna(0), method='ward', 
               vmin=0, vmax=8,
               cmap='vlag', center=0, linewidths=0.5, figsize=(12,18),
              cbar_pos=(0.9, 0.9, 0.05, 0.08),
                   # cbar_pos=None,
                   row_cluster=False,
              dendrogram_ratio=.05,)
# plt.xticks(rotation=90)
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xticklabels(), rotation=90)
h.savefig(os.path.join(FLD,'Reactome_clustermap_with_legend.png'), dpi=300, bbox_inches='tight')
