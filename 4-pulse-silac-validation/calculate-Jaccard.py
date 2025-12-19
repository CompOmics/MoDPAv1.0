#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from Jaccard_idx import *


# In[2]:
data_binary = pd.read_csv(
    "binary-matrix.csv.gz",
    # nrows=100, #for testing
    index_col=0
)


# In[16]:
results_df = calculate_jaccard_matrix(data_binary, pvalue_threshold=0.05)


# In[17]:
significant = results_df[results_df['significant']].copy()
significant.drop(columns=['pvalue'], inplace=True)
significant.jaccard = significant.jaccard.round(2)
significant.rename(columns={'obs1':'nodeA', 'obs2':'nodeB'}, inplace=True)
significant.nodeA = significant.nodeA.str.split('[][]').apply(lambda x: x[0]+x[1])
significant.nodeB = significant.nodeB.str.split('[][]').apply(lambda x: x[0]+x[1])
print(f"Significant pairs (p < 0.05): {len(significant)}/{len(results_df)}")
significant.head()


# In[18]:
significant.to_csv('20251021-1101-relaxed_carver/jaccard-similarities.csv', index=False)
