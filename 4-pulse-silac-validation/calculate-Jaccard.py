#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from itertools import combinations
from Jaccard_idx import *


# In[2]:
data_binary = pd.read_csv(
    "C:/Users/Enrico/OneDrive - UGent/FWO_project/MoDPA/pulse-silac-test/binary-matrix.csv.gz",
    # nrows=1000,
    index_col=0
)
# data_binary


# In[16]:
results_df = calculate_jaccard_matrix(data_binary, pvalue_threshold=0.05)


# In[17]:
significant = results_df[results_df['significant']]
significant = significant.rename(columns={'obs1':'nodeA', 'obs2':'nodeB'})
significant.nodeA = significant.nodeA.str.split('[][]').apply(lambda x: x[0]+x[1])
significant.nodeB = significant.nodeB.str.split('[][]').apply(lambda x: x[0]+x[1])
print(f"Significant pairs (p < 0.05): {len(significant)}/{len(results_df)}")
significant.head()


# In[18]:
significant.to_csv('20251021-1101-relaxed_carver/jaccard-similarities.csv.gz', index=False, compression='gzip')
