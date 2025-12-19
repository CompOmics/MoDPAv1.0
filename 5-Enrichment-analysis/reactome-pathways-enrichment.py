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
clustering_col = '__leidenCluster'
data = pd.read_csv(data_path)


for cluster,df in data.groupby(clustering_col).__iter__():
    outpath = os.path.join(FLD, f'cluster{cluster}-reac.csv')
    markers = list(set(df.UniAcc))
    
    if len(markers)<20:
        continue

    print('>> Cluster', cluster, '~', len(markers), 'proteins')
    result = analysis.identifiers(ids=','.join(markers))
    time.sleep(5)
    token = result['summary']['token']
    token_result = analysis.token(token, species='Homo sapiens', page_size='-1', page='-1', 
                                  # sort_by='ENTITIES_FDR', order='ASC', include_disease=True, 
                                  resource='TOTAL', p_value='0.05', 
                                  min_entities=5, max_entities=None)
    enriched_pathways = pd.DataFrame(token_result['pathways'])
    
    enriched_pathways['species_id'] = enriched_pathways.species.apply(lambda x: x.get('taxId',np.nan))
    enriched_pathways['FDR'] = enriched_pathways.entities.apply(lambda x: x.get('fdr',1))
    enriched_pathways['entities_tot'] = enriched_pathways.entities.apply(lambda x: x.get('total',0))
    enriched_pathways['entities_found'] = enriched_pathways.entities.apply(lambda x: x.get('found',0))
    enriched_pathways['reactions_tot'] = enriched_pathways.reactions.apply(lambda x: x.get('total',0))
    enriched_pathways['reactions_found'] = enriched_pathways.reactions.apply(lambda x: x.get('found',0))
    
    enriched_pathways.sort_values('FDR', inplace=True)
    enriched_pathways = enriched_pathways[enriched_pathways.FDR < .05].copy(deep=True)
    print(len(enriched_pathways), 'pathway(s) enriched with FDR < 5%')
    if len(enriched_pathways) > 0:
        enriched_pathways[[
                'species_id','stId','name','FDR',
                'entities_tot','entities_found',
                'reactions_tot','reactions_found'
            ]].to_csv(outpath, index=False)
    print('https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + token)


# In[26]:
combo = []
for cluster,df in data.groupby(clustering_col).__iter__():
    try:
        tmp = pd.read_csv(os.path.join(FLD, f'cluster{cluster}-reac.csv')).head(8)
    except:
        continue
    tmp['cluster'] = f"C{cluster:02}"
    combo.append(tmp)
    del tmp
combo = pd.concat(combo, ignore_index=True)
combo.name = '[' + combo.stId + '] ' + combo.name.apply(lambda x: x[:70])
combo.to_csv(os.path.join(FLD,'results_combo.csv'), index=False)



