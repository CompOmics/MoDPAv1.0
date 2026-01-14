#!/usr/bin/env python
# coding: utf-8
from reactome2py import content, analysis
import os, time, sys, argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("data_path", type=str, help="Path to the clustered nodes CSV file.")
    p.add_argument('clustering_col', type=str, help="Column with cluster names/indices.")
    return p.parse_args()

args = parse_cli()
data_path = args.data_path
clustering_col = args.clustering_col

if not os.path.isfile(data_path):
    sys.exit(f"Error: data file not found: {data_path}")

FLD = os.path.split(data_path)[0]
data = pd.read_csv(data_path)

data.rename(columns={
    'shared name':'PTM_ID', 
    'mod-label':'modification',
    'RES':'residue',
    'POS':'position',
    'MOD':'Unimod_ID'
    }, inplace=True)
print(data.shape)
data[ # saves the nodes in a slightly better format
    ['PTM_ID',clustering_col,'UniAcc','residue','position','modification','Unimod_ID']
].to_csv(os.path.join(FLD, 'nicely-formatted-nodes.csv'))

for cluster,df in data.groupby(clustering_col).__iter__():
    outpath  = os.path.join(FLD, f'cluster{cluster}-reac.csv')
    outpath2 = os.path.join(FLD, f'cluster{cluster}.csv')
    markers = list(set(df.UniAcc))
    
    df[ # save the nodes of the clusters into separate .csv files
        ['PTM_ID',clustering_col,'UniAcc','residue','position','modification','Unimod_ID']
    ].to_csv(outpath2, index=False)

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



