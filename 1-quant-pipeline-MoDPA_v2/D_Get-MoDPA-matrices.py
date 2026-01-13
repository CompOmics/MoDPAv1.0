#!/usr/bin/env python
# coding: utf-8
import polars as pl
import pandas as pd
import numpy as np
import argparse, re
from datetime import datetime


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("data_folder", type=str, help="Path to processed folder with the data.")
    p.add_argument("date", type=str, help="Date in YYYY-MM-DD format.")
    p.add_argument('--myptms', dest='myptms', type=str, default='./PTMs-of-interest.csv', 
               help="Path to list of PTMs to analyze. (default: ./PTMs-of-interest.csv)")
    return p.parse_args()

# In[6]:
def get_MoDPA_matrix_newer(relcounts_, myptms_):
    RES,ID,MOD = myptms_
    mod_of_interest = f"{RES}|[{ID}]{MOD}"
    
    relcounts_ = relcounts_[relcounts_.PTM_type==mod_of_interest].copy(deep=True)
    relcounts_.sort_values('ptm_name', inplace=True)
    
    try:
        return relcounts_.pivot(columns='file_name',index='PTM_ID',values='relative_psm_counts')
    except:
        pass
        
    relcounts_matrix = pd.DataFrame(relcounts_.groupby(['PTM_ID','file_name'])['relative_psm_counts'].max())
    ptms = list(relcounts_matrix.index.get_level_values(0).unique())
    len(ptms), ptms
    
    tmp = pd.DataFrame() 
    for i,j in enumerate(ptms):
        # clear_output(wait=True)
        tmp2 = relcounts_matrix.loc[j]
        if len(tmp2)>=10:
            tmp2.columns = [j]
            # print(f"{mod_of_interest} ({i+1}/{len(ptms)})")
            tmp = pd.concat((tmp,tmp2.T))
    
    tmp.index.names = ['PTM_ID']
    tmp = tmp.reset_index().sort_values('PTM_ID').set_index('PTM_ID').T
    tmp = tmp.reset_index().sort_values('file_name').set_index('file_name').T
    return tmp


# In[2]:
def main():
    args = parse_cli()
    
    relcounts = pl.read_csv(
        f"./{args.data_folder}/{args.date}_PTMs_counts_relative_prefiltered.csv.gz", 
        columns=['file_name','UniAcc','ptm_loc','ptm_res','ptm_name','classification','relative_psm_counts']
    )
    relcounts = relcounts.with_columns(
        pl.struct(['UniAcc','ptm_loc','ptm_res','ptm_name']).map_elements(
            lambda row: f"{row['UniAcc']}|{row['ptm_loc']}|{row['ptm_res']}|{row['ptm_name']}"
        ).alias('PTM_ID'),
        pl.col('ptm_name').map_elements(
            lambda x: int(re.match(r'\[(\d+)\]',x).groups()[0]),
            return_dtype=pl.Int32
        ).alias('unimod_id'),
        (pl.col('ptm_res') + '|' + pl.col('ptm_name')).alias('PTM_type')
    )
    print('#Raw files =', relcounts.select('file_name').n_unique() )
    print('#Modifications =', relcounts.select('PTM_ID').n_unique() )
    relcounts = relcounts.to_pandas()
    print(relcounts.ptm_name.value_counts())
    
    
    # In[5]:
    myptms = pl.read_csv(args.myptms) # csv with aminoacid, unimod id, and name
    myptms = myptms.rows()
    print(myptms)
    
    
    # In[7]:
    for n,mod in enumerate(myptms):
        savepath = f'{args.data_folder}/MoDPA_Rel_[{mod[1]}]{mod[2]}_{mod[0]}.pkl.gz'
        
        relcounts_matrix = get_MoDPA_matrix_newer(relcounts, mod)
        
        if len(relcounts_matrix)>0:
           relcounts_matrix.T.to_pickle(savepath, compression='gzip')
           print(n+1, '/', len(myptms), '\t', 
                 savepath, relcounts_matrix.shape, relcounts_matrix.notna().sum().sum())
        else:
            print(n+1, '/', len(myptms), '\t', 
                  mod, relcounts_matrix.shape, relcounts_matrix.notna().sum().sum())
        

if __name__ == "__main__":
    START = datetime.now()

    main()
    
    END = datetime.now()
    print("Done!!")
    print("Started: ", START.isoformat())
    print("Finished:", END.isoformat(), '\n')