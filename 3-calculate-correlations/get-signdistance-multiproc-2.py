#!/usr/bin/env python
# coding: utf-8
import concurrent.futures
import multiprocessing as mp
import pandas as pd
import numpy as np
import scipy, os, sys, time, dcor, math
from more_itertools import batched 
start = time.perf_counter()

FLD = sys.argv[1]
try:
    latent = pd.read_pickle(os.path.join(FLD,'Latent-space.pkl.gz'))
except:
    latent = pd.read_csv(os.path.join(FLD,'Latent-space.csv.gz'), index_col=0)
# print(latent.dtypes)

def calculate_correlations_w_pval(i, latent_space=latent):
    correlations = []
    print(i+1,'/',len(latent_space)-1)
    ptmA = latent_space.index[i]
    for j in range(i+1, len(latent_space)):
        ptmB = latent.index[j]
        arrayA, arrayB = latent.iloc[i], latent.iloc[j]
        with np.errstate(divide='ignore'):
            dpval, _ = dcor.independence.distance_correlation_t_test(arrayA, arrayB)
        if dpval < .05:
            pearson, _ = scipy.stats.pearsonr(arrayA, arrayB)
            d = dcor.distance_correlation(arrayA, arrayB)
            signd = math.copysign(d,pearson)
            correlations.append([ptmA, ptmB, f'{d:.2f}', f'{signd:.2f}', f'{dpval:.4f}'])
    return pd.DataFrame(correlations, columns=['nodeA','nodeB','distance','Score','pval'])

if __name__=='__main__':
    # mp.freeze_support()
    
    batches = batched(range(len(latent)-1), n=1000)

    for batch_id,batch in enumerate(batches):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(calculate_correlations_w_pval, batch)
            correlations = pd.concat(results, ignore_index=True)
            
        savepath = os.path.join(FLD,f"{FLD.split('-')[-1]}-signed-distances-{batch_id}-partial.csv.gz")
        correlations.to_csv(savepath, compression='gzip', index=False)
        print(savepath)
        del correlations

finish = time.perf_counter()
print(f'Finished in {round(finish-start, 2)} second(s)')
