#!/usr/bin/env python
# coding: utf-8
import concurrent.futures
import pandas as pd
import polars as pl
import numpy as np
import scipy, os, time, dcor, math, argparse
from more_itertools import batched
from tqdm import tqdm
from datetime import datetime

worker_latent = None

def init_worker(latent_space):
    global worker_latent
    worker_latent = latent_space

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("model", type=str, help="Path to VAE folder.")
    return p.parse_args()

def calculate_correlations_w_pval(i):
    correlations = []
    print(i+1, '/', len(worker_latent)-1)
    ptmA = worker_latent.index[i]
    for j in range(i+1, len(worker_latent)):
        ptmB = worker_latent.index[j]
        arrayA, arrayB = worker_latent.iloc[i], worker_latent.iloc[j]
        with np.errstate(divide='ignore'):
            dpval, _ = dcor.independence.distance_correlation_t_test(arrayA, arrayB)
        pearson, _ = scipy.stats.pearsonr(arrayA, arrayB)
        d = dcor.distance_correlation(arrayA, arrayB)
        signd = math.copysign(d, pearson)
        correlations.append([ptmA, ptmB, signd, dpval, d, pearson])
    return pd.DataFrame(correlations, columns=['nodeA', 'nodeB', 'Score', 'pval', 'distance', 'PCC'])

def clean_tmp_files(list_of_files):
    for path in list_of_files:
        if os.path.exists(path):
            os.remove(path)
            # print(f"Removed {path}")


def main():
    # --- Setup ---
    start = time.perf_counter()
    date_and_time = datetime.today().strftime("%Y%m%d-%H%M")
    print(date_and_time)
    args = parse_cli()
    model_fld = args.model.rstrip('/')
    pickle_path = os.path.join(model_fld, 'Latent-space.pkl.gz')
    csv_path = os.path.join(model_fld, 'Latent-space.csv.gz')
    if os.path.exists(pickle_path):
        latent = pd.read_pickle(pickle_path)
    else:
        latent = pd.read_csv(csv_path, index_col=0)
    
    # --- Step 1: compute signed distances in batches ---
    partial_files = []
    try:
        batches = batched(range(len(latent)-1), n=1000)

        with concurrent.futures.ProcessPoolExecutor(
            initializer=init_worker,
            initargs=(latent,)
        ) as executor:
            for batch_id, batch in enumerate(batches):
                results = [_ for _ in executor.map(calculate_correlations_w_pval, batch) if len(_) > 0]
                if results:
                    correlations = pd.concat(results, ignore_index=True)
                else:
                    correlations = pd.DataFrame(columns=['nodeA', 'nodeB', 'Score', 'pval', 'distance', 'PCC'])

                savepath = os.path.join(model_fld, f"{model_fld.split('-')[-1]}-signed-distances-{batch_id}-partial.csv.gz")
                partial_files.append(savepath)
                correlations.to_csv(savepath, compression='gzip', index=False, encoding='utf-8')
                print(savepath)
                del correlations

        # --- Step 2: combine partial results ---
        combo = []
        for path in tqdm(partial_files):
            tmp = pl.read_csv(path, encoding='utf8')
            if len(tmp) > 0:
                combo.append(tmp)

        if not combo:
            raise ValueError('No correlations were produced; all partial result files were empty.')

        combo = pl.concat(combo)
        combo = combo.sort('pval')
        print(combo)

        combo = combo.to_pandas()
        combo['adj_pval'] = scipy.stats.false_discovery_control(combo.pval)

        combined_savepath = os.path.join(model_fld, f"{date_and_time}-{model_fld.split('-')[-1]}-signed-distances.csv.gz")
        combo[combo.adj_pval < .01].to_csv(combined_savepath, index=False, compression='gzip', encoding='utf8')
        print(f"Saved combined signed distances to {combined_savepath}")

        finish = time.perf_counter()
        print(f'Finished in {round(finish-start)} second(s)')

    finally:
        clean_tmp_files(partial_files)

if __name__ == '__main__':
    main()