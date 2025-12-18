#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import polars as pl 
import numpy as np
import scipy
import os, sys
from tqdm import tqdm
from datetime import datetime
date_and_time = datetime.today().strftime("%Y%m%d-%H%M")
print(date_and_time)

FLD = sys.argv[1]

# def combine_results(FLD):
to_be_combined = [
    _ for _ in os.scandir(FLD) if 'signed-distances' in _.name and _.name.endswith('partial.csv.gz')
    ]
combo = []
for _ in tqdm(to_be_combined):
    tmp = pl.read_csv(_.path, encoding='utf8')
    combo.append(tmp)
del tmp

combo = pl.concat(combo)
combo = combo.sort('pval')
print(combo)

combo=combo.to_pandas()
combo['adj_pval'] = scipy.stats.false_discovery_control(combo.pval)

savepath = os.path.join(FLD, f"{date_and_time}-{FLD.split('-')[-1]}-signed-distances.csv.gz")
combo[combo.adj_pval<.01].to_csv(savepath, index=False, compression='gzip', encoding='utf8')
print(f"Saved combined signed distances to {savepath}")

for _ in to_be_combined:
    os.remove(_.path)