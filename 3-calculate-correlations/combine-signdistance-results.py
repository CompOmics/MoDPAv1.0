#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import polars as pl 
import numpy as np
import scipy
import os, argparse
from tqdm import tqdm
from datetime import datetime
date_and_time = datetime.today().strftime("%Y%m%d-%H%M")
print(date_and_time)

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("model", type=str, help="Path to VAE folder.")
    return p.parse_args()

args = parse_cli()

to_be_combined = [
    _ for _ in os.scandir(args.model) if 'signed-distances' in _.name and _.name.endswith('partial.csv.gz')
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

savepath = os.path.join(args.model, f"{date_and_time}-{args.model.split('-')[-1]}-signed-distances.csv.gz")
combo[combo.adj_pval<.01].to_csv(savepath, index=False, compression='gzip', encoding='utf8')
print(f"Saved combined signed distances to {savepath}")

for _ in to_be_combined:
    os.remove(_.path)