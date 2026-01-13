#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import polars as pl
import numpy as np
import argparse, os, shutil
from string import ascii_uppercase
from datetime import datetime

# In[]: def functions
def ptm_dictionary(row):
    x = {}
    if row['is_modified']:
        for p,m,r,c in zip(row['ptm_loc'],row['ptm_name'],row['ptm_res'],row['classification']):
            if p==p:
                x[int(float(p))] = {'mod': m, 'res': r, 'cla': c}
            else:
                print(row)
    return x
    
def map_modifications_to_peptide(row, unmod_label='[0]Unmod'):
    mods,classes = [],[]
    ptm_dict = ptm_dictionary(row)
    # return ptm_dict
    for _ in row['ptm_loc_2']:
        try:
            mods.append(ptm_dict[_]['mod'])
            classes.append(ptm_dict[_]['cla'])
        except KeyError:
            mods.append(unmod_label)
            classes.append(unmod_label)
    return mods,classes

def zip_modifications(row):
    iterator = zip(row['ptm_loc'], row['ptm_name'], row['ptm_res'], row['classification'])
    prot = row['LeadProt']
    return [f'{prot}|{p}|{r}|{m}|{c}' for p,m,r,c in iterator if r in list(ascii_uppercase) and c not in ['ragging','semi_tryptic']]


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("data_folder", type=str, help="Path to processed folder with the data.")
    p.add_argument("date", type=str, help="Date in YYYY-MM-DD format.")
    return p.parse_args()

def process_chunk(data: pd.DataFrame, chunk_path: str) -> None:
    lf = data.lazy()

    # print("Step 1: Filtering...")
    lf = lf.with_columns(
        pl.col('ptm_name').str.split(';'),
        pl.col('ptm_res').str.split(';'),
        pl.col('classification').str.split(';'),
        pl.col('ptm_loc').str.split(';'),
        )

    # print("Step 2: Getting total counts...")
    # lf = lf.with_columns(
    #     pl.col('psm_counts').sum().over(['file_name','sequence']).alias('TOT_seq_counts')
    # )

    # print("Step 3: Filling unmodified residues...")
    lf = lf.with_columns(
        pl.col("sequence").map_elements(
            list, 
            return_dtype=pl.List(pl.String)
        ).alias("ptm_res_2"),
        pl.col("sequence").map_elements(
            lambda x: list(range(len(x))),
            return_dtype=pl.List(pl.Int32)
        ).alias("ptm_loc_2"),
    )
    # ptm_loc_2: add pep_start to each position
    lf = lf.with_columns(
        pl.struct(["ptm_loc_2", "pep_start"]).map_elements(
            lambda row: [i + row["pep_start"] for i in row["ptm_loc_2"]]
        ).alias("ptm_loc_2")
    )
    # ptm_name_2, classification_2: map modifications to peptide
    lf = lf.with_columns(
        pl.struct(["is_modified","ptm_loc","ptm_name","ptm_res","classification","ptm_loc_2"]).map_elements(
            map_modifications_to_peptide,
        ).alias("tmp")
    )
    lf = lf.with_columns([
        # Split tuple output into two columns
        pl.col("tmp").list.get(0).alias("ptm_name_2"),
        pl.col("tmp").list.get(1).alias("classification_2"),
    ]).drop("tmp")
    # Overwrite original columns
    for i in ['ptm_name','ptm_loc','ptm_res','classification']:
        lf = lf.with_columns([pl.col(f"{i}_2").alias(i)])
    # Drop temporary columns
    lf = lf.drop(['pep_start','ptm_res_2','ptm_name_2','ptm_loc_2','classification_2'])
    lf = lf.explode(['ptm_loc','ptm_name','ptm_res','classification'])


    # print("...Calculating relative counts...")
    # A: Group by all columns and sum psm_counts
    grouped_counts = lf.group_by(['file_name', 'LeadProt', 'ptm_loc', 'ptm_res', 'ptm_name', 'classification']).agg([
        pl.col('psm_counts').sum().alias('ptm_psm_counts')
    ]).filter(
        pl.col('ptm_psm_counts') >= 3
    )
    # B: Calculate total psm_counts per file-LeadProt-ptm_loc-ptm_res group
    totals_per_location = lf.group_by(['file_name', 'LeadProt', 'ptm_loc', 'ptm_res']).agg([
        pl.col('psm_counts').sum().alias('total_site_psm_counts')
    ])
    # C: Join back to get relative counts
    grouped_counts = grouped_counts.collect()
    totals_per_location = totals_per_location.collect()
    rel_counts = grouped_counts.join(
        totals_per_location, 
        on=['file_name', 'LeadProt', 'ptm_loc', 'ptm_res']
    ).with_columns([
        (pl.col('ptm_psm_counts') / pl.col('total_site_psm_counts')).alias('relative_psm_counts')
    ])

    rel_counts.write_csv(chunk_path, float_precision=3)

# In[]: main
def main():
    args = parse_cli() 
    in_path = f"{args.data_folder}/{args.date}_Peptidoforms_counts_mapped.csv.gz"
    tmp_dir_path = f"{args.data_folder}/counts-per-msrun"

    lf = pl.scan_csv(
        in_path,
        encoding='utf8'
    ).filter(
        (~pl.col('ptm_res').str.contains("N-TERM", literal=True))&
        (~pl.col('classification').str.contains("semi_tryptic", literal=True))&
        (~pl.col('classification').str.contains("ragging", literal=True))&
        (pl.col('psm_counts') >= 2)
    )
    data_chunks = lf.collect().partition_by("file_name")

    os.makedirs(tmp_dir_path, exist_ok=False)
    n_chunks = len(data_chunks)

    chunk_path_list = []
    for i in range(n_chunks):
        print(i+1, "/", n_chunks)
        chunk_path = f"{tmp_dir_path}/msrun-{i+1}.csv"
        data = data_chunks.pop()
        process_chunk(data, chunk_path)
        chunk_path_list.append(chunk_path)

    final_data = []
    for _ in chunk_path_list:
        final_data.append(pl.read_csv(_))
    final_data = pl.concat(final_data)
    final_data = final_data.to_pandas()
    
    final_data.rename(columns={'LeadProt':'UniAcc'}, inplace=True)
    final_data.to_csv(
        f"{args.data_folder}/{args.date}_PTMs_counts_relative.csv.gz",
        index=False, 
        compression='gzip',
        encoding='utf-8'
    )
    try:
        shutil.rmtree(tmp_dir_path)
    except:
        print('!! could not remove temporary folder counts-per-msrun')


if __name__ == "__main__":
    START = datetime.now()
    print(START.isoformat())

    main()
    
    END = datetime.now()
    print("Done!!")
    print("Started: ", START.isoformat())
    print("Finished:", END.isoformat(), '\n')