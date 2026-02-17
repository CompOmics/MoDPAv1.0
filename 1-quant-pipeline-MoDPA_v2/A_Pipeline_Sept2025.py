#!/usr/bin/env python
# coding: utf-8
import os, sys
import pandas as pd
import polars as pl
import numpy as np
from PTMmap import Fasta
import argparse, re, time
import pyteomics.parser
from tqdm import tqdm

from datetime import date, datetime
TODAY = date.today().isoformat() 
print('The date is:',TODAY)


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('peptidoform_ids', type=existing_file, help="Path to 'Peptidoforms_IDs' file")
    p.add_argument('peptidoform_counts', type=existing_file, help="Path to 'peptidoform_counts' file")
    p.add_argument('fasta', type=existing_file, help='Path to FASTA file')
    p.add_argument("-t", "--threads", type=int, default=16, help="Number of threads to use (Default: 16)")
    return p.parse_args()

def existing_file(path: str) -> str:
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found! --> {path}")
    else:
        return path


# ----------
# define functions
# ----------
# New Polars implementation
def getClass(acc, fasta) -> str:
    '''to sort canonical prots, isoforms, ORFs, etc.'''
    try:
        return fasta[acc]['Class']
    except:
        return 'zzz'

def unique_counts_per_UniAcc(prot_list, counts):
    """Map a sequence of UniProt accessions to their unique-peptide counts, defaulting to 0."""
    return [counts.get(acc, 0) for acc in prot_list]

def read_peptide_to_protein_mappings(fasta_path, pep_set):
    ## CONSTANTS
    MIN_PMZ = 500
    MAX_PMZ = 6000
    MIN_LENGTH = 7
    MAX_LENGTH = 30
    AMINO_MASSES = {
        "A": 71.037114,
        "C": 103.009185,
        "D": 115.026943,
        "E": 129.042593,
        "F": 147.068414,
        "G": 57.021464,
        "H": 137.058912,
        "I": 113.084064,
        "K": 128.094963,
        "L": 113.084064,
        "M": 131.040485,
        "N": 114.042927,
        "P": 97.052764,
        "Q": 128.058578,
        "R": 156.101111,
        "S": 87.032028,
        "T": 101.047679,
        "V": 99.068414,
        "W": 186.079313,
        "Y": 163.063329,
    }
    AMINO_MAP = {
        a: i for i, a in enumerate(sorted((a for a in AMINO_MASSES.keys() if a != "L")))
    }
    AMINO_MAP["L"] = AMINO_MAP["I"]
    non_aa_characters = re.compile(
        r'[^%s]' % ("".join(list(AMINO_MASSES.keys())))
    ) #Matches anything that is NOT a valid AA
    
    # read Fasta to classify proteins into canonical, isoforms, ecc
    fasta = {}
    Fasta.getFasta(fasta_path, fasta)
    Fasta.addClassification(fasta)
    proteins = list(fasta.keys())
    pep_df = []
    print("Parsing .fasta file...")
    for UniAcc in tqdm(proteins, total=len(proteins)):
        sequence = fasta[UniAcc]['Seq']
    
        for peptide in pyteomics.parser.cleave(sequence, 'Trypsin/P', 2):
            # 'Trypsin/P' cleaves after K/R, ignores Proline
            if len(peptide) < MIN_LENGTH or len(peptide) > MAX_LENGTH or non_aa_characters.findall(peptide):
                continue
            # if peptide not in pep_set:
            #     continue
            if peptide.startswith("M") and sequence.startswith(peptide) and len(peptide)-1 >= MIN_LENGTH:
                pep_df.append((peptide[1:], UniAcc, fasta[UniAcc]['Entry'], sequence.find(peptide[1:]) + 1))
            pep_df.append((peptide, UniAcc, fasta[UniAcc]['Entry'], sequence.find(peptide) + 1))
    
    pep_df = pl.DataFrame(
        pep_df, orient="row",
        schema={'sequence':pl.String, 'UniAcc':pl.String, 'entry':pl.String, 'start':pl.Int64}
        )
    pep_df = pep_df.filter(
        [True if pfid in pep_set else False for pfid in pep_df["sequence"]]
    ).lazy()
    pep_df = pep_df.filter(pl.col('sequence').is_in(pep_set))
    pep_df = pep_df.sort('UniAcc')
    pep_df = pep_df.sort(
        pl.col('UniAcc').map_elements(lambda x: getClass(x,fasta))
        )
    pep_df = pep_df.group_by('sequence').agg(
        [pl.col('start'), pl.col('UniAcc'), pl.col('entry')]
        )
    pep_df = pep_df.with_columns(
        (pl.col('UniAcc').list.len() > 1).alias('ambiguous_map')
        )
    return pep_df.collect()
    

def get_ambiguous_peptides(fasta_path, pep_set):
    maps_ambig_partial = read_peptide_to_protein_mappings(fasta_path, pep_set)
    # maps_ambig_partial = maps_ambig_partial.filter(pl.col('sequence').is_in(pep_set))
    unique_mappings = maps_ambig_partial.filter(~pl.col('ambiguous_map'))

    print("Calculating unique peptide counts per protein...")
    
    unique_mappings = unique_mappings.with_columns(
        pl.col('UniAcc').list.get(0).alias('unique_protein')
    )

    unique_counts = unique_mappings.group_by("unique_protein").len()
    unique_counts = dict(zip(unique_counts['unique_protein'], unique_counts['len']))
    
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('UniAcc').map_elements(lambda x: unique_counts_per_UniAcc(x,unique_counts)).alias('unique_counts')
    )

    print("Selecting leading proteins...")
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('unique_counts').map_elements(np.argmax).alias('max_idx')
    ).with_columns(
        pl.col('UniAcc').list.get(pl.col('max_idx')).alias('LeadProt'),
        pl.col('entry').list.get(pl.col('max_idx')).alias('LeadEntry'),
        pl.col('start').list.get(pl.col('max_idx')).alias('pep_start') 
    )
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('UniAcc').map_elements(lambda x: '||'.join(list(x))),
        pl.col('unique_counts').map_elements(lambda x: '_'.join([str(y) for y in x]))
    )
    maps_ambig_partial = maps_ambig_partial.sort('ambiguous_map')
    return maps_ambig_partial.select(['sequence','pep_start','LeadProt','LeadEntry','UniAcc'])

def map_ionbot_IDs(IDs_path, psm_counts_path, fasta_path):
    my_ids = pl.read_csv(psm_counts_path,columns=['peptidoform_id'])
    my_ids = set(my_ids['peptidoform_id'].to_list())
    ids = pl.read_csv(IDs_path)
    ids = ids.filter(
        [True if pfid in my_ids else False for pfid in ids["peptidoform_id"]]
        )
    del my_ids

    print("Mapping peptides to proteins...")
    maps_ambig = get_ambiguous_peptides(
        fasta_path, 
        set(ids['sequence'].to_list())
    ) 
    maps_ambig.columns = ['sequence','pep_start','LeadProt','LeadEntry','all_UniAcc']
    
    print("IDs table size:", ids.shape)
    mapped_ids = ids.join(maps_ambig, on='sequence', validate='m:1')
    print("Mapped_IDs table size:",mapped_ids.shape)
    print(f"({ ids.shape[0]-mapped_ids.shape[0] } IDs discarded)")
    mapped_ids = mapped_ids.with_columns(
        pl.col('ptm_loc') + pl.col('pep_start') - 1
    )
    return mapped_ids


## PANDAS IMPLEMENTATION ##
def group_IDs_into_peptidoforms(mapped_ids):
    mapped_ids_grouped = []
    iterator = mapped_ids.groupby(['peptidoform_id','peptide_id','is_modified',
                                   'sequence','LeadProt','LeadEntry','all_UniAcc','pep_start'])
    for (peptidoform_id,peptide_id,is_modified,seq,leadprot,leadentry,allprots,pepstart),df in iterator.__iter__():
        mapped_ids_grouped.append([
            peptidoform_id,
            peptide_id,
            is_modified=='t',
            seq,
            ';'.join(list(df.ptm_name.apply(str))),
            ';'.join(list(df.ptm_loc.apply(str))),
            ';'.join(list(df.ptm_res.apply(str))),
            ';'.join(list(df.classification.apply(str))),
            pepstart,
            leadprot,
            leadentry,
            allprots,
        ])
    return pd.DataFrame(mapped_ids_grouped, columns=mapped_ids.columns)

def add_psm_counts(mapped_peptidoforms, psm_counts_path):
    counts  = pd.read_csv(psm_counts_path, usecols=['file_name','peptidoform_id','psm_counts'])
    print(counts.shape)
    print(mapped_peptidoforms.shape)
    counts = counts.merge(mapped_peptidoforms, on='peptidoform_id')
    print(counts.shape)
    return counts

def psm_counts_per_PTM(mapped_peptidoforms_counts):
    modcounts = []
    for _,row in mapped_peptidoforms_counts.iterrows():
        iterator = zip(row.ptm_loc, row.ptm_name, row.ptm_res, row.classification)
        for p,m,r,c in iterator:
            if c=='ragging':
                continue
            tmp = [
                row.file_name, row.peptidoform_id, row.psm_counts, row.sequence, row.is_modified,  
                row.LeadProt, row.LeadEntry, row.all_UniAcc,
                p, m, r, c
            ]
            modcounts.append(tmp)
    
    cols = [
        'file_name',
        'peptidoform_id',
        'psm_counts',
        'sequence',
        'is_modified',
        'LeadProt',
        'LeadEntry',
        'all_UniAcc',
        'ptm_loc',
        'ptm_name',
        'ptm_res',
        'ptm_class',
    ]
    modcounts = pd.DataFrame(modcounts, columns=cols)
    
    modcounts = modcounts.groupby(['file_name','LeadProt','LeadEntry',
                                   'ptm_loc','ptm_name','ptm_res','ptm_class']).sum()[['psm_counts']]
    modcounts.columns = ['total_counts']
    print('\n',modcounts.describe(),'\n')
    return modcounts.reset_index()


# ---------------
# MAIN CODE 
# ---------------
start = time.perf_counter()
args = parse_cli()
os.environ["POLARS_MAX_THREADS"] = f"{args.threads}" # to stop polars from using all available memory
mapped_ids = map_ionbot_IDs(
    args.peptidoform_ids, 
    args.peptidoform_counts,
    args.fasta
)
print(mapped_ids.head())
mapped_ids = mapped_ids.to_pandas()
peptidoforms = group_IDs_into_peptidoforms(mapped_ids)
mapped_peptidoforms_counts = add_psm_counts(peptidoforms, args.peptidoform_counts)

print(f"#peptidoforms = {len(set(peptidoforms.peptidoform_id)):,}")
print("Unique peptidoforms:", peptidoforms.shape)

# Safer counting to avoid unpacking errors
vc = peptidoforms.is_modified.value_counts()
mod = int(vc.get(True, 0))
unmod = int(vc.get(False, 0))
print('_')
print(vc)
total = mod + unmod
if total > 0:
    print(f"% Unmodified peptides = {unmod / total:.1%}", '\n')

print("Mapped peptidoforms with counts:", mapped_peptidoforms_counts.shape)
mapped_peptidoforms_counts.to_csv(f'{TODAY}_Peptidoforms_counts_mapped.csv.gz', compression='gzip', index=False, encoding='utf-8')

finish = time.perf_counter()
print(f'Finished in {round(finish-start, 2)} second(s)')