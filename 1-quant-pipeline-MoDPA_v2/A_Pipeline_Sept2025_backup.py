#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import polars as pl
import numpy as np
from PTMmap import Fasta
import argparse
from datetime import date, datetime
TODAY = date.today().isoformat() 
print('The date is:',TODAY)


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('peptidoform_ids', metavar='peptidoform_ids', type=str, 
                        help="Path to 'Peptidoforms' file.")
    p.add_argument('peptidoform_counts', metavar='peptidoform_counts', type=str, 
                        help="Path to 'peptidoform_counts' file.")
    p.add_argument('peptides_mappings', metavar='peptides_mappings', type=str, 
                        help="Path to 'peptides_mappings' file.")
    p.add_argument('fasta', metavar='fasta', type=str, 
                        help='Path to FASTA file.')
    return p.parse_args()


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

def read_peptide_to_protein_mappings(pep_dict_path, pep_set, fasta_path):
    # read Fasta to classify proteins into canonical, isoforms, ecc
    fasta = {}
    Fasta.getFasta(fasta_path, fasta)
    Fasta.addClassification(fasta)
    # read ionbot pep dict file
    pepdict = pl.scan_csv(pep_dict_path)
    pepdict = pepdict.rename({'peptide':'sequence'})
    # print(f"#Unique pep sequences in searchDB = {len(pepdict.select(pl.col('sequence')).unique()):,}")
    pepdict = pepdict.filter(pl.col('sequence').is_in(pep_set))
    # print(f"#Unique pep sequences found by Ionbot = {len(pepdict.select(pl.col('sequence')).unique()):,}")

    pepdict = pepdict.with_columns(pl.col("proteins").str.split(by="||"))
    pepdict = pepdict.explode('proteins')
    pepdict = pepdict.with_columns(
        pl.col("proteins").str.replace_all(r"\(\(|\)\)", "//")
        )
    pepdict = pepdict.with_columns(
        pl.col('proteins').str.split('//').list.get(3).str.replace_all(r"\|", ".").alias('UniAcc'),
        pl.col('proteins').str.split('//').list.get(0).alias('entry'),
        pl.col('proteins').str.split('//').list.get(1).str.split('-').list.get(0).cast(pl.Int32).alias('start')
        )
    pepdict = pepdict.sort(
        pl.col('UniAcc').map_elements(lambda x: getClass(x,fasta))
        )
    pepdict = pepdict.group_by('sequence').agg(
        [pl.col('start'), pl.col('UniAcc'), pl.col('entry')]
        )
    pepdict = pepdict.with_columns(
        (pl.col('UniAcc').list.len() > 1).alias('ambiguous_map')
        )
    return pepdict.collect()
    
# def get_leading(row):
#     counts = list(row['unique_counts'])
#     accs = list(row['UniAcc'])
#     entries = list(row['entry'])
#     starts = list(row['start'])
#     leading_prot = accs[counts.index(max(counts))]
#     leading_entry = entries[counts.index(max(counts))]
#     leading_start = starts[counts.index(max(counts))]
#     return leading_prot, leading_entry, leading_start

def get_ambiguous_peptides(pep_dict_path, pep_set, fasta_path):
    maps_ambig_partial = read_peptide_to_protein_mappings(pep_dict_path, pep_set, fasta_path)    
    unique_mappings = maps_ambig_partial.filter(~pl.col('ambiguous_map'))
    unique_mappings = unique_mappings.with_columns(
        pl.col('UniAcc').list.get(0).alias('unique_protein')
    )

    unique_counts = unique_mappings.group_by("unique_protein").len()
    unique_counts = dict(zip(unique_counts['unique_protein'], unique_counts['len']))
    
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('UniAcc').map_elements(lambda x: unique_counts_per_UniAcc(x,unique_counts)).alias('unique_counts')
    )
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('unique_counts').map_elements(np.argmax).alias('max_idx')
    )
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('UniAcc').list.get(pl.col('max_idx')).alias('LeadProt'),
        pl.col('entry').list.get(pl.col('max_idx')).alias('LeadEntry'),
        pl.col('start').list.get(pl.col('max_idx')).alias('pep_start') 
    )
    maps_ambig_partial = maps_ambig_partial.with_columns(
        pl.col('UniAcc').map_elements(lambda x: '||'.join(list(x))),
        pl.col('unique_counts').map_elements(lambda x: '_'.join([str(y) for y in x]))
    )
    maps_ambig_partial = maps_ambig_partial.sort('ambiguous_map')
    return maps_ambig_partial.select(['sequence','pep_start','LeadProt','LeadEntry','UniAcc','unique_counts','max_idx'])

def map_ionbot_IDs(IDs_path, pep_dict_path, psm_counts_path, fasta_path):
    my_ids = pl.read_csv(psm_counts_path,columns=['peptidoform_id'])
    my_ids = my_ids['peptidoform_id'].to_list()
    ids = pl.scan_csv(IDs_path)
    ids = ids.filter(
        pl.col('peptidoform_id').is_in(my_ids)
    ).collect()
    del my_ids
    
    maps_ambig = get_ambiguous_peptides(
        pep_dict_path, 
        ids['sequence'].unique().to_list(),
        fasta_path
    ) 
    maps_ambig.columns = ['sequence','pep_start','LeadProt','LeadEntry','all_UniAcc','unique_counts','max_idx']
    
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
                                   'sequence','LeadProt','LeadEntry','all_UniAcc','pep_start','unique_counts','max_idx'])
    for (peptidoform_id,peptide_id,is_modified,seq,leadprot,leadentry,allprots,pepstart,unique_counts,max_idx),df in iterator.__iter__():
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
            unique_counts,max_idx
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
# CODE 
# ---------------
START = datetime.now()
print(START.isoformat())
args = parse_cli()
mapped_ids = map_ionbot_IDs(
    args.peptidoform_ids, 
    args.peptides_mappings,
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
mapped_peptidoforms_counts.to_csv(f'{TODAY}_Peptidoforms_counts_mapped_backup.csv.gz', compression='gzip', index=False, encoding='utf-8')

END = datetime.now()
print("Done!!")
print("Started: ", START.isoformat())
print("Finished:", END.isoformat(), '\n')