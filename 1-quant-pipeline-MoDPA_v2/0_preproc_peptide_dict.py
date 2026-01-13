#!/usr/bin/env python
# coding: utf-8
from PTMmap import Fasta
import gzip, re, argparse
import pyteomics.parser
import pandas as pd
from tqdm import tqdm


# In[2]:
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('fasta', metavar='fasta', type=str, help="Path to gzipped .fasta file")
    return p.parse_args() 

def getClass(acc, fasta) -> str:
    '''to sort canonical prots, isoforms, ORFs, etc.'''
    try:
        return fasta[acc]['Class']
    except:
        return 'zzz'


def get_peptide_identifier(prot, uniprot_id, seq, peptide):
    return "%s((%i-%i))((%s))" % (
        prot,
        seq.find(peptide) + 1,
        seq.find(peptide) + len(peptide) + 1,
        uniprot_id,
    )


def add_pep_to_dict(pep_dict, peptide, fasta, UniAcc, sequence):
    if peptide in pep_dict:
        pep_dict[peptide].append(get_peptide_identifier(fasta[UniAcc]['Entry'], UniAcc, sequence, peptide))
    else:
        pep_dict[peptide] = [get_peptide_identifier(fasta[UniAcc]['Entry'], UniAcc, sequence, peptide)]


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


# In[ ]:
args = parse_cli()
fasta = {}
Fasta.getFasta(args.fasta, fasta)
Fasta.addClassification(fasta)


# In[35]:
proteins = list(fasta.keys())
pep_dict = {}

for UniAcc in tqdm(proteins, total=len(proteins)):
    sequence = fasta[UniAcc]['Seq']
    
    for peptide in pyteomics.parser.cleave(sequence, 'Trypsin/P', 2):
        # 'Trypsin/P' cleaves after K/R, ignores Proline
        
        if len(peptide) < MIN_LENGTH or len(peptide) > MAX_LENGTH:
            continue
        if non_aa_characters.findall(peptide):
            continue

        add_pep_to_dict(pep_dict, peptide, fasta, UniAcc, sequence)

        if peptide.startswith("M") and sequence.startswith(peptide) and len(peptide)-1 >= MIN_LENGTH:
            add_pep_to_dict(pep_dict, peptide[1:], fasta, UniAcc, sequence)

pep_dict = pd.DataFrame(
    [(peptide,'||'.join(ids)) for peptide,ids in pep_dict.items()],
    columns=['peptide','proteins']
)
pep_dict.to_csv(
    'peptide_dict_new.csv.gz', index=False, compression='gzip', encoding='utf-8'
)

