#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Map_and_count_peptides import get_peptidoforms_counts_mapped
from Relative_counts import get_relative_ptms
from Generate_PTM_matrices import generate_ptm_submatrices
from Combine_PTM_matrices import combine_ptm_submatrices


# In[2]:


_prefix     = "Validation-20260812"
_fasta      = "./Human_2026_01_canonical.fasta.gz"
_pep_counts = "./v0113-2026/unprocessed/20260220_Peptidoforms_counts_v0113_breast_cancer.csv.gz"
# "./v0113-2026/unprocessed/20260220_Peptidoforms_counts_v0113_HamletAnnotated.csv.gz"
_pep_ids    = "./v0113-2026/unprocessed/20260220_Peptidoforms_IDs_v0113.csv.gz"
_myptms     = "./PTMs-of-interest-submission.csv"


# In[3]:


PATH = get_peptidoforms_counts_mapped(
    _pep_ids,
    _pep_counts,
    _fasta,
    p=_prefix
)


# In[4]:


PATH = get_relative_ptms(PATH)


# In[5]:


PATH = generate_ptm_submatrices(
    PATH,
    target_fasta_path="./Human_2026_01_canonical.fasta.gz",
    contam_fasta_path="./MQcontaminants_2023_11_14.fasta.gz",
    myptms=_myptms
)


# In[6]:


combine_ptm_submatrices(
    PATH,
    prefix=_prefix,
    myptms=_myptms
)


# -------------
