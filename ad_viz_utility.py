import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw
import tarfile
import pdb
import math

# Writing in files
filtered_reactome_gene_sets = pd.read_csv('data/gene_sets/filtered_reactome_gene_sets.csv', index_col=0)

def load_GO_data(cell_type):
    filtered_GO_gene_sets = pd.read_csv(f'data/gene_sets/filtered_GO_{cell_type}.csv', index_col=0)
    return filtered_GO_gene_sets
    
# Functions

# --- DATA LOADING (CACHED) ---
@st.cache_data
def load_data(data_path, samp_meta_path, entry_meta_path, sample_col):
    # load CSVs
    data = pd.read_csv(data_path, index_col=0)
    samp_meta = pd.read_csv(samp_meta_path, index_col=0)
    entry_meta = pd.read_csv(entry_meta_path, index_col=0)
    
    ## formatting
    # transposes the table so that samples are rows
    data = data.T

    # merge with relevant sample metadata columns
    data = data.merge(samp_meta[[
        'Parental.line', 'CRISPR.edits', sample_col]],
        left_index = True, right_index = True)

    return data, samp_meta, entry_meta
    
def read_targz_to_pandas(file, tar_path):
    with tarfile.open(tar_path, "r:gz") as tar:
        file_obj = tar.extractfile(file)
        dataframe = pd.read_csv(file_obj)
    return dataframe

def pathway_to_uniprot(database, db_name, pathway_identifier, pathway_name, protein_identifier, entry_meta):
    if db_name == 'GO':
        filtered = database[database[pathway_identifier] == pathway_name]
        proteins = list(filtered[protein_identifier])
    elif db_name == 'reactome':
        filtered = database[database[pathway_identifier] == pathway_name]
        genes = filtered[protein_identifier]  
        match = entry_meta[entry_meta['id.mapped'].isin(genes)]
        proteins = list(match.index.values)
    return proteins

def lookup_by_protid(protein_id: str, protein_metadata: pd.DataFrame,
                     name_col: str):
    '''
    Uses a protein ID to look up the rows of the `protein_metadata` table
    where that ID matches the index. Returns the contents of `name_col`
    for that row.
    For proteomics, this is the uniprot ID.
    Args:
    - protein_id (str): The protein identifier to look up in the index
    - protein_metadata (DataFrame): the protein metadata dataframe
    - name_col (str): Column name of the identifier/descriptor to look up
    Returns:
    - str
    '''

    # returns np.nan if the protein_id is not found.
    if protein_id not in protein_metadata.index:
        return np.nan
    # otherwise returns the contents of name_col for that row
    # assumes only one row will be returned (index is unique)
    return protein_metadata.loc[protein_id, name_col]

def file_keys(comparison, enrichment, database, cell_type):
    #comparison
    if comparison == "ABCA7LOF vs APOE3":
        comp_key = "ABCA7LOFvnormal"
    elif comparison == "APOE4 vs APOE3":
        comp_key = "APOE4v3"

    #enrichment
    if enrichment == "GSEA":
        csv_ending = "logFC_pval_table"
    elif enrichment == "ORA":
        csv_ending = "enrichment_table"
        
    #database
    if database == "reactome":
        pathway_identifier = "clean_name"
        gene_set = filtered_reactome_gene_sets 
        protein_identifier = "gene"
    elif database == "GO":
        pathway_identifier = "TERM"
        gene_set = load_GO_data(cell_type)
        protein_identifier = "UNIPROT"

    return comp_key, csv_ending, pathway_identifier, gene_set, protein_identifier

def i_to_row_col(i, cols):
    row = i // cols + 1 #returns division
    col = i % cols + 1 #returns remainder
    return row, col
