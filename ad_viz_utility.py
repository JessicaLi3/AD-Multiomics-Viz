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

def file_keys(comparison, enrichment):
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
    elif enrichment == "All":
        csv_ending = None
        
    #database
    #if database == "reactome":
        #pathway_identifier = "clean_name"
        #gene_set = filtered_reactome_gene_sets 
        #protein_identifier = "gene"
    #elif database == "GO":
        #pathway_identifier = "TERM"
        #gene_set = load_GO_data(cell_type)
        #protein_identifier = "UNIPROT"

    return comp_key, csv_ending

def load_enrichment_df(cell_type, comp_key, enrichment, database, selection):
    csv_ending_map = {"GSEA": "logFC_pval_table", "ORA": "enrichment_table"}

    if database == "All" and enrichment == "All":
        go_gsea = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_GSEA_GO_logFC_pval_table.csv', selection=selection)
        go_ora = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_ORA_GO_enrichment_table.csv', selection=selection)
        reactome_gsea = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_GSEA_reactome_logFC_pval_table.csv', selection=selection)
        reactome_ora = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_ORA_reactome_enrichment_table.csv', selection=selection)
        go_gsea["Database"] = "GO"
        go_gsea["enrichment_type"] = "GSEA"
        go_ora["Database"] = "GO" 
        go_ora["enrichment_type"] ="ORA"
        reactome_gsea["Database"] = "reactome" 
        reactome_gsea["enrichment_type"] = "GSEA"
        reactome_ora["Database"] = "reactome" 
        reactome_ora["enrichment_type"] = "ORA"
        return pd.concat([go_gsea, go_ora, reactome_gsea, reactome_ora], ignore_index=True)

    elif database == "All":
        go_df = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_GO_{csv_ending_map[enrichment]}.csv', selection=selection)
        reactome_df = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_reactome_{csv_ending_map[enrichment]}.csv',   
                      selection=selection)
        go_df["Database"] = "GO" 
        go_df["enrichment_type"] = enrichment
        reactome_df["Database"] = "reactome" 
        reactome_df["enrichment_type"] = enrichment
        return pd.concat([go_df, reactome_df], ignore_index=True)

    elif enrichment == "All":
        gsea_df = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_GSEA_{database}_logFC_pval_table.csv', selection=selection)
        ora_df = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_ORA_{database}_enrichment_table.csv', selection=selection)
        gsea_df["Database"] = database
        gsea_df["enrichment_type"] = "GSEA"
        ora_df["Database"] = database
        ora_df["enrichment_type"] = "ORA"
        return pd.concat([gsea_df, ora_df], ignore_index=True)

    else:
        df = check_read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_{database}_{csv_ending_map[enrichment]}.csv', selection=selection)
        df["Database"] = database
        df["enrichment_type"] = enrichment
        return df

def i_to_row_col(i, cols):
    row = i // cols + 1 #returns division
    col = i % cols + 1 #returns remainder
    return row, col

def get_font_size(title_text):
    if len(title_text[0]) > 65: 
        font_size = 19 
    else: 
        font_size = 25
    return font_size

def check_read_csv(path, index_col=None, selection=None):
    try:
        return pd.read_csv(path, index_col=index_col)
    except FileNotFoundError:
        st.error(
            f"""
            No significant results for this comparison.

            Selection:
            - Cell type: {selection['cell_type']}
            - Comparison: {selection['comparison']}
            - Enrichment: {selection['enrichment']}
            - Database: {selection['database']}
            """
        )
        st.stop()

def set_config():
    st.set_page_config(page_title="AD Multiomics Viz", layout="wide")

def sidebar_selectors():
    comparison = st.sidebar.selectbox("Select Comparison", ["ABCA7LOF vs APOE3", "APOE4 vs APOE3"])
    enrichment = st.sidebar.radio("Select Enrichment Type", ["GSEA", "ORA", "All"], index=2)
    database = st.sidebar.radio("Select Database", ["reactome", "GO", "All"], index=2)
    return comparison, enrichment, database

def pathway_filters(enrichment_df, enrichment):
    if enrichment in ("GSEA", "All"):
        nes_min = float(enrichment_df['NES'].min())
        nes_max = float(enrichment_df['NES'].max())
        nes_range = st.sidebar.slider("Normalized Enrichment Score (NES)", nes_min, nes_max, (nes_min, nes_max))
        gsea_rows = enrichment_df[enrichment_df['NES'].between(nes_range[0], nes_range[1])]
        ora_rows = enrichment_df[enrichment_df['NES'].isna()]
        filtered_enrichment_df = pd.concat([gsea_rows, ora_rows], ignore_index=True)
    else:
        filtered_enrichment_df = enrichment_df
        
    minp = filtered_enrichment_df['p.adjust'].min()
    minvalue = int(np.floor(np.log10(minp)))
    maxvalue = int(np.ceil(np.log10(0.05)))
    padj_options = [10 ** e for e in range(minvalue, maxvalue + 1)]
    padj_range = st.sidebar.select_slider("Adjusted p-value", options=padj_options, value=(padj_options[0], padj_options[-1]))
    filtered_enrichment_df = filtered_enrichment_df[filtered_enrichment_df['p.adjust'].between(padj_range[0],padj_range[1])] 
    return filtered_enrichment_df