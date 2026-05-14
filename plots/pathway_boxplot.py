import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import numpy as np
import textwrap as tw
import tarfile
import pdb
import math
import sys
import textwrap as tw
sys.path.append('../')
from ad_viz_utility import *

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")

# # --- SIDEBAR FILTERS --- 
st.sidebar.header("Controls")

cell_type = st.sidebar.selectbox("Select Cell Type", ["astrocytes", "microglia", "pericytes"])
if cell_type == "astrocytes":
    sampcol = 'Samples'
elif cell_type == "microglia":
    sampcol = 'Sample'
elif cell_type == "pericytes":
    sampcol = 'Sample'
with st.spinner("Loading 10,000 proteins..."):
    intensity, samp_meta, entry_meta = load_data(
        f'data/proteomics/{cell_type}_prot_data.csv',
        f'data/proteomics/{cell_type}_prot_samp_meta.csv',
        f'data/proteomics/{cell_type}_prot_meta.csv', 
        sample_col = sampcol)

comparison = st.sidebar.selectbox("Select Comparison", ["ABCA7LOF vs APOE3", "APOE4 vs APOE3"])
enrichment = st.sidebar.selectbox("Select Enrichment Type", ["GSEA", "ORA"])
database = st.sidebar.selectbox("Select Database", ["reactome", "GO"])
 
comp_key, csv_ending, pathway_identifier, gene_set, protein_identifier = file_keys(comparison, enrichment, database, cell_type)

diff = pd.read_csv(f'data/differential_analysis/{cell_type}_proteome/{comp_key}_covDonor_genotype.csv', index_col = 0)
enrichment_df = pd.read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_{database}_{csv_ending}.csv')

#sorting by significance 
diff['log_adjpval'] = diff['adj.P.Val'].apply(lambda x: np.log10(x) * -1)
diff['logFC_logpval'] = np.abs(diff['logFC']*diff['log_adjpval']) # rank proteins by this

# selecting a pathway
pathway_list = enrichment_df['Description']
pathway_list = sorted(pathway_list)
pathway_name = st.sidebar.selectbox("Select or Search Pathways", pathway_list)

# Setting up pages
id_list = pathway_to_uniprot(gene_set, database, pathway_identifier, pathway_name, protein_identifier, entry_meta)
id_list = [p for p in id_list if p in entry_meta.index]
ranked_id_list = (diff.loc[id_list].sort_values(by='logFC_logpval', ascending=False).index.tolist())

plots_per_page = 9

if 'current_page' not in st.session_state:
    st.session_state.current_page = 0

total_pages = math.ceil(len(ranked_id_list) / plots_per_page)
start_idx = st.session_state.current_page * plots_per_page
end_idx = start_idx + plots_per_page
current_proteins = ranked_id_list[start_idx:end_idx]
valid_proteins = [p for p in current_proteins if p in intensity.columns]

st.write(f"### Showing proteins {start_idx + 1} to {min(end_idx, len(ranked_id_list))} "
        f"out of {len(ranked_id_list)} proteins (Page {st.session_state.current_page + 1} of {total_pages})")

with st.spinner("Loading pathway boxplots..."):
    
    pathway_info = enrichment_df[enrichment_df['Description'] == pathway_name].copy()
    if 'setSize' not in pathway_info.columns:
        pathway_info['setSize'] = np.nan
    if 'NES' not in pathway_info.columns:
        pathway_info['NES'] = np.nan

    summary_table = pathway_info[[
        'Description',
        'p.adjust',
        'setSize',
        'NES']]
    st.dataframe(summary_table)

subplot_titles = [
    "<br>".join(
        tw.wrap(
            f"{lookup_by_protid(p, entry_meta, 'id.mapped')} "
            f"({lookup_by_protid(p, entry_meta, 'id.description')})<br>"
            f"{lookup_by_protid(p, entry_meta, 'entry_name')}",
            width=45
        )
    )
    for p in current_proteins
]

# --- PLOTS ---
fig = make_subplots(
    rows=3,
    cols=3,
    subplot_titles=subplot_titles
)

for i, prot_id in enumerate(valid_proteins):
# boxplots
    plot_df = intensity[[prot_id]].copy()
    plot_df.columns = ["value"]
    plot_df = plot_df.merge(samp_meta[['Parental.line', 'CRISPR.edits', sampcol]], left_index = True, right_index = True)
    row, col = i_to_row_col(i, 3)
    
    subfig = px.box(
        plot_df,
        x="Parental.line",
        y="value",
        color="CRISPR.edits",
    )
        
    subfig.update_traces(
        opacity=0.3,
        hoverinfo="skip", # disables hover
        hovertemplate=None, # disables hover
        boxpoints=False,
        selector=dict(type="box")
    )

    # draw stripplot points on top (disable hover)
    # use this instead of boxpoints="all" in the boxplot
    # because it draws the points directly on top of the box
    strip = px.strip(plot_df, x="Parental.line", y="value", color="CRISPR.edits")
    strip.update_traces(
        marker_size=5,
        hoverinfo="skip", # disables hover
        hovertemplate=None, # disables hover
        showlegend=False,
    )

    # draw an invisible barplot
    # the invisible bars go all the way to the top of the plot ("y_max")
    medians = plot_df.groupby(["Parental.line", "CRISPR.edits"],observed=True)["value"].median().reset_index()
    medians["y_max"] = plot_df["value"].max()
    medians["y_min"] = plot_df["value"].min()
    medians["y_span"] = medians["y_max"] - medians["y_min"]
    bar = px.bar(medians, x="Parental.line", y="y_span", base = "y_min", color="CRISPR.edits", custom_data=["value"])
    bar.update_traces(
        marker_color = "rgba(0,0,0,0)",
        marker_line_color = "rgba(0,0,0,0)",
        marker_line_width = 0,
        showlegend=False,
        hovertemplate = "Parental.line: %{x}<br>"
        "CRISPR.edits: %{fullData.name}<br>"
        "Median: %{customdata[0]:.2f}"
        "<extra></extra>"
    )
    
    # adds the title
    prot_name = lookup_by_protid(prot_id, entry_meta, 'id.mapped')
    prot_desc = lookup_by_protid(prot_id, entry_meta, 'entry_name')
    prot_id = lookup_by_protid(prot_id, entry_meta, 'id.description')

    # add the stripplot and barplot to the original figure
    subfig.add_traces(strip.data)
    subfig.add_traces(bar.data)
    fig.add_traces(subfig.data, rows=row, cols=col)

    fig.update_yaxes(
        title_text="Normalized levels",
        showgrid=False,
        tickfont=dict(size=11),
        zeroline=False,
        showline=True,
        linecolor="black",
        linewidth=1,
    )

    fig.update_xaxes(
        title_text="",
        tickfont=dict(size=11),
        showgrid=False,
        zeroline = False,
        showline = True,
        linecolor="black",
        linewidth=1,
    )

    names = set()
    for trace in fig.data:
        # If the legend name is already in the set, hide it
        if trace.name in names:
            trace.showlegend = False
        # Otherwise, show it and add it to the set
        else:
            trace.showlegend = True
            names.add(trace.name)
    
    fig.update_layout(
        height=900,
        width=1200,
        boxmode = 'group'
    )

    fig.update_xaxes(title_font=dict(color="black"))
    fig.update_yaxes(title_font=dict(size=16, color="black"))
    fig.update_xaxes(tickfont=dict(size=16, color="black"))
    fig.update_yaxes(tickfont=dict(size=16, color="black"))
    fig.update_annotations(font=dict(size=16, color="black"))

st.plotly_chart(fig, use_container_width=True)

# Page controls
col1, col2, col3 = st.columns([1,2,1])

with col1:
    if st.button("Previous") and st.session_state.current_page > 0:
        st.session_state.current_page -= 1
        st.rerun()

with col2:
    target_page = st.number_input(
        "Jump to page",
        min_value=1,
        max_value=total_pages,
        value=st.session_state.current_page + 1,
        step=1,
        label_visibility="collapsed"
    )
    if target_page-1 != st.session_state.current_page:
        st.session_state.current_page = target_page-1
        st.rerun()

with col3:
    if st.button("Next") and st.session_state.current_page < total_pages-1:
        st.session_state.current_page += 1
        st.rerun()
