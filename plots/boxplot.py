import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import textwrap as tw
sys.path.append('../')
from ad_viz_utility import *

max_msg_size = st.get_option("server.maxMessageSize")

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")
# add a comment

st.title("🔬 Proteomics Data Explorer v0.1")

# Load Data
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

# --- SIDEBAR FILTERS ---
# Exclude metadata columns so they don't show up in the protein dropdown
meta_cols = ['Parental.line', 'CRISPR.edits', 'Samples']
entry_meta['uniprot_id.mapped_concat'] = entry_meta['id.mapped'] + '_' + entry_meta.index
protein_list = entry_meta['uniprot_id.mapped_concat']
protein_list = sorted(protein_list)

st.sidebar.header("Controls")
id_name = st.sidebar.selectbox("Select or Search Protein", protein_list)
prot_id = entry_meta[entry_meta['uniprot_id.mapped_concat'] == id_name].index[0]

# --- VISUALIZATION ---
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader(f"Raw Values: {prot_id}")
    protein_data = intensity[[prot_id]].copy()
    protein_data.columns = ["Normalized Abundance"]
    st.dataframe(protein_data)

with col2:
    st.subheader("Expression Boxplot")
    if prot_id in intensity.index:
        plot_df = intensity.loc[[prot_id]].T
    elif prot_id in intensity.columns:
        plot_df = intensity[[prot_id]]

    plot_df = plot_df.rename(columns={prot_id: "value"})

    plot_df = plot_df.merge(
        samp_meta[['Parental.line', 'CRISPR.edits']],
        left_index=True,
        right_index=True
    )

    plot_df = plot_df.reset_index().rename(columns={"index": "sample"})
    medians = plot_df.groupby(["Parental.line", "CRISPR.edits"], observed=True)["value"].median().reset_index()

    # adds the title
    prot_name = lookup_by_protid(prot_id, entry_meta, 'id.mapped')
    prot_desc = lookup_by_protid(prot_id, entry_meta, 'entry_name')
    prot_id = lookup_by_protid(prot_id, entry_meta, 'id.description')

    title_text = [
        "<br>".join(
            tw.wrap(
                f"{prot_name} ({prot_id})<br>{prot_desc}",
                width=60
            )
        )
    ]

    font_size = get_font_size(title_text)
 
    # --- PLOTLY ---
    fig = make_subplots(
        rows=1,
        cols=1,
        subplot_titles=title_text,
        vertical_spacing=1)

    subfig = px.box(
        plot_df,
        x="Parental.line",
        y="value",
        color="CRISPR.edits",
    )

    # draw stripplot points on top (disable hover)
    # use this instead of boxpoints="all" in the boxplot
    # because it draws the points directly on top of the box
    subfig.update_traces(
        opacity=0.3,
        hoverinfo="skip",
        hovertemplate=None,
        boxpoints=False,
        selector=dict(type="box")
    )

    strip = px.strip(
        plot_df,
        x="Parental.line",
        y="value",
        color="CRISPR.edits"
    )

    strip.update_traces(
        marker_size=5,
        hoverinfo="skip",
        hovertemplate=None,
        showlegend=False,
    )

    # draw an invisible barplot
    # the invisible bars go all the way to the top of the plot ("y_max")
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

    # add the stripplot and barplot to the figure
    subfig.add_traces(strip.data)
    subfig.add_traces(bar.data)
    fig.add_traces(subfig.data, rows=1, cols=1)

    fig.update_yaxes(
        title_text="Normalized levels",
        title_font=dict(size=18, color="black"),
        tickfont=dict(size=16, color="black"),
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="black",
        linewidth=1,
        ticks="outside",
        tickcolor="black"
    )

    fig.update_xaxes(
        title_text="",
        tickfont=dict(size=18, color="black"),
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="black",
        linewidth=1,
        ticks="outside",
        tickcolor="black"
    )
    
    fig.update_annotations(
        font=dict(size=font_size, color="black"),
        align="center"
    )

    fig.update_layout( 
        height=500,
        width=800,
        hovermode="closest",
        plot_bgcolor="white",
        boxmode="group"
    )
    
    # fig.update_layout(
    #     legend=dict(
    #         x=1,
    #         y=0.98,
    #         xanchor="left",
    #         yanchor="top"
    #     )
    # )
    
    # fig.update_layout(
    #     margin=dict(l=80, r=150, t=80, b=80)
    # )

    config = {
        'toImageButtonOptions': {
        'format': 'svg', 
        'filename': 'newplot',
        'height': 500,
        'width': 1100,
        'scale': 1 # can increase for higher resolution
      }
    }

    st.plotly_chart(fig, config=config, use_container_width=True)