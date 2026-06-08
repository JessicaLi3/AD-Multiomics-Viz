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

# Page Config
set_config()

# Load Data
st.sidebar.header("Controls")
cell_type = st.sidebar.selectbox("Select Cell Type", ["astrocytes", "microglia", "pericytes"])
if cell_type == "astrocytes":
    sampcol = 'Samples'
elif cell_type == "microglia":
    sampcol = 'Sample'
elif cell_type == "pericytes":
    sampcol = 'Sample'
    
intensity, samp_meta, entry_meta = load_data(
f'data/proteomics/{cell_type}_prot_data.csv',
f'data/proteomics/{cell_type}_prot_samp_meta.csv',
f'data/proteomics/{cell_type}_prot_meta.csv', 
sample_col = sampcol)

comparison = st.sidebar.selectbox("Select Comparison", ["ABCA7LOF vs APOE3", "APOE4 vs APOE3"])

if comparison == "ABCA7LOF vs APOE3":
    comp_key = "ABCA7LOFvnormal"
elif comparison == "APOE4 vs APOE3":
    comp_key = "APOE4v3"

diff = check_read_csv(f'data/differential_analysis/{cell_type}_proteome/{comp_key}_covDonor_genotype.csv', index_col = 0)

diff['id.mapped'] = diff.index.map(entry_meta['id.mapped'])
diff['entry_name'] = diff.index.map(entry_meta['entry_name'])
diff['log_adjpval'] = diff['adj.P.Val'].apply(lambda x: np.log10(x) * -1)

diff['color'] = "gray"
diff.loc[(diff['log_adjpval'] > 1) & (diff['logFC'] < 0), 'color'] = "blue"
diff.loc[(diff['log_adjpval'] > 1) & (diff['logFC'] > 0), 'color'] = "red"
    
fig = make_subplots(
    rows=1,
    cols=1,
)

scatter = px.scatter(
    diff,
    x='logFC',
    y='log_adjpval',
    color="color",
    color_discrete_map={
        "gray": "lightgray",
        "blue": "blue",
        "red": "red"
    }
)

scatter.update_traces(
    selector=dict(mode='markers'),
    marker=dict(
        size=6,           
        opacity = 0.5,
    ),
    customdata=np.stack([
        diff['id.mapped'],
        diff['entry_name'],
        diff.index,
    ], axis=-1),
    hovertemplate=(
        "Gene Symbol: %{customdata[0]}<br>"
        "Full name: %{customdata[1]}<br>"
        "UniProt: %{customdata[2]}<br>"
        "logFC: %{x:.2f}<br>"
        "adj p-value: %{y:.2f}<extra></extra>"
    )
)

fig.add_traces(scatter.data)

if comp_key == "ABCA7LOFvnormal":
    xaxis_title = "logFC<br><-- higher in APOE3      higher in ABCA7 -->"
elif comp_key == "APOE4v3":
    xaxis_title = "logFC<br><-- higher in APOE3      higher in APOE4 -->"
fig.update_layout(
    xaxis_title=dict(
        text=xaxis_title,
        font=dict(size=27, color="black")),
    yaxis_title=dict(
        text="log10(adj. pval)",
        font=dict(size=27, color="black")),
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, gridcolor='black', showline=True, zeroline=False),
    yaxis=dict(showgrid=False, gridcolor='black', showline=True, zeroline=False),
    height=600,
    width=700,
    showlegend=False
)

fig.update_xaxes(
    showline=True,
    linecolor="black",
    tickcolor="black",
    tickfont=dict(size=24, color="black"),
    zeroline=False
)

fig.update_yaxes(
    showline=True,
    linecolor="black",
    tickcolor="black",
    tickfont=dict(size=24, color="black"),
    zeroline=False
)

fig.add_hline(y=1)

st.plotly_chart(fig, use_container_width=True)