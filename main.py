import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw 

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="centered")

pages = {
    "Protein-Level Plots": [
        st.Page("plots/boxplot.py", title="Expression Boxplot"), st.Page("plots/proteins.py", title="All Proteins")
    ],
    "Pathway-Level Plots": [
        st.Page("plots/volcano.py", title="Pathway Volcano Plots"), st.Page("plots/pathway_boxplot.py", title="Pathway Boxplots"), st.Page("plots/scatterplot.py",          title = "Pathway Scatterplots")
    ],
}
pg = st.navigation(pages, position = "top") 
pg.run()

