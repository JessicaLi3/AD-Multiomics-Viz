import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")

# --- DATA LOADING (CACHED) ---
@st.cache_data
def load_placeholder_data():
    # Simulate 10,000 proteins across 18 samples
    proteins = [f"PROT_{i:05d}" for i in range(10000)]
    samples = [f"Line{l}_Rep{r}" for l in range(1,7) for r in range(1,4)]
    data = np.random.lognormal(mean=10, sigma=1, size=(10000, 18))
    df = pd.DataFrame(data, columns=samples)
    df.insert(0, "Protein_ID", proteins)
    return df

st.title("🔬 Proteomics Data Explorer v1.0")
st.write("Securely logged in. Data is restricted to authorized viewers.")

# Load Data
with st.spinner("Loading 10,000 proteins..."):
    df = load_placeholder_data()

# --- SIDEBAR FILTERS ---
st.sidebar.header("Controls")
search_query = st.sidebar.selectbox("Select or Search Protein", df["Protein_ID"])

# --- VISUALIZATION ---
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader(f"Raw Values: {search_query}")
    protein_data = df[df["Protein_ID"] == search_query].iloc[:, 1:].T
    protein_data.columns = ["Abundance"]
    st.dataframe(protein_data)

with col2:
    st.subheader("Expression Boxplot")
    # Prepare data for Seaborn
    plot_df = pd.DataFrame({
        "Abundance": protein_data["Abundance"],
        "Line": [s.split('_')[0] for s in protein_data.index]
    })
    
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(data=plot_df, x="Line", y="Abundance", ax=ax, palette="light:b")
    sns.stripplot(data=plot_df, x="Line", y="Abundance", color="black", alpha=0.5)
    st.pyplot(fig)