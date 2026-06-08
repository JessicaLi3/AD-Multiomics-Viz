# AD-Multiomics-Viz

> **Status: Actively in development** — built as part of a UROP in Fraenkel Lab at MIT, February 2025 - Present.

An interactive web app for exploratory analysis of Alzheimer's disease multi-omics data, including **proteomics**, **phosphoproteomics**, and **acetylomics** datasets.

---
## Background

Alzheimer's disease (AD) involves complex molecular dysregulation across multiple layers of the proteome. This tool is designed to support hypothesis generation and visualization for:

- **Proteomics** — global protein abundance changes in AD vs. control brain tissue
- **Phosphoproteomics** — differential phosphorylation events, relevant to kinase signaling and tau pathology
- **Acetylomics** — lysine acetylation changes tied to metabolic and epigenetic regulation
---
## Methods & Features

| Feature | Description |
|---|---|
| Differential analysis | Log2 fold-change + adjusted p-value across conditions |
| Volcano plots | Interactive visualization of significant hits |
| Gene Set Enrichment Analysis (GSEA) | Pre-ranked GSEA using curated pathway databases (e.g. MSigDB, KEGG) |
| Over-representation analysis (ORA) | Fisher's exact test for pathway enrichment in hit lists |
| Multi-omics integration | Cross-layer comparison of proteins, phosphosites, and acetylsites |

---

## Tech Stack

- **Python 3.x**
- **Streamlit** — interactive web interface
- **pandas / numpy** — data wrangling
- **scipy / statsmodels** — statistical testing
- **plotly** — interactive volcano plots and charts
- **gseapy** — GSEA and ORA
---
## Project context

This app is being developed as part of a UROP (Undergraduate Research Opportunities Program) in Fraenkel Lab. The goal is to build a reproducible, user-friendly interface for collaborators to explore omics data without needing to write code.

**Ongoing work includes:**
- [ ] [Adding in phosphoproteomics and acetylomics data]
- [ ] [Adding in heat maps and other plots]
