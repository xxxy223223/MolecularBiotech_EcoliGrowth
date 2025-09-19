# MolecularBiotech_EcoliGrowth
Python-based analysis of E. coli microplate growth data across multiple temperatures. Course project for KETN20 Molecular Biotechnology, Lund University.

**Python-based analysis of *E. coli* microplate growth data across multiple temperatures.**  
Course project for **KETN20 Molecular Biotechnology**, Lund University.  

## 📌 Project Overview
This project investigates how temperature (27–45 °C) affects the growth of *Escherichia coli* K12 NCM3722 in minimal medium with glucose.  
Growth was monitored in microplates and measured as optical density (OD) every 10 minutes.  


## 📊 Dataset
- Source: Katipoglu-Yazan et al., 2023  
- DOI: [10.57745/GCKG7W](https://doi.org/10.57745/GCKG7W)  
- Replicates: 40 per temperature, 48-well plates  
- Measurements: Optical density (OD), every 10 minutes


## 🖥️ Methods & Tools
- **Python** (Jupyter Notebook)  
- Data cleaning & processing with **pandas**  
- Visualization with **matplotlib**  
- Estimation of growth rate (μ) from OD measurements  
- Comparison of growth at different temperatures to identify the **optimal growth temperature**
- Experimental results showed that the maximum growth rate (μ_max) was significantly higher at **36 °C** compared to 37 °C, indicating that under the given conditions, 36 °C may represent the optimal growth temperature for *E. coli* K12 NCM3722.


## 📂 Repository Structure
MolecularBiotech_EcoliGrowth/  
├── Results/               # Processed results and filtered datasets  
├── ecoli_growth_data/     # Raw experimental data  
├── EcoliGrowth.ipynb      # Main Jupyter Notebook for analysis (Python + visualization)  
├── EcoliGrowth.py         # Python script version of the analysis (developed in Spyder)  
├── README.md              # Project description and documentation  
├── LICENSE                # License information (MIT)  
└── .gitignore             # Ignore unnecessary files (Python cache, checkpoints, etc.)


Project completed by **Chengxi**  
as part of the **KETN20 Molecular Biotechnology** course at **Lund University**.
