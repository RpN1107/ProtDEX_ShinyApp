"""# Protein DEX: Shiny App for Proteomics Differential Expression

The **Protein DEX** app is a Shiny-based interactive platform for analyzing differential expression in proteomics datasets.  
It provides preprocessing, normalization, statistical testing, and visualization tools in one place, with easy export of results and plots.

------------------------------------------------------------
🚀 Features
------------------------------------------------------------
- Upload CSV/Excel proteomics data
- Preprocess: impute missing values, remove contaminants
- Normalize data (Median / Quantile)
- Perform differential expression analysis (FC & p-value)
- Generate plots: PCA, volcano, violin/boxplots, heatmaps
- Download normalized datasets and results

------------------------------------------------------------
🛠️ Installation
------------------------------------------------------------
1. Clone this repository:
   git clone https://github.com/<your-username>/protein-dex.git
   cd protein-dex

2. Install required R packages:
   install.packages(c(
     "shiny", "tidyr", "readxl", "writexl", 
     "preprocessCore", "ggplot2", "DT", 
     "ggrepel", "dplyr", "reshape2", "pheatmap"
   ))

------------------------------------------------------------
▶️ Running the App
------------------------------------------------------------
In R or RStudio, run:

   library(shiny)
   runApp("app.R")   # or the filename where you saved the code

The app will open in your browser.

------------------------------------------------------------
📂 File Structure
------------------------------------------------------------
protein-dex/
│── app.R              # Main Shiny application
│── README.txt         # Project documentation
│── ABOUT.md           # About page
|── TestData/          # Sample Datasets
   | ── test.csv
   | ── test.xlsx

------------------------------------------------------------
🧪 Test Data
------------------------------------------------------------
A small test dataset is provided under the `tests/` folder:

- `tests/test_data.csv`
- `tests/test_data.xlsx`

You can upload these files in the Shiny app to explore all features.

------------------------------------------------------------
📊 Example Workflow
------------------------------------------------------------
1. Upload your proteomics dataset (.csv or .xlsx)
2. Choose label column, group1, and group2 samples
3. Normalize (Median/Quantile)
4. Run analysis → view volcano plot, PCA, and heatmap
5. Export results and plots for publication

------------------------------------------------------------
👨‍💻 Author
------------------------------------------------------------
Developed by Rithwik Nambiar (2025)

------------------------------------------------------------
📜 License
------------------------------------------------------------
This project is licensed under the MIT License – see LICENSE for details.
"""
