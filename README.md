# GSE206848-Data-Normalization-and-Subtype-Analysis
GSE206848 Data Normalization and Subtype Analysis
# 📊 Microarray Preprocessing Pipeline: Osteoarthritis, Rheumatoid Arthritis, and Normal Synovium (GSE206848)

This R script automates the essential first steps for analyzing the **GSE206848** microarray dataset, which investigates gene expression differences in synovial fluid cells among **Osteoarthritis (OA)**, **Rheumatoid Arthritis (RA)**, and **Normal** controls.

The pipeline focuses on data quality control (QC), normalization, low-expression filtering, and preliminary visualization to ensure the data is suitable for comparative statistical analysis.

## 🚀 Key Features

* **Automated Data Retrieval:** Fetches expression data and metadata directly from the **GEO database (GSE206848)** using `GEOquery`.
* **Quantile Normalization:** Applies **quantile normalization** (`limma::normalizeBetweenArrays`) to correct technical variation and ensure comparable gene distributions across all samples.
* **Targeted Subtyping:** Extracts and clearly labels the three distinct sample groups: **Normal**, **Osteoarthritis (OA)**, and **Rheumatoid Arthritis (RA)** based on sample titles.
* **QC Visualization:** Generates **boxplots** before and after normalization to visually confirm data distribution consistency.
* **Low-Expression Filtering:** Filters out genes with mean expression in the lowest quartile (25th percentile) to reduce noise.
* **Subtype Clustering:** Performs and visualizes **Principal Component Analysis (PCA)** to check for clear sample separation based on the biological condition.
* **Data Persistence:** Saves the fully processed and filtered expression matrix, metadata, and subtype information to an `.RData` file.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE206848 | Gene expression from synovial fluid cells in joint disease. |
| **Normalization** | Quantile Normalization | Standardizes gene expression distributions across all arrays/samples. |
| **Filtering** | Interquartile Range (IQR) Filtering | Removes non-informative, low-expressed genes for improved statistical power. |
| **Clustering** | PCA | Assesses global sample similarity and confirms separation between Normal, OA, and RA groups. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script loads the following essential packages:
* `GEOquery` (For data retrieval)
* `limma` (For normalization)
* `dplyr` (For data manipulation and subtype extraction)
* `ggplot2` (For PCA visualization)

### ⚙️ Execution

1.  **Download** the `GSE206848 Data Normalization and Subtype Analysis.R` file.
2.  **Execute** the script in your R environment:
    ```R
    source("GSE206848 Data Normalization and Subtype Analysis.R")
    ```
    *Note: All output files are saved to the current working directory where the script is executed.*

---

## 📁 Output Files (3 Plots + 1 Data File)

| Filename | Type | Description |
| :--- | :--- | :--- |
| `GSE206848_processed_data_step1.RData` | R Binary Data | Contains the final, filtered, and normalized `expression_data`, `metadata`, and `subtype` factor for downstream DGE analysis. |
| `GSE206848_boxplot_before_normalization.png` | QC | Boxplot illustrating raw data distributions. |
| `GSE206848_boxplot_after_normalization.png` | QC | Boxplot confirming uniform data distributions across all samples after quantile normalization. |
| `GSE206848_pca_plot.png` | Clustering | **Principal Component Analysis (PCA)** plot, showing sample clustering colored by condition (Normal, OA, RA). |
