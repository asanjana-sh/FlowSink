# An optimal transport-based low-dimensional visualization framework for high-parameter flow cytometry

## Overview
This paper provides a comprehensive framework for analyzing and visualizing high-parameter flow cytometry data, and demonstrates the methodology using a **Malignant Peritoneal Mesothelioma (MPM) immunotherapy dataset**. The main dataset consists of PBMC samples collected from 14 patients at three time points: baseline (before vaccination), after the first vaccine, and after the third vaccine. The dataset was collected under six panels, of which three T-cell panels (Co-inhibition, Co-stimulation, Cytokine) are the primary focus of analysis in this paper.

The paper implements a **low-dimensional visualization framework based on optimal transport (Sinkhorn distance) and Graph Edit Distance (GED)**, enabling comparative analysis of samples over time and correlation with patient outcomes (OS, PFS, and response to therapy).  

---

## Data Details

The data files are available [here](https://drive.google.com/file/d/1YHeYZlL0vsViHzcQzffB1q7HxTgDtr-0/view?usp=sharing). Downloading the zip folder and extracting it under the `FlowSink/` directory should work with the relative paths in the code files.

### 1. Raw Data (`FlowSink/Data/Raw Data`)
- The primary data in `.fcs` format are provided by the original study authors. 
- The `.fcs` files are preprocessed in R using Bioconductor packages (`PeacoQC`, `flowCore`, `CytoNorm2.0`) for margin event removal, compensation, transformation, scaling, and quality control. `.fcs` files are stored after every step of the preprocessing pipeline.
- Important metadata files (`_cell_type.xlsx`) are stored in each panel subfolder under `Normalized Data`:
  - **Channel:** Marker-fluorochrome mapping.
  - **Phenotype:** Targeted cell population definitions for semi-automated cell assignment.
  - **Marker_threshold:** Expert-curated gating thresholds.
  - **Test_edge_filter:** Expert-curated adjacency matrix for graph visualization.
- The final cell type assigned data are in `.xlsx` format. These files are used in the code.
- **File naming convention:**  
  `NormalizationStep_PanelName_GroupID_TubeID_PatientCode_TimePoint_QC_ParentPopulation_CellClassification.xlsx`  
  Example: `Norm_Norm_coinhib_g1_Tube_001_MCV001_VAC_1_QC_Lymph_CellType.xlsx`
- Common **cell type markers**: `CD56, CD3, CD4, FoxP3, CD8, CD45RA, CCR7, Live/Dead`  
- Panel-specific **cell state markers**:
  - Co-inhibition: `LAG3, PD1, TIM3, CD39, KI67, CTLA-4`
  - Co-stimulation: `CD28, CD137, PD1, HLA-DR, ICOS, KI67`
  - Cytokine: `PD1, TBET, IL-10, TNF-a, IL2, IFN-y`
- Cell type assignment is semi-automated, guided by expert manual gating. Outlier cells labeled as `"Other cell"` are excluded during experiments.

### 2. Data for Downstream Tasks (`FlowSink/Data/Malignant Peritoneal Mesothelioma Final Dataset`)
- Files are stored as `.xlsx` in separate folders for the three T-cell panels (`Co-inhibition`, `Co-stimulation`, `Cytokine`) and for three time points (`TP1`, `TP2`, `TP3`).

### 3. Metadata
- `metadata.xlsx` contains the mapping of original and renamed files, and patient clinical outcomes (OS, PFS, response).
- OS and PFS values are provided in months.

---

## Code Organization

### 1. Preprocessing 
- Preprocessing pipeline for original `.fcs` files:
  1. Rename raw `.fcs` files and organize by time point. (`RenameFCSFiles.R`)
  2. Margin event removal, compensation, transformation, scaling, QC. (`PeacoQC.R`, `MPM…Panel.R`)
  3. Lymphocyte gating and merging of gated `.fcs` files. (`MPMGate1Merge.R`)
  4. Two-step normalization (cell type markers, then cell state markers). (`MPMTcellNormalize.R`, `MPMTcellNormalize_2.R`)
  5. Cell type assignment (semi-automated). (`MPM…CellAssign.R`)  

*Note:* The preprocessed `.xlsx` files provided can be used directly; re-running preprocessing is generally unnecessary.

### 2. Main Analysis
The `MPM...ComparePlot.R` scripts can be used as the main demonstration code.

- **Panel-specific workflows:** Separate folders for `CoInhibition`, `CoStimulation`, `Cytokine`.
- **Core utilities:** `MPMUtilities.R` (shared computational functions), `MPMMarkerChannelMap.R` (marker-fluorochrome mapping), and result summarization scripts (`MPMDimRedExample.R`, `MPMTcellGEDDistribution.R`).
- **Analysis workflow:**
  1. Compute Sinkhorn distance matrices for each sample (Python via R interface: `TestOTSinkhorn.py`).
  2. Construct filtered phenotype graphs and compute GED between samples.
  3. Visualize single-sample graphs and comparative graphs across time points.
  4. Generate low-dimensional MDS visualizations and correlate GED with clinical outcomes (OS, PFS, response).

*Note:* Python environment setup is required for Sinkhorn and GED computations. Python scripts `TestOTSinkhorn.py` and `TestGED.py` are called from R.

---

## Key Outputs
- **Filtered Graphs:** Single-sample and comparative graphs per panel and time point.
- **GED Matrices:** Pairwise Graph Edit Distances between samples across time points.
- **MDS Plots:** Low-dimensional embeddings of samples based on GED.
- **Survival Analysis:** Visualizations correlating GED with patient response, OS, and PFS.

---

## References
1. Shemonti AS, Gmyrek GB, Quintelier KL, Van Gassen S, Saeys Y, Willemsen M, Aerts JG, Madsen EV, Robinson JP, Pothen A, Rajwa B. *[An optimal transport-based low-dimensional visualization framework for high-parameter flow cytometry](https://www.biorxiv.org/content/10.1101/2025.08.19.670604v1.full.pdf)*. bioRxiv. 2025:2025-08. 
2. Quintelier KL, Willemsen M, Bosteels V, Aerts JG, Saeys Y, Van Gassen S. *CytoNorm 2.0: A flexible normalization framework for cytometry data without requiring dedicated controls*. Cytometry Part A. 2025 Feb;107(2):69-87.
3. Dietz MV, Quintelier KL, van Kooten JP, de Boer NL, Vink M, Brandt-Kerkhof AR, Verhoef C, Saeys Y, Aerts JG, Willemsen M, Van Gassen S. *Adjuvant dendritic cell-based immunotherapy after cytoreductive surgery and hyperthermic intraperitoneal chemotherapy in patients with malignant peritoneal mesothelioma: a phase II clinical trial*. Journal for Immunotherapy of Cancer. 2023 Aug 3;11(8):e007070.

---

## Notes
- The pipeline is modular: preprocessing is separate from analysis, and panel-specific differences are handled in dedicated folders.
- Paths in code may be system-dependent; update carefully before running.
- The dataset and code can serve as a reference for analyzing high-dimensional flow cytometry data using optimal transport and graph-based methods.

Correspondence: Abida Sanjana Shemonti (asanjana@miftek.com)
