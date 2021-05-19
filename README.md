

SPmarker
===
SPmarker is a machine learning based approach for identification of marker genes and classification of cells in plant tissues

# Introduction
An essential facet of the single-cell RNA sequencing (scRNA-seq) studies is to use marker genes that distinguish heterogeneous cell populations and dissect biological functions of each cell. In plants, the scRNA-seq focused mostly on the understood Arabidopsis root system. However, few suitable computational methodologies aid the identification of the novel marker genes in the root system. Here, we introduce SPmarker, a machine learning based method to identify the marker genes via identifying their feature importance from Random Forests (RF) model using the SHapley Additive exPlanations (SHAP) package (SHAP markers).

# Dependence and requirements
SPmarker is developed in Python with modules and external tools.

Before running this pipeline, a dependency check should be performed first to make sure every dependency is correctly installed.

For information about installing the dependencies, please see below. The version numbers listed below represents the version this pipeline is developed with, and using the newest version is recommended.

## Requirements
**Python** (v3.0 or more; Developed and tested with version 3.7.1)  
**pandas** (python package; Developed and tested with version 0.24.2)  
**numpy** (python package; Developed and tested with version 1.16.4)  
**itertools** (python package; Developed and tested with version 2.3)  
**sklearn** (python package; Developed and tested with version 0.22.2)  
**shap** (python package; Developed and tested with version 0.35.0)  
#**Seurat** (R package; Developed and tested with version 3.0)

# Quick Start
## Installation
1.	Download SPmarker  
> git clone https://github.com/LiLabAtVT/SPMarker.git
2.	Use scripts wrapped in the SPmarker to run pipeline

## Required input data
**1. (mandatory) Gene expression matrix file (.csv)**  
The matrix row names are features/genes and column names are cell names  

|  | cell1 | cell2 | cell3 |
| -------- | ----- | ----- | ----- |
| gene1    | 2   | 1     | 2
| gene2    | 3     | 5     |3
| gene3    | 3     | 5     |0

**2. (optional) cell meta file (.csv)**
This file contains three columns. prob means probability that the cell will be assigned to the cell type. this prob value is obtained from ICI method. If it is not from the ICI method or other methods that could give a probability, we can use 1 to represent prop.  
|  | cell_type | prob |
| -------- | -------- | -------- |
| cell1    | celltype1     | 1     |
| cell2    | celltype1     | 0.8     |
| cell3     | celltype2     | 0.6     |

**3. (optional) unknown cell matrix (.csv)**  
UN means unknown
|  | UNcell1 | UNcell2 | UNcell3 |
| -------- | ----- | ----- | ----- |
| gene1    | 2   | 1     | 2
| gene2    | 3     | 5     |3
| gene3    | 3     | 5     |0

**Note**: Users must provide **optional 2** or **optional 3**.

## Example Run 1
### Input data
**1. Gene expression matrix file**  
(ipt_test1_exp.csv)
**2. cell meta file**  
(ipt_test1_meta.csv)
**3. unknown cell matrix**  
(ipt_test1_unknown_exp.csv)

### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_exp.csv \\  
-meta ipt_test1_meta.csv \\  
-ukn_mtx ipt_test1_unknown_exp.csv  



## Example Run 2
### Input data
**1. Gene expression matrix file**  
(ipt_test2_exp.csv)
**2. GFP marker gene name**  
(eg. GFP_marker, the name 'GFP_marker' should be one of feature names in the ipt_test2_exp.csv)

### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test2_exp.csv \\  
-m GFP_marker \\  
-ukn_mtx ipt_test2_unknown_exp.csv  

## Outputs
**1. opt_SHAP_markers_dir and opt_SVM_markers_dir**  
a. opt_top_20_novel_known_marker.txt  
b. opt_top_20_novel_marker.txt  
c. opt_top_20_summary_marker_composition.txt  
**Note:**  
'**top_20**' means users define they want to select the top 20 markers based on the feature importance values.  
'**novel**' means the predicted markers are novel markers. If users do not provide known marker file using '-kmar_fl', all the markers will be labeled as novel.  
'**novel_known**' means the output contains the novel and known markers at the same time.  
'**marker_composition**' means among the 20 markers, how many of them are novel and how many of them are known markers.  

**2.opt_prediction_dir**
a. opt_prediction_RFindetest.txt  
b. opt_prediction_SVMindetest.txt  
**Note:**  
'**RF**' means the prediction of cell types on unknown cells are based on Random Forest model.  
'**SVM**' means the prediction of cell types on unknown cells are based on Support Vector Machine model.  



## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::


