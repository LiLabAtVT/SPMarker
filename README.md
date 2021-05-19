

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

## Example Run 1
## Input data
**1. Gene expression matrix file (.csv)**  
**(ipt_test1_exp.csv)**
The matrix row names are features/genes and column names are cell names  
|  | cell1 | cell2 | cell3 |
| -------- | ----- | ----- | ----- |
| gene1    | 2   | 1     | 2
| gene2    | 3     | 5     |3
| gene3    | 3     | 5     |0

**2. cell meta file (.csv)**  
**(ipt_test1_meta.csv)**
This file contains three columns. prob means probability that the cell will be assigned to the cell type. this prob value is obtained from ICI method. If it is not from the ICI method or other methods that could give a probability, we can use 1 to represent prop.  
|  | cell_type | prob |
| -------- | -------- | -------- |
| cell1    | celltype1     | 1     |
| cell2    | celltype1     | 0.8     |
| cell3     | celltype2     | 0.6     |

**3. unknown cell matrix (.csv)**  
**(ipt_test1_unknown_exp.csv)**
UN means unknown
|  | UNcell1 | UNcell2 | UNcell3 |
| -------- | ----- | ----- | ----- |
| gene1    | 2   | 1     | 2
| gene2    | 3     | 5     |3
| gene3    | 3     | 5     |0

## Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_exp.csv \\  
-meta ipt_test1_meta.csv \\  
-ukn_mtx ipt_test1_unknown_exp.csv 

**Outputs:**  

## Example Run 2
## Input data





## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::


