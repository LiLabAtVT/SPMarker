

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
**1. Gene expression matrix file**  
(ipt_test1_exp.csv)  
**2. cell meta file**  
(ipt_test1_meta.csv)  
**3. unknown cell matrix**  
(ipt_test1_unknown_exp.csv)  
**4. selected_features.csv**  
(ipt_test1_selected_features.csv)  

### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_exp.csv \\  
-meta ipt_test1_meta.csv \\  
-ukn_mtx ipt_test1_unknown_exp.csv \\  
-feat_fl ipt_test1_selected_features.csv



## Example Run 3
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
d. opt_all_novel_known_marker.txt  
e. opt_all_novel_marker.txt  
f. opt_all_summary_marker_composition.txt  
**Note:**  
'**top_20**' means users define they want to select the top 20 markers based on the feature importance values.  
'**novel**' means the predicted markers are novel markers. If users do not provide known marker file using '-kmar_fl', all the markers will be labeled as novel.  
'**novel_known**' means the output contains the novel and known markers at the same time.  
'**marker_composition**' means among the 20 markers, how many of them are novel and how many of them are known markers.  
'**all**' means the results report all markers intead of top markers.

**2.opt_prediction_dir**  
a. opt_prediction_RFindetest.txt  
b. opt_prediction_SVMindetest.txt  
**Note:**  
'**RF**' means the prediction of cell types on unknown cells are based on Random Forest model.  
'**SVM**' means the prediction of cell types on unknown cells are based on Support Vector Machine model.  

# Usage
```
usage:
**SPmarker**  
SPmarker.py [-h] required: [-d working_dir][-o output_dir]
                           [-mtx expression_matrix_file]
                           [-ukn_mtx unknown_expression_matrix_file]
                           ([-m marker_name]|[-meta meta_file])
                 optional: [-bns no][-bns_ratio 1:1][-cv_num 5]
                           [-indep_ratio 0.1][-eval_score MCC]
                           [-mar_num 20][-kmar_fl known_marker_file]
                           [-SVM yes][-feat_fl feature_file]

arguments:
-h, --help        Show this help message and exit.

-d                Working directory to store intermediate files of each step. 
                  Default: ./ .

-o                Output directory to store the output files. 
                  Default: ./ .

-mtx              Training expression matrix file.
                  Rowname is gene, and column name is cell barcode.
                  Please make sure the gene name do not contain space. Otherwise, the gene name will be transfered to a name with "_" connected.

-ukn_mtx          An expression matrix file with cells that need to be annotated. 
                  Please keep same format as 'Expression matrix file'.
                  If the genes have different orders as training matrix, the SPmarker will automatically keep the same feature orders between the unknown matrix and training matrix.
                  If there are some genes missing in the unknown matrix compared to the training matrix, this tool will assign '0' across all cells for this gene.

-m                Provide a marker name such as a name form an internal GFP marker. 
                  This marker will assign the GFP-related cell identity to a cell where reads can map to this GFP marker.
                  If users provide '-m', '-meta' cannot be provided.

-bns              Balance the matrix of cell identities. This option works only when -m is initiated.

-bns_ratio        Provide a ratio of cell number from different identities.
                  For example, if users set the ratio to be 1:1, they must initiate the '-m' to assign a GFP-related cell identity to cells (eg. 500) where reads can map to this GFP marker. 
                  If there are 2000 cells not be assigned, SPmarker will sample 500 cells from these 1000 cells to allow the ratio to be 1:1. 
                  If users set ratio to be 1:2, SPmarker will sample 1000 cells.
                  Default: 1:1.

-meta             Provide a meta that contains known cell identity for all cells in the training matrix.
                  If the '-meta' is initiated, we should not provide '-m'.

-cv_num           Initiate x fold cross validation.
                  Default: 5.

-indep_ratio      Provide ratio of cells from independent dataset to all cells.
                  For example. If users provide 1000 cells in the '-mtx', SPmarker will sample 100 cells to be independent dataset (default).
                  Default: 0.1.

-eval_score       Provide a type of evaluation score to decide the best model that will be used for marker identification.
                  Default: MCC.

-mar_num          Provide the number of top candidate marker users want to extract as output markers from each cell type.
                  Default: 20.
                  If the feature number is below 20, we will extract all the features under the cell type.

-kmar_fl          Provide the known marker gene list file. 
                  Once users provide this file, they will obtain a file that only contains novel marker genes.
                  
-SVM              Decide to generate the SVM markers.
                  Default: -SVM yes.

-feat_fl          Provide the features that will be kept in the expression file that is used for the training.
                  If users do not provide the argument, we will use all the features from the training matrix (-mtx).


                      



## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::


