
SPmarker
===
SPmarker is a machine learning based approach for identification of marker genes and classification of cells in plant tissues

Current release: 07/09/23 v0.3  
Current release: 10/17/21 v0.2  

If you use SPmarker in your own study, please consider citing the following article:  

Haidong Yan,Jiyoung Lee,Qi Song,Qi Li,John Schiefelbein,Bingyu Zhao,Song Li 2022 Identification of new marker genes from plant single-cell RNA-seq data using interpretable machine learning methods. New phytologist. https://doi.org/10.1111/nph.18053.

# Introduction
In order to dissect the biological functions of each individual cells, an essential step in the analysis of single-cell RNA sequencing data is to classify specific cell types with marker genes. In this study, we have developed a machine learning pipeline called Single cell Predictive markers (SPmarker) to assign cell types and to identify novel cell-type marker genes in the Arabidopsis root. Our method can (1) assign cell types based on cells that were labeled using published methods, (2) project cell types identified by trajectory analysis from one dataset to other datasets, and (3) assign cell types based on internal GFP markers. Using SPmarker, we have identified hundreds of new marker genes and majority of these machine learning-derived marker genes were not identified before. As compared to known marker genes, we have found more orthologous genes of these new marker genes in corresponding rice single cell clusters. We have also found 172 new marker genes for Trichoblast in five non-Arabidopsis species, which expands the number of marker genes for this cell type by 35-154%. Our results represent a new approach to identify cell-type marker genes from scRNA-seq data and pave the way for cross-species mapping of scRNA-seq data in plants. 

![](https://i.imgur.com/tUBBGIg.png)





# Dependence and requirements
SPmarker is developed in Python with modules and external tools.

Before running this pipeline, a dependency check should be performed first to make sure every dependency is correctly installed.

For information about installing the dependencies, please see below. The version numbers listed below represents the version this pipeline is developed with, and using the newest version is recommended.

## Requirements
**Python** (v3.0 or more; Developed and tested with version 3.7.1)  
**pandas** (python package; Developed and tested with version 1.1.1)  
**sklearn** (python package; Developed and tested with version 0.24.2)  
**shap** (python package; Developed and tested with version 0.39.0)  
**keras** (python package; Developed and tested with version 2.4.3)  

## Use conda to install required packages (Recommend)  
conda create -n py37 python=3.7  
conda activate py37  
conda install pandas  
conda install scikit-learn  
conda install -c conda-forge shap  
conda install keras  




# Quick Start
## Installation
1.	Download SPmarker  
> git clone https://github.com/LiLabAtVT/SPMarker.git
2.	Use scripts wrapped in the SPmarker to run pipeline

## Required input data
**1. (mandatory) Gene expression matrix file (.csv)**  
Note: The matrix row names are features/genes and column names are cell names. The value here shows normalized read count of each cell mapping to the genes. It is recommended for users to use ‘SCTransform’ in Seurat tool (https://satijalab.org/seurat/articles/sctransform_vignette.html) to normalize the gene expression for each cell.   

|  | cell1 | cell2 | cell3 |
| -------- | ----- | ----- | ----- |
| gene1    | 2.1   | 1.2     | 2.7
| gene2    | 3.3     | 5.7     |3.2
| gene3    | 3.6     | 5.2     |0  


**2. (mandatory) Provide cell meta file OR marker gene list file OR both**  

**Optional a.** cell meta file (.csv)  
Note: This file contains three columns. prob means probability that the cell will be assigned to the cell type. this prob value is obtained from ICI method. If it is not from the ICI method or other methods that could give a probability, we can use 1 to represent prop.  
|  | cell_type | prob |
| -------- | -------- | -------- |
| cell1    | celltype1     | 1     |
| cell2    | celltype1     | 0.8     |
| cell3     | celltype2     | 0.6     |

**Optional b.** marker gene list file (.txt)  
Note: SPmarker will utilize a correlation-based method to predict cell identities that will be used for generate a meta file.
| gene | cell_type |
| -------- | -------- |
| gene1    | celltype1     |
| gene1    | celltype1     | 
| gene2     | celltype2     | 
| gene3     | celltype3     | 

**Optional c.** cell meta file (.csv) and marker gene list file (.txt)
Note: If users provide both files, the final output will return novel markers that do not contain the markers provided in the marker gene list file.

**3. (optional) unknown cell matrix (.csv)**  
UN means unknown
|  | UNcell1 | UNcell2 | UNcell3 |
| -------- | ----- | ----- | ----- |
| gene1    | 3.6   | 4.5     | 1.1
| gene2    | 7.3     | 0     |1.8
| gene3    | 2.1     | 0     |8.1

## **Example Run 1, Identify novel cell-type marker genes**
### Input data
**1. Gene expression matrix file**  
(ipt_test1_gene_cell_mtx.csv)  
**2. cell meta file**  
(ipt_test1_meta.csv)  
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_gene_cell_mtx.csv \\  
-meta ipt_test1_meta.csv  

## Other options:  
**Case 1:** if users do not have meta file and instead they have a marker gene list file (ipt_test1_marker_list.txt).  
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_gene_cell_mtx.csv \\  
-mlist ipt_test1_marker_list.txt  

**Case 2:** Change independent testing data size and cross-validation setting. SPmarker will divide the provided cells to training and independent datasets based on the ‘-indep_ratio’. The default setting is ‘0.1’ which means if users provide 1000 cells in the ‘-mtx’, SPmarker will sample 100 cells to be independent dataset. Also, the SPmarker automatically uses five-fold cross validation to do the training, and users can change it using ‘-cv_num’. If users want to use a set of genes/features other than all features during training, they can use ‘-feat_fl’.  
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_gene_cell_mtx.csv \\  
-meta ipt_test1_meta.csv \\  
-indep_ratio 0.1 \\  
-cv_num 5 \\  
-feat_fl ipt_test1_selected_features.csv  

**Case 3:** Identify novel markers. The above command lines could also generate a candidate new marker list with default top 20 markers for each cell type. Users can change the number by setting ‘-mar_num’. If users use ‘-SVM’, the SPmarker will return markers based on Support Vector Machine (SVM) approach. If users provide a known marker list using ‘-kmar_fl’, SPmarker will return a novel marker list that does not include the provided known markers.
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_gene_cell_mtx.csv \\  
-meta ipt_test1_meta.csv \\  
-indep_ratio 0.1 \\  
-cv_num 5 \\  
-mar_num 20 \\  
-SVM yes \\  
-kmar_fl ipt_test1_known_marker_list.txt  


## Example Run 2: Identify novel cell-type marker genes based on a GFP marker
When users provide a GFP marker name presented in the gene expression matrix, SPmarker can label cells by identifying cells where the GFP marker expressed. Briefly, for GFP-tagged cells, the SPmarker will label cells as ‘positive’ because these cells contain reads that can map to the GFP gene, and other cells without GFP reads as ‘negative’ examples. 
### Input data
**1. Gene expression matrix file**  
(ipt_test2_gene_cell_mtx.csv)  
**2. GFP marker gene name**  
(eg. GFP_marker, the name 'GFP_marker' should be one of feature names in the ipt_test2_exp.csv)
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test2_gene_cell_mtx.csv \\  
-m GFP_marker_name  
## Other options:  
**Case 1:** For assigning the ‘negative’ cells, SPmarker will keep balance of positive and negative cells by setting ‘-bns_ratio’. If users want to allow number of ‘negative’ cells to have two times of ‘positive’ cells, they can use ‘-bns_ratio 1:2’.   
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test2_gene_cell_mtx.csv \\  
-m GFP_marker_name \\  
-bns yes \\  
-bns_ratio 1:1  

## Example Run 3: Assign unknown cells  
After building the classifiers, SPmarker can assign the cell identities if users provide a gene expression matrix including unknown cells. If the genes have different orders as training matrix, the SPmarker will automatically keep the same feature orders between the unknown matrix and training matrix. If there are some genes missing in the unknown matrix compared to the training matrix, this tool will assign '0' across all cells for this gene.  
### Input data
**1. Gene expression matrix file**  
(ipt_test1_gene_cell_mtx.csv)  
**2. cell meta file**  
(ipt_test1_meta.csv)  
**3. unknown cell matrix**  
(ipt_test1_unknown_exp.csv)  
### Running
* > mkdir output_dir working_dir  
* > python SPmarker.py \\  
-d working_dir/ -o output_dir/ \\  
-mtx ipt_test1_gene_cell_mtx.csv \\  
-meta ipt_test1_meta.csv \\  
-ukn_mtx ipt_test1_unknown_exp.csv  


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
                 optional: [-lowestCellN 50][-bns no][-bns_ratio 1:1][-cv_num 5]
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

-lowestCellN      If users provide a marker list, please provide an additional cell number cutoff, 
                  which enables the predicted cell type having number of cells above this cutoff.
                  Default: 50.

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

-eval_score       Bascially, we will use the all the training dataset to identify markers. If users provide a specific evaluation score such as MCC in this argument, SPmarker will use this score to decide the best model that will be used for marker identification.
                  Default: All.

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


