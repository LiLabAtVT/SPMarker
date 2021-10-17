#!/usr/bin/env python

##this script is to use the marker list to predict the cell identity based on the correlation methods
##we will generate a file meta file
import re
import glob
import sys
import subprocess
import os

input_mtx_fl = sys.argv[1] ##04_generate_WER_training_012921/generate_mtx_013021/output_dir
input_output_dir = sys.argv[2]
marker_gene_list_fl = sys.argv[3]
input_correlation_Rscript_fl = sys.argv[4]

def prepare_datasets (input_mtx_fl,input_output_dir,marker_gene_list_fl):


    ##create a 1 0 matrix
    store_celltype_gene_dic = {}
    store_gene_dic = {}
    count = 0
    with open (marker_gene_list_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            celltype = col[1]
            temp_gene = col[0]

            count += 1
            if count != 1:

                if ' ' in temp_gene:
                    gene = col[0].replace(' ', '_')
                else:
                    gene = temp_gene

                if celltype in store_celltype_gene_dic:
                    store_celltype_gene_dic[celltype][gene] = 1
                else:
                    store_celltype_gene_dic[celltype] = {}
                    store_celltype_gene_dic[celltype][gene] = 1

                store_gene_dic[gene] = 1

    store_final_line_list = []
    first_line = 'Gene'
    for eachcelltype in store_celltype_gene_dic:
        first_line = first_line + '\t' + eachcelltype
    store_final_line_list.append(first_line)

    for eachgene in store_gene_dic:

        score_string = eachgene
        for eachcelltype in store_celltype_gene_dic:
            celltype_gene_dic = store_celltype_gene_dic[eachcelltype]

            if eachgene in celltype_gene_dic:
                score = '1'
            else:
                score = '0'

            score_string = score_string + '\t' + score

        store_final_line_list.append(score_string)

    with open (input_output_dir + '/temp_marker_celltype_mtx.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##prepare the gene by cell matrix
    if '/' in input_mtx_fl:
        mt = re.match('.+/(.+)',input_mtx_fl)
        flnm = mt.group(1)
    else:
        flnm = input_mtx_fl

    store_final_line_list = []
    count = 0
    with open (input_mtx_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')

            if '.csv' not in flnm:
                print("Please provide the 'csv' format matrix file")
                return

            col = eachline.strip().split(',')

            count += 1
            if count == 1:
                store_final_line_list.append(eachline)
            else:
                genenm = col[0]
                if genenm in store_gene_dic:
                    store_final_line_list.append(eachline)

    with open (input_output_dir + '/temp_gene_cell_exp_mtx.csv','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def make_correlation (input_output_dir,input_correlation_Rscript_fl):

    cmd = 'Rscript ' + input_correlation_Rscript_fl + \
          ' ' + input_output_dir + '/temp_gene_cell_exp_mtx.csv' + \
          ' ' + input_output_dir + '/temp_marker_celltype_mtx.txt' + \
          ' ' + input_output_dir
    subprocess.call(cmd,shell=True)


prepare_datasets (input_mtx_fl,input_output_dir,marker_gene_list_fl)
make_correlation (input_output_dir,input_correlation_Rscript_fl)







