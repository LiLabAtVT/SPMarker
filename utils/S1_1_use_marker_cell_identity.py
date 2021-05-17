#!/usr/bin/env python

##updating 020121 set the target gene as labled gene
##this script we will divide the three replicates to pos and neg
import sys
import re
import glob
import os
import subprocess

input_mtx_fl = sys.argv[1] ##04_generate_WER_training_012921/generate_mtx_013021/output_dir
input_output_dir = sys.argv[2]
target_labled_gene = sys.argv[3]
#opt_nm = sys.argv[4] ##WERGFP

#target_labled_gene = 'AT5G14750.Araport11.447' ##previous is GFP

##the aim is to generate the meta file first
def divide_pos_neg_cells_meta (input_mtx_fl,input_output_dir,target_labled_gene):

    store_final_line_list = []
    first_line = ',cell_type,prob'
    store_final_line_list.append(first_line)
    #fl_list = glob.glob(input_mtx_dir + '/*mtx.csv')
    #for eachfl in fl_list:
    #mt = re.match('.+/(.+)',eachfl)
    #fl_nm = mt.group(1)
    #mt = re.match('.+_(\d)_mtx.csv',fl_nm)
    #rep_nm = mt.group(1)

    store_all_line_dic = {} ##key is the first item and value should be the whole line
    count = 0
    with open (input_mtx_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            count += 1
            if count != 1:
                store_col_id_dic = {}##key is the order id and value is the count
                for i in range(1,len(col)):
                    store_col_id_dic[str(i)] = col[i]

                if re.match('"(.+)"',col[0]):
                    mt = re.match('"(.+)"',col[0])
                    gene_nm = mt.group(1)
                else:
                    gene_nm = col[0]
                store_all_line_dic[gene_nm] = store_col_id_dic
            else:
                store_cellnm_id_dic = {}
                for i in range(1,len(col)):

                    if re.match('"(.+)"',col[i]):
                        mt = re.match('"(.+)"',col[i])
                        cell_nm = mt.group(1)
                    else:
                        cell_nm = col[i]
                    store_cellnm_id_dic[str(i)] = cell_nm
                store_all_line_dic['cellnm'] = store_cellnm_id_dic

    gpf_line_dic = store_all_line_dic[target_labled_gene]

    store_pos_order_list = []
    store_neg_order_list = []
    store_pos_cell_nm_list = []
    store_neg_cell_nm_list = []

    for eachid in gpf_line_dic:
        cellnm = store_all_line_dic['cellnm'][eachid]
        gpf_count = gpf_line_dic[eachid]
        if float(gpf_count) == 0:
            store_neg_order_list.append(eachid)
            store_neg_cell_nm_list.append(cellnm)
        else:
            store_pos_order_list.append(eachid)
            store_pos_cell_nm_list.append(cellnm)

    ##obtain the name of cell nm
    for eachcell in store_pos_cell_nm_list:
        final_cellnm = eachcell
        final_line = final_cellnm + ',Positive,1'
        store_final_line_list.append(final_line)

    for eachcell in store_neg_cell_nm_list:
        final_cellnm = eachcell
        final_line = final_cellnm + ',Negative,1'
        store_final_line_list.append(final_line)


    with open (input_output_dir + '/opt_all_meta.csv','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

divide_pos_neg_cells_meta (input_mtx_fl,input_output_dir,target_labled_gene)
