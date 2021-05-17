#!/usr/bin/env python

##updating 013021 we need to keep balance for the data negative and postivie
##this script will correspond the exp mtx and meta
##generate a final exp and meta files with same genes and features
##the meta file has the largest number of cells so we only subset cells in the exp

import sys
import re
import glob
import os
import subprocess
import random

input_mtx_fl = sys.argv[1] ##00_data_integration_Seurat_013021/output_dir/exp_data_integrated_21264genes_12282cells_mtx.csv
input_meta_fl = sys.argv[2] ##01_divide_pos_neg_cells_012921/output_dir/opt_all_meta.csv
input_ratio = sys.argv[3] ##default is 1:1
input_output_dir = sys.argv[4] ##


def correspond_meta_exp (input_mtx_fl,input_meta_fl,input_output_dir):

    store_new_meta_line_list = []

    ##store_meta_infor
    store_meta_cells_dic = {}
    count = 0
    with open (input_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            count += 1
            if count != 1:
                store_meta_cells_dic[col[0]] = eachline
            else:
                store_new_meta_line_list.append(eachline)

    store_cells_dic = {}
    count = 0
    with open (input_mtx_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            count += 1
            if count == 1:
                for i in range(1,len(col)):
                    if '"' in col[i]:
                        mt = re.match('\"(.+)\"',col[i])
                        new_cellnm = mt.group(1)
                    else:
                        new_cellnm = col[i]

                    meta_line = store_meta_cells_dic[new_cellnm]
                    store_new_meta_line_list.append(meta_line)
                    store_cells_dic[new_cellnm] = 1

    cell_num = str(len(list(store_cells_dic.keys())))


    with open (input_output_dir + '/opt_all_meta_' + cell_num + 'cells.csv','w+') as opt:
        for eachline in store_new_meta_line_list:
            opt.write(eachline + '\n')

##select the balance information
def keep_balance (new_meta_fl,input_mtx_fl,input_output_dir,input_ratio):

    mt = re.match('(.+):(.+)',input_ratio)
    left_prop = mt.group(1)
    right_prop = mt.group(2)


    ########
    store_final_meta_line_list = []
    ########

    store_meta_cells_line_dic = {}
    count = 0
    pos_cell_list = []
    neg_cell_list = []
    with open (new_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            count += 1
            if count != 1:
                store_meta_cells_line_dic[col[0]] = eachline
                if col[1] == 'Positive':
                    pos_cell_list.append(col[0])
                else:
                    neg_cell_list.append(col[0])
            else:
                store_final_meta_line_list.append(eachline)

    ##compare number of pos cells and neg cells
    store_total_cellnm_dic = {}
    if len(pos_cell_list) <= len(neg_cell_list):

        ratio = float(right_prop) / float(left_prop)

        ##random the neg_cell_list
        random_neg_list = random.sample(neg_cell_list, len(neg_cell_list))
        pos_cell_num = len(pos_cell_list)

        subset_rd_neg_cellnm_list = []
        for i in range(int(pos_cell_num*ratio)):
            subset_rd_neg_cellnm_list.append(random_neg_list[i])

        for eachcell in pos_cell_list:
            store_total_cellnm_dic[eachcell] = 1
        for eachcell in subset_rd_neg_cellnm_list:
            store_total_cellnm_dic[eachcell] = 1

    else:

        ratio = float(left_prop) / float(right_prop)

        random_pos_list = random.sample(pos_cell_list,len(pos_cell_list))
        neg_cell_num = len(neg_cell_list)

        subset_rd_pos_cellnm_list = []
        for i in range(int(neg_cell_num*ratio)):
            subset_rd_pos_cellnm_list.append(random_pos_list[i])

        for eachcell in subset_rd_pos_cellnm_list:
            store_total_cellnm_dic[eachcell] = 1
        for eachcell in neg_cell_list:
            store_total_cellnm_dic[eachcell] = 1

    ##subset the neg cell exp file
    ##by store_total_cellnm_dic
    ########
    store_final_exp_line_list = []
    ########

    store_select_order_id_list = []
    store_select_cellnm_list = []
    count = 0
    with open (input_mtx_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            count += 1

            if count == 1:
                subset_first_line_str = ''
                for i in range(1,len(col)):
                    if col[i] in store_total_cellnm_dic:
                        store_select_order_id_list.append(i)
                        subset_first_line_str = subset_first_line_str + ',' + col[i]
                        store_select_cellnm_list.append(col[i])

                store_final_exp_line_list.append(subset_first_line_str)

            else:
                subset_col_str = col[0]
                for eachselect_id in store_select_order_id_list:
                    subset_col_str = subset_col_str + ',' + col[eachselect_id]
                store_final_exp_line_list.append(subset_col_str)

    for eachcell in store_select_cellnm_list:
        meta_line = store_meta_cells_line_dic[eachcell]
        store_final_meta_line_list.append(meta_line)

    final_meta_cell_num = len(store_final_meta_line_list) - 1
    with open (input_output_dir + '/opt_balance_meta_' + str(final_meta_cell_num) + 'cells.csv','w+') as opt:
        for eachline in store_final_meta_line_list:
            opt.write(eachline + '\n')

    final_exp_cell_num = len(store_select_cellnm_list)
    with open (input_output_dir + '/opt_balance_exp_' + str(final_exp_cell_num) + 'cells.csv','w+') as opt:
        for eachline in store_final_exp_line_list:
            opt.write(eachline + '\n')



correspond_meta_exp (input_mtx_fl,input_meta_fl,input_output_dir)
opt_fl_list = glob.glob(input_output_dir + '/opt_all_meta_*')
new_meta_fl = opt_fl_list[0]
keep_balance (new_meta_fl,input_mtx_fl,input_output_dir,input_ratio)
