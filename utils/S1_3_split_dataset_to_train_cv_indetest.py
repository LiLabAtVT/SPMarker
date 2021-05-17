#!/usr/bin/env python

####do not know why did this

##updating 021121 add a new independent labeled by GFP this performance is the same as the previous one without labeled by GFP
##

##updating 020521 generate more kinds of dataset to train
##updating 020221 we need to set the list to remove the AT or GFP or both genes
##updating 013121 add a function to remove the GFP gene and AT5G14750
##updating 013021 this script aims to split the exp and meta to independent testing trainign and test as well as cross validation
##set to the remain number
##updation 041120 add all meta output
##updation 040120
# this script will randomly remove the 50% selected single cell

##updation 012720 decrease the number of category with large number of single number
##first we sort the value of the category and extract the first xx part of data

##this script is to divide the dataset into test and train dataset
## 5% 95%
##divide each group into 5% and 95%
##generate 5% random id and 95% remained id

import os
import pandas as pd
import random
import numpy as np
from itertools import islice
import sys
import subprocess
import glob
import re

input_dataset_exp_fl = sys.argv[1] ##00_data_integration_Seurat_013021/output_dir/exp_data_integrated_21264genes_12282cells_mtx.csv
input_meta_fl = sys.argv[2] ##02_correspond_meta_exp_genes_cells_013021/output_dir/opt_all_meta_12282cells.csv
fold_num = sys.argv[3] ##default is 5
input_thr_par = sys.argv[4] ##default is 0.1
input_output_dir = sys.argv[5]


input_fl_out_fam_list = []
cell_type_dic = {}

def split_train_test(input_dataset_exp_fl, input_meta_fl, input_thr_par, input_opt_dir, input_fl_out_fam_list,
                     cell_type_dic):
    ##generate a temp dir to store the temp meta data for the cell type whose size will be decreased
    opt_temp_store_meta_dir = input_opt_dir + '/opt_temp_store_meta_dir'
    if not os.path.exists(opt_temp_store_meta_dir):
        os.makedirs(opt_temp_store_meta_dir)

    ##store sg cell name in each family in a dic
    ##key is the cell name
    fam_nm_dic = {}
    count = 0
    with open(input_meta_fl, 'r') as ipt:
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            if count != 1:
                col = eachline.strip().split(',')
                fam_nm = col[1]
                fam_nm_dic[fam_nm] = 1

    ##value is the list of sg cell name
    total_sg_cell_list = []

    train_sg_cell_list = []
    test_sg_cell_list = []
    for eachfamnm in fam_nm_dic:

        if eachfamnm not in input_fl_out_fam_list:  ##it does not work

            total_fam_count = 0
            sg_cell_list = []
            ##updaiton 012720
            ##check if the eachfamnm is in the cell_type_dic
            ##if this generate a temp dir to sort the data from this cell type
            if eachfamnm in cell_type_dic:

                raw_meta_fl_dt = pd.read_csv(input_meta_fl, header=0, index_col=0)
                ##select the specific fam name
                famn_dt = raw_meta_fl_dt.loc[raw_meta_fl_dt['cell_type'] == eachfamnm]

                ##updation 040120 shuffle the dataframe
                shf_famn_dt = famn_dt.iloc[np.random.permutation(len(famn_dt))]

                ##sort the value in the famn_dt
                # sort_famn_dt = famn_dt.sort_values(by='prob',ascending=False)

                ##save the dt
                shf_famn_dt.to_csv(opt_temp_store_meta_dir + '/' + eachfamnm + '.csv', index=True)

                ##decrease the size of specific family
                remain_cell_num = cell_type_dic[eachfamnm]
                #remain_cell_num = de_por * int(shf_famn_dt.shape[0])

                de_shf_famn_dt = shf_famn_dt.head(int(remain_cell_num))
                ##save to the dir
                de_shf_famn_dt.to_csv(opt_temp_store_meta_dir + '/' + eachfamnm + '_final.csv', index=True)

                ##read the new meta fl
                count = 0
                with open(opt_temp_store_meta_dir + '/' + eachfamnm + '_final.csv', 'r') as ipt:
                    for eachline in ipt:
                        count += 1
                        eachline = eachline.strip('\n')
                        if count != 1:
                            col = eachline.strip().split(',')
                            if eachfamnm == col[1]:
                                total_fam_count += 1
                                sg_cell_nm = col[0]
                                sg_cell_list.append(sg_cell_nm)

            else:
                ##this type will not be analyzed
                count = 0
                with open(input_meta_fl, 'r') as ipt:
                    for eachline in ipt:
                        count += 1
                        eachline = eachline.strip('\n')
                        if count != 1:
                            col = eachline.strip().split(',')
                            if eachfamnm == col[1]:
                                total_fam_count += 1
                                sg_cell_nm = col[0]
                                sg_cell_list.append(sg_cell_nm)

            ##shuffle the list
            random.shuffle(sg_cell_list)
            fam_count = 0

            for eachid in sg_cell_list:
                fam_count += 1
                if fam_count <= int(total_fam_count * float(input_thr_par)):
                    test_sg_cell_list.append(eachid)
                else:
                    train_sg_cell_list.append(eachid)

                total_sg_cell_list.append(eachid)

    ##target on the meta and exp data
    ##for the exp data
    raw_exp_fl_dt = pd.read_csv(input_dataset_exp_fl, header=0, index_col=0)
    exp_test_id_dt = raw_exp_fl_dt[test_sg_cell_list]
    exp_test_id_dt.to_csv(input_opt_dir + '/opt_exp_test.csv', index=True)
    exp_train_id_dt = raw_exp_fl_dt[train_sg_cell_list]
    exp_train_id_dt.to_csv(input_opt_dir + '/opt_exp_train.csv', index=True)

    ##updaton 040120 set a total sg cell list to do the feature selection
    exp_id_dt = raw_exp_fl_dt[total_sg_cell_list]
    exp_id_dt.to_csv(input_opt_dir + '/opt_exp_all.csv', index=True)

    ##for the meta data
    raw_meta_fl_dt = pd.read_csv(input_meta_fl, header=0, index_col=0)
    meta_test_id_dt = raw_meta_fl_dt.loc[test_sg_cell_list]
    meta_test_id_dt.to_csv(input_opt_dir + '/opt_meta_test.csv', index=True)
    meta_train_id_dt = raw_meta_fl_dt.loc[train_sg_cell_list]
    meta_train_id_dt.to_csv(input_opt_dir + '/opt_meta_train.csv', index=True)

    ##updation 041120 add opt_meta_all.csv
    meta_id_dt = raw_meta_fl_dt.loc[total_sg_cell_list]
    meta_id_dt.to_csv(input_opt_dir + '/opt_meta_all.csv', index=True)

    ##check the number of meta data
    store_meta_test_nm_dic = {}
    store_meta_train_nm_dic = {}
    with open(input_opt_dir + '/opt_meta_test.csv', 'r') as ipt:
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            if count != 1:
                col = eachline.strip().split(',')
                if col[1] in store_meta_test_nm_dic:
                    store_meta_test_nm_dic[col[1]] += 1
                else:
                    store_meta_test_nm_dic[col[1]] = 1

    with open(input_opt_dir + '/opt_meta_train.csv', 'r') as ipt:
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            if count != 1:
                col = eachline.strip().split(',')
                if col[1] in store_meta_train_nm_dic:
                    store_meta_train_nm_dic[col[1]] += 1
                else:
                    store_meta_train_nm_dic[col[1]] = 1

    ##change the col name of the opt_fam_nm_num.txt
    with open(input_opt_dir + '/opt_fam_nm_num.txt', 'w+') as opt:
        for nm in store_meta_test_nm_dic:
            if nm == 'cell_type':
                opt.write('Cell_type' + '\t' + 'Validation' + '\t' + 'Training' + '\n')
            else:
                opt.write(nm + '\t' + str(store_meta_test_nm_dic[nm]) + '\t' + \
                      str(store_meta_train_nm_dic[nm]) + '\n')

def split_data_five_fold_cross_val (input_meta_fl,input_exp_fl,input_working_dir,input_output_dir,fold_num):

    ##store sg cell name in each family in a dic
    ##key is the cell name
    fam_nm_dic = {}
    count = 0
    with open(input_meta_fl, 'r') as ipt:
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            if count != 1:
                col = eachline.strip().split(',')
                fam_nm = col[1]
                fam_nm_dic[fam_nm] = 1

    ##mk a dir to store the replicates from each single cell type
    divide_family_output_dir = input_working_dir + '/divide_family_output'
    if not os.path.exists(divide_family_output_dir):
        os.makedirs(divide_family_output_dir)

    #############################
    ##store all the path together
    store_famnm_split_cnt_dic = {}  ##key is the famnm value is another dic whose key is str(count) and value is the path of csv
    for eachfamnm in fam_nm_dic:

        ##create a fam dir
        fam_divide_dir = input_working_dir + '/' + eachfamnm + '_dir'
        if not os.path.exists(fam_divide_dir):
            os.makedirs(fam_divide_dir)


        ##work in the fam_divide_dir
        raw_meta_fl_dt = pd.read_csv(input_meta_fl, header=0, index_col=0)
        ##select the specific fam name
        famn_dt = raw_meta_fl_dt.loc[raw_meta_fl_dt['cell_type'] == eachfamnm]
        ##generate shuffle dt
        shf_famn_dt = famn_dt.iloc[np.random.permutation(len(famn_dt))]
        ##save the dt
        shf_famn_dt.to_csv(fam_divide_dir + '/shuffle_' + eachfamnm + '.csv', index=True)

        ##collect the index of shf_famn_dt to one list
        index_list = list(shf_famn_dt.index.values)

        num, div = len(index_list), int(fold_num) ##five fold cross validation
        length_to_split = [num // div + (1 if x < num % div else 0) for x in range(div)]

        Inputt = iter(index_list)
        split_index_list = [list(islice(Inputt, elem)) for elem in length_to_split]
        ##this split_index_list contains three lists of each contains

        ##split the dataframes into five dataframes
        split_count = 0
        ##store the split dataframe path into a dic
        store_split_path_dic = {}
        for eachindex_list in split_index_list:
            split_shf_famn_dt = shf_famn_dt.loc[eachindex_list]
            split_count += 1
            ##write to csv file
            split_shf_famn_dt.to_csv(fam_divide_dir + '/split_' + str(split_count) + '_' + eachfamnm + '.csv',index=True)
            store_split_path_dic[str(split_count)] = fam_divide_dir + '/split_' + str(split_count) + '_' + eachfamnm + '.csv'

        ##store the store_split_path_dic to the store_famnm_split_cnt_dic
        store_famnm_split_cnt_dic[eachfamnm] = store_split_path_dic

    ########################
    ##generate five datasets
    letter_list = ['a', 'b', 'c', 'd', 'e','f','g','h','i','j']
    let_count = 0
    split_count = 0
    for each_let_nm in letter_list:
        let_count += 1

        ##do not exceed the fold number
        if let_count <= int(fold_num):

            split_count += 1

            ##generate test dt
            test_dt_frames = []
            train_dt_frames = []
            for eachfamnm in store_famnm_split_cnt_dic:
                store_split_path_dic = store_famnm_split_cnt_dic[eachfamnm]
                for eachsplit_str in store_split_path_dic:
                    if eachsplit_str == str(split_count):
                        fam_split_path = store_split_path_dic[eachsplit_str]
                        fam_split_dt = pd.read_csv(fam_split_path, header=0, index_col=0)
                        ##combine all the fam split path together
                        ##this is for the test
                        test_dt_frames.append(fam_split_dt)

                    ##it means b c d e if the previous one is a
                    else:
                        fam_split_path = store_split_path_dic[eachsplit_str]
                        fam_split_dt = pd.read_csv(fam_split_path, header=0, index_col=0)
                        train_dt_frames.append(fam_split_dt)

            combine_test_frame_dt = pd.concat(test_dt_frames)
            combine_train_frame_dt = pd.concat(train_dt_frames)
            ##write out frame dataframe
            ##mk a dir to store the data from combining family
            combine_family_output_dir = input_output_dir + '/' + each_let_nm
            if not os.path.exists(combine_family_output_dir):
                os.makedirs(combine_family_output_dir)

            combine_test_frame_dt.to_csv(combine_family_output_dir + '/opt_meta_test.csv', index=True)
            combine_train_frame_dt.to_csv(combine_family_output_dir + '/opt_meta_train.csv',index=True)

            ##generate expression
            ##extract index list
            ##for the testing
            exp_fl_dt = pd.read_csv(input_exp_fl, header=0, index_col=0)
            combine_test_dt_index_list = list(combine_test_frame_dt.index.values)
            exp_test_dt = exp_fl_dt[combine_test_dt_index_list]
            exp_test_dt.to_csv(combine_family_output_dir + '/opt_exp_test.csv', index=True)

            ##for the training
            combine_train_dt_index_list = list(combine_train_frame_dt.index.values)
            exp_train_dt = exp_fl_dt[combine_train_dt_index_list]
            exp_train_dt.to_csv(combine_family_output_dir + '/opt_exp_train.csv', index=True)


def split_dataset (input_dataset_exp_fl,input_meta_fl,input_output_dir):
    ########################
    ##begin to split dataset
    step1_split_train_indetest_opt_dir = input_output_dir + '/step1_split_train_indetest_opt_dir'
    if not os.path.exists(step1_split_train_indetest_opt_dir):
        os.makedirs(step1_split_train_indetest_opt_dir)
    step2_split_train_cross_val_opt_dir = input_output_dir + '/step2_split_train_cross_val_opt_dir'
    if not os.path.exists(step2_split_train_cross_val_opt_dir):
        os.makedirs(step2_split_train_cross_val_opt_dir)
    ##generate working and output directories
    S2_split_cross_val_w_dir = step2_split_train_cross_val_opt_dir + '/working_dir'
    if not os.path.exists(S2_split_cross_val_w_dir):
        os.makedirs(S2_split_cross_val_w_dir)
    S2_split_cross_val_o_dir = step2_split_train_cross_val_opt_dir + '/output_dir'
    if not os.path.exists(S2_split_cross_val_o_dir):
        os.makedirs(S2_split_cross_val_o_dir)

    #####################################
    ##use the generated opt_exp_train.csv to conduct the cross validation
    split_train_test(input_dataset_exp_fl, input_meta_fl, input_thr_par, step1_split_train_indetest_opt_dir, input_fl_out_fam_list,cell_type_dic)
    cmd = 'cp ' + step1_split_train_indetest_opt_dir + '/opt_fam_nm_num.txt ' + input_output_dir + '/opt_train_validation_cell_type_number.txt'
    subprocess.call(cmd, shell=True)
    input_exp_fl = step1_split_train_indetest_opt_dir + '/opt_exp_train.csv'
    input_meta_fl = step1_split_train_indetest_opt_dir + '/opt_meta_train.csv'
    input_inde_exp_fl = step1_split_train_indetest_opt_dir + '/opt_exp_test.csv'
    input_inde_meta_fl = step1_split_train_indetest_opt_dir + '/opt_meta_test.csv'

    ##################################
    ##cross validation dataset prepare
    split_data_five_fold_cross_val(input_meta_fl, input_exp_fl, S2_split_cross_val_w_dir, S2_split_cross_val_o_dir,int(fold_num))
    ##copy indep testing data to output dir. the data will be used for detecting shap marker
    cmd = 'cp ' + input_inde_exp_fl + ' ' + input_output_dir + '/opt_exp_indep_test.csv'
    subprocess.call(cmd, shell=True)
    cmd = 'cp ' + input_inde_meta_fl + ' ' + input_output_dir + '/opt_meta_indep_test.csv'
    subprocess.call(cmd, shell=True)
    ##cp the independent exp and meta fl to the output dir from the the split cross validation dir
    cross_vali_dir_list = glob.glob(S2_split_cross_val_o_dir + '/*')
    for eachcross_dir in cross_vali_dir_list:
        ##eachcross is a b c d e dir
        cmd = 'cp ' + input_inde_exp_fl + ' ' + eachcross_dir + '/opt_exp_indep_test.csv'
        subprocess.call(cmd, shell=True)
        cmd = 'cp ' + input_inde_meta_fl + ' ' + eachcross_dir + '/opt_meta_indep_test.csv'
        subprocess.call(cmd, shell=True)


split_dataset (input_dataset_exp_fl,input_meta_fl,input_output_dir)