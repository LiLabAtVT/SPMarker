#!/usr/bin/env python

##updating 051821 collect the SVM markers output
##updation060420 add other scripts to generate other datasets
##eg. the number of fts in each cell type how many cells are come from ORIM
##updation060320 debug
##this script will generate average fts from the previous 45 combinations of feats per class
from pandas import read_csv
import sys
import glob
import re
import pickle
from statistics import mean
import os
import pandas as pd

input_imp_five_cross_dir = sys.argv[1] ##this is the output dir from the generate_imp_from_built_SVM_m_060120
##store the imp from the five cross validation
##06_SVMM_markers_identify_013121/01_generate_imp_from_built_SVM_m_013121/output_dir

input_marker_gene_fl = sys.argv[2]
input_select_top_feat_num = sys.argv[3] ##select top 20

input_working_dir = sys.argv[4]
input_output_dir = sys.argv[5]

check_known_marker = sys.argv[6] ##yes or no


##define a function that calculate the average of imp for each svm in the cross validation
def average_imp_cross (input_imp_five_cross_dir,input_working_dir):

    ##generate_another dir to store the output that contains each cell type file
    ##this file is used to assign the cell type to the features
    temp_sort_svm_dir = input_working_dir + '/temp_sort_svm_dir'
    if not os.path.exists(temp_sort_svm_dir):  ##if dir not exit, it will create the dir
        os.makedirs(temp_sort_svm_dir)

    store_cell_type_dic = {} ##key is '1' and value is cell name
    store_all_feat_dic = {}
    col_len = 0
    opt_imp_svm_fl_list = glob.glob(input_imp_five_cross_dir + '/opt_imp_svm_*')
    for eachimp_fl in opt_imp_svm_fl_list:
        count = 0
        with open (eachimp_fl ,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                count += 1
                col = eachline.strip().split('\t')
                if count != 1:
                    store_all_feat_dic[col[0]] = 1
                    col_len = len(col)
                else:
                    ##for the title line
                    for i in range(1 ,len(col)):
                        ##store the cell name into a dic
                        store_cell_type_dic[str(i)] = col[i]

    ##for eachcell type
    for i in range(1 ,col_len): ##col_len is 11 since we have ten cell type

        store_final_line_list = []
        first_line = 'cell_name' + '\t' + 'feature' + '\t' + 'svm_value'
        store_final_line_list.append(first_line)

        ##cell name
        cell_name = store_cell_type_dic[str(i)]
        for eachft in store_all_feat_dic:
            ##for each ft
            store_cv_ft_val_list = []  ##store the imp for each ft from the five cross validation models
            for eachimp_fl in opt_imp_svm_fl_list:
                count = 0
                with open(eachimp_fl, 'r') as ipt:
                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        count += 1
                        if count != 1:
                            col = eachline.strip().split('\t')
                            if eachft == col[0]:
                                class_imp = col[i]
                                store_cv_ft_val_list.append(float(class_imp))

            ft_aver_imp = mean(store_cv_ft_val_list)
            other_line = cell_name + '\t' + eachft + '\t' + str(ft_aver_imp)
            store_final_line_list.append(other_line)

        with open (temp_sort_svm_dir + '/opt_' + cell_name + '_svm.txt' ,'w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')


        ##sort the svm file
        ##sort the svm value file
        svm_fl_dt = pd.read_table(temp_sort_svm_dir + '/opt_' + cell_name + '_svm.txt',
                                  delimiter=r"\s+")  ##has header so we use header=None shut down
        ##sort the value
        cell_type_svm_fl_sort_dt = svm_fl_dt.sort_values(by='svm_value', ascending=False)
        ##save sorted dt to output
        cell_type_svm_fl_sort_dt.to_csv(temp_sort_svm_dir + '/opt_' + cell_name + '_sort.txt',
                                        index=False,
                                        sep='\t')  ##no index

    return (temp_sort_svm_dir,store_all_feat_dic,store_cell_type_dic)

##step 2
def assign_cell_type_to_feature (temp_sort_svm_dir ,store_all_feat_dic ,input_working_dir):

    store_cell_type_name_dic = {}
    store_final_assign_cell_type_line_list = []
    ft_id = 0
    for eachft in store_all_feat_dic:

        ft_id += 1
        print('analyze ft id is ' + str(ft_id))

        store_cell_type_svm_dic = {} ##each cell type has one svm value

        ##extract the svm for each cell type
        fl_list = glob.glob(temp_sort_svm_dir + '/*_sort.txt')
        for eachfl in fl_list:
            # print(eachfl)
            mt = re.match('.+/(.+)', eachfl)
            fl_nm = mt.group(1)
            mt = re.match('opt_(.+)_sort\.txt', fl_nm)
            cell_type_nm = mt.group(1)
            store_cell_type_name_dic[cell_type_nm] = 1

            svm_value = ''
            count = 0
            with open(eachfl, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split('\t')

                    count += 1
                    if count != 1:  ##no column
                        gene_nm = col[1]

                        if eachft == gene_nm:
                            ##collect the svm

                            if len(col) == 3:
                                svm_value = col[2]
                                # print(svm_value)
                            else:
                                print('the line is wrong ' + eachline)

            if svm_value != '':
                store_cell_type_svm_dic[cell_type_nm] = float(svm_value)

            else:
                print('the svm value is empty ' + eachft)


        if len(list(store_cell_type_svm_dic.keys())) != 0:
            ##selec the key max as the cell type
            Keymax = max(store_cell_type_svm_dic, key=store_cell_type_svm_dic.get)

            final_line = Keymax + '\t' + eachft + '\t' + str(store_cell_type_svm_dic[Keymax])
            store_final_assign_cell_type_line_list.append(final_line)

    ##store the file a working_dir and divide the file into multiple files each files contains a cell type
    with open (input_working_dir + '/temp_assign_cell_type_fl.txt' ,'w+') as opt:
        for eachline in store_final_assign_cell_type_line_list:
            opt.write(eachline + '\n')


def create_sort_svmvalue (input_working_dir ,temp_assign_cell_type_fl):

    ##create a dir in the working_dir to store the sorted unique fts in it.
    select_top_cell_type_sort_uni_fs_dir = input_working_dir + '/select_top_cell_type_sort_uni_fs_dir'
    if not os.path.exists(select_top_cell_type_sort_uni_fs_dir):
        os.makedirs(select_top_cell_type_sort_uni_fs_dir)

    ##generate cell type dic
    store_cell_type_true_dic = {}
    with open(temp_assign_cell_type_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            cell_type_nm = col[0]
            store_cell_type_true_dic[cell_type_nm] = 1

    ##divide the file to multiple
    for eachcell_type in store_cell_type_true_dic:

        store_final_line_list = []
        with open (temp_assign_cell_type_fl ,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')
                cell_type_nm = col[0]

                if eachcell_type == cell_type_nm:
                    store_final_line_list.append(eachline)

        with open (input_working_dir + '/temp_' + eachcell_type + '_svm.txt' ,'w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        ##sort the temp svm file
        svm_fl_dt = pd.read_table(input_working_dir + '/temp_' + eachcell_type + '_svm.txt', header=None, delimiter=r"\s+")
        svm_fl_dt.columns = ['cell_name', 'feature', 'svm_value']
        cell_type_svm_fl_sort_dt = svm_fl_dt.sort_values(by='svm_value', ascending=False)
        cell_type_svm_fl_sort_dt.to_csv \
            (select_top_cell_type_sort_uni_fs_dir + '/opt_' + eachcell_type + '_uni_fts_sort.txt', index=False
            ,sep='\t')  ##no index


    return (select_top_cell_type_sort_uni_fs_dir)

##udpating 051921
####this function does not work now we will generate
##updation 051720 unique fts version
def generate_dir_to_store_sort_for_top_uni_fts (select_top_cell_type_sort_uni_fs_dir ,input_select_top_feat_num
                                               ,input_marker_gene_fl,check_known_marker):

    ##generate dir to store the specific top feat dir
    cell_top_feat_dir = input_output_dir + '/select_svm_top_feat_sort_uni_fs_dir'
    if not os.path.exists(cell_top_feat_dir):
        os.makedirs(cell_top_feat_dir)

    cell_top_feat_svm_dir = input_output_dir + '/select_svm_top_feat_sort_svm_uni_fs_dir'
    if not os.path.exists(cell_top_feat_svm_dir):
        os.makedirs(cell_top_feat_svm_dir)

    cell_top_feat_orim_dir = input_output_dir + '/select_svm_top_feat_sort_orim_uni_fs_dir'
    if not os.path.exists(cell_top_feat_orim_dir):
        os.makedirs(cell_top_feat_orim_dir)

    ##updation 052420 create a dir to store the shap markers that are up to required number
    cell_top_feat_svm_uptop_num_dir = input_output_dir + '/select_top_feat_sort_svm_uptop_num_uni_fs_dir'
    if not os.path.exists(cell_top_feat_svm_uptop_num_dir):
        os.makedirs(cell_top_feat_svm_uptop_num_dir)

    ##store the gene information
    store_marker_gene_dic = {}

    if check_known_marker == 'yes':
        count = 0
        with open (input_marker_gene_fl ,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                count += 1
                if count != 1:
                    col = eachline.strip().split('\t')
                    store_marker_gene_dic[col[0]] = 1


    ##updation 052420
    collect_all_feat_svm_uptop_num_dic = {}

    sort_fl_list = glob.glob(select_top_cell_type_sort_uni_fs_dir + '/*')
    for eachsort_fl in sort_fl_list:
        mt = re.match('.+/(.+)', eachsort_fl)
        fl_nm = mt.group(1)
        mt = re.match('opt_(.+)_uni_fts_sort\.txt', fl_nm)
        cell_type_nm = mt.group(1)

        store_top_feat_line_list = []
        store_top_feat_svm_line_list = []
        store_top_feat_orim_line_list = []

        ##updation 0524
        store_top_feat_svm_uptop_num_line_list = []

        ##updation 0524
        svm_count = -1

        count = -1
        with open (eachsort_fl ,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')

                count += 1
                if count != 0:
                    if count <= int(input_select_top_feat_num):
                        store_top_feat_line_list.append(eachline)

                        if col[1] in store_marker_gene_dic:
                            store_top_feat_orim_line_list.append(eachline)
                        else:
                            store_top_feat_svm_line_list.append(eachline)

                ##updation 052420
                ##do not contain any marker gene and select the top 20
                if col[1] not in store_marker_gene_dic:
                    svm_count += 1
                    if svm_count != 0:
                        if svm_count <= int(input_select_top_feat_num):
                            store_top_feat_svm_uptop_num_line_list.append(eachline)
                            collect_all_feat_svm_uptop_num_dic[col[1]] = 1




        with open (cell_top_feat_dir + '/opt_top_' + input_select_top_feat_num + '_' + cell_type_nm + '_sort.txt'
                  ,'w+') as opt:
            for eachline in store_top_feat_line_list:
                opt.write(eachline + '\n')

        with open (cell_top_feat_svm_dir + '/opt_top_' + input_select_top_feat_num + '_' + cell_type_nm + '_sort_svm.txt',
                'w+') as opt:
            for eachline in store_top_feat_svm_line_list:
                opt.write(eachline + '\n')

        with open(cell_top_feat_orim_dir + '/opt_top_' + input_select_top_feat_num + '_' + cell_type_nm + '_sort_orim.txt',
                'w+') as opt:
            for eachline in store_top_feat_orim_line_list:
                opt.write(eachline + '\n')

        ##updation 052420
        with open(cell_top_feat_svm_uptop_num_dir + '/opt_top_' + input_select_top_feat_num + '_' + cell_type_nm + '_sort_svm_uptop_num.txt',
                'w+') as opt:
            for eachline in store_top_feat_svm_uptop_num_line_list:
                opt.write(eachline + '\n')

    ##updation 052420
    with open(input_output_dir + '/opt_top_' + input_select_top_feat_num + '_svm.txt', 'w+') as opt:
        for eachft in collect_all_feat_svm_uptop_num_dic:
            opt.write(eachft + '\n')

##updating 051921 use the same scripts as SHAP to generate output
def generate_dir_to_store_sort_for_top_uni_fts_add_markerres (select_top_cell_type_sort_uni_fs_dir,
                                                input_select_top_feat_num,input_marker_gene_fl,check_known_marker):


    ##store the gene information
    store_marker_gene_dic = {}
    if check_known_marker == 'yes':
        count = 0
        with open (input_marker_gene_fl ,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                count += 1
                if count != 1:
                    col = eachline.strip().split('\t')
                    store_marker_gene_dic[col[0]] = 1

    ##updation 052620
    store_final_unis_top_feat_count_list = []

    ##store final top marker without known markers  eg. top 20 markers are all new
    store_final_top_marker_no_known_line_list = []
    ##store final top marker with known markers
    store_final_top_marker_with_known_line_list = []

    sort_fl_list = glob.glob(select_top_cell_type_sort_uni_fs_dir + '/*')
    for eachsort_fl in sort_fl_list:
        mt = re.match('.+/(.+)', eachsort_fl)
        fl_nm = mt.group(1)
        mt = re.match('opt_(.+)_uni_fts_sort\.txt', fl_nm)
        cell_type = mt.group(1)

        ##updation 052620
        cell_type_nm = ''
        if cell_type == 'Meri..Xylem':
            cell_type_nm = 'Meri_Xylem'
        if cell_type == 'Phloem..CC.':
            cell_type_nm = 'Phloem_CC'
        if cell_type == 'Cortext':
            cell_type_nm = 'Cortex'
        if cell_type != 'Meri..Xylem' and cell_type != 'Phloem..CC.' and cell_type != 'Cortext':
            cell_type_nm = cell_type

        ##updation 052620
        top_ori_marker_num = 0
        top_new_marker_num = 0

        ##updation 052420
        store_top_feat_shap_uptop_num_line_list = []

        ##updation 052420
        shap_count = -1

        count = -1
        with open(eachsort_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                count += 1
                if count != 0:
                    if count <= int(input_select_top_feat_num):

                        if col[1] in store_marker_gene_dic:

                            ##updation 052620
                            top_ori_marker_num += 1

                            with_known_marker_line = col[1] + '\t' + cell_type_nm + '\t' + 'Known' + '\t' + col[2]
                            store_final_top_marker_with_known_line_list.append(with_known_marker_line)

                        else:
                            ##updation 052620
                            top_new_marker_num += 1

                            with_known_marker_line = col[1] + '\t' + cell_type_nm + '\t' + 'Novel' + '\t' + col[2]
                            store_final_top_marker_with_known_line_list.append(with_known_marker_line)

                ##updation 052420
                ##do not contain any marker gene and select the top 20
                if col[1] not in store_marker_gene_dic:
                    shap_count += 1
                    if shap_count != 0:
                        if shap_count <= int(input_select_top_feat_num):
                            store_top_feat_shap_uptop_num_line_list.append(eachline)

                            ##updating 051821 add the value
                            no_known_marker_line = col[1] + '\t' + cell_type_nm + '\t' + 'Novel' + '\t' + col[2]
                            store_final_top_marker_no_known_line_list.append(no_known_marker_line)

        ##updation 052620 calculate the unis shap and orim markers
        ##store top marker gene information
        final_line = cell_type_nm + '\t' + 'Known_maker_feature' + '\t' + str(top_ori_marker_num) + '\n' + \
                     cell_type_nm + '\t' + 'Novel_marker_feature' + '\t' + str(top_new_marker_num)
        store_final_unis_top_feat_count_list.append(final_line)


    return (store_final_top_marker_no_known_line_list,store_final_top_marker_with_known_line_list,store_final_unis_top_feat_count_list)



temp_sort_svm_dir,store_all_feat_dic,store_cell_type_dic = average_imp_cross (input_imp_five_cross_dir,input_working_dir)

assign_cell_type_to_feature(temp_sort_svm_dir, store_all_feat_dic, input_working_dir)
select_top_cell_type_sort_uni_fs_dir = create_sort_svmvalue(input_working_dir,input_working_dir + '/temp_assign_cell_type_fl.txt')

#generate_dir_to_store_sort_for_top_uni_fts(select_top_cell_type_sort_uni_fs_dir, input_select_top_feat_num,
#                                           input_marker_gene_fl,check_known_marker)

##udpating 051921 write resutls same as SHAP markers
store_final_top_marker_no_known_line_list, \
    store_final_top_marker_with_known_line_list, \
    store_final_unis_top_feat_count_list = generate_dir_to_store_sort_for_top_uni_fts_add_markerres (select_top_cell_type_sort_uni_fs_dir,
                                                input_select_top_feat_num,input_marker_gene_fl,check_known_marker)

##write out results to the output_dir
with open (input_output_dir + '/opt_top_' + str(input_select_top_feat_num) + '_novel_marker.txt','w+') as opt:
    for eachline in store_final_top_marker_no_known_line_list:
        opt.write(eachline + '\n')

with open (input_output_dir + '/opt_top_' + str(input_select_top_feat_num) + '_novel_known_marker.txt', 'w+') as opt:
    for eachline in store_final_top_marker_with_known_line_list:
        opt.write(eachline + '\n')

with open(input_output_dir + '/opt_top_' + str(input_select_top_feat_num) + '_summary_marker_composition.txt', 'w+') as opt:
    for eachline in store_final_unis_top_feat_count_list:
        opt.write(eachline + '\n')




