#!/usr/bin/env python

##updating 051821 add the SHAP value to the output fl
##updating 031521 check string to be float
##updation 052620 change all the cell type name
##updation 052620 calculate the how many features belong to shap and how many features belong to known markers
##updation 052420 select the top number of shap that do not contain any original marker
##and randomly select marker to generate another files in a dir that will be used to compare performance by retraining models

##updation 052420 use all the collected markers
##updation 051720 rm the duplicate fts by calculating the proportion of the SHAP value in each cell type
##for each feature
##comapre shap
##assign feature to cell type
##select the top features
##generate a dir to contain different files


##updation 051720 add another fts detected from Seurat
##updation 042620 generate a dir to store the specific top number of features
##updation 042320 generate features specific to each cell type
##updation 042220 generate a pan features under the threshold mentioned in this script
##this script will compare number of markers and none markers
##generate a format for the R
##cell_type_nm \t type_marker \t number

import re
import glob
import os
import pandas as pd


def isfloat(value):
  try:
    float(value)
    return 'True'
  except ValueError:
    return 'False'


def assign_cell_type_to_feature (input_shap_dir):

    ##collect all the fts
    store_ft_dic = {}
    fl_list = glob.glob(input_shap_dir + '/*_sort.txt')
    count = 0
    with open(fl_list[0], 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                store_ft_dic[col[1]] = 1

    store_cell_type_name_dic = {}
    store_final_assign_cell_type_line_list = []
    ft_id = 0
    for eachft in store_ft_dic:

        ft_id += 1
        print('analyze ft id is ' + str(ft_id))

        store_cell_type_shap_dic = {} ##each cell type has one shap value

        ##extract the shap for each cell type
        fl_list = glob.glob(input_shap_dir + '/*_sort.txt')
        for eachfl in fl_list:
            #print(eachfl)
            mt = re.match('.+/(.+)', eachfl)
            fl_nm = mt.group(1)
            mt = re.match('opt_(.+)_sort\.txt', fl_nm)
            cell_type = mt.group(1)

            ##updation 052620
            cell_type_nm = ''
            if cell_type == 'Meri..Xylem' or cell_type == 'Meri.Xylem':
                cell_type_nm = 'Meri_Xylem'
            if cell_type == 'Phloem..CC.':
                cell_type_nm = 'Phloem_CC'
            if cell_type == 'Cortext':
                cell_type_nm = 'Cortex'
            if cell_type != 'Meri..Xylem' and cell_type != 'Meri.Xylem' and cell_type != 'Phloem..CC.' and cell_type != 'Cortext':
                cell_type_nm = cell_type

            store_cell_type_name_dic[cell_type_nm] = 1

            shap_value = ''
            count = 0
            with open(eachfl, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split()

                    count += 1
                    if count != 1:  ##no column
                        gene_nm = col[1]

                        if eachft == gene_nm:
                            ##collect the shap
                            shap_value = col[2]
                            #print(shap_value)

            ##updating 031521
            check_res = isfloat(shap_value)
            if check_res == 'True':
                store_cell_type_shap_dic[cell_type_nm] = float(shap_value)
            else:
                store_cell_type_shap_dic[cell_type_nm] = float('0')


        ##selec the key max as the cell type
        Keymax = max(store_cell_type_shap_dic, key=store_cell_type_shap_dic.get)

        final_line = Keymax + '\t' + eachft + '\t' + str(store_cell_type_shap_dic[Keymax])
        store_final_assign_cell_type_line_list.append(final_line)

    return (store_final_assign_cell_type_line_list)


def create_sort_shapvalue (input_shap_dir,input_working_dir,temp_assign_cell_type_fl):

    store_cell_type_name_dic = {}
    fl_list = glob.glob(input_shap_dir + '/*_sort.txt')
    for eachfl in fl_list:
        # print(eachfl)
        mt = re.match('.+/(.+)', eachfl)
        fl_nm = mt.group(1)
        mt = re.match('opt_(.+)_sort\.txt', fl_nm)
        cell_type = mt.group(1)

        ##updation 052620
        cell_type_nm = ''
        if cell_type == 'Meri..Xylem' or cell_type == 'Meri.Xylem':
            cell_type_nm = 'Meri_Xylem'
        if cell_type == 'Phloem..CC.':
            cell_type_nm = 'Phloem_CC'
        if cell_type == 'Cortext':
            cell_type_nm = 'Cortex'
        if cell_type != 'Meri..Xylem' and cell_type != 'Meri.Xylem' and cell_type != 'Phloem..CC.' and cell_type != 'Cortext':
            cell_type_nm = cell_type

        store_cell_type_name_dic[cell_type_nm] = 1

    ##create a dir in the working_dir to store the sorted unique fts in it.
    select_top_cell_type_sort_uni_fs_dir = input_working_dir + '/select_top_cell_type_sort_uni_fs_dir'
    if not os.path.exists(select_top_cell_type_sort_uni_fs_dir):
        os.makedirs(select_top_cell_type_sort_uni_fs_dir)

    ##divide the file to multiple
    for eachcell_type in store_cell_type_name_dic:

        store_final_line_list = []
        with open (temp_assign_cell_type_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                cell_type = col[0]

                ##updation 052620
                cell_type_nm = ''
                if cell_type == 'Meri..Xylem' or cell_type == 'Meri.Xylem':
                    cell_type_nm = 'Meri_Xylem'
                if cell_type == 'Phloem..CC.':
                    cell_type_nm = 'Phloem_CC'
                if cell_type == 'Cortext':
                    cell_type_nm = 'Cortex'
                if cell_type != 'Meri..Xylem' and cell_type != 'Meri.Xylem' and cell_type != 'Phloem..CC.' and cell_type != 'Cortext':
                    cell_type_nm = cell_type

                if eachcell_type == cell_type_nm:
                    store_final_line_list.append(eachline)

        with open (input_working_dir + '/temp_' + eachcell_type + '_shap.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        ##sort the temp shap file
        shap_fl_dt = pd.read_csv(input_working_dir + '/temp_' + eachcell_type + '_shap.txt', header=None, delimiter=r"\s+")
        shap_fl_dt.columns = ['cell_name', 'feature', 'shap_value']
        cell_type_shap_fl_sort_dt = shap_fl_dt.sort_values(by='shap_value', ascending=False)
        cell_type_shap_fl_sort_dt.to_csv(select_top_cell_type_sort_uni_fs_dir + '/opt_' + eachcell_type + '_uni_fts_sort.txt', index=False,sep='\t')  ##no index


    return (select_top_cell_type_sort_uni_fs_dir)


##updation 051720 unique fts version
def generate_dir_to_store_sort_for_top_uni_fts (select_top_cell_type_sort_uni_fs_dir,
                                                input_select_top_feat_num,input_marker_gene_fl):

    ##store the gene information
    store_marker_gene_dic = {}
    if input_marker_gene_fl != '':
        count = 0
        with open(input_marker_gene_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                count += 1
                if count != 1:
                    col = eachline.strip().split()
                    store_marker_gene_dic[col[0]] = 1

    ##updation 052620
    store_final_unis_top_feat_count_list = []
    ##store final top marker without known markers  eg. top 20 markers are all new
    store_final_top_marker_no_known_line_list = []
    ##store final top marker with known markers
    store_final_top_marker_with_known_line_list = []

    ##updating 051921
    ##save all the marker instead of the top markers
    store_final_unis_feat_count_list = []
    store_final_marker_no_known_line_list = []
    store_final_marker_with_known_line_list = []


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

        ##updating 051921
        ori_marker_num = 0
        new_marker_num = 0

        ##updation 052420
        store_top_feat_shap_uptop_num_line_list = []

        ##updating 051921
        store_feat_shap_uptop_num_line_list = []

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

                    #################
                    ##udpating 051921
                    if col[1] in store_marker_gene_dic:
                        ##updation 052620
                        ori_marker_num += 1
                        with_known_marker_line = col[1] + '\t' + cell_type_nm + '\t' + 'Known' + '\t' + col[2]
                        store_final_marker_with_known_line_list.append(with_known_marker_line)
                    else:
                        ##updation 052620
                        new_marker_num += 1
                        with_known_marker_line = col[1] + '\t' + cell_type_nm + '\t' + 'Novel' + '\t' + col[2]
                        store_final_marker_with_known_line_list.append(with_known_marker_line)

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

                #################
                ##updating 051921
                ##consider all features
                if col[1] not in store_marker_gene_dic:
                    shap_count += 1
                    if shap_count != 0:
                        store_feat_shap_uptop_num_line_list.append(eachline)
                        ##updating 051821 add the value
                        no_known_marker_line = col[1] + '\t' + cell_type_nm + '\t' + 'Novel' + '\t' + col[2]
                        store_final_marker_no_known_line_list.append(no_known_marker_line)


        ##updation 052620 calculate the unis shap and orim markers
        ##store top marker gene information
        final_line = cell_type_nm + '\t' + 'Known_maker_feature' + '\t' + str(top_ori_marker_num) + '\n' + \
                     cell_type_nm + '\t' + 'Novel_marker_feature' + '\t' + str(top_new_marker_num)
        store_final_unis_top_feat_count_list.append(final_line)

        ##updating 051921
        ##consider all features
        final_line = cell_type_nm + '\t' + 'Known_maker_feature' + '\t' + str(ori_marker_num) + '\n' + \
                     cell_type_nm + '\t' + 'Novel_marker_feature' + '\t' + str(new_marker_num)
        store_final_unis_feat_count_list.append(final_line)

    return (store_final_top_marker_no_known_line_list,store_final_top_marker_with_known_line_list,store_final_unis_top_feat_count_list,
            store_final_marker_no_known_line_list,store_final_marker_with_known_line_list,store_final_unis_feat_count_list)


