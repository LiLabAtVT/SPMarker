#!/usr/bin/env python

##this script will select the marker genes with most contribution for each specific cell type


import pandas as pd


def select_contribution (input_shap_value_fl,input_working_dir,input_output_dir):

    ##store the name
    store_cell_type_nm_dic = {}
    with open (input_shap_value_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_cell_type_nm_dic[col[0]] = 1

    #shap_fl_dt = pd.read_table(input_shap_value_fl, header=None, delimiter=r"\s+")
    shap_fl_dt = pd.read_csv(input_shap_value_fl, header=None, sep='\t')

    shap_fl_dt.columns = ['cell_name','feature','shap_value']
    for eachcell_nm in store_cell_type_nm_dic:

        ##select the feature
        cell_type_shap_fl_dt = shap_fl_dt.loc[shap_fl_dt['cell_name'] == eachcell_nm]
        ##sort the value
        cell_type_shap_fl_sort_dt = cell_type_shap_fl_dt.sort_values(by='shap_value', ascending=False)
        ##save sorted dt to output
        cell_type_shap_fl_sort_dt.to_csv(input_output_dir + '/opt_' + eachcell_nm + '_sort.txt', index=False, sep='\t') ##no index

        ##cumulative shap
        cell_type_shap_fl_sort_cum_dt = cell_type_shap_fl_sort_dt['shap_value'].cumsum(axis=0)
        #total_shap_value = cell_type_shap_fl_sort_dt['shap_value'].sum()

        ##save to a file
        cell_type_shap_fl_sort_cum_dt.to_csv(input_working_dir + '/' + eachcell_nm + '_cum.txt', index=True, sep='\t')

        ##read a file
        #cell_type_shap_fl_sort_cum_dt = pd.read_table(input_working_dir + '/' + eachcell_nm + '_cum.txt', delimiter=r"\s+")
        #cell_type_shap_fl_sort_cum_dt = pd.read_csv(input_working_dir + '/' + eachcell_nm + '_cum.txt',
        #                                              sep="\t")

        ##divided by total_shape_value
        #cell_type_shap_fl_sort_cum_dt['shap_value'] = cell_type_shap_fl_sort_cum_dt['shap_value'].div(total_shap_value)

        ##save to a file
        ##use this output to draw figures
        ##we should have index that indicates the feature number
        #cell_type_shap_fl_sort_cum_dt.to_csv(input_working_dir + '/opt_' + eachcell_nm + '_cum_div.txt', index=True, sep='\t')

        ##modify the cum_div.txt that will fit the format in R
        #store_final_line_list = []
        #count = 0
        #feat_count = 0
        #with open (input_working_dir + '/opt_' + eachcell_nm + '_cum_div.txt','r') as ipt:
        #    for eachline in ipt:
        #        eachline = eachline.strip('\n')
        #        count += 1
        #        if count != 1:
        #            feat_count += 1
        #            col = eachline.strip().split()
        #            final_line = str(feat_count) + '\t' + col[1]
        #            store_final_line_list.append(final_line)

        #with open (input_output_dir + '/opt_' + eachcell_nm + '_cum_div.txt','w+') as opt:
        #    for eachline in store_final_line_list:
        #        opt.write(eachline + '\n')




