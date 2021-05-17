#!/usr/bin/env python

##updating 031821 we need to make sure we have filter out cells from cookie's data
##also we need to keep the WER cells expressed in the ori five datasets

import re
import glob
import sys
import subprocess
import os
import pandas as pd
import numpy as np

input_unknown_mtx_fl = sys.argv[1]

input_pre_exp_fl = sys.argv[2]
##look for the features
##Step1_3_split_data_o_dir/opt_exp_indep_test.csv

input_output_dir = sys.argv[3]


def remove_additional_features(input_unknown_mtx_fl, input_pre_exp_fl,input_output_dir):

    store_target_fts_list = []
    count = 0
    with open(input_pre_exp_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            count += 1
            if count != 1:
                store_target_fts_list.append(col[0])

    store_ft_exp_line_dic = {}
    first_line = ''
    count = 0
    cellnum = 0
    with open(input_unknown_mtx_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split(',')
            cellnum = len(col) - 1
            count += 1
            if count == 1:
                first_line = eachline
            else:
                store_ft_exp_line_dic[col[0]] = eachline

    ##generate a final exp file
    store_final_line_list = []
    store_final_line_list.append(first_line)
    no_ft_count = 0
    for eachft in store_target_fts_list:

        if eachft in store_ft_exp_line_dic:
            final_line = store_ft_exp_line_dic[eachft]
            store_final_line_list.append(final_line)
        else:
            ##simulate 0 for all the cells
            no_ft_count += 1
            other_line = eachft
            for i in range(cellnum):
                other_line = other_line + ',' + '0'
            store_final_line_list.append(other_line)
    print('after checking unknown matrix features')
    print('a total of ' + str(no_ft_count) + ' fts contain all 0 information in the previous training matrix')

    with open(input_output_dir + '/opt_final_testing_mtx.csv', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


remove_additional_features(input_unknown_mtx_fl, input_pre_exp_fl,input_output_dir)










