#!/usr/bin/env python

##this script will generate average fts from the previous 45 combinations of feats per class
from pandas import read_csv
import sys
import glob
import re
import pickle
from statistics import mean

input_build_output_dir = sys.argv[1] ##04_train_model_balanced_nR5000_17kcellver_ninte_013021/working_dir/opt_store_svm_models_dir
input_five_cross_ipt_dir = sys.argv[2]  ##extract the train fl information
##03_split_dataset_to_train_cv_indetest_balanced_nR5000_17kcellver_ninte_013021/output_dir/step2_split_train_cross_val_opt_dir/output_dir/
class_num = sys.argv[3] ##we need to change as we want

input_output_dir = sys.argv[4]

#class_num = 10 ##we need to change as we want

def generate_imp (input_build_output_dir,input_five_cross_ipt_dir,input_output_dir,class_num):

    rep_let_dir_list = glob.glob(input_build_output_dir + '/*')
    for each_rep_let_dir in rep_let_dir_list:
        mt = re.match('.+/(.+)',each_rep_let_dir)
        rep_let = mt.group(1)
        svm_model_fl = each_rep_let_dir + '/svm_model.pkl'


        ##extract train file information
        train_fl = input_five_cross_ipt_dir + '/' + rep_let + '/opt_exp_test.csv'
        input_meta_train_file = input_five_cross_ipt_dir + '/' + rep_let + '/opt_meta_test.csv'

        ##extract label_names
        meta_train = read_csv(input_meta_train_file, header=0, index_col=0)
        label_names = meta_train.loc[:, "cell_type"].unique()

        ##extract feature importance from all the features
        store_final_weight_line_list = []
        ##generate the first line
        store_new_names_to_id_dic = {str(val): key for val, key in enumerate(label_names)}
        first_line = 'Gene'
        for i in range(0,int(class_num)):
            ##extract class name
            class_nm = store_new_names_to_id_dic[str(i)]
            first_line = first_line + '\t' + class_nm
        store_final_weight_line_list.append(first_line)


        ##store the gene name information
        exp_train = read_csv(train_fl, index_col=0)

        ##load the model
        with open(svm_model_fl, 'rb') as file:
            svm_load_model = pickle.load(file)

        ##extract the weights and calculate the average weights for each class
        ##1, store fts to each feature by constructing a list
        for i in range(0, svm_load_model.coef_.shape[1]):

            print('analyze ft count is ' + str(i))

            imp_list = []
            for j in range(0,svm_load_model.coef_.shape[0]):
                imp_list.append(abs(svm_load_model.coef_[j][i]))

            ##2, calculate the average of fts
            store_class_str_list = []
            for x in range(0, int(class_num)):
                for y in range(0, int(class_num)):
                    if x < y:
                        order_str = str(x) + str(y)
                        store_class_str_list.append(order_str)

            store_order_imp_dic = {}
            for p in range(0, int(int(class_num)*(int(class_num)-1)/2)): ##if class_mum is equal to 10 it will be 45
                store_order_imp_dic[store_class_str_list[p]] = imp_list[p] ##eg. store_order_imp_dic['01'] = 0.156

            ##3 calculate imp for each class
            ft_line = exp_train.index[i]
            ##use the this step to keep the same order of ft importance
            for n in range(0,int(class_num)):

                store_class_imp_list = [] ##store imp for each class
                for eachp in store_order_imp_dic: ##eachp is '01' '02' ..
                    if str(n) in eachp:
                        store_class_imp_list.append(float(store_order_imp_dic[eachp]))
                imp_class = mean(store_class_imp_list)
                ft_line = ft_line + '\t' + str(imp_class)

            store_final_weight_line_list.append(ft_line)

        ##write out final_weight_line_information
        with open(input_output_dir + '/opt_imp_svm_' + rep_let + '.txt', 'w+') as opt:
            for eachline in store_final_weight_line_list:
                opt.write(eachline + '\n')






generate_imp (input_build_output_dir,input_five_cross_ipt_dir,input_output_dir,class_num)