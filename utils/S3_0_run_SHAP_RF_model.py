#!/usr/bin/env python

##updation 060520 generate a function to run the process
##updation 041920 this script will only run shap for RF since RF and SVM have different scripts to run the shap
##updation 032120 use shap to call feature importance
##updation 032120 change the train to test and only for the svm
##updation 030820 change the test dataset to the train dataset when conducting permulation
##this script will directly load the model and run the permutation

##this script is to train the PCA SVM RF model using python script

##import modules
##utils is another script used for evaluation

import numpy as np
from pandas import read_csv
import pickle
import shap

def run_SHAP (input_rf_model_fl,input_test_file,input_meta_test_file):

    ##step 0 load the data
    #features_train = read_csv(input_train_file, index_col = 0).values.T
    features_test =  read_csv(input_test_file, index_col = 0).values.T
    #meta_train = read_csv(input_meta_train_file,header = 0, index_col = 0)
    meta_test = read_csv(input_meta_test_file,header = 0, index_col = 0)
    label_names = meta_test.loc[:,"cell_type"].unique()
    names_to_id = { key:val for val,key in enumerate(label_names)}
    #labels_train = meta_train.loc[:,"cell_type"].replace(names_to_id).values
    labels_test = meta_test.loc[:,"cell_type"].replace(names_to_id).values

    ##step 1 load the model
    with open (input_rf_model_fl,'rb') as file:
        rf_model = pickle.load(file)

    ##step 2 conduct the permutation
    exp_test = read_csv(input_test_file, index_col = 0)
    feature_name_dic = {} ##key is string number from 0 and value is the feature name
    for i in range(len(exp_test.index)):
        feature_name_dic[str(i)] = exp_test.index[i]

    ##for the randomforest
    explainer = shap.TreeExplainer(rf_model)
    shap_values = explainer.shap_values(features_test,check_additivity=False) ##use test data
    ##store cell number and its relative cell type name
    label_names = meta_test.loc[:,"cell_type"].unique()
    names_to_id = { key:val for val,key in enumerate(label_names)}
    id_nm_dic = {}
    for eachnm in names_to_id:
        id_nm_dic[str(names_to_id[eachnm])] = eachnm
    ##run shap
    class_num = len(label_names)
    store_shap_value_dic = {}
    for i in range(class_num):
        class_nm = id_nm_dic[str(i)]  ##class name could be called by single_cell_name = id_nm_dic[str(single_cell_id)]
        store_shap_value_dic[class_nm] = {}
        shap_class_array = np.absolute(shap_values[i]).mean(axis=0)
        feature_id = -1
        for eachshap_val in shap_class_array:
            feature_id += 1
            feature_nm = feature_name_dic[str(feature_id)]
            store_shap_value_dic[class_nm][feature_nm] = eachshap_val


    return (store_shap_value_dic)



