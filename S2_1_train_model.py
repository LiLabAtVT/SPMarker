#!/usr/bin/env python

##updating 020521 add the output comparing list file
##updating 020221 we need to add the auprc and auroc
##this script will train SVM and RF model and identify the markers
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from pandas import read_csv
from sklearn import svm
import pickle
from sklearn.metrics import classification_report
from sklearn.metrics import multilabel_confusion_matrix
import math
import os
import pandas as pd
import random
import numpy as np
from itertools import islice
import sys
import subprocess
import glob
import re

##updating 020221 add the auprc and auroc
from keras.utils import np_utils
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import auc
from sklearn.metrics import roc_curve, auc


input_S2_split_cross_val_o_dir = sys.argv[1] ##03_split_dataset_to_train_cv_indetest_013021/output_dir/step2_split_train_cross_val_opt_dir/output_dir
evaluation_score = sys.argv[2] ##default: All. We can also use the MCC
input_output_dir = sys.argv[3]
input_working_dir = sys.argv[4]
##updating 033022
step1_split_train_indetest_opt_dir = sys.argv[5]

##updating 020521 add a function to write out cells
def writeout_compare_file (labels_pred,labels_test,meta_test,id_nm_dic,input_opt_dir,testing_type):

    ##testing_type is SVMtest or SVMindetest RFtest RFindetest

    #################
    ##updating 020521
    pred_nm_list = []
    for eachpred_id in labels_pred[0]:
        #print(eachpred_id)
        pred_nm = id_nm_dic[str(eachpred_id)]
        pred_nm_list.append(pred_nm)

    ori_nm_list = []
    for eachori_id in labels_test:
        ori_nm = id_nm_dic[str(eachori_id)]
        ori_nm_list.append(ori_nm)

    store_final_line_list = []
    list_count = -1
    for eachcellnm in meta_test.index:
        list_count += 1
        pred_celltypenm = pred_nm_list[list_count]
        ori_celltypenm = ori_nm_list[list_count]
        final_line = eachcellnm + ' \t' + ori_celltypenm + '\t' + pred_celltypenm
        store_final_line_list.append(final_line)

    with open (input_opt_dir + '/opt_compare_prediction_' + testing_type + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



def mean_average_precision(test_labels, pred_labels, pred_decision_val):
    order = np.argsort(-pred_decision_val)
    test_labels_reord = test_labels[order]
    pred_labels_reord = pred_labels[order]
    precision_list = np.zeros((test_labels.shape[0]))
    for i in range(test_labels.shape[0]):
        precision_list[i] = np.sum(test_labels_reord[:i + 1] == pred_labels_reord[:i + 1]) / np.float32(i + 1)

    return np.mean(precision_list)

# Compute cell-type-specific MAP
def get_cell_type_map(n_cell_types, test_labels, pred_labels, pred_decision_val):
    cell_type_map = np.zeros((n_cell_types))
    for i, label in enumerate(np.unique(test_labels)):
        mask = (test_labels == label)
        cell_type_map[i] = mean_average_precision(test_labels[mask], pred_labels[mask], pred_decision_val[mask])
    return (cell_type_map)

##define a function to evaluate
def evaluation (labels_pred,prob_pred,labels_test,label_names,test_type,id_nm_dic,model_name,input_opt_dir,prob_pred_all):

    ##model_name is svm or rf

    ##test_type is independent or validate

    ##evaluate predicted cell types for svm, rf and knn
    acc = []
    overall_map = []
    cell_type_map = []
    for y_pred, y_prob_pred, y_test in zip(labels_pred, prob_pred, np.tile(labels_test,(3,1))):
        acc.append(accuracy_score(y_test, y_pred))
        overall_map.append(mean_average_precision(y_test, y_pred, y_prob_pred))
        cell_type_map.append(get_cell_type_map(label_names.shape[0], y_test, y_pred, y_prob_pred))

    ###########################
    ##write map results for knn
    ##initiate a dic to store method id and its name
    store_method_id_nm = {'1':model_name}
    meth_id = 0

    for eachpred_method in cell_type_map:
        meth_id += 1
        ##write name of method
        meth_nm = store_method_id_nm[str(meth_id)]
        single_cell_id = -1
        with open(input_opt_dir + '/opt_' + meth_nm + '_' + test_type + '_results_multiple_class.txt', 'w+') as opt:
            for eachvalue in eachpred_method:
                single_cell_id += 1
                single_cell_name = id_nm_dic[str(single_cell_id)]
                final_line = str(single_cell_id) + '\t' + single_cell_name + '\t' + str(eachvalue)
                opt.write(final_line + '\n')

    ################################
    ##write other evaluation reports
    target_names = []
    for i in range(len(list(id_nm_dic.keys()))):
        target_names.append(id_nm_dic[str(i)])
    with open(input_opt_dir + '/opt_' + model_name + '_' + test_type + '_class_report.txt', 'w+') as opt:
        opt.write(classification_report(labels_test, labels_pred[0],
                                        target_names=target_names) + '\n')  ##labels_pred is a list

    ##calculate TP TN FP FN
    mcm = multilabel_confusion_matrix(labels_test, labels_pred[0])
    tn = mcm[:, 0, 0]
    tp = mcm[:, 1, 1]
    fn = mcm[:, 1, 0]
    fp = mcm[:, 0, 1]

    SE = tp / (tp + fn)
    AC = (tp + tn) / (tp + tn + fp + fn)
    PR = tp / (tp + fp)
    SP = tn / (fp + tn)
    ##MCC and GM
    MCC_list = []
    GM_list = []
    for i in range(len(AC)):
        s_tp = tp[i]  ##single tp
        s_tn = tn[i]
        s_fp = fp[i]
        s_fn = fn[i]

        if (s_tp + s_fp) != 0 and (s_tp + s_fn) != 0 and (s_tn + s_fp) != 0 and (s_tn + s_fn) != 0:
            MCC = (s_tp * s_tn - s_fp * s_fn) / math.sqrt((s_tp + s_fp) * (s_tp + s_fn) * (s_tn + s_fp) * (s_tn + s_fn))
        else:
            MCC = 0

        if (s_tp + s_fn) != 0 and (s_fp + s_tn) != 0:
            GM = math.sqrt((s_tp / (s_tp + s_fn)) * (s_tn / (s_fp + s_tn)))
        else:
            GM = 0
        MCC_list.append(str(MCC))
        GM_list.append(str(GM))

    ##write out results
    ##all the evaluaton scores are array
    store_final_line_list = []
    ##generate the first line
    first_line = 'Cell_type' + '\t' + 'AC' + '\t' + 'SE' + '\t' + 'PR' + '\t' + \
                 'SP' + '\t' + 'MCC' + '\t' + 'GM' + '\t' + 'MAP'
    store_final_line_list.append(first_line)

    for i in range(len(list(id_nm_dic.keys()))):
        final_line = id_nm_dic[str(i)] + '\t' + str(AC[i]) + '\t' + str(SE[i]) + '\t' + \
                     str(PR[i]) + '\t' + str(SP[i]) + '\t' + str(MCC_list[i]) + '\t' + \
                     str(GM_list[i]) + '\t' + str(cell_type_map[0][i])
        store_final_line_list.append(final_line)

    with open(input_opt_dir + '/opt_' + model_name + '_' + test_type + '_class_evaluation_score.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    #########################################
    ##updation 020221 add the prauc and roauc
    labels_test_one_hot = np_utils.to_categorical(labels_test,len(list(id_nm_dic.keys())))  ##has 10 classes
    ##generate the AuPRC for each class
    precision = {}
    recall = {}
    average_precision = {}
    PRauc_all_dic = {}
    n_classes = labels_test_one_hot.shape[1]
    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(labels_test_one_hot[:, i],
                                                            prob_pred_all[:, i])
        average_precision[i] = average_precision_score(labels_test_one_hot[:, i], prob_pred_all[:, i])
        PRauc_all_dic[i] = auc(recall[i], precision[i])
    ##store the PRauc
    model_nm = model_name
    #model_nm = 'svm'
    with open (input_opt_dir + '/opt_auPRC_' + model_nm + '.txt','w+') as opt:
        for eachclass in PRauc_all_dic:
            final_line = str(eachclass) + '\t' + str(id_nm_dic[str(eachclass)]) + '\t' + str(PRauc_all_dic[eachclass])
            opt.write(final_line + '\n')

    ##generate the AuROC for each class
    # Compute ROC curve and ROC area for each class
    fpr = {}
    tpr = {}
    roc_auc_dic = {}
    n_classes = int(labels_test_one_hot.shape[1])
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(labels_test_one_hot[:, i], prob_pred_all[:, i], drop_intermediate=False)
        roc_auc_dic[i] = auc(fpr[i], tpr[i])
    ##store the roc
    #model_nm = 'svm'
    model_nm = model_name
    with open(input_opt_dir + '/opt_auROC_' + model_nm + '.txt', 'w+') as opt:
        for eachclass in roc_auc_dic:
            final_line = str(eachclass) + '\t' + str(id_nm_dic[str(eachclass)]) + '\t' + str(roc_auc_dic[eachclass])
            opt.write(final_line + '\n')




def train_evaluation (input_train_file,input_test_file,input_meta_train_file,input_meta_test_file,
                      input_test_indep_file,input_meta_test_indep_file,input_opt_dir):

    ##run script
    ##load import file
    features_train = read_csv(input_train_file, index_col = 0).values.T
    features_test =  read_csv(input_test_file, index_col = 0).values.T
    meta_train = read_csv(input_meta_train_file,header = 0, index_col = 0)
    meta_test = read_csv(input_meta_test_file,header = 0, index_col = 0)

    label_names = meta_train.loc[:,"cell_type"].unique()
    names_to_id = { key:val for val,key in enumerate(label_names)}
    labels_train = meta_train.loc[:,"cell_type"].replace(names_to_id).values
    labels_test = meta_test.loc[:,"cell_type"].replace(names_to_id).values

    ##initiate a function to store single cell id and its real name
    ##change the order of id and name
    id_nm_dic = {}
    for eachnm in names_to_id:
        id_nm_dic[str(names_to_id[eachnm])] = eachnm

    #########
    ##For svm
    clf_svm = svm.SVC(kernel='linear',probability=True)
    clf_svm.fit(features_train, labels_train)
    ##save svm model
    pkl_filename = input_opt_dir + "/svm_model.pkl"
    with open(pkl_filename, 'wb') as file:
        pickle.dump(clf_svm, file)

    ##for the validation data
    clfs = [clf_svm]
    labels_pred = []
    prob_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_test)
        labels_pred.append(np.argsort(-prob_pred_all,axis = 1)[:,0])
        prob_pred.append(np.array([prob_pred_all[i,labels_pred[-1][i]] for i in range(labels_pred[-1].shape[0])]))
    ##evaluation

    ##updating 020221
    prob_pred_all = ''
    labels_pred = []
    prob_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_test)
        labels_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
        prob_pred.append(
            np.array([prob_pred_all[i, labels_pred[-1][i]] for i in range(labels_pred[-1].shape[0])]))

    ##updating 020521
    writeout_compare_file(labels_pred, labels_test, meta_test, id_nm_dic, input_opt_dir, 'SVMtest')

    evaluation (labels_pred,prob_pred,labels_test,label_names,'validate',id_nm_dic,'svm',input_opt_dir,prob_pred_all)

    ##for the independent data
    features_indep_test =  read_csv(input_test_indep_file, index_col = 0).values.T
    meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
    labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
    labels_indep_pred = []
    prob_indep_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_indep_test)
        labels_indep_pred.append(np.argsort(-prob_pred_all,axis = 1)[:,0])
        prob_indep_pred.append(np.array([prob_pred_all[i,labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))
    ##evaluation
    ##updating 020221
    prob_pred_all = ''
    labels_indep_pred = []
    prob_indep_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_indep_test)
        labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
        prob_indep_pred.append(
            np.array([prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

    ##updating 020521
    writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'SVMindetest')

    evaluation (labels_indep_pred,prob_indep_pred,labels_indep_test,label_names,'independent',id_nm_dic,'svm',input_opt_dir,prob_pred_all)

    ########
    ##For RF
    clf_rf = RandomForestClassifier(n_estimators = 500, n_jobs = 42).fit(features_train, labels_train)
    ##save rf model
    pkl_filename = input_opt_dir + "/rf_model.pkl"
    with open(pkl_filename, 'wb') as file:
        pickle.dump(clf_rf, file)

    ##for the validation data
    clfs = [clf_rf]
    labels_pred = []
    prob_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_test)
        labels_pred.append(np.argsort(-prob_pred_all,axis = 1)[:,0])
        prob_pred.append(np.array([prob_pred_all[i,labels_pred[-1][i]] for i in range(labels_pred[-1].shape[0])]))

    ##updating 020221
    prob_pred_all = ''
    labels_pred = []
    prob_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_test)
        labels_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
        prob_pred.append(
            np.array([prob_pred_all[i, labels_pred[-1][i]] for i in range(labels_pred[-1].shape[0])]))

    ##updating 020521
    writeout_compare_file(labels_pred, labels_test, meta_test, id_nm_dic, input_opt_dir, 'RFtest')

    ##evaluation
    evaluation (labels_pred,prob_pred,labels_test,label_names,'validate',id_nm_dic,'rf',input_opt_dir,prob_pred_all)

    ##for the independent data
    features_indep_test =  read_csv(input_test_indep_file, index_col = 0).values.T
    meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
    labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
    labels_indep_pred = []
    prob_indep_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_indep_test)
        labels_indep_pred.append(np.argsort(-prob_pred_all,axis = 1)[:,0])
        prob_indep_pred.append(np.array([prob_pred_all[i,labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

    ##updating 020221
    prob_pred_all = ''
    labels_indep_pred = []
    prob_indep_pred = []
    for clf in clfs:
        prob_pred_all = clf.predict_proba(features_indep_test)
        labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
        prob_indep_pred.append(
            np.array([prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

    ##updating 020521
    writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'RFindetest')

    ##evaluation
    evaluation (labels_indep_pred,prob_indep_pred,labels_indep_test,label_names,'independent',id_nm_dic,'rf',input_opt_dir,prob_pred_all)

def extract_output (input_opt_dir,output_dir,method_list):

    for eachmethod_nm in method_list:

        opt_method_dir = output_dir + '/' + eachmethod_nm
        if not os.path.exists(opt_method_dir):
            os.makedirs(opt_method_dir)

        ##generate two output dir testing and validation
        opt_testing_dir = opt_method_dir + '/testing'
        if not os.path.exists(opt_testing_dir):
            os.makedirs(opt_testing_dir)

        opt_validation_dir = opt_method_dir + '/validation'
        if not os.path.exists(opt_validation_dir):
            os.makedirs(opt_validation_dir)

        rep_dir_list = glob.glob(input_opt_dir + '/*')
        for eachrep_dir in rep_dir_list:
            mt = re.match('.+/(.+)',eachrep_dir)
            rep_nm = mt.group(1) ##rep_nm is a b c d e

            fl_list = glob.glob(eachrep_dir + '/*class_evaluation_score.txt')
            for eachfl in fl_list:

                mt = re.match('.+/(.+)',eachfl)
                fl_nm = mt.group(1)

                mt = re.match('opt_(.+)_(.+)_class_evaluation_score\.txt',fl_nm)
                method_nm = mt.group(1)
                test_type = mt.group(2)

                if method_nm == eachmethod_nm:
                    ##cp the fl to the output file
                    ##consider two situations
                    if test_type == 'independent':
                        cmd = 'cp ' + eachfl + ' ' + opt_testing_dir + '/opt_' + rep_nm + '_class_evaluation_score.txt'
                        print(cmd)
                        subprocess.call(cmd,shell=True)
                    if test_type == 'validate':
                        cmd = 'cp ' + eachfl + ' ' + opt_validation_dir + '/opt_' + rep_nm + '_class_evaluation_score.txt'
                        print(cmd)
                        subprocess.call(cmd,shell=True)


def select_best_model (input_opt_store_model_dir,final_output_dir,collect_opt_dir,eval_score_nm,model_nm,
                       step1_split_train_indetest_opt_dir,input_working_dir):


    ##updating 033022 decide whether we will use specific score to decide the best model
    if eval_score_nm != 'All':

        store_rep_nm_dic = {} ##key is the rep_nm and value is the eval score
        opt_eval_score_fl_list = glob.glob(collect_opt_dir + '/' + model_nm + '/testing/*')
        for eacheval_score_fl in opt_eval_score_fl_list:

            mt = re.match('.+/(.+)',eacheval_score_fl)
            fl_nm = mt.group(1)

            mt = re.match('opt_(.+)_class_evaluation_score\.txt',fl_nm)
            rep_nm = mt.group(1)

            dt = pd.read_csv(eacheval_score_fl,sep='\t')
            ##updation 060620 change the mean_eval_score to the sum_eval_score since MCC has negative value
            sum_eval_score = dt[eval_score_nm].sum(axis=0)

            store_rep_nm_dic[rep_nm] = float(sum_eval_score)

        largest_rep_nm = max(store_rep_nm_dic, key=store_rep_nm_dic.get)

        ##cp the best model the final_output_dir
        cmd = 'cp ' +  input_opt_store_model_dir + '/' + largest_rep_nm + '/' + model_nm + '_model.pkl ' + final_output_dir
        subprocess.call(cmd,shell=True)

    else:
        ##if the eval_score_nm == All
        ##we will use all the training dataset to train
        ##create a dir to store the new training information
        opt_use_all_train_data_dir = input_working_dir + '/opt_use_all_train_data_dir'
        if not os.path.exists(opt_use_all_train_data_dir):
            os.makedirs(opt_use_all_train_data_dir)

        ipt_exp_train_fl = step1_split_train_indetest_opt_dir + '/opt_exp_train.csv'
        ipt_meta_train_fl = step1_split_train_indetest_opt_dir + '/opt_meta_train.csv'
        ipt_exp_test_fl = step1_split_train_indetest_opt_dir + '/opt_exp_test.csv'
        ipt_meta_test_fl = step1_split_train_indetest_opt_dir + '/opt_meta_test.csv'

        ipt_inde_exp_test_fl = step1_split_train_indetest_opt_dir + '/opt_exp_indep_test.csv'
        ipt_inde_meta_test_fl = step1_split_train_indetest_opt_dir + '/opt_meta_indep_test.csv'

        train_evaluation(ipt_exp_train_fl, ipt_exp_test_fl, ipt_meta_train_fl, ipt_meta_test_fl,
                         ipt_inde_exp_test_fl, ipt_inde_meta_test_fl, opt_use_all_train_data_dir)

        ##collect the models to final_output_dir
        cmd = 'cp ' + opt_use_all_train_data_dir + '/' + model_nm + '_model.pkl ' + final_output_dir
        subprocess.call(cmd,shell=True)



##################
##train the models
#print('Step 4 Model training and evaluation')
def train_models (input_S2_split_cross_val_o_dir,input_output_dir):
    S3_train_model_eval_o_dir = input_output_dir + '/S3_train_model_eval_o_dir'
    if not os.path.exists(S3_train_model_eval_o_dir):
        os.makedirs(S3_train_model_eval_o_dir)

    ipt_dt_list = glob.glob(input_S2_split_cross_val_o_dir + '/*')
    for eachipt_dir in ipt_dt_list:
        mt = re.match('.+/(.+)',eachipt_dir)
        rep_nm = mt.group(1) ##rep_nm is a, b, c, d, e

        ##generate output dir
        rep_output_dir = S3_train_model_eval_o_dir + '/' + rep_nm
        if not os.path.exists(rep_output_dir):
            os.makedirs(rep_output_dir)

        ipt_exp_train_fl = eachipt_dir + '/opt_exp_train.csv'
        ipt_meta_train_fl = eachipt_dir + '/opt_meta_train.csv'
        ipt_exp_test_fl = eachipt_dir + '/opt_exp_test.csv'
        ipt_meta_test_fl = eachipt_dir + '/opt_meta_test.csv'

        ipt_inde_exp_test_fl = eachipt_dir + '/opt_exp_indep_test.csv'
        ipt_inde_meta_test_fl = eachipt_dir + '/opt_meta_indep_test.csv'

        ##updating we will cp the exp and meta to output dir
        cmd = 'cp ' + ipt_inde_exp_test_fl + ' ' + input_output_dir
        subprocess.call(cmd,shell=True)
        cmd = 'cp ' + ipt_inde_meta_test_fl + ' ' + input_output_dir
        subprocess.call(cmd,shell=True)

        ##run the training process
        ##run script
        ##train both models since we also generate the svm and rf or shap markers at the same time
        train_evaluation(ipt_exp_train_fl, ipt_exp_test_fl, ipt_meta_train_fl, ipt_meta_test_fl,
                         ipt_inde_exp_test_fl, ipt_inde_meta_test_fl, rep_output_dir)


###########################
##collect model performance
#print('Step 5 Collect model performance')
def collect_model_performance (input_output_dir,S3_train_model_eval_o_dir,input_working_dir,evaluation_score,step1_split_train_indetest_opt_dir):

    S4_collect_res_o_dir = input_output_dir + '/S4_collect_res_o_dir'
    if not os.path.exists(S4_collect_res_o_dir):
        os.makedirs(S4_collect_res_o_dir)

    method_list = ['svm', 'rf']
    extract_output (S3_train_model_eval_o_dir,S4_collect_res_o_dir,method_list)

    ########
    ##Step 6 model select for the rf models
    print('Step 6 select the best rf model and collect all svm models')
    ##the output contains two dirs: rf and svm
    ##rf is used to compare to select the best model for detecting the SHAP markers
    ##svm model is used to extract feature importance to generate SVM markers
    ##we need to compare performance of each model using the testing evaluation
    select_best_model(S3_train_model_eval_o_dir, input_output_dir, S4_collect_res_o_dir, evaluation_score,'rf',step1_split_train_indetest_opt_dir,input_working_dir)

    ##updating also select best for the svm
    select_best_model(S3_train_model_eval_o_dir, input_output_dir, S4_collect_res_o_dir, evaluation_score, 'svm',step1_split_train_indetest_opt_dir,input_working_dir)

    ##we need to cp all the svm models to the output_dir
    ##generate a output in the final output dir
    opt_store_svm_models_dir = input_working_dir + '/opt_store_svm_models_dir'
    if not os.path.exists(opt_store_svm_models_dir):
        os.makedirs(opt_store_svm_models_dir)

    opt_rep_dir_list = glob.glob(S3_train_model_eval_o_dir + '/*')
    for eachrep_dir in opt_rep_dir_list:

        mt = re.match('.+/(.+)',eachrep_dir)
        rep_nm = mt.group(1)

        opt_store_svm_model_rep_nm_dir = opt_store_svm_models_dir + '/' + rep_nm
        if not os.path.exists(opt_store_svm_model_rep_nm_dir):
            os.makedirs(opt_store_svm_model_rep_nm_dir)

        cmd = 'cp ' + eachrep_dir + '/svm_model.pkl ' + opt_store_svm_model_rep_nm_dir
        subprocess.call(cmd,shell=True)

    ##the opt_store_svm_models_dir will be used for the SVM marker detection


train_models (input_S2_split_cross_val_o_dir,input_output_dir)
collect_model_performance (input_output_dir,input_output_dir + '/S3_train_model_eval_o_dir',input_working_dir,evaluation_score,step1_split_train_indetest_opt_dir)