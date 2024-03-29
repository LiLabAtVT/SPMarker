#!/usr/bin/env python

##updating 070923 set an option to filter the cell types with less than a specified cell number

##this script will use the Rscript to
##conduct analysis on the Step 1 to 3.
##1. prepare data and generate meta file
##2. train models and make prediction on independent datasets
##3. identify markers

import os
import argparse
import sys
import subprocess
import re
import glob

def get_parsed_args():

    parser = argparse.ArgumentParser(description="SPmarker prepare data")

    ##require files
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store intermediate files of "
                                                                     "each step. Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    parser.add_argument('-mtx', dest='exp_matrix', help='Provide expression matrix. Rowname is gene, and column name is cell.'
                                                        'Please make sure the gene name do not contain space. '
                                                        'Otherwise, the gene name will be transfered to a name with "_" connected')

    ##optional parameters
    parser.add_argument('-mlist',dest='marker_list', help='Provide marker list with two columns seperated by space or tab.'
                                                           'First column is geneID,'
                                                           'Second column is celltype.')

    parser.add_argument('-m', dest='marker', help='Provide a marker that would help to define cell identity.')

    ##updating 070923
    parser.add_argument('-lowestCellN', dest='lowest_cell_num_per_celltype', default = '50',help='If users provide a marker list, '
                                                                                                 'please provide an additional cell number cutoff,'
                                                                                                 ' which enables the predicted cell type having number of cells above this cutoff.'
                                                                                                 'Default: 50')

    parser.add_argument('-meta',dest='meta',help= 'Provide a meta that contains known cell identity.'
                                                  'If the -meta is initiated, we should not provide -m')

    parser.add_argument('-ukn_mtx', dest='unknown_cell_fl', help='Provide unknown cell matrix file that is need to be assigned with cell type.')

    parser.add_argument('-feat_fl',dest="feature_file",help="Provide the features that will be kept in the expression file that is used for the training."
                                                            "If users do not provide the argument, we will use all the features.")

    ##other parameters
    parser.add_argument('-bns', dest='keep_balance', help='Balance the matrix of cell identities. This option works only -m is initiated.'
                                                          'Default: -bns no')

    parser.add_argument('-bns_ratio', dest='ratio_of_balance',default='1:1',help='Provide ratio of different cell identities.'
                                                                                 'If users set 1:1, and if number of marker labeled cells have less cells than the non-marker labeled cells,'
                                                                                 'it will sample same number of non-marker labeled cells as the marker labeled cells.'
                                                                                 'Default: 1:1. Left 1 is marker labeled cell identity')

    parser.add_argument('-cv_num',dest='cross_vali_num',help='Initiate x fold cross validation.'
                                                             'Default: 5')

    parser.add_argument('-indep_ratio',dest='indep_ratio',help='Provide ratio of independent dataset.'
                                                               'Default: 0.1')

    parser.add_argument('-eval_score',dest='type_eval_score',help='Bascially, we will use the all the training dataset to identify markers'
                                                                  'If users provide a specific evaluation score such as MCC in this argument, SPmarker will use this score to decide the best model that will be used for marker identification.'
                                                                  'Default: All')

    parser.add_argument("-mar_num", dest="marker_number", help="Provide the candidate marker number users want to extract from each cell type."
                                                                "Default is 20."
                                                                "If the feature number is below 20, we will extract all the features under the cell type.")

    ##updating 101221 this marker fl is replaced by -mlist
    #parser.add_argument("-kmar_fl", dest="known_marker_fl", help="Provide the known marker gene list file. Once users provide this file, "
    #                                                             "they will obtain a file that contains novel marker genes.")


    parser.add_argument('-SVM', dest='SVM_marker',help='Decide to generate the SVM markers.'
                                                       'Default: -SVM yes')

    ##updating 052121
    #parser.add_argument('-feat_fl',dest="feature_file",help="Provide the features that will be kept in the expression file that is used for the training."
    #                                                        "If users do not provide the argument, we will use all the features.")

    #parser.add_argument("-SPmarker_dir" ,dest="SPmarker_directory",help="Provide the path to the SPmarker_directory")

    #parser.add_argument("-merged_obj", dest="merged_object", help="Provide merged object generated from Seurat.")

    #parser.add_argument("-R_p", dest="R_path", help="Provide Rscript path."
    #                                                "Default: /usr/bin/Rscript.")

    ##Optional parameters
    #parser.add_argument("-kmar_fl", dest="known_marker_fl", help="Provide the known marker gene list file. Once users provide this file, "
    #                                                             "they will obtain a file that contains novel marker genes.")

    ##parse of parameters
    args = parser.parse_args()
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    #######################################
    ##check the required software and files
    ##for the input files
    ##the exp_matrix must be provided
    if args.exp_matrix is None:
        print('Cannot find expression matrix, please provide it')
        return
    else:
        try:
            file = open(args.exp_matrix, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the matrix file!')
            return

    ##three options to provide the meta file
    if args.meta is not None:
        print('A meta file is provided that contains cell annotation.')
        try:
            file = open(args.meta, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the matrix file!')
            return

        if args.marker_list is not None:

            print('A known marker list is provided, and SPmarker will return novel candidate markers')

            try:
                file = open(args.marker_list, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the matrix file!')
                return

            if args.marker is not None:
                print('Do not provide marker once marker list has been provided')
                return

    else:

        if args.marker_list is not None:

            print('A known marker list is provided that will be used to create a candidate meta file with cell annotation.')

            try:
                file = open(args.marker_list, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the matrix file!')
                return

            if args.marker is not None:
                print('Do not provide marker once marker list has been provided')
                return

        else:

            if args.marker is not None:
                print('Single marker is provided that will be used to create a candidate meta file.')

            else:
                print('Please provide meta, marker list or single marker information.')
                return

    ##updating 101221 set the unknown cell fl to be optional choice
    if args.unknown_cell_fl is not None:
        try:
            file = open(args.unknown_cell_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the matrix file!')
            return

    if args.ratio_of_balance is None:
        ratio_of_balance = '1:1'
    else:
        if not re.match('\d+:\d+',args.ratio_of_balance):
            print ('Please provide right format of ratio of balance')
            return
        else:
            ratio_of_balance = args.ratio_of_balance

    if args.cross_vali_num is None:
        cross_vali_num = '5'
    else:
        cross_vali_num = args.cross_vali_num

    if args.indep_ratio is None:
        indep_ratio = '0.1'
    else:
        indep_ratio = args.indep_ratio

    ##updating 033022 change the MCC to All
    if args.type_eval_score is None:
        type_eval_score = 'All'
    else:
        type_eval_score = args.type_eval_score

    ###parameters
    if args.marker_number is not None:
        marker_number = args.marker_number
    else:
        marker_number = '20'

    if args.SVM_marker is None:
        SVM_marker = 'yes'
    else:
        SVM_marker = args.SVM_marker

        if SVM_marker != 'yes' and SVM_marker != 'no':
            print ('Please use yes or no to open or close the identification of SVM markers')
            return

    if args.feature_file is not None:
        try:
            file = open(args.feature_file ,'r')
        except IOError:
            print('There was an error opening the feature_file!')
            return

    ###########################################
    ##create the working and output directories
    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    #################
    ##Run the process
    print ('Begin to run the whole process')
    print ('Step 1 prepare data')

    Step1_prepare_data_dir = working_dir + '/Step1_prepare_data_dir'
    if not os.path.exists(Step1_prepare_data_dir):
        os.makedirs(Step1_prepare_data_dir)

    ##obtain the path of utils
    run_script_path = __file__
    if '/' in run_script_path:
        mt = re.match('(.+)/.+',run_script_path)
        run_script_dir = mt.group(1)
        utils_dir = run_script_dir + '/utils'
    else:
        utils_dir = './utils'

    ############
    input_mtx_fl = args.exp_matrix
    ############

    ###############
    print('Check and change genes in matrix with space')
    store_final_line_list = []
    count = 0
    with open(input_mtx_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split(',')
                if ' ' in col[0]:
                    new_name = col[0].replace(' ', '_')

                    new_line = new_name
                    for i in range(1, len(col)):
                        new_line = new_line + ',' + col[i]
                    store_final_line_list.append(new_line)
                else:
                    store_final_line_list.append(eachline)
            else:
                store_final_line_list.append(eachline)

    with open(Step1_prepare_data_dir + '/temp_modified_gene_matrix.csv','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##now the new mtx is temp_modified_gene_matrix.csv
    input_mtx_fl = Step1_prepare_data_dir + '/temp_modified_gene_matrix.csv'


    ##start to prepare the training dataset
    if args.meta is not None:

        print ('Users choose to provide meta file')

        ############
        meta_fl_path = args.meta
        ############

    else:

        if args.marker_list is not None:
            print('A known marker list is provided that will be used to create a candidate meta file with cell annotation.')

            marker_list_fl = args.marker_list
            s1_0_use_markerlist_cell_identity_script = utils_dir + '/S1_0_use_markerlist_cell_identity.py'
            S1_0_use_markerlist_cell_identity_R_script = utils_dir + '/S1_0_use_markerlist_cell_identity.R'

            Step1_0_generate_meta_dir = Step1_prepare_data_dir + '/Step1_0_generate_meta_dir'
            if not os.path.exists(Step1_0_generate_meta_dir):
                os.makedirs(Step1_0_generate_meta_dir)
            Step1_0_generate_meta_o_dir = Step1_0_generate_meta_dir + '/output_dir'
            if not os.path.exists(Step1_0_generate_meta_o_dir):
                os.makedirs(Step1_0_generate_meta_o_dir)

            cmd = 'python ' + s1_0_use_markerlist_cell_identity_script + \
                  ' ' + input_mtx_fl + \
                  ' ' + Step1_0_generate_meta_o_dir + \
                  ' ' + marker_list_fl + \
                  ' ' + S1_0_use_markerlist_cell_identity_R_script
            subprocess.call(cmd,shell=True)

            ##put the file of meta fl here
            meta_fl_path = Step1_0_generate_meta_o_dir + '/opt_meta.csv'



            ##updating 070923
            print('Check if the cell number per cell type is above the cell number cutoff')
            if args.lowest_cell_num_per_celltype is not None:
                cellnum_cutoff = args.lowest_cell_num_per_celltype
                print('Users provide a lowest cell num cutoff: ' + cellnum_cutoff)
            else:
                cellnum_cutoff = '50'
                print('Users use the default lowest cell num cutoff: ' + '50')

            store_celltype_cellnum_dic = {}
            count = 0
            with open(meta_fl_path, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split(',')
                    count += 1
                    if count != 1:
                        celltype = col[1]
                        if celltype in store_celltype_cellnum_dic:
                            store_celltype_cellnum_dic[celltype] += 1
                        else:
                            store_celltype_cellnum_dic[celltype] = 1

            store_celltype_pass_cellnumcutoff_dic = {}
            for eachcelltype in store_celltype_cellnum_dic:
                cellnum = store_celltype_cellnum_dic[eachcelltype]
                if cellnum >= int(cellnum_cutoff):
                    store_celltype_pass_cellnumcutoff_dic[eachcelltype] = 1

            store_final_line_list = []
            store_target_cell_dic = {}
            count = 0
            with open(meta_fl_path, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split(',')
                    count += 1
                    if count != 1:
                        celltype = col[1]
                        cellnm = col[0]
                        if celltype in store_celltype_pass_cellnumcutoff_dic:
                            store_final_line_list.append(eachline)
                            store_target_cell_dic[cellnm] = 1
                    else:
                        store_final_line_list.append(eachline)

            with open(Step1_0_generate_meta_o_dir + '/opt_all_meta_flt_celltype.csv', 'w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            ##assign new meta fl path
            meta_fl_path = Step1_0_generate_meta_o_dir + '/opt_all_meta_flt_celltype.csv'

            ##allow the matrix has the same cell name
            store_final_line_list = []
            store_order_list = []
            count = 0
            with open(input_mtx_fl, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split(',')
                    count += 1
                    if count == 1:

                        store_new_header_sample_list = []
                        for i in range(1, len(col)):
                            if col[i] in store_target_cell_dic:
                                store_order_list.append(i)
                                store_new_header_sample_list.append(col[i])

                        final_line = ',' + ','.join(store_new_header_sample_list)
                        store_final_line_list.append(final_line)

                    else:

                        store_new_val_list = []
                        for eachorder in store_order_list:
                            val = col[eachorder]
                            store_new_val_list.append(val)

                        final_line = col[0] + ',' + ','.join(store_new_val_list)
                        store_final_line_list.append(final_line)

            with open(Step1_0_generate_meta_o_dir + '/temp_gene_cell_exp_mtx_flt_celltype.csv', 'w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            input_mtx_fl = Step1_0_generate_meta_o_dir + '/temp_gene_cell_exp_mtx_flt_celltype.csv'


        else:

            if args.marker is not None:
                print('Single marker is provided that will be used to create a candidate meta file.')

                target_marker = args.marker
                s1_use_marker_cell_idenity_script = utils_dir + '/S1_1_use_marker_cell_identity.py'

                Step1_1_generate_meta_dir = Step1_prepare_data_dir + '/Step1_1_generate_meta_dir'
                if not os.path.exists(Step1_1_generate_meta_dir):
                    os.makedirs(Step1_1_generate_meta_dir)
                Step1_1_generate_meta_o_dir = Step1_1_generate_meta_dir + '/output_dir'
                if not os.path.exists(Step1_1_generate_meta_o_dir):
                    os.makedirs(Step1_1_generate_meta_o_dir)

                cmd = 'python ' + s1_use_marker_cell_idenity_script + \
                      ' ' + input_mtx_fl + \
                      ' ' + Step1_1_generate_meta_o_dir + \
                      ' ' + target_marker
                subprocess.call(cmd, shell=True)
                # print ('The meta file has been created with three columns: cellnames,identity,probability')

                meta_fl_path = Step1_1_generate_meta_o_dir + '/opt_all_meta.csv'

                ##updating 070923
                print('Check if the cell number per cell type is above the cell number cutoff')
                if args.lowest_cell_num_per_celltype is not None:
                    cellnum_cutoff = args.lowest_cell_num_per_celltype
                    print('Users provide a lowest cell num cutoff: ' + cellnum_cutoff)
                else:
                    cellnum_cutoff = '50'
                    print('Users use the default lowest cell num cutoff: ' + '50')

                store_celltype_cellnum_dic = {}
                count = 0
                with open (meta_fl_path,'r') as ipt:
                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        col = eachline.strip().split(',')
                        count += 1
                        if count != 1:
                            celltype = col[1]
                            if celltype in store_celltype_cellnum_dic:
                                store_celltype_cellnum_dic[celltype] += 1
                            else:
                                store_celltype_cellnum_dic[celltype] = 1

                store_celltype_pass_cellnumcutoff_dic = {}
                for eachcelltype in store_celltype_cellnum_dic:
                    cellnum = store_celltype_cellnum_dic[eachcelltype]
                    if cellnum >= int(cellnum_cutoff):
                        store_celltype_pass_cellnumcutoff_dic[eachcelltype] = 1

                store_final_line_list = []
                store_target_cell_dic = {}
                count = 0
                with open (meta_fl_path,'r') as ipt:
                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        col = eachline.strip().split(',')
                        count += 1
                        if count != 1:
                            celltype = col[1]
                            cellnm = col[0]
                            if celltype in store_celltype_pass_cellnumcutoff_dic:
                                store_final_line_list.append(eachline)
                                store_target_cell_dic[cellnm] = 1
                        else:
                            store_final_line_list.append(eachline)

                with open (Step1_1_generate_meta_o_dir + '/opt_all_meta_flt_celltype.csv','w+') as opt:
                    for eachline in store_final_line_list:
                        opt.write(eachline + '\n')

                ##assign new meta fl path
                meta_fl_path = Step1_1_generate_meta_o_dir + '/opt_all_meta_flt_celltype.csv'

                ##allow the matrix has the same cell name
                store_final_line_list = []
                store_order_list = []
                count = 0
                with open (input_mtx_fl,'r') as ipt:
                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        col = eachline.strip().split(',')
                        count += 1
                        if count == 1:

                            store_new_header_sample_list = []
                            for i in range(1,len(col)):
                                if col[i] in store_target_cell_dic:
                                    store_order_list.append(i)
                                    store_new_header_sample_list.append(col[i])

                            final_line = ',' + ','.join(store_new_header_sample_list)
                            store_final_line_list.append(final_line)

                        else:

                            store_new_val_list = []
                            for eachorder in store_order_list:
                                val = col[eachorder]
                                store_new_val_list.append(val)

                            final_line = col[0] + ',' + ','.join(store_new_val_list)
                            store_final_line_list.append(final_line)

                with open (Step1_1_generate_meta_o_dir + '/temp_gene_cell_exp_mtx_flt_celltype.csv','w+') as opt:
                    for eachline in store_final_line_list:
                        opt.write(eachline + '\n')

                input_mtx_fl = Step1_1_generate_meta_o_dir + '/temp_gene_cell_exp_mtx_flt_celltype.csv'




                if args.keep_balance == 'yes':
                    print('Users choose to keep balance of cell identities of the meta file')

                    s1_keep_balance_of_meta_script = utils_dir + '/S1_2_keep_balance_of_meta.py'

                    Step1_2_keep_balance_of_meta_dir = Step1_prepare_data_dir + '/Step1_2_keep_balance_of_meta_dir'
                    if not os.path.exists(Step1_2_keep_balance_of_meta_dir):
                        os.makedirs(Step1_2_keep_balance_of_meta_dir)

                    Step1_2_keep_balance_of_meta_o_dir = Step1_2_keep_balance_of_meta_dir + '/output_dir'
                    if not os.path.exists(Step1_2_keep_balance_of_meta_o_dir):
                        os.makedirs(Step1_2_keep_balance_of_meta_o_dir)

                    cmd = 'python ' + s1_keep_balance_of_meta_script + \
                          ' ' + input_mtx_fl + \
                          ' ' + meta_fl_path + \
                          ' ' + ratio_of_balance + \
                          ' ' + Step1_2_keep_balance_of_meta_o_dir
                    subprocess.call(cmd, shell=True)

                    ##we need to update the meta_fl_path and input_mtx_fl
                    opt_fl_list = glob.glob(Step1_2_keep_balance_of_meta_o_dir + '/*')
                    for eachfl in opt_fl_list:
                        mt = re.match('.+/(.+)', eachfl)
                        flnm = mt.group(1)
                        if 'balance_meta' in flnm:
                            meta_fl_path = eachfl
                        if 'balance_exp' in flnm:
                            input_mtx_fl = eachfl

            else:
                print('Please provide meta, marker list or single marker information.')
                return

    ##udpating 052121
    ##check whether we will select a part of features to be training
    if args.feature_file is not None:

        print ('Users choose to use picked features to do the training')

        s1_select_feature_script = utils_dir + '/S1_2_select_feature.py'

        ipt_feature_file = args.feature_file
        ipt_expression_data = input_mtx_fl
        ##create a dir under the working_dir
        S1_2_select_feature_dir = working_dir + '/S1_2_select_feature_dir'
        if not os.path.exists(S1_2_select_feature_dir):
            os.makedirs(S1_2_select_feature_dir)

        cmd = 'python ' + s1_select_feature_script + \
              ' ' + ipt_expression_data + \
              ' ' + ipt_feature_file + \
              ' ' + S1_2_select_feature_dir
        subprocess.call(cmd,shell=True)

        input_mtx_fl = S1_2_select_feature_dir + '/opt_select_feat_exp.csv'


    ##splite the dataset
    s1_split_dataset_script = utils_dir + '/S1_3_split_dataset_to_train_cv_indetest.py'

    ##use the meta_fl_path to generate testing training and independent testing dataset
    Step1_3_split_data_dir = Step1_prepare_data_dir + '/Step1_3_split_data_dir'
    if not os.path.exists(Step1_3_split_data_dir):
        os.makedirs(Step1_3_split_data_dir)
    Step1_3_split_data_o_dir = Step1_3_split_data_dir + '/output_dir'
    if not os.path.exists(Step1_3_split_data_o_dir):
        os.makedirs(Step1_3_split_data_o_dir)

    cmd = 'python ' + s1_split_dataset_script + \
          ' ' + input_mtx_fl + \
          ' ' + meta_fl_path + \
          ' ' + cross_vali_num + \
          ' ' + indep_ratio + \
          ' ' + Step1_3_split_data_o_dir
    subprocess.call(cmd,shell=True)

    #####################
    ##Step 2 train models
    #####################
    print('Step 2 train models')

    s2_train_model_script = utils_dir + '/S2_1_train_model.py'

    Step2_train_models_dir = working_dir + '/Step2_train_models_dir'
    if not os.path.exists(Step2_train_models_dir):
        os.makedirs(Step2_train_models_dir)

    Step2_train_models_w_dir = Step2_train_models_dir + '/Step2_train_models_w_dir'
    if not os.path.exists(Step2_train_models_w_dir):
        os.makedirs(Step2_train_models_w_dir)

    Step2_train_models_o_dir = Step2_train_models_dir + '/Step2_train_models_o_dir'
    if not os.path.exists(Step2_train_models_o_dir):
        os.makedirs(Step2_train_models_o_dir)

    ipt_cross_dataset_dir = Step1_3_split_data_o_dir + '/step2_split_train_cross_val_opt_dir/output_dir'

    ##updating 033022 add s1 input data
    cmd = 'python ' + s2_train_model_script + \
          ' ' + ipt_cross_dataset_dir + \
          ' ' + type_eval_score + \
          ' ' + Step2_train_models_o_dir + \
          ' ' + Step2_train_models_w_dir + \
          ' ' + Step1_3_split_data_o_dir + '/step1_split_train_indetest_opt_dir'
    subprocess.call(cmd,shell=True)

    ##we need to generate an output to collect all the prediction from the independent datasets

    ##collect the model to the output dir

    #############################
    ##Step 3 identify SHAP marker
    #############################
    ##this step will identify SHAP and SVM markers at same time
    ##identify the SHAP markers
    print('Step 3 identify SHAP markers')

    s3_pipeline_identify_SHAP_marker_script = utils_dir + '/S3_pipeline_identify_SHAP_marker.py'

    Step3_identify_marker_dir = working_dir + '/Step3_identify_SHAP_marker_dir'
    if not os.path.exists(Step3_identify_marker_dir):
        os.makedirs(Step3_identify_marker_dir)

    Step3_identify_marker_w_dir = Step3_identify_marker_dir + '/Step3_identify_marker_w_dir'
    if not os.path.exists(Step3_identify_marker_w_dir):
        os.makedirs(Step3_identify_marker_w_dir)

    ##create a dir in the major output_dir to store the SHAP markers output
    opt_SHAP_markers_dir = output_dir + '/opt_SHAP_markers_dir'
    if not os.path.exists(opt_SHAP_markers_dir):
        os.makedirs(opt_SHAP_markers_dir)

    if args.marker_list is not None:
        known_marker_fl = args.marker_list
        cmd = 'python ' + s3_pipeline_identify_SHAP_marker_script + \
              ' -d ' + Step3_identify_marker_w_dir + \
              ' -o ' + opt_SHAP_markers_dir + \
              ' -m ' + Step2_train_models_o_dir + '/rf_model.pkl' + \
              ' -exp_fl ' + Step2_train_models_o_dir + '/opt_exp_indep_test.csv' + \
              ' -meta_fl ' + Step2_train_models_o_dir + '/opt_meta_indep_test.csv' + \
              ' -kmar_fl ' + known_marker_fl + \
              ' -mar_num ' + marker_number
        print(cmd)
        subprocess.call(cmd,shell=True)

    else:
        cmd = 'python ' + s3_pipeline_identify_SHAP_marker_script + \
              ' -d ' + Step3_identify_marker_w_dir + \
              ' -o ' + opt_SHAP_markers_dir + \
              ' -m ' + Step2_train_models_o_dir + '/rf_model.pkl' + \
              ' -exp_fl ' + Step2_train_models_o_dir + '/opt_exp_indep_test.csv' + \
              ' -meta_fl ' + Step2_train_models_o_dir + '/opt_meta_indep_test.csv' + \
              ' -mar_num ' + marker_number
        print(cmd)
        subprocess.call(cmd,shell=True)

    ############################
    ##Step 4 identify SVM marker
    ############################
    if SVM_marker == 'yes':

        print('Step 4 identify SVM markers')

        s4_generate_imp_for_built_SVM_script = utils_dir + '/S4_1_generate_imp_for_built_SVM.py'

        Step4_identify_SVM_marker_dir = working_dir + '/Step4_identify_SVM_marker_dir'
        if not os.path.exists(Step4_identify_SVM_marker_dir):
            os.makedirs(Step4_identify_SVM_marker_dir)

        S4_1_generate_imp_for_built_SVM = Step4_identify_SVM_marker_dir + '/S4_1_generate_imp_for_built_SVM'
        if not os.path.exists(S4_1_generate_imp_for_built_SVM):
            os.makedirs(S4_1_generate_imp_for_built_SVM)

        S4_1_generate_imp_for_built_o_SVM = S4_1_generate_imp_for_built_SVM + '/S4_1_generate_imp_for_built_o_SVM'
        if not os.path.exists(S4_1_generate_imp_for_built_o_SVM):
            os.makedirs(S4_1_generate_imp_for_built_o_SVM)

        ##extract number of cell type
        store_celltype_num_dic = {}
        count = 0
        with open (meta_fl_path,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                count += 1
                if count != 1:
                    col = eachline.strip().split(',')
                    store_celltype_num_dic[col[1]] = 1
        celltype_num = str(len(list(store_celltype_num_dic.keys())))

        ipt_cross_dataset_dir = Step1_3_split_data_o_dir + '/step2_split_train_cross_val_opt_dir/output_dir'

        cmd = 'python ' + s4_generate_imp_for_built_SVM_script + \
              ' ' + Step2_train_models_w_dir + '/opt_store_svm_models_dir' + \
              ' ' + ipt_cross_dataset_dir + \
              ' ' + celltype_num + \
              ' ' + S4_1_generate_imp_for_built_o_SVM
        subprocess.call(cmd,shell=True)

        ##generate SVM markers
        s4_generate_SVM_markers_script = utils_dir + '/S4_2_generate_SVM_markers.py'

        Step4_2_identify_SVM_marker_dir = Step4_identify_SVM_marker_dir + '/Step4_2_identify_SVM_marker_dir'
        if not os.path.exists(Step4_2_identify_SVM_marker_dir):
            os.makedirs(Step4_2_identify_SVM_marker_dir)

        Step4_2_identify_SVM_marker_w_dir = Step4_2_identify_SVM_marker_dir + '/working_dir'
        if not os.path.exists(Step4_2_identify_SVM_marker_w_dir):
            os.makedirs(Step4_2_identify_SVM_marker_w_dir)

        Step4_2_identify_SVM_marker_o_dir = Step4_2_identify_SVM_marker_dir + '/output_dir'
        if not os.path.exists(Step4_2_identify_SVM_marker_o_dir):
            os.makedirs(Step4_2_identify_SVM_marker_o_dir)

        ##create a dir to store the SVM markers
        opt_SVM_markers_dir = output_dir + '/opt_SVM_markers_dir'
        if not os.path.exists(opt_SVM_markers_dir):
            os.makedirs(opt_SVM_markers_dir)

        if args.marker_list is not None:
            known_marker_fl = args.marker_list

            cmd = 'python ' + s4_generate_SVM_markers_script + \
                  ' ' + S4_1_generate_imp_for_built_o_SVM + \
                  ' ' + known_marker_fl + \
                  ' ' + marker_number + \
                  ' ' + Step4_2_identify_SVM_marker_w_dir + \
                  ' ' + opt_SVM_markers_dir + \
                  ' ' + 'yes'
            subprocess.call(cmd,shell=True)

        else:
            cmd = 'python ' + s4_generate_SVM_markers_script + \
                  ' ' + S4_1_generate_imp_for_built_o_SVM + \
                  ' ' + 'no_known_marker_provided' + \
                  ' ' + marker_number + \
                  ' ' + Step4_2_identify_SVM_marker_w_dir + \
                  ' ' + opt_SVM_markers_dir + \
                  ' ' + 'no'
            subprocess.call(cmd,shell=True)

    ####################################
    ##Step 5 prediction of unknown cells
    ####################################
    if args.unknown_cell_fl is not None:

        print ('Step 5 start to predict identities of unknown cells')

        unknown_cell_fl = args.unknown_cell_fl

        s5_keep_same_feature_as_training = utils_dir + '/S5_1_keep_same_feature_as_training.py'

        Step5_predict_unknown_cells_dir = working_dir + '/Step5_predict_unknown_cells_dir'
        if not os.path.exists(Step5_predict_unknown_cells_dir):
            os.makedirs(Step5_predict_unknown_cells_dir)

        Step5_1_keep_same_feature_as_training_dir = Step5_predict_unknown_cells_dir + '/Step5_1_keep_same_feature_as_training_dir'
        if not os.path.exists(Step5_1_keep_same_feature_as_training_dir):
            os.makedirs(Step5_1_keep_same_feature_as_training_dir)

        Step5_1_keep_same_feature_as_training_o_dir = Step5_1_keep_same_feature_as_training_dir + '/Step5_1_keep_same_feature_as_training_o_dir'
        if not os.path.exists(Step5_1_keep_same_feature_as_training_o_dir):
            os.makedirs(Step5_1_keep_same_feature_as_training_o_dir)

        ##create a dir to store the prediction results
        opt_prediction_dir = output_dir + '/opt_prediction_dir'
        if not os.path.exists(opt_prediction_dir):
            os.makedirs(opt_prediction_dir)

        cmd = 'python ' + s5_keep_same_feature_as_training + \
              ' ' + unknown_cell_fl + \
              ' ' + Step1_3_split_data_o_dir + '/opt_exp_indep_test.csv' + \
              ' ' + Step5_1_keep_same_feature_as_training_o_dir
        subprocess.call(cmd,shell=True)

        ## we need to modify it
        ##make a prediction
        Step5_2_make_prediction_dir = Step5_predict_unknown_cells_dir + '/Step5_2_make_prediction_dir'
        if not os.path.exists(Step5_2_make_prediction_dir):
            os.makedirs(Step5_2_make_prediction_dir)

        #Step5_2_make_prediction_o_dir = Step5_2_make_prediction_dir + '/Step5_2_make_prediction_o_dir'
        #if not os.path.exists(Step5_2_make_prediction_o_dir):
        #    os.makedirs(Step5_2_make_prediction_o_dir)

        opt_exp_indep_test_path = Step5_1_keep_same_feature_as_training_o_dir + '/opt_final_testing_mtx.csv'
        opt_meta_train_path = ipt_cross_dataset_dir + '/a/opt_meta_train.csv'

        rf_model_path = Step2_train_models_o_dir + '/rf_model.pkl'
        svm_model_path = Step2_train_models_o_dir + '/svm_model.pkl'

        s5_make_prediction_script = utils_dir + '/S5_2_make_prediction.py'

        cmd = 'python ' + s5_make_prediction_script + ' ' + \
              opt_exp_indep_test_path + ' ' + \
              rf_model_path + ' ' + \
              svm_model_path + ' ' + \
              opt_meta_train_path + ' ' + \
              opt_prediction_dir
        subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()



