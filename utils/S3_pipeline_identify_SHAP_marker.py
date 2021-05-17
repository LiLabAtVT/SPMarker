#!/usr/bin/env python

##Step 3 identification of markers from rf and svm

##BUILT-IN MODULES
import os
import argparse
import sys

import S3_0_run_SHAP_RF_model as S3_0_run_SHAP
import S3_1_generate_sort_SHAP_val as S3_1_sort_SHAP
import S3_2_assign_celltype_to_feature as S3_2_assign_feature

def get_parsed_args():

    parser = argparse.ArgumentParser(description="SPmarker identify markers")

    ##require files
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store intermediate files of "
                                                                     "each step. Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    parser.add_argument("-m", dest='RF_model_fl',help="Provide RF model file that is used to detect feature importance from the SHAP package."
                                                      "Users can find this model from the output directory of Step 2.")

    parser.add_argument("-exp_fl", dest="test_exp_file", help="Provide independent testing exp file from the output directory of Step 2.")

    parser.add_argument("-meta_fl", dest="test_meta_file", help="Provide independent testing meta file from the output directory of Step 2.")


    ##optional parameters
    parser.add_argument("-mar_num", dest="marker_number", help="Provide the candidate marker number users want to extract from each cell type."
                                                                "Default is 20."
                                                                "If the feature number is below 20, we will extract all the features under the cell type.")

    parser.add_argument("-kmar_fl", dest="known_marker_fl", help="Provide the known marker gene list file. Once users provide this file, "
                                                                 "they will obtain a file that contains novel marker genes.")

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
    if args.RF_model_fl is None:
        print('Cannot find model file, please provide it')
        return
    else:
        try:
            file = open(args.RF_model_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the model data!')
            return

    if args.test_exp_file is None:
        print('Cannot find testing expression file, please provide it')
        return
    else:
        try:
            file = open(args.test_exp_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the testing expression data!')
            return

    if args.test_meta_file is None:
        print('Cannot find testing meta file, please provide it')
        return
    else:
        try:
            file = open(args.test_meta_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the testing meta data!')
            return

    if args.known_marker_fl is not None:
        try:
            file = open(args.known_marker_fl ,'r')
        except IOError:
            print('There was an error opening the known marker file!')
            return


    ###parameters
    if args.marker_number is not None:
        marker_number = args.marker_number
    else:
        marker_number = '20'

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

    ########
    ##Step 0 run the SHAP package
    print('Step 0 Run the SHAP')
    RF_model_fl = args.RF_model_fl
    test_exp_file = args.test_exp_file
    test_meta_file = args.test_meta_file

    store_shap_value_dic = S3_0_run_SHAP.run_SHAP(RF_model_fl, test_exp_file, test_meta_file)

    with open(working_dir + '/opt_S0_shap_value_rf.txt', 'w+') as opt:
        for eachclass in store_shap_value_dic:
            for eachfeature in store_shap_value_dic[eachclass]:
                final_line = eachclass + '\t' + eachfeature + '\t' + str(store_shap_value_dic[eachclass][eachfeature])
                opt.write(final_line + '\n')

    ########
    ##Step 1 generate markers for each cell type
    print('Step 1 Generate markers for each cell type')
    ##generate a dir in the working_dir
    S1_sort_shap_cell_type_dir = working_dir + '/S1_sort_shap_cell_type_dir'
    if not os.path.exists(S1_sort_shap_cell_type_dir):
        os.makedirs(S1_sort_shap_cell_type_dir)

    S1_sort_shap_cell_type_w_dir = S1_sort_shap_cell_type_dir + '/S1_sort_shap_cell_type_w_dir'
    if not os.path.exists(S1_sort_shap_cell_type_w_dir):
        os.makedirs(S1_sort_shap_cell_type_w_dir)

    S1_sort_shap_cell_type_o_dir = S1_sort_shap_cell_type_dir + '/S1_sort_shap_cell_type_o_dir'
    if not os.path.exists(S1_sort_shap_cell_type_o_dir):
        os.makedirs(S1_sort_shap_cell_type_o_dir)

    input_shap_value_fl = working_dir + '/opt_S0_shap_value_rf.txt'

    S3_1_sort_SHAP.select_contribution (input_shap_value_fl,S1_sort_shap_cell_type_w_dir,S1_sort_shap_cell_type_o_dir)

    ########
    ##Step 2 assign the cell type to the cells
    print('Step 2 Assign cell type to features')
    ##this step will give several outputs
    ##1) top 20 novel markers that do not contain marker gene (one file, annotated with novel)
    ##2) top 20 markers that contain novel and known markers (one file, annoated with known or novel)
    ##3) summary of the top 20 marker count that includes known or novel markers
    S2_assign_cell_type_dir = working_dir + '/S2_assign_cell_type_dir'
    if not os.path.exists(S2_assign_cell_type_dir):
        os.makedirs(S2_assign_cell_type_dir)

    S2_assign_cell_type_w_dir = S2_assign_cell_type_dir + '/S2_assign_cell_type_w_dir'
    if not os.path.exists(S2_assign_cell_type_w_dir):
        os.makedirs(S2_assign_cell_type_w_dir)

    known_marker_fl = ''
    if args.known_marker_fl is not None:
        known_marker_fl = args.known_marker_fl

    store_final_assign_cell_type_line_list = S3_2_assign_feature.assign_cell_type_to_feature(S1_sort_shap_cell_type_o_dir)

    ##store the file a working_dir and divide the file into multiple files each files contains a cell type
    with open(S2_assign_cell_type_w_dir + '/temp_assign_cell_type_fl.txt', 'w+') as opt:
        for eachline in store_final_assign_cell_type_line_list:
            opt.write(eachline + '\n')

    ##generate a dir that contains the sorted value
    select_top_cell_type_sort_uni_fs_dir = S3_2_assign_feature.create_sort_shapvalue(S1_sort_shap_cell_type_o_dir,S2_assign_cell_type_w_dir,
                                                                                     S2_assign_cell_type_w_dir + '/temp_assign_cell_type_fl.txt')

    ##this case has known marker fl
    store_final_top_marker_no_known_line_list, \
    store_final_top_marker_with_known_line_list, \
    store_final_unis_top_feat_count_list = S3_2_assign_feature.generate_dir_to_store_sort_for_top_uni_fts(select_top_cell_type_sort_uni_fs_dir,
                                                                                                          marker_number,
                                                                                                          known_marker_fl)

    ##write out results to the output_dir
    with open (output_dir + '/opt_top_' + marker_number + '_novel_marker.txt','w+') as opt:
        for eachline in store_final_top_marker_no_known_line_list:
            opt.write(eachline + '\n')

    with open (output_dir + '/opt_top_' + marker_number + '_novel_known_marker.txt', 'w+') as opt:
        for eachline in store_final_top_marker_with_known_line_list:
            opt.write(eachline + '\n')

    with open(output_dir + '/opt_top_' + marker_number + '_summary_marker_composition.txt', 'w+') as opt:
        for eachline in store_final_unis_top_feat_count_list:
            opt.write(eachline + '\n')


if __name__ == "__main__":
    main()

















