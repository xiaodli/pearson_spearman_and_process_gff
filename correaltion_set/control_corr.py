import read_file_shuffle as read_file
import process_file
import pandas as pd
# import sys
import corr_analysis


def block_total_corr(parameter):
    expression_file = parameter.expression
    wgdi_collinearity = parameter.WGDI
    corr = parameter.corr
    method = parameter.method
    column = parameter.column
    bool_log = parameter.bool_log
    length = parameter.Block_min
    sheet_name = parameter.sheet_name
    drop_method = parameter.drop_method
    bootstrap = parameter.bootstrap
    size = parameter.sample_size

    df_list = []
    query_expr_df, ref_expr_df = read_file.read_expression_matrix_xlsx(expression_file, bool_log, sheet_name)
    ref_df_columns = ref_expr_df.columns.copy()
    query_expr_dict, ref_expr_dict = process_file.get_expr_dict(query_expr_df, ref_expr_df, column, ref_df_columns)

    if wgdi_collinearity != "":
        _, total_pre_df_core_record_list = read_file.read_wgdi_collinearity(wgdi_collinearity, length)
        if drop_method == "drop_all":
            total_pre_df_core_record_list = process_file.drop_duplicate_all(total_pre_df_core_record_list)
        if drop_method == "retain_one":
            total_pre_df_core_record_list = process_file.drop_duplicate_core_record(total_pre_df_core_record_list)
        same, inv = process_file.split_inv_normal_to_pair(total_pre_df_core_record_list)
        same_list_ele_tuple = process_file.get_corr_data(query_expr_dict, ref_expr_dict, same, bootstrap, size)
        same_dict_corr = corr_analysis.size_corr(same_list_ele_tuple, method)
        corr_df1 = pd.DataFrame(list(same_dict_corr.items()), columns=["Index1", "Corr_Value1"])
        df_list.append(corr_df1)

        inv_list_ele_tuple = process_file.get_corr_data(query_expr_dict, ref_expr_dict, inv, bootstrap, size)
        inv_dict_corr = corr_analysis.size_corr(inv_list_ele_tuple, method)
        corr_df2 = pd.DataFrame(list(inv_dict_corr.items()), columns=["Index2", "Corr_Value2"])
        df_list.append(corr_df2)

    corr_df = pd.concat(df_list, axis=1, ignore_index=False)
    # df_no_na = corr_df.replace({pd.NA: ''})
    df_no_na = corr_df.fillna('')
    df_no_na.to_csv(corr, header=True, index=False)


def three_types_sample_corr(parameter):
    expression_file = parameter.expression
    anchorwave1 = parameter.AnchorWave1
    anchorwave2 = parameter.AnchorWave2
    blast = parameter.blast
    corr = parameter.corr
    method = parameter.method
    column = parameter.column
    bool_log = parameter.bool_log
    length = parameter.Block_min
    sheet_name = parameter.sheet_name
    bootstrap = parameter.bootstrap
    size = parameter.sample_size

    df_list = []
    query_expr_df, ref_expr_df = read_file.read_expression_matrix_xlsx(expression_file, bool_log, sheet_name)
    ref_df_columns = ref_expr_df.columns.copy()
    query_expr_dict, ref_expr_dict = process_file.get_expr_dict(query_expr_df, ref_expr_df, column, ref_df_columns)

    if anchorwave1 != "" and anchorwave2 != "" and blast != "":
        recent_wgd_set, pre_wgd_set, homologous_pair_set = read_file.read_three_file(anchorwave1, anchorwave2, blast, length)

        recent_wgd_list_ele_tuple = process_file.get_corr_data(query_expr_dict, ref_expr_dict, recent_wgd_set, bootstrap, size)
        recent_wgd_dict_corr = corr_analysis.size_corr(recent_wgd_list_ele_tuple, method)
        corr_df1 = pd.DataFrame(list(recent_wgd_dict_corr.items()), columns=["Index3", "Corr_Value3"])
        df_list.append(corr_df1)

        pre_wgd_list_ele_tuple = process_file.get_corr_data(query_expr_dict, ref_expr_dict, pre_wgd_set, bootstrap, size)
        pre_wgd_dict_corr = corr_analysis.size_corr(pre_wgd_list_ele_tuple, method)
        corr_df2 = pd.DataFrame(list(pre_wgd_dict_corr.items()), columns=["Index4", "Corr_Value4"])
        df_list.append(corr_df2)

        homologous_list_ele_tuple = process_file.get_corr_data(query_expr_dict, ref_expr_dict, homologous_pair_set, bootstrap, size)
        homologous_dict_corr = corr_analysis.size_corr(homologous_list_ele_tuple, method)
        corr_df3 = pd.DataFrame(list(homologous_dict_corr.items()), columns=["Index5", "Corr_Value5"])
        df_list.append(corr_df3)
    corr_df = pd.concat(df_list, axis=1, ignore_index=False)
    # df_no_na = corr_df.replace({pd.NA: ''})
    df_no_na = corr_df.fillna('')
    df_no_na.to_csv(corr, header=True, index=False)
