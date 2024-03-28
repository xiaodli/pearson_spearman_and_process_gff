import read_file
import corr_analysis
import pandas as pd


def row_corr(parameter):
    expression_file = parameter.expression
    collinearity = parameter.AnchorWave
    wgdi_collinearity = parameter.WGDI
    mcscanx_collinearity = parameter.MCScanX
    jcvi_anchors = parameter.jcvi
    corr = parameter.corr
    method = parameter.method
    bool_log = parameter.bool_log
    # read four files
    query_df, ref_df = read_file.read_expression_matrix_xlsx(expression_file, bool_log)

    df_list = []
    if collinearity != "":
        query_to_ref, _ = read_file.read_anchor_wave_collinearity(collinearity)
        query_df_new, ref_df_new = read_file.ref_equal_query_line(query_df, ref_df, query_to_ref)
        dict_corr = corr_analysis.row_corr_analysis(query_df_new, ref_df_new, method)
        corr_df1 = pd.DataFrame(list(dict_corr.items()), columns=["Gene_name1", "Corr_Value1"])
        corr_df1 = corr_df1.dropna(subset=['Corr_Value1'])
        corr_df1.drop_duplicates()
        df_list.append(corr_df1)
    if wgdi_collinearity != "":
        query_to_ref_wgdi, _ = read_file.read_wgdi_collinearity(wgdi_collinearity)
        query_df_new_wgdi, ref_df_new_wgdi = read_file.ref_equal_query_line(query_df, ref_df, query_to_ref_wgdi)
        dict_corr_wgdi = corr_analysis.row_corr_analysis(query_df_new_wgdi, ref_df_new_wgdi, method)
        corr_df2 = pd.DataFrame(list(dict_corr_wgdi.items()), columns=["Gene_name2", "Corr_Value2"])
        corr_df2 = corr_df2.dropna(subset=["Corr_Value2"])
        corr_df2.drop_duplicates()
        df_list.append(corr_df2)
    if mcscanx_collinearity != "":
        query_to_ref_mcsc, _ = read_file.read_mc_scanx_collinearity(mcscanx_collinearity)
        query_df_new_mcsc, ref_df_new_mcsc = read_file.ref_equal_query_line(query_df, ref_df, query_to_ref_mcsc)
        dict_corr_mcsc = corr_analysis.row_corr_analysis(query_df_new_mcsc, ref_df_new_mcsc, method)
        corr_df3 = pd.DataFrame(list(dict_corr_mcsc.items()), columns=["Gene_name3", "Corr_Value3"])
        corr_df3 = corr_df3.dropna(subset=["Corr_Value3"])
        corr_df3.drop_duplicates()
        df_list.append(corr_df3)
    if jcvi_anchors != "":
        query_to_ref_jcvi, _ = read_file.read_jcvi_collinearity(jcvi_anchors)
        query_df_new_jcvi, ref_df_new_jcvi = read_file.ref_equal_query_line(query_df, ref_df, query_to_ref_jcvi)
        dict_corr_jcvi = corr_analysis.row_corr_analysis(query_df_new_jcvi, ref_df_new_jcvi, method)
        corr_df4 = pd.DataFrame(list(dict_corr_jcvi.items()), columns=["Gene_name4", "Corr_Value4"])
        corr_df4 = corr_df4.dropna(subset=["Corr_Value4"])
        corr_df4.drop_duplicates()
        df_list.append(corr_df4)

    corr_df = pd.concat(df_list, axis=1, ignore_index=False)
    df_no_na = corr_df.fillna('')
    # df_no_na = corr_df.replace({pd.NA: '', np.nan: ''})
    df_no_na.to_csv(corr, header=True, index=False)


def block_total_corr(parameter):
    expression_file = parameter.expression
    collinearity = parameter.AnchorWave
    wgdi_collinearity = parameter.WGDI
    mcscanx_collinearity = parameter.MCScanX
    jcvi_anchors = parameter.jcvi
    corr = parameter.corr
    method = parameter.method
    column = parameter.column
    bool_log = parameter.bool_log
    length = parameter.Block_min
    is_random = parameter.is_random
    two_block = parameter.two_block

    df_list = []
    query_df, ref_df = read_file.read_expression_matrix_xlsx(expression_file, bool_log)
    ref_df_columns = ref_df.columns.copy()
    if collinearity != "":
        query_to_ref, block_len1 = read_file.read_anchor_wave_collinearity(collinearity)
        query_df_new, ref_df_new = read_file.total_block(query_df, ref_df, query_to_ref, column, ref_df_columns, is_random)
        dict_corr = corr_analysis.total_block_corr(query_df_new, ref_df_new, block_len1, method, length, two_block)
        corr_df_1 = pd.DataFrame(list(dict_corr.items()), columns=["Index1", "Corr_Value1"])
        df_list.append(corr_df_1)
    if wgdi_collinearity != "":
        query_to_ref_wgdi, block_len2 = read_file.read_wgdi_collinearity(wgdi_collinearity)
        query_df_new_wgdi, ref_df_new_wgdi = read_file.total_block(query_df, ref_df, query_to_ref_wgdi, column, ref_df_columns, is_random)
        dict_corr_wgdi = corr_analysis.total_block_corr(query_df_new_wgdi, ref_df_new_wgdi, block_len2, method, length, two_block)
        corr_df_2 = pd.DataFrame(list(dict_corr_wgdi.items()), columns=["Index2", "Corr_Value2"])
        df_list.append(corr_df_2)
    if mcscanx_collinearity != "":
        query_to_ref_mcsc, block_len3 = read_file.read_mc_scanx_collinearity(mcscanx_collinearity)
        query_df_new_mcsc, ref_df_new_mcsc = read_file.total_block(query_df, ref_df, query_to_ref_mcsc, column, ref_df_columns, is_random)
        dict_corr_mcsc = corr_analysis.total_block_corr(query_df_new_mcsc, ref_df_new_mcsc, block_len3, method, length, two_block)
        corr_df_3 = pd.DataFrame(list(dict_corr_mcsc.items()), columns=["Index3", "Corr_Value3"])
        df_list.append(corr_df_3)
    if jcvi_anchors != "":
        query_to_ref_jcvi, block_len4 = read_file.read_jcvi_collinearity(jcvi_anchors)
        query_df_new_jcvi, ref_df_new_jcvi = read_file.total_block(query_df, ref_df, query_to_ref_jcvi, column, ref_df_columns, is_random)
        dict_corr_jcvi = corr_analysis.total_block_corr(query_df_new_jcvi, ref_df_new_jcvi, block_len4, method, length, two_block)
        corr_df_4 = pd.DataFrame(list(dict_corr_jcvi.items()), columns=["Index4", "Corr_Value4"])
        df_list.append(corr_df_4)
    corr_df = pd.concat(df_list, axis=1, ignore_index=False)
    # df_no_na = corr_df.replace({pd.NA: ''})
    df_no_na = corr_df.fillna('')
    df_no_na.to_csv(corr, header=True, index=False)
