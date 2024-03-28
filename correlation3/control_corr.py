# output four collinearity expression files
# import numpy as np

import base
import corr_analysis
import pandas as pd
# from argparse import ArgumentParser
# import sys


def row_corr(parameter):
    tassel = parameter.tassel
    ear = parameter.ear
    sorghum = parameter.sorghum

    collinearity = parameter.AnchorWave
    wgdi_collinearity = parameter.WGDI
    mcscanx_collinearity = parameter.MCScanX
    jcvi_anchors = parameter.jcvi

    corr = parameter.corr
    method = parameter.method
    bool_log = parameter.bool_log

    df_list = []
    tassel_ref_df, sor_query_df1 = base.read_expression_matrix(tassel, sorghum, bool_log)
    ear_ref_df, sor_query_df1 = base.read_expression_matrix(ear, sorghum, bool_log)

    if collinearity != "":
        query_to_ref, _ = base.read_collinearity(collinearity)
        sor_query_tassel, sor_ref_tassel = base.ref_equal_query_line(sor_query_df1, tassel_ref_df, query_to_ref)
        sor_query_ear, sor_ref_ear = base.ref_equal_query_line(sor_query_df1, ear_ref_df, query_to_ref)

        tassel_sor_anchorwave = corr_analysis.row_corr_analysis(sor_query_tassel, sor_ref_tassel, method)
        ear_sor_anchorwave = corr_analysis.row_corr_analysis(sor_query_ear, sor_ref_ear, method)

        corr_df_anchorwave_1 = pd.DataFrame(list(tassel_sor_anchorwave.items()), columns=["tassel_gene_name1", "tassel_corr_value1"]).drop_duplicates()
        corr_df_anchorwave_2 = pd.DataFrame(list(ear_sor_anchorwave.items()), columns=["ear_gene_name1", "ear_corr_value1"]).drop_duplicates()
        # corr_df = pd.concat([corr_df_anchorwave_1, corr_df_anchorwave_2], axis=1, ignore_index=False)
        # df_no_na_1 = corr_df.replace({pd.NA: ''})
        df_list.extend([corr_df_anchorwave_1, corr_df_anchorwave_2])
    if wgdi_collinearity != "":
        query_to_ref_wgdi, _ = base.read_wgdi_collinearity(wgdi_collinearity)
        sor_query_tassel_wgdi, sor_ref_tassel_wgdi = base.ref_equal_query_line(sor_query_df1, tassel_ref_df, query_to_ref_wgdi)
        sor_query_ear_wgdi, sor_ref_ear_wgdi = base.ref_equal_query_line(sor_query_df1, ear_ref_df, query_to_ref_wgdi)

        tassel_sor_wgdi = corr_analysis.row_corr_analysis(sor_query_tassel_wgdi, sor_ref_tassel_wgdi, method)
        ear_sor_wgdi = corr_analysis.row_corr_analysis(sor_query_ear_wgdi, sor_ref_ear_wgdi, method)

        corr_df_wgdi_1 = pd.DataFrame(list(tassel_sor_wgdi.items()), columns=["tassel_gene_name2", "tassel_corr_value2"]).drop_duplicates()
        corr_df_wgdi_2 = pd.DataFrame(list(ear_sor_wgdi.items()), columns=["ear_gene_name2", "ear_corr_value2"]).drop_duplicates()
        # corr_df = pd.concat([corr_df_wgdi_1, corr_df_wgdi_2], axis=1, ignore_index=False)
        # df_no_na_2 = corr_df.replace({pd.NA: ''})
        df_list.extend([corr_df_wgdi_1, corr_df_wgdi_2])
    if mcscanx_collinearity != "":
        query_to_ref_mcscan, _ = base.read_mc_scanx_collinearity(mcscanx_collinearity)
        sor_query_tassel_mcscan, sor_ref_tassel_mcscan = base.ref_equal_query_line(sor_query_df1, tassel_ref_df, query_to_ref_mcscan)
        sor_query_ear_mcscan, sor_ref_ear_mcscan = base.ref_equal_query_line(sor_query_df1, ear_ref_df, query_to_ref_mcscan)

        tassel_sor_mcscan = corr_analysis.row_corr_analysis(sor_query_tassel_mcscan, sor_ref_tassel_mcscan, method)
        ear_sor_mcscan = corr_analysis.row_corr_analysis(sor_query_ear_mcscan, sor_ref_ear_mcscan, method)

        corr_df_mcscan_1 = pd.DataFrame(list(tassel_sor_mcscan.items()), columns=["tassel_gene_name3", "tassel_corr_value3"]).drop_duplicates()
        corr_df_mcscan_2 = pd.DataFrame(list(ear_sor_mcscan.items()), columns=["ear_gene_name3", "ear_corr_value3"]).drop_duplicates()
        # corr_df = pd.concat([corr_df_mcscan_1, corr_df_mcscan_2], axis=1, ignore_index=False)
        # df_no_na_3 = corr_df.replace({pd.NA: ''})
        df_list.extend([corr_df_mcscan_1, corr_df_mcscan_2])
    if jcvi_anchors != "":
        query_to_ref_jcvi, _ = base.read_jcvi_collinearity(jcvi_anchors)
        sor_query_tassel_jcvi, sor_ref_tassel_jcvi = base.ref_equal_query_line(sor_query_df1, tassel_ref_df, query_to_ref_jcvi)
        sor_query_ear_jcvi, sor_ref_ear_jcvi = base.ref_equal_query_line(sor_query_df1, ear_ref_df, query_to_ref_jcvi)

        tassel_sor_jcvi = corr_analysis.row_corr_analysis(sor_query_tassel_jcvi, sor_ref_tassel_jcvi, method)
        ear_sor_jcvi = corr_analysis.row_corr_analysis(sor_query_ear_jcvi, sor_ref_ear_jcvi, method)

        corr_df_jcvi_1 = pd.DataFrame(list(tassel_sor_jcvi.items()), columns=["tassel_gene_name4", "tassel_corr_value4"]).drop_duplicates()
        corr_df_jcvi_2 = pd.DataFrame(list(ear_sor_jcvi.items()), columns=["ear_gene_name4", "ear_corr_value4"]).drop_duplicates()
        # corr_df = pd.concat([corr_df_jcvi_1, corr_df_jcvi_2], axis=1, ignore_index=False)
        # df_no_na_1 = corr_df.replace({pd.NA: ''})
        df_list.extend([corr_df_jcvi_1, corr_df_jcvi_2])

    # output row corr
    # corr_df_1 = pd.DataFrame(list(dict1.items()), columns=["tassel_gene_name", "corr_value1"])
    # corr_df_2 = pd.DataFrame(list(dict2.items()), columns=["ear_gene_name", "corr_value2"])
    corr_df = pd.concat(df_list, axis=1, ignore_index=False)
    # df_no_na = corr_df.replace({pd.NA: '', np.nan: ''})
    df_no_na = corr_df.fillna('')
    df_no_na.to_csv(corr, header=True, index=False)
