import pandas as pd
import numpy as np
from itertools import islice
import random


# read AnchorWave collinearity file
def read_anchor_wave_collinearity(file):
    block_len = []
    with open(file) as f:
        _ = next(f)
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[2].split(sep="N=")[1]))
    # read anchorWave output collinearity file as a dataframe, every line of the result occurs once different from wgdi.
    collinearity_df = pd.read_table(file, comment="#", header=0, low_memory=False)

    # get two gene lists and delete "gene-" string
    ref_gene_list = list(collinearity_df.loc[:, "refGene"])   # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.loc[:, "queryGene"])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    for i in range(len(query_gene_list)):
        query_gene_list[i] = query_gene_list[i].split(sep=":")[1]
    for i in range(len(ref_gene_list)):
        ref_gene_list[i] = ref_gene_list[i].split(sep=":")[1]
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    return query_to_ref, block_len


def read_wgdi_collinearity(file):
    block_len = []
    with open(file) as f:
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[5].split(sep="N=")[1]))
    # read wgdi output collinearity file as a dataframe.
    collinearity_df = pd.read_table(file, comment="#", header=None, sep=r'\s+', low_memory=False)
    collinearity_df.columns = ["queryGene", "order_1", "refGene", "order_2", "direction"]
    # get two gene lists and delete "gene:" string
    ref_gene_list = list(collinearity_df.loc[:, "refGene"])   # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.loc[:, "queryGene"])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    for i in range(len(query_gene_list)):
        query_gene_list[i] = query_gene_list[i].split(sep=":")[1]
    for i in range(len(ref_gene_list)):
        ref_gene_list[i] = ref_gene_list[i].split(sep=":")[1]
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    return query_to_ref, block_len


def read_mc_scanx_collinearity(file):
    block_len = []
    with open(file) as f:
        _ = next(islice(f, 10, 11))
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[5].split(sep="N=")[1]))
    # read MCScanX output collinearity file as a dataframe.
    collinearity_df = pd.read_table(file, comment="#", header=None, sep='\t', low_memory=False)
    # print(collinearity_df)
    # collinearity_df.columns = ["queryGene", "order_1", "refGene", "order_2", "direction"]
    # get two gene lists and delete "gene:" string
    ref_gene_list = list(collinearity_df.iloc[:, 2])   # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.iloc[:, 1])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    # print(query_to_ref)
    return query_to_ref, block_len


# mcscan(python version) + quota_align(python version)
def read_jcvi_collinearity(file):
    block_len = []
    length = 0
    with open(file) as f:
        for line in f:
            # if re.match(r"###", line):
            if line.startswith('#') & length != 0:
                block_len.append(length)
                length = 0
            else:
                length += 1
        block_len.append(length)

    # read mcscan and quota_align output anchors file as a dataframe.
    collinearity_df = pd.read_table(file, comment="#", header=None, sep='\t', low_memory=False)
    # print(collinearity_df)
    # collinearity_df.columns = ["queryGene", "order_1", "refGene", "order_2", "direction"]
    # get two gene lists and delete "gene:" string
    ref_gene_list = list(collinearity_df.iloc[:, 1])   # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.iloc[:, 0])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    for i in range(len(query_gene_list)):
        query_gene_list[i] = query_gene_list[i].split(sep=":")[1]
    for i in range(len(ref_gene_list)):
        ref_gene_list[i] = ref_gene_list[i].split(sep=":")[1]
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    # print(query_to_ref)
    return query_to_ref, block_len


def read_expression_matrix_xlsx(file, bool_log):
    maize_root = pd.read_excel(file, header=0, sheet_name="YM_root")
    if bool_log == "T":
        numeric_columns = maize_root.select_dtypes(include=np.number).columns
        maize_root[numeric_columns] = maize_root[numeric_columns].applymap(lambda x: np.log(x+1))
    maize_root.columns.values[0] = "gene_name"

    sorghum_root = pd.read_excel(file, header=0, sheet_name="GL_root")
    if bool_log == "T":
        numeric_columns = sorghum_root.select_dtypes(include=np.number).columns
        sorghum_root[numeric_columns] = sorghum_root[numeric_columns].applymap(lambda x: np.log(x+1))
    sorghum_root.columns.values[0] = "gene_name"
    return sorghum_root, maize_root

#####################################################################################################################
#                                         len(query_df_o) = len(ref_df_o)                                           #
#####################################################################################################################


def ref_equal_query_line(query_df, ref_df, query_to_ref):
    query_series_list = list(query_df.iloc[:, 0])
    ref_series_list = list(ref_df.iloc[:, 0])
    headers = [column_name[2:] if column_name.startswith("Y") else column_name for column_name in ref_df.columns]
    query_df.columns = headers
    ref_df.columns = headers
    query_df_o = pd.DataFrame(columns=headers)
    ref_df_o = pd.DataFrame(columns=headers)
    my_series = pd.Series([0] * 16, index=headers, dtype=object)
    for query, ref in query_to_ref:
        if (query in query_series_list) or (ref in ref_series_list):
            if query in query_series_list:
                index1 = query_series_list.index(query)
                query_df_o.loc[len(query_df_o)] = query_df.iloc[index1]
#                query_df_o._append(query_df.iloc[index1], ignore_index=True, inplace=True)
            else:
                my_series.iloc[0] = query
                query_df_o.loc[len(query_df_o)] = my_series
#                query_df_o._append(my_series, ignore_index=True, inplace=True)
            if ref in ref_series_list:
                index2 = ref_series_list.index(ref)
                ref_df_o.loc[len(ref_df_o)] = ref_df.iloc[index2]
#                ref_df_o._append(ref_df.iloc[index2], ignore_index=True, inplace=True)
            else:
                my_series.iloc[0] = ref
                ref_df_o.loc[len(ref_df_o)] = my_series
#                ref_df_o._append(ref, ignore_index=True, inplace=True)
#         query_df_o["sum"] = query_df_o.iloc[:, 1:-1].sum(axis=1)
#         ref_df_o["sum"] = ref_df_o.iloc[:, 1:-1].sum(axis=1)
    del query_df, ref_df
    return query_df_o, ref_df_o


#####################################################################################################################
#              collinearity file has block1 which includes 365 gene pairs.                                          #
#####################################################################################################################
def total_block(query_df, ref_df, query_to_ref, column, ref_df_columns, is_random):
    if is_random.upper() == "T":
        random.shuffle(query_to_ref)
    # by set 1,2,3...,15 to using 15 different column
    headers = [column_name[2:] if not column_name.startswith("g") else column_name for column_name in ref_df_columns]
    # delete YM or GL
    query_df.columns = headers
    ref_df.columns = headers
    query_df = query_df.loc[:, ["gene_name", str(column)]]
    ref_df = ref_df.loc[:, ["gene_name", str(column)]]
    # print(query_df)

    query_df_o = pd.DataFrame(columns=["gene_name", str(column)])
    ref_df_o = pd.DataFrame(columns=["gene_name", str(column)])
    my_series = pd.Series([0] * 2, index=["gene_name", str(column)], dtype=object)
    query_df_series = list(query_df.iloc[:, 0])
    ref_df_series = list(ref_df.iloc[:, 0])
    # print(ref_df_series)
    for query, ref in query_to_ref:
        # print(query)
        # print(query_df_series)
        if query in query_df_series:
            index1 = query_df_series.index(query)
            query_df_o.loc[len(query_df_o)] = query_df.iloc[index1]
            # print(query_df_o)
            # print(query_df_o)
        else:
            my_series.iloc[0] = query
            query_df_o.loc[len(query_df_o)] = my_series

        if ref in ref_df_series:
            index2 = ref_df_series.index(ref)
            ref_df_o.loc[len(ref_df_o)] = ref_df.iloc[index2]
        else:
            my_series.iloc[0] = ref
            ref_df_o.loc[len(ref_df_o)] = my_series
    # query_df_o["sum"] = query_df_o.iloc[:, 1:-1].sum(axis=1)
    # ref_df_o["sum"] = ref_df_o.iloc[:, 1:-1].sum(axis=1)
    del query_df, ref_df
    return query_df_o, ref_df_o
