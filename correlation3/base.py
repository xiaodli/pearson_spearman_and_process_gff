# tassel and ear of maize gene expression matrix and sorghum gene expression matrix to collinearity expression matrix
# from argparse import ArgumentParser
import pandas as pd
from itertools import islice
import numpy as np


def read_collinearity(collinearity_result):
    block_len = []
    with open(collinearity_result) as f:
        _ = next(f)
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[2].split(sep="N=")[1]))
    # read anchorWave output collinearity file as a dataframe, every line of the result occurs once different from wgdi.
    collinearity_df = pd.read_table(collinearity_result, comment="#", header=0)

    # get two gene lists and delete "gene-" string
    ref_gene_list = list(collinearity_df.loc[:, "refGene"])   # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.loc[:, "queryGene"])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    for i in range(len(query_gene_list)):
        query_gene_list[i] = query_gene_list[i].split(sep="-")[1]
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    return query_to_ref, block_len


def read_wgdi_collinearity(collinearity_result):
    block_len = []
    with open(collinearity_result) as f:
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[5].split(sep="N=")[1]))
    # read wgdi output collinearity file as a dataframe.
    collinearity_df = pd.read_table(collinearity_result, comment="#", header=None, sep=r'\s+')
    # print(collinearity_df.columns)
    collinearity_df.columns = ["queryGene", "order_1", "refGene", "order_2", "direction"]
    # print(collinearity_df)
    # get two gene lists and delete "gene-" string
    ref_gene_list = list(collinearity_df.loc[:, "refGene"])   # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.loc[:, "queryGene"])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    for i in range(len(query_gene_list)):
        query_gene_list[i] = query_gene_list[i].split(sep="-")[1]
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
    for i in range(len(query_gene_list)):
        query_gene_list[i] = query_gene_list[i].split(sep="-")[1]
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
        query_gene_list[i] = query_gene_list[i].split(sep="-")[1]
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    # print(query_to_ref)
    return query_to_ref, block_len


def read_expression_matrix(tassel_or_ear_expression_matrix, sorghum_expression_matrix, bool_log):
    tassel_or_ear_df = pd.read_csv(tassel_or_ear_expression_matrix, header=0)
    sorghum_df = pd.read_csv(sorghum_expression_matrix, header=0)
    tassel_or_ear_df.columns.values[0] = "gene_name"
    sorghum_df.columns.values[0] = "gene_name"
    if bool_log == "T":
        numeric_columns = tassel_or_ear_df.select_dtypes(include=np.number).columns
        tassel_or_ear_df[numeric_columns] = tassel_or_ear_df[numeric_columns].applymap(lambda x: np.log(x + 1))

    if bool_log == "T":
        numeric_columns = sorghum_df.select_dtypes(include=np.number).columns
        sorghum_df[numeric_columns] = sorghum_df[numeric_columns].applymap(lambda x: np.log(x + 1))
    # print(sorghum_df)
    return tassel_or_ear_df, sorghum_df


# collinearity file has block1 which includes 365 gene pairs, removing gene pair which don't express. In other words, one gene express at least.
def ref_equal_query_line(query_df, ref_df, query_to_ref):
    query_series_list = list(query_df.iloc[:, 0])
    ref_series_list = list(ref_df.iloc[:, 0])
    headers = ref_df.columns.tolist()
    query_df_o = pd.DataFrame(columns=headers)
    ref_df_o = pd.DataFrame(columns=headers)
    my_series = pd.Series([0] * 20, index=headers, dtype=object)
    # print(my_series)
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
