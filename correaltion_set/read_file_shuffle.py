import pandas as pd
import numpy as np


def trans_dtype(record_list):
    new_record_list = []
    query_name = str(record_list[0])
    query_order = int(record_list[1])
    ref_name = str(record_list[2])
    ref_order = int(record_list[3])
    gene_pair_direction = str(record_list[4])
    new_record_list.extend([query_name, query_order, ref_name, ref_order, gene_pair_direction])
    return new_record_list


# step(1)选择block大于一定数量的block
# step(2)计算基因对距离block边缘的距离
def read_wgdi_collinearity(wgdi_file, min_block_length):
    total_info = []
    block_info = []
    total_pre_df_core_record_list = []
    flag = True
    with open(wgdi_file) as f:
        block_index = 1
        record_index = 1
        for line in f:
            if line.startswith('#'):
                # append recent block_info
                if block_info:
                    record_index = 1
                    total_info.append(block_info)
                    block_info = []

                header_record_list = line.split()
                # get block length
                block_length = int(header_record_list[5].split("=")[1])
                if block_length < int(min_block_length):
                    flag = False
                    continue
                else:
                    # [core_record_list1, core_record_list2] --> block_pre_df_core_record_list
                    # [index, direction, length]  --> block_info
                    flag = True
                    # get block direction
                    block_direction = header_record_list[7]
                    if block_direction == "plus":
                        block_direction = "1"
                    if block_direction == "minus":
                        block_direction = "-1"
                    block_info.append(block_index)
                    block_info.append(block_direction)
                    block_info.append(block_length)
                    # TODO
                    block_info.append([])
                    block_index += 1

            else:
                # startswith not "#"
                if flag:
                    # length >= min_block_length
                    core_record_list = line.split()
                    record_index_pair = [record_index, block_length+1-record_index]
                    distance_edge = min(record_index_pair)
                    new_core_record_list = trans_dtype(core_record_list)
                    if new_core_record_list[-1] == block_direction:
                        judge_direction = "same"
                    else:
                        judge_direction = "inverse"
                    new_core_record_list[-1] = judge_direction
                    new_core_record_list.append(distance_edge)
                    new_core_record_list.append(block_index)
                    total_pre_df_core_record_list.append(new_core_record_list)
                    record_index += 1
        if block_info:
            total_info.append(block_info)
    return total_info, total_pre_df_core_record_list


def read_anchor_wave_collinearity(file, length):
    query_to_ref = set()
    flag = True
    with open(file) as f:
        _ = next(f)
        _ = next(f)
        for line in f:
            if line.startswith('#'):
                bk_length = int(line.split()[2].split(sep="N=")[1])
                if bk_length < length:
                    flag = False
                    continue
                else:
                    flag = True
                    continue
            else:
                if flag:
                    record_list = line.split()
                    query_to_ref.add((record_list[5], record_list[0]))

    return query_to_ref


# recent wgd. pre-wgd, homologous pair
def read_three_file(col_file1, col_file2, blast_file, length):
    recent_wgd_set = read_anchor_wave_collinearity(col_file1, length)

    pre_wgd_total_set = read_anchor_wave_collinearity(col_file2, length)
    pre_wgd_set = pre_wgd_total_set.difference(recent_wgd_set)

    homologous = pd.read_csv(blast_file, comment="#", header=0, index_col=None, sep="\t")
    homologous.iloc[:, 0] = homologous.iloc[:, 0].astype(str)
    homologous.iloc[:, 1] = homologous.iloc[:, 1].astype(str)
    homologous = homologous[homologous.iloc[:, 0] != homologous.iloc[:, 1]]
    homologous_total = set(zip(homologous.iloc[:, 0], homologous.iloc[:, 1]))
    homologous_set = homologous_total.difference(pre_wgd_total_set)

    return recent_wgd_set, pre_wgd_set, homologous_set


def read_expression_matrix_xlsx(file, bool_log, sheet_name):
    maize_root = pd.read_excel(file, header=0, sheet_name="YM_"+sheet_name)
    maize_root.iloc[:, 0] = maize_root.iloc[:, 0].astype(str)
    if bool_log == "T":
        numeric_columns = maize_root.select_dtypes(include=np.number).columns
        maize_root[numeric_columns] = maize_root[numeric_columns].applymap(lambda x: np.log2(x+1))
    maize_root.columns.values[0] = "gene_name"

    sorghum_root = pd.read_excel(file, header=0, sheet_name="GL_"+sheet_name)
    sorghum_root.iloc[:, 0] = sorghum_root.iloc[:, 0].astype(str)
    if bool_log == "T":
        numeric_columns = sorghum_root.select_dtypes(include=np.number).columns
        sorghum_root[numeric_columns] = sorghum_root[numeric_columns].applymap(lambda x: np.log2(x+1))
    sorghum_root.columns.values[0] = "gene_name"

    return sorghum_root, maize_root
