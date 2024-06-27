import pandas as pd
import random


# step(2)去除重复，保留离block边缘比较远的基因对
# 对每一个group(by=[0,1,2,3])选择距离最大的基因对
# 如果距离相同,选择same
# 然后再次去重df
def drop_duplicate_core_record(block_list):
    record_df = pd.DataFrame(block_list)

    # distance to edge(select max)
    max_values = record_df.groupby(by=[0, 1, 2, 3])[5].transform('max')
    df_filtered = record_df[max_values == record_df[5]]

    # TODO
    # rm all same record may be useful for difference expression-->(def drop_duplicate_all(block_list):)
    # consider remove same distance(rm this part)
    # distance to edge(select same) if same max value
    first_filtered_to_grouped = df_filtered.groupby(by=[0, 1, 2, 3])

    need_to_select_same_direction_groups = first_filtered_to_grouped.filter(lambda x: len(x) > 1)
    short_groups = first_filtered_to_grouped.filter(lambda x: len(x) == 1)
    filtered_long_groups = need_to_select_same_direction_groups.groupby(by=[0, 1, 2, 3]).apply(
        lambda gop: gop[gop[5] == 'same']
    )
    result = pd.concat([filtered_long_groups, short_groups])
    result.reset_index(drop=True, inplace=True)
    result.drop_duplicates(subset=[0, 1, 2, 3], keep=False)
    result.columns = ["query_name", "query_order", "ref_name", "ref_order", "same_or_inverse", "distance_block_edge", "block_index"]
    return result


def drop_duplicate_all(block_list):
    record_df = pd.DataFrame(block_list)
    result = record_df.drop_duplicates(subset=[0, 1, 2, 3], keep=False)
    result.columns = ["query_name", "query_order", "ref_name", "ref_order", "same_or_inverse", "distance_block_edge", "block_index"]
    return result


def split_inv_normal_to_pair(df):
    same_direction_short_groups = df.query('`same_or_inverse` == "same"')
    same_query_sorghum = same_direction_short_groups.loc[:, "query_name"].to_list()
    same_ref_maize = same_direction_short_groups.loc[:, "ref_name"].to_list()
    same = list(zip(same_query_sorghum, same_ref_maize))

    inv_direction = df.query('`same_or_inverse` == "inverse"')
    inv_query_sorghum = inv_direction.loc[:, "query_name"].to_list()
    inv_ref_maize = inv_direction.loc[:, "ref_name"].to_list()
    inv = list(zip(inv_query_sorghum, inv_ref_maize))
    random.shuffle(same)
    random.shuffle(inv)
    return same, inv


def get_expr_dict(query_df, ref_df, column, ref_df_columns):
    headers = [column_name[2:] if not column_name.startswith("g") else column_name for column_name in ref_df_columns]
    # delete YM or GL
    query_df.columns = headers
    ref_df.columns = headers
    new_column_list = ["gene_name", "R1-" + column, "R2-" + column, "R3-" + column]
    mean_list_column = ["R1-" + column, "R2-" + column, "R3-" + column]
    query_df = query_df.loc[:, new_column_list]
    query_df[column] = query_df[mean_list_column].mean(axis=1)
    query_df = query_df.loc[:, ["gene_name", column]]
    query_df.set_index("gene_name", inplace=True)
    query_dict = query_df[column].to_dict()

    ref_df = ref_df.loc[:, new_column_list]
    ref_df[column] = ref_df[mean_list_column].mean(axis=1)
    ref_df = ref_df.loc[:, ["gene_name", column]]
    ref_df.set_index("gene_name", inplace=True)
    ref_dict = ref_df[column].to_dict()
    return query_dict, ref_dict


# same inv can be replaced by recent or pre or homologous
def get_corr_data(query_expr_dict, ref_expr_dict, same_or_inv, bootstrap, size):
    new_same_or_inv_expr_list = []
    for query, ref in same_or_inv:
        if query[5:] not in query_expr_dict or ref[5:] not in ref_expr_dict:

            continue
        else:
            query_expr_value = query_expr_dict[query[5:]]
            ref_expr_value = ref_expr_dict[ref[5:]]
            new_same_or_inv_expr_list.append((query_expr_value, ref_expr_value))
    total_info = [random.sample(new_same_or_inv_expr_list, size) for _ in range(bootstrap)]
    return total_info


def tuple_ele_list_to_two_list(sub_list):
    query_list = []
    ref_list = []
    for query, ref in sub_list:
        query_list.append(query)
        ref_list.append(ref)
    return query_list, ref_list
