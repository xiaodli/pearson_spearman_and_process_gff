from scipy.stats import pearsonr, ConstantInputWarning, spearmanr, PermutationMethod
import warnings
# import numpy as np
# from scipy import stats
warnings.filterwarnings('error', category=ConstantInputWarning)


def row_corr_analysis(query_df, ref_df, method):
    dict_corr = {}
    for i in range(len(query_df)):
        # if (query_df.loc[i, "sum"] != 0) and (ref_df.loc[i, "sum"] != 0):
        list1 = query_df.iloc[i, 1:]
        list2 = ref_df.iloc[i, 1:]
        name_1 = query_df.iloc[i, 0]
        name_2 = ref_df.iloc[i, 0]
        if method == "pearson":
            try:
                # correlation_coefficient, p_value = pearsonr(list1, list2, alternative="greater", method=perm_method)
                correlation_coefficient = pearsonr(list1, list2, alternative="greater").statistic
                dict_corr[name_1 + "_" + name_2] = correlation_coefficient
            except ConstantInputWarning as e:
                print(e)
                print("ref or query gene don't express in the early period of the plant ")
        if method == "spearman":
            try:
                correlation_coefficient = spearmanr(list1, list2, alternative="greater").statistic
                # ref = stats.permutation_test((list1, list2), statistic, alternative='greater', permutation_type='pairings')
                dict_corr[name_1 + "_" + name_2] = correlation_coefficient
            except ConstantInputWarning as e:
                print(e)
                print("ref or query gene don't express in the early period of the plant ")
    return dict_corr


def total_block_corr(query_df, ref_df, block_len, method, length, two_block):
    if two_block.upper() == "T":
        dict_corr = {}
        block_index = 0
        start = 0
        end = 0
        for size in block_len:
            end = end + size
            if size == end:
                continue
            if (end - start) >= int(length):
                block_index += 1
                list1 = query_df.iloc[start:end, 1]
                list2 = ref_df.iloc[start:end, 1]
                assert len(list2) == len(list1)
                if method == "pearson":
                    try:
                        correlation_coefficient, _ = pearsonr(list1, list2)
                        dict_corr[block_index] = correlation_coefficient
                    except ConstantInputWarning:
                        print("pearsonr don't agree with constant input,in other words query block or ref block expression sum is zero.")
                else:
                    try:
                        correlation_coefficient, _ = spearmanr(list1, list2)
                        dict_corr[block_index] = correlation_coefficient
                    except ConstantInputWarning:
                        print("pearsonr don't agree with constant input,in other words query block or ref block expression sum is zero.")
            start = start + block_len[block_index - 1]
        return dict_corr
    else:
        dict_corr = {}
        block_index = 1
        start = 0
        for size in block_len:
            end = start + size
            if size >= int(length):
                list1 = query_df.iloc[start:end, 1]
                list2 = ref_df.iloc[start:end, 1]
                assert len(list2) == len(list1)
                if method == "pearson":
                    try:
                        correlation_coefficient, _ = pearsonr(list1, list2)
                        dict_corr[block_index] = correlation_coefficient
                    except ConstantInputWarning:
                        print("pearsonr don't agree with constant input,in other words query block or ref block expression sum is zero.")
                else:
                    try:
                        correlation_coefficient, _ = spearmanr(list1, list2)
                        dict_corr[block_index] = correlation_coefficient
                    except ConstantInputWarning:
                        print("pearsonr don't agree with constant input,in other words query block or ref block expression sum is zero.")
            block_index += 1
            start = end
        return dict_corr
