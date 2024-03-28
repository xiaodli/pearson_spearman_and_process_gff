# import pandas as pd
from scipy.stats import pearsonr, ConstantInputWarning, spearmanr
import warnings
# import numpy as np
# import base
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
        # correlation_coefficient, _ = pearsonr(list1, list2)
        # dict_corr[name_1 + name_2] = correlation_coefficient
    return dict_corr
