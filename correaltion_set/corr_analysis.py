from process_file import tuple_ele_list_to_two_list
from scipy.stats import pearsonr, ConstantInputWarning, spearmanr, PermutationMethod
import warnings
warnings.filterwarnings('error', category=ConstantInputWarning)


def size_corr(same_or_inv_ele_tuple, method):
    dict_corr = {}
    index = 1
    for corr in same_or_inv_ele_tuple:
        query_list, ref_list = tuple_ele_list_to_two_list(corr)
        assert len(query_list) == len(ref_list)
        if method == "pearson":
            try:
                correlation_coefficient, _ = pearsonr(query_list, ref_list)
                dict_corr[index] = correlation_coefficient
            except ConstantInputWarning:
                print("pearsonr don't agree with constant input,in other words query block or ref block expression sum is zero.")
        else:
            try:
                correlation_coefficient, _ = spearmanr(query_list, ref_list)
                dict_corr[index] = correlation_coefficient
            except ConstantInputWarning:
                print("pearsonr don't agree with constant input,in other words query block or ref block expression sum is zero.")
        index += 1
    return dict_corr
