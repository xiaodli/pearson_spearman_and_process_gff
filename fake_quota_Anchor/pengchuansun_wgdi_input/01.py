import pandas as pd
import sys
data = pd.read_csv(sys.argv[1], header=None, sep="\t", comment="#")
data = data[(data[2] == 'mRNA') | (data[2] =='V_gene_segment')]
data = data.loc[:, [0, 8, 3, 4, 6]]
data[8], data[1] = data[8].str.split(';|:|=', expand=True)[5], data[8].str.split(';|:|=', expand=True)[2]   # 看样子是  08346是列名
data.to_csv(sys.argv[2], sep="\t", header=None, index=False)

