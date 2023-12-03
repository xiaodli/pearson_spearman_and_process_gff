# Combine pengchuan_sun 's 01.py and this file using gff file to get wgdi input (6 columns gff)
# I have had blasp result in the anchorwave analysis process.

import pandas as pd
import sys

data = pd.read_csv(sys.argv[1], header=None, sep="\t")
order = ""
data[2].astype(int)
data[3].astype(int)
data["dif"] = data[3]-data[2]
for name, group in data.groupby([1]):
    nu = len(group)
    if nu == 1:
        continue
    ind = group.sort_values(by="dif", ascending=False).index[1:].values
    data.drop(index=ind, inplace=True)
for ch, group in data.groupby(data[0]):
    num = len(group)
    group.sort_values(by=[2])
    data.loc[group.index, "order"] = list(range(1, num+1))
    data.loc[group.index, "newname"] = list(["zm" + str(chr) + "g" + str(i).zfill(5) for i in range(1, num+1)])

data['order'] = data['order'].astype('int')
data = data[[0, 1, 2, 3, 4, 'order']]
data.to_csv(sys.argv[2], sep="\t", index=False, header=None)
lens = data.groupby([0]).max()[[3, 'order']]
lens.to_csv(sys.argv[3], sep="\t", header=None)
