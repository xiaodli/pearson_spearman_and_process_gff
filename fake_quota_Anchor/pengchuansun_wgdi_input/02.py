import pandas as pd 
import sys

data = pd.read_csv(sys.argv[1], header=None, sep="\t")
order = ""
newname = ""
#print(type(data[0][0]))
#print(data[0][0])
data[2].astype(int)
data[3].astype(int)
data["dif"] = data[3]-data[2]
for name,group in data.groupby([5]):
    nu = len(group)
    if nu == 1:
        continue
    ind = group.sort_values(by="dif").index[1:].values
    data.drop(index=ind, inplace=True)
for chr, group in data.groupby(data[0]):
    num = len(group)
    group.sort_values(by=[2])
    data.loc[group.index, "order"] = list(range(1,num+1))
    data.loc[group.index, "newname"] = list(["zm" + str(chr) + "g" + str(i).zfill(5) for i in range(1,num+1)])
    #print(type(chr))
    #ind.append(chr)
data['order'] = data['order'].astype('int')
data = data[[0,'newname',2,3,4,'order',1]]
data.to_csv(sys.argv[2],sep="\t",index=False, header=None)
lens = data.groupby([0]).max()[[3, 'order']]
#lens.index = ind
lens.to_csv(sys.argv[3], sep="\t",header=None)
#lens = data.groupby([0]).max()[[3, 'order']]
#lens.to_csv(sys.argv[3], sep="\t",header=None)
