import pandas as pd
import sys


input_file = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]

df = pd.read_csv(input_file, header=0, sep=",")
df["random"] = ["random"] * len(df)
columns_list = ["Software", "Corr_Value", "status"]
df .columns = columns_list

df1 = pd.read_csv(input_file2, header=0, sep=",")
df1["predictable"] = ["predictable"] * len(df1)
columns_list = ["Software", "Corr_Value", "status"]
df1 .columns = columns_list

merged_df1 = pd.concat([df, df1])
merged_df1 = merged_df1.dropna(subset=['Corr_Value'])
merged_df1.to_csv(output_file, header=True, index=False)
