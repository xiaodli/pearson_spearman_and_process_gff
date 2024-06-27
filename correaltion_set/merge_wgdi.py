import pandas as pd
import sys


input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
input_file3 = sys.argv[3]
input_file4 = sys.argv[4]
input_file5 = sys.argv[5]
input_file6 = sys.argv[6]

output_file = sys.argv[7]
columns_list = ["Software", "Corr_Value", "status"]

df1 = pd.read_csv(input_file1, header=0, sep=",")
df1["status"] = ["shoot_0_h"] * len(df1)
df1.columns = columns_list

df2 = pd.read_csv(input_file2, header=0, sep=",")
df2["status"] = ["shoot_2_h"] * len(df2)
df2.columns = columns_list

df3 = pd.read_csv(input_file3, header=0, sep=",")
df3["status"] = ["shoot_72_h"] * len(df3)
df3.columns = columns_list

df4 = pd.read_csv(input_file4, header=0, sep=",")
df4["status"] = ["root_0_h"] * len(df4)
df4.columns = columns_list

df5 = pd.read_csv(input_file5, header=0, sep=",")
df5["status"] = ["root_2_h"] * len(df5)
df5.columns = columns_list

df6 = pd.read_csv(input_file6, header=0, sep=",")
df6["status"] = ["root_72_h"] * len(df6)
df6.columns = columns_list


merged_df1 = pd.concat([df1, df2, df3, df4, df5, df6])
merged_df1 = merged_df1.dropna(subset=['Corr_Value'])
merged_df1.to_csv(output_file, header=True, index=False)