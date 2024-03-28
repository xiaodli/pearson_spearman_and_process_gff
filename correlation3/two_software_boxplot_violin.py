import pandas as pd


def process(parameter):
    input_file = parameter.input
    output_file1 = parameter.output1
    output_file2 = parameter.output2
    output_file3 = parameter.output3

    df = pd.read_csv(input_file, header=0, sep=",")
    df = df.iloc[:, [1, 3, 5, 7, 9, 11, 13, 15]]
    # tassel_anchor ear_anchor tassel_wgdi ear_wgdi ......
    df.columns = ["t_a", "e_a", "t_w", "e_w", "t_m", "e_m", "t_j", "e_j"]
    # tassel_gene_name4 tassel_corr_value4
    length = len(df)

    list_1 = ["tassel"] * length
    list_2 = ["ear"] * length
    list_3 = ["AnchorWave"] * length
    list_4 = ["WGDI"] * length
    list_5 = ["MCScanX"] * length
    list_6 = ["JCVI"] * length

    df["tassel"] = list_1
    df["ear"] = list_2
    df["AnchorWave"] = list_3
    df["WGDI"] = list_4
    df["MCScanX"] = list_5
    df["JCVI"] = list_6
    # AnchorWave WGDI
    columns_list = ["corr_value", "tas_ear", "software"]
    df1 = df.loc[:, ["t_a", "tassel", "AnchorWave"]]
    df1.drop_duplicates(keep='first')
    df1.columns = columns_list
    df2 = df.loc[:, ["t_w", "tassel", "WGDI"]]
    df2.columns = columns_list
    df2.drop_duplicates(keep='first')
    merged_df1 = pd.concat([df1, df2])
    merged_df1 = merged_df1.dropna(subset=['corr_value'])

    df3 = df.loc[:, ["e_a", "ear", "AnchorWave"]]
    df3.columns = columns_list
    df3.drop_duplicates(keep='first')
    df4 = df.loc[:, ["e_w", "ear", "WGDI"]]
    df4.columns = columns_list
    df4.drop_duplicates(keep='first')
    merged_df2 = pd.concat([df3, df4])
    merged_df2 = merged_df2.dropna(subset=['corr_value'])

    merged_df_o1 = pd.concat([merged_df1, merged_df2])
    merged_df_o1.fillna('')
    merged_df_o1.to_csv(output_file1, header=True, index=False)

    # AnchorWave MCScanX
    df1 = df.loc[:, ["t_a", "tassel", "AnchorWave"]]
    df1.columns = columns_list
    df1.drop_duplicates(keep='first')
    df2 = df.loc[:, ["t_m", "tassel", "MCScanX"]]
    df2.columns = columns_list
    df2.drop_duplicates(keep='first')
    merged_df1 = pd.concat([df1, df2])
    merged_df1 = merged_df1.dropna(subset=['corr_value'])

    df3 = df.loc[:, ["e_a", "ear", "AnchorWave"]]
    df3.columns = columns_list
    df3.drop_duplicates(keep='first')
    df4 = df.loc[:, ["e_m", "ear", "MCScanX"]]
    df4.columns = columns_list
    df4.drop_duplicates(keep='first')
    merged_df2 = pd.concat([df3, df4])
    merged_df2 = merged_df2.dropna(subset=['corr_value'])

    merged_df_o2 = pd.concat([merged_df1, merged_df2])
    merged_df_o2.fillna('')
    merged_df_o2.to_csv(output_file2, header=True, index=False)

    # AnchorWave jcvi
    df1 = df.loc[:, ["t_a", "tassel", "AnchorWave"]]
    df1.columns = columns_list
    df1.drop_duplicates(keep='first')
    df2 = df.loc[:, ["t_j", "tassel", "JCVI"]]
    df2.columns = columns_list
    df2.drop_duplicates(keep='first')
    merged_df1 = pd.concat([df1, df2])
    merged_df1 = merged_df1.dropna(subset=['corr_value'])

    df3 = df.loc[:, ["e_a", "ear", "AnchorWave"]]
    df3.columns = columns_list
    df3.drop_duplicates(keep='first')
    df4 = df.loc[:, ["e_j", "ear", "JCVI"]]
    df4.columns = columns_list
    df4.drop_duplicates(keep='first')
    merged_df2 = pd.concat([df3, df4])
    merged_df2 = merged_df2.dropna(subset=['corr_value'])

    merged_df_o3 = pd.concat([merged_df1, merged_df2])
    merged_df_o3.fillna('')
    merged_df_o3.to_csv(output_file3, header=True, index=False)
