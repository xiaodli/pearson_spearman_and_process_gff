import pandas as pd


def process(parameter):
    input_file = parameter.input
    output_file1 = parameter.output1
    output_file2 = parameter.output2
    output_file3 = parameter.output3
    output_file4 = parameter.output4
    output_file5 = parameter.output5
    output_file6 = parameter.output6

    df = pd.read_csv(input_file, header=0, sep=",")
    df = df.iloc[:, [0, 2, 4, 6, 8, 10, 12, 14]]
    # tassel_anchor ear_anchor tassel_wgdi ear_wgdi ......
    df.columns = ["tassel_anchor", "ear_anchor", "tassel_wgdi", "ear_wgdi", "tassel_mcscanx", "ear_mcscanx", "tassel_jcvi", "ear_jcvi"]
    # tassel_gene_name4 tassel_corr_value4

    df1_tassel = df.loc[:, ["tassel_anchor", "tassel_wgdi"]]
    df1_ear = df.loc[:, ["ear_anchor", "ear_wgdi"]]

    df2_tassel = df.loc[:, ["tassel_anchor", "tassel_mcscanx"]]
    df2_ear = df.loc[:, ["ear_anchor", "ear_mcscanx"]]

    df3_tassel = df.loc[:, ["tassel_anchor", "tassel_jcvi"]]
    df3_ear = df.loc[:, ["ear_anchor", "ear_jcvi"]]

    df1_tassel = df1_tassel.drop_duplicates(keep='first')
    df1_ear = df1_ear.drop_duplicates(keep='first')
    df2_tassel = df2_tassel.drop_duplicates(keep='first')
    df2_ear = df2_ear.drop_duplicates(keep='first')
    df3_tassel = df3_tassel.drop_duplicates(keep='first')
    df3_ear = df3_ear.drop_duplicates(keep='first')

    df1_tassel.to_csv(output_file1, header=True, index=False)
    df1_ear.to_csv(output_file2, header=True, index=False)

    df2_tassel.to_csv(output_file3, header=True, index=False)
    df2_ear.to_csv(output_file4, header=True, index=False)

    df3_tassel.to_csv(output_file5, header=True, index=False)
    df3_ear.to_csv(output_file6, header=True, index=False)
