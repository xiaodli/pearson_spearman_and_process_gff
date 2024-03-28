import pandas as pd


def process(parameter):
    input_file = parameter.input
    output_file = parameter.output

    df = pd.read_csv(input_file, header=0, sep=",")
    columns_list = ["Software", "Corr_Value"]

    df_list = []
    length = len(df)
    column = df.columns
    width = len(column)

    for i in range(width):
        if i % 2 == 0:
            a = column[i]
            str_number = str(a)[-1]
            if str_number == "1":
                list_2 = ["AnchorWave"] * length
                df["AnchorWave"] = list_2
                df1 = df.loc[:, ["AnchorWave", "Corr_Value1"]]
                df1.columns = columns_list
                df_list.append(df1)
            if str_number == "2":
                list_4 = ["WGDI"] * length
                df["WGDI"] = list_4
                df2 = df.loc[:, ["WGDI", "Corr_Value2"]]
                df2.columns = columns_list
                df_list.append(df2)
            if str_number == "3":
                list_6 = ["MCScanX"] * length
                df["MCScanX"] = list_6
                df3 = df.loc[:, ["MCScanX", "Corr_Value3"]]
                df3.columns = columns_list
                df_list.append(df3)
            if str_number == "4":
                list_8 = ["JCVI"] * length
                df["JCVI"] = list_8
                df4 = df.loc[:, ["JCVI", "Corr_Value4"]]
                df4.columns = columns_list
                df_list.append(df4)

    merged_df1 = pd.concat(df_list)
    merged_df1 = merged_df1.dropna(subset=['Corr_Value'])
    merged_df1.to_csv(output_file, header=True, index=False)
