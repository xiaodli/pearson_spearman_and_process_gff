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
                list_2 = ["Same_Direction_Set"] * length
                df["Same_Direction_Set"] = list_2
                df1 = df.loc[:, ["Same_Direction_Set", "Corr_Value1"]]
                df1.columns = columns_list
                df_list.append(df1)
            if str_number == "2":
                list_4 = ["Inverse_Direction_Set"] * length
                df["Inverse_Direction_Set"] = list_4
                df2 = df.loc[:, ["Inverse_Direction_Set", "Corr_Value2"]]
                df2.columns = columns_list
                df_list.append(df2)
            if str_number == "3":
                list_6 = ["Recent_WGD_Pair"] * length
                df["Recent_WGD_pair"] = list_6
                df3 = df.loc[:, ["Recent_WGD_pair", "Corr_Value3"]]
                df3.columns = columns_list
                df_list.append(df3)
            if str_number == "4":
                list_8 = ["Pre_WGD_Pair"] * length
                df["Pre_WGD_Pair"] = list_8
                df4 = df.loc[:, ["Pre_WGD_Pair", "Corr_Value4"]]
                df4.columns = columns_list
                df_list.append(df4)
            if str_number == "5":
                list_10 = ["Homologous_Pair"] * length
                df["Homologous_Pair"] = list_10
                df5 = df.loc[:, ["Homologous_Pair", "Corr_Value5"]]
                df5.columns = columns_list
                df_list.append(df5)

    merged_df1 = pd.concat(df_list)
    merged_df1 = merged_df1.dropna(subset=['Corr_Value'])
    merged_df1.to_csv(output_file, header=True, index=False)
