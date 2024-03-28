import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os


def all_value_df_dotplot1(file):
    file_name_without_extension = os.path.splitext(os.path.basename(file))[0]
    path = os.path.dirname(file)
    df = pd.read_csv(file, header=0, index_col=0)
    numeric_columns = df.select_dtypes(include=np.number).columns
    df[numeric_columns] = df[numeric_columns].applymap(lambda x: np.log(x+1) if x != 0 else 0)
    for column in df.columns:
        my_list = df[column]
        my_list = [x for x in my_list if x != 0]
        plt.figure(figsize=(8, 6))
        sns.histplot(my_list, kde=True, color='blue', stat='density')
        plt.title(f'Histogram for {column}')
        plt.xlabel(column)
        plt.ylabel('Density')

        plt.xlim(0, 4)
        # print(my_list)
        if file_name_without_extension.startswith("s"):
            plt.savefig(f"{path}/result/normal_density/log+1/sorghum/{file_name_without_extension}_{column}.png")
        if file_name_without_extension.startswith("e"):
            plt.savefig(f"{path}/result/normal_density/log+1/ear/{file_name_without_extension}_{column}.png")
        if file_name_without_extension.startswith("t"):
            plt.savefig(f"{path}/result/normal_density/log+1/tassel/{file_name_without_extension}_{column}.png")
        # plt.show()


all_value_df_dotplot1(sys.argv[1])
all_value_df_dotplot1(sys.argv[2])
all_value_df_dotplot1(sys.argv[3])


def all_value_df_dotplot2(file):
    file_name_without_extension = os.path.splitext(os.path.basename(file))[0]
    path = os.path.dirname(file)
    df = pd.read_csv(file, header=0, index_col=0)
    numeric_columns = df.select_dtypes(include=np.number).columns
    df[numeric_columns] = df[numeric_columns].applymap(lambda x: np.log(x) if x != 0 else 0)
    for column in df.columns:
        my_list = df[column]
        my_list = [x for x in my_list if x != 0]
        plt.figure(figsize=(8, 6))
        sns.histplot(my_list, kde=True, color='blue', stat='density')
        plt.title(f'Histogram for {column}')
        plt.xlabel(column)
        plt.ylabel('Density')

        plt.xlim(-6, 4)
        # print(my_list)
        if file_name_without_extension.startswith("s"):
            plt.savefig(f"{path}/result/normal_density/log/sorghum/{file_name_without_extension}_{column}.png")
        if file_name_without_extension.startswith("e"):
            plt.savefig(f"{path}/result/normal_density/log/ear/{file_name_without_extension}_{column}.png")
        if file_name_without_extension.startswith("t"):
            plt.savefig(f"{path}/result/normal_density/log/tassel/{file_name_without_extension}_{column}.png")
        # plt.show()


all_value_df_dotplot2(sys.argv[1])
all_value_df_dotplot2(sys.argv[2])
all_value_df_dotplot2(sys.argv[3])


def all_value_df_dotplot3(file):
    file_name_without_extension = os.path.splitext(os.path.basename(file))[0]
    path = os.path.dirname(file)
    df = pd.read_csv(file, header=0, index_col=0)
    # numeric_columns = df.select_dtypes(include=np.number).columns
    # df[numeric_columns] = df[numeric_columns].applymap(lambda x: np.log(x) if x != 0 else 0)
    for column in df.columns:
        my_list = df[column]
        my_list = [x for x in my_list if x != 0]
        plt.figure(figsize=(8, 6))
        sns.histplot(my_list, kde=True, color='blue', stat='density')
        plt.title(f'Histogram for {column}')
        plt.xlabel(column)
        plt.ylabel('Density')

        plt.xlim(-6, 4)
        # print(my_list)
        if file_name_without_extension.startswith("s"):
            plt.savefig(f"{path}/result/normal_density/tpm/sorghum/{file_name_without_extension}_{column}.png")
        if file_name_without_extension.startswith("e"):
            plt.savefig(f"{path}/result/normal_density/tpm/ear/{file_name_without_extension}_{column}.png")
        if file_name_without_extension.startswith("t"):
            plt.savefig(f"{path}/result/normal_density/tpm/tassel/{file_name_without_extension}_{column}.png")


all_value_df_dotplot3(sys.argv[1])
all_value_df_dotplot3(sys.argv[2])
all_value_df_dotplot3(sys.argv[3])
