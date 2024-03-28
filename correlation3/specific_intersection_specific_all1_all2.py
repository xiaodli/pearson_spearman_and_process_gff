import pandas as pd


# intersection wgdi_specific anchor_specifc anchor wgdi
def specific_all_inner(parameter):
    file = parameter.input
    output1 = parameter.output1
    output2 = parameter.output2
    output3 = parameter.output3
    df = pd.read_csv(file, header=0)
    tassel_anchor = df.iloc[:, 0:2]
    tassel_anchor = tassel_anchor.dropna(subset=['tassel_corr_value1'])
    ear_anchor = df.iloc[:, 2:4]
    ear_anchor = ear_anchor.dropna(subset=['ear_corr_value1'])

    tassel_wgdi = df.iloc[:, 4:6]
    tassel_wgdi = tassel_wgdi.dropna(subset=['tassel_corr_value2'])
    ear_wgdi = df.iloc[:, 6:8]
    ear_wgdi = ear_wgdi.dropna(subset=['ear_corr_value2'])

    tassel_mc = df.iloc[:, 8:10]
    tassel_mc = tassel_mc.dropna(subset=['tassel_corr_value3'])
    ear_mc = df.iloc[:, 10:12]
    ear_mc = ear_mc.dropna(subset=['ear_corr_value3'])

    tassel_jcvi = df.iloc[:, 12:14]
    tassel_jcvi = tassel_jcvi.dropna(subset=['tassel_corr_value4'])
    ear_jcvi = df.iloc[:, 14:16]
    ear_jcvi = ear_jcvi.dropna(subset=['ear_corr_value4'])

    # all gene pair
    tassel_anchor = tassel_anchor.drop_duplicates(keep='first')
    tassel_wgdi = tassel_wgdi.drop_duplicates(keep='first')
    tassel_mc = tassel_mc.drop_duplicates(keep='first')
    tassel_jcvi = tassel_jcvi.drop_duplicates(keep='first')

    ear_anchor = ear_anchor.drop_duplicates(keep='first')
    ear_wgdi = ear_wgdi.drop_duplicates(keep='first')
    ear_mc = ear_mc.drop_duplicates(keep='first')
    ear_jcvi = ear_jcvi.drop_duplicates(keep='first')

    for df in [tassel_anchor, tassel_wgdi, tassel_mc, tassel_jcvi]:
        # df = reduce_mem_usage(df)
        df.columns = ["tassel_gene_name", "tassel_corr_value"]
    for df in [ear_anchor, ear_wgdi, ear_mc, ear_jcvi]:
        # df = reduce_mem_usage(df)
        df.columns = ["ear_gene_name", "ear_corr_value"]
    # intersection
    tassel_inner_anchor_wgdi = pd.merge(tassel_anchor, tassel_wgdi, how="inner",
                                        on=["tassel_gene_name", "tassel_corr_value"]).drop_duplicates(keep='first')
    tassel_inner_anchor_mc = (pd.merge(tassel_anchor, tassel_mc, how="inner",
                                       on=["tassel_gene_name", "tassel_corr_value"]).drop_duplicates(keep='first'))
    tassel_inner_anchor_jcvi = pd.merge(tassel_anchor, tassel_jcvi, how="inner",
                                        on=["tassel_gene_name", "tassel_corr_value"]).drop_duplicates(keep='first')
    ear_inner_anchor_wgdi = pd.merge(ear_anchor, ear_wgdi, how="inner",
                                     on=["ear_gene_name", "ear_corr_value"]).drop_duplicates(keep='first')
    ear_inner_anchor_mc = pd.merge(ear_anchor, ear_mc, how="inner",
                                   on=["ear_gene_name", "ear_corr_value"]).drop_duplicates(keep='first')
    ear_inner_anchor_jcvi = pd.merge(ear_anchor, ear_jcvi, how="inner",
                                     on=["ear_gene_name", "ear_corr_value"]).drop_duplicates(keep='first')

    # specific
    tassel_only_anchor1 = pd.concat([tassel_anchor, tassel_inner_anchor_wgdi]).drop_duplicates(keep='first')
    tassel_only_wgdi = pd.concat([tassel_wgdi, tassel_inner_anchor_wgdi]).drop_duplicates(keep='first')
    ear_only_anchor1 = pd.concat([ear_anchor, ear_inner_anchor_wgdi]).drop_duplicates(keep='first')
    ear_only_wgdi = pd.concat([ear_wgdi, ear_inner_anchor_wgdi]).drop_duplicates(keep='first')

    tassel_only_anchor2 = pd.concat([tassel_anchor, tassel_inner_anchor_mc]).drop_duplicates(keep='first')
    tassel_only_mc = pd.concat([tassel_wgdi, tassel_inner_anchor_mc]).drop_duplicates(keep='first')
    ear_only_anchor2 = pd.concat([ear_anchor, ear_inner_anchor_mc]).drop_duplicates(keep='first')
    ear_only_mc = pd.concat([ear_mc, ear_inner_anchor_mc]).drop_duplicates(keep='first')

    tassel_only_anchor3 = pd.concat([tassel_anchor, tassel_inner_anchor_jcvi]).drop_duplicates(keep='first')
    tassel_only_jcvi = pd.concat([tassel_jcvi, tassel_inner_anchor_jcvi]).drop_duplicates(keep='first')
    ear_only_anchor3 = pd.concat([ear_anchor, ear_inner_anchor_jcvi]).drop_duplicates(keep='first')
    ear_only_jcvi = pd.concat([ear_jcvi, ear_inner_anchor_wgdi]).drop_duplicates(keep='first')

    list1 = [tassel_anchor, tassel_wgdi, tassel_inner_anchor_wgdi, tassel_only_anchor1, tassel_only_wgdi,
             ear_anchor, ear_wgdi, ear_inner_anchor_wgdi, ear_only_anchor1, ear_only_wgdi]
    list2 = [tassel_mc, tassel_inner_anchor_mc, tassel_only_anchor2, tassel_only_mc,
             ear_mc, ear_inner_anchor_mc, ear_only_anchor2, ear_only_mc]
    list3 = [tassel_jcvi, tassel_inner_anchor_jcvi, tassel_only_anchor3, tassel_only_jcvi,
             ear_jcvi, ear_inner_anchor_jcvi, ear_only_anchor3, ear_only_jcvi]
    for lt in [list1, list2, list3]:
        for df in lt:
            # print(df.columns)
            df.drop(df.columns[0], axis=1, inplace=True)
            df.columns = ["Corr_value"]
    tassel_wgdi["type"] = ["WGDI"] * len(tassel_wgdi)
    tassel_wgdi["tas_ear"] = ["Tassel"] * len(tassel_wgdi)
    tassel_anchor["type"] = ["AnchorWave"] * len(tassel_anchor)
    tassel_anchor["tas_ear"] = ["Tassel"] * len(tassel_anchor)
    tassel_mc["type"] = ["MCScanX"] * len(tassel_mc)
    tassel_mc["tas_ear"] = ["Tassel"] * len(tassel_mc)
    tassel_jcvi["type"] = ["JCVI"] * len(tassel_jcvi)
    tassel_jcvi["tas_ear"] = ["Tassel"] * len(tassel_jcvi)
    ear_anchor["type"] = ["AnchorWave"] * len(ear_anchor)
    ear_anchor["tas_ear"] = ["Ear"] * len(ear_anchor)
    ear_wgdi["type"] = ["WGDI"] * len(ear_wgdi)
    ear_wgdi["tas_ear"] = ["Ear"] * len(ear_wgdi)
    ear_mc["type"] = ["MCScanX"] * len(ear_mc)
    ear_mc["tas_ear"] = ["Ear"] * len(ear_mc)
    ear_jcvi["type"] = ["JCVI"] * len(ear_jcvi)
    ear_jcvi["tas_ear"] = ["Ear"] * len(ear_jcvi)
    tassel_inner_anchor_wgdi['type'] = ["Intersection"] * len(tassel_inner_anchor_wgdi)
    tassel_inner_anchor_wgdi['tas_ear'] = ["Tassel"] * len(tassel_inner_anchor_wgdi)
    tassel_inner_anchor_mc['type'] = ["Intersection"] * len(tassel_inner_anchor_mc)
    tassel_inner_anchor_mc['tas_ear'] = ["Tassel"] * len(tassel_inner_anchor_mc)
    tassel_inner_anchor_jcvi['type'] = ["Intersection"] * len(tassel_inner_anchor_jcvi)
    tassel_inner_anchor_jcvi['tas_ear'] = ["Tassel"] * len(tassel_inner_anchor_jcvi)
    ear_inner_anchor_wgdi['type'] = ["Intersection"] * len(ear_inner_anchor_wgdi)
    ear_inner_anchor_wgdi['tas_ear'] = ["Ear"] * len(ear_inner_anchor_wgdi)
    ear_inner_anchor_mc['type'] = ["Intersection"] * len(ear_inner_anchor_mc)
    ear_inner_anchor_mc['tas_ear'] = ["Ear"] * len(ear_inner_anchor_mc)
    ear_inner_anchor_jcvi['type'] = ["Intersection"] * len(ear_inner_anchor_jcvi)
    ear_inner_anchor_jcvi['tas_ear'] = ["Ear"] * len(ear_inner_anchor_jcvi)
    tassel_only_anchor1["type"] = ["AnchorWave_specific"] * len(tassel_only_anchor1)
    tassel_only_anchor1["tas_ear"] = ["Tassel"] * len(tassel_only_anchor1)
    tassel_only_wgdi["type"] = ["WGDI_specific"] * len(tassel_only_wgdi)
    tassel_only_wgdi["tas_ear"] = ["Tassel"] * len(tassel_only_wgdi)
    ear_only_anchor1["type"] = ["AnchorWave_specific"] * len(ear_only_anchor1)
    ear_only_anchor1["tas_ear"] = ["Ear"] * len(ear_only_anchor1)
    ear_only_wgdi["type"] = ["WGDI_specific"] * len(ear_only_wgdi)
    ear_only_wgdi["tas_ear"] = ["Ear"] * len(ear_only_wgdi)
    tassel_only_anchor2["type"] = ["AnchorWave_specific"] * len(tassel_only_anchor2)
    tassel_only_anchor2["tas_ear"] = ["Tassel"] * len(tassel_only_anchor2)
    tassel_only_mc["type"] = ["MCScanX_specific"] * len(tassel_only_mc)
    tassel_only_mc["tas_ear"] = ["Tassel"] * len(tassel_only_mc)
    ear_only_anchor2["type"] = ["AnchorWave_specific"] * len(ear_only_anchor2)
    ear_only_anchor2["tas_ear"] = ["Ear"] * len(ear_only_anchor2)
    ear_only_mc["type"] = ["MCScanX_specific"] * len(ear_only_mc)
    ear_only_mc["tas_ear"] = ["Ear"] * len(ear_only_mc)
    tassel_only_anchor3["type"] = ["AnchorWave_specific"] * len(tassel_only_anchor3)
    tassel_only_anchor3["tas_ear"] = ["Tassel"] * len(tassel_only_anchor3)
    tassel_only_jcvi["type"] = ["JCVI_specific"] * len(tassel_only_jcvi)
    tassel_only_jcvi["tas_ear"] = ["Tassel"] * len(tassel_only_jcvi)
    ear_only_anchor3["type"] = ["AnchorWave_specific"] * len(ear_only_anchor3)
    ear_only_anchor3["tas_ear"] = ["Ear"] * len(ear_only_anchor3)
    ear_only_jcvi["type"] = ["JCVI_specific"] * len(ear_only_jcvi)
    ear_only_jcvi["tas_ear"] = ["Ear"] * len(ear_only_jcvi)

    merge_df1 = pd.concat([tassel_anchor, tassel_wgdi, tassel_inner_anchor_wgdi, tassel_only_anchor1, tassel_only_wgdi,
                           ear_anchor, ear_wgdi, ear_inner_anchor_wgdi, ear_only_anchor1, ear_only_wgdi])
    merge_df1.to_csv(output1, header=True, index=False)
    merge_df2 = pd.concat([tassel_anchor, tassel_mc, tassel_inner_anchor_mc, tassel_only_anchor2, tassel_only_mc,
                           ear_anchor, ear_mc, ear_inner_anchor_mc, ear_only_anchor2, ear_only_mc])
    merge_df2.to_csv(output2, header=True, index=False)
    merge_df3 = pd.concat([tassel_anchor, tassel_jcvi, tassel_inner_anchor_jcvi, tassel_only_anchor3, tassel_only_jcvi,
                           ear_anchor, ear_jcvi, ear_inner_anchor_jcvi, ear_only_anchor3, ear_only_jcvi])
    merge_df3.to_csv(output3, header=True, index=False)
