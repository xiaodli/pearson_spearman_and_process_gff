#!/bin/bash
input_dir="/home/xiaodong/Desktop/recent/strand_WGD/gene_pair_shuffle_correlation"
constant="/three_result/four_software_"
folder_array=("shoot_u1" "shoot_u2" "shoot_u3" "root_d1" "root_d2" "root_d3")
shoot_or_root_array=("shoot" "shoot" "shoot" "root" "root" "root")
column_array=("u1" "u2" "u3" "d1" "d2" "d3")

for folder_index in ${!folder_array[@]}; do {
        column=${column_array[folder_index]}
        shoot_or_root=${shoot_or_root_array[folder_index]}
        out_path1=${input_dir}${constant}${folder_array[folder_index]}/block_10/raw
        out_path2=${input_dir}${constant}${folder_array[folder_index]}/block_10
        mkdir -p "${out_path1}"
        

        python main.py three_bootstrap_corr -A zm1_sb2.collinearity -AT zm20_sb20.collinearity -BT maize.sorghum.blastp -e all_special_tpm.xlsx -co ${out_path1}/${column}_pearson_correlation_log.csv -m pearson -c ${column} -b T -B 0 -s ${shoot_or_root} -boot 10000 -sample_size 50
        python main.py process -i ${out_path1}/${column}_pearson_correlation_log.csv -o ${out_path2}/${column}_pearson_correlation_log.R.csv

        python main.py three_bootstrap_corr -A zm1_sb2.collinearity -AT zm20_sb20.collinearity -BT maize.sorghum.blastp -e all_special_tpm.xlsx -co ${out_path1}/${column}_spearman_correlation_log.csv -m spearman -c ${column} -b T -B 0 -s ${shoot_or_root} -boot 10000 -sample_size 50
        python main.py process -i ${out_path1}/${column}_spearman_correlation_log.csv -o ${out_path2}/${column}_spearman_correlation_log.R.csv
        
        python main.py three_bootstrap_corr -A zm1_sb2.collinearity -AT zm20_sb20.collinearity -BT maize.sorghum.blastp -e all_special_tpm.xlsx -co ${out_path1}/${column}_pearson_correlation.csv -m pearson -c ${column} -b F -B 0 -s ${shoot_or_root} -boot 10000 -sample_size 50
        python main.py process -i ${out_path1}/${column}_pearson_correlation.csv -o ${out_path2}/${column}_pearson_correlation.R.csv
        
        python main.py three_bootstrap_corr -A zm1_sb2.collinearity -AT zm20_sb20.collinearity -BT maize.sorghum.blastp -e all_special_tpm.xlsx -co ${out_path1}/${column}_spearman_correlation.csv -m spearman -c ${column} -b F -B 0 -s ${shoot_or_root} -boot 10000 -sample_size 50
        python main.py process -i ${out_path1}/${column}_spearman_correlation.csv -o ${out_path2}/${column}_spearman_correlation.R.csv
}&
done
wait

