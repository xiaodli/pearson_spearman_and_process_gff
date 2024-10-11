#!/bin/bash
input_dir=/home/xiaodong/Desktop/recent/correlation3/result/gene_pair/raw        
output_dir=/home/xiaodong/Desktop/recent/correlation3/result/gene_pair   
for input_file in $input_dir/*.csv; do {
    filename=$(basename $input_file)
    out_path=$output_dir/venn_${filename%.csv}
    python main.py process_venn -i $input_file -o1 $out_path/o1.csv -o2 $out_path/o2.csv -o3 $out_path/o3.csv -o4 $out_path/o4.csv -o5 $out_path/o5.csv -o6 $out_path/o6.csv 
}&
done
# venn plot 4 * 6
