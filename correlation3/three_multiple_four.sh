#!/bin/bash
input_dir=/home/xiaodong/Desktop/recent/correlation3/result/gene_pair/raw        
output_dir=/home/xiaodong/Desktop/recent/correlation3/result/gene_pair   
for input_file in $input_dir/*.csv; do {
    filename=$(basename $input_file)
    out_path=$output_dir/five_${filename%.csv}
    python main.py process_five -i $input_file -o1 $out_path/one.csv -o2 $out_path/two.csv -o3 $out_path/three.csv 
}&
done
# boxplot violin 3 * 4
