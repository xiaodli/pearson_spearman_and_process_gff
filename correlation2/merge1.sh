#!/bin/bash
input_dir=/home/xiaodong/Desktop/recent/correlation2/result/four_software/block_30
input_dir1=/home/xiaodong/Desktop/recent/correlation2/result/four_software/block_30_two
output_dir=/home/xiaodong/Desktop/recent/correlation2/result/four_software/merge_block30_two
for input_file in $input_dir/*.csv; do {
    filename=$(basename $input_file)
    out_path=$output_dir/merge_${filename}
    python merge.py ${input_dir}/${filename} ${input_dir1}/${filename} ${out_path}
}&
done
