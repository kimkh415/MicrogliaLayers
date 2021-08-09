#! /bin/bash


set -e

meta_file=$1
count_file=$2
out_path=$3

cellphonedb method statistical_analysis --threads=8 --subsampling --subsampling-log true ${meta_file} ${count_file} --output-path ${out_path}

echo "DONE"
