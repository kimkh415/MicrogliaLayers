#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=4g
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -binding linear:8
#$ -R y

source /broad/software/scripts/useuse
# anaconda3 configuration
# source /stanley/levin_dr/kwanho/anaconda3/etc/profile.d/conda.sh
export PATH=/stanley/levin_dr/kwanho/anaconda3/bin:$PATH:$HOME/bin
use UGER
use .r-3.6.0
use .samtools-1.8
use .hdf5-1.8.16
use .boost-1.70.0
use .java-jdk-1.8.0_181-x86-64
export R_LIBS=/stanley/levin_dr/kwanho/R-3.6_libs

set -e

source /home/unix/kwanho/kwanho/cpdb-venv/bin/activate

meta_file=$1
count_file=$2
out_path=$3

cellphonedb method statistical_analysis --threads=8 --subsampling --subsampling-log true ${meta_file} ${count_file} --output-path ${out_path}

echo "DONE"
