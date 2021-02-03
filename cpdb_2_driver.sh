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
export PATH=/stanley/levin_dr/kwanho/anaconda3/bin:$PATH:$HOME/bin
use UGER
use .hdf5-1.8.16
use .boost-1.70.0

set -e

source /home/unix/kwanho/kwanho/cpdb-venv/bin/activate

cd /stanley/levin_dr/kwanho/projects/microglia/fezf2_wt_ko/ligand_receptor/cpdb_20200530

cellphonedb method statistical_analysis --threads=8 --subsampling --subsampling-log true combined_meta.tsv combined_count.tsv &> stat_analysis.log

echo "DONE"
