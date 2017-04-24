#!/bin/bash
##SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --mem=10MB
#SBATCH --time=5:00
#SBATCH --mail-type=END
#SBATCH --mail-user=tem298@hpc.nyu.edu
#SBATCH --job-name=ChIPseq_closing

module purge

# This function is designed to merge all of the error/output files into a single file
# this function wants two variables: $IDS and $OUT

IFS=":" read -ra ID_LIST <<< $IDS

OUTNAME="ChIPseq_Pipeline_${ID_LIST[0]}.output"

cd /scratch/$USER

mv $OUT $OUTNAME

for i in `echo "${ID_LIST[@]}"`;
do
    cat *$i*.out >> $OUTNAME
    rm *$i*.out
done

mv ChIPseq_Pipeline_${ID_LIST[0]}.output ChIPseq_Pipeline_${ID_LIST[0]}.out




