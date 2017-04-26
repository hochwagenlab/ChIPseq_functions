#!/bin/bash
##SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --mem=10MB
#SBATCH --time=5:00
#SBATCH --mail-type=END
#SBATCH --mail-user=tem298@nyu.edu
#SBATCH --job-name=closing_ChIPseq
#SBATCH --output=/scratch/%u/%x_%j.out

module purge

# This function is designed to merge all of the error/output files into a single file
# this function wants three variables: $IDS, $TMPID, and $TAGC

IFS=":" read -ra ID_LIST <<< $IDS

#OUTNAME="ChIPseq_Pipeline_${ID_LIST[0]}.output"

cd /scratch/$USER

if [[ $TAGC != "" ]]; then
    OUTNAME="${TAGC}_${ID_LIST[0]}.output"
else
    OUTNAME="${TMPID}_${ID_LIST[0]}.output"
fi

mv MACS2_pipeline_$TMPID.txt $OUTNAME

echo "moving these files:" 
for i in `echo "${ID_LIST[@]}"`;
do
    ls *$i*.out 
    cat *$i*.out >> $OUTNAME
    rm *$i*.out
done

if [[ $TAGC != "" ]]; then
    mv ${TAGC}_${ID_LIST[0]}.output ${TAGC}_${ID_LIST[0]}.out
else
    mv ${TMPID}_${ID_LIST[0]}.output ${TMPID}_${ID_LIST[0]}.out
fi




