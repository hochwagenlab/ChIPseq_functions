#!/bin/bash

# To use this script for paired-end sequencing data, 
# fill in the next six lines and
# submit as `source ChIPseq_PE_example.sh`

# Set PE="False" if you want to treat data as single-end reads for pileups
# You will also want to adjust M2FILEPE to say "SE" but DON'T touch TREAT or CONTROL

# PE sequencing produces two FASTQ files. You only need to list one for CHIP and INPUT. The scripts will find the paired FASTQ.

CHIP=
TAGC=""
INPUT=
TAGI=""
FLDR=""
GENOME=""

TREAT="Bowtie/$TAGC-$GENOME-PE_B4.sam"
CONTROL="Bowtie/$TAGI-$GENOME-PE_B4.sam"
FOLDER="$TAGC-$GENOME-B4W4"
M2FILEPE="$TAGC-$GENOME-PE_B4W4"

ji1=$(sbatch --export FLDR=$FLDR,FASTQIN=$CHIP,TAG=$TAGC,GENOME=$GENOME \
    ~/ChIPseq_Pipeline_v3/Bowtie_PE_B4.sbatch | awk '{print $NF}')

ji2=$(sbatch --export FLDR=$FLDR,FASTQIN=$INPUT,TAG=$TAGI,GENOME=$GENOME \
    ~/ChIPseq_Pipeline_v3/Bowtie_PE_B4.sbatch | awk '{print $NF}')

sbatch --dependency=afterok:$ji1:$ji2 --kill-on-invalid-dep=yes \
    --export FLDR=$FLDR,TREAT=$TREAT,CONTROL=$CONTROL,FOLDER=$FOLDER,M2FILE=$M2FILEPE,BDG="T",PEAK="BOTH",PE="True"  \
    ~/ChIPseq_Pipeline_v3/MACS2_FE_W4.sbatch

