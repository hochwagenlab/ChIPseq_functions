#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --time=1:00:00
#SBATCH --job-name=MACS2_FE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tem298@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out

### 6GB sam files need 7GB RAM
### 17GB sam files need 7GB RAM

module purge
module load macs2/intel/2.1.1
module load samtools/intel/1.3.1

echo "#####################################################"
echo ""
echo "Running: $SLURM_JOB_NAME"
echo "JobID is: $SLURM_JOBID"
echo ""

####################################
## PURPOSE ##
# To use MACS2 to call peaks and/or make normalized wiggle
# files.
# MACS2 settings: no model, 200bp extension size, SPMR,
#    keep duplicates based upon binomial distribution
#    (or all reads for rDNA),
#    FE normalization of wiggle files
# For Paired-end: MACS2 only accepts BAM files
#    No model or extension involved
#
## VARIABLES ##
#
# TREAT:   experiment SAM/BAM file, can be space delimited for replicates
# CONTROL: control SAM/BAM file, can be space delimited for replicates
# FOLDER:  the name of the folder where all of the MACS2 files will go
# M2FILE:  the tagname of the files created with this script
# BDG:     [optional] whether or not to make FE-normalized bedgraph files
#               default is F [ "TRUE", "FALSE", "RDNA" ]
# WIG:     [optional] whether or not to convert bedgraphs into wiggle files
#               BDG is assumed when WIG is not FALSE
#               default is F [ "TRUE", "FALSE", "RDNA" ]
# PEAK:    [optional] type of peaks to be called
#               default is "" (none) [ "", "NARROW", "BROAD", "BOTH" ]
# FLDR:    [optional] working directory
#               default is working directory when submitting
# PE:      [optional] whether input data should be treated as paired-end
#               sequencing data
#               default is F ["TRUE", "FALSE"]
#
####################################

# Go to the currect starting folder
cd /scratch/$USER
if [[ $FLDR != "" ]]; then
    cd $FLDR
else
    cd $SLURM_SUBMIT_DIR
fi
DIR=$PWD

echo "working in $DIR"
echo "Treatment file is: $TREAT"
echo "Control file is: $CONTROL"
echo "Output folder is: $FOLDER"
echo "Experiment ID is: $M2FILE"
echo ""

# ensure PEAK and NORM is capital
PEAK=$(echo $PEAK | tr '[:lower:]' '[:upper:]')

GENSIZE=1.2e7

# ensuring that BDG and WIG are defined
if [[ $BDG == "" ]]; then
    if [[ $WIG == "" ]]; then
	BDG="FALSE"
    else
	BDG=$WIG
    fi
fi

if [[ $BDG =~ ^[Rr] ]]; then
    DUP="all"
else
    DUP="auto"
fi

FOLDER=$FOLDER-MACS2

if [[ ! $M2FILE =~ "W4" ]]; then
    M2FILE="${M2FILE}_W4_MACS2"
else
    M2FILE="${M2FILE}_MACS2"
fi

# to deal with paired-end sequencing data
# MACS2 requires that data is in BAM format
# Paired-end ignores both --extsize, and all parameters relating to model-building
if [[ $PE =~ ^[Tt] ]]; then

    if [[ $TREAT == *.bam ]]; then
	echo "Warning: unable to ensure that this is a paired-end file."
    elif [[ $TREAT == *.sam ]]; then
	if [[ $TREAT =~ " " ]]; then
	    echo "Convert paired-end files into BAM format before running replicates through MACS2."
	    exit 1
	'''fi
	READ1=`tail -n 1 $TREAT | cut -f 1 | cut -f 1 -d " "`
	READ2=`tail -n 2 $TREAT | head -n 1 | cut -f 1 | cut -f 1 -d " "`
	if [[ $READ1 != $READ2 ]]; then
	    echo "ChIP file $TREAT does not contain paired-end data"
	    exit 1'''
	fi
    else
	echo "ChIP file $TREAT cannot be used for paired-end analysis"
	exit 1
    fi

    if [[ $CONTROL == *.bam ]]; then
	echo "Warning: unable to ensure that this is a paired-end file."
    elif [[ $CONTROL == *.sam ]]; then
	if [[ $CONTROL =~ " " ]]; then
	    echo "Convert paired-end files into BAM format before running replicates through MACS2."
	    exit 1
	'''fi
	READ3=`tail -n 1 $CONTROL | cut -f 1 | cut -f 1 -d " "`
	READ4=`tail -n 2 $CONTROL | head -n 1 | cut -f 1 | cut -f 1 -d " "`
	if [[ $READ3 != $READ4 ]]; then
	    echo "Input file $CONTROL does not contain paired-end data"
	    exit 1'''
	fi
    else
	echo "Input file $CONTROL cannot be used for paired-end analysis"
	exit 1
    fi

RE='(.+).sam'
if [[ $TREAT =~ $RE ]]; then
    TREAT1=${BASH_REMATCH[1]}.bam
    echo "Converting ChIP file to BAM"
    samtools view -b -h $TREAT > $TREAT1
    samtools sort -o $TREAT1 $TREAT1
    samtools index $TREAT1
fi

if [[ $CONTROL =~ $RE ]]; then
    echo "Converting Input file to BAM"
    CONTROL1=${BASH_REMATCH[1]}.bam
    samtools view -b -h $CONTROL > $CONTROL1
    samtools sort -o $CONTROL1 $CONTROL1
    samtools index $CONTROL1
fi

    FORMAT="BAMPE"
else
    FORMAT="AUTO"
    TREAT1=$TREAT
    CONTROL1=$CONTROL
fi

#### Run MACS2 and downstream

# use MACS2 to call peaks and create SPMR-normalized bedgraph files
if [[ $BDG =~ [fF] ]]; then
    # if wiggle files do not need to be created, don't make bedgraphs
    # currently not ensuring that files are not already created before running functions
    if [[ $PEAK == "NARROW" || $PEAK == "BOTH" ]]; then
	    macs2 callpeak -t $TREAT1 -c $CONTROL1 -n $M2FILE \
		--outdir $FOLDER -g $GENSIZE --extsize 200 \
		--SPMR --keep-dup=$DUP --nomodel -f $FORMAT
    fi
    if [[ $PEAK == "BROAD" || $PEAK == "BOTH" ]]; then
	macs2 callpeak -t $TREAT1 -c $CONTROL1 -n $M2FILE \
	    --outdir $FOLDER -g $GENSIZE --extsize 200 \
	    --SPMR --keep-dup=$DUP --nomodel --broad -f $FORMAT
    fi
else
    if [[ $PEAK == "NARROW" || $PEAK == "" ]]; then
	macs2 callpeak -t $TREAT1 -c $CONTROL1 -n $M2FILE \
	    --outdir $FOLDER -B -g $GENSIZE --extsize 200 \
	    --SPMR --keep-dup=$DUP --nomodel -f $FORMAT
    elif [[ $PEAK == "BROAD" ]]; then
	macs2 callpeak -t $TREAT1 -c $CONTROL1 -n $M2FILE \
	    --outdir $FOLDER -B -g $GENSIZE --extsize 200 \
	    --SPMR --keep-dup=$DUP --nomodel --broad -f $FORMAT
    elif [[ $PEAK == "BOTH" ]]; then
	macs2 callpeak -t $TREAT1 -c $CONTROL1 -n $M2FILE \
	    --outdir $FOLDER -B -g $GENSIZE --extsize 200 \
	    --SPMR --keep-dup=$DUP --nomodel -f $FORMAT
	macs2 callpeak -t $TREAT1 -c $CONTROL1 -n $M2FILE \
	    --outdir $FOLDER -g $GENSIZE --extsize 200 \
	    --SPMR --keep-dup=$DUP --nomodel --broad -f $FORMAT
    else
	echo "Do not recognize method chosen for PEAK"
	exit 1
    fi

    echo ""
    echo "working in $FOLDER"
    cd $FOLDER

    # get fold-enrichment normalized bedgraphs
    macs2 bdgcmp -t ${M2FILE}_treat_pileup.bdg \
	-c ${M2FILE}_control_lambda.bdg \
	-o ${M2FILE}_FE.bdg -m FE

    if [[ $WIG =~ ^[TtRr] ]]; then
	# convert bedgraphs created by MACS2 into wiggle files
	if [ ! -f ~/ChIPseq_Pipeline_v3/Bedgraph2VariableStepWiggle.py ]; then
	    echo "Cannot find Bedgraph2VariableStepWiggle.py"
	    echo "This file should be located in ~/ChIPseq_Pipeline_v3"
	    exit 1
	else
	    if [ -d ${M2FILE}_FE ]; then
		echo "Already converted ${M2FILE}_FE.bdg into wiggle"
	    else
		module purge
		module load python3/intel/3.5.3

		echo "Converting ${M2FILE}_FE.bdg into wiggle"
		python3 ~/ChIPseq_Pipeline_v3/Bedgraph2VariableStepWiggle.py \
		    -b ${M2FILE}_FE.bdg

		cd $DIR/$FOLDER
		if [ -d ${M2FILE}_FE ]; then
		    gzip ${M2FILE}_FE/*.wig
		else
		    gzip $M2FILE*.wig
		fi
	    fi
	fi
    fi
    gzip $M2FILE*.bdg
fi

echo "complete"
date
echo ""
