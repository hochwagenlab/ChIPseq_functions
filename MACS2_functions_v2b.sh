#!/bin/bash

# written by: Tovah Markowitz
# Date: 4/11/16
# updated: 4/20/16 to change bdg/wig normalization from SPMR to median
# updated: 9/27/16 to include checks to ensure dependency files are found
# updated: 10/4/16 to fix issues with double normalization
# updated: 11/13/16 to work with files on github

##########################################
### FUNCTION1
### For replicates
# Inputs: INPUT1, INPUT2, OUTROOT

function replicate_merge {
    INPUT1=$1
    INPUT2=$2
    OUTROOT=$3

    module load samtools/intel/1.3
    
    test -d RepBam || mkdir RepBam

    echo "Running replicate_merge on $INPUT1 and $INPUT2"
    echo "Outroot name is $OUTROOT"
    date

# split INPUT1 and INPUT2 file names into pieces
## Pattern is looking to split into tagname, rootname and file type
## Tagname should be specific to a specific library
## Rootname must begin with -S and should be identical between files
## Files must be sam or bam files
if [[ $INPUT1 == */* ]]; then
    RE1='.*/(.+)-(S.+).(sam|bam)'
else
    RE1='(.+)-(S.+).(sam|bam)'
fi
if [[ $INPUT2 == */* ]]; then
    RE2='.*/(.+)-(S.+).(sam|bam)'
else
    RE2='(.+)-(S.+).(sam|bam)'
fi

    if [[ $INPUT1 =~ $RE1 ]]; then
	IN1=${BASH_REMATCH[1]}
	ROOT1=${BASH_REMATCH[2]}
	TYPE1=${BASH_REMATCH[3]}
    else
	echo "Error: First input name ($INPUT1) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

if [[ $INPUT2 =~ $RE2 ]]; then
	IN2=${BASH_REMATCH[1]}
	ROOT2=${BASH_REMATCH[2]}
	TYPE2=${BASH_REMATCH[3]}
    else
	echo "Error: Second input name ($INPUT2) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $ROOT1 != $ROOT2 ]]; then
	echo "Filenames ($INPUT1 and $INPUT2) do not have the same root names and may not have been compared to the same genome scaffold or software."
	exit 1
    fi

# Before running this function, make sure output bam file doesn't already exist
if [ -f RepBam/$OUTROOT-reps-${ROOT1}.bam ]; then
    echo "RepBam/$OUTROOT-reps-${ROOT1}.bam already exists. Continuing."
else

# determine if first file is sam or bam, then convert into sorted bam
    if [[ $TYPE1 == sam ]]; then
	samtools view -b -h $INPUT1 > RepBam/$IN1-$ROOT1.bam
	echo "sorting RepBam/$IN1-$ROOT1.bam"
	samtools sort -o RepBam/$IN1-$ROOT1.bam RepBam/$IN1-$ROOT1.bam
    else
	echo "Sorting $INPUT1"
	samtools sort -o RepBam/$IN1-$ROOT1.bam $INPUT1
    fi

# determine if second file is sam or bam, then convert into sorted bam
    if [[ $TYPE2 == sam ]]; then
	samtools view -b -h $INPUT2 > RepBam/$IN2-$ROOT1.bam
	echo "sorting RepBam/$IN2-$ROOT2.bam"
	samtools sort -o RepBam/$IN2-$ROOT1.bam RepBam/$IN2-$ROOT1.bam
    else
	echo "Sorting $INPUT2"
	samtools sort -o RepBam/$IN2-$ROOT1.bam $INPUT2
    fi

# merge files
    if [[ "RepBam/$OUTROOT-reps-${ROOT1}.bam" == "RepBam/$IN2-$ROOT1.bam" ]]; then
	echo "temporarily adjusting files"
	TMPID=`date +%N`
	mv RepBam/$IN2-$ROOT1.bam RepBam/tmp-$TMPID.bam
	BAM2="RepBam/tmp-$TMPID.bam"
    else
	BAM2="RepBam/$IN2-$ROOT1.bam"
    fi

    samtools merge RepBam/$OUTROOT-reps-${ROOT1}.bam RepBam/$IN1-$ROOT1.bam $BAM2
    echo "wrote RepBam/$OUTROOT-reps-${ROOT1}.bam"

    if [ -e RepBam/tmp*.bam ]; then
	rm RepBam/tmp*.bam
    fi
fi
    echo "complete"
    date
    echo ""
}

#########################################
### FUNCTION2
### get normalized wiggle from bedgraph
# Inputs: TREAT, CONTROL, MLOW
# TREAT/CONTROL are SAM or BAM files

function create_bedgraph {
    TREAT=$1
    CONTROL=$2
    MLOW=$3

    module load macs2/intel/2.1.0
    module load r/intel/3.2.2
    echo "running create_bedgraph on $TREAT and $CONTROL"
    date

# split TREAT and CONTROL file names into pieces
## Pattern is looking to split into tagname, rootname and file type
## Tagname should be specific to a specific library
## Rootname must begin with -S and should be identical between files
## Files must be sam or bam files
if [[ $TREAT == */* ]]; then
    RE1='/(.+)-(S.+).(sam|bam)'
else
    RE1='(.+)-(S.+).(sam|bam)'
fi
if [[ $CONTROL == */* ]]; then
    RE2='/(.+)-(S.+).(sam|bam)'
else
    RE2='(.+)-(S.+).(sam|bam)'
fi

    if [[ $TREAT =~ $RE1 ]]; then
	IN1=${BASH_REMATCH[1]}
	ROOT1=${BASH_REMATCH[2]}
    else
	echo "Error: Treatment filename ($TREAT) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $CONTROL =~ $RE2 ]]; then
	IN2=${BASH_REMATCH[1]}
	ROOT2=${BASH_REMATCH[2]}
    else
	echo "Error: Control filename ($CONTROL) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $ROOT1 != $ROOT2 ]]; then
	echo "Filenames ($TREAT and $CONTROL) do not have the same root names and may not have been compared to the same genome scaffold or software."
	exit 1
    else
	TAGT=$IN1-$ROOT1
    fi

# check MLOW
    if [[ $MLOW == "" ]]; then
	echo "Using default mfold values."
	MLOW=5
    elif [[ $MLOW =~ [0-9]* ]]; then
	if [[ $MLOW -gt 49 ]]; then
	    echo "MLOW ($MLOW) is higher than the maximum mfold value [50]."
	    exit 1
	fi
    else
	echo "MLOW ($MLOW) is not numeric."
	exit 1
    fi
    
    if [[ $ROOT1 =~ "reps" ]]; then
	IN=$IN1-reps
    else
	IN=$IN1
    fi
    
    STDIR=$PWD
    
# make bedgraphs for individual files
    if [ -e $IN-MACS2/$TAGT-M${MLOW}_MACS2_treat_pileup.bdg ]; then
	echo "already created necessary bedgraph in MACS2"
    else
	macs2 callpeak -t $TREAT -c $CONTROL -n $TAGT-M${MLOW}_MACS2 --outdir $IN-MACS2 -g 1.2e7 --bw=350 --keep-dup="auto" --mfold $MLOW 50 -B --fix-bimodal
    fi

if [ ! -f ~/Pipeline/Bedgraph2VariableStepWiggle.py ]; then
    echo "Cannot find Bedgraph2VariableStepWiggle.py. Should be in ~/Pipeline/"
    exit 1
fi

# convert to wiggle
    cd $IN-MACS2
    if [ -d $TAGT-M${MLOW}_MACS2_control_lambda ]; then
	echo "already converted $TAGT-M${MLOW}_MACS2_control_lambda.bdg into wiggle"
    else
	echo "converting $TAGT-M${MLOW}_MACS2_control_lambda.bdg to wiggle"
	python ~/Pipeline/Bedgraph2VariableStepWiggle.py -b $TAGT-M${MLOW}_MACS2_control_lambda.bdg
    fi
    if [ -d $TAGT-M${MLOW}_MACS2_treat_pileup ]; then
	echo "already converted $TAGT-M${MLOW}_MACS2_treat_pileup.bdg into wiggle"
    else
	echo "converting $TAGT-M${MLOW}_MACS2_treat_pileup.bdg to wiggle"
	python ~/Pipeline/Bedgraph2VariableStepWiggle.py -b $TAGT-M${MLOW}_MACS2_treat_pileup.bdg
    fi

if [ ! -f ~/Pipeline/NormalizeWiggle.R ]; then
    echo "Cannot find NormalizeWiggle.R. Should be in ~/Pipeline/"
    exit 1
fi

    OUTFILE="$TAGT-M${MLOW}_MACS2_wiggle_norm"
    if [ -d  $OUTFILE ]; then
	echo "already normalized $OUTFILE"
    else
	R CMD BATCH "--args $TAGT-M${MLOW}_MACS2_treat_pileup $TAGT-M${MLOW}_MACS2_control_lambda" ~/Pipeline/NormalizeWiggle.R
    fi

    cd $STDIR
    mv $IN-MACS2/NormalizeWiggle.Rout $IN-MACS2/$TAGT-M${MLOW}_MACS2_NormalizeWiggle.Rout
    gzip $IN-MACS2/*/*.wig
    
    echo "Extracting median values from 'NormalizeWiggle.Rout' and saving in ~/Genomewide_medians_MACS2.txt'"
    perl ~/Pipeline/ExtractMedian_MACS2_pipeline.pl -f $IN-MACS2 -u $USER
    echo "complete"
    date
    echo ""
}

#########################################
### FUNCTION3
### get median normalized wiggle file from bedgraph (nomodel)
# Inputs: TREAT, CONTROL
# TREAT/CONTROL are SAM or BAM files

function create_bedgraph_nomodel {
    TREAT=$1
    CONTROL=$2
   
    module load macs2/intel/2.1.0
    module load r/intel/3.2.2
    echo "running create_bedgraph_nomodel on $TREAT and $CONTROL"
    date
    
 # split TREAT and CONTROL file names into pieces
## Pattern is looking to split into tagname, rootname and file type
## Tagname should be specific to a specific library
## Rootname must begin with -S and should be identical between files
## Files must be sam or bam files
if [[ $TREAT == */* ]]; then
    RE1='/(.+)-(S.+).(sam|bam)'
else
    RE1='(.+)-(S.+).(sam|bam)'
fi
if [[ $CONTROL == */* ]]; then
    RE2='/(.+)-(S.+).(sam|bam)'
else
    RE2='(.+)-(S.+).(sam|bam)'
fi

    if [[ $TREAT =~ $RE1 ]]; then
	IN1=${BASH_REMATCH[1]}
	ROOT1=${BASH_REMATCH[2]}
    else
	echo "Error: Treatment filename ($TREAT) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $CONTROL =~ $RE2 ]]; then
	IN2=${BASH_REMATCH[1]}
	ROOT2=${BASH_REMATCH[2]}
    else
	echo "Error: Control filename ($CONTROL) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $ROOT1 != $ROOT2 ]]; then
	echo "Filenames ($TREAT and $CONTROL) do not have the same root names and may not have been compared to the same genome scaffold or software."
	exit 1
    else
	TAGT=$IN1-$ROOT1
    fi

    if [[ $ROOT1 =~ "reps" ]]; then
	IN=$IN1-reps
    else
	IN=$IN1
    fi
    
    STDIR=$PWD
# make bedgraphs for individual files
    if [ -e $IN-MACS2/${TAGT}_MACS2_treat_pileup.bdg ]; then
	echo "already created necessary bedgraph in MACS2"
    else
	macs2 callpeak -t $TREAT -c $CONTROL -n $TAGT_MACS2 --outdir $IN-MACS2 -g 1.2e7 --bw=350 --keep-dup="auto" -B --nomodel
    fi

if [ ! -f ~/Pipeline/Bedgraph2VariableStepWiggle.py ]; then
    echo "Cannot find Bedgraph2VariableStepWiggle.py. Should be in ~/Pipeline/"
    exit 1
fi

# convert to wiggle
    cd $IN-MACS2
    if [ -d ${TAGT}_MACS2_control_lambda ]; then
	echo "already converted ${TAGT}_MACS2_control_lambda.bdg into wiggle"
    else
	echo "converting ${TAGT}_MACS2_control_lambda.bdg to wiggle"
	python ~/Pipeline/Bedgraph2VariableStepWiggle.py -b ${TAGT}_MACS2_control_lambda.bdg
    fi
    if [ -d ${TAGT}_MACS2_treat_pileup ]; then
	echo "already converted ${TAGT}_MACS2_treat_pileup.bdg into wiggle"
    else
	echo "converting ${TAGT}_MACS2_treat_pileup.bdg to wiggle"
	python ~/Pipeline/Bedgraph2VariableStepWiggle.py -b ${TAGT}_MACS2_treat_pileup.bdg
    fi

if [ ! -f ~/Pipeline/NormalizeWiggle.R ]; then
    echo "Cannot find NormalizeWiggle.R. Should be in ~/Pipeline/"
    exit 1
fi

    OUTFILE="${TAGT}_MACS2_wiggle_norm"
    if [ -d  $OUTFILE ]; then
	echo "already normalized $OUTFILE"
    else
	R CMD BATCH "--args ${TAGT}_MACS2_treat_pileup ${TAGT}_MACS2_control_lambda" ~/Pipeline/NormalizeWiggle.R
    fi

    cd $STDIR
    mv $IN-MACS2/NormalizeWiggle.Rout $IN-MACS2/${TAGT}_MACS2_NormalizeWiggle.Rout
    gzip $IN-MACS2/*/*.wig
    perl ~/Pipeline/ExtractMedian_MACS2_pipeline.pl -f $IN-MACS2 -u $USER


    echo "complete"
    date
    echo ""
}

#########################################
### FUNCTION4
### call narrow and broad peaks
# Inputs: TREAT, CONTROL, MLOW
# TREAT/CONTROL are SAM or BAM files

function peak_call {
    TREAT=$1
    CONTROL=$2
    MLOW=$3
    DNORM=$4
    
    module load macs2/intel/2.1.0
    
    echo "running peak_call on $TREAT and $CONTROL"
    date
    
# split TREAT and CONTROL file names into pieces
## Pattern is looking to split into tagname, rootname and file type
## Tagname should be specific to a specific library
## Rootname must begin with -S and should be identical between files
## Files must be sam or bam files
if [[ $TREAT == */* ]]; then
    RE1='/(.+)-(S.+).(sam|bam)'
else
    RE1='(.+)-(S.+).(sam|bam)'
fi
if [[ $CONTROL == */* ]]; then
    RE2='/(.+)-(S.+).(sam|bam)'
else
    RE2='(.+)-(S.+).(sam|bam)'
fi
    
    if [[ $TREAT =~ $RE1 ]]; then
	IN1=${BASH_REMATCH[1]}
	ROOT1=${BASH_REMATCH[2]}
    else
	echo "Error: Treatment filename ($TREAT) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi
    
    if [[ $CONTROL =~ $RE2 ]]; then
	IN2=${BASH_REMATCH[1]}
	ROOT2=${BASH_REMATCH[2]}
    else
	echo "Error: Control filename ($CONTROL) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi
    
    if [[ $ROOT1 != $ROOT2 ]]; then
	echo "Filenames ($TREAT and $CONTROL) do not have the same root names and may not have been compared to the same genome scaffold or software."
	exit 1
    else
	TAGT=$IN1-$ROOT1
    fi
    
# check MLOW
    if [[ $MLOW == "" ]]; then
	echo "Using default mfold values."
	MLOW=5
    elif [[ $MLOW =~ [0-9]* ]]; then
	if [[ $MLOW -gt 49 ]]; then
	    echo "MLOW ($MLOW) is higher than the maximum mfold value [50]."
	    exit 1
	fi
    else
	echo "MLOW ($MLOW) is not numeric."
	exit 1
    fi
    
    if [[ $DNORM == "" ]]; then
	NAME=$TAGT-M${MLOW}_MACS2
    else
	NAME=$TAGT-TagNorm-M${MLOW}_MACS2
    fi

    if [[ $ROOT1 =~ "reps" ]]; then
	IN=$IN1-reps
    else
	IN=$IN1
    fi
    
    macs2 callpeak -t $TREAT -c $CONTROL -n $NAME --outdir $IN-MACS2 -g 1.2e7 --bw=350 --keep-dup="auto" --mfold $MLOW 50 --fix-bimodal
    
    macs2 callpeak -t $TREAT -c $CONTROL -n $NAME --outdir $IN-MACS2 -g 1.2e7 --bw=350 --keep-dup="auto" --mfold $MLOW 50 --fix-bimodal --broad
    
    echo "complete"
    date
    echo ""
}

##########################################
### FUNCTION5
### call narrow and broad peaks (nomodel)
# Inputs: TREAT, CONTROL
# TREAT/CONTROL are SAM or BAM files

function peak_call_nomodel {
    TREAT=$1
    CONTROL=$2
    DNORM=$3

    module load macs2/intel/2.1.0
    
    echo "running peak_call_nomodel on $TREAT and $CONTROL"
    date

    # split TREAT and CONTROL file names into pieces
## Pattern is looking to split into tagname, rootname and file type
## Tagname should be specific to a specific library
## Rootname must begin with -S and should be identical between files
## Files must be sam or bam files
if [[ $TREAT == */* ]]; then
    RE1='/(.+)-(S.+).(sam|bam)'
else
    RE1='(.+)-(S.+).(sam|bam)'
fi
if [[ $CONTROL == */* ]]; then
    RE2='/(.+)-(S.+).(sam|bam)'
else
    RE2='(.+)-(S.+).(sam|bam)'
fi

    if [[ $TREAT =~ $RE1 ]]; then
	IN1=${BASH_REMATCH[1]}
	ROOT1=${BASH_REMATCH[2]}
    else
	echo "Error: Treatment filename ($TREAT) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $CONTROL =~ $RE2 ]]; then
	IN2=${BASH_REMATCH[1]}
	ROOT2=${BASH_REMATCH[2]}
    else
	echo "Error: Control filename ($CONTROL) is incorrect."
	echo "Either it is not labelled as a .sam or .bam or it is not in the correct format for this command."
	exit 1
    fi

    if [[ $ROOT1 != $ROOT2 ]]; then
	echo "Filenames ($TREAT and $CONTROL) do not have the same root names and may not have been compared to the same genome scaffold or software."
	exit 1
    else
	TAGT=$IN1-$ROOT1
    fi

    if [[ $DNORM == "" ]]; then
	NAME=${TAGT}_MACS2
    else
	NAME=$TAGT-TagNorm_MACS2
    fi

    if [[ $ROOT1 =~ "reps" ]]; then
	IN=$IN1-reps
    else
	IN=$IN1
    fi

    macs2 callpeak -t $TREAT -c $CONTROL -n $NAME --outdir $IN-MACS2 -g 1.2e7 --bw=350 --keep-dup="auto" --no-model

    macs2 callpeak -t $TREAT -c $CONTROL -n $NAME --outdir $IN-MACS2 -g 1.2e7 --bw=350 --keep-dup="auto" --no-model --broad

    echo "complete"
    date
    echo ""
}

##########################################
### FUNCTION6
### double normalization of wiggle files
# Inputs: TREAT, CONTROL
# Assumes both wiggle files came from this pipeline and have already been normalized once

function wiggle_dnorm {
    TREAT=$1
    CONTROL=$2
    FLMKRIDS=$3
    
    module load macs2/intel/2.1.0

    echo "running bedgraph_norm on $TREAT and $CONTROL"
    date

# split TREAT and CONTROL folder names into two pieces to work with normalizewig_batch.R output
    RE_wig="(.+)-S(.+)_MACS2_wiggle_norm"
    if [[ $TREAT =~ $RE_wig ]]; then
	IN1=${BASH_REMATCH[1]}
    else
	echo "Error with double norm file name $TREAT"
	exit 1
    fi
    if [[ $CONTROL =~ $RE_wig ]]; then
	IN2=${BASH_REMATCH[1]}
    else
	echo "Error with double norm file name $CONTROL"
	exit 1
    fi
   
    STDIR=$PWD
    
    echo "ChIP folder is $TREAT"
    echo "Control folder is $CONTROL"
    echo "Filemaker IDs are $FLMKRIDS"

    gzip -d $TREAT/*.wig.gz
    gzip -d $CONTROL/*.wig.gz

if [ ! -f ~/Pipeline/NormalizeWiggle.R ]; then
    echo "Cannot find NormalizeWiggle.R. Should be in ~/Pipeline/"
    exit 1
fi

    R CMD BATCH "--args $TREAT $CONTROL double_norm NA NA $FLMKRIDS" ~/Pipeline/NormalizeWiggle.R

    cd $STDIR
    gzip $CONTROL/*.wig
    gzip $TREAT/*.wig
#    gzip *double_norm/*.wig

    mv NormalizeWiggle.Rout ${TREAT}_doublenorm_NormalizeWiggle.Rout

    echo "complete"
    date
    echo ""

}

##########################################
### FUNCTION7
### bowtie alignment- SK1
# Inputs: FASTQ, TAG

function bowtie_SK1 {
    FASTQIN=$1
    TAG=$2

    module load bowtie/intel/1.1.1
    echo "running bowtie_SK1 on $FASTQIN with tagname $TAG"
    date
    
    RE="(.+)\.gz"

    if [[ $FASTQIN =~ $RE ]]; then
	gzip -d $FASTQIN
	FASTQ=${BASH_REMATCH[1]}
    else
	FASTQ=$FASTQIN
    fi
	
    test -d Bowtie || mkdir Bowtie
    test -d Unaligned || mkdir Unaligned
    
    # only perfect matches
    bowtie -q -m 1 -v 0 -S --un Unaligned/$TAG-SK1K-Unaligned.fastq --max Unaligned/$TAG-SK1K-Max.fastq ~/Library/sk1_MvO_V1 $FASTQ Bowtie/$TAG-SK1K-PM.sam

    gzip Unaligned/*.fastq
    gzip $FASTQ

    echo "complete"
    date
    echo ""
}

##########################################
### FUNCTION8
### bowtie alignment- SacCer3
# Inputs: FASTQ, TAG

function bowtie_SacCer3 {
    FASTQIN=$1
    TAG=$2
    
    module load bowtie/intel/1.1.1
    echo "running bowtie_SacCer3 on $FASTQIN with tagname $TAG"
    date
    
    RE="(.+)\.gz"

    if [[ $FASTQIN =~ $RE ]]; then
	gzip -d $FASTQIN
	FASTQ=${BASH_REMATCH[1]}
    else
	FASTQ=$FASTQIN
    fi
    
    test -d Bowtie || mkdir Bowtie
    test -d Unaligned || mkdir Unaligned

    # allow 2 mismatches, only map to one location
    bowtie -q -m 1 -S --un Unaligned/$TAG-SacCer3-2mis-Unaligned.fastq --max Unaligned/$TAG-SacCer3-2mis-Max.fastq ~/Library/SacCer3 $FASTQ Bowtie/$TAG-SacCer3-2mis.sam

    gzip Unaligned/*.fastq
    gzip $FASTQ
    
    echo "complete"
    date
    echo ""
}
