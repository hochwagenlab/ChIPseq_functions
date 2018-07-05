#!/bin/bash
# written by: Tovah Markowitz

# known issues:
# 1) no checks for the following: 
# -- do all SAM files have the same root?
# -- have wiggle files/peak files already been created?
# -- is filemaker ID already in the name? (only works if they are identical)
# 2) memory issue
# 3) three error/output files created at the end of each job

###########################################################
### FUNCTION1

function run_closing {
# run sbatch script ChIPseq_closing.sbatch to organize error/output files
    TMPID=$1
    JIDS=$2
    TAGOUT=$3

    if [[ $JIDS != "" ]]; then
	sbatch --dependency=afterany:$JIDS \
	    --export IDS=$JIDS,TMPID=$TMPID,TAGOUT=$TAGOUT \
            ~/ChIPseq_Pipeline_v3/ChIPseq_closing.sbatch
    fi
}

###########################################################
### FUNCTION2

function check_fastq {
    # check input variables to determine if they are to be run through Bowtie
    # purpose: to count how many Bowtie scripts to run
    # POTFQ = potential fastq [file to check]
    # returns: FQTF [FastqT/F] as "TRUE/FALSE"
    TMPID=$1
    TAGOUT=$2
    POTFQ=$3
    TAGN=$4
    JIDS=$5

    if [[ $POTFQ =~ .fastq|.fq ]]; then
	if [[ $TAGN == "" ]]; then
	    echo "$POTFQ was not given a tagname" >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
	    run_closing $TMPID $JIDS $TAGOUT
	    exit 1
	else
	    FQTF="TRUE"
	fi
    else
	FQTF="FALSE"
    fi
}

###########################################################
### FUNCTION3

function assign_genroot {
    # adjusting output names to fit with genome requested
    # returns: GENROOT and GENROOT2
    TMPID=$1
    TAGOUT=$2
    GEN=$3
    PE=$4
    JIDS=$5

    if [[ $PE =~ [Tt] ]]; then
	case $GEN in
            SK1K ) GENROOT="SK1K-PM-PE_B4"
		   GENROOT2="SK1K-rDNA-PE_B4"
		   ;;
            SacCer3 ) GENROOT="SacCer3-2mis-PE_B4"
	              GENROOT2="SacCer3-rDNA-PE_B4"
		      ;;
	    SK1Yue-PM ) GENROOT="SK1Yue-PM-PE_B4"
			GENROOT2="SK1Yue-rDNA-PE_B4"
			;;
	    SK1Yue-2mis ) GENROOT="SK1Yue-2mis-PE_B4"
			  GENROOT2="SK1Yue-rDNA-PE_B4"
			  ;;
            * ) echo "Unknown genome listed." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
		run_closing $TMPID  $JIDS $TAGOUT
		exit 1
	esac
    else
        case $GEN in
            SK1K ) GENROOT="SK1K-PM_B3"
		   GENROOT2="SK1K-rDNA_B3"
		   ;;
            SacCer3 ) GENROOT="SacCer3-2mis_B3"
	              GENROOT2="SacCer3-rDNA_B3"
		      ;;
	    SK1Yue-PM ) GENROOT="SK1Yue-PM_B3"
			GENROOT2="SK1Yue-rDNA_B3"
			;;
	    SK1Yue-2mis ) GENROOT="SK1Yue-2mis_B3"
			  GENROOT2="SK1Yue-rDNA_B3"
			  ;;
            * ) echo "Unknown genome listed." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
		run_closing $TMPID $JIDS $TAGOUT
		exit 1
	esac
    fi
}

###########################################################
### FUNCTION4

function parse_inCHIPname {
    # if running MACS and not Bowtie and need Genroot from files
    TMPID=$1
    TAGOUT=$2
    TREATMENT1=$3
    JIDS=$4

    # if there are spaces in CHIP, it must be replicates
    # continue to parse first file for ROOT
    RESPACE="(^[^[:space:]]*)[[:space:]].*"
    if [[ $TREATMENT1 =~ $RESPACE ]]; then
        TREATMENT=${BASH_REMATCH[1]}
    else
        TREATMENT=$TREATMENT1
    fi

    # parse filename of first ChIP file
    if [[ $TREATMENT == */* ]]; then
        RE='.*/(.+)-(S.+).(sam|bam)'
    else
        RE='(.+)-(S.+).(sam|bam)'
    fi
    if [[ $TREATMENT =~ $RE ]]; then
        TAGC=${BASH_REMATCH[1]}
        GENROOT=${BASH_REMATCH[2]}
    else
        echo "Error: ChIP filename ($TREATMENT) is incorrect." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
        echo "Either it is not labelled as a .sam or .bam " >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
        echo "or it is not in the correct format for this command." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
        echo "Correct pattern is (.+)-(S.+).(sam|bam)." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
        echo "Make sure the file has one and only one -S." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt$OUT
        echo "Both relative and complete paths are accepted." >> /scratch/$USER/MACS2_pipeline_$TMPID.txt
        run_closing $TMPID $JIDS $TAGOUT
        exit 1
    fi
}

###########################################################
### FUNCTION5

function change_genroot {
    # adjust genroot for MACS2
    # create ROOTF (for the files) and ROOTFR (for the outer folder)
    # shortening folder names for MACS2 outputs
    ROOT=$1

if [[ $ROOT =~ "B3"$ ]]; then
    ROOTF=$ROOT
    VER="-B3"
elif [[ $ROOT =~ "B4"$ ]]; then
    ROOTF=$ROOT
    VER="-B4"
elif [[ $ROOT =~ "2mis-PM" ]]; then
    ROOTF="SacCer3-2mis_"
else
    ROOTF="${ROOT}_"
fi

if [[ $ROOT =~ "SK1K" ]]; then 
    ROOTFR="SK1K$VER"
elif [[ $ROOT =~ "SacCer3" ]]; then
    ROOTFR="SacCer3$VER"
elif [[ $ROOT =~ "SK1Yue" ]]; then
    ROOTFR="SK1Yue$VER"
fi
}

###########################################################
### FUNCTION6

function define_MACS2_filenames {
    IN=$1
    ROOTF=$2
    ROOTFR=$3
    TYPE=$4
    FLMKR=$5

    if [[ $FLMKR != "" ]]; then
	if [[ ! $IN =~ $FLMKR ]]; then
	    FLMKR2="-$FLMKR"
	fi
    fi

    TP2=""
    if [[ $TYPE == "NORMAL" ]]; then
	TP2=""
    elif [[ $TYPE =~ "REPS" ]]; then
	TP2="-Reps"
    fi

   
    if [[ $TYPE =~ "TAGNORM" ]]; then
	M2FILE1=$IN$TP2-ChvCh$FLMKR2-${ROOTF}W4
	FOLDER=$IN$TP2-ChvCh$FLMKR2-${ROOTFR}W4
    else
	M2FILE=$IN$FLMKR2$TP2-${ROOTF}W4
	FOLDER=$IN$FLMKR2$TP2-${ROOTFR}W4
    fi
}

