#! /bin/bash
##SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --mem=10MB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tem298@nyu.edu
#SBATCH --time=5:00
#SBATCH --job-name=ChIPseq_Pipeline_prince

module purge

#################################################################################################################
## VARIABLES                                                                                                    #
#                                                                                                               #
# MAIN VARIABLES                                                                                                #
#   CHIP     -- Only required feature. May be SAM, BAM, or FASTQ file. Relative or absolute path accepted.      #
#               See below for specialized cases (bowtie alignment, replicates, or ChIP-vs-ChIP normalization )  #
#   INPUT    -- May be SAM, BAM, or FASTQ file. See below for specialized cases.                                #
#                                                                                                               #
# GENERIC VARIABLES                                                                                             #
## To run Bowtie without running MACS2, set PEAK="" and WIG="F"                                                 #
#   FLDR      -- Input folder name [ default: working directory when submitting job ]                           #
#                Location where output files will reside                                                        #
#   PEAK      -- Type of peaks to be called. Options are: "", "NARROW", "BROAD", or "BOTH" [ default: "" ]      #
#   WIG       -- Whether or not to create wiggle files. Options are: "TRUE" or "FALSE" [ default: "T" ]         #
#   FLMKR     -- filemaker IDs from Datasets of ChIP samples in "-"-delimited and quoted list                   #
#                Should not be strain IDS. Will be incorporated into all MACS2-related files.                   #
#                REQUIRED for replicates and ChIP-vs-ChIP normalization.                                        #
#                                                                                                               #
# CHIP vs CHIP NORMALIZATION                                                                                    #
## For mock controls, no tag controls, or histones                                                              #
## Note: for MACS2 analysis, CHIP,INPUT,INCON, and CHCON must be included to run                                #
#   INCON     -- control input file (see INPUT or CHIP for details)                                             #
#   CHCON     -- control ChIP file (see INPUT or CHIP for details)                                              #
#                                                                                                               #
# REPLICATES                                                                                                    #
#   REP       -- tagname for MACS2 outputs of replicates                                                        #
#   CHIP(etc) -- all included files must be SAM/BAM files [currently no error checking].                        #
#                CHIP,INPUT,INCON, and CHCON will all accept replicates as a " "-delimited list                 #
#                                                                                                               #
# BOWTIE ALIGNMENT                                                                                              #
## Available for CHIP,INPUT,INCON, and CHCON [Fastq files can, but do not have to be zipped.]                   #
## Cannot be combined with replicates                                                                           #
#   TAGI      -- Tagname for INPUT fastq                                                                        #
#   TAGC      -- Tagname for CHIP fastq                                                                         #
#   TAGICON   -- Tagname for INCON fastq                                                                        #
#   TAGCCON   -- Tagname for CHCON fastq                                                                        #
#   GEN       -- Genome to align to. Current options: SK1 or SacCer3                                            #
#                                                                                                               #
# SOME EXAMPLES                                                                                                 #
## Example 1: Start with input file A.fastq and ChIP file B.fastq, map to SacCer3, and create wiggle files      #
#             and both narrow and broad peak files (will also create wiggle plot of the rDNA)                   #
#      sbatch --export INPUT=A.fastq,CHIP=B.fastq,GEN="SacCer3",TAGI="A",TAGC="B",PEAK="BOTH" \                 #
#             ~/ChIPseq_Pipeline/ChIPseq-pipeline_v3.sbatch                                                     #
#                                                                                                               #
## Example 2: Get broad peak files for replicate input files A_1.sam and A_2.sam and                            #
#             ChIP files B_1.sam and B_2.sam                                                                    #
#      sbatch --export INPUT="A_1.sam A_2.sam",CHIP="B_1.sam B_2.sam",REP="B",FLMKR="1-3",PEAK="BROAD",WIG="F" \#
#             ~/ChIPseq_Pipeline/ChIPseq-pipeline_v3.sbatch                                                     #
#                                                                                                               #
## Example 3: Get normalized wiggle files (and narrow peaks) with ChIP file A1.sam, input file A2.sam,          #
#             mock ChIP file B1.sam, and mock input file B2.sam                                                 #
#      sbatch --export INPUT=A2.sam,CHIP=A1.sam,INCON=B2.sam,CHCON=B1.sam,FLMKR="1-2" \                         #
#             ~/ChIPseq_Pipeline/ChIPseq-pipeline_v3.sbatch                                                     #
#                                                                                                               #
#################################################################################################################               

# check for shell functions
if [ -f ~/ChIPseq_Pipeline_v3/ChIPseq_functions_prince.sh ]; then
    source ~/ChIPseq_Pipeline_v3/ChIPseq_functions_prince.sh
else
    echo "Cannot find MACS2_functions_prince.sh"
    exit 3
fi

# TMPID should be the root of the ChIP file (without any periods)
TMP=${CHIP##*/}
TMPID=${TMP%%.*}
OUT=/scratch/$USER/MACS2_pipeline_$TMPID.txt
touch $OUT

# change to folder
if [[ $FLDR == "" ]]; then
    cd $PWD
else
    cd /scratch/${USER}
    cd $FLDR
fi
echo "Working directory is specified as: $PWD" >> $OUT
jids="" # all jobs
jidBs="" # Bowtie jobs
jidRs="" # rDNA Bowtie jobs

FLMKR=$(echo $FLMKR | tr ':' '-')

############## BOWTIE ALIGNMENT ################
############### preparation ####################
if [[ $GEN != "" ]]; then
    echo "Bowtie mapping to $GEN requested." >> $OUT

    # ensure that there are no replicates
    if [[ $CHIP =~ " " || $INPUT =~ " " || \
	$CHCON =~ " " || $INCON =~ " " ]]; then
	echo "Error: spaces included for variable CHIP." >> $OUT
	echo "Cannot have spaces and run Bowtie alignment for the same variable." >> $OUT
	echo "If file has spaces, use rm_space_filename.py" >> $OUT
	run_closing $TMPID $jids $TAGC
	exit 1
    fi

    # fqCount is an array which denotes which of the four variables 
    # have fastq files (CHIP, INPUT, INCON, CHCON)
    # make sure that replicates are not involved (check files for spaces)
    fqCount=()

    check_fastq $OUT $CHIP $TAGC $jids
    if [[ $FQTF == "TRUE" ]]; then fqCount+=(1); fi

    check_fastq $OUT $INPUT $TAGI $jids
    if [[ $FQTF == "TRUE" ]]; then fqCount+=(2); fi

    check_fastq $OUT $CHCON $TAGCCON $jids
    if [[ $FQTF == "TRUE" ]]; then fqCount+=(3); fi

    check_fastq $OUT $INCON $TAGICON $jids
    if [[ $FQTF == "TRUE" ]]; then fqCount+=(4); fi

    # adjusting output names to fit with genome requested
    assign_genroot $OUT $GEN $jids

################# setting up bowtie sbatch ###############
    for fq in ${fqCount[@]}
    do
	if [[ $fq == 1 ]]; then
	    FQIN=$CHIP
	    TAGIN=$TAGC
	    TREAT1=Bowtie/$TAGC-$GENROOT.sam
	    TREAT1R=rDNA/$TAGC-$GENROOT2.sam
	elif [[ $fq == 2 ]]; then
	    FQIN=$INPUT
	    TAGIN=$TAGI
	    CONTROL1=Bowtie/$TAGI-$GENROOT.sam
	    CONTROL1R=rDNA/$TAGI-$GENROOT2.sam
	elif [[ $fq == 3 ]]; then
	    FQIN=$CHCON
	    TAGIN=$TAGCCON
	    TREAT2=Bowtie/$TAGCCON-$GENROOT.sam
	    TREAT2R=rDNA/$TAGCCON-$GENROOT2.sam
	else
	    FQIN=$INCON
	    TAGIN=$TAGICON
	    CONTROL2=Bowtie/$TAGICON-$GENROOT.sam
	    CONTROL2R=rDNA/$TAGICON-$GENROOT2.sam
	fi

	if [[ $jidBs != "" ]]; then jidBs+=":"; fi
	jidBs+=$(sbatch --export FLDR=$FLDR,FASTQIN=$FQIN,TAG=$TAGIN,GENOME=$GEN \
	    ~/ChIPseq_Pipeline_v3/Bowtie_B3.sbatch | awk '{print $NF}')

	#for rDNA mapping
	if [[ " ${fqCount[@]} " == " 1 2 " || \
	    " ${fqCount[@]} " == " 1 2 3 4 " ]]; then
	    if [[ $jidRs != "" ]]; then jidRs+=":"; fi
	    if [[ $GEN == "SacCer3" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
		    --export FLDR=$FLDR,FASTQIN=$FQIN,TAG=$TAGIN,GENOME="SacCer3_rDNA" \
		    ~/ChIPseq_Pipeline_v3/Bowtie_B3.sbatch | awk '{print $NF}')
	    else
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
		    --export FLDR=$FLDR,FASTQIN=$FQIN,TAG=$TAGIN,GENOME="SK1_rDNA" \
		    ~/ChIPseq_Pipeline_v3/Bowtie_B3.sbatch | awk '{print $NF}')
	    fi
	fi
    done

    if [[ $jidBs != "" ]]; then
	if [[ $jids != "" ]]; then jids+=":"; fi
	jids+="$jidBs"
    fi

    if [[ $jidRs != "" ]]; then
	if [[ $jids != "" ]]; then jids+=":"; fi
	jids+="$jidRs"
    fi


# if GEN is empty, check for any fastq files and exit if TRUE
else
    if [[ $CHIP =~ .fastq|.fq || $INPUT  =~ .fastq|.fq || \
	$INCON =~ .fastq|.fq || $CHCON =~ .fastq|.fq ]]; then
	echo "Variable GEN is empty." >> $OUT
	echo "Please include genome type to analyze FASTQ files." >> $OUT
	run_closing $TMPID $jids $TAGC
	exit 1
    fi
fi

############### MACS2 ANALYSIS #################
if [[ $WIG == "" ]]; then WIG="TRUE"; fi 
if [[ $WIG =~ [Tt] || $PEAK != "" ]]; then
    echo "wiggle files or peak files requested" >> $OUT

    # some file name parsing for output filename information
    # only necessary if Bowtie alignment was not performed
    if [[ $TAGC == "" ]]; then
	parse_inCHIPname $OUT $CHIP $jids
    fi
    
    change_genroot $GENROOT
	
############### ChIP vs ChIP Normalization ###############
    # if working with ChIP vs ChIP normalization, run MACS2_FE.sbatch twice
    if [[ $INCON != "" && $CHCON != "" ]]; then
	if [[ $INPUT == "" || $CHIP == "" ]]; then
	    echo "For MACS2 analysis, this pipeline requires an input and a ChIP file" >> $OUT
	    run_closing $TMPID $jids $TAGC
	    exit 1
	else
	    # where roundA is the ChIP vs ChIP normalization and roundB is the input vs input normalization
	    if [[ $TREAT1 != "" ]]; then TREATA=$TREAT1; else TREATA=$CHIP; fi
	    if [[ $CONTROL1 != "" ]]; then TREATB=$CONTROL1; else TREATB=$INPUT; fi
	    if [[ $TREAT2 != "" ]]; then CONTROLA=$TREAT2; else CONTROLA=$CHCON; fi
	    if [[ $CONTROL2 != "" ]]; then CONTROLB=$CONTROL2; else CONTROLB=$INCON; fi
	    
	    if [[ $FLMKR == "" ]]; then
		    echo "For normalization against another ChIP sample," >> $OUT
		    echo "FLMKR must be defined" >> $OUT
		    run_closing $TMPID $jids $TAGC
		    exit 1
	    fi

	    # define output filenames here
	    # if no replicates expected
	    if [[ ! $TREATA =~ " " && ! $CONTROLA =~ " " && \
		! $TREATB =~ " " && ! $CONTROLB =~ " " ]]; then
		define_MACS2_filenames $TAGC $ROOTF $ROOTFR "TAGNORM" $FLMKR
	    # if replicates
	    else
		if [[ $REP != "" ]]; then
		    define_MACS2_filenames $REP $ROOTF $ROOTFR "TAGNORM_REPS" $FLMKR
		else
		    echo "Either CHIP or INPUT has multiple files." >> $OUT
		    echo "REP required" >> $OUT
		    run_closing $TMPID $jids $TAGC
		    exit 1
		fi    
	    fi

	    if [[ $jids != "" ]]; then jids+=":"; fi
	    if [[ $jidBs != "" ]]; then DEPB="--dependency=afterok:$jidBs"; fi
	    jids+=$(sbatch $DEPB \
		--export FLDR=$FLDR,WIG=$WIG,PEAK=$PEAK,TREAT="$TREATA",CONTROL="$CONTROLA",FOLDER=$FOLDER,M2FILE=$M2FILE1 \
		~/ChIPseq_Pipeline_v3/MACS2_FE_W3.sbatch | awk '{print $NF}')
	    jids+=":"
	    jids+=$(sbatch $DEPB \
		--export FLDR=$FLDR,WIG=$WIG,PEAK=$PEAK,TREAT="$TREATB",CONTROL="$CONTROLB",FOLDER=$FOLDER,M2FILE=$M2FILE2 \
		~/ChIPseq_Pipeline_v3/MACS2_FE_W3.sbatch | awk '{print $NF}')

	fi
    elif [[ $INCON != "" || $CHCON != "" ]]; then
	echo "For MACS2 analysis, this pipeline requires that" >> $OUT
	echo "if you define a control ChIP file, " >> $OUT
	echo "you must also include a control input file" >> $OUT
	run_closing $TMPID $jids $TAGC
	exit 1

################ ChIP vs Input Normalization ###############3
    # if doing a more traditional normalization
    else
	if [[ $INPUT == "" || $CHIP == "" ]]; then
	    echo "For MACS2 analysis, this pipeline requires an input and a ChIP file" >> $OUT
	    run_closing $TMPID $jids $TAGC
	    exit 1
	else
	    if [[ $TREAT1 != "" ]]; then TREAT=$TREAT1; else TREAT=$CHIP; fi
	    if [[ $CONTROL1 != "" ]]; then CONTROL=$CONTROL1; else CONTROL=$INPUT; fi

	    # define output filenames here
	    # if no replicates expected
	    if [[ ! $TREAT =~ " " && ! $CONTROL =~ " " ]]; then
		define_MACS2_filenames $TAGC $ROOTF $ROOTFR "NORMAL" $FLMKR
	    # if replicates
	    else
		if [[ $REP != "" ]]; then
		    define_MACS2_filenames $REP $ROOTF $ROOTFR "REPS" $FLMKR
		else
		    echo "Either CHIP or INPUT has multiple files." >> $OUT
		    echo "REP required" >> $OUT
		    run_closing $TMPID $jids $TAGC
		    exit 1
		fi    
	    fi
	    
	    if [[ $jids != "" ]]; then jids+=":"; fi
	    if [[ $jidBs != "" ]]; then DEPB="--dependency=afterok:$jidBs"; fi
	    jids+=$(sbatch $DEPB \
		--export FLDR=$FLDR,WIG=$WIG,PEAK=$PEAK,TREAT="$TREAT",CONTROL="$CONTROL",FOLDER=$FOLDER,M2FILE=$M2FILE \
		~/ChIPseq_Pipeline_v3/MACS2_FE_W3.sbatch | awk '{print $NF}')
	fi
    fi
fi

#################### rDNA Pipeline ################
# only runs if using the entire pipeline and genome is SacCer3
if [[ $GENROOT2 != "" ]]; then
    change_genroot $GENROOT2
    if [[ $INCON != "" ]]; then
	if [[ $REP != "" ]]; then
	    define_MACS2_filenames $REP $ROOTF $ROOTFR "TAGNORM_REPS" $FLMKR
	else
	    define_MACS2_filenames $TAGC $ROOTF $ROOTFR "TAGNORM" $FLMKR
	fi
	jids+=":"
	if [[ $jidRs != "" ]]; then DEPR="--dependency=afterok:$jidRs"; fi
	jids+=$(sbatch $DEPR \
	    --export FLDR=$FLDR,WIG="RDNA",PEAK=$PEAK,TREAT="$TREAT1R",CONTROL="$TREAT2R",FOLDER=$FOLDER,M2FILE=$M2FILE1 \
	    ~/ChIPseq_Pipeline_v3/MACS2_FE_W3.sbatch | awk '{print $NF}')
	jids+=":"
	jids+=$(sbatch $DEPR \
	    --export FLDR=$FLDR,WIG="RDNA",PEAK=$PEAK,TREAT="$CONTROL1R",CONTROL="$CONTROL2R",FOLDER=$FOLDER,M2FILE=$M2FILE2 \
	    ~/ChIPseq_Pipeline_v3/MACS2_FE_W3.sbatch | awk '{print $NF}')
    else
	if [[ $REP == "" ]]; then
	    define_MACS2_filenames $TAGC $ROOTF $ROOTFR "NORMAL" $FLMKR
	else
	    define_MACS2_filenames $REP $ROOTF $ROOTFR "REPS" $FLMKR
	fi
	jids+=":"
	if [[ $jidRs != "" ]]; then DEPR="--dependency=afterok:$jidRs"; fi
	jids+=$(sbatch $DEPR \
	    --export FLDR=$FLDR,WIG="RDNA",PEAK=$PEAK,TREAT="$TREAT1R",CONTROL="$CONTROL1R",FOLDER=$FOLDER,M2FILE=$M2FILE \
	    ~/ChIPseq_Pipeline_v3/MACS2_FE_W3.sbatch | awk '{print $NF}')
    fi
fi

################### CLOSING FUNCTION ################
run_closing $TMPID $jids $TAGC

echo "Job ids are: $jids" >> $OUT
echo "" >> $OUT





