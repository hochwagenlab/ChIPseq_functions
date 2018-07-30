#! /bin/bash
##SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --mem=10GB
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
#   BDG       -- Whether or not to create bedgraphs. Options are: "TRUE" or "FALSE" [ default: "T" ]            #
#                Not needed when WIG is TRUE.                                                                   #
#   WIG       -- Whether or not to create wiggle files. Options are: "TRUE" or "FALSE" [ default: "F" ]         #
#   FLMKR     -- filemaker IDs from Datasets of ChIP samples in "-"-delimited and quoted list                   #
#                Should not be strain IDS. Will be incorporated into all MACS2-related files.                   #
#                REQUIRED for replicates and ChIP-vs-ChIP normalization.                                        #
#   PE        -- Whether or not the samples are paired-end and should be treated as such [ default: "F" ]       #
#                                                                                                               #
# CHIP vs CHIP NORMALIZATION                                                                                    #
## For mock controls, no tag controls, or histones                                                              #
#   CHCON     -- control ChIP file. Must be matched with a CHIP file, must be SAM files.                        #
#                                                                                                               #
# REPLICATES                                                                                                    #
#   REP       -- tagname for MACS2 outputs of replicates                                                        #
#   CHIP(etc) -- all included files must be SAM/BAM files [currently no error checking].                        #
#                CHIP,INPUT, and CHCON will all accept replicates as a " "-delimited list                       #
#                                                                                                               #
# BOWTIE ALIGNMENT                                                                                              #
## Available for CHIP,INPUT,INCON, and CHCON [Fastq files can, but do not have to be zipped.]                   #
## Cannot be combined with replicates                                                                           #
#   TAGI      -- Tagname for INPUT fastq                                                                        #
#   TAGC      -- Tagname for CHIP fastq                                                                         #
#   GEN       -- Genome to align to. Current options: SK1K, SacCer3, SK1Yue-PM, SK1Yue-2mis                     #
#   CLIP3P    -- Number of bases to clip from the 3' end of reads before aligning [ default: 0 ]                #
#                                                                                                               #
# SOME EXAMPLES                                                                                                 #
## Example 1: Start with input file A.fastq and ChIP file B.fastq, map to SacCer3, and create wiggle files      #
#             and both narrow and broad peak files (will also create wiggle plot of the rDNA)                   #
#      sbatch --export INPUT=A.fastq,CHIP=B.fastq,GEN="SacCer3",TAGI="A",TAGC="B",PEAK="BOTH" \                 #
#             ~/ChIPseq_Pipeline_v4/ChIPseq-pipeline_v4.sbatch                                                  #
#                                                                                                               #
## Example 2: Get broad peak files for replicate input files A_1.sam and A_2.sam and                            #
#             ChIP files B_1.sam and B_2.sam                                                                    #
#      sbatch --export INPUT="A_1.sam A_2.sam",CHIP="B_1.sam B_2.sam",REP="B",FLMKR="1-3",PEAK="BROAD",BDG="F" \#
#             ~/ChIPseq_Pipeline_v4/ChIPseq-pipeline_v4.sbatch                                                  #
#                                                                                                               #
## Example 3: Get normalized bedgraph files (and narrow peaks) with ChIP file A1.sam and mock ChIP file B1.sam. #
#      sbatch --export CHIP=A1.sam,CHCON=B1.sam,FLMKR="1-2" \                                                   #
#             ~/ChIPseq_Pipeline_v4/ChIPseq-pipeline_v4.sbatch                                                  #
#                                                                                                               #
## Example 4: Start with input file A.fastq and ChIP file B.fastq, map to SK1Yue with no mismatches, clipping   #
#             50bp off the 3' end, and create bedgraph files andnarrow file.                                    #
#      sbatch --export INPUT=A.fastq,CHIP=B.fastq,GEN="SK1Yue-PM",TAGI="A",TAGC="B",CLIP3P=50 \                 #
#             ~/ChIPseq_Pipeline_v4/ChIPseq-pipeline_v4.sbatch                                                  #
#                                                                                                               #
#################################################################################################################               

# check for shell functions
if [ -f ~/ChIPseq_Pipeline_v4/ChIPseq_functions_prince.sh ]; then
    source ~/ChIPseq_Pipeline_v4/ChIPseq_functions_prince.sh
else
    echo "Cannot find MACS2_functions_prince.sh"
    exit 3
fi

# TMPID should be the root of the ChIP file (without any periods)
# naming of out files, if no other sbatch scripts run
TMP=${CHIP##*/}
TMPID1=${TMP%%.*}
if [[ $TMPID =~ "*" ]]; then
    IFS="*" read -ra TMPID <<< $TMPID1
else
    TMPID=$TMPID1
fi    

OUT=/scratch/$USER/MACS2_pipeline_$TMPID.txt
touch $OUT

# if other sbatch scripts run, use this notation for out file
# new system for naming file output file
if [[ $REP != "" ]]; then
    TAGOUT=$REP
elif [[ $TAGC != "" ]]; then
    TAGOUT=$TAGC
else
    TAGOUT=$TMPID
fi

# change to folder
if [[ $FLDR == "" ]]; then
    cd $PWD
else
    cd /scratch/${USER}
    cd $FLDR
fi
echo "Working directory is specified as: $PWD" >> $OUT

# preallocating for jobid lists
jids="" # all jobs
jidBs="" # Bowtie jobs
jidRs="" # rDNA Bowtie jobs

FLMKR=$(echo $FLMKR | tr ':' '-')

############## BOWTIE ALIGNMENT ################
############### preparation ####################
if [[ $GEN != "" ]]; then
    echo "Bowtie mapping to $GEN requested." >> $OUT

    if [[ $CLIP3P == "" ]]; then
	CLIP3P=0
    fi
    
    # ensure that there are no replicates
    if [[ $CHIP =~ " " || $INPUT =~ " " || \
	$CHCON =~ " " || $INCON =~ " " ]]; then
	echo "Error: spaces included for variable CHIP." >> $OUT
	echo "Cannot have spaces and run Bowtie alignment for the same variable." >> $OUT
	echo "If file has spaces, use rm_space_filename.py" >> $OUT
	run_closing $TMPID $jids $TAGOUT
	exit 1
    fi

    # fqCount is an array which denotes which of the four variables 
    # have fastq files (CHIP, INPUT, INCON, CHCON)
    # make sure that replicates are not involved (check files for spaces)
    fqCount=()

    check_fastq $TMPID $TAGOUT $CHIP $TAGC $jids
    if [[ $FQTF == "TRUE" ]]; then fqCount+=(1); fi

    check_fastq $TMPID $TAGOUT $INPUT $TAGI $jids
    if [[ $FQTF == "TRUE" ]]; then fqCount+=(2); fi

    # adjusting output names to fit with genome requested
    assign_genroot $TMPID $TAGOUT $GEN $PE $jids


###### dealing with TAGNAME aka mapped read length #######
    
# first check if file can be found (either zipped or not) and unzip if necessary
RE1="(.+)\.gz"
if [[ $CHIP =~ $RE1 ]]; then
   FASTQ1=${BASH_REMATCH[1]}
else
    FASTQ1=$CHIP
fi

if [ -f $FASTQ1.gz ]; then
    if [ ! -f $FASTQ1 ]; then
	gzip -d $FASTQ1.gz
    fi
elif [ ! -f $FASTQ1 ]; then
     echo "Cannot find first file $FASTQ1. Quitting."
     exit 1
fi    

# to adjust for different sized reads to clip t0 50bp for alignment
# aka to define -3
# defines number of bases to be clipped off
if [[ $CLIP3P == "" ]]; then
    Prime3=0
else
    Prime3=$CLIP3P
fi

# calculate length of the read
READLEN=$(head -n 2 $FASTQ1 | tail -n 1 | wc -m)
if [[ $READLEN == 51 || $READLEN == 76 ]]; then
    MAPLENGTH=$(($READLEN-$Prime3-1))
else
    MAPLENGTH=$(($READLEN-$Prime3-2))
fi

# new variable for adjusting read name
if [[ $MAPLENGTH != 50 ]]; then
    MAPLEN="-${MAPLENGTH}bp"
else
     MAPLEN=""
fi
    
################# setting up bowtie sbatch ###############
    for fq in ${fqCount[@]}
    do
	if [[ $fq == 1 ]]; then
	    FQIN=$CHIP
	    TAGIN=$TAGC
	    TAGMACS=$TAGC$MAPLEN
	    TREAT1=Bowtie/$TAGC$MAPLEN-$GENROOT.sam
	    TREAT1R=rDNA/$TAGC$MAPLEN-$GENROOT2.sam
	elif [[ $fq == 2 ]]; then
	    FQIN=$INPUT
	    TAGIN=$TAGI
	    CONTROL1=Bowtie/$TAGI$MAPLEN-$GENROOT.sam
	    CONTROL1R=rDNA/$TAGI$MAPLEN-$GENROOT2.sam
	fi

	if [[ $PE =~ [Tt] ]]; then
	    if [[ $jidBs != "" ]]; then jidBs+=":"; fi
	    jidBs+=$(sbatch --export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME=$GEN,PE="TRUE",CLIP3P=$CLIP3P \
			    ~/ChIPseq_Pipeline_v4/Bowtie_PE_B4.sbatch | awk '{print $NF}')

	    #for rDNA mapping
	    if [[ $jidRs != "" ]]; then jidRs+=":"; fi
	    if [[ $GEN == "SacCer3" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
				--export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME="SacCer3_rDNA",PE="TRUE",CLIP3P=$CLIP3P \
				~/ChIPseq_Pipeline_v4/Bowtie_PE_B4.sbatch | awk '{print $NF}')
	    elif [[ $GEN == "SK1K" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
				--export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME="SK1K_rDNA",PE="TRUE",CLIP3P=$CLIP3P \
				~/ChIPseq_Pipeline_v4/Bowtie_PE_B4.sbatch | awk '{print $NF}')
	    elif [[ $GEN =~ "SK1Yue" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
				--export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME="SK1Yue_rDNA",PE="TRUE",CLIP3P=$CLIP3P \
				~/ChIPseq_Pipeline_v4/Bowtie_PE_B4.sbatch | awk '{print $NF}')		
	    fi
	else
	    if [[ $jidBs != "" ]]; then jidBs+=":"; fi
	    jidBs+=$(sbatch --export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME=$GEN,CLIP3P=$CLIP3P \
			    ~/ChIPseq_Pipeline_v4/Bowtie_B3.sbatch | awk '{print $NF}')

	    #for rDNA mapping
	    if [[ $jidRs != "" ]]; then jidRs+=":"; fi
	    if [[ $GEN == "SacCer3" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
				--export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME="SacCer3_rDNA",CLIP3P=$CLIP3P \
				~/ChIPseq_Pipeline_v4/Bowtie_B3.sbatch | awk '{print $NF}')
	    elif [[ $GEN == "SK1K" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
				--export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME="SK1K_rDNA",CLIP3P=$CLIP3P \
				~/ChIPseq_Pipeline_v4/Bowtie_B3.sbatch | awk '{print $NF}')
	    elif [[ $GEN =~ "SK1Yue" ]]; then
		jidRs+=$(sbatch --dependency=afterany:$jidBs \
				--export FLDR=$FLDR,FASTQIN=$FQIN,TAGNAME=$TAGIN,GENOME="SK1Yue_rDNA",CLIP3P=$CLIP3P \
				~/ChIPseq_Pipeline_v4/Bowtie_B3.sbatch | awk '{print $NF}')		
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
	run_closing $TMPID $TAGOUT $jids
	exit 1
    fi
fi

############### MACS2 ANALYSIS #################
if [[ $BDG == "" ]]; then BDG="TRUE"; fi
if [[ $WIG == "" ]]; then WIG="FALSE"; fi
if [[ $BDG =~ [Tt] || $PEAK != "" || $WIG =~ [Tt] ]]; then
    echo "bedgraphs, wiggle files, or peak files requested" >> $OUT

    # some file name parsing for output filename information
    # only necessary if Bowtie alignment was not performed
    if [[ $TAGC == "" ]]; then
	parse_inCHIPname $TMPID $TAGOUT $CHIP $jids
    fi
    
    change_genroot $GENROOT
	
############### ChIP vs ChIP Normalization ###############
    if [[ $CHCON != "" ]]; then
	if [[ $CHIP == "" ]]; then
	    echo "For MACS2 analysis, this pipeline requires an control and a ChIP file" >> $OUT
	    run_closing $TMPID $TAGOUT $jids 
	    exit 1
	else
	    TREATA=$CHIP
	    CONTROLA=$CHCON
	    
	    if [[ $FLMKR == "" ]]; then
		    echo "For normalization against another ChIP sample," >> $OUT
		    echo "FLMKR must be defined" >> $OUT
		    run_closing $TMPID $TAGOUT $jids
		    exit 1
	    fi

	    # define output filenames here
	    # if no replicates expected
	    if [[ ! $TREATA =~ " " && ! $CONTROLA =~ " " ]]; then
		    define_MACS2_filenames $TAGC $ROOTF $ROOTFR "TAGNORM" $FLMKR
	    # if replicates
	    else
		if [[ $REP != "" ]]; then
		    define_MACS2_filenames $REP $ROOTF $ROOTFR "TAGNORM_REPS" $FLMKR
		else
		    echo "Either CHIP or CONTROL has multiple files." >> $OUT
		    echo "REP required" >> $OUT
		    run_closing $TMPID $jids $TAGOUT
		    exit 1
		fi    
	    fi

	    if [[ $jids != "" ]]; then jids+=":"; fi
	    jids+=$(sbatch $DEPB \
		--export FLDR=$FLDR,BDG=$BDG,WIG=$WIG,PEAK=$PEAK,TREAT="$TREATA",CONTROL="$CONTROLA",FOLDER=$FOLDER,M2FILE=$M2FILE1,PE=$PE \
		~/ChIPseq_Pipeline_v4/MACS2_FE_W4.sbatch | awk '{print $NF}')
	fi
    elif [[ $CHCON != "" ]]; then
	echo "For MACS2 analysis, this pipeline requires that" >> $OUT
	echo "if you define a control ChIP file, " >> $OUT
	echo "you must also include a control input file" >> $OUT
	run_closing $TMPID $jids $TAGOUT $jids

################ ChIP vs Input Normalization ###############3
    # if doing a more traditional normalization
    else
	if [[ $INPUT == "" || $CHIP == "" ]]; then
	    echo "For MACS2 analysis, this pipeline requires an input and a ChIP file" >> $OUT
	    run_closing $TMPID $jids $TAGOUT 
	    exit 1
	else
	    if [[ $TREAT1 != "" ]]; then TREAT=$TREAT1; else TREAT=$CHIP; fi
	    if [[ $CONTROL1 != "" ]]; then CONTROL=$CONTROL1; else CONTROL=$INPUT; fi

	    # define output filenames here
	    # if no replicates expected
	    if [[ ! $TREAT =~ " " && ! $CONTROL =~ " " ]]; then
		if [[ $CHIP =~ ".sam" ]]; then
		    define_MACS2_filenames $TAGC $ROOTF $ROOTFR "NORMAL" $FLMKR
		    else
		    define_MACS2_filenames $TAGMACS $ROOTF $ROOTFR "NORMAL" $FLMKR
		fi
	    # if replicates
	    else
		if [[ $REP != "" ]]; then
		    define_MACS2_filenames $REP $ROOTF $ROOTFR "REPS" $FLMKR
		else
		    echo "Either CHIP or INPUT has multiple files." >> $OUT
		    echo "REP required" >> $OUT
		    run_closing $TMPID $jids $TAGOUT
		    exit 1
		fi    
	    fi
	    
	    if [[ $jids != "" ]]; then jids+=":"; fi
	    if [[ $jidBs != "" ]]; then DEPB="--dependency=afterok:$jidBs --kill-on-invalid-dep=yes"; fi
	    jids+=$(sbatch $DEPB \
		--export FLDR=$FLDR,BDG=$BDG,WIG=$WIG,PEAK=$PEAK,TREAT="$TREAT",CONTROL="$CONTROL",FOLDER=$FOLDER,M2FILE=$M2FILE,PE=$PE \
		~/ChIPseq_Pipeline_v4/MACS2_FE_W4.sbatch | awk '{print $NF}')
	fi
    fi

#################### rDNA Pipeline ################
# only runs if using the entire pipeline [Bowtie included]
    if [[ $GENROOT2 != "" ]]; then
	change_genroot $GENROOT2
	if [[ $REP == "" ]]; then
	    define_MACS2_filenames $TAGC $ROOTF $ROOTFR "NORMAL" $FLMKR
	else
	    define_MACS2_filenames $REP $ROOTF $ROOTFR "REPS" $FLMKR
	fi
	jids+=":"
	if [[ $jidRs != "" ]]; then DEPR="--dependency=afterok:$jidRs --kill-on-invalid-dep=yes"; fi
	jids+=$(sbatch $DEPR \
	    --export FLDR=$FLDR,WIG="RDNA",PEAK=$PEAK,TREAT="$TREAT1R",CONTROL="$CONTROL1R",FOLDER=$FOLDER,M2FILE=$M2FILE,PE=$PE \
	    ~/ChIPseq_Pipeline_v4/MACS2_FE_W4.sbatch | awk '{print $NF}')
    fi
fi

################### CLOSING FUNCTION ################
run_closing $TMPID $jids $TAGOUT

echo "Job ids are: $jids" >> $OUT
echo "" >> $OUT





