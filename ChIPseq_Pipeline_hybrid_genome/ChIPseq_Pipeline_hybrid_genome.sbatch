#!/bin/bash
#
#SBATCH --nodes=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --time=2:00:00
#SBATCH --mem=26GB
#SBATCH --job-name=Hybrid_genome_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lv38@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out
#SBATCH --error=/scratch/%u/%x_%j.err

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#
# 1. If GENNAME is provided:
#       . Maps FASTQ data to indicated genome (allowing no mismatches)
#       . Converts output from SAM to sorted and indexed BAM format
# 2. Normalizes by library size using MACS2 SPMR
# 3. Removes noise using MACS2 fold enrichment compared to input

### Argument options:
# EXPID     Custom ID for output files.
# RUNDIR    Path to directory to run script and save output in.
# CHIP      Absolute path(s) to ChIP (treatment) sample file(s).
# INPUT     Absolute path(s) to input (control) sample file(s).
#           NOTE: CHIP and INPUT file formats can be either:
#                   - FASTQ         if GENNAME is provided; pipeline will start
#                                   with Bowtie mapping.
#                   - SAM or BAM    if GENNAME is not provided; pipeline will
#                                   skip Bowtie mapping.
# GENNAME   Basename of reference genome FASTA file preceded by absolute path to
#           directory containing the file. Must be provided if data files are
#           in FASTQ format; must not be provided if data files are alignments
#           (BAM or SAM format).
#           An existing Bowtie index with a basename ("ebwt") matching the
#           file's is used if found in the same directory; otherwise a new index
#           is built.

# WITH REPLICATES:
# The pipeline can be used with alignment maps (BAM or SAM files) of replicate
# samples. If starting from FASTQ files they must be run individually.
# Simply add file names in the appropriate variable separated by spaces:
# CHIP="path/to/chip1.bam path/to/chip2.bam"

### EXAMPLES: 
# NO REPLICATES
# sbatch --export EXPID="AH119spikein-060717_YueSK1_S288C_PM_SPMR",\
# RUNDIR="/scratch/lv38",\
# CHIP="/scratch/lv38/HLYHHAFXX_n01_ah119spikea-062817.fastq.gz",\
# INPUT="/scratch/lv38/HLYHHAFXX_n01_ah119spikeb-062817.fastq.gz",\
# GENNAME="/home/lv38/Library/S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue" \
# ~/Pipeline/ChIPseq_Pipeline_hybrid_genome.sbatch

# sbatch --export EXPID="AH119spikein-060717_YueSK1_S288C_PM_SPMR",\
# RUNDIR="/scratch/lv38",\
# CHIP="/scratch/lv38/ah119spikeb-062817_S288C_SK1_Yue_PM.sam",\
# INPUT="/scratch/lv38/ah119spikea-062817_S288C_SK1_Yue_PM.sam",\
# ~/Pipeline/ChIPseq_Pipeline_hybrid_genome.sbatch

# REPLICATES
# sbatch --export EXPID="Red1-WT-155-175-reps_S288C_SK1_Yue_PM_SPMR",\
# RUNDIR="/scratch/lv38",\
# CHIP="/scratch/lv38/ah119spikeb-062817_S288C_SK1_Yue_PM.sam \
# /scratch/lv38/ah119spiked-01012018_S288C_SK1_Yue_PM.sam",\
# INPUT="/scratch/lv38/ah119spikea-062817_S288C_SK1_Yue_PM.sam \
# /scratch/lv38/ah119spiked-01012018_S288C_SK1_Yue_PM.sam" \
# ~/Pipeline/ChIPseq_Pipeline_hybrid_genome.sbatch

#------------------------------------------------------------------------------#
#                                  Functions                                   #
#------------------------------------------------------------------------------#

function elapsed_time() # Time elapsed since a STARTTIME=$(date +%s")
 {
     local ENDTIME=$(date +%s)

     local TIME=$(($ENDTIME - $1))
     if [ $TIME -lt 60 ]
     then
         echo "$TIME second(s)"
     elif [ $TIME -ge 60 ]  && [ $TIME -lt 3600 ]
     then
         echo "$(($TIME / 60)) minute(s)"
     else
         echo "$(($TIME / 60 / 60)) hour(s)"
     fi
 }

function check_arg() # Check if a variable is available (and not empty)
 {
     if [ -z "$1" ]
     then
         echo "Please provide values for all required arguments"
         exit 1
     fi
 }

function sam_to_sorted_and_indexed_bam() # Convert to BAM file, sort, and index
 {
     local BASE=${1%.sam}
     local BAM=${BASE}.bam
     local S_BAM=${BASE}_sorted.bam

     samtools view -bS $1 > $BAM
     samtools sort -o $S_BAM $BAM
     samtools index $S_BAM

     # Clean up
     rm $1
     rm $BAM

     echo "$S_BAM"
}

#------------------------------------------------------------------------------#

# Check arguments
check_arg $EXPID
check_arg $RUNDIR 
check_arg $CHIP
check_arg $INPUT

# Abort if: GENNAME is provided with alignment files
#           GENNAME is not provided with FASTQ files
if [ ! -z "$GENNAME" ]
then
    # Make sure the files are not alignments (GENNAME mistakenly provided)
    case "${CHIP: -4}" in .sam|.SAM|.bam|.BAM)
        echo "GENNAME provided: expected FASTQ files"; exit 2 ;;
    esac
    case "${INPUT: -4}" in .sam|.SAM|.bam|.BAM)
        echo "GENNAME provided: expected FASTQ files"; exit 2 ;;
    esac
else
    # Make sure the files are alignments (GENNAME mistakenly not provided)
    case "${CHIP: -4}" in .sam|.SAM|.bam|.BAM)
        true ;;
    *)
        echo "GENNAME not provided: expected alignment files"; exit 2 ;;
    esac
    case "${INPUT: -4}" in .sam|.SAM|.bam|.BAM)
        true ;;
    *)
        echo "GENNAME not provided: expected alignment files"; exit 2 ;;
    esac
fi

# Abort if output directory already exists
if [ -d "$EXPID" ]
then
    echo "Output directory already exists"
    exit 2
fi

# Check input files / dirs
[ -f $CHIP ] || { echo "Could not find file $CHIP"; exit 2; }
[ -f $INPUT ] || { echo "Could not find file $INPUT"; exit 2; }
[ -d $RUNDIR ] || { echo "Could not find directory $RUNDIR"; exit 2; }

# Was GENNAME provided?
if [ ! -z "$GENNAME" ]
then
    GENDIR=$(dirname "$GENNAME")
    GENNAME=$(basename "$GENNAME")
    [ -d $GENDIR ] || { echo "Could not find directory $GENDIR"; exit 2; }

    # Search for reference genome file; exit if not found
    FA=$(find $GENDIR -iname "${GENNAME}.fa*")
    if [ -z "$FA" ]
    then
        echo "Could not find reference genome FASTA file"
        echo "(Searching for file named ${GENNAME}.fa*)"
        exit 2
    fi

    # Search for Bowtie index (a file named as "fasta_base_name.1.ebwt")
    # Build index if not found
    # The following code will return the full basename
    # (the provided $GENNAME might not include it in full)
    IX=$(basename $FA)                               # Drop path to file
    IX=${IX%.*}                                      # Drop extension
    CHECKIX=$(find $GENDIR -iname "${IX}.1.ebwt")    # Search file
    #IX=$(basename $IX | cut -d '.' -f 1)

    if [ -z "$CHECKIX" ]
    then
        echo ">>>>> Building Bowtie index..."
        module purge
        module load bowtie/gnu/1.2.0
        # Build index
        cd $GENDIR
        bowtie-build -f $FA $IX
    else
        echo ">>>>> Found Bowtie index:"
        echo "      $GENDIR/$IX"
    fi
fi

#------------------------------------------------------------------------------#
#                                                                              #
#                                Run pipeline                                  #
#                                                                              #
#------------------------------------------------------------------------------#

STARTTIME=$(date +%s)
echo ">>>>> Started pipeline: $EXPID"
echo "      `date`"

cd $RUNDIR
mkdir ${EXPID}/
cd ${EXPID}/

#------------------------------------------------------------------------------#
#                           If GENNAME is provided:                            #
#                 Align reads to reference genome with Bowtie                  #
#               Convert output SAM files to BAM, sort and index                #
#------------------------------------------------------------------------------#

if [ ! -z "$GENNAME" ]
then
    # Bowtie 1.2.0 does not accept gzipped files
    # Unzip them if gzipped (set empty variable to later know to rezip them)
    if [ ${CHIP: -3} == ".gz" ]
    then
        echo ">>>>> Unzip $CHIP"
        gzip -d $CHIP
        CHIP=${CHIP%.gz}
        ZIP_CHIP=
    fi

    if [ ${INPUT: -3} == ".gz" ]
    then
        echo ">>>>> Unzip $INPUT"
        gzip -d $INPUT
        INPUT=${INPUT%.gz}
        ZIP_INPUT=
    fi

    echo ">>>>> Align reads with Bowtie:"
    [ -d $RUNDIR/Bowtie_hybrids ] || mkdir $RUNDIR/Bowtie_hybrids
    [ -d $RUNDIR/Unaligned_hybrids ] || mkdir $RUNDIR/Unaligned_hybrids
    
    module purge
    module load bowtie/gnu/1.2.0
    
    # Get input file names (without path to file or extension)
    # Build output file names
    CHIP_NAME=$(basename ${CHIP%.f*})
    INPUT_NAME=$(basename ${INPUT%.f*})
    CHIP_OUT=$RUNDIR/Bowtie_hybrids/${CHIP_NAME}_${IX}-PM.sam
    INPUT_OUT=$RUNDIR/Bowtie_hybrids/${INPUT_NAME}_${IX}-PM.sam

    echo "      $CHIP"
    bowtie -q -m 1 -v 0 -p 8 -S \
    --un $RUNDIR/Unaligned_hybrids/${CHIP_NAME}_${IX}_PM-Unaligned.fastq \
    --max $RUNDIR/Unaligned_hybrids/${CHIP_NAME}_${IX}_PM-Max.fastq \
    $GENDIR/$IX $CHIP $CHIP_OUT
    
    echo "      $INPUT"
    bowtie -q -m 1 -v 0 -p 8 -S \
    --un $RUNDIR/Unaligned_hybrids/${INPUT_NAME}_${IX}_PM-Unaligned.fastq \
    --max $RUNDIR/Unaligned_hybrids/${INPUT_NAME}_${IX}_PM-Max.fastq \
    $GENDIR/$IX $INPUT $INPUT_OUT

    # Were the files originally gzipped?
    [ ! -z ${ZIP_CHIP+x} ] && echo ">>>>> gzip $CHIP" && gzip $CHIP
    [ ! -z ${ZIP_INPUT+x} ] && echo ">>>>> gzip $INPUT" && gzip $INPUT

    # convert to BAM, sort and index
    # In the process, reset file variables to map files (to use in next step)
    echo ">>>>> Convert SAM to sorted and indexed BAM:"
    module purge
    module load samtools/intel/1.3.1

    CHIP=$(sam_to_sorted_and_indexed_bam $CHIP_OUT)
    echo "      $CHIP"
    INPUT=$(sam_to_sorted_and_indexed_bam $INPUT_OUT)
    echo "      $INPUT"
    
    # Reset file variables to map files (to use in the next step)
    #CHIP=$CHIP_OUT
    #INPUT=$INPUT_OUT
fi

#------------------------------------------------------------------------------#
#                          MACS2 callpeak with SPMR                            #
#------------------------------------------------------------------------------#

echo ">>>>> Run macs2 callpeak..."
module purge
module load macs2/intel/2.1.1

macs2 callpeak -t $CHIP -c $INPUT \
    --keep-dup="auto" -B --nomodel \
    --extsize 200 --SPMR -g 2.4e7 \
    -n $EXPID

#------------------------------------------------------------------------------#
#                            MACS2 bdgcmp with FE                              #
#------------------------------------------------------------------------------#

echo ">>>>> Run macs2 bdgcmp with FE..."
macs2 bdgcmp \
    -t ${EXPID}_treat_pileup.bdg \
    -c ${EXPID}_control_lambda.bdg \
    -o ${EXPID}_FE.bdg -m FE

#------------------------------------------------------------------------------#
ELAPSEDTIME=$(elapsed_time $STARTTIME)
echo "-----"
echo "-----"
echo "Completed pipeline in $ELAPSEDTIME"
