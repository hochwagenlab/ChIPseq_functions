#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --time=0:30:00
#SBATCH --job-name=Bowtie_PE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tem298@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out

module purge
module load bowtie/gnu/1.2.0
module load samtools/intel/1.3.1
module load perl/intel/5.24.0

echo "#############################################################"
echo ""
echo "Running: $SLURM_JOB_NAME"
echo "JobID is: $SLURM_JOBID "

####################################
## PURPOSE ##
# To use Bowtie and a few associated tools to create SAM files
# Designed specifically for Paired-End sequencing
# Can be used to map to: SacCer3, SK1, or SacCer3 rDNA (single repeat)
# SacCer3: 2 mismatches allowed, read must only map to one position
# SK1: only perfect matches, read must only map to one position
# rDNA: maps to a single rDNA repeat from SK1 and adjusted to align with
#         either SK1 or SacCer3 reference genome
#       2 mismatches allowed, keep all alignments [ file has 50bp repeating sequence ]

## VARIABLES ##
#
# FLDR:         [optional] working directory
#               default is working directory when submitting
# FASTQIN:      the fastq file to be aligned to a reference genome (zipped or unzipped)
#               only one file from the pair is needed
# TAG:          the tag name for the output sam file
# GENOME:       the reference genome to align to
#               Options include: SacCer3, SK1K, SK1Yue, SacCer3_rDNA, SK1K_rDNA, or SK1Yue_rDNA

####################################
# adjusting file name information is as easy as adjusting $ROOT within each section

# Define -X parameter of Bowtie for maximum distance between paired reads after alignment
X=1000

# Go to the currect starting folder
cd /scratch/$USER
if [[ $FLDR != "" ]]; then
    cd $FLDR
else
    cd $SLURM_SUBMIT_DIR
fi
DIR=$PWD

echo "working in folder: $DIR"
echo "running Bowtie on $FASTQIN with tagname $TAG to $GENOME"
date

# first check if file can be found (either zipped or not) and unzip if necessary
RE1="(.+)\.gz"
if [[ $FASTQIN =~ $RE1 ]]; then
   FASTQ1=${BASH_REMATCH[1]}
else
    FASTQ1=$FASTQIN
fi

if [ -f $FASTQ1.gz ]; then
    if [ ! -f $FASTQ1 ]; then
	gzip -d $FASTQ1.gz
    fi
elif [ ! -f $FASTQ1 ]; then
     echo "Cannot find first file $FASTQ1. Quitting."
     exit 1
fi

# second ensure that there are two matching paired end sequencing files
RE2="(.+)(R1|R2|n01|n02)(.+)"
if [[ $FASTQ1 =~ $RE2 ]]; then
    P1=${BASH_REMATCH[1]}
    P2=${BASH_REMATCH[2]}
    P3=${BASH_REMATCH[3]}
fi

if [[ $P2 =~ 1 ]]; then
    P2="${P2%?}2"
else
    P2="${P2%?}1"
fi

FASTQ2=$P1$P2$P3

# repeat searching for file
if [ -f $FASTQ2.gz ]; then
    if [ ! -f $FASTQ2 ]; then
	gzip -d $FASTQ2.gz
    fi
elif [ ! -f $FASTQ2 ]; then
     gzip $FASTQ2
     echo "Cannot find second file $FASTQ2. Quitting."
     exit 1
fi

test -d Bowtie || mkdir Bowtie
test -d Unaligned || mkdir Unaligned
test -d ~/Library/Bowtie1_2 || mkdir ~/Library/Bowtie1_2

########## SK1_KEENEY ###########
# outer if statement to check genome requested
if [[ $GENOME == "SK1K" ]]; then
    ROOT="SK1K-PM-PE_B4"
    # check if sam file exists in either of the expected spots
   if [ -f $TAG-$ROOT.sam ]; then
      echo "Sam file $TAG-$ROOT.sam already exists. Continuing."
   elif [ -f Bowtie/$TAG-$ROOT.sam ]; then
      echo "Sam file Bowtie/$TAG-$ROOT.sam already exists. Continuing."
   else
      # check if database files exist
      if [ ! -f ~/Library/Bowtie1_2/sk1_MvO_V1.1.ebwt ]; then
	  if [ ! -f ~/Library/sk1_MvO_V1.fasta ]; then
	      echo "Put sk1_MvO_V1.fasta in ~/Library/ to align files to SK1K"
	      exit 1
	  else
	      cd ~/Library/Bowtie1_2
	      bowtie-build -f ../sk1_MvO_V1.fasta sk1_MvO_V1
	      cd $DIR
	  fi
      fi
      bowtie -p 8 -q -m 1 -v 0 -X $X -S --un Unaligned/$TAG-${ROOT}_Unaligned.fastq \
	  --max Unaligned/$TAG-${ROOT}_Max.fastq \
	  ~/Library/Bowtie1_2/sk1_MvO_V1 -1 $FASTQ1 -2 $FASTQ2 \
	  Bowtie/$TAG-$ROOT.sam --verbose
      gzip Unaligned/$TAG*.fastq
   fi

elif [[ $GENOME == "SK1Yue-PM" ]]; then
    ROOT="SK1Yue-PM-PE_B4"
    # check if sam file exists in either of the expected spots
   if [ -f $TAG-$ROOT.sam ]; then
      echo "Sam file $TAG-$ROOT.sam already exists. Continuing."
   elif [ -f Bowtie/$TAG-$ROOT.sam ]; then
      echo "Sam file Bowtie/$TAG-$ROOT.sam already exists. Continuing."
   else
      # check if database files exist
      if [ ! -f ~/Library/Bowtie1_2/SK1Yue.1.ebwt ]; then
	  if [ ! -f ~/Library/Yue.SK1.genome.fa ]; then
	      echo "Put Yue.SK1.genome.fa in ~/Library/ to align files to SK1Yue"
	      exit 1
	  else
	      cd ~/Library/Bowtie1_2
	      bowtie-build -f ../Yue.SK1.genome.fa SK1Yue
	      cd $DIR
	  fi
      fi
      bowtie -p 8 -q -m 1 -v 0 -S -X $X --un Unaligned/$TAG-${ROOT}_Unaligned.fastq \
	  --max Unaligned/$TAG-${ROOT}_Max.fastq \
	  ~/Library/Bowtie1_2/SK1Yue -1 $FASTQ1 -2 $FASTQ2 \
	  Bowtie/$TAG-$ROOT.sam --verbose
      gzip Unaligned/$TAG*.fastq
   fi

elif [[ $GENOME == "SK1Yue-2mis" ]]; then
   ROOT="SK1Yue-2mis-PE_B4"
    # check if sam file exists in either of the expected spots
   if [ -f $TAG-$ROOT.sam ]; then
      echo "Sam file $TAG-$ROOT.sam already exists. Continuing."
   elif [ -f Bowtie/$TAG-$ROOT.sam ]; then
      echo "Sam file Bowtie/$TAG-$ROOT.sam already exists. Continuing."
   else
      # check if database files exist
      if [ ! -f ~/Library/Bowtie1_2/SK1Yue.1.ebwt ]; then
	  if [ ! -f ~/Library/Yue.SK1.genome.fa ]; then
	      echo "Put Yue.SK1.genome.fa in ~/Library/ to align files to SK1Yue"
	      exit 1
	  else
	      cd ~/Library/Bowtie1_2
	      bowtie-build -f ../Yue.SK1.genome.fa SK1Yue
	      cd $DIR
	  fi
      fi
      bowtie -p 8 -q -m 1 -S -X $X --un Unaligned/$TAG-${ROOT}_Unaligned.fastq \
	  --max Unaligned/$TAG-${ROOT}_Max.fastq \
	  ~/Library/Bowtie1_2/SK1Yue -1 $FASTQ1 -2 $FASTQ2 \
	  Bowtie/$TAG-$ROOT.sam --verbose
      gzip Unaligned/$TAG*.fastq
   fi

########### SacCer3 ##########
elif [[ $GENOME == "SacCer3" ]]; then
    ROOT="SacCer3-2mis-PE_B4"
    # check if sam file exists in either of the expected spots              
   if [ -f $TAG-$ROOT.sam ]; then
      echo "Sam file $TAG-$ROOT.sam already exists. Continuing."
   elif [ -f Bowtie/$TAG-$ROOT.sam ]; then
      echo "Sam file Bowtie/$TAG-$ROOT.sam already exists. Continuing."
   else
      # check if database files exist 
      if [ ! -f ~/Library/Bowtie1_2/SacCer3.1.ebwt ]; then
	  if [ ! -f ~/Library/S288C_R64.fasta ]; then
	      echo "Put S288C_R64.fasta into ~/Library/ to align files to SacCer3"
	      exit 1
	  else
	      cd ~/Library/Bowtie1_2
	      bowtie-build -f ../S288C_R64.fasta SacCer3
	      cd $DIR
	  fi
      fi
      bowtie -p 8 -q -m 1 -X $X -S --un Unaligned/$TAG-${ROOT}_Unaligned.fastq \
	  --max Unaligned/$TAG-${ROOT}_Max.fastq \
	  ~/Library/Bowtie1_2/SacCer3 -1 $FASTQ1 -2 $FASTQ2 \
	  Bowtie/$TAG-$ROOT.sam --verbose
      gzip Unaligned/$TAG*.fastq
   fi

######### SacCer3 rDNA #########
elif [[ $GENOME =~ "rDNA" ]]; then
    if [[ $GENOME =~ "SacCer3" ]]; then
	ROOT="SacCer3-rDNA-PE_B4"
    elif [[ $GENOME =~ "SK1K" ]]; then
	ROOT="SK1K-rDNA-PE_B4"
    elif [[ $GENOME =~ "SK1Yue" ]]; then
	ROOT="SK1Yue-rDNA-PE_B4"
    fi
   test -d rDNA || mkdir rDNA

   # check that file doesn't already exist
   if [ -f $TAG-$ROOT.sam ]; then
      echo "Sam file $TAG-$ROOT.sam already exists. Continuing."
   elif [ -f Bowtie/$TAG-$ROOT.sam ]; then
      echo "Sam file Bowtie/$TAG-$ROOT.sam already exists. Continuing."
   elif [ -f rDNA/$TAG-$ROOT.sam ]; then
      echo "Sam file rDNA/$TAG-$ROOT.sam already exists. Continuing."
   else
      # check if database files exist 
      if [ ! -f ~/Library/Bowtie1_2/rDNA_SK1.1.ebwt ]; then
	  if [ ! -f ~/Library/SK1_rDNA_50bp_overhang.fasta ]; then
	      echo "Put SK1_rDNA_50bp_overhang.fasta into ~/Library/ to align rDNA files"
	      exit 1
	  else
	      cd ~/Library/Bowtie1_2
	      bowtie-build -f ../SK1_rDNA_50bp_overhang.fasta rDNA_SK1
	      cd $DIR
	  fi
      fi
      # align to a single rDNA repeat 
      # (note: SK1_rDNA_50bp_overhang.fastA has 50bp of overlapping sequence 
      # so multiple mappings are allowed)
      bowtie -p 8 -a -q -X $X -S ~/Library/Bowtie1_2/rDNA_SK1 \
	 -1 $FASTQ1 -2 $FASTQ2 rDNA/$TAG-rDNA-raw.sam --verbose
      
      cd rDNA
      # to limit memory costs
      samtools view -F 4 -h $TAG-rDNA-raw.sam > $TAG-rDNA-raw2.sam
      rm $TAG-rDNA-raw.sam

      if [ ! -f ~/ChIPseq_Pipeline_v3/Position_shift_SAM.pl ]; then
	  echo "Cannot find Position_shift_SAM.pl. Should be in ~/ChIPseq_Pipeline_v3/"
	  exit 1
      fi
      # put sam outputs on chrXII
      if [[ $GENOME =~ "SacCer3" ]]; then
	  perl ~/ChIPseq_Pipeline_v3/Position_shift_SAM.pl -i $TAG-rDNA-raw2.sam \
	      -c chrXII -p 451418 -o $TAG-$ROOT.sam
      elif [[ $GENOME =~ "SK1K" ]]; then
	  perl ~/ChIPseq_Pipeline_v3/Position_shift_SAM.pl -i $TAG-rDNA-raw2.sam \
	      -c chr12 -p 433028 -o $TAG-$ROOT.sam
      elif [[ $GENOME =~ "SK1Yue" ]]; then
	  perl ~/ChIPseq_Pipeline_v3/Position_shift_SAM.pl -i $TAG-rDNA-raw2.sam \
	      -c chrXII -p 447013 -o $TAG-$ROOT.sam
      fi

      rm $TAG-rDNA-raw2.sam

   fi

else
	echo "Unknown genome listed."
	exit 1
fi

cd $DIR
test -f $FASTQ1.gz || gzip $FASTQ1
test -f $FASTQ2.gz || gzip $FASTQ2

echo "complete"
date
echo ""



