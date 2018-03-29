# Hybrid genome ChIP-seq pipeline

ChIP-seq experiment analysis code for samples containing hybrid S288C and SK1
genomes:

* SK1 samples spiked with SK288C cells
* S288C-SK1 hybrid strains

## Analysis pipeline Slurm job:

__`ChIPseq_Pipeline_hybrid_genome.sbatch`__

1. If GENNAME is provided:
    * Maps `FASTQ` data to indicated genome
    * Converts output from `SAM` to sorted and indexed `BAM` format
(allowing no mismatches)
2. Normalizes by library size using `MACS2` SPMR
3. Removes noise using `MACS2` fold enrichment compared to input

#### Argument options:

* __EXPID__     Custom ID for output files.
* __RUNDIR__    Path to directory to run script and save output in.
* __CHIP__      Absolute path(s) to ChIP (treatment) sample file(s).
* __INPUT__     Absolute path(s) to input (control) sample file(s).

    * NOTE: CHIP and INPUT file formats can be either:
        - __`FASTQ`__ if `GENNAME` is provided; pipeline will start with `Bowtie` mapping.
        - __`SAM`__ or __`BAM`__ if `GENNAME` is not provided; pipeline will skip `Bowtie` mapping.
 
* __GENNAME__   Basename of reference genome `FASTA` file preceded by absolute path to
directory containing the file. Must be provided if data files are in `FASTQ` format;
must not be provided if data files are alignments (`BAM` or `SAM` format).
An existing `Bowtie` index with a basename (`ebwt`) matching the file's is used if
found in the same directory; otherwise a new index is built.

### Running the pipeline:

To prepare the pipeline for running, follow the following steps:

* Clone the `ChIPseq_functions` GitHub repository into your preferred location:
    * `git clone https://github.com/hochwagenlab/ChIPseq_functions.git`
* Move into the hybrid genome pipeline folder:
    * `cd ChIPseq_Pipeline_hybrid_genome`
* Change the existing email address to your own in all files (replace "X" by
your user name):
    * `find . -type f | xargs sed -i 's/lv38@nyu.edu/X@nyu.edu/g'`
* Point to the location of the `sbatch` file when submittting the job (as in
the examples below) 

#### Example job submission:

The pipeline can be run on the raw `FASTQ` output of the sequencing run or on
previously made alignment maps (`BAM` or `SAM` files). In the latter case,
replicate sample files can be run simply by adding file names in the appropriate
variable separated by spaces, as in: `CHIP="path/to/chip1.bam path/to/chip2.bam"`.
If starting from `FASTQ` files, they must be run individually.

Starting from `FASTQ` files; no replicates:

```
sbatch --export EXPID="AH119spikein-060717_YueSK1_S288C_PM_SPMR",\
RUNDIR="/scratch/lv38",\
CHIP="/scratch/lv38/HLYHHAFXX_n01_ah119spikeb-062817.fastq.gz",\
INPUT="/scratch/lv38/HLYHHAFXX_n01_ah119spikea-062817.fastq.gz",\
GENNAME="/home/lv38/Library/S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue" \
~/ChIPseq_functions/ChIPseq_Pipeline_hybrid_genome/ChIPseq_Pipeline_hybrid_genome.sbatch
```

Starting from alignment files; no replicates:

```
sbatch --export EXPID="AH119spikein-060717_YueSK1_S288C_PM_SPMR",\
RUNDIR="/scratch/lv38",\
CHIP="/scratch/lv38/ah119spikeb-062817_S288C_SK1_Yue_PM.sam",\
INPUT="/scratch/lv38/ah119spikea-062817_S288C_SK1_Yue_PM.sam" \
~/ChIPseq_functions/ChIPseq_Pipeline_hybrid_genome/ChIPseq_Pipeline_hybrid_genome.sbatch
```

Starting from alignment files with replicates:

```
sbatch --export EXPID="Red1-WT-155-175-reps_S288C_SK1_Yue_PM_SPMR",\
RUNDIR="/scratch/lv38",\
CHIP="/scratch/lv38/ah119spikeb-062817_S288C_SK1_Yue_PM.sam \
/scratch/lv38/ah119spiked-01012018_S288C_SK1_Yue_PM.sam",\
INPUT="/scratch/lv38/ah119spikea-062817_S288C_SK1_Yue_PM.sam \
/scratch/lv38/ah119spiked-01012018_S288C_SK1_Yue_PM.sam" \
~/ChIPseq_functions/ChIPseq_Pipeline_hybrid_genome/ChIPseq_Pipeline_hybrid_genome.sbatch
```

## Spike-in normalization factor calculation

If you are running this pipeline for a spike-in experiment, you can find details
about how to calculate the normalization fator in the 
[dedicated page](https://github.com/hochwagenlab/ChIPseq_functions/blob/master/ChIPseq_Pipeline_hybrid_genome/Spike-in_normalization_factor.md).
