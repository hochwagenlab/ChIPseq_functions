# README
# ChIP-seq Pipeline version 3

## PREPARATION
1. Move entire folder into your home on prince (keep folder name as  
   "ChIPseq_Pipeline_v3")
2. Move the 3 fasta files into a folder called "Library" in your home on prince
3. Open each sbatch script and adjust the email address
4. Start running jobs (sbatch ~/ChIPseq_Pipeline/ChIPseq-pipeline_v3.sbatch)
5. Keep an eye out for the following error:  
   **slurmstepd: error: Exceeded step memory at some point.**  
   If you see this in association with a Bowtie job, adjust the memory  
   requirements for Bowtie_B3.sbatch and rerun pipeline. You should be  
   able to ignore it for other portions of the pipeline.

## SOME RUN EXAMPLES
Example 1: Start with input file A.fastq and ChIP file B.fastq, map to SacCer3, 
and create wiggle files and both narrow and broad peak files 
(will also create wiggle plot of the rDNA)
```Bash
sbatch --export INPUT=A.fastq,CHIP=B.fastq,GEN="SacCer3",TAGI="A",TAGC="B",PEAK="BOTH" \ 
~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

Example 2: Get broad peak files for replicate input files A_1.sam and A_2.sam, and
ChIP files B_1.sam and B_2.sam
```Bash
sbatch --export INPUT="A_1.sam A_2.sam",CHIP="B_1.sam B_2.sam",REP="B",FLMKR="1-3",PEAK="BROAD",WIG="F" \
~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

Example 3: Get normalized wiggle files (and narrow peaks) with ChIP file A1.sam, input file A2.sam,
mock ChIP file B1.sam, and mock input file B2.sam
```Bash
sbatch --export INPUT=A2.sam,CHIP=A1.sam,INCON=B2.sam,CHCON=B1.sam,FLMKR="1-2" \
~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

## OUTPUT FILE NAMES
All outputs will get a version name identifier:
  - B3 means that the Sam file was created with this pipeline version
  - W3 means that the peak/bedgraph/wiggle files were made with  
    this pipeline version
  - B3W3 means that the everything was made with this pipeline version

### SAM:
All new sam files will be given the following names: $TAG-$GENROOT.sam
where $TAG is the user defined name and $GENROOT is one of the following:
  - SacCer3-2mis_B3
  - SK1K-PM_B3
  - SacCer3-rDNA_B3
  - SK1K-rDNA_B3

### MACS2 FOLDER:
Macs2 files will now go in folders with the notation: $TAGC-$GEN-$VER-MACS2 where  
  $TAGC is either the tagname of the ChIP file or the user-defined name $REP,  
  $GEN is either SacCer3 or SK1K, and  
  $VER will currently be either B3W3, 2mis_W3, or PM_W3.

$FLMKR will be incorporated into all MACS2 output files. 
The term "Reps" will be incorporated into all MACS2 output files when applicable.

### MACS2 FILES:
See MACS2 FOLDER with the following changes:
  - All files will include the method of Bowtie analysis: 2mis, PM, or rDNA.
  - Wiggle and bedgraph files will include the normalization method: FE.
  - ChIP vs ChIP analysis is identified as either: ChvCh or InvIn.

## ERROR/OUTPUT FILES
- The pipeline automatically makes two slurm_$JOBID.out files that I  
  haven't figured out how to remove yet. Both should be essentially  
  empty at all times. One will go to $FLDR and the other will go to the  
  working directory where you submit the job.  
  	  - I have changed one of these files to ChIPseq_closing_$JOBID.sbatch   
	  for troubleshooting purposes. It will now go to /scratch/$USER  
	  instead of $FOLDER
- The pipeline will initially make a /scratch/$USER/MACS2_pipeline_$CHIP.txt.  
  If none of the bowtie or macs2 jobs are initiated, this file will include  
  all necessary error messages.
- While running, the pipeline will create many .out files in  
  /scratch/$USER. These will eventually be consolidated along with  
  the above .txt into ChIPseq_Pipeline_$JOBID1.out where $JOBID1 is the  
  ID of the first job initiated from within the job.  
     	- This output file has been changed to $CHIP_$JOBID1.out.

## FOR TROUBLESHOOTING
If you get a job that has issues running, it may not consolidate all 
output files into a single document. If this is the case, follow these 
rules to begin troubleshooting:
1) Check the file /scratch/$USER/MACS2_pipeline_$CHIP.txt. It will tell   
   you all of the jobs that ran for this pipeline run as a colon-separated list.
2) For each job, just type: less /scratch/$USER/*[insert JobID here]* to   
   see what might have gone wrong. These may be called:  
       - Bowtie_$JOBID.out  
       - MACS2_FE_$JOBID.out  
       - ChIPseq_closing_$JOBID.out  

## NEW FEATURES
- Bowtie_B3.sbatch automatically checks for a file in both zipped and  
  unzipped forms, and proceeds accordingly. Therefore, inputs into the  
  pipeline can be written in either format.
- Bowtie_B3.sbatch includes a check for the databases needed, and if  
  not found, attempts to build them. No need for a separate  
  bowtie-build.
- When starting with two (or four for ChIP vs ChIP normalization)  
  fastq files, the pipeline automatically processes the rDNA along with  
  the rest of the pipeline.
- For details on the Bowtie or MACS2 conditions, see Bowtie_B3.sbatch  
  or MACS2_FE_W3.sbatch.

## KNOWN ISSUES
- There are currently no error checks for the following:
  - Do all SAM files have the same root?
  - Have wiggle or peak files already been created?
  - Is the FileMaker ID currently in the name? (Currently  
    there is only a superficial search).
  - Bowtie_B3.sbatch still takes 5 minutes even if it should only  
    be checking file names and unzipping/zipping a single fastq file.
  - Two extra slurm output files are created
  - No simple method to create an rDNA wiggle file from two  
       (or more) sam files.
