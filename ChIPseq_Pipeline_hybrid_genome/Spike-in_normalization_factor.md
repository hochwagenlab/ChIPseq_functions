# Calculating spike-in normalization factors

Once you've sequenced your SK1 samples spiked with an S288C background-based strain and aligned the resulting reads to the SK1-S288C hybrid genome (using the "Hybrid genome ChIP-seq pipeline" here), you are ready to compute the spike-in normalization factors. The first step is to get counts of reads aligned to each genome.

Let's take two example conditions, wild type strain AH119 and low-Red1 strain AH9048, with input and IP samples for each. Running the hybrid genome pipeline produced BAM alignment maps, with their respective index files.

```bash
ah119_spike_chip_S288c_SK1_Yue-PM_sorted.bam
ah119_spike_chip_S288c_SK1_Yue-PM_sorted.bam.bai
ah119_spike_input_S288c_SK1_Yue-PM_sorted.bam
ah119_spike_input_S288c_SK1_Yue-PM_sorted.bam.bai
ah9048_spike_chip_S288c_SK1_Yue-PM_sorted.bam
ah9048_spike_chip_S288c_SK1_Yue-PM_sorted.bam.bai
ah9048_spike_input_S288c_SK1_Yue-PM_sorted.bam
ah9048_spike_input_S288c_SK1_Yue-PM_sorted.bam.bai
```

We get aligned read counts by obtaining counts of reads aligning to each chromosome, using `samtools idxstats` with the produced alignment map BAM files as the input.


```bash
for ALN in *_sorted.bam
do
    echo ">>>>> $ALN"
    samtools idxstats $ALN | cut -f 1,3 > stats_${ALN%_sorted.bam}.txt
done
```

This produces the following text files:

```bash
stats_ah119_spike_chip_S288c_SK1_Yue-PM.txt
stats_ah119_spike_input_S288c_SK1_Yue-PM.txt
stats_ah9048_spike_chip_S288c_SK1_Yue-PM.txt
stats_ah9048_spike_input_S288c_SK1_Yue-PM.txt
```

Here's an example of what they look like.

```bash
$ more stats_ah119_spike_chip_S288c_SK1_Yue-PM.txt
chrI_S288C      28799
chrII_S288C     58882
chrIII_S288C    30758
chrIV_S288C     109397
chrV_S288C      45905
chrVI_S288C     35852
chrVII_S288C    80876
chrVIII_S288C   47754
chrIX_S288C     33365
chrX_S288C      53175
chrXI_S288C     44937
chrXII_S288C    140746
chrXIII_S288C   63952
chrXIV_S288C    51483
chrXV_S288C     81037
chrXVI_S288C    75778
chrI_SK1        119295
chrII_SK1       286200
chrIII_SK1      125943
chrIV_SK1       530850
chrV_SK1        221427
chrVI_SK1       191777
chrVII_SK1      390217
chrVIII_SK1     216634
chrIX_SK1       156239
chrX_SK1        254171
chrXI_SK1       224275
chrXII_SK1      340318
chrXIII_SK1     319701
chrXIV_SK1      261328
chrXV_SK1       383405
chrXVI_SK1      357341
*       0
```

We can now use these files in R to calculate a reference condition-centric normalization factor, which will be `1.0` for the reference condition (the wild type) and the computed normalization factor for the test condition. The following code is an example implementation.

```r
library(stringr)
library(readr)

#' Compute spike-in normalization factor from total read counts
#'
#' Computes spike-in normalization factor between two spiked-in samples using
#' total counts of aligned reads. Inputs paths to text files containing counts
#' of aligned reads per chromosome of a hybrid SK1:S288C genome.
#' @param ref_chip_counts Either a single or a list of paths to reference ChIP
#' samples' read counts file. No default.
#' @param ref_input_counts Either a single or a list of paths to reference input
#' samples' read counts file. No default.
#' @param test_chip_counts Either a single or a list of paths to test ChIP
#' samples' read counts file. No default.
#' @param test_input_counts Either a single or a list of paths to test input
#' samples' read counts file. No default.
#' @param return_counts Logical indicating whether to return the computed read
#' counts instead of the normalization factor. Defaults to \code{FALSE}.
#' @return Numeric normalization factor.
#' @examples
#' \dontrun{
#' spikein_normalization_factor_from_counts(
#'      ref_chip_counts='Counts_AH119_chip.txt',
#'      ref_input_counts='Counts_AH119_input.txt',
#'      test_chip_counts='Counts_AH8104_chip.txt',
#'      test_input_counts='Counts_AH8104_input.txt')
#'
#' spikein_normalization_factor_from_counts(
#'     ref_chip_counts=list('Counts_AH119_chip_1.txt',
#'                          'Counts_AH119_chip_2.txt',
#'                          'Counts_AH119_chip_3.txt'),
#'     ref_input_counts=list('Counts_AH119_inp_1.txt',
#'                           'Counts_AH119_inp_2.txt',
#'                           'Counts_AH119_inp_3.txt'),
#'     test_chip_counts='Counts_AH8104_chip.txt',
#'     test_input_counts='Counts_AH8104_input.txt')
#' }
#' @export
spikein_normalization_factor_from_counts <- function(
  ref_chip_counts, ref_input_counts, test_chip_counts, test_input_counts,
  return_counts=FALSE) {

  # Put paths in list
  files <- list(ref_chip=ref_chip_counts, ref_input=ref_input_counts,
                test_chip=test_chip_counts, test_input=test_input_counts)

  # Convert each element into list, if not one already
  for (i in seq_along(files)) {
    if (!is.list(files[[i]])) files[[i]] <- list(files[[i]])
  }

  # Print files to read to console
  message('>>> Read alignment count files:')
  for (i in seq_along(files)) {
    for (file in files[[i]]) {
      message('   ', basename(file))
    }
  }    

  message()
  # Read files into tibble in list
  tables <- list()
  for (i in seq_along(files)) {
    tables[[i]] <- sapply(files[[i]], FUN=read_tsv, col_names=F,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(tables) <- names(files)

  message()
  # Get read counts per chromosome
  message('>>> Count reads per genome:')
  counts <- list()
  for (i in seq_along(tables)) {
    counts[[i]] <- sapply(tables[[i]], FUN=sum_per_genome,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(counts) <- names(tables)

  # Add-up counts for replicates (results in nested lists)
  for (i in seq_along(counts)) {
    if (length(counts[[i]]) > 1) {
      total <- counts[[i]][[1]]
      for (j in 2:length(counts[[i]])) {
        total <- total + counts[[i]][[j]]
      }
      counts[[i]] <- total
    } else counts[[i]] <- unlist(counts[[i]])
  }

  if (return_counts) {
    message('---')
    message('Done!')
    return(counts)
  }

  # Compute normalization factor
  result <- normalization_factor(ctrl_input=counts$ref_input,
                                 ctrl_chip=counts$ref_chip,
                                 test_input=counts$test_input,
                                 test_chip=counts$test_chip)

  message('---')
  message('Done!')

  return(result)
}

# Helper functions
sum_per_genome <- function(df) {
  # Compute sum of reads aligned to each genome
  S288C <- sum(
    df[apply(df, 1, function(x) str_detect(x[1],'_S288C')), 2])
  SK1 <- sum(
    df[apply(df, 1, function(x) str_detect(x[1], '_SK1')), 2])

  # Print result to console
  message('  S288C: ', formatC(S288C, big.mark=",",
                               drop0trailing=TRUE, format="f"))
  message('  SK1: ', formatC(SK1, big.mark=",",
                             drop0trailing=TRUE, format="f"))
  message('      ', round(S288C * 100 / (SK1 + S288C), 1), '% spike-in reads')

  # Return result as named vector
  c('S288C'=S288C, 'SK1'=SK1)
}


normalization_factor <- function(ctrl_input, ctrl_chip,
                                 test_input, test_chip) {
  # Compute Q values
  Q_ctrl_input <- ctrl_input['S288C'] / ctrl_input['SK1']
  Q_ctrl_chip <- ctrl_chip['S288C'] / ctrl_chip['SK1']

  Q_test_input <- test_input['S288C'] / test_input['SK1']
  Q_test_chip <- test_chip['S288C'] / test_chip['SK1']

  # Compute normalization factors
  a_ctrl <- Q_ctrl_input / Q_ctrl_chip
  a_test <- Q_test_input / Q_test_chip

  # Return reference strain-centric normalization factor
  a_test/ a_ctrl
}
```

We can use this code to calculate the result.

```r
> read_counts <- data.frame(
  Condition=c('WT', 'Low Red1'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_ah119_spike_chip_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_ah119_spike_input_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_ah9048_spike_chip_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_ah9048_spike_input_S288c_SK1_Yue-PM.txt')
     )
)
```

```
>>> Read alignment count files:
   stats_ah119_spike_chip_S288c_SK1_Yue-PM.txt
   stats_ah119_spike_input_S288c_SK1_Yue-PM.txt
   stats_ah9048_spike_chip_S288c_SK1_Yue-PM.txt
   stats_ah9048_spike_input_S288c_SK1_Yue-PM.txt

Parsed with column specification:
cols(
  X1 = col_character(),
  X2 = col_integer()
)
Parsed with column specification:
cols(
  X1 = col_character(),
  X2 = col_integer()
)
Parsed with column specification:
cols(
  X1 = col_character(),
  X2 = col_integer()
)
Parsed with column specification:
cols(
  X1 = col_character(),
  X2 = col_integer()
)

>>> Count reads per genome:
  S288C: 982,696
  SK1: 4,379,121
      18.3% spike-in reads
  S288C: 115,044
  SK1: 441,130
      20.7% spike-in reads
  S288C: 1,036,418
  SK1: 1,470,564
      41.3% spike-in reads
  S288C: 400,313
  SK1: 1,574,041
      20.3% spike-in reads
---
Done!
```

This produces the following `data.frame`.

```r
> read_counts

Condition        NF
       WT 1.0000000
 Low Red1 0.3105042
```

The resulting normalization factor reflects the percentage of target protein in that condition relative to the wild type, 31.05%. To compare any data for the two conditions directly, data (that has not been normalized by sequencing depth) for the test condition should be multiplied by the calculated normalization factor. We will typically __normalize the data by genome-wide average__, to center it at `1.0` for all conditions, and then multiply by the appropriate normalization factor.
