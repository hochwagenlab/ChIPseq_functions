#!/usr/bin/perl

# Name:ExtractMedian_MACS2_pipeline.pl
# Author: Tovah Markowitz
# Date: 9/30/16
# Updated: 11/13/16 to work with NormalizeWiggle.R
# Input: folder where NormalizeWiggle.Rout file(s) is located
# Output: a tab delimited file in home with folder name and genome-wide ChIP and input medians

use strict;
use Data::Dumper;
use Getopt::Long;                                                                             
use Pod::Usage;

my $help=0;
my ($folderName, $user);

# set commandline options.                                        
GetOptions ("folderName=s" => \$folderName,
	    "user=s" => \$user,
            "help" => \$help,) or pod2usage(2);
pod2usage(1) if $help;


################################################################

# find all NormalizeWiggle.Rout files in directory
opendir (DIR, $folderName), or die "Cannot open directory: $!";
my @files = grep { /NormalizeWiggle.Rout$/ && -f "$folderName/$_" } readdir(DIR);
closedir DIR;

my $inFile;
my $fName;
my $tmp;

# for each file found
for (my $count=0; $count < scalar(@files); $count++) {
  if ( $folderName =~ m(/$) ) {
      $tmp = $`;
      if ($files[$count] =~ /(.*)_NormalizeWiggle.Rout/ ) {
	  $fName = $1;
      } else {
	  $fName = $tmp;
      }
      $inFile = $folderName . $files[$count];
  } else {
    if ($files[$count] =~ /(.*)_NormalizeWiggle.Rout/ ) {
	$fName = $1;
    } else {
	$fName = $folderName;
    }
    $inFile = $folderName . "/" . $files[$count];
  }

# ensure input file exists else die
  open (IN, "<$inFile") or die ("Can not find your input file $inFile: No such file or directory\n");

  my $controlMedian;
  my $treatMedian;

  # for as long as the file has more data:
  while (<IN>) {
    chomp;
    
    if ($_ =~ /"Median value for control \w* is: ([0-9.]*)"/ ) {
      $controlMedian =  $1;
    } elsif ($_ =~ /"Median value for experiment \w* is: ([0-9.]*)"/ ) {
      $treatMedian =  $1;
    }
    
  }
  
  close IN;
  
  my $outFile = "/home/". $user ."/Genomewide_medians_MACS2.txt";

  # if no output file exists, write header to file
  unless (-e $outFile) {
    open (OUT, ">$outFile") or die ("Can not find your output file $outFile: No such file or directory\n");
    print OUT "SourceFolder" . "\t" . "InputMedian" . "\t" . "ChIPMedian" . "\n";
    close OUT;
  }

  open (OUT, ">>$outFile") or die ("Can not find your output file $outFile: No such file or directory\n");
  print OUT $fName . "\t" . $controlMedian . "\t" . $treatMedian . "\n";
  close OUT;

}

################################################################

__END__

=head1 SYNOPSIS

ExtractMedian_MACS2_pipeline.pl

Purpose: To get genome-wide median values for ChIP-seq data after being run through
and before normalization.

Output: A tab-delimited file in home with folder name and genome-wide ChIP 
and input medians

Options:

    -h, --help        brief help message
    -f, --folderName  name of folder with Rout files of interest
    -u, --user        User name

=cut
