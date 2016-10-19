#########################################################################
#
# File name: NormalizeWiggle.R
# Written by: Tovah Markowitz
# Creation date: 10/18/16
# Purpose: To median normalize or double normalize wiggle data
# Functions included: readAll.tab, medianCalc, normalize, doubleNormalize,
# writeWiggle, and main
# Dependencies:
# Latest Update:
# By:
#
## ARGUMENTS
# Note: All wiggle files must be unzipped on linux computers
# experiment           - wiggle file of experimental data to be normalized (if single chromosome) or
#                        folder containing 16 or 17 wiggle files [if 17 assumes first file is all
#                              chromosomes] REQUIRED
# control              - wiggle file of experimental data to be normalized (if single chromosome) or
#                        folder containing 16 or 17 wiggle files [if 17 assumes first file is all
#                              chromosomes] REQUIRED
# normalizationType    - Can be any name that might be added to file name
#                        Also used to control some features of function
#                        If blank, genome-wide median will be calculated and used for normalization.
#                              Files will gain the words "wiggle_norm"
#                        To only calculate median values and then exit, type: "MedianOnly"
#                        To double normalize data, type: "double_norm"
# exptMedian           - To use an externally provided median on the experimental wiggle
# controlMedian        - To use an externally provided median on the control wiggle
# filemakerIDs         - To incorporate filemakerID numbers of experiments into the output files
#                               Must be a dash-delimited string. If a number is already within name,
#                               it is used to position new IDs within the name.
#
#
## EXAMPLES
# Example 1: MACS2 basic median normalization
# R CMD BATCH "--args AH119B-053013-SacCer3-2mis-PM-M5_MACS2_treat_pileup AH119B-053013-SacCer3-2mis-PM-M5_MACS2_control_lambda" NormalizeWiggle.R
# Example 2: MACS14 basic median normalization
# R CMD BATCH "--args AH119B-053013-SacCer3-2mis_P15_MACS2_wiggle/treat AH119B-053013-SacCer3-2mis_P15_MACS2_wiggle/control" NormalizeWiggle.R
# Example 3: MACS2 rDNA normalized with genome-wide medians
# R CMD BATCH "--args AH119B-053013-SacCer3-2mis-rDNA_MACS2_treat_pileup_chrXII.wig AH119B-053013-SacCer3-2mis-rDNA_MACS2_control_lambda_chrXII.wig rDNA 150 141" NormalizeWiggle.R
# Example 4: MACS14 rDNA normalized with genome-wide medians
# R CMD BATCH "--args AH119B-053013-SacCer3-2mis_P15_MACS2_wiggle/treat/AH119B-053013-SacCer3-2mis_P15_treat_afterfiting_chrXII.wig AH119B-053013-SacCer3-2mis_P15_MACS2_wiggle/control/AH119B-053013-SacCer3-2mis_P15_control_afterfiting_chrXII.wig rDNA 150 141" NormalizeWiggle.R
# Example 4: just calculating median values
# R CMD BATCH "--args AH119B-053013-SacCer3-2mis-PM-M5_MACS2_treat_pileup AH119B-053013-SacCer3-2mis-PM-M5_MACS2_control_lambda MedianOnly" NormalizeWiggle.R
# Example 5: double normalization
# R CMD BATCH "--args AH6408I-144-183-reps-SacCer3-2mis-PM-M5_MACS2_wiggle_norm AH7797I-148-185-reps-SacCer3-2mis-PM-M5_MACS2_wiggle_norm double_norm NA NA 144-183-148-185" ../NormalizeWiggle.R  

#########################################################################
### FUNCTION
#
# Function name: readAll.tab
# Purpose: read in whole genome wiggle data
# Arguments/Inputs: fileLocation
# Outputs: list of wiggle file data, excluding headers
# Assumptions: 16 or 17 .wig files in folder and nothing else
#              variable step wiggle with two columns
#              header information is only within the first 2 columns
# Example usage:

readAll.tab <- function( fileLocation ) {
    print("Reading in data")
    filenames <- list.files( fileLocation, full=T )
    # check if 'all chromosomes' file exists
    if ( length(filenames) == 17 ) {
        filenames <- filenames[2:17]
    }
    if ( sum( grepl( "\\.wig", filenames) ) == 16 ) {
        alldata <- lapply(filenames, read.table, skip=2, sep="\t" )
        names(alldata) = filenames
    } else {
        stop( "At least one file is not a wiggle file." )
    }
    return (alldata)
}

#########################################################################
### FUNCTION
#
# Function name: medianCalc
# Purpose: median calculation of wiggle data
# Arguments/Inputs: list of 1 or 16 chromosomes (such as output from readAll.tab)
# Outputs: the median value
# Assumptions:
# Example usage:

medianCalc <- function( data ) {
    print("calculating mean")
    if (length(data) == 1) { ## only one chromosome
        med <- median( data[,2] )
    } else if ( length(data) == 16 ) {  ## all chromosomes
        med <- median( c( data[[1]][,2], data[[2]][,2], data[[3]][,2], data[[4]][,2],
                         data[[5]][,2], data[[6]][,2], data[[7]][,2], data[[8]][,2],
                         data[[9]][,2], data[[10]][,2], data[[11]][,2], data[[12]][,2],
                         data[[13]][,2], data[[14]][,2], data[[15]][,2], data[[16]][,2]) )
    } else {
        stop( paste0( "Cannot calculate median value. Expecting either 1 or 16 chromosomes. Got ",
                     length(data), " chromosomes." ) )
    }
    return( med )
}

#########################################################################
### FUNCTION
#
# Function name: normalize
# Purpose: to median normalize data
# Formula at each position is: ChIP/(ChIP median) - Input/(Input median)
# Arguments/Inputs: two lists of wiggle data and two median values
# Outputs: normalized wiggle data as list
# Assumptions: assumes every base in ChIP is found in Input and vice versa
# Example usage:

normalize <- function( ChIP, Input, ChIPmedian, Inputmedian ) {
    print("median normalizing")
    # Formula at each position is: ChIP/(ChIP median) - Input/(Input median)
    out <- list()
    for ( i in 1:length(ChIP) ) {
	ChIP[[i]][,2] <- ChIP[[i]][,2] / ChIPmedian
	Input[[i]][,2] <- Input[[i]][,2] / Inputmedian
	out[[i]] <- data.frame( matrix( NA, nrow=dim(ChIP[[i]][1]), ncol=2 ) )
	out[[i]][,1] <- ChIP[[i]][,1]
	out[[i]][,2] <- ChIP[[i]][,2] - Input[[i]][,2]
    }
    return (out)
}

#########################################################################
### FUNCTION
#
# Function name: doubleNormalize
# Purpose: to normalize data that has already been median normalized by subtraction
# Arguments/Inputs: two lists of wiggle data
# Outputs: list of normalized wiggle data
# Assumptions: assumes every base in experiment is found in control and vice versa
# Example usage:

doubleNormalize <- function( experiment, control ) {
    print("double normalizing")
    out <- experiment
    for ( i in 1:length(experiment) ) {
        out[[i]][,2] <- experiment[[i]][,2] - control[[i]][,2]
    }
    return (out)
}

#########################################################################
### FUNCTION
#
# Function name: adjustHeaders
# Purpose: to adjust headers to include control file (if applicable), as well as
# normalization type, and median values
# Arguments/Inputs: many, including headers from wiggle files as list with each line being a row
# Outputs: new and improved headers in list format
# Assumptions: 
# Example usage:

adjustHeaders <- function( headers, experiment, control, normalizationType, MACS, exptMedian, controlMedian ) {
    # to create a new descriptor that describes most aspects of this function
    part1 <- strsplit( tail( strsplit( experiment, split="/" )[[1]], 1), split="_" )[[1]][1]
    part2 <- strsplit( tail( strsplit( control, split="/" )[[1]], 1), split="_" )[[1]][1]
    if ( part1 == part2 ) {
        exptName <- part1
    } else {
        exptName <- paste(part1, part2, sep="_")
    }
    
    if ( is.na( exptMedian ) ) {
        newdesc <- paste( exptName, MACS, normalizationType, sep="_" )
    } else {
        newdesc <- paste( exptName, MACS, normalizationType, "exptMedian", exptMedian,
                         "controlMedian", controlMedian, sep="_" )
    }

    for( i in 1:length(headers) ) {
    # to replace old header with new header
        tmp <- headers[[i]][1,]
        tmp <- strsplit( tmp, "=" )[[1]]
        if ( MACS == "MACS2" ) {
        # for each header replace description with new description while leaving name unchanged
            tmp[4] <- newdesc 
        } else if ( MACS == "MACS14" ) {
            tmp[3] <- newdesc
        }
        headers[[i]][1,] <- paste( tmp, collapse= "=" )
    }
    return (headers)
}
    

#########################################################################
### FUNCTION
#
# Function name: writeWiggle
# Purpose: to write output wiggle data to files
# Arguments/Inputs: list of normalized wiggle data, list of header data, vector of output filenames
# Outputs: .wig files for every chromosome, plus a single file with all chromosomes (if applicable)
# Assumptions:
# Example usage:

writeWiggle <- function( normData, headers, fileNames ) {
    print("writing output")
    if ( length(normData) == 1 ) {
        if ( length(headers) != 1 || length(fileNames) != 1 ) {
            stop( paste0("Normalized data has ", length(normData), " chromosomes./n",
                         "Headers have ", length(headers), " chromosomes./n",
                         "Output filenames have ", length(fileNames), " chromosomes./n",
                         "All should have had the same length." ) )
        } else {
            write.table( headers[[1]], file=fileNames, sep="\n", quote=F, row.names=F, col.names=F )
            write.table( normData[[1]], file=fileNames, sep="\t", quote=F, row.names=F, col.names=F, append=T)
        }
    } else {
        if ( length(normData) != length(headers) || length(normData) != ( length(fileNames) - 1 ) ) {
            stop( paste0("Normalized data has ", length(normData), " chromosomes./n",
                         "Headers have ", length(headers), " chromosomes./n",
                         "Output filenames have ", (length(fileNames) - 1), " chromosomes./n",
                         "All should have had the same length." ) )
        } else {
            for ( i in 1:length(normData) ) {
                write.table( headers[[i]], file=fileNames[i], sep="\n", quote=F, row.names=F, col.names=F )
                write.table( normData[[i]], file=fileNames[i], sep="\t", quote=F, row.names=F, col.names=F, append=T)
                if (i==1) {
                    write.table( headers[[i]], file=fileNames[17], sep="\n", quote=F, row.names=F, col.names=F)
                } else {
                    write.table(headers[[i]], file=fileNames[17], sep="\n", quote=F, row.names=F, col.names=F, append=T)
                }
                write.table( normData[[i]],file=fileNames[17], sep="\t", quote=F, row.names=F, col.names=F, append=T)
            }
        }
    }
}

# Found online
# Purpose: to convert any "NA"s given as arguments into NA without affecting anything else
make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
                                  is.na(x) <- x=="NA"; x} else {
                                                            x}



#########################################################################
### MAIN

main <- function( experiment, control, normalizationType=NA, exptMedian=NA, controlMedian=NA, filemakerIDs=NA ) {

    normalizationType <- make.true.NA( normalizationType )
    exptMedian <- make.true.NA( exptMedian )
    controlMedian <- make.true.NA( controlMedian )
    filemakerIDs <- make.true.NA (filemakerIDs )
    
    # check if experiment or control are files
    test1 <- grepl( "\\.", experiment )
    test2 <- grepl( "\\.", control )
    if ( test1 != test2 ) {
        stop( "Experiment and control are not the same file type. They need to either both be files or both be folders." )

    } else if ( test1 == TRUE ) {
        # if files
        if ( sum( grepl( "\\.wig", c(experiment, control) ) ) == 2 ) {
            # if both files are wiggle files
            exptData <- read.table( experiment, skip=2 )
            controlData <- read.table( control, skip=2 )
            headers <- read.table( experiment, nrows=2, sep="\n", stringsAsFactors= FALSE )
             # To remove all positions which do not have at least one read aligned for BOTH treatment and control
            controlData <- list( na.omit( controlData[match( exptData[,1], controlData[,1] ),] ) )
            exptData <- list( na.omit( exptData[match( controlData[,1], exptData[,1] ),] ) )
        } else {
            stop( "At least one file is not a wiggle file." )
        }

    } else {
        # if folders
        exptData <- readAll.tab( experiment )
        controlData <- readAll.tab( control )
        exptInFiles <- list.files( experiment, full=T )
        # check if 'all chromosomes' file exists
        if ( length(exptInFiles) == 17 ) {
            exptInFiles <- exptInFiles[2:17]
        }
        headers <- lapply( exptInFiles, read.table, nrows=2, sep="\n", stringsAsFactors= FALSE )
        # To remove all positions which do not have at least one read aligned for BOTH treatment and control
        for ( i in 1:length(controlData) ) {
            controlData[[i]] <- na.omit( controlData[[i]][match( exptData[[i]][,1], controlData[[i]][,1] ),] )
            exptData[[i]] <- na.omit( exptData[[i]][match( controlData[[i]][,1], exptData[[i]][,1] ),] )
        }
    }
    
    if (is.na(normalizationType) | normalizationType != "double_norm" ) {
        # if median values were not already included
        if ( is.na(exptMedian) & is.na(controlMedian) ) {
            controlMedian <- medianCalc( controlData )
            exptMedian <- medianCalc( exptData )
            print ( paste0( "Median value for control (", control, ") is: ", controlMedian ) )
            print ( paste0( "Median value for experiment (", experiment, ") is: ", exptMedian ) )
        } else if (  is.na(exptMedian) | is.na(controlMedian) ) {
            stop( "Be careful: Either include two median values or none." )
        } else {
            exptMedian <- as.numeric(exptMedian)
            controlMedian <- as.numeric(controlMedian)
            if ( is.na(normalizationType) ) {
                normalizationType <- "genomeMedian_norm"
            } else {
                normalizationType <- paste0( "genomeMedian_", normalizationType )
            }
        }
    }

    # normalize data if requested
    if ( !( is.na(normalizationType) ) ) {
        if ( normalizationType == "MedianOnly" ) {
            stop( "Only median values requested. Finished." )
        } else if ( normalizationType == "double_norm" ) {
            outData <- doubleNormalize( exptData, controlData )
        } else {
            outData <- normalize( exptData, controlData, exptMedian, controlMedian )
        }
    } else {
        outData <- normalize( exptData, controlData, exptMedian, controlMedian )
        normalizationType <- "wiggle_norm"
    }
    

    # determining names of output files/folders (part1)
    outname.base <- strsplit( tail( strsplit( experiment, split="/" )[[1]], 1), split="_" )[[1]][1]
    pwd <- getwd()
    # if MACS14
    if ( outname.base == "treat" ) {
        outname.base <- strsplit( tail( strsplit( experiment, split="/" )[[1]], 2)[1], split="_" )[[1]][1]
        MACS <- "MACS14"
    } else {
        MACS <- "MACS2"
    }

    # dealing with filemakerIDs
    if ( !( is.na(filemakerIDs) ) ) {
        indIDs <- strsplit( filemakerIDs, "-" )[[1]]
        knownID <- ""
        unknownID <- ""
        for (i in 1:length(indIDs) ) {
           tmp <- paste0( "-", indIDs[i], "-" )
            if ( grepl( tmp, outname.base ) ) {
                knownID <- c( knownID, indIDs[i] )
            } else {
                unknownID <- c( unknownID, indIDs[i] )
            }
        }
        if ( length(unknownID) != 1 ) {
            if ( length(knownID) != 1 ) {
                outname.base.part <- strsplit( outname.base, knownID[2] )[[1]]
                newIDs <- paste( c(knownID[2], unknownID[2:length(unknownID)]), collapse= "-" )
            } else {
                outname.base.part <- strsplit( outname.base, "-S")
                tmp <- regexpr("-S\\w+", outname.base, perl=T )
                newIDs <- paste( c(unknownID[2:length(unknownID)], regmatches( outname.base, tmp )), collapse= "-" )
            }
            if ( length(outname.base.part) == 2 ) {
                outname.base <- paste( c( outname.base.part[1], newIDs, outname.base.part[2]), collapse="" )
            } else {
                outname.base <- paste0( c(outname.base.part[1], newIDs, outname.base.part[2:length(outname.base.part)]),
                                       collapse="")
            }
        }                                                
    }
    
    if ( test1 == TRUE ) {
    # to get final file name for output when single file
        chr <- strsplit( tail( strsplit( experiment, split="_" )[[1]], n=1), split="\\." )[[1]][1]
        fileNames <- paste0( outname.base, "_", MACS, normalizationType, "_", chr, ".wig" )
    } else {
        outRoot <- paste0( outname.base, "_", MACS, "_", normalizationType )
        if ( normalizationType != "wiggle_norm" ) {
            outFolder <- paste0( outname.base, "_", MACS, "_wiggle_", normalizationType )
        } else {
            outFolder = outRoot
        }
        dir.create( outFolder )
        setwd( paste0( pwd, "/", outFolder ) )
        for ( i in 1:length(exptInFiles) ) {
            chr <- strsplit( tail( strsplit( exptInFiles[i], split="_" )[[1]], n=1), split="\\." )[[1]][1]
            if ( i == 1 ) {
                fileNames <- paste0( outRoot, "_", chr, ".wig" )
            } else {
                fileNames <- c( fileNames, paste0( outRoot, "_", chr, ".wig" ) )
            }
        }
        fileNames <- c( fileNames, paste0( outRoot, "_all.wig" ) )
    }
    
    headers <- adjustHeaders(headers, experiment, control, normalizationType, MACS, exptMedian, controlMedian )
    # write output wiggle files
    writeWiggle( outData, headers, fileNames )
    setwd( pwd )
}

args<-commandArgs(TRUE)

main(args[1],args[2],args[3],args[4],args[5],args[6])
# args[1]: experiment
# args[2]: control
# args[3]: normalizationType
# args[4]: exptMedian
# args[5]: controlMedian
# args[6]: filemakerIDs
