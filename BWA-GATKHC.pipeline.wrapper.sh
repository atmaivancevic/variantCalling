#!/bin/bash

# A script to set off a series of parallel mapping and calling from a folder of fastq files.

# Example usage:
#./BWA-GATKHC.pipeline.wrapper.sh \
#-s /data/neurogenetics/sequences/Illumina/WES/KruerMarch2017 \
#-p /fast/users/a1211880/references/idlists/KruerMarch2017.txt \
#-o /fast/users/a1211880/outputs/variantCalling

usage()
{
echo "# A script to set off a series of parallel mapping and calling from a folder of fastq files.
# Requires: BWA-GATKHC.q and dependencies.  
#
# Usage $0 -s /path/to/sequences [-p /path/to/prefix-file.txt][ -o /path/to/output] | [ - h | --help ]
#
# Options
# -s	Path to the sequence files
# -p	Path and file name of a list of prefixes for the sequence files to run if not specified all sequences in the folder will be run
# -o	Path to where you want to find your file output (if not specified current directory is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 27/03/2014
# Modified: Atma Ivancevic, 9/5/2017, translating to SLURM
# Modified: (Date; Name; Description)
#
"
}
## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-s )			shift
					SEQPATH=$1
					;;
		-o )			shift
					WORKDIR=$1
					;;
		-p )			shift
					prefixFile=$1
					;;
		-h | --help )		usage
					exit 1
					;;
		* )			usage
					exit 1
	esac
	shift
done

if [ -z "$SEQPATH" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "#ERROR: You need to specify the path to your sequence files"
	exit 1
fi
if [ -z "$WORKDIR" ]; then # If no output directory then use current directory
	WORKDIR=$(pwd)
	echo "Using current directory as the working directory"
fi

if [ -z "$prefixFile" ]; then # If a list of sequences was not supplied then extract from seq folder
	prefixList=$(ls $SEQPATH | awk -F "_" '{print $1}' | sort | uniq) 
else
	prefixList=$(cat $prefixFile)
fi

for p in $prefixList
do
	 OUTPREFIX=$p SEQPATH=$SEQPATH WORKDIR=$WORKDIR sbatch BWA-GATKHC.q -J BWA-GATKHC.$p
done
