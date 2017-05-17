#!/bin/bash
# A script to map reads and then call variants using the GATK v3.x best practices designed for the eResearchSA supercomputer but will work on stand alone machines too

# Variables that usually don't need changing once set for your system
export FASTDIR=/fast/users/$USER

gVcfFolder=$FASTDIR/vcf/gVcfDumpingGround/Exomes # A place to dump gVCFs for later genotyping
BWAINDEXPATH=/data/neurogenetics/RefSeq/BWA/hg19_1stM_unmask_ran_all # Your genome reference path for BWA
BWAINDEX=hg19_1stM_unmask_ran_all.fa # name of the genome reference
GATKPATH=$FASTDIR/executables/GenomeAnalysisTK-3.7 # Where the GATK program.  Be mindful that GATK is under rapid development so things may change over time!
GATKREFPATH=/data/neurogenetics/RefSeq/GATK #Refseq index library locations
GATKINDEX=$BWAINDEX # Base name of GATK indexes (usually the same as the $BWAINDEX)
ChrIndexPath=$GATKREFPATH/$BWAINDEX.chridx #Location of index bed files
IndexBedFiles=01.hg19-M1.bed,02.hg19-2-3.bed,03.hg19-4-5.bed,04.hg19-6-7.bed,05.hg19-8-10.bed,06.hg19-11-13.bed,07.hg19-14-17.bed,08.hg19-18etc.bed # A comma separated array of names of index files
arrIndexBedFiles=$(echo $IndexBedFiles | tr "," "\n")
PICARDPATH=$FASTDIR/executables/Picard-2.9.2 # Where the picard program is.  Picard is also under rapid development so may change over time.
DBSNP=dbsnp_138.hg19.vcf
BUILD=$(echo $BWAINDEX | awk '{print substr($1, 1, length($1) - 3)}') # Genome build used = $BWAINDEX less the .fa, this will be incorporated into file names.

usage()
{
echo "# Script for processing and mapping Illumina 100bp pair-end sequence data and optionally plotting coverage for an interval
# Requires: BWA 0.7.x, Picard, samtools, GATKv3.x, BWA-Picard-GATK-CleanUp.sh.  
# This script assumes your sequence files are gzipped
#
# Usage $0 -p file_prefix -s /path/to/sequences [ -o /path/to/output] [-L LibraryName][-I ID][-i /path/to/bedfile.bed] | [ - h | --help ]
#
# Options
# -p	A prefix to your sequence files of the form PREFIX_R1.fastq.gz
# -s	Path to the sequence files
# -o	Path to where you want to find your file output (if not specified current directory is used)
# -S	Sample name if not specified then will be set the same as -p
# -L	Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20140317). Default \"IlluminaExome\"
# -I	ID for the sequence library (to go into the @RG line). If not specified script will make one up from the first read header
# --emu	Set parameters specific for the emu cluster in the cloud
# -i	Path and file name of .bed file with intervals of interest for coverage calculations (optional)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# System variables currently set:
# gVcfFolder=$gVcfFolder
# BWAINDEXPATH=$BWAINDEXPATH
# BWAINDEX=$BWAINDEX
# ChrIndexPath=$ChrIndexPath
# IndexBedFiles=$IndexBedFiles
# GATKPATH=$GATKPATH
# GATKREFPATH=$GATKREFPATH
# GATKINDEX=$GATKINDEX
# PICARDPATH=$PICARDPATH
# SCRIPTPATH=$SCRIPTPATH
# BUILD=$BUILD
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Contact: mark.corbett@adelaide.edu.au
# Modified: (Date; Name; Description)
# 18/09/2014; Mark Corbett; Update GATK version 3.2-2
# 21/11/2014; Mark Corbett; Incorporate GATKv3.x.HC.by.Quarters.HPC.sh 
# 03/12/2014; Mark Corbett; Update to use tizard Picard installation
# 15/01/2015; Mark Corbett; Put in temp file saves on GATK error
# 21/04/2016; Mark Corbett; Bring up to date with GATK current best practices
# 09/05/2017; Atma Ivancevic; Translating for Phoenix
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-p )			shift						
					OUTPREFIX=$1
					;;
		-s )			shift
					SEQPATH=$1
					;;
		-S )			shift
					SAMPLE=$1
					;;
		-o )			shift
					WORKDIR=$1
					;;
		-L )			shift
					LB=$1
					;;
		-I )			shift
					ID=$1
					;;
		-i )			shift
					INTERVAL=$1
					;;
		--emu )			shift
					gVcfFolder=/scratch/mcorbett/VCF/gVcfFileLibrary
					;;
		-h | --help )		usage
					exit 1
					;;
		* )			usage
					exit 1
	esac
	shift
done

if [ -z "$OUTPREFIX" ]; then # If no file prefix specified then do not proceed
	usage
	echo "#ERROR: You need to specify a file prefix (PREFIX) referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$SEQPATH" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "#ERROR: You need to specify the path to your sequence files"
	exit 1
fi
if [ -z "$WORKDIR" ]; then # If no output directory then use current directory
	WORKDIR=$(pwd)
	echo "Using current directory as the working directory"
fi
if [ -z "$SAMPLE" ]; then # If sample name not specified then use "OUTPREFIX"
	SAMPLE=$OUTPREFIX
	echo "Using $OUTPREFIX for sample name"
fi
if [ -z "$LB" ]; then # If library not specified then use "IlluminaExome"
	LB=IlluminaExome
	echo "Using \"IlluminaExome\" for library name"
fi
if [ -z "$INTERVAL" ]; then # If no interval set then 
	IntervalCoverage=false # Don't run coverage calculation
	else
	IntervalCoverage=true # If it is set then run the coverage calculation
fi
if [ ! -d $gVcfFolder ]; then
        mkdir -p $gVcfFolder
fi
tmpDir=/tmp/$USER/$OUTPREFIX # Use a tmp directory in tmp for all of the GATK and samtools temp files
if [ ! -d $tmpDir ]; then
	mkdir -p $tmpDir # -p means the parent directories will also be created if they don't exist
fi
	
# Locate sequence file names.
# This is a bit awkward and prone to errors since relies on only a few file naming conventions and assumes how they will line up after ls of files
cd $SEQPATH
SEQFILE1=$(ls | grep $OUTPREFIX\_ | head -n 1) # Assume sequence files are some form of $OUTPREFIX_fastq.gz
if [ -f $SEQFILE1 ]; then		 
	fileCount=$(ls | grep $OUTPREFIX\_ | wc -l | sed 's/[^0-9]*//g')
	if [ $fileCount -ne "2" ]; then
		echo "Sorry I've found the wrong number of sequence files and there's a risk I will map the wrong ones!"
		exit 1
	fi
	SEQFILE2=$(ls | grep $OUTPREFIX\_ | tail -n 1)
else
	fileCount=$(ls | grep -w $OUTPREFIX | wc -l | sed 's/[^0-9]*//g')
	if [ $fileCount -ne "2" ]; then
		echo "Sorry I've found the wrong number of sequence files and there's a risk I will map the wrong ones!"
		exit 1
	fi
	SEQFILE1=$(ls | grep -w $OUTPREFIX | head -n 1) 
	SEQFILE2=$(ls | grep -w $OUTPREFIX | tail -n 1)
fi
if [ ! -f $SEQFILE1 ]; then
	echo "Sorry I can't find your sequence files! I'm using $OUTPREFIX as the start of the filename to locate them"
	exit 1
fi
if [ -z $ID ]; then
	ID=$(zcat $SEQFILE1 | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}').$OUTPREFIX # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq
fi

## Start of the script ##
# Map reads to genome using BWA-MEM 
cd $tmpDir
bwa mem -M -t 8 -R "@RG\tID:$ID\tLB:$LB\tPL:ILLUMINA\tSM:$SAMPLE" $BWAINDEXPATH/$BWAINDEX $SEQPATH/$SEQFILE1 $SEQPATH/$SEQFILE2 |\
samtools view -bT $GATKREFPATH/$GATKINDEX - |\
samtools sort -l 5 -m 4G -@8 -T $OUTPREFIX -o $OUTPREFIX.samsort.bwa.$BUILD.bam -

# Mark duplicates
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $PICARDPATH/picard.jar MarkDuplicates \
INPUT=$OUTPREFIX.samsort.bwa.$BUILD.bam \
OUTPUT=$OUTPREFIX.marked.sort.bwa.$BUILD.bam \
METRICS_FILE=$OUTPREFIX.md.metrics \
MAX_RECORDS_IN_RAM=3000000 \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

rm $OUTPREFIX.samsort.bwa.$BUILD.bam

echo "# Flagstats" > $WORKDIR/$OUTPREFIX.Stat_Summary.txt
samtools flagstat $OUTPREFIX.marked.sort.bwa.$BUILD.bam >> $WORKDIR/$OUTPREFIX.Stat_Summary.txt
echo "
# MarkDuplicates metrics" >> $WORKDIR/$OUTPREFIX.Stat_Summary.txt
cat $OUTPREFIX.md.metrics >> $WORKDIR/$OUTPREFIX.Stat_Summary.txt 
rm $tmpDir/$OUTPREFIX.md.metrics

# Polish the BAM using GATK best practices
# Realign indels
# Create a file of intervals to realign to
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar \
-I $OUTPREFIX.marked.sort.bwa.$BUILD.bam \
-R $GATKREFPATH/$GATKINDEX \
-T RealignerTargetCreator \
-o $OUTPREFIX.forIndelRealigner.intervals \
-known $GATKREFPATH/$DBSNP \
-known $GATKREFPATH/hapmap_3.3.hg19.vcf \
-known $GATKREFPATH/1000G_omni2.5.hg19.vcf \
-known $GATKREFPATH/Mills_and_1000G_gold_standard.indels.hg19.vcf \
-rf BadCigar \
-nt 8 >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

# For the next steps split the bams into bits based on the IndexBedFiles
# First make tmp dirs
for bed in $arrIndexBedFiles; do
	mkdir -p $tmpDir/$bed
done

# Perform the local realignment around indels
for bed in $arrIndexBedFiles; do
	(
	java -Xmx6g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.jar \
	-I $OUTPREFIX.marked.sort.bwa.$BUILD.bam \
	-R $GATKREFPATH/$GATKINDEX \
	-T IndelRealigner \
	-L $ChrIndexPath/$bed \
	-known $GATKREFPATH/$DBSNP \
	-known $GATKREFPATH/hapmap_3.3.hg19.vcf \
	-known $GATKREFPATH/1000G_omni2.5.hg19.vcf \
	-known $GATKREFPATH/Mills_and_1000G_gold_standard.indels.hg19.vcf \
	-targetIntervals $tmpDir/$OUTPREFIX.forIndelRealigner.intervals \
	-rf BadCigar \
	-o $bed.$OUTPREFIX.realigned.bwa.$BUILD.bam >> $tmpDir/$bed.$OUTPREFIX.pipeline.log 2>&1
	) &
done
wait

cat $tmpDir/*.$OUTPREFIX.pipeline.log >> $WORKDIR/$OUTPREFIX.pipeline.log

ls | grep $OUTPREFIX.realigned.bwa.$BUILD.bam$ - > $tmpDir/$OUTPREFIX.bam.list.txt
sed 's,^,-I '"$tmpDir"'\/,g' $tmpDir/$OUTPREFIX.bam.list.txt > $tmpDir/$OUTPREFIX.inputBAM.txt

# Base quality score recalibration
java -Xmx30g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $GATKREFPATH/$GATKINDEX \
$(cat $tmpDir/$OUTPREFIX.inputBAM.txt) \
-knownSites $GATKREFPATH/$DBSNP \
-o $tmpDir/$OUTPREFIX.recal.grp \
-rf BadCigar \
-nct 8 >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

echo "
# BQSR metrics" >> $WORKDIR/$OUTPREFIX.Stat_Summary.txt
cat $tmpDir/$OUTPREFIX.recal.grp >> $WORKDIR/$OUTPREFIX.Stat_Summary.txt

# PrintReads
# Dumps the final .bam file in the project working directory.
# None of the previous .bam files will be kept!
for bed in $arrIndexBedFiles; do
	(
	java -Xmx6g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R $GATKREFPATH/$GATKINDEX \
	-I $bed.$OUTPREFIX.realigned.bwa.$BUILD.bam \
	-L $ChrIndexPath/$bed \
	-BQSR $tmpDir/$OUTPREFIX.recal.grp \
	-rf BadCigar \
	-o $bed.$OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bam > $tmpDir/$bed.$OUTPREFIX.pipeline.log 2>&1
	) &
done
wait

rm *$OUTPREFIX.realigned.bwa.$BUILD.bam
rm *$OUTPREFIX.realigned.bwa.$BUILD.bai
cat $tmpDir/*.$OUTPREFIX.pipeline.log >> $WORKDIR/$OUTPREFIX.pipeline.log

# Merge and sort with Picard
ls | grep $OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bam$ - > $tmpDir/$OUTPREFIX.bam.list.txt
sed 's,^,INPUT='"$tmpDir"'\/,g' $tmpDir/$OUTPREFIX.bam.list.txt > $tmpDir/$OUTPREFIX.inputBAM.txt

java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $PICARDPATH/picard.jar MergeSamFiles \
$(cat $tmpDir/$OUTPREFIX.inputBAM.txt) \
OUTPUT=$WORKDIR/$OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
USE_THREADING=true \
MERGE_SEQUENCE_DICTIONARIES=true \
MAX_RECORDS_IN_RAM=3000000 \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

# Clean up temporary .bam files
for bed in $arrIndexBedFiles; do
	(
	rm $bed.$OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bam
	rm $bed.$OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bai
	) &
done

wait

# As of GATK v3.x you can now run the haplotype caller directly on a single bam
# Run haplotype caller in gVCF mode

for bed in $arrIndexBedFiles; do
	java -Xmx4g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.jar \
	-I $WORKDIR/$OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bam \
	-R $GATKREFPATH/$GATKINDEX \
	-T HaplotypeCaller \
	-L $ChrIndexPath/$bed \
	--dbsnp $GATKREFPATH/$DBSNP \
	--min_base_quality_score 20 \
	--emitRefConfidence GVCF \
	-o $tmpDir/$bed.$OUTPREFIX.snps.g.vcf > $tmpDir/$bed.$OUTPREFIX.pipeline.log 2>&1 &
done
wait

cat *.$OUTPREFIX.pipeline.log >> $WORKDIR/$OUTPREFIX.pipeline.log
ls | grep $OUTPREFIX.snps.g.vcf$ > $OUTPREFIX.gvcf.list.txt
sed 's,^,-V '"$tmpDir"'\/,g' $OUTPREFIX.gvcf.list.txt > $OUTPREFIX.inputGVCF.txt

java -cp $GATKPATH/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R $GATKREFPATH/$GATKINDEX \
-out $OUTPREFIX.snps.g.vcf \
$(cat $OUTPREFIX.inputGVCF.txt) \
--assumeSorted >> $WORKDIR/$OUTPREFIX.pipeline.log  2>&1

bgzip $OUTPREFIX.snps.g.vcf
tabix $OUTPREFIX.snps.g.vcf.gz

mv $OUTPREFIX.snps.g.vcf.gz $gVcfFolder/$OUTPREFIX.snps.g.vcf.gz
mv $OUTPREFIX.snps.g.vcf.idx $gVcfFolder/$OUTPREFIX.snps.g.vcf.idx
mv $OUTPREFIX.snps.g.vcf.gz.tbi $gVcfFolder/$OUTPREFIX.snps.g.vcf.gz.tbi

## Other optional metrics ##
# Coverage for an interval
if $IntervalCoverage ; then
	samtools depth -b $INTERVAL $WORKDIR/$OUTPREFIX.realigned.recal.sorted.bwa.$BUILD.bam > $OUTPREFIX.Coverage.$BUILD.txt
	gzip $OUTPREFIX.Coverage.$BUILD.txt > $WORKDIR/$OUTPREFIX.Coverage.$BUILD.txt.gz
fi

## Check for bad things and clean up
grep ERROR $WORKDIR/$OUTPREFIX.pipeline.log > $WORKDIR/$OUTPREFIX.pipeline.ERROR.log
if [ -z $(cat $WORKDIR/$OUTPREFIX.pipeline.ERROR.log) ]; then
	rm $WORKDIR/$OUTPREFIX.pipeline.ERROR.log $OUTPREFIX.marked.sort.bwa.$BUILD.bam $OUTPREFIX.marked.sort.bwa.$BUILD.bai
	rm -r $tmpDir
else 
	echo "Some bad things went down while this script was running please see $OUTPREFIX.pipeline.ERROR.log and prepare for disappointment."
fi
