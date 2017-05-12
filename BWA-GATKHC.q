#!/bin/bash

# Example usage:
# <variables> sbatch BWA-GATKHC.q

#SBATCH -J BWA-GATKHC
#SBATCH -o /fast/users/a1211880/slurmOUT/slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3-00:00 # change this to 3 days for real set
#SBATCH --mem=50GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au

# load modules
module load BWA/0.7.15-foss-2017a
module load Java/1.8.0_121
module load HTSlib/1.3.1-GCC-5.3.0-binutils-2.25
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
### module load GATK 3.7
### module load picard/2.6.0 or higher

# run the executable
./Illumina-Phred33-PE-BWA-Picard-GATKv3.x.HPC.sh -p $OUTPREFIX -s $SEQPATH -o $WORKDIR
