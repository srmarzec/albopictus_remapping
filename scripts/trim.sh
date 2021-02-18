#!/bin/bash
#SBATCH --job-name=bedtools_getfasta --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=[Insert user]@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs bedtools getfasta to pull sequence from full assembly   #
#-----------------------------------------------------------------------------#

## module load bedtools




#- RUN bedtools getfasta----------------------------------------------------------------#

java -Xmx2G -jar /share/BIO379/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
          WTC2_1.fq.gz \
          WTC2_2.fq.gz \
          WTC2_1_trPE.fq.gz WTC2_1_trSE.fq.gz \
          WTC2_2_trPE.fq.gz WTC2_2_trSE.fq.gz \
          LEADING:10 \
          TRAILING:10 \
          SLIDINGWINDOW:4:15 \
          ILLUMINACLIP:/share/BIO379/TruSeq3-PE-2.fa:2:30:10 \
          MINLEN:50


#- FIN -----------------------------------------------------------------------#
