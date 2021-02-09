#!/bin/bash
#SBATCH --job-name=trim --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=[netID]@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs trimmomatic to clean up raw reads #
#-----------------------------------------------------------------------------#


#- Set variables ----------------------------------------------------------------#

raw_dir=[location of raw reads]

trim_dir=[output directory for trimmed reads]

trim=[path to trimmomatic]/Trimmomatic-0.39/trimmomatic-0.39.jar

adapter=[path to trimmomatic]/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10

#- RUN Trimmomatic----------------------------------------------------------------#

files=(${raw_dir}/*_1.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _1.fastq.gz`
java -Xmx2G -jar ${trim} PE \
          ${raw_dir}/${base}_1.fastq.gz \
          ${raw_dir}/${base}_2.fastq.gz \
          ${trim_dir}/${base}_1_PE.fastq.gz ${trim_dir}/${base}_1_SE.fastq.gz \
          ${trim_dir}/${base}_2_PE.fastq.gz ${trim_dir}/${base}_2_SE.fastq.gz \
          ILLUMINACLIP:${adapter} \
          LEADING:6 \
          TRAILING:6 \
          SLIDINGWINDOW:4:15 \
          MINLEN:50
done

#- FIN -----------------------------------------------------------------------#
