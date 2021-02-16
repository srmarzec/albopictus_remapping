#!/bin/bash
#SBATCH --job-name=fastqcembryo --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mlp134@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script gives quality control of embryo fastq files #
#-----------------------------------------------------------------------------#


#- RUN fastqc ----------------------------------------------------------------#

/home/mlp134/FastQC/fastqc -o /home/mlp134/RNAseqproject/Fastqoutput /home/mlp134/RNAseqproject/trim_embryo_output/*PE.fastq.gz

#- FIN -----------------------------------------------------------------------#
