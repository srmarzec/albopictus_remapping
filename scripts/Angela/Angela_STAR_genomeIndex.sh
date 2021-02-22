#!/bin/bash
#SBATCH --job-name=STAR_genomeIndex --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=zz220@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=100G

#-----------------------------------------------------------------------------#
# This script makes index of the ref genome using STAR #--------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

refgen_dir=/home/zz220/albopictus_remap/albopictus_genome
refgen_index=/home/zz220/albopictus_remap/albopictus_genome_index


#- RUN STAR----------------------------------------------------------------#

STAR --runMode genomeGenerate \
        --genomeDir ${refgen_index} \
        --genomeFastaFiles ${refgen_dir}/aedes_albopictus_AalbF3.fa \
        --sjdbGTFfile ${refgen_dir}/aedes_albopictus_AalbF3.gff3

#- FIN -----------------------------------------------------------------------#
