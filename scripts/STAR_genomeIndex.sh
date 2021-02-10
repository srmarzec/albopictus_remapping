#!/bin/bash
#SBATCH --job-name=STAR_genomeIndex --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=[netID]@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=100G

#-----------------------------------------------------------------------------#
# This script makes index of the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

refgen_dir=[path to ref genome and gff3 file]
refgen_index=[path to directory of new indexed genome]


#- RUN STAR----------------------------------------------------------------#

STAR --runMode genomeGenerate \
        --genomeDir ${refgen_index} \
        --genomeFastaFiles ${refgen_dir}/aedes_albopictus_AalbF3.fa \
        --sjdbGTFfile ${refgen_dir}/aedes_albopictus_AalbF3.gff3

#- FIN -----------------------------------------------------------------------#
