#!/bin/bash
#SBATCH --job-name=InterProScan_multi --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=[netID]@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs InterPro on multiple protein fasta files #
#-----------------------------------------------------------------------------#

module load java/11.0.10+9

#- Set variables ----------------------------------------------------------------#

InterPro=/home/sm3679/bin/my_interproscan/interproscan-5.50-84.0/interproscan.sh

input_dir=/home/sm3679/albopictus_remap/albopictus_genome/broken_protein_fasta

output_dir=/home/sm3679/albopictus_remap/interpro_dir


#- RUN InterPro ----------------------------------------------------------------#

files=(${input_dir}/*.fa)
for file in ${files[@]}
do
${InterPro} -i ${file} -d ${output_dir} -f TSV -goterms
done

#- FIN -----------------------------------------------------------------------#
