#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=zz220@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G
#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#
source /home/zz220/python-env/bin/activate
#- Set variables ----------------------------------------------------------------#
bam_dir=/home/zz220/albopictus_remap/bam_dir
count_dir=/home/zz220/albopictus_remap/htseqcounts_dir
htseq=/home/zz220/albopictus_remap/scripts/htseqPhython.SBATCH 
ref=/home/zz220/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.gff3
#- RUN fastqc ----------------------------------------------------------------#
files=(${bam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do`
base=`basename ${file} _Aligned.out.bam`
${htseq} -f bam -r pos -s no -t exon -i gene ${bam_dir}/${base}_Aligned.out.bam ${ref} > ${count_dir}/${base}_htseqCount
done
#- FIN -----------------------------------------------------------------------#
