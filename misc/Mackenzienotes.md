# Workflow Record
## February 4
Made a directory for this workflow titled "RNAseqproject"
Downloaded Trimmomatic using following commands: 
``` 
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/
Trimmomatic-0.39.zip 
unzip Trimmomatic-0.39.zip 
rm Trimmomatic-0.39.zip
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version
```
Installed SRA Toolkit (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) 
Ran into issues with step 3, Quick Toolkit Configuration: can't seem to find fastq-dump even though it appears under ls -l

Used format ./vdb-config, worked; Confirmed functioning with 
```fastq-dump --stdout SRR390728 | head -n 8```

Assigned embryo data. 
```
	SRR458463
	SRR458464
	SRR458465
	SRR458466
	SRR458467
	SRR458468
	SRR458469
	SRR458470
	SRR458471
	SRR458472
	SRR458473
  ```
Put accession numbers into text file entitled: SRAAccList.txt 
```
prefetch --option-file SraAccList.txt
fasterq-dump sra_directory/sra/*.sra
gzip *.fastq
```
## February 11 
Renaming the downloaded files to have more info in the names 
```
mv SRR458462_1.fastq.gz DI_72h_rep1_1.fastq.gz
mv SRR458463_1.fastq.gz DI_72h_rep2_1.fastq.gz
mv SRR458464_1.fastq.gz DI_72h_rep3_1.fastq.gz
mv SRR458465_1.fastq.gz DI_135h_rep1_1.fastq.gz
mv SRR458466_1.fastq.gz DI_135h_rep2_1.fastq.gz
mv SRR458467_1.fastq.gz DI_135h_rep3_1.fastq.gz
mv SRR458468_1.fastq.gz NDI_72h_rep1_1.fastq.gz
mv SRR458469_1.fastq.gz NDI_72h_rep2_1.fastq.gz
mv SRR458470_1.fastq.gz NDI_72h_rep3_1.fastq.gz
mv SRR458471_1.fastq.gz NDI_135h_rep1_1.fastq.gz
mv SRR458472_1.fastq.gz NDI_135h_rep2_1.fastq.gz
mv SRR458473_1.fastq.gz NDI_135h_rep3_1.fastq.gz
mv SRR458462_2.fastq.gz DI_72h_rep1_2.fastq.gz
mv SRR458463_2.fastq.gz DI_72h_rep2_2.fastq.gz
mv SRR458464_2.fastq.gz DI_72h_rep3_2.fastq.gz
mv SRR458465_2.fastq.gz DI_135h_rep1_2.fastq.gz
mv SRR458466_2.fastq.gz DI_135h_rep2_2.fastq.gz
mv SRR458467_2.fastq.gz DI_135h_rep3_2.fastq.gz
mv SRR458468_2.fastq.gz NDI_72h_rep1_2.fastq.gz
mv SRR458469_2.fastq.gz NDI_72h_rep2_2.fastq.gz
mv SRR458470_2.fastq.gz NDI_72h_rep3_2.fastq.gz
mv SRR458471_2.fastq.gz NDI_135h_rep1_2.fastq.gz
mv SRR458472_2.fastq.gz NDI_135h_rep2_2.fastq.gz
mv SRR458473_2.fastq.gz NDI_135h_rep3_2.fastq.gz
```
## February 14th 
Carrying out Trimmomatic read cleaning script. Ran the following script job. Output files all there after about four hours. Did not receive email confirming end job- unsure if this is an issue. 
```
#!/bin/bash
#SBATCH --job-name=trim_embryo --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mlp134@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G
#-----------------------------------------------------------------------------#
# This script runs trimmomatic to clean up raw reads for the embryo sequences #
#-----------------------------------------------------------------------------#
#- Set variables ----------------------------------------------------------------#
raw_dir=/home/mlp134/RNAseqproject/raw_embryo_data
trim_dir=/home/mlp134/RNAseqproject/trim_embryo_output
trim=/home/mlp134/RNAseqproject/Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=/home/mlp134/RNAseqproject/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10
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
          HEADCROP:15 \
          SLIDINGWINDOW:4:15 \
          MINLEN:50
done
#- FIN -----------------------------------------------------------------------#
``` 
## February 15th 
Download Win/Linux zip file (https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) to local machine. Had issue where local device no longer had the Google SDK software. Downloaded using these commands. 
```
curl https://sdk.cloud.google.com | bash
exec -l $SHELL
gcloud init
``` 
Successfully fixed issue. Brought file up to bananas controller. 
``` 
gcloud compute scp /Users/mackenzieparsons/Downloads/fastqc_v0.11.9.zip bananas-controller:.
```
