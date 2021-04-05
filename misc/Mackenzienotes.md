# Workflow Record
Multiple projects and publications have been completed and generated understanding expression difference of diapause and non-diapause at multiple timepoints during the Aedes albopictus lifecycle. However, now a better quality reference assembly exists and old data needs to be reanalyzed based on this reference. We will re-analyze the various diapause transcriptomes. Potentially, once all analyses are done against common reference, we may also look for clusters of genes that change across timepoints in a coordinated way (look into Kostal et al 2017 and Dowle et al 2020).
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
Installed [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) 
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
Successfully fixed issue. Brought file up to bananas controller using this command on local device. 
``` 
gcloud compute scp /Users/mackenzieparsons/Downloads/fastqc_v0.11.9.zip bananas-controller:.
```
Unzip file and change permissions of fastqc within the FastQC folder. Unsure if I interpreted this step correctly? Used following command. Turned file from red to green color.
```
chmod u+x fastqc_v0.11.9.zip 
unzip fastqc_v0.11.9.zip 
```
Running FastQC using the following script, named fastqc_embryoscript.SBATCH in RNAseqproject directory. 
```
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
```
## February 16th 
FastQC job did not run. Looked in job output file- access issue. Cd into FastQC directory and did command within directory.
```chmod u+x fastqc```

Reran script- worked successfully. Logged onto Terminal on local device and ran ```gcloud compute scp bananas-controller:/home/mlp134/RNAseqproject/Fastqoutput/*_fastqc.html .``` Files successfully downloaded.

Notes on FastQC output files (downloaded to local device). 24 output files. Files look good- all only have red for per sequence GC content. Yellow seen for Sequence Length Distribution and Sequence Duplication levels for most files. Yellow seen for Overrepresented sequences on a few files. 

## February 22th: Mapping 
First step: had to reset the authorization from my accounnt to google bucket. https://wiki.uis.georgetown.edu/pages/viewpage.action?pageId=87097411

``` gcloud auth login ``` 

Made a new folder titled mlp134 within paa9 bucket at link: https://console.cloud.google.com/storage/browser/gu-biology-pi-paa9;tab=objects?forceOnBucketsSortingFiltering=false&authuser=0&project=gcp-gu-hpc-medusa&prefix=&forceOnObjectsSortingFiltering=false. 

Moved raw data into folder within my google bucket folder 

``` gsutil -m mv *.fastq.gz gs://gu-biology-pi-paa9/mlp134/embryo_rawdata ``` 

Make directory within RNAseqproject titled mapping. (/home/mlp134/RNAseqproject/mapping)
Bring down genome fasta file. 

``` gsutil cp gs://gu-biology-pi-paa9/aedes_albopictus_AalbF3.fa . ``` 

Bring down genome annotation file. 

``` gsutil cp gs://gu-biology-pi-paa9/aedes_albopictus_AalbF3.gff3 . ``` 

run STAR index script: in RNAseqproject folder titled STAR_genomeIndex.SBATCH  

``` 
#!/bin/bash
#SBATCH --job-name=STAR_genomeIndex --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mlp134@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=100G

#-----------------------------------------------------------------------------#
# This script makes index of the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

refgen_dir=/home/mlp134/RNAseqproject/mapping
refgen_index=/home/mlp134/RNAseqproject/indexed_genome


#- RUN STAR----------------------------------------------------------------#

STAR --runMode genomeGenerate \
        --genomeDir ${refgen_index} \
        --genomeFastaFiles ${refgen_dir}/aedes_albopictus_AalbF3.fa \
        --sjdbGTFfile ${refgen_dir}/aedes_albopictus_AalbF3.gff3

#- FIN -----------------------------------------------------------------------#
``` 
Took about little over one hour to run. Next ran the STAR mapping script (STARmap.SBATCH). Took about eight hours and ran successfully. Output files placed into sam_dir (/home/mlp134/RNAseqproject/sam_dir)

``` 
#!/bin/bash
#SBATCH --job-name=STAR --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mlp134@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=100G
#-----------------------------------------------------------------------------#
# This script maps reads to the ref genome using STAR #
#-----------------------------------------------------------------------------#
module load star/2.7.1a
#- Set variables ----------------------------------------------------------------#
trim_dir=/home/mlp134/RNAseqproject/trim_embryo_output
out_dir=/home/mlp134/RNAseqproject/sam_dir
refgen_dir=/home/mlp134/RNAseqproject/indexed_genome
#- RUN STAR----------------------------------------------------------------#
files=(${trim_dir}/*_1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _1_PE.fastq.gz`
STAR --genomeDir ${refgen_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${trim_dir}/${base}_1_PE.fastq.gz ${trim_dir}/${base}_2_PE.fastq.gz \
        --outFileNamePrefix ${out_dir}/${base}_  
done
#- FIN -----------------------------------------------------------------------#
``` 
# February 23: Examine mapping output files 
Used ```less *Log.final.out``` command on the output mapped files (within sam_dir folder) and :n to click through each file. Results below. 
```
Sample		Read Align Rate (Uniquely mapped reads?)
SRR458462	81.11%
SRR458463	79.95%
SRR458464	77.56%
SRR458465	78.60%
SRR458466	79.56%
SRR458467	80.60%
SRR458468	80.39%
SRR458469	80.27%
SRR458470	80.15%
SRR458471	79.25%
SRR458472	81.01%
SRR458473	79.75%
```
Now running the script to convert sam to bam files and indexing for future work. (samtobam.SBATCH) 
``` 
#!/bin/bash
#SBATCH --job-name=sam2bam --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mlp134@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script converts sam to bam, sorts/indexes the bam files #
#-----------------------------------------------------------------------------#

module load samtools/1.9

#- Set variables ----------------------------------------------------------------#

sam_dir=/home/mlp134/RNAseqproject/sam_dir
bam_dir=/home/mlp134/RNAseqproject/bam_dir

#- RUN fastqc ----------------------------------------------------------------#

files=(${sam_dir}/*_Aligned.out.sam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.sam`
samtools view -b -S ${sam_dir}/${base}_Aligned.out.sam | samtools sort -o ${bam_dir}/${base}_Aligned.out.bam
rm -f ${sam_dir}/${base}_Aligned.out.sam
samtools index ${bam_dir}/${base}_Aligned.out.bam

done

#- FIN -----------------------------------------------------------------------#
``` 
Installing HITSeq. Make a virtual environment called python-ev (```python3 -m venv python-env```) 
Activate the environment from within the environment. (```source python-env/bin/activate```) 
download HTSeq (and its prerequisites): 
```pip install numpy
pip install pysam
pip install HTSeq 
``` 

Exit virtual environment (```deactivate```) 

## March 3 
Tried running HTSeq count script, failed and received the following error message. 

<img width="928" alt="Screen Shot 2021-03-03 at 9 43 34 PM" src="https://user-images.githubusercontent.com/78369439/110047352-bb9ab780-7d1b-11eb-861c-921df042d6d6.png">

Trying these commands to uninstall the numpy, pysam, HTSeq because issue seems to be with this software. Possibly got interrupted when downloading- similar issues happened to Angela and Sarah as well. 
``` python3 -m venv python-env 
pip uninstall numpy
pip uninstall pysam
pip uninstall HTSeq 
``` 
Did not work. Still getting same error from failed script. Also getting the same error when 
I check for strandedness. 
``` 
source python-env/bin/activate
htseq-count -f bam -a 20 -r pos -s no -t gene -i gene /home/mlp134/RNAseqproject/bam_di
r/DI_135h_rep1_Aligned.out.bam /home/mlp134/RNAseqproject/mapping/aedes_albopictus_AalbF3.gff3 > /home/mlp134/RNAseqpr
oject/counts_dir/DI_135h_rep1_htseq_gff_gene_no 
``` 
<img width="1071" alt="Screen Shot 2021-03-04 at 7 25 15 PM" src="https://user-images.githubusercontent.com/78369439/110049209-69f42c00-7d1f-11eb-884f-e8c029fa490c.png">

## March 4 
Repeated process and eventually worked by updating and rerunning the script. 
 
Reran the test command.
``` 
htseq-count -f bam -a 20 -r pos -s no -t gene -i gene /home/mlp134/RNAseqproject/bam_dir/DI_135h_rep1_Aligned.out.bam /home/mlp134/RNAseqproject/mapping/aedes_albopictus_AalbF3.gff3 > /home/mlp134/RNAseqproject/counts_dir/DI_135h_rep1_htseq_gff_gene_no - 8832673 alignment pairs processed 

htseq-count -f bam -a 20 -r pos -s yes -t gene -i gene /home/mlp134/RNAseqproject/bam_dir/DI_135h_rep1_Aligned.out.bam /home/mlp134/RNAseqproject/mapping/aedes_albopictus_AalbF3.gff3 > /home/mlp134/RNAseqproject/counts_dir/DI_135h_rep1_htseq_gff_gene_yes

htseq-count -f bam -a 20 -r pos -s reverse -t gene -i gene /home/mlp134/RNAseqproject/bam_dir/DI_135h_rep1_Aligned.out.bam /home/mlp134/RNAseqproject/mapping/aedes_albopictus_AalbF3.gff3 > /home/mlp134/RNAseqproject/counts_dir/DI_135h_rep1_htseq_gff_gene_reverse
```
## March 8 
Testing for strandedness by using tail command on each file to look at unassigned reads. Looks like data is unstranded since 1/4 roughly of the number of reads assigned no feature in the -no file compared to the -reverse and -yes files.

| Sample_and_info | # Reads assigned no feature|
| ---: | ---: |
| DI_135h_rep1_htseq_gff_gene_no | __no_feature    1,351,474 | 
| DI_135h_rep1_htseq_gff_gene_reverse | __no_feature	4,517,024 |
| DI_135h_rep1_htseq_gff_gene_yes | __no_feature 4,541,824 |

Ran htseq.count script (htseq_count.SBATCH)
``` 
#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mlp134@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G
#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#
source [path to python env]/bin/activate
#- Set variables ----------------------------------------------------------------#
bam_dir=/home/mlp134/RNAseqproject/bam_dir
count_dir=/home/mlp134/RNAseqproject/counts_dir
htseq=/home/mlp134/RNAseqproject/python-env/bin/htseq-count
ref=/home/mlp134/RNAseqproject/mapping/aedes_albopictus_AalbF3.gff3
#- RUN fastqc ----------------------------------------------------------------#
files=(${bam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.bam`
${htseq} -f bam -r pos -s no -t exon -i gene ${bam_dir}/${base}_Aligned.out.bam ${ref} > ${count_dir}/${base}_htseqCount
done
#- FIN -----------------------------------------------------------------------#
```
Got twelve output files in the counts_dir directory (/home/mlp134/RNAseqproject/counts_dir). Moved files down to my local device on the Terminal app using ```gcloud compute scp bananas-controller:/home/mlp134/RNAseqproject/counts_dir/*htseqCount /Users/mackenzieparsons/Downloads/HTSeqfiles``` 

Once files were downloaded to my computer, I ran the script (Mackenzie_CountReadsAssignedToGenes.R) to count reads assigned to genes in R. 

<img width="258" alt="Screen Shot 2021-03-09 at 8 53 45 PM" src="https://user-images.githubusercontent.com/78369439/110563594-8c1fec80-8119-11eb-8e90-227eff1a40de.png">

## March 15: DESeq
Ran script Mackenzie_DESeq.R in R. Put results files into the folder Downloads/misc. 

Barplot of library sizes for 72h: 
![Barplot of library sizes (line 47) for 72h](https://user-images.githubusercontent.com/78369439/111237258-1e693a00-85cb-11eb-9e76-80178db15a25.png)

Box plot of per gene counts for 72h: 
![Boxplot of per gene counts for sample groups, 72 hours](https://user-images.githubusercontent.com/78369439/111237352-46589d80-85cb-11eb-9504-9f9571646d42.png)

PCA graph for 72 hours:  
![PCA for 72 hours](https://user-images.githubusercontent.com/78369439/111237438-756f0f00-85cb-11eb-9d19-cf3b684dd5e7.png)

PCA graph with transformed normalized counts for 72 hours: 
![PCA with transformed normalized counts, 72 hr](https://user-images.githubusercontent.com/78369439/111237495-96cffb00-85cb-11eb-8604-b1610873302e.png)

Heatmap Pairwise Correlation Values: 
![Heatmap Pairwise Correlation Values, 72 hr](https://user-images.githubusercontent.com/78369439/111237557-b6672380-85cb-11eb-9053-3edf33f69aa4.png)

PlotDispEsts (dds) for 72 hours 
![plotDispEsts (dds), 72 hr](https://user-images.githubusercontent.com/78369439/111237579-caab2080-85cb-11eb-90e6-f40ea5227732.png)

MA plot for  72 hours 

![MA plot, 72 hr](https://user-images.githubusercontent.com/78369439/111237708-0219cd00-85cc-11eb-8dae-2be5a614f32e.png)

plotMA (res_LFC) for 72 hours 
![plotMA (res_LFC), 72 hr](https://user-images.githubusercontent.com/78369439/111237742-12ca4300-85cc-11eb-8d88-3918b0109758.png)

Volcano plot for 72 hours
![Volcano plot, 72 hr](https://user-images.githubusercontent.com/78369439/111237767-1eb60500-85cc-11eb-9905-95661101626a.png)

Volcano plot for 135 hours 

![Volcano plot, 135h](https://user-images.githubusercontent.com/78369439/111237792-2aa1c700-85cc-11eb-902c-066f8a3a03a0.png)

Venn diagram 

![venn diagram](https://user-images.githubusercontent.com/78369439/111237886-5c1a9280-85cc-11eb-9981-e8d0d2c31bfd.png)

## March 24: KEGG analysis 
Worked on running final lines of updated DESeq script. 

## April 4 
Finished working on final lines of DESeq script. Ran Kegg script for both 72h and 135h data sets. The output file contains many more enriched pathways for 135h than for 72h- unsure the reason for this? Recieved the same error message for running the script both times aroune line 200 but rest of script ran fine.
