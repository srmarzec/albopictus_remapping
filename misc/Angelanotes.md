## Jan-Feb 10
### Downloading trimmomatic into home directory
Enter directory you want to store trimmomatic in: /home/zz220
```
$ curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
$ unzip Trimmomatic-0.39.zip
$ rm Trimmomatic-0.39.zip
```
### You can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version

### Downloading the SRA toolkit
Chose Ubuntu: $ wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```
$ tar -vxzf sratoolkit.tar.gz
$ export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin
```
test: 
```
$ which fastq-dump
~/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump
```
test worked

completed https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration\
test
`$ fastq-dump --stdout SRR390728 | head -n 8` 
```
@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
fastq-dump was killed (signal 13 SIGPIPE)
```

Notes: 
1. you must reset path everytime you use the SRA toolkit: `$ export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin`
1. You must cd to the directory you set as your repository directory during the config; in this case the directory is: albopictus_remap

### Downloading SRA adult files
Now that SRA toolkit is downloaded, we want to use it do download the (Sequence Read Archive (SRA) files for the *adult* Ae albo the list of accession numbers can be found here: albopictus_remapping/misc/sra_accession/SRR_Acc_List_Adult.txt on github. Also listed here: 
```
SRR1663689
SRR1663685
SRR1663687
SRR1663697
SRR1663700
SRR1663703
SRR1663707
SRR1663709
SRR1663754
SRR1663769
SRR1663843
SRR1663911
SRR1663913
SRR1663916
SRR1664190
SRR1664192
```
Following the method Sarah outlined as the most efficient in master notes: working within Remapping directory

1. put all the numbers in a nano text editor file: $ nano SraAccList.txt
1. $ srun --pty bash
1. prefetch --option-file SraAccList.txt *note: need to be in the same folder as the SraAccList file
1. [zz220@bananas-controller remapping_RNAseq]$ fasterq-dump sra/*.sra
1. gzip *.fastq

*the txt file was really helpful because prefetching the files took a few hours*
Result: 32 fastq.gz files within `/home/zz220/albopictus_remap/adult_rawdata`



Renamed the files:
```
#!/bin/bash
#This is just a bunch of renaming commands to change sra accession numbers into actual relevant info

mv SRR1663685_1.fastq.gz SD_BM_rep1_1.fastq.gz
mv SRR1663687_1.fastq.gz SD_BM_rep2_1.fastq.gz
mv SRR1663689_1.fastq.gz SD_BM_rep3_1.fastq.gz
mv SRR1663697_1.fastq.gz SD_BM_rep4_1.fastq.gz
mv SRR1663700_1.fastq.gz SD_NB_rep1_1.fastq.gz
mv SRR1663703_1.fastq.gz SD_NB_rep2_1.fastq.gz
mv SRR1663707_1.fastq.gz SD_NB_rep3_1.fastq.gz
mv SRR1663709_1.fastq.gz SD_NB_rep4_1.fastq.gz
mv SRR1663754_1.fastq.gz LD_BM_rep1_1.fastq.gz
mv SRR1663769_1.fastq.gz LD_BM_rep2_1.fastq.gz
mv SRR1663843_1.fastq.gz LD_BM_rep3_1.fastq.gz
mv SRR1663911_1.fastq.gz LD_BM_rep4_1.fastq.gz
mv SRR1663913_1.fastq.gz LD_NB_rep1_1.fastq.gz
mv SRR1663916_1.fastq.gz LD_NB_rep2_1.fastq.gz
mv SRR1664190_1.fastq.gz LD_NB_rep3_1.fastq.gz
mv SRR1664192_1.fastq.gz LD_NB_rep4_1.fastq.gz
mv SRR1663685_2.fastq.gz SD_BM_rep1_2.fastq.gz
mv SRR1663687_2.fastq.gz SD_BM_rep2_2.fastq.gz
mv SRR1663689_2.fastq.gz SD_BM_rep3_2.fastq.gz
mv SRR1663697_2.fastq.gz SD_BM_rep4_2.fastq.gz
mv SRR1663700_2.fastq.gz SD_NB_rep1_2.fastq.gz
mv SRR1663703_2.fastq.gz SD_NB_rep2_2.fastq.gz
mv SRR1663707_2.fastq.gz SD_NB_rep3_2.fastq.gz
mv SRR1663709_2.fastq.gz SD_NB_rep4_2.fastq.gz
mv SRR1663754_2.fastq.gz LD_BM_rep1_2.fastq.gz
mv SRR1663769_2.fastq.gz LD_BM_rep2_2.fastq.gz
mv SRR1663843_2.fastq.gz LD_BM_rep3_2.fastq.gz
mv SRR1663911_2.fastq.gz LD_BM_rep4_2.fastq.gz
mv SRR1663913_2.fastq.gz LD_NB_rep1_2.fastq.gz
mv SRR1663916_2.fastq.gz LD_NB_rep2_2.fastq.gz
mv SRR1664190_2.fastq.gz LD_NB_rep3_2.fastq.gz
mv SRR1664192_2.fastq.gz LD_NB_rep4_2.fastq.gz
```

## Feb 13
Reorganized my directory:
- made directory: `zz220@bananas-controller albopictus_remap` This will be the overarching directory with the files related to this remapping project
- within `albopictus_remap` made a directory called `adult_rawdata` Moved my \*.gz i.e. raw adult albo data files to this dir
- within `albopictus_remap` made dir called `scripts`; will contain my trimmomatic script
- within `albopictus_remap` made dir called `trim_ouput`

Edited [trim.sh](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/trim.sh) into [Angela_trim.sh](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_trim.sh); this is the script I will be using to clean the raw SRA adult reads\
- in my directory, it is:`/home/zz220/albopictus_remap/scripts/trim_adult_script.SBATCH`
- files appeared in the trim_out directory, but did not recieve an email confirming job success- unsure if this is an issue.
- result: 32 parired end files, moved to `/home/zz220/albopictus_remap/trim_output/pairedendfiles` and 32 SE files

## Feb 15
### Downloading and running FastQC on bananas controller
- downloaded Win/Linux zip file from [here](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) to local machine. 
- brought up to google cloud with 
```
gcloud compute scp /Users//Users/cottonellezhou/Downloads/fastqc_v0.11.9.zip bananas-controller:.
```
- unzipped fastqc_v0.11.9.zip with: `$ unzip fastqc_v0.11.9.zip` This created a folder called FastQC
- modified permissions with `$ chmod u+x FastQC` Must do this within the FastQC folder.
- modified Sarah's [script for FastQC](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/fastqc.sh) to [fastqc_script.SBATCH](https://github.com/srmarzec/albopictus_remapping/tree/main/scripts)
- after running script for FastQC, result is 25 html files and 25 gz files: 7/32 were lost during this process
- when checking the FastQC job report file, fastqc.96531.out, some areas said failed to process file -----. For example:
```
Failed to process file LD_BM_rep2_2_PE.fastq.gz
uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Ran out of data in the middle of a fastq entry.  Your file is probably truncated
        at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:179)
        at uk.ac.babraham.FastQC.Sequence.FastQFile.next(FastQFile.java:125)
        at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:77)
        at java.lang.Thread.run(Thread.java:748)
```

Bring the html files back to local to veiw in web format (run command from local)
```
Cottonelles-MBP:~ cottonellezhou$ gcloud compute scp bananas-controller:/home/zz220/albopictus_remap/fastqc_out/*_fastqc.html .
```
## Feb 17
Will redownload the raw data for the files that did not succesfully go through FastQC:
```
mv SRR1663769_2.fastq.gz LD_BM_rep2_2.fastq.gz
mv SRR1663916_1.fastq.gz LD_NB_rep2_1.fastq.gz
mv SRR1663916_2.fastq.gz LD_NB_rep2_2.fastq.gz
mv SRR1663685_1.fastq.gz SD_BM_rep1_1.fastq.gz
mv SRR1663685_2.fastq.gz SD_BM_rep1_2.fastq.gz
mv SRR1663697_1.fastq.gz SD_BM_rep4_1.fastq.gz
mv SRR1663697_2.fastq.gz SD_BM_rep4_2.fastq.gz
 ```
Repeated the trimmomatic and FastQC workflow for the above files:
You can see the difference between the size of raw files that were truncated (thus did not go through FastQC succesfully) and the redownloaded version of these files. For example:\
```
new file sizes:
-rw-r--r--. 1 zz220 users 2.3G Feb 17 19:18 LD_BM_rep2_1.fastq.gz
-rw-r--r--. 1 zz220 users 2.3G Feb 17 19:18 LD_BM_rep2_2.fastq.gz
-rw-r--r--. 1 zz220 users 2.5G Feb 17 19:25 LD_NB_rep2_1.fastq.gz
-rw-r--r--. 1 zz220 users 2.4G Feb 17 19:25 LD_NB_rep2_2.fastq.gz
-rw-r--r--. 1 zz220 users 2.1G Feb 17 19:00 SD_BM_rep1_1.fastq.gz
-rw-r--r--. 1 zz220 users 2.0G Feb 17 19:00 SD_BM_rep1_2.fastq.gz
-rw-r--r--. 1 zz220 users 3.7G Feb 17 19:11 SD_BM_rep4_1.fastq.gz
-rw-r--r--. 1 zz220 users 3.6G Feb 17 19:11 SD_BM_rep4_2.fastq.gz

vs

old file sizes:
-rw-------. 1 zz220 users 716M Feb  8 20:49 LD_BM_rep2_2.fastq.gz
-rw-------. 1 zz220 users 2.1G Feb  9 00:02 LD_NB_rep2_1.fastq.gz
-rw-r--r--. 1 zz220 users 2.4G Feb  8 13:08 LD_NB_rep2_2.fastq.gz
```

*note: the trim.out file tells you the number of raw reads and the number that survived in the paired end file. fastqc reports tell you the number of reads that survived in the paired end file.* You can find the data [here](https://docs.google.com/spreadsheets/d/1hIqqMIk8ZVw56BJ8_YN_OnwuzdWBSH7bujwDMTqsCKs/edit#gid=0)

### Mapping
Downloading the reference genome:

```
#From within the genome directory you make: /home/zz220/albopictus_remap/albopictus_genome
#genome fasta file
gsutil cp gs://gu-biology-pi-paa9/aedes_albopictus_AalbF3.fa .
#genome annotation file
gsutil cp gs://gu-biology-pi-paa9/aedes_albopictus_AalbF3.gff3 .
```
result: one fasta file and one annotation file downloaded in genome folder\

Making the index:
Ran the following script from: [scripts/Angela/Angela_STAR_genomeIndex.sh](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_STAR_genomeIndex.sh)

Result: 15 files in /home/zz220/albopictus_remap/albopictus_genome_index

Mapping: ran the following script: /home/zz220/albopictus_remap/scripts/STAR_map.SBATCH\
Can also find it here on gitbhub: [scripts/Angela/Angela_STAR_map.sh](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_STAR_map.sh)

Took around 11 hours to run
Result: 6 files for each of the 16 genomes in `/home/zz220/albopictus_remap/sam_dir` ending in
- Aligned.out.sam
- Log.final.out
- SD_BM_rep3_Log.out
- Log.progress.out
- SJ.out.tab
- STARtmp

less *Log.final.out and typed :n to go to the next file each time (typing :p takes you to the previous file). And of course, typing q lets you escape less.

## Feb 24
### Deleted raw data files to save space\
### Uploading the trimmed data files into the google bucket to save space\
1. created a [trimmed file folder](https://console.cloud.google.com/storage/browser/gu-biology-pi-paa9/zz220/trimmed_files?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22) on google bucket within my netID directory
2. from within my bananas controller, trim_output directory: `$ gcloud auth login`
3. `$ gsutil -m mv *.fastq.gz gs://gu-biology-pi-paa9/zz220/trimmed_files

Result: succesfully uploaded 32 PE files and 32 SE onto the google bucket. Files are no longer in my bananas controller directory
### Converted sam to bam: ran [Angela_sam2bam.sh](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_sam2bam.sh)
Result: removed the Aligned.out.sam files from the sam folder; created a .bam and .bai for each of the 16 genomes

## Feb 27
counted read align rate from the sam files by using less *Log.final.out and typed :n to go to the next file each time (typing :p takes you to the previous file)*
put the Uniquely mapped reads % on shared google sheets

### Installing HTSeq:
- Make a virtual environment called python-ev in home directory (```python3 -m venv python-env```) 
- Activate the virtual environment (```source python-env/bin/activate```) 
- Download HTSeq (and its prerequisites): 

```
pip install numpy
pip install pysam
pip install HTSeq 
``` 

Note: to exit virtual environment: `deactivate`

## March 2
Ran into issues trying to use HTSeq, so reinstalled and updated pip:
``` 
python3 -m venv python-env 
pip uninstall numpy
pip uninstall pysam
pip uninstall HTSeq 
``` 

```
pip install numpy
pip install pysam
pip install HTSeq 
```

Testing for strandedness using HTSeq:

Information from [here](https://chipster.csc.fi/manual/library-type-summary.html) and other options to test for strandedness [here](https://github.com/srmarzec/albopictus_remapping/blob/main/misc/SarahNotes.md). I ran can HTSeq on a single sample using either yes `-s yes` or reverse `-s reverse` for the stranded flag and then checked the output files using `last file_name` to find the number of reads that could not be assigned to any gene:

```
htseq-count -f bam -a 20 -r pos -s no -t gene -i gene /home/zz220/albopictus_remap/bam_dir/LD_NB_rep2_Aligned.out.bam /home/zz220/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.gff3 > /home/zz220/albopictus_remap/counts_dir/LD_NB_rep2_htseq_gff_gene_no

htseq-count -f bam -a 20 -r pos -s yes -t gene -i gene /home/zz220/albopictus_remap/bam_dir/LD_NB_rep2_Aligned.out.bam /home/zz220/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.gff3 > /home/zz220/albopictus_remap/counts_dir/LD_NB_rep2_htseq_gff_gene_yes

htseq-count -f bam -a 20 -r pos -s reverse -t gene -i gene /home/zz220/albopictus_remap/bam_dir/LD_NB_rep2_Aligned.out.bam /home/zz220/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.gff3 > /home/zz220/albopictus_remap/counts_dir/LD_NB_rep2_htseq_gff_gene_reverse
```
Found that using the unstranded setting, resulted in 4X fewer unassigned reads so I then decided to also run the stranded flag with the no setting `-s no`:
| Sample_and_info | # Reads assigned no feature|
| ---: | ---: |
| LD_NB_rep2_htseq_gff_gene_no | __no_feature	2,286,514 | 
| LD_NB_rep2_htseq_gff_gene_reverse | __no_feature	9,994,135 |
| LD_NB_rep2_htseq_gff_gene_yes | __no_feature	9,964,263 |

## March 3rd
Ran [script for STAR](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_STAR_map.sh)
Result: took around 11 hours to run; output was 16 \*htseqCount files, 1 for each sample

## March 8th
### Counting Reads Assigned To Genes on R
Created a folder on my laptop for DE analysis
Downloaded HTSeq files onto my local desktop (will move to the DE analysis data folder) for differential analysis on RStudio: `gcloud compute scp bananas-controller: /home/zz220/albopictus_remap/htseqcounts_dir/*htseqCount .`

Ran [script](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_CountReadsAssignedToGenes.R) and put results in shared google sheet

## March 15th
### DESeq
Ran [Angela_DESeq.R script](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/Angela/Angela_DESeq.R)\
I had 16 data files in total, but we decided to 
Results:
For NB (non-blood fed, 8 files):

![barplot of library sizes](https://user-images.githubusercontent.com/78465068/111240013-c5040980-85d0-11eb-9b40-95842d5cebad.png)
Figure 1. Barplot of library sizes. The bars show the distribution of the count data across the 8 samples. The bars are similar in height and our data looks normal.\


![Logcounts](https://user-images.githubusercontent.com/78465068/111241101-d4845200-85d2-11eb-8cb1-42c9d15ca159.png)
Figure 2. Boxplot of per gene counts for each of the sample groups. The median log counts are similar across samples.

![PCA counts_percentage of the variance can be explained by PC1 (LD vs SD](https://user-images.githubusercontent.com/78465068/111241126-e36b0480-85d2-11eb-9802-4131c0ea632c.png)
Figure 3. PCA plot for the nonblood fed samples. Most of the variance (67%) can be explained by PC1 which is short-day vs long day conditions.

![heatmap copy](https://user-images.githubusercontent.com/78465068/111241135-e6fe8b80-85d2-11eb-83ee-10ba9b35e0b5.png)
Figure 7. Heat map- we want to look out for a row or block that show the same underexpression or over expression pattern (i.e. are either all red or all blue)

![MA plot NB](https://user-images.githubusercontent.com/78465068/111242570-cf74d200-85d5-11eb-9301-cbced2a69dbd.png)
Figure 4. MA plot: an application of a Bland–Altman plot for visual representation of genomic data. The plot visualizes the differences between measurements taken in two samples, by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values. In this case, we are plotting the log fold change against the mean of normalized counts. 

![MA LFC](https://user-images.githubusercontent.com/78465068/111242578-d0a5ff00-85d5-11eb-98a2-ff38fa749034.png)
Figure 4. MA plot after cleaning up the data: excluded data points with low mean counts and a large log fold change.

![plot of variation of data](https://user-images.githubusercontent.com/78465068/111241130-e534c800-85d2-11eb-96fc-b7e0fcb5e9d2.png)
Figure 5. Plot of the variation or dispersion in the data

<img width="689" alt="Volcano plot_NB" src="https://user-images.githubusercontent.com/78465068/111241861-56c14600-85d4-11eb-8fb8-2547afc01519.png">
Figure 6. A volcano plot showing the statistical significance (P value) versus magnitude of change (fold change). Points in red are genes with large fold changes that are also statistically significant. These may be the most biologically significant genes. 

Table of significant differentially expressed genes:
[NB_SDvLD_test.txt](https://github.com/srmarzec/albopictus_remapping/files/6145371/NB_SDvLD_test.txt)

Table with only the gene names (to be used later in KEGG analysis):
[NB_SDvLD_LFCshrink_padj.txt](https://github.com/srmarzec/albopictus_remapping/files/6145372/NB_SDvLD_LFCshrink_padj.txt)

For BM (8 files):
![Volcano plot BM](https://user-images.githubusercontent.com/78465068/111241276-2fb64480-85d3-11eb-903c-5ca323b5580f.png)

[BM_SDvLD_LFCshrink_padj.txt](https://github.com/srmarzec/albopictus_remapping/files/6145374/BM_SDvLD_LFCshrink_padj.txt)

Venn Diagram
![Venn Diagram](https://user-images.githubusercontent.com/78465068/111241469-95a2cc00-85d3-11eb-86d7-e6efbe6239d3.png)

## March 25
