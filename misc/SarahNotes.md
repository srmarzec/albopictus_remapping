# Sarah's Notes
Things Sarah has tried, wants to keep track of. A guide and the reasoning behind the main workflow.

## SRA Accession
I had to continually set the path each time I stearted to use the sra-toolkit on the hpc: `export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin`
Once I had prefetched the files using my accession list, I was able convert them to fastq files using fasterq-dump. This is much faster then using fasterq-dump to access files directly as prefetch has copy to use in machine now.
 ```
 fasterq-dump /home/sm3679/sra_directory/sra/*.sra -O /home/sm3679/albopictus_remap/albopictus_larva_rawData
 ```
Note that fasterq-dump does not work with lists (another reason to do prefetch since that will do lists) and additionally it will not output a zipped file so you will need to zip the files later. Also remember to remove the sra files you have in the sra-directory.
[Notes I made before](https://github.com/srmarzec/albopictus_remapping/blob/main/misc/sra_accession/sraRetrievalTips.md)

### Renaming files after pulling them from SRA
I knew that when I downloaded the sra files (and converted them into fastq.gz) that the names would be the accession numbers, which isn't as helpful as having the sample names. I thus also downloaded the metadata to my local computer and deleted everythign but two columns (the accession # and sample names). I looked into a bunch of way to rename files based on a csv and I'm sure somehow this method would have worked but I couldn't make it work in a day and I could easily write out a command in excel to print lines pasted with mv and oldname, newname. Maybe figure out how to rename based on a csv later...

## Cleaning Data
In the past, I've used trimmomatic for trim low quality and adapters and have mapped remaining reads to the reference (assuming leftover reads are from extraneous sequences). This has generally worked in the past and I think will continue to work in the future. However, in PI's past work, steps have been taken for pre-processing of RNAseq data to reduce the amount of rRNA and vector contaminants. Generally preprocessing encompasses removing all of: Unwanted sequence (Ex. polyA tails in RNAseq data), artificially added onto sequence of primary interest (vectors, adapters, primers), join short overlapping paired-end reads, low quality bases, originate from PCR duplication. [This](https://ucdavis-bioinformatics-training.github.io/2019_March_UCSF_mRNAseq_Workshop/data_reduction/preproc_htstream.html) seems to have a nice summary of why to preprocess and also references using [HTStream](https://s4hts.github.io/HTStream/) which could be looked into more, but seems promising. Back to PI, PI's old student used ssaha to identify contaminants and rRNA (?). However, newer than that is using [sortmeRNA](https://bioinfo.lifl.fr/sortmerna/sortmerna.php) which does remove specified reads. That said, I'm not entirely sure how to create the database for this yet. I pulled all the rRNA reads for albopitus from [SILVA](https://www.arb-silva.de/) but I'm not sure if you index these fasta files against the package's existing database or somehow make your own... Additionally I would also want to remove vector contaminants probablly based on [UniVec](ftp://ftp.ncbi.nih.gov/pub/UniVec/). 

All of that said, there may be no reason to really clean the data beforehand. Trimmomatic will remove adapters and low quality reads which is important. After that, mapping to a good reference will mean any unmapped reads are extraneous. But PI brings up an important point: if for some reason we have low read-alignment rates, is this due to contaminants or is it due to bad alignment (based on ref or parameters)? (In this case low read alignment is around 50%, and high read alignment is 70-80%). We've decided to go ahead with trimming a subset of reads and apping them to check read-alignment rates.

### Downloading the right version of trimmomatic
We will be using version 0.39
```
# enter directory you want to store trimmomatic in
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

# you can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version
```
Using trimmomatic: LEADING/TRAILING 6 is a light touch. We went with MINLEN 50 which is a general default. PI has used higher MINLEN in the past (75) for albopictus but lots of data removed. I've used lower MINLEN in the past but albopictus has lots of repeat regions.

UPDATE: We decided that instead of the LEADING/TRAILING flags (which when removed,did little) that we would crop the first 15 bases off using HEADCROP. This is because the "per sequence base content" failed during fastqc and this is common in RNAseq analysis due to 'random' hexamer priming during library prep. Since this was common in all my datasets, we decided to make this a lateral move for all three (adult, embryo, larva) datasets.

### Downloading/Using FastQC
Download Win/Linux zip file from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to local machine. 
Bring up to google cloud with 
```
gcloud compute scp /Users/sarah/Downloads/fastqc_v0.11.9.zip bananas-controller:.
```
Unzip file and change permissions of fastqc within the FastQC folder with `chmod u+x fastqc` in order to run interactively. I would likely run this as a script (even though it takes such a short time) but on test files I run it interactively:
```
/home/sm3679/FastQC/fastqc -o /home/sm3679/albopictus_remap/test /home/sm3679/albopictus_remap/test/*PE.fastq.gz
```
Bring the html files back to local to veiw in web format (run command from local)
```
gcloud compute scp bananas-controller:/home/sm3679/albopictus_remap/test/*_fastqc.html .
```
Notes: Generally my fastqc seem to fail in ceratin areas. I've read that often with RNAseq data the 'per base sequence content' fails [link](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html) (this we fixed by indiscriminantly cropping off the first 15 bases, see update for trimmomatic above). I've also read that 'sequence duplication levels' often fails in RNAseq but this could be due to either technical or biological duplicates [link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html); however, deduplication steps are not common in RNAseq analysis [link](https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/) so I figure we acknowledge this but continue on normally. Lastly, 'per sequence GC content' fails but this could be relatively normal; not sure what it typically looks like for albopictus.

Update: When we ran all samples through trimmomatic and veiwed fastqc, it turns out we get a few overrepresented sequences, but these all fall in albopictus (potentially are biologically relevant because they're RNA data and not technical artefacts) and so we don't worry. Additionally, we found that the adult data failed at times in the 'per sequence tile score' and this could be because there were problems with the flow cell (most common cause of this failure is air bubble in the flow cell). That said, the data has passed trim and seems good quality so we'll keep it. This seems most common step. Although if a good reference is not available then further steps could be taken to sort based on tile quality (see [here](https://www.biostars.org/p/228762/)). 

___

After we know we have the trimmed files we want (all the fastqc looks good for each dataset), we moved raw reads up to the google bucket for saving space.
```
#from within my folder containing raw reads
gsutil -m mv *.fastq.gz gs://gu-biology-pi-paa9/sm3679/albopictus_remap/larva_rawData/
```

## Mapping
Note that when generating the genome index using STAR, you need a bunch of memory (I put a 100G), or it will fail. We used a new reference genome generated by PI's group so it is not publically accessible yet. PI put a copy up for us to access on lab google bucket
```
#From within the genome directory you make
#genome fasta file
gsutil cp gs://gu-biology-pi-paa9/aedes_albopictus_AalbF3.fa .
#genome annotation file
gsutil cp gs://gu-biology-pi-paa9/aedes_albopictus_AalbF3.gff3 .
```
I put these in one directory separate from what will be the genome index directory for STAR. This is mostly because I didn't want any issues as STAR commands call for the directory name and not specific files (I didn't want file confusion).

___
When finding the read mapping rates after using STAR, I went through the files one after another using `less`. I'm sure you can cat the files together but the problem is that they don't include the file name in the final output file (although for the sake of saying it, my files and how I record are both in the same order, alphabetical, and thus I could have tried this). Instead I used `less *Log.final.out` and typed `:n` to go to the next file each time (typing `:p` takes you to the previous file). And of course, typing `q` lets you escape `less`.
___
For some reason, using the option to get sorted-by-coordinate bam files as an output from star resulted in quite variable bam files sizes (one of them was even listed as 1.9M although the mapping rate was close to 80%). I wasn't sure why this was from the log files or looking online. I decided to try leaving the default output as a sam file, and then sorting and indexing later using samtools. This is an extra step (although not completely troublesome as I could spare the space for temporary sam storage and combined the sorting and indexing commands in one script), but it resulted in much more consistent bam file sizes. I'm still not sure of the initial issue but I'm going to go with the sam default output, which I then converted to bam and sorted/indexed with samtools.

Table: File sizes with different output methods. All sizes listed are for the bam files made either initially by STAR or which were converted to bam using samtools.
|	bamOutput_sortedByCoordinateSTAR	|	samOutput_ConvertedSortedSamtools	|	FileName	|
|	---:	|	---:	|	:---	|
|	2.0G	|	1.7G	|	D_11d_rep1_Aligned.out.bam	|
|	1.9M	|	2.3G	|	D_11d_rep2_Aligned.out.bam	|
|	2.5G	|	2.1G	|	D_11d_rep4_Aligned.out.bam	|
|	3.5G	|	2.9G	|	D_21d_rep1_Aligned.out.bam	|
|	1.6G	|	1.3G	|	D_21d_rep2_Aligned.out.bam	|
|	2.4G	|	2.0G	|	D_21d_rep4_Aligned.out.bam	|
|	1.9G	|	1.6G	|	D_40d_rep2_Aligned.out.bam	|
|	1.7G	|	1.4G	|	D_40d_rep3_Aligned.out.bam	|
|	3.7G	|	3.1G	|	D_40d_rep4_Aligned.out.bam	|
|	2.8G	|	2.3G	|	ND_11d_rep1_Aligned.out.bam	|
|	1.1G	|	872M	|	ND_11d_rep2_Aligned.out.bam	|
|	1.4G	|	1.2G	|	ND_11d_rep4_Aligned.out.bam	|
|	1.3G	|	1.1G	|	ND_21d_rep1_Aligned.out.bam	|
|	3.0G	|	2.5G	|	ND_21d_rep2_Aligned.out.bam	|
|	1.2G	|	1000M	|	ND_21d_rep3_Aligned.out.bam	|
|	2.4G	|	2.0G	|	ND_40d_rep3_Aligned.out.bam	|
|	3.4G	|	2.8G	|	ND_40d_rep4_Aligned.out.bam	|

## Generating count matrix with HTSeq (htseq-count)
Future analysis will be done in DESeq and HTSeq generated count data is an acceptable input. I think I'm using HTSeq mostly because I have seen it used before and we have decided to map to a reference genome (as compared to a transcriptome which seems favored with other packages).

### Installing HTSeq 
So to install HTSeq on the hpc, you will have to make a virtual environment. This is because based on the [GU HPC wiki](https://wiki.uis.georgetown.edu/display/HPCGCP/Python) personalized python installtion need to be done in [virtual environments](https://docs.python.org/3.6/tutorial/venv.html).

Within the home directory (or whever you want the virtual environment to go), we can make a virtual environment called "python-env": `python3 -m venv python-env` \
We can then activate the environment: `source python-env/bin/activate` \
Within the environment, we can download HTSeq (and its prerequisites):
```
pip install numpy
pip install pysam
pip install HTSeq
``` 
To exit the virtual environment: `deactivate`

We should be able to run scripts with HTSeq as long as we activate the virtual environment within the script (see [example](https://implement.pt/2018/09/using-python-virtual-environments-with-slurm/))

### HTSeq parameters
#### Stranded `-s`
So I think one of the biggest issues with using HTSeq is determining strandedness which we have to specify when running ([htseq-count manual](https://htseq.readthedocs.io/en/master/count.html)). I ran a couple of [tests](https://github.com/srmarzec/albopictus_remapping/blob/main/misc/TestingStrandedness.md) to check for my data, and found my data was not strand specific.

#### Type `-t`
I also tested the type input: either exon or gene. (This next part was doen in R)
```
type_gene <- read.table("D_11d_rep1_htseq_gff_gene_no_bamV2", header = F) #ran with htseq -t gene
type_exon <- read.table("D_11d_rep1_htseq_gff_exon_no_bamV2", header = F) #ran with htseq -t exon (which is the default)

merged <- merge(type_gene, type_exon, by = "V1")

head(merged)
```
```
##                       V1    V2.x    V2.y
## 1 __alignment_not_unique       0       0
## 2            __ambiguous  435151  102258
## 3           __no_feature 1637117 1284132
## 4          __not_aligned       0       0
## 5        __too_low_aQual       0       0
## 6           LOC109396977      62      62
```
Note that the above files were generated with a minimum quality flag of 20, `-a 20`, and since no reads were unasigned due to quality, it doesn't matter if I choose this or the default of a min score of 10.
```
#getting rid of anything that wasn't assigned a gene
merged_v2 <- merged[c(-1,-2,-3,-4,-5),]
#counting up everything assigned a gene
colSums(merged_v2[,c(2,3)])
```
```
##    V2.x    V2.y 
## 9527473 9648216 
```
Looks like type exon assigns more reads to genes, so I will be going with type exon instead of gene, `-t exon`.

#### ID attribute `-i`
This is what the counts are sorted under. The default is gene_id if you are using a gtf file, and in fact HTSeq expects you to use a gtf file instead of a gff file. I actually tried converting my gff file to a gtf file using [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread). However, some of the the lines in the gtf file that are "exons" are missing a "gene_id". I wasn't sure why this was. Apparently you can sort by "transcript_id" which by cursory glance seems to be present at many of the exon lines, however I read online that this leads to a lot more 'ambiguous reads' since transcripts could map to several genes. 

But HTSeq runs well if I give it the gff3 file and specify a different ID attribute. I decided to go with "gene" and this works well. So I think converting it to a gtf file isn't necessary.

# Differential Expression analysis with DESeq
[Vingette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

[Helpful post for QC](https://hbctraining.github.io/DGE_workshop_salmon/lessons/03_DGE_QC_analysis.html)
