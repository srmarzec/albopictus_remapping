# Sarah's Notes
Things Sarah has tried, wants to keep track of

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
Notes: Generally my fastqc seem to fail in ceratin areas. I've read that often with RNAseq data the 'per base sequence content' fails [link](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html). I've also read that 'sequence duplication levels' often fails in RNAseq but this could be due to either technical or biological duplicates [link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html); however, deduplication steps are not common in RNAseq analysis [link](https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/) so I figure we acknowledge this but continue on normally. Lastly, 'per sequence GC content' fails but this could be relatively normal; not sure what it typically looks like for albopictus.

## Mapping
Note that when generating the genome index using STAR, you need a bunch of memory (I put a 100G), or it will fail.

