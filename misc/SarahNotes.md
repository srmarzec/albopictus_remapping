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

_________________

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

### Using HTSeq
So I think one of the biggest issues with using HTSeq is determining strandedness which we have to specify when running ([htseq-count manual](https://htseq.readthedocs.io/en/master/count.html)). Apparently this can be confusing because this is done during library prep, and it's possible that we won't have this info as we will be using data accessed from the SRA. That said, there are ways to determine this after. [This](https://github.com/mcadamme/Culex_RNAseq_Chemosensory/blob/master/Upstream_processing/strandedness_and_htseq.md) gives an example of why they chose the 'stranded=reverse' flag based on the [Srinivasan et al. 2020](https://academic.oup.com/bfg/article/19/5-6/339/5837822?login=true). I'm going to check the SRA database for what information we have for the datasets (whether they were based on cDNA which, based on Srinivasan et al., gives the strandedness), otherwise I may run two files with either strandedness to cinfirm what my data looks like. If I have to do the test (may be good to do anyways), I will need to communicate these steps with Angela and Mackenzie.
