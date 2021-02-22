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
Chose Ubuntu: $ wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz\
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
