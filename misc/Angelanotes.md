## Downloading trimmomatic into home directory
### enter directory you want to store trimmomatic in: /home/zz220
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

### you can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version

### downloading the SRA toolkit
Chose Ubuntu: $ wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz\
```
$ tar -vxzf sratoolkit.tar.gz\
$ export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin\
```
test: 
```
$ which fastq-dump\
~/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump\
```
test worked

completed https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration\
test
`$ fastq-dump --stdout SRR390728 | head -n 8`\
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
1. you must reset path everytime you use the SRA toolkit: $ export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin

1. You must cd to the directory you set as your repository directory during the config; in this case the directory is: remapping_RNAseq

Now that SRA toolkit is downloaded, we want to use it do download the (Sequence Read Archive (SRA) files for the *adult* Ae albo the list of accession numbers can be found here: albopictus_remapping/misc/sra_accession/SRR_Acc_List_Adult.txt on github. Also listed here: 

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

Following the method Sarah outlined as the most efficient in master notes: working within Remapping directory

1. put all the numbers in a nano text editor file: $ nano SraAccList.txt

1. $ srun --pty bash

1. prefetch --option-file SraAccList.txt

1. [zz220@bananas-controller remapping_RNAseq]$ fasterq-dump sra/*.sra

1. gzip *.fastq

*the txt file was really helpful because prefetching the files took a few hours*


