## Downloading trimmomatic into home directory
### enter directory you want to store trimmomatic in: /home/zz220
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

### you can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version

### downloading the SRA toolkit
Chose Ubuntu: $ wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
$ tar -vxzf sratoolkit.tar.gz

$ export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin

test: $ which fastq-dump

~/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump

test worked

completed https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration

test: $ fastq-dump --stdout SRR390728 | head -n 8

@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
fastq-dump was killed (signal 13 SIGPIPE)

#### not sure if the last line is ok

Note: you must cd to the directory you set as your repository directory during the config


