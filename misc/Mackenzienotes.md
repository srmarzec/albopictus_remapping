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
