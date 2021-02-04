# Workflow Record
## February 4
Make a directory for this workflow titled "RNAseqproject"
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

