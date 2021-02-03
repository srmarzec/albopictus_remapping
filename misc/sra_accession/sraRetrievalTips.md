# Accessing data from the SRA
Sending information [directly to a google bucket](https://www.ncbi.nlm.nih.gov/sra/docs/data-delivery/) requires google bucket admin permission which I don't have. Instead I [downlowded the SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to my home directory on the GU hpc (Note: I used the ubuntu version of the toolkit and it works fine).

Within the sra-toolkit I have a particular problem where I must [append the PATH](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit#3-for-convenience-and-to-show-you-where-the-binaries-are-append-the-path-to-the-binaries-to-your-path-environment-variable) each time I start work on the hpc again. Not sure how much this impacts anyone who completes steps in one go.

Tried fasterq-dump to bring a file down (converts sra file into fastq file for you). Two things to note: 1) Took ~13 minutes for transfer of two 11G files and 2) files are named with accession number and no identifying data. Additionally, fasterq-dump can only do a single file at a time... Outputs files into working directory
```
fasterq-dump SRR1663689
```

Tried prefetch which brings down the sra file and deposits file into the sra folder of the chosen local directory (set when configuring the sratoolkit). Bringing the sra file takes just over a minute. You can also use prfetch with an accession list. After you have the sra file you can use fasterq-dump to get the \*.fastq files (took ~6 minutes to convert 2 7.8G files)
```
prefetch SRR1663685
fasterq-dump sra_directory/sra/SRR1663685.sra
```

Considering you can only use fasterq-dump one file at a time, I'm not really sure what the most effective way to bring fastq files down from the sra site...

### Most efficient way
Apparently you can run fasterq-dump with a wildcard if you have prefetched the sra files to a local directory. Need to gzip as well. Pefetching each file can take quite a bit of time. I ran this interactively and found that for 17 sra files, it took me 2 hours. However, the fasterq-dump comand takes about 5 minutes per sra file since they are prefetched in a local directory. 
```
prefetch --option-file SraAccList.txt
fasterq-dump sra_directory/sra/*.sra
gzip *.fastq
```
