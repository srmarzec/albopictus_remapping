# Sarah's Notes
Things Sarah has tried, wants to keep track of

## Renaming files after pulling them from SRA
I knew that when I downloaded the sra files (and converted them into fastq.gz) that the names would be th accession numbers, which isn't as helpful as having the sample names. I thus also downloaded the metadata to my local computer and deleted everythign but two columns (the accession # and sample names). I looked into a bunch of way to rename files based on a csv and I'm sure somehow this method would have worked but I couldn't make it work in a day and I could easily write out a command in excel to print lines pasted with mv and oldname, newname. Maybe figure out how to rename based on a csv later...

## Downloading the right version of trimmomatic
We will be using version 0.39
```
# enter directory you want to store trimmomatic in
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

# you can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version
```
