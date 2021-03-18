# GO term annotation for *Aedes albopictus*

One issue with working with non-model organisms (especially one with a new reference genome) is that is is hard to do canned GO enrichment analyses that everyone suggests simply because we don't have readily available GO terms. One way we are fortunate is that the new reference genome we have (AalbF3) actually has the exact same annotation as the previous reference (AalbF2). The difference between these two genome is a change of coordinates.

I was able to get the old protein fasta so that I could use that to search databases and hopefully get GO terms. It turns out the [InterPro](http://www.ebi.ac.uk/interpro/) was truly great at this and there is a command line version in [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/Introduction.html) so I was able to query the protein sequences for *A. albopictus* to find GO terms.

# Obtain protein fasta
Took protein fasta and annotation file from ncbi for AAlbF2 (which is Aedes albopictus annotation release 102; accession number GCF_006496715.1) (upload date stated is 2019-08-06 20:54). We do not have a version of this for our genome since it's new.
```
curl -LO https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/7160/102/GCF_006496715.1_Aalbo_primary.1/GCF_006496715.1_Aalbo_primary.1_protein.faa.gz
gunzip GCF_006496715.1_Aalbo_primary.1_protein.faa.gz
rm GCF_006496715.1_Aalbo_primary.1_protein.faa.gz
```

# Split the protein fasta
The complete protein fasta (AalbF2) is way too large to run in a single go. I first had to split the files.

Split a large fasta file into multiple smaller files based on ">". The number of files made is based on how large to original file is. Each file is named "myseq*.fa" where the * is a number
```
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ../GCF_006496715.1_Aalbo_primary.1_protein.faa
```

# Run InterProScan
After I had the fasta in manageable portions, I ran a [script](https://github.com/srmarzec/albopictus_remapping/blob/main/scripts/InterProScan.sh) made to run over all the files

# Manipulate outputs to get a set or GO annotations for AalbF2 

## Combine the broken up results tsv files into one large file 
Since we ran this on broken up fasta files, there are multiple output tsv files.
```
cat myseq*.fa.tsv >> FULL.tsv
```

## Subset the big tsv by just protein ID and GO-terms
This command prints the first column (protein name) and the GO column
```
awk -v OFS='\t' '{for(i=11;i<=NF;i++){if($i~/^GO:/){a=$i}} print $1,a}' FULL.tsv > FULL_protein_GOterm.tsv
```
`for(...)` loops through all fields, starting with field 11 (i=11; this is because I want to save time instead of going through all columns and I think GO terms will always be after column 11). `if($i~/^GO:/)` checks if the field starts with "GO:". `a=$i` if yes, set variable a to that value. `-v OFS='\t'` sets the "output file separator" as a tab (instead of space which is default)

## Replace protein ID with gene ID
```
awk -F'\t' 'NR==FNR{a[$2]=$1} NR>FNR{$1=a[$1];print}' OFS='\t' "../albopictus_genome/gene_proteinID_list.txt" "FULL_protein_GOterm.tsv" > FULL_gene_GOterm.tsv
```
The first file ("gene_proteinID_list.txt") has column 1 as gene IDs and column 2 as protein IDs

The second file ("FULL_protein_GOterm.tsv") has column 1 as protein IDs and column 2 as GO annotations

`-F'\t'` says the input file is tab delimited. `OFS='\t'` says the output file is tab delimited (default is space delimited)

`NR==FNR{a[$2]=$3}` - NR is an awk internal variable which keeps track of the total number of rows read since the program began. FNR is similar, but keeps track of the number of rows of the current file which have been read. So `NR==FNR` is an awk idiom which means "if this is the first file to be read", and the associated action is `a[$2]=$1` which saves the value of field 1 (gene IDs) in the array a, with the string index being set to the value of field 2 (protein IDs).

`NR>FNR{$1=a[$1];print}'` - similar to the previous, but this time operates only on files other than the first to be read. For each line, we use the value of field 1 (protein IDs) as the index to look up the value in the array, then re-assign field 1 to the array value. Finally, the whole line is printed.

## Clean up the file a little bit 
I have an issue where some of the lines in my file are missing gene IDs. This is an issue from the "replacing protein ID" code above. I checked a few examples by printing out lines in the protein and gene files (e.g. `sed -n '45,53p' FULL_protein_GOterm.tsv` & `sed -n '45,53p' FULL_gene_GOterm.tsv`). Anyways, once I had protein IDs I checked if they were in the list ("gene_proteinID_list.txt") using `grep` and I guess they aren't...

Anyways, this next line subset original file based on unique combinations in column 1 and 2 (unique pairs or gene id and GO-terms) and removes will also remove empty field 1
```
awk -F"\t" '!seen[$1, $2]++' FULL_gene_GOterm.tsv | awk -F"\t" '$1!=""' > FULL_gene_GOterm_unique.tsv
```
Crazy, the above line turns my 18M file into a 1.1M file. This is because InterPro searches a bunch of databases to find GO-terms and so there's a lot of duplicates

`-F"\t"` Must specify that the separator is a tab (otherwise it hates me)

`'!seen[$1, $2]++'` keep lines with unique combinations of fields 1 & 2 (i.e. gene ID and GO-terms) and gets rid of the duplicates. (It says the combo of [$1, $2] has to not been seen before to print)

`'$1!=""'` keeps only lines where the field 1 ($1) isn't (!=) empty ("")

## Rearrange all the GO-terms to be on one line
So we now have many rows with replicate gene IDs that have lists of unique combos of GO-terms. Because InterPro searches different databases, we get different GO-terms, or different combos.

For instance:
```
LOC109400251    GO:0005515
LOC109400251    GO:0005524|GO:0016887
LOC109400064    GO:0003677|GO:0003899|GO:0006351
LOC109400064    GO:0003899|GO:0006351|GO:0032549
LOC109400064    GO:0006351
LOC109400064    GO:0003899
```
So I need to split the GO-terms into separate rows, find unique combos again, and try to print all the unique GO-terms into a single line again nicely. 
```
tr '|' '$' < FULL_gene_GOterm_unique.tsv | tr '\t' ',' | awk '{v1=$0; gsub(/\$[^,]+/,""); gsub (/,[^,]+\$/,",",v1); print $0; if (v1!=$0) print v1}' | tr "," "\t" | awk '!seen[$1, $2]++' | awk -F'\t' -v OFS='\t' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' > FULL_gene_GOterm_single.tsv
```
`tr '|' '$'` Changes the pipe separators to dollar signs cause pipes scare me as separators

`tr '\t' ','` Truthfully only changed the tab to a comma cause the next line of code I found uses comma as a separator and I copied it but don't care enough to change the commas in it cause my file is small enough to run with this extra command and I don't understand the next command fully but it works, so yeah...

`awk '{v1=$0; gsub(/\$[^,]+/,""); gsub (/,[^,]+\$/,",",v1); print $0; if (v1!=$0) print v1}'` The "next command" that works. This searches for multiple variables in a field that are separated by a $ and then should print new rows with everything else in the row the same, except now each variable that was originally in the multiple variable field will have it's own row. Exactly what we need to separate these GO terms onto new rows.

`tr "," "\t"` Change the commas back to tabs cause I can't be bother to change the above awk command...

`awk '!seen[$1, $2]++'` Find unique rows

`awk -F'\t' -v OFS='\t' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}'` This awk looks through all the rows and finds all the rows with the same field 1, then it prints all the field 2 as a new field within one row with the matching field 1 

The result file of this will have one copy of each unique GO term associated with a protein ID. The format is one protein ID per row with GO terms as fields in the rows. For example:
```
LOC109417152            GO:0003824
LOC109417153            GO:0005515      GO:0006886      GO:0031267
LOC115266828            GO:0003714      GO:0006351
LOC109418934            GO:0005249      GO:0008076      GO:0003677      GO:0032508      GO:0016791
LOC109404534            GO:0000786      GO:0003677      GO:0046982
LOC109418935            GO:0003756      GO:0003854      GO:0016616
LOC109422451            GO:0046856
LOC109417156            GO:0003677
```
Of course, now we have an uneven number of fields/columns for each row...

## Combine GO terms on single line
First we need to count how many columns are in the file
```
awk '{print NF}' FULL_gene_GOterm_single.tsv | sort -nu | tail -n 1
```
The output of this is 44 columns... Which means there's no way I'm going to individually code the separators in order to get the right format

The GO-terms annotation file needed has a specific format. Each line has a gene ID followed by a tab, then all the GO-terms for the gene in a single column separated by a comma and space
```
tr '\t' ',' < FULL_gene_GOterm_single.tsv | sed -e "s/,,/\t/g" -e "s/,/, /g" > FULL_GOannotations.txt
```
`tr '\t' ','` Replace all the tabs with commas

`sed -e "s/,,/\t/g"` Replace the double column after the gene ID with a tab and `-e "s/,/, /g"` changes all the commas to commas plus a space. Couldn't get that to work with tr....

Now we should have a GO annotations file!!! Example:
```
LOC109417152    GO:0003824
LOC109417153    GO:0005515, GO:0006886, GO:0031267
LOC115266828    GO:0003714, GO:0006351
LOC109418934    GO:0005249, GO:0008076, GO:0003677, GO:0032508, GO:0016791
LOC109404534    GO:0000786, GO:0003677, GO:0046982
LOC109418935    GO:0003756, GO:0003854, GO:0016616
LOC109422451    GO:0046856
LOC109417156    GO:0003677
LOC109399785    GO:0003824, GO:0050660
LOC109418936    GO:0016787
```


