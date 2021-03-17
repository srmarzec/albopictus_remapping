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



