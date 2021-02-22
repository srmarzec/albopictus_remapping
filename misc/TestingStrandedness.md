# Figuring out strandedness on data for using HTSeq
So I think one of the biggest issues with using HTSeq is determining strandedness which we have to specify when running ([htseq-count manual](https://htseq.readthedocs.io/en/master/count.html)). Apparently this can be confusing because this is done during library prep, and it's possible that we won't have this info as we will be using data accessed from the SRA. That said, there are ways to determine this after. [This](https://github.com/mcadamme/Culex_RNAseq_Chemosensory/blob/master/Upstream_processing/strandedness_and_htseq.md) gives an example of why they chose the 'stranded=reverse' flag based on the [Srinivasan et al. 2020](https://academic.oup.com/bfg/article/19/5-6/339/5837822?login=true).

## Using HTSeq to test strandedness

Information from [here](https://chipster.csc.fi/manual/library-type-summary.html) (which also privides a nice stranded RNA graphic) basically suggests we can check the output files from HTSeq to see which flag we should use. I decided to run htseq-count on a single sample but using either yes `-s yes` or reverse `-s reverse` for the stranded flag. I was quite surprised when in the end, both had almost equivalent amount of reads that could not be assigned to any gene. I then decided to also run the stranded flag with the no setting `-s no`. This showed that I had much fewer reads that could not be assigned.

Output for each trial. These were run with the gff3 files looking for type gene and id attribute as gene. An example command for the code is listed below the table.
| Sample_and_info | # Reads assigned no feature|
| ---: | ---: |
| D_11d_rep1_htseq_gff_gene_no | __no_feature	1637117 | 
| D_11d_rep1_htseq_gff_gene_reverse | __no_feature	6488707 |
| D_11d_rep1_htseq_gff_gene_yes | __no_feature	6513039 |
```
htseq-count -f bam -a 20 -r pos -s no -t gene -i gene /home/sm3679/albopictus_remap/bam_dir/D_11d_rep1_Aligned.sortedByCoord.out.bam /home/sm3679/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.gff3 > /home/sm3679/albopictus_remap/counts_dir/D_11d_rep1_htseq_gff_gene_no
```
From the outputs, it seemed my data was not from a strand-specific assay (use `-s no` for stranded flag).

## Using RSeQC to test strandedness
I decided to double-check the strandedness of the data using RSeQC which can infer strandedness for you. I figured this could be a good sanity check. One slight hurdle is that RSeQC needs a bed file as an input (instead of a gff file).

### Convert gff to bed
I used bedops to convert my gff file to bed.
[Download bedops](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html)

Upload and unzip. As a note, unzipping this makes a directory called 'bin' which was unfortunate cause I had already put all my other binaries in a folder called bin... So then I had to go a pull all those other binaries out... It is probably best to move the zipped file into the bin directory and unzip in there where it will make its own bin directory
```
gcloud compute scp /Users/sarah/Downloads/bedops_linux_x86_64-v2.4.39.tar.bz2 bananas-controller:.
tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2
```
In order to run this you must set the path to find it (I did this from $HOME) `export PATH=$PATH:$PWD/bin/bedops/bin`

Ran for the gff file
```
/home/sm3679/bin/bedops/bin/gff2bed < /home/sm3679/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.gff3 > /home/sm3679/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.bed
```

### Run RSeQC
I looked up information about RSeQC on the [official website](http://rseqc.sourceforge.net/#infer-experiment-py) and from this [page](https://chipster.csc.fi/manual/rseqc_infer_rnaseq_experiment.html).

I had already created a virtual environment to use python in so I did not need to recreate one but the steps I used are listed in my notes for [installing HTSeq](https://github.com/srmarzec/albopictus_remapping/blob/main/misc/SarahNotes.md#installing-htseq). Installing RSeQC:
```
cd $HOME
source python-env/bin/activate
cd python-env
pip install RSeQC
```
After I installed RSeQC I was able to test my file and see the stranded output. 
```
bin/infer_experiment.py -r /home/sm3679/albopictus_remap/albopictus_genome/aedes_albopictus_AalbF3.bed -i /home/sm3679/albopictus_remap/bam_dir/D_11d_rep1_Aligned.sortedByCoord.out.bam
```
The output looked equivalent to "Example 1: Pair-end non strand specific" from the [infer_experiment.py page](http://rseqc.sourceforge.net/#infer-experiment-py), thus confirming what I saw when looking at the htseq-count outputs.
