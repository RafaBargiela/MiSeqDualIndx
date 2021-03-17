# MiSeqDualIndx
## Pipeline to process  ILLUMINA MISEQ results produced by improved DUAL-INDEXING methodology

This protocol is based on the pipeline mentioned by Fadrosh et al. (Microbiome, 2014) at github https://github.com/igsbma/MiSeq16S and also implemented by Alesei Korzshenkov at https://github.com/laxeye/BFXscripts/blob/master/16S/16S-DI-workflow.sh. In this case, dual-indexing sequencing methodology is combined with the use of an heterogeneity spacer on the primer design to improve the quality of the reads. 

Steps from 1 to 4 complete the first phase of the pipeline, where the pair-end reads are joint, trimmed and cleaned, and barcodes are separated from the reads. Further steps are performed under Qiime to final demultiplexing, OTU picking and taxonomy assignation. Currently is prepared to use Qiime 1 but there will be an update in the future to performed the same steps under Qiime 2.

### Initial files:
- R1.fastq (Fastq file with the forward reads).
- R2.fastq (Fastq file with the reverse reads).
- sample.sheet.data.csv -> File with the SampleIDs and barcodes used in the MiSeq run.

### Dependencies:
We need to have installed the following tools. A good options could be using anaconda and create an environment with all the programs installed on it, to make sure that everything is in the right place.

- fastx-toolkit v0.0.14 (http://hannonlab.cshl.edu/fastx_toolkit/)
- seqtk v1.3 (https://github.com/lh3/seqtk)
- SeqPrep v0.1 (https://github.com/jstjohn/SeqPrep)
- tagcleaner v0.16 (http://tagcleaner.sourceforge.net/)
- Qiime 1.9.1 (http://qiime.org/) or
- Qiime 2021.2 (https://docs.qiime2.org/)
- scripts from:  https://github.com/igsbma/MiSeq16S

Some of the scripts also require the use of biopython and python2, so be aware that some warnings or errors can prompt because of this issue.

## 1.PREPROCESSING OF THE SEQUENCING RESULTS

### 1.1.trimming the barcodes and generating barcodes.fastq files and reads files with barcodes trimmed

a) Trimming barcodes from fastq files and generating barcodes.fastq files
```shell
fastx_trimmer -i R1.fastq -f 1 -l 12 -Q 33 -o R1_barcode.fastq
fastx_trimmer -i R2.fastq -f 1 -l 12 -Q 33 -o R2_barcode.fastq

cat R1_barcode.fastq | fq_mergelines.pl > R1_barcode_temp
cat R2_barcode.fastq | fq_mergelines.pl > R2_barcode_temp

paste R1_barcode_temp R2_barcode_temp | awk -F"\t" '{print $5"\t"$2$6"\t"$3"\t"$4$8}' | fq_splitlines.pl > R1R2_barcode.fastq
```
b) Trimming the barcodes from original files and generating trimmed.fastq files
```shell
seqtk trimfq -b 12 R1.fastq > R1_trimmed_seq.fastq
seqtk trimfq -b 12 R2.fastq > R2_trimmed_seq.fastq
```

### 1.2.Assembling paired-end reads using SeqRep
```shell
SeqPrep -f R1_trimmed_seq.fastq -r R2_trimmed_seq.fastq -1 R1_trimmed_SeqRep.fastq.gz -2 R2_trimmed_SeqRep.fastq.gz -s merged.fastq.gz -m 0.15 -n 0.8 -o 5
```
### 1.3.Cleaning and trimming adapter/primers sequences

In this step we need to use tagcleaner. The adapters correspond to the last part of the forward (tag5) and reverse (tag3) oligo sequence in the sample_sheet.csv. First command is to get statistics, and the second to really trim the sequences:
```shell
tagcleaner.pl -fastq merged.fastq -verbose -stats -tag5 GTGBCAGCMGCCGCGGTAA  -tag3 GGACTACHVGGGTWTCTAAT
```
You will be prompted for some statistics which could help you to choose the values for maximum mismatch on for the next step. Values between 0-2 are OK for $mm3 and $mm5.
**NOTE:** If the statistics show that your reads start to group from values higher than 2, you can SKIP this step. Probably that means that the adapter has been removed during the sequencing process.
```shell
tagcleaner.pl -fastq merged.fastq  -tag5 GTGBCAGCMGCCGCGGTAA  -tag3 GGACTACHVGGGTWTCTAAT -mm3 $mm3 -mm5 $mm5 -out Merged.clean.fastq -nomatch 3
```
### 1.4.Macthing up barcodes and merged reads

Note: Python2 and biopython are required in this step. Would be usefull to create a conda environment under python2.7.5 and install there biopython.
```shell
fq_getPairAndOrphan1.8.py Merged.clean.fastq R1R2_barcode.fastq Reads.ready.fastq Barcodes.ready.fastq Orphan.fastq
```

## 2.PROCESSING SEQUENCING READS USING QIIME 1.

This pipeline is thought to be used under Qiime1. There is also the possibility to use Qiime2, though.  You just have to change the scripts and adapt the input and mapping file to use the new version of Qiime.
As recommendation, is better to install qiime1 through anaconda in a specific environment in order to avoid any problem with the numpy version of your system:
```shell
conda create -n qiime1 qiime matplotlib=1.4.3 mock nose -c bioconda
```
Then activate the environment and run the rest of the pipeline with the envionment activated:
```shell
source activate qiime1
```
Once you finish all the process you can deactivate the environment using:
```shell
source deactivate
```
### 2.1.Create a mappting file for qiime. 

For this step it is possible to use the script by Aleksei Korzhenkov provided in https://github.com/laxeye/mapgen.

In this protocol we also provide a script to create a mapping file from the SampleSheet file created to run the MiSeq protocol. An example of mapping file is also provided.

Basically, the file has four basic columns (tab separated):
 #SampleID[tab]BarcodeSequence[tab]LinkerPrimerSequence[tab]Description.

In the first one, with the name of our samples, be careful to not use forbidding characters. In the BarcodeSequence column, we add a 24 barcode sequence, which will be the result of joining both Index sequences. LinkerPrimerSequence can go empty (should be removed in previous steps). Description column can be empty as well.

Use validation tool from qiime to check the file:
```shell
validate_mapping_file.py -m mappingFile.txt -p -o validation-results/
```
It can return some errors, but will be probably because the lack of Linker sequence and the empty field on Description. Anyway, the program returns a corrected file on the validadtion-results dir. You can use the corrected version provided by the script in de output directory.

### 2.2.Demultiplexing the sequences

Now we going to run split_libraries_fastq.py, we must remember to use the option –barcode_type and set it in 24, which is the length of our current barcodes:
```shell
split_libraries_fastq.py -i Reads.ready.fastq -b Barcodes.ready.fastq -m mappingFile_corrected.txt --barcode-type 24 -o split_output_dir
```
In terms of a future submission to any database of the reads in fastq QIIME offers an interesting option called store_demultiplexed_fastq. With this flag an aditional fastq file is generated with the same reads contained in the seqs.fna generated as default. This is specially interesting when several samples from different experiments are analyzed in the same sequencing run and we need a specific fastq file for each experiment samples (i.e. to submit to Genebank):
```shell
split_libraries_fastq.py -i Reads.ready.fastq -b Barcodes.ready.fastq -m mappingFile_corrected.txt --barcode-type 24 -o split_output_dir --store_demultiplexed_fastq
```
It is possible to separate the reads belonging to each group of samples using filter_fasta.py. First, you will need to add the list of samples of each group in a text list and then run it like this:
```shell
filter_fasta.py -f seqs.fna -o Group.seqs.fna --sample_id_fp List.txt
```


### 2.3.Picking OTUS and classify sequences

Once the demultiplexing process has ended we will get in the OUT directory a fasta file named seqs.fna, among other files. This is the file which we will use to get the OTUS and assign them to a taxonomy. In order to do this we will need to download as well the database from the ARB-SILVA website. 

Otu picking will be performed against the SILVA database as reference, using the 97% of confidence at all type of taxa (meaning Eukarya, Bacteria and Archaea). The picking method chosen is open against a reference database, which allows to pick OTUs taking the SILVA database as reference but also try to pick de novo OTUs with those reads which couldn’t be aligned against the database. The only problem is that this method consume more CPU resources. Also, this method assigns taxonomy by default, using the script assign_taxonomy.py after picking OTUs. We can change setting to run assign_taxonomy.py by our own after the picking OTUs process, using the option –suppress_taxonomy_assignment. Otherwise, we can use a parameters files ($PARAMETERS) to directly set specific parameters to pass to assign_taxonomy.py and run it with  pick_open_reference_otus.py. In this parameters file can be specified, for example, the taxonomy file to use as reference in the assignation (in our case, 97% of confidence for all type of taxa, for all levels. An example is attached). Check http://qiime.org/documentation/qiime_parameters_files.html. 
```shell
pick_open_reference_otus.py -i seqs.fna -r silva132_97.fna -t  -o pickOtus_output_dir -a  -O 40 -p $PARAMETERS
```
Among all the output files created in the output directory, the most important is the table named as otu_table_mc2_w_tax_no_pynast_failures.biom, which contains final OTUs selection after some filterings. Also, rep_set.fna is also important, so it has the reprsentative sequences for each OTU.

In case you want to run separately the taxonomy assignation you must run:
```shell
pick_open_reference_otus.py -i seqs.fna -r silva132_97.fna -t  -o pickOtus_output_dir -a  -O 40 –suppress_taxonomy_assignment
```
Then, you can independently classify the OTU sequences:
```shell
parallel_assign_taxonomy_blast.py -i rep_set.fna -r  silva132_97.fna -t consensus_taxonomy_all_levels.txt -o Taxonomy_assignation_output_dir -O 40
```
Regard that in both scripts we use the option -O, to set the number of CPU cores to use and make faster the process. In case of assign_taxonomy.py, if we want to use this option we have to run the alternative script parallel_assign_taxonomy_blast.py (there are also other option which use other algorithms like uclust or vsearch).

The most important output is the *otu_table_mc2_w_tax_no_pynast_failures.biom* file. This table provide the abundance of reads for each otu on each sample, with the related taxonomy. *mc2* means minimum count of 2, so it only shows those OTUs with a number of reads >1. It is very convenient to convert this file in a most readable format, using the following command:
```shell
biom convert -i otu_table_mc2_w_tax_no_pynast_failures.biom -o Final_otu_table.txt --to-tsv --header-key taxonomy
```
Option *--to-tsv* is to convert it into a tab-delimited table, while *--header-key taxonomy* is needed to add the taxonomy assigned to each otu at the end of the table.

Once you get the OTU table, it is possible to get specific OTU sequences belonging to specific sample groups, using filter_fasta.py. First, you need to make list of OTUs ids on a text file:
```shell
filter_fasta.py -f OTUs_rep_set.fna -o OTUs.out.fna -s OTUs.List.txt
```

### 2.8.Diversity and taxonomy distribution analysis

Final step will be to analyze the final results form picking otus and taxonomy assignation. To do this we can use two different scripts: summarize_taxa_through_plots.py, which just analyzes the taxonomy distribution; and core_diversity_analyses.py, which analyzes the taxonomy distribution and also diversity parameters like rarefaction curves and diversity indexes.
```shell
summarize_taxa_through_plots.py -i otu_table_mc2_w_tax_no_pynast_failures.biom -o Summary_output_dir
```
Inside the output directory, the most important files are the .html files into the taxa_summary_plots directory. These files allow you to check all the summary charts in the web browser and also download summary tables in text format.
```shell
core_diversity_analyses.py -i  otu_table_mc2_w_tax_no_pynast_failures.biom -o core_diversity_output_dir -m mappingFile_corrected.txt -e 20 --nonphylogenetic_diversity
```
In the case of core_diversity_analyses.py, the main output file is the index.html file, which address you through the web browser into the different taxonomy and diversity charts generated by the analysis. This script provides a deeper study, but both should return similar results in the taxonomy distribution charts. In this option you need to add additional parameters like the mapping file and the sampling depth to use for even sub-sampling and maximum rarefaction depth. On Qiime website they set the -e parameter on 20, but you could use the command biom summarize-table to decide value.

## 3. PIPELINE USING QIIME 2

For the analysis using QIIME2 we start from Step 4, once all preprocessing is done. So, we have know all our reads and barcodes in the files **Reads.ready.fastq** and **Barcodes.ready.fastq**. To avoid any possible issues with names and location using QiIME2 we will:
- compress these files with *gzip*
- rename them to **sequences.fastq.gz** and **barcodes.fastq.gz**
- Create a new folder with only these two files on it, which we will call *RAW* on this tutorial

Also, we need to create a *sample-metadata.tsv* file to use with QIIME2, a bit different from the *mappingFile.txt* used with QIIME1. This file has the following tab-separated format:
`
 sample-id barcode-sequence description
 #q2:types categorical categorical
 SAMPLE1 BARCODE_SEQUENCE OPTIONAL_DATA
 `
 It is very similar to the mapping file used on QIIME1, with the 24bp barcodes sequences (Reverse-forward joined) in the second column. Be careful with the first column name, QIIME2 only accepts some different variants.
 As it was with QIIME1, it is also recommended to install QIMME2 under a conda environment. You can do it as follows:
 ```shell
 wget https://data.qiime2.org/distro/core/qiime2-2021.2-py36-linux-conda.yml
conda env create -n qiime2-2021.2 --file qiime2-2021.2-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2021.2-py36-linux-conda.yml
```
 Then, activate the environment:
 ```shell
 conda activate qiime2-2021.2
 ```
 
 ### 3.1.Importing files to QIIME2 environment
 
 Working with QIIME2 is completely different than QIIME1. It works with specific _artifacts_ that are modified and transformed inside its own environment. So, firt we need to create a specific _.qza_ file, which will include our sequences and barcodes, for further processing on QIIME2.
 
 ```shell
 qiime tools import \
  --type EMPSingleEndSequences \
  --input-path RAW \
  --output-path sequences.qza
 ```
 
 Note that on the _--input-path_ option we write only the FULL PATH to the directory where are located the files, here we will assume that it is inside the folder where we are located. Also, regard that _--type_ is set with _EMPSingleEndSequences_, this is because we will follow the pipeline as if sequences were Single End sequences. We know that is not the case, but R1 and R2 have been assembled and merged on steps 1.2-1.4.
 
 
 
Will be edited in the future...
