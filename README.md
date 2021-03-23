# MiSeqDualIndx
## Pipeline to process  ILLUMINA MISEQ results produced by improved DUAL-INDEXING methodology

This protocol is based on the pipeline mentioned by Fadrosh et al. (Microbiome, 2014) at github https://github.com/igsbma/MiSeq16S and also implemented by Alesei Korzshenkov at https://github.com/laxeye/BFXscripts/blob/master/16S/16S-DI-workflow.sh. In this case, dual-indexing sequencing methodology is combined with the use of an heterogeneity spacer on the primer design to improve the quality of the reads. 

Currently, a preprocessing of the initial files is needed, where the pair-end reads are joint, trimmed and cleaned, and barcodes are separated from the reads. Further steps are performed under Qiime for final demultiplexing of the samples, picking of representative sequences (OTUs using Qiime1 or ASVs using Qiime2 with DADA2) and final taxonomic assignation.

## Processing index
[1.PREPROCESSING OF THE SEQUENCING RESULTS](#1.PREPROCESSING-OF-THE-SEQUENCING-RESULTS)
 [1.1.trimming the barcodes and generating barcodes.fastq files and reads files with barcodes trimmed](#1.1.trimming-the-barcodes-and-generating-barcodes.fastq-files-and-reads-files-with-barcodes-trimmed)
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
 
 Imput directory: _RAW/_ (containing sequences.fastq.gz and barcodes.fastq.gz)
 ```shell
 qiime tools import \
  --type EMPSingleEndSequences \
  --input-path RAW \
  --output-path sequences.qza
 ```
  Files created: _sequences.qza_
  
 Note that on the _--input-path_ option we write only the FULL PATH to the directory where are located the files, here we will assume that it is inside the folder where we are located. Also, regard that _--type_ is set with _EMPSingleEndSequences_, this is because we will follow the pipeline as if sequences were Single End sequences. We know that is not the case, but R1 and R2 have been assembled and merged on steps 1.2-1.4.
 
  ### 3.2.Demultiplexing
  
  Here we will separate every read to each sample, creating two new files, one containing the reads per sample and another one with details about the demultiplexing process.
  
  Imput files: _sequences.qza_
  ```shell
  qiime demux emp-single \
  --i-seqs sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --p-no-golay-error-correction \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza
  ```
 Files created: _demux.qza_ and _demux-details.qza_

Regard on the _--p-no-golay-error-correction_ option, which must be set on to **avoid golay correction**, since we have 24nt barcodes and this option performs a correction over 12nt barcodes. This could produce an error in the demultiplexing process.
Addtionally, you could need a summary of the demultiplexing process, which you can get as follows:

```shell
  qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
  mkdir DEMUX_SUMMARY
  qiime tools export --input-path demux.qzv --output-path DEMUX_SUMMARY
```
First command produce a _.qzv_ file which you could visualize using _Qiime view_. Second command uses this file to return different summary tables which will be send the the _output-path_ directory.
Also, you could need to get the individual fastq file, which you could get into a common directory with this command:

```shell
mkdir FASTQ_FILES
 qiime tools export --input-path demux.qza --output-path FASTQ_FILES
```

 ### 3.3.Quality control, denoising and picking representative sequences using DADA2

During this chimeras will be removed and reads will be filtered by length and quality. Also, dereplication will be performed, resulting in a set of uique Amplicon Sequence Variants (ASVs), representing the unique nucleotides sequence variants among the whole dataset of reads.

```shell
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len 0 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza
```
Files created: re-seqs-dada2.qza, table-dada2.qza and stats-dada2.qza

Regards here that we set truncation length by 0, meaning there won't be any truncation by length. This can be changed of course and also complete it with a filtration by quality value, but we have also take in mind that we perform some filtration and truncation processes during the preprocessing step, so maybe is not necessary. 
Three files are created in the process. _ stats-dada2.qza_ is an artifact contaning information about the filtering and read selection process. It could be visualize on _Qiime view_ or exported to get the information on manageble formats. But the most important outputs are the other two. _rep-seqs-dada2.qza_ contains the selected unique read sequences or ASVs and _table-dada2.qza_ is the summary of the abundance of these ASV along the samples. We will export them to the same directory:

```shell
   mkdir DADA2
    qiime tools export --input-path table-dada2.qza --output-path DADA2
    qiime tools export --input-path rep-seqs-dada2.qza --output-path DADA2
    
    biom convert -i DADA2/feature-table.biom -o DADA2/feature-table.tsv --to-tsv 
```
Files created: **_dna-sequences.fasta_**, _feature-table.biom_, **_feature-table.tsv_**

Hence, on _DADA2_ directory we will get all ASVs sequences in fasta format and also a _.tsv_ table showing the abundances by sample (equvalent to the previous OTU table).

### 3.4.Taxonomic classification of reads

We will assume that our sequencing product come from a 16S rRNA metabarcoding analysis, so we are interest on knowing the abundace in Bacteria and Archaea in our samples. Therefore, we will perform a taxonomic assignation using last version of SILVA database, for what we need a trained classifier Qiime2 artifact to run with classification command. There are precomputed classifiers available here: https://docs.qiime2.org/2019.10/data-resources. However, we would recommend to train previously our own classifier, so we could also use the very last version of the database we are interested. Also, maybe you want to classify using your custom set of sequences as database or to use a different database without precomputed classifier, so you would need also to train your own one.

#### 3.4.1. Training you own classifier

For this tutorial we will train a classifier to use the last available version of SILVA. To do this we will need two files: 1) DNA sequences on fasta format and 2) a taxonomy table file. Accession number must be the same in both files. In our case, both files are available on ARB-SILVA website. It could be possible that sequences are in RNA instead of DNA, so you will need to translate to DNA your sequences. Also, taxonomy file is just a tab separated file with the accessions on the first column/place and the taxonomy assigned in the second one. 
First, we need to import to Qiime2 both files (which we will call _SILVA.seqs.fasta_ and _SILVA.taxonomy.txt_):

```shell
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path SILVA.seqs.fasta \
  --output-path SILVA.seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path SILVA.taxonomy.txt \
  --output-path SILVA.taxonomy.qza
```
Files created: _SILVA.seqs.qza__ and _SILVA.taxonomy.qza_

Optionally, we can adapt our classifier to the type of sequences we usually analyse, extracting the reference sequences corresponding the specific length we usually manage in our reads and setting the primers that have use during the sequencing process:

```shell
    qiime feature-classifier extract-reads \
      --i-sequences SILVA.seqs.qza \
      --p-f-primer GTGBCAGCMGCCGCGGTAA \
      --p-r-primer GGACTACHVGGGTWTCTAAT \
      --p-min-length 150 \
      --p-max-length 450 \
      --o-reads SILVA.ref-seqs.qza
```
Files created: _SILVA.ref-seqs-qza_

The advantage of doing this is a classifier adapted to our reads, providing a faster classification over the specific amplified region of the 16S rRNA gene.

Finally, we run the training process of the classifier:

```shell
      qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads SILVA.ref-seqs.qza\
      --i-reference-taxonomy SILVA.axonomy.qza \
      --o-classifier SILVA.classifier.qza
  ```
  File created: **_SILVA.classifier.qza_**
  
  Regard that if we don't want to extract reference reads we just run the training process straight with the _SILVA.seqs.qza_ file rather than _SILVA.ref-seqs.qza_.

#### 3.4.2.Running the taxonomic classification

Once we have our classifier artifact, we can run the classification command to assign a taxonomy to the representative (ASVs) sequences from our samples selected by DADA2.

```shell
qiime feature-classifier classify-sklearn \
  --i-classifier  SILVA.classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.results.qza
  
    mkdir TAX
  qiime metadata tabulate --m-input-file taxonomy.results.qza --o-visualization "TAX/taxonomy.results.qzv
  qiime tools export --input-path taxonomy.results.qza --output-path TAX
```
File created: _taxonomy.results.qza_, _taxonomy.ressults.qzv, _taxonomy.tsv_

With the first command we run the classification itself, creating the artifact _.qza__ with the results. Second command creates a _.qzv_ file to check results on _qiime view_. Finally, _.qza_ file is exported to get a tab-separated text file with the taxonomy in _.tsv_ format. Regard that we previously have created a directory named _TAX_ to store there the results.

### 3.5.Summarizing results from taxonomic classification

Now we want to our results in a readable format, exportable to keep working with it out of _Qiime_. In the following command we will combine some of the files created in previous steps to, first, create a visualization artifact to use in _qiime view_ and check barplots and tables with the results at different levels of the classification done over the microbiota inhabiting our samples.

```shell
  qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.results.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```
File created: _taxa-bar-plots.qzv_

Finally, we will export the _.qzv_ artifact created to a directory, storing all the tables with taxonomies at different levels and the abundances of each of them in our samples.

```shell
 mkdir SUMMARY
    qiime tools export --input-path taxa-bar-plots.qzv --output-path SUMMARY
```

### 3.6.Joining everything together

Among all the steps, three are the most important files created: _dna-sequences.fasta_, _feature-table.tsv_ and _taxonomy.tsv_. These last two can be joined together in excel, as they should be in the same order with the same number of rows (ASVs IDs), so taxonomy can directly copy/paste as an additional column of the feature-table, getting a nice tipical _OTU_ table with the taxonomy added. Also, the _fasta_ file is useful if you need to check the sequence of any of the reads and _BLAST_ it or to do any other process.

Another useful files are those created inside the _SUMMARY_ directory, particularly the _level.6.csv_ file. Here, the taxonomy abundance described on the _taxonomy.tsv_ file is collapsed on unique taxonomies. In this case, unique taxonomies until genus level (6). Thiss is a good complement to stick next to the raw _otu table_, to check directly how many reads of a specific taxon are in our samples.



