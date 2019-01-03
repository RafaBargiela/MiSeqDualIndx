# MiSeqDualIndx
## Pipeline to process  ILLUMINA MISEQ results produced by DUAL-INDEXING methodology

Following protocol is a mixture of the steps mentioned by Fadrosh et al. (Microbiome, 2014) at github https://github.com/igsbma/MiSeq16S and those mentioned in https://github.com/laxeye/BFXscripts/blob/master/16S/16S-DI-workflow.sh by Aleksei  Korzhenkov.

Steps from 1 to 4 complete the first phase of the pipeline, where the pair-end reads are joint, trimmed and cleaned, and barcodes are separated from the reads. Further steps are performed under Qiime to final demultiplexing, OTU picking and taxonomy assignation. Currently is prepared to use Qiime 1 but there will be an update in the future to performed the same steps under Qiime 2.

### Initial files:
-R1.fastq (Fastq file with the forward reads).
-R2.fastq (Fastq file with the reverse reads).
- sample.sheet.data.csv -> File with the SampleIDs and barcodes used in the MiSeq run.

### Dependencies:
We need to have installed the following tools. A good options could be using anaconda and create an environment with all the programs installed on it, to make sure that everything is in the right place.

- fastx toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
- seqtk (https://github.com/lh3/seqtk)
- SeqPrep (https://github.com/jstjohn/SeqPrep)
- tagcleaner (http://tagcleaner.sourceforge.net/)
- Qiime (http://qiime.org/)
- scripts from:  https://github.com/igsbma/MiSeq16S

Some of the scripts also require the use of biopython and python2, so be aware that some warnings or errors can prompt because of this issue.

### Step 1: trimming the barcodes and generating barcodes.fastq files and reads files with barcodes trimmed

a) Trimming barcodes from fastq files and generating barcodes.fastq files

> fastx_trimmer -i R1.fastq -f 1 -l 12 -Q 33 -o R1_barcode.fastq
> fastx_trimmer -i R2.fastq -f 1 -l 12 -Q 33 -o R2_barcode.fastq

> cat R1_barcode.fastq | fq_mergelines.pl > R1_barcode_temp
> cat R2_barcode.fastq | fq_mergelines.pl > R2_barcode_temp

> paste R1_barcode_temp R2_barcode_temp | awk -F"\t" '{print $5"\t"$2$6"\t"$3"\t"$4$8}' | fq_splitlines.pl > R1R2_barcode.fastq

b) Trimming the barcodes from original files and generating trimmed.fastq files

> seqtk trimfq -b 12 R1.fastq > R1_trimmed_seq.fastq
> seqtk trimfq -b 12 R2.fastq > R2_trimmed_seq.fastq


### Step2: assemble paired end read using SeqRep

> SeqPrep -r R1_trimmed_seq.fastq -f R2_trimmed_seq.fastq -1 R1_trimmed_SeqRep.fastq -2 R2_trimmed_SeqRep.fastq -s merged.fastq -A GATCGGAAGAGCACACG -B AGATCGGAAGAGCGTCGT -m 0.15 -n 0.8 -o 5

Step3: Clean and trim sequences from adapter/primers

In this step we need to use tagcleaner. The adapters correspond to the last part of the forward (tag5) and reverse (tag3) oligo sequence in the sample_sheet.csv. First command is to get statistics, and the second to really trim the sequences:

> tagcleaner.pl -fastq merged.fastq -verbose -stats -tag5 GGACTACHVGGGTWTCTAAT -tag3 TTACCGCGGCKGCTGVCAC

According to the stats displayed, you can choose values for mismatches.

> tagcleaner.pl -fastq merged.fastq  -tag5 GTGBCAGCMGCCGCGGTAA  -tag3 GGACTACHVGGGTWTCTAAT -mm3 XX -mm5 XX -out Merged.clean.fastq -nomatch 3

### Step4: Macthing up barcodes and merged reads

Note: Python2 and biopython are required in this step. Would be usefull to create a conda environment under python2.7.5 and install there biopython.

> fq_getPairAndOrphan1.8.py Merged.clean.fastq R1R2_barcode.fastq Reads.ready.fastq Barcodes.ready.fastq Orphan.fastq


## PIPELINE USING QIIME 1.

This pipeline is though to be used under Qiime1. There is also the possibility to use Qiime2, though.  You just have to change the scripts and adapt the input and mapping file to use the new version of Qiime.
As recommendation, is better to install qiime1 trough anaconda in a specific environment in order to avoid any problem with the numpy version of your system:

>  conda create -n qiime1 qiime matplotlib=1.4.3 mock nose -c bioconda

Then activate the environment and run the rest of the pipeline with the envionment activated:

> source activate qiime1

Once you finish all the process you can deactivate the environment using:

> source deactivate

### Step5: Create a mappting file for qiime. 

For this step it is possible to use the script by Aleksei Korzhenkov provided in https://github.com/laxeye/mapgen.

In this protocol we also provide a script to create a mapping file from the SampleSheet file created to run the MiSeq protocol. An example of mapping file is also provided.

Basically, the file has four basic columns (tab separated):
 #SampleID[tab]BarcodeSequence[tab]LinkerPrimerSequence[tab]Description.

In the first one with indicate the name of our samples, been careful to not use forbidding characters. In the BarcodeSequence column, we add a 24 barcode sequence, which will be the result of joining both Index sequences. LinkerPrimerSequence can go empty (should be removed in previous steps). Description column can be empty as well.

Use validation tool from qiime to check the file:

> validate_mapping_file.py -m mappingFile.txt -p -o validation-results/

It can return some errors, but will be probably because the lack of Linker sequence and the empty field on Description. Anyway, the program returns a corrected file on the validadtion-results dir. You can use the corrected version provided by the script in de output directory.

### Step 6: Demultiplexing the sequences

Now we going to run split_libraries_fastq.py, we must remember to use the option –barcode_type and set it in 24, which is the length of our current barcodes:

> split_libraries_fastq.py -i Reads.ready.fastq -b Barcodes.ready.fastq -m mappingFile_corrected.txt --barcode-type 24 -o split_output_dir

### Step 7: Picking OTUS and classify sequences

Once the demultiplexing process has ended we will get in the OUT directory a fasta file named seqs.fna, among other files. This is the file which we will use to get the OTUS and assign them to a taxonomy. In order to do this we will need to download as well the database from the ARB-SILVA website. 

Otu picking will be performed against the SILVA database as reference, using the 97% of confidence at all type of taxa (meaning Eukarya, Bacteria and Archaea). The picking method chosen is open against a reference database, which allows to pick OTUs taking the SILVA database as reference but also try to pick de novo OTUs with those reads which couldn’t be aligned against the database. The only problem is that this method consume more CPU resources. Also, this method assigns taxonomy by default, using the script assign_taxonomy.py after picking OTUs. We can change setting to run assign_taxonomy.py by our own after the picking OTUs process, using the option –suppress_taxonomy_assignment. Otherwise, we can use a parameters files ($PARAMETERS) to directly set specific parameters to pass to assign_taxonomy.py and run it with  pick_open_reference_otus.py. In this parameters file can be specified, for example, the taxonomy file to use as reference in the assignation (in our case, 97% of confidence for all type of taxa, for all levels. An example is attached). Check http://qiime.org/documentation/qiime_parameters_files.html. 

> pick_open_reference_otus.py -i seqs.fna -r silva132_97.fna -t  -o pickOtus_output_dir -a  -O 40 -p $PARAMETERS

Among all the output file created in the output directory, the most important is the out table named as otu_table_mc2_w_tax_no_pynast_failures.biom, which contains final OTUs selection after some filterings. Also, rep_set.fna is also important, so it has the reprsentative sequences for each OTU.

In case you want to run separately the taxonomy assignation you must run:

> pick_open_reference_otus.py -i seqs.fna -r silva132_97.fna -t  -o pickOtus_output_dir -a  -O 40 –suppress_taxonomy_assignment
Then, you can independently classify the OTU sequences:

> parallel_assign_taxonomy_blast.py -i rep_set.fna -r  silva132_97.fna -t consensus_taxonomy_all_levels.txt -o Taxonomy_assignation_output_dir -O 40

Regard that in both scripts we use the option -O, to set the number of CPU cores to use and make faster the process. In case of assign_taxonomy.py, if we want to use this option we have to run the alternative script parallel_assign_taxonomy_blast.py (there are also other option which use other algorithms like uclust or vsearch).

### Step 8: Diversity and taxonomy distribution analysis

Final step will be to analyze the final results form picking otus and taxonomy assignation. To do this we can use two different scripts: summarize_taxa_through_plots.py, which just analyzes the taxonomy distribution; and core_diversity_analyses.py, which analyzes the taxonomy distribution and also diversity parameters like rarefaction curves and diversity indexes.

- summarize_taxa_through_plots.py -i otu_table_mc2_w_tax_no_pynast_failures.biom -o Summary_output_dir

Inside the output directory, the most important files are the .html files into the taxa_summary_plots directory. These files allow you to check all the summary charts in the web browser and also download summary tables in text format.

- core_diversity_analyses.py -i  otu_table_mc2_w_tax_no_pynast_failures.biom -o core_diversity_output_dir

In the case of core_diversity_analyses.py, the main output file is the index.html file, which address you to through the web browser into the different taxonomy and diversity charts generated by the different analysis. This script provides a deeper analysis, but both should return similar results in the taxonomy distribution charts.

## PIPELINE USING QIIME 2

Will be edited in the future...
