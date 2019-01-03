#!/usr/bin/perl 
####################################################################################################
# Name: Make_MappingFile_Qiime.pl
#
# Author: Rafael Bargiela, PhD. Bangor University. 2018.
#
# Resume: Perl script to make a mapping file suitable to be used with Qiime. The program uses 
# the sampleSheet file (in txt format) produced to run the Illumina MiSeq and convert it in a 
# mapping file joining the sequences of the two indexes used by the Dual-index sequencing protocol.

# Usage: perl Make_MappingFile_Qiime.pl <SAMPLESHEET.txt> > MappingFile.txt
####################################################################################################
use strict;
####################################################################################################
my %hash;
my %hash2;

print "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n";

my $data=0;
open(FILE,"<$ARGV[0]") || subdie($ARGV[0]);
  while(my $l=<FILE>){
    chomp($l);
    if($l=~/^Sample_ID/){
      $data=1;
    }
    else{
        if($data==1){
          my @a=split("\t",$l);
          $a[0]=~s/_//;
          my $sample=$a[0];
          my $barcode24=$a[5].$a[7];
          my $Des=$a[9];
          print "$sample\t$barcode24\t\t$Des\n"

        }
        else{
          next;
        }
    }
  }
close(FILE);

exit;
##############0######################################################################################
sub subdie{

  die"
####################################################################################################
# Name: Make_MappingFile_Qiime.pl
#
# Author: Rafael Bargiela, PhD. Bangor University. 2018.
#
# Resume: Perl script to make a mapping file suitable to be used with Qiime. The program uses 
# the sampleSheet file (in txt format) produced to run the Illumina MiSeq and converts it in a 
# mapping file joining the sequences of the two indexes used by the Dual-index sequencing protocol.

# Usage: perl Make_MappingFile_Qiime.pl <SAMPLESHEET.txt> > MappingFile.txt
####################################################################################################
  ";
}
