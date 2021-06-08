#!/usr/bin/perl 
####################################################################################################
# Name: Make_SampleMetadata_Qiime2.pl
#
# Author: Rafael Bargiela, PhD. Bangor University. 2020.
#
# Resume: Perl script to make the sample-metadata.tsv needed to run classification through Qiime2. 
# The program uses the sampleSheet file (in txt format) produced after running the Illumina MiSeq
# and converts it in a sample-metadata file.

# Usage: perl Make_SampleMetadata_Qiime2.pl <SAMPLESHEET.txt> > sample-metadata.tsv

# NOTE: there are barcodes sequences added on the bottom with their cooresponding ID. If they
#       don't correspond to those you are using, just add them on the hash made on the bottom by
#       'indexes' subroutine.
####################################################################################################
use strict;
####################################################################################################
my %in=indexes();


print "sample-id\tbarcode-sequence\tbarcodes-order\n";
print "#q2:types\tcategorical\tcategorical\n";

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
            my $order="$a[2] $a[1]";
           my $barcode24=$in{$a[2]}.$in{$a[1]}; # Reverse barcode first
          print "$sample\t$barcode24\t$order\n";

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
# Name: Make_SampleMetadata_Qiime2.pl
#
# Author: Rafael Bargiela, PhD. Bangor University. 2020.
#
# Resume: Perl script to make the sample-metadata.tsv needed to run classification through Qiime2. 
# The program uses the sampleSheet file (in txt format) produced after running the Illumina MiSeq
# and converts it in a sample-metadata file.

# Usage: perl Make_SampleMetadata_Qiime2.pl <SAMPLESHEET.txt> > sample-metadata.tsv
####################################################################################################
  ";
}

sub indexes{
  my %hash=(
    "F1"=>"CCTAAACTACGG",
    "F2"=>"TGCAGATCCAAC",
    "F3"=>"CCATCACATAGG",
    "F4"=>"GTGGTATGGGAG",
    "F5"=>"ACTTTAAGGGTG",
    "F6"=>"GAGCAACATCCT",
    "F7"=>"TGTTGCGTTTCT",
    "F8"=>"ATGTCCGACCAA",
    "F9"=>"AGGTACGCAATT",
    "F10"=>"ACAGCCACCCAT",
    "F11"=>"TGTCTCGCAAGC",
    "F12"=>"GAGGAGTAAAGC",
    "F13"=>"GTTACGTGGTTG",
    "F14"=>"TACCGCCTCGGA",
    "F15"=>"CGTAAGATGCCT",
    "F16"=>"TACCGGCTTGCA",
    "F17"=>"ATCTAGTGGCAA",
    "F18"=>"CCAGGGACTTCT",
    "F19"=>"CACCTTACCTTA",
    "F20"=>"ATAGTTAGGGCT",
    "F21"=>"GCACTTCATTTC",
    "F22"=>"TTAACTGGAAGC",
    "F23"=>"CGCGGTTACTAA",
    "F24"=>"GAGACTATATGC",
    "R1"=>"CCTAAACTACGG",
    "R2"=>"TGCAGATCCAAC",
    "R3"=>"CCATCACATAGG",
    "R4"=>"GTGGTATGGGAG",
    "R5"=>"ACTTTAAGGGTG",
    "R6"=>"GAGCAACATCCT",
    "R7"=>"TGTTGCGTTTCT",
    "R8"=>"ATGTCCGACCAA",
    "R9"=>"AGGTACGCAATT",
    "R10"=>"ACAGCCACCCAT",
    "R11"=>"TGTCTCGCAAGC",
    "R12"=>"GAGGAGTAAAGC",
    "R13"=>"GTTACGTGGTTG",
    "R14"=>"TACCGCCTCGGA",
    "R15"=>"CGTAAGATGCCT",
    "R16"=>"TACCGGCTTGCA",
    "R17"=>"ATCTAGTGGCAA",
    "R18"=>"CCAGGGACTTCT",
    "R19"=>"CACCTTACCTTA",
    "R20"=>"ATAGTTAGGGCT",
    "R21"=>"GCACTTCATTTC",
    "R22"=>"TTAACTGGAAGC",
    "R23"=>"CGCGGTTACTAA",
    "R24"=>"GAGACTATATGC"
      );
  return(%hash);
}
