#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

###split sequences into 4fold and 0fold sites
# Four-fold degenerate only at 3rd postion if CT-, GT-, -C-, CG-, GG-
# Zero-fold at all 3 positions for TGG and ATG
# Zero-fold degenerate at 2nd position always
# Watch out for STOP codons at TAA, TAG and TGA

# Create hash each for 4-fold, 0-fold 1st position, and 0-fold 3rd position.

my %four = (
#CT-
   CTA => '',
   CTC => '',
   CTG => '',
   CTT => '',
#GT-
   GTA => '',
   GTC => '',
   GTG => '',
   GTT => '',
#CG-
   CGA => '',
   CGC => '',
   CGG => '',
   CGT => '',
#GG-
   GGA => '',
   GGC => '',
   GGG => '',
   GGT => '',
#-C-
   CCA => '',
   CCC => '',
   CCG => '',
   CCT => '',
   GCA => '',
   GCC => '',
   GCG => '',
   GCT => '',
   ACA => '',
   ACC => '',
   ACG => '',
   ACT => '',
   TCA => '',
   TCC => '',
   TCG => '',
   TCT => '');
# Zero-fold degenerate at 1st position if: G--, AT-, AC-, AA-, AG(T/C), CC-, CA-, CT(T/C), CG(T/C), TC-, TT(T/C), TA(T/C), TG(T/C/G)
my %zero = (
#1st codon position zero-degenerate
#G--
   GAA => '',
   GAC => '',
   GAG => '',
   GAT => '',
   GCA => '',
   GCC => '',
   GCG => '',
   GCT => '',
   GGA => '',
   GGC => '',
   GGG => '',
   GGT => '',
   GTA => '',
   GTC => '',
   GTG => '',
   GTT => '',
#AT-
   ATA => '',
   ATC => '',
   ATG => '',
   ATT => '',
#AC-
   ACA => '',
   ACC => '',
   ACG => '',
   ACT => '',
#AA-
   AAA => '',
   AAC => '',
   AAG => '',
   AAT => '',
#AG
   AGT => '',
   AGC => '',
#CT
   CTT => '',
   CTC => '',
#CC-
   CCA => '',
   CCC => '',
   CCG => '',
   CCT => '',
#CA-
   CAA => '',
   CAC => '',
   CAG => '',
   CAT => '',
#CG
   CGT => '',
   CGC => '',
#TT-
   TTT => '',
   TTC => '',
#TC-
   TCA => '',
   TCC => '',
   TCG => '',
   TCT => '',
#TA-
   TAT => '',
   TAC => '',
#TG-
   TGC => '',
   TGG => '',
   TGT => ''
);

# Push codons into an array. then go through the array in a foreach loop, counting as you go.
# if the codon matches one of the four-fold codons record number of codon x 3 for position.
# same for 1st position for zero-fold, all 2nd positions will be kept. if codon matches TGG and ATG record 3rd position.
# If it matches a stop codon, stop and write warning of where it is.


open GENES, '<TAIR10_cds.fa';
chomp(my @lines = <GENES>);

open FOUR, '>four_fold';
open ZERO, '>zero_fold';

my $gene_string = join('',@lines);
my @genes = split (/>/, $gene_string);

my $gene;
my $name;
my @codons;
my $codon;
my $n;
my @four_positions;
my @zero_positions;
my $position;
my @second_positions;
my $length;

open STOP, ">prem_stops"; #print names of genes with prem stop codons

foreach $gene(@genes){
#match sequence except first (START) and last (STOP) codon
#>AT1G01010.1 gene=AT1G01010
     if($gene =~ /^(AT\S+)\s+gene\=AT\dG\d+ATG([ACTG-]+)[ACTG]{3}$/i){
#print "$1\n$2\n";
      $n = 1; #to account for start codon
      $name = $1;
      $length = length($2) + 3;
      @codons = ($2 =~ m/.../g);
      @four_positions = ();
      @zero_positions = ();
#now loop through the codons looking for 4-fold sites
      foreach $codon(@codons){
         $n += 1; #starting with second codon
         if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA"){
            print STOP "$name position $n\n";
	    #need to break the loop here so this is not included in the output files
         }
         elsif(exists $four{$codon}){
            $position = $n * 3;#3rd position in codon
            push @four_positions, $position; 
         }
         elsif(exists $zero{$codon}){
            $position = ($n * 3) - 2;#first position in codon
            push @zero_positions, $position; 
         }
         elsif($codon eq "TGG" || $codon eq "ATG"){#all three positions are 0-fold
            $position = ($n * 3) - 2;
            push @zero_positions, $position; 
            $position = ($n * 3) - 1;
            push @zero_positions, $position; 
            $position = $n * 3;
            push @zero_positions, $position; 
         }
      }
#now to get all 2nd positions
      @second_positions = ();
      for(my $i = 5; $i <= $length; $i += 3 ){
      push @second_positions, $i;
      }
      push @zero_positions, @second_positions;
      @zero_positions = sort { $a <=> $b } @zero_positions;
      @zero_positions = uniq @zero_positions;
      print FOUR "$name\t@four_positions\n";
      print ZERO "$name\t@zero_positions\n";
   }

}




