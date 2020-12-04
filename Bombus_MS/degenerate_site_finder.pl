#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

###split sequences into 4fold and 0fold sites
# Four-fold degenerate only at 3rd postion if CT-, GT-, -C-, CG-, GG-
# Zero-fold at all 3 positions for TGG and ATG only
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
   TCT => ''
);
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
#TC-
   TCA => '',
   TCC => '',
   TCG => '',
   TCT => '',
#AG-
   AGT => '',
   AGC => '',
#CT-
   CTT => '',
   CTC => '',
#CG-
   CGT => '',
   CGC => '',
#TT-
   TTT => '',
   TTC => '',
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


open GENES, '<Bter_cds_loIso.fa';
chomp(my @lines = <GENES>);

open FOUR, '>four_fold';
open ZERO, '>zero_fold';
open STOP, '>prem_stops';

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
my %stops;
my %n_codon;

MAIN: foreach $gene(@genes){
#match sequence except first (START) and last (STOP) codon
   #if($gene =~ /^rna\d+\s+gene=(\S+\d)ATG(\S+)[ACTGactg]{3}$/ || $gene =~ /^rna\d+\s+gene=(\S+?[ePfkmm])ATG(\S+)[ACTGactg]{3}$/){
   if($gene =~ /^(\S*?)[Aa][Tt][Gg](\S+)[ACTGactg]{3}$/ ){
#print "$1\n$2\n";
      $name = $1;
      if(length($2) % 3 != 0){#only if full codons
print "$name not full codons\n";
         next MAIN;
      } 
      $n = 1; #to account for start codon
      $length = length($2) + 3;
      @codons = ($2 =~ m/.../g);
      @four_positions = ();
      @zero_positions = ();
      %n_codon = ();
#now loop through the codons looking for 4-fold sites
      foreach $codon(@codons){
         $n += 1; #starting with second codon
#if premature stop codon, stop loop
         if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA"){
            $stops{$name} = ();#record gene with premature stop codon to exclude it
         }
         else{
            if(exists $four{$codon}){
               $position = $n * 3;#3rd position in codon
               push @four_positions, $position; 
            }
            if(exists $zero{$codon}){
               $position = ($n * 3) - 2;#first position in codon
               push @zero_positions, $position; 
            }
            if($codon eq "TGG" || $codon eq "ATG"){#all three positions
               $position = ($n * 3) - 2;
               push @zero_positions, $position; 
               $position = ($n * 3) - 1;
               push @zero_positions, $position; 
               $position = $n * 3;
               push @zero_positions, $position; 
            }
            if($codon =~ /N/){#if it contains an N the codon must be removed, i.e. 2nd positions not included below
               $position = ($n * 3) - 1;
               $n_codon{$position} = ();
#print "$name\tN \= $position\n";
            }
         }
      }
#now to get all 2nd positions
      unless(exists $stops{$name}){
         @second_positions = ();
         for(my $i = 5; $i <= $length; $i += 3 ){
            unless(exists $n_codon{$i}){
               push @second_positions, $i;
            }
         }
         push @zero_positions, @second_positions;
         @zero_positions = sort { $a <=> $b } @zero_positions;
         @zero_positions = uniq @zero_positions;
         print FOUR "$name\t@four_positions\n";
         print ZERO "$name\t@zero_positions\n";
      }
   }

}

foreach(keys %stops){
print STOP "$_\n";
}


