#!/usr/bin/perl
use strict;
use warnings;
use IO::Handle;

#count number of SNPs in each alignment
open OUT, ">zerofold_snps";
print OUT "Gene\tnumber of sites\tnumber of SNPs\n";
my @files=`ls *.zero.fas`;
foreach my $file (@files){
   chomp($file);
   $file =~ m/(\S+).zero.fas/;
   my $gene = $1;

   open (IN, "<$file") or die "Can't open file\n";
   chomp(my @lines = <IN>);

   my $line = join('',@lines);
   my @sequences = split(/>/,$line);
   shift @sequences; #get rid of leading \n
   my $length = 0;
   my @thal_bases;
   my @lyr_bases;
#only continue if there are 2 sequences
   if(scalar @sequences == 2){

      if($sequences[0] =~ /^\S+\d([A-Za-z-]+)$/){
         @thal_bases = ( $1 =~ m/./g );
         $length = scalar @thal_bases;
      }

      if($sequences[1] =~ /^\S+?-PA([A-Za-z-]+)$/){
         @lyr_bases = ( $1 =~ m/./g );
      }
      my $count = 0;
      for(my $i = 0; $i <= scalar @thal_bases - 1; $i ++){
         if($thal_bases[$i] ne $lyr_bases[$i]){
            $count += 1;
         }
      }
      print OUT "$gene\t$length\t$count\n";

   }

}
