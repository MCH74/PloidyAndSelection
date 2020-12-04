#!/usr/bin/perl
use strict;
use warnings;
use IO::Handle;

open OUTPUT, '>', "output.txt" or die $!;
open ERROR,  '>', "error.txt"  or die $!;

STDOUT->fdopen( \*OUTPUT, 'w' ) or die $!;
STDERR->fdopen( \*ERROR,  'w' ) or die $!;

open MUT, ">mutations";

open FOUR, "<../four_fold";
chomp(my @four = <FOUR>);
open ZERO, "<../zero_fold";
chomp(my @zero = <ZERO>);

open STOPS, "<../prem_stops";
chomp(my @prems = <STOPS>);

my %prems;

foreach(@prems){
   $prems{$_} = ();
}

#split each gene sequence into four-fold and zero-fold sites only
#first gaps in the alignment need to be processed. We need to remove all gaps from terrestris sequence and the corresponding bases in impatiens sequence for the 4-fold and 0-fold positions to be correct

#open alignment, record name of terrestris gene and push sequences into array
my @files=`ls *cds_aln.fa`;
MAIN: foreach my $file (@files){
   chomp($file);

   open (IN, "<$file") or die "Can't open file\n";
   chomp(my @lines = <IN>);
   my $line = join('',@lines);
   my @sequences = split(/>/,$line);
   shift @sequences; #get rid of leading \n
#only continue if there are 2 sequences
   if(scalar @sequences != 2){
print "$file not 2 seqs\n";
      next MAIN;
   }
#check for gaps within sequences and if so whether they are full codons with the reading frame
   my $current_position = 0; #keep track of current base position within the sequence
   my $gap_length = 0; #measure length of gap to check for divisibility by 3
#split terrestris sequence into bases
   my %frameshift;
   my @terr_bases = ();
   my @terr_gaps = ();
   my $gene;
   if($sequences[0] =~ /^(\S+\d)([A-Za-z-]+)$/ || $sequences[0] =~ /(\S+?[ePfkmm])([A-Za-z-]+)$/){
      $gene = $1;
      @terr_bases = ( $2 =~ m/./g );
   }
   else{
      print "match didn't work on terrestris seq $file\n";
      next MAIN;#if any other bases are encountered or no gene name, go to next file
   }
#only continue if gene isn't in the list of genes with prem stop codons
   if (exists $prems{$gene}){
print "$gene in prems\n";
      next MAIN;
   }

#check for gaps at end --> stop codon mutation $sequences[0] is terrestris and [1] is impatiens sequence
   my %stop;
   
   if($sequences[0] =~ /-$/){
      $stop{$gene} = "late SC";
   }
   elsif($sequences[1] =~ /-$/){
      $stop{$gene} = "premature SC";
   }
   else{
      $stop{$gene} = "no SC";
   }

   foreach(@terr_bases){
      if($_ =~ /\w/){   
#now test if previous gap is divisible by 3
         if($gap_length > 0){
            if($gap_length % 3 != 0 || $current_position % 3 != 0){
            $frameshift{$gene} = "frameshift";
            }
         }
         $current_position += 1; 
         $gap_length = 0;
      }
      elsif($_ eq "-"){
#         if($current_position % 3 != 0){#check the gap is within reading frame
#            $frameshift{$gene} = "frameshift";      
#         }
         $current_position += 1;
         push @terr_gaps, $current_position;#record gap positions in array
         $gap_length += 1;
      }
   }

#and now the same for impatiens

   my @imp_bases;
   my @imp_gaps;
   my $imp_gene;

   $current_position = 0;
   $gap_length = 0;

   if($sequences[1] =~ /^(\S+-PA)([A-Z-]+)$/){
      $imp_gene = $1;
      @imp_bases = ( $2 =~ m/./g );
   }
   else{
      print "$file impatiens match\n";
      next MAIN;
   }
   foreach(@imp_bases){
      if($_ =~ /\w/){
#now test if previous gap is divisible by 3
         if($gap_length > 0){
            if($gap_length % 3 != 0 || $current_position % 3 != 0){
            $frameshift{$gene} = "frameshift";
            }
         }
         $current_position += 1; 
         $gap_length = 0;
      }
      elsif($_ eq "-"){
#         if($current_position % 3 != 0){#check the gap is within reading frame
#            $frameshift{$gene} = "frameshift";        
#         }
         $current_position += 1;
         push @imp_gaps, $current_position;#record gap positions in array
         $gap_length += 1;
      }
   }

#print out mutations if they exist
   if(exists $frameshift{$gene}){
      print MUT "$gene\t$stop{$gene}\tFS\n";
   }
   else{
      print MUT "$gene\t$stop{$gene}\tno FS\n";
   }


#if gaps exist, record sequences without gaps for performing 4-fold/0-fold split
   my $index;
#unless(exists $frameshift{$gene}){
#remove terrestris gaps from both sequences. For the first removal we need to subtract 1 to get correct index, but after that the array becomes smaller and we need to keep subtracting one more
   my $n = 1;
   my $gap;
   if(scalar @terr_gaps > 1){
      foreach $gap(@terr_gaps){
         $index = $gap - $n;
         splice @terr_bases, $index, 1;
         $n += 1;	
      }
      $n = 1;
      foreach $gap(@terr_gaps){
         $index = $gap - $n;
         splice @imp_bases, $index, 1;
         $n += 1;	
      }
   }

#before splitting into 4- and 0-fold sites, check for codons with more than one SNP. Split reduced sequences (without terrestris gaps) into chunks of 3 and compare each codon between terrestris and impatiens. Loop through all elements of the codon arrays and split each one into three bases for the comparison. $a != $b = $n; if $n > 1, push that codon position into an array. Later when sorting out 4- and 0-fold positions, check positions are not in this array
   my $terr_bases = join('',@terr_bases);
   my $imp_bases = join('',@imp_bases);

   my @terr_codons = ($terr_bases =~ m/.../g);
   my @imp_codons = ($imp_bases =~ m/.../g);
   my @terr_codon_split;
   my @imp_codon_split;
   my $count;
   my %double_snps = ();
   for(my $i = 0; $i <= scalar @terr_codons - 1; $i ++){
      @terr_codon_split = ($terr_codons[$i] =~ m/./g);
      @imp_codon_split = ($imp_codons[$i] =~ m/./g);
      $count = 0;
      for(my $n = 0; $n <= 2; $n ++){
         if($terr_codon_split[$n] ne $imp_codon_split[$n]){
            $count += 1;
         }
      }
      if($count > 1){
         $double_snps{$i*3+1} = ();
         $double_snps{$i*3+2} = ();
         $double_snps{$i*3+3} = ();#record all 3 bases positions of codon
      }
   }

#now get 4- and 0-fold sites
#first open files containing positions of 4-fold and 0-fold sites and push them into a hash of arrays, with gene name as key and positions as array. Remove all positions which are included in @double_snps. These positions will be used as array indices on the sequences
   my @four_pos = ();#array of four-fold positions
   my @zero_pos = ();#array of zero-fold positions
   my @temp = ();

   foreach(@four){
      if($_ =~ /^$gene\s+(\S+.+)$/){
         @temp = split ' ' , $1;
      }
   }
   foreach(@temp){
      unless(exists $double_snps{$_}){
         push @four_pos, $_;
      }
   }

   @temp = ();
   foreach(@zero){
      if($_ =~ /^$gene\s+(\S+.+)$/){
         @temp = split ' ' , $1;
      }
   }
   foreach(@temp){
      unless(exists $double_snps{$_}){
         push @zero_pos, $_;
      }
   }


#create and open out files
   my $four_file="$gene".".four.fas";
   my $zero_file="$gene".".zero.fas";

   open FOUR_OUT, ">$four_file";
   open ZERO_OUT, ">$zero_file";

   my @four_seq;
   my @zero_seq;
   my $four_seq;
   my $zero_seq;

##we need to add random element to each sequence array so that the indices coincide with real base positions
#terrestris
   unshift @terr_bases, "X";
   @four_seq = ();#reset each loop
   @zero_seq = ();#reset each loop   
#foreach element (=position) of gene array push that index of @seq into new array
   foreach (@four_pos){
      push @four_seq, $terr_bases[$_];
   }
   $four_seq = join('',@four_seq);
   print FOUR_OUT ">$gene\n$four_seq\n";
   foreach (@zero_pos){
      push @zero_seq, $terr_bases[$_];
   }
   $zero_seq = join('',@zero_seq);
   print ZERO_OUT ">$gene\n$zero_seq\n";
   
#impatiens
   unshift @imp_bases, "X";
   @four_seq = ();#reset each loop
   @zero_seq = ();#reset each loop
#foreach element (=position) of gene array push that index of @seq into new array
   foreach (@four_pos){
      push @four_seq, $imp_bases[$_];
   }
   $four_seq = join('',@four_seq);
   print FOUR_OUT ">$imp_gene\n$four_seq\n";
   foreach (@zero_pos){
      push @zero_seq, $imp_bases[$_];
   }
   $zero_seq = join('',@zero_seq);
   print ZERO_OUT ">$imp_gene\n$zero_seq\n";
   
#}#end of unless command



}#end of MAIN loop
