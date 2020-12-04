#!/usr/bin/perl
use strict; use warnings;

#open TEST, ">test.out";
#split each gene sequence into four-fold and zero-fold sites only

#first open files containing positions of 4-fold and 0-fold sites and push them into arrays. These positions will be used as array indices on the sequences
open FOUR, "<four_fold";
chomp(my @four = <FOUR>);
open ZERO, "<zero_fold";
chomp(my @zero = <ZERO>);


#open fasta files, one at a time
my @files=`ls ../*_cds.aln`;
#my @files = `ls ../AT1G01180.1_cds.aln`;
foreach my $file (@files){
	chomp($file);
	$file =~ m/\.\.\/(\S+)_cds.aln/;
	my $gene = $1;

my $four_file="$gene".".four.fas";
my $zero_file="$gene".".zero.fas";


open (IN, "<$file") or die "Can't open file\n";
chomp(my @lines = <IN>);
#print TEST "@lines\n";
open FOUR_OUT, ">$four_file";
open ZERO_OUT, ">$zero_file";

my @seq;
my @four_seq;
my @zero_seq;
my $four_seq;
my $zero_seq;

###########################################################
######first test for presence of more than one SNP per codon. All codons with > 1 SNP will be recorded in a hash and not included when 4- and 0-fold sites are printed.

my %codons;#store array of codons for each allele name
my @test;
my $test;
my @test_bases;
my %first;
my %second;
my %third;
my $snp;
my %double_snps;
my $num_codons;
my $allele_name;
my @allele_codons;

#store codons and allele names in %codons
my $string = join '', @lines;
my @sequences = split />/, $string;
#print TEST "@sequences\n";
foreach(@sequences){
   if($_ =~ /^(variant_\d+)(atg\S+)$/){
	   #      print TEST "$1\n$2\n";
      $allele_name = $1;#get allele name
      @allele_codons = ($2 =~ m/.../g);#split sequence into codons
      #print TEST "$allele_name\t@allele_codons\n";
      $codons{$allele_name} = [@allele_codons];
   }
}
$num_codons = scalar @allele_codons;
#print TEST "@allele_codons\n";
for(my $i = 0; $i < $num_codons; $i ++){#go through all codons
      @test = ();# reset @test;
   foreach my $key (keys %codons) {#for each of the alleles (keys)
      push @test, $codons{$key}[$i];# save codon from position $i
   }
   #print TEST "@test\n";
   #for $i = 0, we have @test = (ATG, ATG, ATG, ATG, ...), we need to compare 1st, 2nd and 3rd bases across the alleles
   %first = ();#reset hashes
   %second = ();
   %third = ();
   @test_bases = ();
   foreach $test(@test){
      @test_bases = ($test =~ m/./g);#foreach allele split codon intro three bases
      $first{$test_bases[0]} = ();#store each position in a separate hash
      $second{$test_bases[1]} = ();
      $third{$test_bases[2]} = ();
   }
#now test if more than one of the hashes has more than 1 key
   $snp = 0;
   if (keys %first > 1) {#if there is more than one key in this hash there is a snp
      $snp += 1;
   }
   if (keys %second > 1) {
      $snp += 1;
   }
   if (keys %third > 1) {
      $snp += 1;
   }
   if($snp > 1){#if there is more than one snp we need to record the codon position
      $double_snps{$i*3+1} = ();
      $double_snps{$i*3+2} = ();
      $double_snps{$i*3+3} = ();#record all 3 bases positions of codon
   }
}
#print "$gene:\n";
#foreach (sort {$a<=>$b} keys %double_snps){
#print "$_\n";
#}

###################################################################
##now get 4- and 0-fold positions for this current gene if they are not in the double_snps hash
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
###########################################
#now print relevant bases for each allele

foreach(@sequences){
   if($_ =~ /^(variant_\d+)(atg\S+)$/){
      print FOUR_OUT ">$1\n";
      print ZERO_OUT ">$1\n";
      @seq = ( $2 =~ m/./g );
#we need to add random element to the array so that the indices coincide with real base positions
      unshift @seq, "X";
      @four_seq = ();#reset each loop
      @zero_seq = ();#reset each loop
#foreach element (=position) of gene array push that index of @seq into new array
      foreach (@four_pos){
         push @four_seq, $seq[$_];
      }
      $four_seq = join('',@four_seq);
      print FOUR_OUT "$four_seq\n";
      foreach (@zero_pos){
         push @zero_seq, $seq[$_];
      }
      $zero_seq = join('',@zero_seq);
      print ZERO_OUT "$zero_seq\n";
   }
}


}

