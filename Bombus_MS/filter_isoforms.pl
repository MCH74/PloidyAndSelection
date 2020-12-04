#!/usr/bin/perl
use strict; use warnings;

my($fasta,$gff,$out)=@ARGV;

#find longest isoforms:
open GFF, "<$gff";
chomp(my @gff = <GFF>);

open ISO, ">longest_isoforms";

my $length;
my $locus;
my $rna;
my %rna_length;
my %isoforms;

foreach my $gff(@gff){
   if($gff =~ /^.+CDS\s+(\d+)\s+(\d+)\s+.+Parent=(rna\d+);.+gene=(\S+?);/){
	$length = $2-$1 + 1;
	$rna = $3;
	$locus = $4;
#	$rna_length{$mrna} = $length;
	if(exists $rna_length{$rna}){
        $rna_length{$rna} = $rna_length{$rna} + $length;
	}
	else{
        $rna_length{$rna} = $length;
	}
	if(exists $isoforms{$locus}){
		if($rna_length{$rna} > $rna_length{$isoforms{$locus}}){
			$isoforms{$locus} = $rna;
		}
	}
	else{
		$isoforms{$locus} = $rna;
	}
   }
}

foreach my $key(sort keys %isoforms){
    print ISO "$key\t$isoforms{$key}\n";
}

#filter fasta files to contain only longest isoforms



open IN, "<$fasta";
chomp(my @lines = <IN>);

my $line = join '', @lines;
my @seqs = split />/, $line;

open OUT, ">$out";

my %printed;

foreach my $seq(@seqs){
    if($seq =~ /^(rna\d+)\s+gene\=(\S+?\d)([A-Za-z]+)$/){  #for cds file
#    if($seq =~ /^(rna\d+)\s+gene\=(\S+?\d)([A-Z\.]+)$/){  #for protein file
            if($isoforms{$2} eq $1 ){
#               print OUT ">$1$2\n$3\n";   #for protein file
                print OUT ">$2\n$3\n";  #for cds file
        }
    }
    elsif($seq =~ /^(rna\d+)\s+gene\=(\S+?[ePfkmm])([A-Za-z]+)$/){  #for cds file
	    #    elsif($seq =~ /^(rna\d+)\s+gene\=(\S+?[ePfkmm])([A-Z\.]+)$/){  #for protein file
#print "$1\t$2\n";
    if($isoforms{$2} eq $1 ){
#               print OUT ">$1$2\n$3\n";   #for protein file
                print OUT ">$2\n$3\n";  #for cds file
        }
    }
}
