#!/usr/bin/perl
use strict; use warnings;

#find longest isoforms:
open GFF, "<GCF_000214255.1_Bter_1.0_genomic.gff";
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

