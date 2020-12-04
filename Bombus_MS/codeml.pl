#!/usr/bin/perl
use strict; use warnings;
use IO::Handle;

my @files=`ls *cds_aln.fa`;
foreach my $file (@files){
	chomp($file);
	$file =~ m/(\S+)cds_aln\.fa/;
	my $outfile = "$1".".out";
        my $ctlfile = "$1".".ctl";

open CONTROL, ">$ctlfile";
print CONTROL
	"seqfile = $file\n",
	"outfile = $outfile\n",
	"noisy = 0\n",
      	"verbose = 0\n",
      	"runmode = -2\n",
    	"cleandata = 1\n",
      	"seqtype = 1\n",
    	"CodonFreq = 2\n",
     	"model = 2\n",
	"NSsites = 0\n",
     	"icode = 0\n",
     	"Mgene = 0\n",
    	"fix_kappa = 0\n",
     	"kappa = 2\n",
    	"fix_omega = 0\n",
     	"omega = 1\n",
    	"fix_alpha = 1\n",
     	"alpha = .0\n",
       	"Malpha = 0\n",
     	"ncatG = 4\n",
       	"clock = 0\n",
     	"getSE = 0\n",
 	"RateAncestor = 0\n",
       	"method = 0\n";

}
