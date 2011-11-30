#!/usr/bin/env perl

=head1 NAME
 
 run_gubbins.pl
 
=head1 SYNOPSIS
 
 run_gubbins.pl 
 
=head1 DESCRIPTION
 
 Assumes you have RaXML in your PATH.
 
=head1 CONTACT
=head1 METHODS
 
=cut

use strict;
use warnings;
use File::Basename;

my $tree_builder_exec = 'raxmlHPC -f d  -m GTRGAMMA';
#-n parisomony.tree.phylip.ra5 -t RAxML_result.sl.FINAL.aln.phylip.ra -s parisomony.tree.phylip
my $number_of_iterations = 5;

my $input_multi_fasta_alignment_file = $ARGV[0];

my($filename, $directories, $suffix) = fileparse($input_multi_fasta_alignment_file,  qr/\.[^.]*/);

# Initial step to find SNPs
system("gubbins -s $ARGV[0]");

for(my $i = 1; $i <= $number_of_iterations; $i++)
{
	my $previous_tree = '';
	if($i > 1)
	{
		$previous_tree  = "-t RAxML_result.iteration_".($i-1);
	}
	system("$tree_builder_exec -s $input_multi_fasta_alignment_file.phylip -n iteration_1 $previous_tree");
	
	# delete old temp files
	# run RaXML 
	system("gubbins -r $input_multi_fasta_alignment_file $input_multi_fasta_alignment_file.vcf RAxML_result.iteration_$i  $input_multi_fasta_alignment_file.phylip");
}
