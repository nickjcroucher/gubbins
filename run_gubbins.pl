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
my $gubbins_exec = 'gubbins';
my $number_of_iterations = 5;

my $input_multi_fasta_alignment_file = $ARGV[0];

my($filename, $directories, $suffix) = fileparse($input_multi_fasta_alignment_file,  qr/\.[^.]*/);

# Initial step to find SNPs
system("$gubbins_exec -s $ARGV[0]");

my $current_time = time();

my $base_filename = $filename.$suffix;

for(my $i = 1; $i <= $number_of_iterations; $i++)
{
  my $current_tree = "RAxML_result.$filename.$current_time.iteration_$i";
  my $previous_tree_name = $base_filename;
  my $previous_tree = "";
  if($i > 1)
  { 
    $previous_tree_name = "RAxML_result.$filename.$current_time.iteration_".($i-1);
    $previous_tree = "-t $previous_tree_name";
    $base_filename = $current_tree;
   }
  system("$tree_builder_exec -s $previous_tree_name.phylip -n $filename.$current_time.iteration_$i $previous_tree");
  system("$gubbins_exec -r $input_multi_fasta_alignment_file $base_filename.vcf $current_tree $base_filename.phylip");
}

