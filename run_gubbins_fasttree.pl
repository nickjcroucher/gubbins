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
no warnings 'uninitialized';
use File::Basename;
use Getopt::Long;

my $tree_building_exec = '/software/pathogen/external/apps/usr/local/fasttree/FastTree -gtr -gamma -nt ';
my $gubbins_exec = 'gubbins';


my($number_of_iterations, $alignment_file, $starting_tree, $help );

GetOptions(
   'i|iterations=s'            => \$number_of_iterations,
   'a|alignment_file=s'       => \$alignment_file,
   't|starting_tree=s'         => \$starting_tree,
   'h|help'                    => \$help,
    );

($alignment_file) or die <<USAGE;

Usage: $0
  -a|alignment_file     <Multi fasta alignment file>
  -i|iterations         <Iterations to run for, defaults to 5>
  -t|starting_tree      <Optional starting tree in Newick format>
  -h|help               <print this message>

This script takes a multifasta alignment file, finds all SNP sites, then iteratively eliminates recombinations, and rebuilds the tree.
USAGE

$number_of_iterations ||= 5;

# find SNP sites
#system("$gubbins_exec -s $alignment_file");

my($filename, $directories, $suffix) = fileparse($alignment_file,  qr/\.[^.]*/);
my $base_filename = $filename.$suffix;

my $iteration_base_name = $base_filename.".iteration_0";
rename("$base_filename.vcf", "$iteration_base_name.vcf");
rename("$base_filename.phylip", "$iteration_base_name.phylip");
rename("$base_filename.snp_sites.aln", "$iteration_base_name.snp_sites.aln");

for(my $i = 1; $i <= $number_of_iterations; $i++)
{
  $iteration_base_name = $base_filename.".iteration_$i";
  my $previous_iteration_base_name = $base_filename.".iteration_".($i-1)."";

 print "$tree_building_exec $previous_iteration_base_name.snp_sites.aln > $iteration_base_name.tre\n";
 print "$gubbins_exec -r $alignment_file $previous_iteration_base_name.vcf $iteration_base_name.tre $previous_iteration_base_name.phylip\n\n";
#  system("$tree_building_exec $previous_iteration_base_name.phylip > $iteration_base_name.tre");
#  system("$gubbins_exec -r $alignment_file $previous_iteration_base_name.vcf $iteration_base_name.tre $previous_iteration_base_name.phylip");
}

