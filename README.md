Please see our website for more information on [Gubbins](http://sanger-pathogens.github.io/gubbins/).


Gubbins
=======
Since the introduction of high-throughput, second-generation DNA sequencing technologies, there has been an enormous increase in the size of datasets being used for estimating bacterial population phylodynamics. Although many phylogenetic techniques are scalable to hundreds of bacterial genomes, methods which have been used for mitigating the effect of mechanisms of horizontal sequence transfer on phylogenetic reconstructions cannot cope with these new datasets. Gubbins (Genealogies Unbiased By recomBinations In Nucleotide Sequences) is an algorithm that iteratively identifies loci containing elevated densities of base substitutions while concurrently constructing a phylogeny based on the putative point mutations outside of these regions. Simulations demonstrate the algorithm generates highly accurate reconstructions under realistic models of short-term bacterial evolution, and can be run in only a few hours on alignments of hundreds of bacterial genome sequences.

Install
=======
Please see the INSTALL file for detailed instructions.

Running Gubbins
===============
To run Gubbins with default settings:

    run_gubbins.py [FASTA alignment]
    
Input options:

    --outgroup, -o	

The name of a sequence in the alignment on which to root the tree

    --starting_tree, -s	

A Newick-format starting tree on which to perform the first iteration analysis. The default is to compute a starting tree using RAxML

    --filter_percentage -f	

Filter out taxa with more than this percentage of missing data. Default is 25%
    
Processing options:

    --tree_builder, -t	
    
The algorithm to use in the construction of phylogenies in the analysis; can be ‘raxml’, to use RAxML, ‘fasttree’, to use Fasttree, or ‘hybrid’, to use Fasttree for the first iteration and RAxML in all subsequent iterations. Default is raxml

    --iterations, -i	
    
The maximum number of iterations to perform; the algorithm will stop earlier than this if it converges on the same tree in two successive iterations. Default is 5.

    --min_snps, -m	
The minimum number of base substitutions required to identify a recombination. Default is 3.

    --converge_method, -z
Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]. Default is weighted_robinson_foulds.
    
Output options:

    --use_time_stamp, -u	
    
Include a time stamp in the name of output files to avoid overwriting previous runs on the same input file. Default is to not include a time stamp.

    --prefix, -p	
    
Specifiy a prefix for output files. If none is provided it defaults to the name of the input FASTA alignment

    --verbose, -v	
    
Print debugging messages. Default is off.

    --no_cleanup, -n
    
Do not remove files from intermediate iterations. This option will also keep other files created by RAxML, fastml and fasttree, which would otherwise be deleted. Default is to only keep files from the final iteration.
    
Output files    
==========

Prefix
------

If a prefix is not defined with the –prefix option, the default prefix of the output files is:
X.Y

where:
X = Prefix taken from the input fasta file
Y = Time stamp. NOTE: This will only be included in the output file prefix if the –u flag has been selected



Output file suffices:
---------------------

    .recombination_predictions.embl

Recombination predictions in EMBL tab file format.

    .recombination_predictions.gff	

Recombination predictions in GFF3 format

    .branch_base_reconstruction.embl	

Base substitution reconstruction in EMBL tab format.

    .summary_of_snp_distribution.vcf	

VCF file summarising the distribution of SNPs

    .per_branch_statistics.csv	

Per branch reporting of the base substitutions inside and outside recombinations events.

    .filtered_polymorphic_sites.fasta	

FASTA format alignment of filtered polymorphic sites used to generate the phylogeny in the final iteration.

    .filtered_polymorphic_sites.phylip	

Phylip format alignment of filtered polymorphic sites used to generate the phylogeny in the final iteration.

    .final_tree.tre	

Final phylogenetic tree in newick format.

