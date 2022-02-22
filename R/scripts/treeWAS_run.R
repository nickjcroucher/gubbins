###############################################################################
## Running treeWAS on gubbins results #########################################
###############################################################################

require(argparse, quietly = TRUE, warn.conflicts = FALSE)
if(!("treeWAS" %in% rownames(installed.packages()))){
  require(devtools, quietly = TRUE, warn.conflicts = FALSE)
  devtools::install_github("caitiecollins/treeWAS", build_vignettes = TRUE)
}
require(adegenet, quietly = TRUE, warn.conflicts = FALSE)
require(treeWAS, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(ape, quietly = TRUE, warn.conflicts = FALSE)

###############################################################################
## Functions ##################################################################
###############################################################################

get_input <- function(){
  parser <- ArgumentParser(description='Perform treeWAS on gubbins results')
  parser$add_argument('--aln', type="character", required = TRUE,
                       help='alignment of isolates (initial input to GUBBINS) ')
  parser$add_argument('--tree', dest='tree', required = TRUE,
                      help='gubbins node labelled tree')
  parser$add_argument('--phen',dest = "phen", type = "character", required = TRUE,
                      help = "two column csv, first column names of isolates, second column phenotype of interest, discrete or binary")
  parser$add_argument('--out', dest = 'out', type = "character", required = TRUE,
                      help = "Location for out pdf and RData objects")
  
  return(parser$parse_args())
}

main <- function(input_args){
  
  cat("\n", "Reading in Alignment ", "\n")
  input_fasta <- ape::read.dna(input_args$aln,
                                 format = "fasta")
  input_mat <- DNAbin2genind(input_fasta)@tab
  cat("Done ", "\n")
  ## Going through the biallelic loci test 
  cat("Subsetting biallelic loci", "\n")
  suffixes <- keepLastN(colnames(input_mat), n = 2)
  suffixes <- unique(suffixes)
  
  if(all(suffixes %in% c(".a",".t",".c",".g"))){
    snps_input <- get.binary.snps(input_mat)
  }
  cat("Done ", "\n")
   
  ## So thats the snp data now in the snps matrix
  ## I'll get the tree and the phenotype data
  cat("Reading in phenotypic data ", "\n")
  input_micro <- read.csv(input_args$phen,
                           stringsAsFactors = FALSE) 
  
  input_disease_phen <- input_micro[,2]
  input_phen <- as.vector(unlist(input_disease_phen))
  names(input_phen) <- input_micro[,1]
  
  ## check if phen names in alignment names
  if(!(all(names(input_phen) %in% rownames(snps_input)))){
    cat("Names of phenotypic data don't match alignment names, please reconfigure names", "\n")
    stop()
  }
  
  cat("Done", "\n")
  
  ## Get the tree loaded
  cat("Loading Tree", "\n")
  input_tree <- read.tree(input_args$tree)
  cat("Done", "\n")
  ## Ok lets try out a treeWAS run!
  cat("Beginning treeWAS run", "\n")
  pdf_name <- paste(input_args$out, ".pdf", sep = "")
  out_disease <- treeWAS(snps = snps_input,
                         phen = input_phen,
                         tree = input_tree, 
                         seed = 1,
                         phen.reconstruction = "ML",
                         snps.reconstruction = "ML",
                         filename.plot = pdf_name)
  
  cat("\n", "Done","\n")
  significant_snps <- length(out_disease$treeWAS.combined$treeWAS.combined)
  cat(sprintf("There are %s significant SNPs identified ", significant_snps), "\n")
  out_name <- paste(input_args$out, ".RData", sep = "")
  cat("Saving out_data to: ", out_name, "\n")
  save(out_disease, file =  out_name)
  cat("Done", "\n")
  
}

###############################################################################
## Main run ###################################################################
###############################################################################

input_args <- get_input()
main(input_args)


