## Mona Arabzadeh
## Test file for IBAid: Interval-Based Allelic Imbalance detection method

library(devtools)
devtools::install_github("marabzadeh/IBAid", force = T)

library(IBAid)

# Example usage
samples <- c("sample1", "sample2")  # Replace with actual sample names
input_rna <- "~/inps/" #path to the test input folder
input_dna <- "~/inps/" #path to the test input folder
output_dir <- "~/outs/" #path to the test output folder

IBAid::process_samples(samples, input_rna, input_dna, output_dir,
                       z=1.96, ABminCut=1.2, ABminCut_noCI=0.6, Threshold = 10 )
## This will generate an output in the outs folder for each sample

## Output parameters (for each gene in each samples):
# Gene == Name of the gene
# AB == The relative expression of Allele A counts on Allele B counts
# BI == The output of the method :: 0 balance, 1 imbalance, -1/-2 not informative
# VRna == Allele frequency of RNA
# A-Rna == Number of reads for Allele A in RNA
# B-Rna == Number of reads for Allele B in RNA
# DRna == Total number of reads for both alleles in RNA (A-Rna + B-Rna)
# VDna == Allele frequency of DNA
# A-Dna == Number of reads for Allele A in DNA
# B-Dna == Number of reads for Allele B in DNA
# DDna == Total number of reads for both alleles in DNA (A-Dna + B-Dna)
