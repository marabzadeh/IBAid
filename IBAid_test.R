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

## The input parameters:
# Threshold == The min number of reads considered as informative (for DNA and RNA counts) 
# ABminCut_noCI == This parameter defines the minimum (and maximum = 1/minimum) interval space based on which if the calculated upper/lower interval (SE in the output file) pass the minimum/maximum bound, the releative calculated AB will be considered as imbalance. This parameter can be defined based on each dataset to clean the output. The tighter the interval the more samples/genes pass as imbalance. 
# ABminCut == This parameter defines the minimum (and maximum = 1/minimum) interval cutt off which if the AB (The relative expression of Allele A counts on Allele B counts) pass that cut off, the gene/sample will be considered as imbalance.  

IBAid::process_samples(samples, input_rna, input_dna, output_dir,
                       z=1.96, ABminCut=1, ABminCut_noCI=0.6, Threshold = 10 )
## This will generate an output in the outs folder for each sample

## Output parameters (for each gene in each sample):
# Gene == Name of the gene
# AB == The relative expression of Allele A counts on Allele B counts
# SE == The confidence interval calculated by the method 
# BI == The output of the method :: 0 balance, 1 imbalance, -1/-2 not informative
# VRna == Allele frequency of RNA
# A-Rna == Number of reads for Allele A in RNA
# B-Rna == Number of reads for Allele B in RNA
# DRna == Total number of reads for both alleles in RNA (A-Rna + B-Rna)
# VDna == Allele frequency of DNA
# A-Dna == Number of reads for Allele A in DNA
# B-Dna == Number of reads for Allele B in DNA
# DDna == Total number of reads for both alleles in DNA (A-Dna + B-Dna)
