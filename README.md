Interval-Based Allelic Imbalance Detection (IB-Aid) is a quantitative framework that uses interval arithmetic to distinguish monoallelic from biallelic expression while normalizing for copy number and tumor purity. 
IB-Aid computes confidence intervals based on depth measurement uncertainty (for both DNA and RNA) to ensure robust allelic ratio estimation. 
To use the method, please download the git package as a library in R with the following commands:

library(devtools)
devtools::install_github("marabzadeh/IBAid", force = T)

IBAid_test.R includes guides to how to use the method and read the output. The input sample files are in inps folder. To use the method, we need the list of genes with allele counts A (a_count) and B (b_count), for both DNA and RNA, each in a separate file. If no counts for DNA is available, any constant value for both counts to show the balance for DNA works as input. 

![image](https://github.com/user-attachments/assets/8451d832-72d7-4833-84f3-de6868006fab)

