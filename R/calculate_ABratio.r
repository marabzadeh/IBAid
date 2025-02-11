#' calculate_ABratio
#'
#' Function to compute AI ratio and confidence intervals
#'
#' @param df_rna rna allele counts
#' @param df_dna dna allele counts
#' @param z
#' @param ABminCut
#' @param ABminCut_noCI
#' @param Threshold
#' @return abratio and confidence interval and the result of genes being balance/imbalance
#' @examples
#' calculate_ABratio(df_rna, df_dna, z, ABminCut, ABminCut_noCI)
calculate_ABratio <- function(df_rna, df_dna, z, ABminCut, ABminCut_noCI, Threshold) {

  if (is.null(df_rna) | is.null(df_dna)) return(NULL)

  ABmaxCut <- 1 / ABminCut
  ABmaxCut_noCI <- 1 / ABminCut_noCI

  a_count_rna <- df_rna$a_count
  b_count_rna <- df_rna$b_count
  a_count_dna <- df_dna$a_count
  b_count_dna <- df_dna$b_count

  total_rna <- a_count_rna + b_count_rna
  total_dna <- a_count_dna + b_count_dna
  VAF_RNA <- vector (mode = "integer", length = length(total_rna))
  VAF_DNA <- vector (mode = "integer", length = length(total_rna))

  altAllele_rna <- vector (mode = "integer", length = length(total_rna))
  altAllele_dna <- vector (mode = "integer", length = length(total_rna))

  SE <- vector (mode = "integer", length = length(total_rna))
  AB_Ratio <- vector (mode = "integer", length = length(total_rna))
  IB <- vector (mode = "integer", length = length(total_rna))
  IB <- unlist(lapply(IB, function(x) x<--2))

  for (i in 1:length(a_count_rna))
  {
    altAllele_rna[i] <-  a_count_rna[i]
    if (a_count_rna[i] < b_count_rna[i]){
      altAllele_dna[i] <- min(a_count_dna[i],b_count_dna[i]) }
    else {
      altAllele_dna[i] <- max(a_count_dna[i],b_count_dna[i]) }

    #A = Vr/Vd
    #B = (1-Vr)/(1-Vd)
    #A/B = Vr(1-Vd) / Vd(1-Vr)

    VAF_RNA[i] = altAllele_rna[i] / total_rna[i]
    VAF_DNA[i] = altAllele_dna[i] / total_dna[i]

    if (((VAF_RNA[i] == 1) & (VAF_DNA[i]==1)) | ((VAF_RNA[i] == 0) & (VAF_DNA[i]==0))) {          #For now consider if one
      AB_Ratio[i] <- 100
      SE[i] <- 0
    }# End if main

    else {
      #A = Vr/Vd
      #B = (1-Vr)/(1-Vd)

      if ((VAF_DNA[i]==0) | (VAF_RNA[i]==1)){
        AB_Ratio[i] <- 100#np.random.normal(100,0.1)
      }else{
        AB_Ratio[i] <- (VAF_RNA[i]*(1-VAF_DNA[i]))/(VAF_DNA[i]*(1-VAF_RNA[i]))
      }

      SE2_DNA <- (z/1+(z^2/total_dna[i]))^2*((VAF_DNA[i]*(VAF_DNA[i])/total_dna[i])+z^2/(4*total_dna[i]^2))
      SE2_RNA <- (z/1+(z^2/total_rna[i]))^2*((VAF_RNA[i]*(1-VAF_RNA[i])/total_rna[i])+z^2/(4*total_rna[i]^2))

      if ((VAF_DNA[i]==0) | (VAF_RNA[i]==1)){
        SE[i] <- 0.01#np.random.normal(1e-4,0.1)
      }else{
        SE[i] <- sqrt(((1-VAF_DNA[i])/VAF_DNA[i])^2 * (1/(1-VAF_RNA[i])^4) * SE2_RNA +
                        (VAF_RNA[i]/(1-VAF_RNA[i]))^2 * (1/VAF_DNA[i]^4) * SE2_DNA)
      }#end else
   }# End else main

    ## IB decision ==>  1/0/-1/-2

    if ((total_rna[i] <Threshold) | (total_dna[i] <Threshold) ){
      IB[i] <- -1
    }else{
      if (  ( ( ((AB_Ratio[i]>ABmaxCut) &  ((AB_Ratio[i]-SE[i])>1)) | (AB_Ratio[i]>ABmaxCut_noCI) )  |
              ( ((AB_Ratio[i]<ABminCut) &  ((AB_Ratio[i]+SE[i])<1)) | (AB_Ratio[i]<ABminCut_noCI) ) ) ){

          IB[i] <- 1 } #Imbalance
      else{IB[i] <- 0} # Balance
    }#end else -1
  } ## End for i

  ## MAKE THE RESULT FRAME
  gene_id <- df_rna$gene_id
  SE <- round(SE,2)
  VAF_RNA <- round(VAF_RNA,2)
  VAF_DNA <- round(VAF_DNA,2)

  results <- data.frame(gene_id, AB_Ratio, SE, IB,
                        VAF_RNA, a_count_rna, b_count_rna, total_rna,
                        VAF_DNA, a_count_dna, b_count_dna, total_dna  )

  colnames(results) <- c("Gene", 	"AB", 	"SE", 	"BI",
                         "VRna", 	"A-Rna", 	"B-Rna",	"DRna",
                         "VDna", 	"A-Dna", 	"B-Dna", "DDna")

  return(results)
}

