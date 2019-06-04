# 3. Utils
#   3.1 COUNT2TPM
#   3.2 qqnorm
#   3.3 DESeq2

#' Title Transform a count table to TPM matrix
#'
#' @param COUNT a dataframe with each row represents a gene and each column represents a cell
#'
#' @return TPM Dataframe
#' @export
#'
#' @examples tpm <- COUNT2TPM(count)
COUNT2TPM <- function(COUNT){
  stopifnot(is.data.frame(COUNT))
  stopifnot(any(is.na(COUNT) | is.infinite(COUNT)))

  eff <- read.table("SuppData/Kallisto_genelen_mm20.txt",stringsAsFactors = F,header=T,row.names = 1)
  genes <- rownames(COUNT)
  eff_TPM <- sweep(COUNT,1,eff[genes,]/1000,"/")
  libsize <- apply(eff_TPM,2,sum)
  eff_TPM <- sweep(eff_TPM,2,libsize,"/")*1e6
  return(eff_TPM)
  print("eff_TPM is done!")
}

#' Title Quantile normalization for single cell expression profile
#'
#' @param x expression matrix row is gene col is cell
#'
#' @return Normalized expression profile
#' @export
#'
#' @examples qqnorm(TPM_matrix)
qqnorm <- function(x){

  stopifnot(is.data.frame(x))
  stopifnot(any(is.na(x) | is.infinite(x)))

  norm_matrix <- matrix(0,nrow = nrow(x),ncol=ncol(x))
  rownames(norm_matrix) <- rownames(x)
  colnames(norm_matrix) <- colnames(x)
  gene_num_dis <- apply(x, 2, function(x){sum(x>0)})
  max_ori <- which.max(gene_num_dis)

  for (i in 1:ncol(x) ){
    ratio <- x[,i]/x[,max_ori]
    ratio <- ratio[!is.infinite(ratio)] # remove infinite
    ratio <- ratio[!is.na(ratio)]
    # sizefactor <- median(ratio[ratio!=0])
    sizefactor <- median(ratio[ratio>10])
    norm_matrix[,i] <- x[,i]*sizefactor
  }

  norm_matrix <- norm_matrix[,!is.na(colSums(norm_matrix))]
  norm_matrix <- log2(norm_matrix+1)
  mean_ <- apply(norm_matrix,1,mean)
  norm_matrix <- sweep(norm_matrix,1,mean_,'-')

  norm_matrix <- norm_matrix[!duplicated(norm_matrix),]
  tnm <- t(norm_matrix)
  tnm <- tnm[!duplicated(tnm),]
  norm_matrix <- t(tnm)
  return(norm_matrix)
}

