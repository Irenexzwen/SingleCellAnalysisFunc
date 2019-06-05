# 1. Summary
#   1.1 TopGenes
#   1.2 TopCells
#   1.3 Grouped dot+box plot (continue)

#' Title Returen Top expressed genes with highest sum(expr)
#'
#' @param COUNT a dataframe with each row represents a gene and each column represents a cell
#' @param top int(default 10) number of top genes to inspect
#'
#' @return DataFrame with genenames mean(expr),variance and dropout rate *expr==0/cell_number*
#' @export
#'
#' @examples tg <- TopGenes(ExprMatx)
TopGenes <- function(COUNT,top=10){
  stopifnot(is.data.frame(COUNT))
  stopifnot(any(is.na(COUNT) | is.infinite(COUNT)))

  #sort the dataframe by the row sums, and return the complete ordered dataframe
  names <- head(names(sort(apply(COUNT,1,sum),decreasing = T)),n=top)
  mean_ <- apply(COUNT[names,],1,mean)
  var_ <- apply(COUNT[names,],1,var)
  DropOutRate_ <- apply(COUNT[names,],1,function(x){sum(x==0)/length(x)})

  out <- data.frame('mean'=mean_,'var'=var_,'DropOutRate'=DropOutRate_,row.names = names)
  print(head(out))
  return(out)
}


#' Title Returen Cells expressed highest number of genes
#'
#' @param COUNT a dataframe with each row represents a gene and each column represents a cell
#' @param top int(default 10) number of top cells to inspect
#'
#' @return DataFrame with Cell names, expressed gene numbers, gene dropout rate
#' @export
#'
#' @examples tc <- TopCells(ExprMatx)
TopCells <- function(COUNT,top=10){
  stopifnot(is.data.frame(COUNT))
  stopifnot(any(is.na(COUNT) | is.infinite(COUNT)))

  #sort the dataframe by the row sums, and return the complete ordered dataframe
  cells <- head(names(sort(apply(COUNT,2,function(x){sum(x>0)}),decreasing = T)),n=top)
  genenum_ <- apply(COUNT[,cells],2,function(x){sum(x>0)})
  DropOutRate_ <- apply(COUNT[,cells],2,function(x){sum(x==0)/length(x)})

  out <- data.frame('GeneNum'=genenum_,'DropOutRate'=DropOutRate_,row.names = cells)
  print(head(out))
  return(out)
}

