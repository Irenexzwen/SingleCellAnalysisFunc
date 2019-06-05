# 2.Cell and gene cleaning
#   2.1 Filter genes by min cells
#   2.2 Filter cell by summary Cutoff - FilterCellBySumCut

#' Title Delete genes expressed lower than min_cell.
#' if min_cell=0, then it will delete genes with zero expression level across all cells.
#' @param count_table A dataframe with each row represents a gene and each column represents a cell
#' @param min_cell int (default 3), minimum number of cells that a gene manifests.
#'
#' @return Filtered single cell expressed profile.
#' @export
#'
#' @examples
FilterGenesMinCells <- function(count_table,min_cell=3){
  # delete genes expressed lower than min_cell.
  name <- apply(count_table,1,function(x){sum(x>0)<min_cell})
  tmp <- count_table[!name,]
  return(tmp)
}



#' Title Quality Control on cells by summary matrix
#'
#' @param summary_ Str, path to summary file (with header and rownames)
#' @param uniqmp list of length two. u[1] is the lower bound, and u[2] is the upper bound.
#' @param prgene list of length two.
#' @param mt list of length two.
#' @param ercc list of length two.
#' @param rRNA list of length two.
#'
#' @return list, list of filtered cell names.
#' @export
#'
#' @examples filtered_cells <- FilterCellBySumCut(sum,uniqmp=c(60,100),prgene=c(4000,12000),mt=c(0,10),ercc=c(0,15,rRNA=c(0,20)))
FilterCellBySumCut <- function(summary_,uniqmp=u,prgene=p,mt=m,ercc=e,rRNA=r){
  library(dplyr)
  if(!is.list(u)) stop("uniqmp must be a positive list with length 2")
  if(!is.list(p)) stop("prgene must be a positive list with length 2")
  if(!is.list(m)) stop("mt must be a positive list with length 2")
  if(!is.list(e)) stop("ercc must be a positive list with length 2")
  if(!is.list(r)) stop("rRNA must be a positive list with length 2")

  summary <- read.table(summary_,header = T,row.names = 1,stringsAsFactors = F)
  stopifnot(is.data.frame(summary))

  summary['names'] <- rownames(summary)
  pr_out <- filter(summary,(ProteinCodingGenes < p[1] | ProteinCodingGenes > p[2]))$names
  uniqmap_out <- filter(summary,(UniqueMapRatio... < u[1] | UniqueMapRatio... > u[2]))$names
  mt_name <- filter(summary, (MT.expr...<m[1] | MT.expr...>m[2]))$names
  ercc_name <- filter(summary,(ERCC... >e[2] | ERCC... <e[1]))$names
  rrna <- filter(summary,(rRNA.. >r[2] | rRNA... <r[1]))$names
  filtered_cells <- list(pr=pr_out,uniqmap=uniqmap_out,mt=mt_name,ercc=ercc_name,rRNA=rrna)
  return(filtered_cells)
}


#' Title Delete rows by names of a dataframe
#'
#' @param DF dataframe
#' @param row.names.remove names of rows to be removed
#'
#' @return dataframe
#' @export
#'
#' @examples clean_cexpr <- DeleteRowsByNames(expr,c("Lars2","Gm10800"))
DeleteRowsByNames <- function(DF,row.names.remove){
  stopifnot(is.data.frame(DF))
  stopifnot(row.names.remove %in% row.names(DF))
  DF <- DF[!(row.names(DF) %in% row.names.remove), ]
  return(DF)
}


#' Title Delete rows by names of a dataframe
#'
#' @param DF dataframe
#' @param col.names.remove names of col to be removed
#'
#' @return dataframe
#' @export
#'
#' @examples clean_cexpr <- DeleteColsByNames(expr,c("Lars2","Gm10800"))
DeleteColsByNames <- function(DF,col.names.remove){
  stopifnot(is.data.frame(DF))
  stopifnot(col.names.remove %in% colnames(DF))
  DF <- DF[,!(colnames(DF) %in% col.names.remove)]
  return(DF)
}


#' Title Filter out genes with expression level below threshold (thre)
#' if the max expression level of gene is below thre, then drop the gene.
#' @param DF dataframe with each row represents a gene and each column represents a cell
#' @param thre expression levele threshold
#'
#' @return A dataframe
#' @export
#'
#' @examples high_gene_expr_mtx <- FilterGenesExpr(matx,thre=5)
FilterGenesExpr <- function(DF,thre=3){
  stopifnot(is.data.frame(DF))
  stopifnot(thre>0)
  DF[!apply(DF,1,function(x){max(x)<thre}),] # delete row if all items are zero
}
