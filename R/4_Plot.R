# 4. Plot
#   4.1 DimReduction use tsne or umap to a embedding dataframe
#   4.2 Use one gene expression level as color scale
#	4.3 Use category vector level as color map


#' Title Dimension reduction plot (2D) with tsne/umap
#'
#' @param exprmatx Dataframe dataframe with each row represents a gene and each column represents a cell
#' @param method str (default tsne). "umap" / "tsne"
#' @param rdm int. (default 123) random seed number for umap (https://umap-learn.readthedocs.io/en/latest/parameters.html).
#' @param n_neighours int. (default=20) param for umap. Higher, more global structure .
#' @param min_dist float. (default=0.1) between 0-1, Higher, more global structure.
#' @param perplex int. default = 25 tsne para.
#' @param iter int. default = 8000 tsne para.
#'
#' @return a ggplot object
#' @export
#'
#' @examples p <- Embedding(expr,"umap")
Embedding <- function(exprmatx,method="umap",rdm=123,n_neighours=20,min_dist=0.5,perplex=25,iter=8000){
  stopifnot(is.data.frame(exprmatx))
  
  if(method=="tsne"){
    matx <- exprmatx %>% as.matrix() %>% t() %>% Rtsne::normalize_input()
    tsne_out <- Rtsne(matx,dims=2,perplexity =perplex,max_iter=iter,verbose=FALSE,is_distance=FALSE,pca=TRUE,check_duplicates = F) # Run TSNE
    embd <- data.frame(tsne_out$Y)
  }
  else if (method=="umap") {
    matx <- exprmatx %>% as.matrix() %>% t()
    custom.config = umap::umap.defaults
    custom.config$random_state = rdm
    custom.config$n_neighbors = n_neighours
    custom.config$min_dist = min_dist
    embd <- data.frame(umap::umap(matx,random_state=rdm)$layout)
  }
  
  return(embd)
}

#' Title Plot 2D embedding with specific gene expression level as color scale 
#'
#' @param exprmatx Dataframe dataframe with each row represents a gene and each column represents a cell
#' @param embd datafram. By calling Embedding(exprmatx). 
#' @param genename Gene Symbol, eg "Lars2"
#' @param title title of plot
#'
#' @return a ggplot object
#' @export
#'
#' @examples lars2 <- Plot_Embed_Continous(exprmatx,embd,"Lars2")
Plot_Embed_Continous <- function(exprmatx,embd,genename,title=""){
  
  stopifnot(is.data.frame(exprmatx))
  stopifnot(genename!="")
  stopifnot(genename %in% rownames(exprmatx))
  
  library(ggplot2)
  
  if(title==""){title <- paste0(genename," Log expression level")}
  
  embd['co'] <- log(as.numeric(exprmatx[genename,]))
  p <- ggplot2::ggplot(embd)+ggplot2::geom_point(aes(x=X1,y=X2,color=co),size=4)+
    ggplot2::ggtitle(title)+ggplot2::xlab("Embedding 1")+ggplot2::ylab("Embedding 2")+
    theme_bw()+ggplot2::scale_color_distiller(palette = "RdBu")
}


#' Title Plot 2D embedding with factor group as color scale
#'
#' @param exprmatx Dataframe dataframe with each row represents a gene and each column represents a cell
#' @param embd datafram. By calling Embedding(exprmatx). 
#' @param group list. (default="red") corresponding to colnames of exprmatx
#' @param title str. Fig title.
#'
#' @return a ggplot object
#' @export
#'
#' @examples sti_vs_ctrl <- Plot_Embed_Category(exprmatx,embd,roup=c(rep("sti",3),rep("ctrl",5)))
Plot_Embed_Category <- function(exprmatx,embd,group="red",title=""){
  stopifnot(is.data.frame(exprmatx))
  stopifnot(length(group)==ncol(exprmatx))
  
  embd['co'] <- factor(group)
  p <- ggplot2::ggplot(embd)+ggplot2::geom_point(aes(x=X1,y=X2,color=co),size=4)+
    ggplot2::ggtitle(title)+ggplot2::xlab("Embedding 1")+ggplot2::ylab("Embedding 2")+
    theme_bw()+ggsci::scale_colour_simpsons()
}



