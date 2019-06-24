# 4. Plot
#   4.1 DimReduction use tsne or umap to a embedding dataframe
#   4.2 Use one gene expression level as color scale
# 	4.3 Use category vector level as color map
# 5. QC result
#   5.1 Clean/Raw Reads 
#   5.2 Map Ratio / Uniq map ratio / rRNA ratio / Procoding Num / ERCC ratio / mt Ratio
#   5.3 Cell-cell correlation


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
  library(dplyr)
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
  
  library(ggplot2)
  embd['co'] <- factor(group)
  p <- ggplot2::ggplot(embd)+ggplot2::geom_point(aes(x=X1,y=X2,color=co),size=4)+
    ggplot2::ggtitle(title)+ggplot2::xlab("Embedding 1")+ggplot2::ylab("Embedding 2")+
    theme_bw()+ggsci::scale_colour_simpsons()
}

# Rawreads/cleanReads

#' Title Barplot of Rawreads and Cleanreads
#'
#' @param sum Dataframe. Summary table. CellxFeature,c("RawReads","CleanReads"...)
#' @param col list (default col=c(1,2)). Indicates which columns represents c("RawReads","CleanReads")
#' @param srt Bool (default srt=TRUE). Whether sort the bar by Clean Reads. 
#'
#' @return a ggplot object
#' @export
#'
#' @examples sum <- read.table("/Analysis/bioinfo/wenxingzhao/project/21_RAW_time_series/2019_06_Plate1_Blank_0h_24h/summary/QC_summary/Merge_Summary.txt",header = T, row.names = 1,stringsAsFactors = F) \cr p <- Plot_Raw_Clean_Reads(sum = sum,srt = F)
Plot_Raw_Clean_Reads <- function(sum,col=c(1,2),srt=T){
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)

# by default, col1 of sum is Raw reads, col2 is Clean Reads.
Reads <- data.frame(t(data.frame(sum[col[1]]-sum[col[2]],sum[col[2]])))
rownames(Reads) <- c("TrimmedReads","CleanReads")
Reads <- data.frame(t(Reads))
Reads['name'] <- rownames(sum)

if(srt==TRUE){
  level <- arrange(Reads,CleanReads)
} else {  level <- Reads}

Reads <- melt(level)

fill <- c("#5F9EA0", "#E1B378")
p <- ggplot()+geom_bar(data=Reads,aes(x=factor(name,levels = level$name),fill=variable,weight=value),width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values=fill)+ xlab("88 single cells")+ylab("Reads number")+
  scale_y_continuous(labels = comma)

return(p)
}



#' Title Box and dot plot of QC features. 
#'
#' @param sum Dataframe. Summary table. CellxFeature,c("RawReads","CleanReads"...)
#' @param group List. pattern of group identifier (which will be used in group)
#' @param features List. Pattern of selected feature, allow fuzzy search. 
#'
#' @return ggplot2 object list. equal length as features list.
#' @export
#'
#' @examples sum <- read.table("/Analysis/bioinfo/wenxingzhao/project/21_RAW_time_series/2019_06_Plate1_Blank_0h_24h/summary/QC_summary/Merge_Summary.txt",header = T, row.names = 1,stringsAsFactors = F) \cr plot_list <- Plot_QC_Features(sum,group = c("blank","0h","24h"),features = c("Clean","Protein"))
Plot_QC_Features <- function(sum,group,features){
names <- rep(group[1],nrow(sum))
message("Please make sure group is valid!")
for(i in 1:length(group)){
  idx <- grep(group[i],rownames(sum),ignore.case = T)
  stopifnot(length(idx)>0)
  names[idx] <- group[i]
}

message("Please make sure items in features could be matched to colnames of sum!")

plot_list <- list()
for(i in seq_along(features)){
feature <- dplyr::select(sum,grep(features[i],colnames(sum),ignore.case = T))
stopifnot(ncol(feature)>0)
cn <- colnames(feature)
feature <- reshape2::melt(feature)
feature["category"]=names

library(ggplot2)
p <- ggplot(feature,aes(y=value,x=category,fill=category))+geom_boxplot(position=position_dodge(width=0.5),width=0.4,alpha=0.8)+theme_classic()+
  xlab("Experiement design")+ylab(cn)+geom_jitter(position=position_dodge(width=0.5),aes(x=category,y=value,shape=category))+
  scale_fill_brewer(palette = 1)

plot_list[[i]] <- p }

return(plot_list)
}


sum <- read.table("/Analysis/bioinfo/wenxingzhao/project/21_RAW_time_series/2019_06_Plate1_Blank_0h_24h/summary/QC_summary/Merge_Summary.txt",header = T, row.names = 1,stringsAsFactors = F)
plot_list <- Plot_QC_Features(sum,group = c("blank","0h","24h"),features = c("Clean","Protein"))

