#Function for Silhouette scores
compute_silhouette_scores <- function(seurat,res.from=.02, res.to=.4,
                                      by=.01,plot=T, plot.prefix=""){
  
  library(cluster)
  library(ggplot2)
  library(cowplot)
  
  cat(paste0("Using assay: ",DefaultAssay(seurat),"\nUpdate your object
before applying the function if you want to use another one."))
  
  sil.results <- data.frame()
  
  for(resolution in seq(res.from,res.to,by = by)){
    seurat<-FindClusters(seurat,resolution=resolution)
    dists <- dist(Embeddings(seurat,reduction = "umap"))
    combined.fac <- factor(paste0(seurat$seurat_clusters))
    clusters<-as.integer(seurat$seurat_clusters)
    sil <- silhouette(clusters, dist = dists)
    if(plot){
      filename=paste0(plot.prefix,"silhouette.res=",resolution,".png")
      png(filename)
      plot(sil, border = NA)
      dev.off()
      cat("Plotted to: ",filename,"\n")
    }
    sil.results <-
      rbind(sil.results,c(resolution,mean(sil[,3]),length(levels(seurat$seurat_clusters))))
    
  }
  
  colnames(sil.results)<-c("Resolution","Average Si","nClusters")
  
  #plot_grid(
    #ggplot(sil.results,aes(x=factor(Resolution),y=`Average Si`))+geom_bar(stat="identity")+theme_cowplot()+xlab("")+geom_bar(data=subset(sil.results,`Average Si`==max(`Average Si`)), aes(factor(Resolution),`Average Si`),fill="red", stat="identity"),
    #ggplot(sil.results,aes(x=factor(Resolution),y=nClusters))+geom_bar(stat="identity")+theme_cowplot()+xlab("Resolution"),ncol=1)
  #ggsave(paste0(plot.prefix,"silhouette_results.pdf"),w=8,h=6)
  #cat("Combined results plotted to silhouette_results.pdf.\n")
  
  return(sil.results)
}




