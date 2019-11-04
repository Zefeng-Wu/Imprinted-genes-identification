library(pheatmap)
library(GenomicFeatures)
tr<-makeTxDbFromGFF("gff/modified.gff")
genes<-genes(tr)

### inspect cluster
data <-read.table("6all_allele_counts.txt",header = TRUE)
pheatmap(cor(data[,5:ncol(data)]))

## annotation each snp with gene name
snp_gr<-makeGRangesFromDataFrame(data,
                                 keep.extra.columns = FALSE,
                                 seqnames.field = "chr",
                                 start.field = "pos",
                                 end.field = "pos")
overlap_with_genes<-findOverlaps(snp_gr,genes)
snp_gr$gene_id<-names(genes[subjectHits(overlap_with_genes)])
data$genes_id <-snp_gr$gene_id

### sum reads for sam genes
library(tidyverse)
data<-as.data.frame(data%>%group_by(genes_id)%>%summarise_at(grep("rep",colnames(data)),sum))

##### plot maternal and paternal bias
bias_plot<-function(df){
  ggplot(df, aes(x=log(df[,1]+df[,4],2), 
                 y=log(df[,2]+df[,3],2))) + 
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    theme(legend.position='none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(size = 20)) +
    #geom_point(aes(shape = rownames(df)[1]),size = 1) + scale_shape(solid = FALSE)+
    ylab("Paternal expression (log2 reads count)")+
    xlab("Maternal expression (log2 reads count)")+
    geom_abline(intercept = 0, slope = 1,col="black",size=1.5,linetype="dashed")+
    geom_abline(intercept = -1, slope = 1,col="black",size=1.5)+
    #geom_abline(intercept = -2, slope = 1,col="orange",size=1.5,linetype="dashed")+
    facet_grid(date~replicate)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))
}
bias_plot_direct1<-function(df){
  ggplot(df, aes(x=log(df[,1],2), 
                 y=log(df[,2],2))) + 
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    theme(legend.position='none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(size = 20)) +
    #geom_point(aes(shape = rownames(df)[1]),size = 1) + scale_shape(solid = FALSE)+
    ylab("Paternal expression (log2 reads count)")+
    xlab("Maternal expression (log2 reads count)")+
    geom_abline(intercept = 0, slope = 1,col="black",size=1.5,linetype="dashed")+
    geom_abline(intercept = -1, slope = 1,col="black",size=1.5)+
    #geom_abline(intercept = -2, slope = 1,col="orange",size=1.5,linetype="dashed")+
    facet_grid(date~replicate)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))
}
bias_plot_direct2<-function(df){
  ggplot(df, aes(x=log(df[,4],2), 
                 y=log(df[,3],2))) + 
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    theme(legend.position='none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(size = 20)) +
    #geom_point(aes(shape = rownames(df)[1]),size = 1) + scale_shape(solid = FALSE)+
    ylab("Paternal expression (log2 reads count)")+
    xlab("Maternal expression (log2 reads count)")+
    geom_abline(intercept = 0, slope = 1,col="black",size=1.5,linetype="dashed")+
    geom_abline(intercept = -1, slope = 1,col="black",size=1.5)+
    #geom_abline(intercept = -2, slope = 1,col="orange",size=1.5,linetype="dashed")+
    facet_grid(date~replicate)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))
}

# 10days ##################################################################################################
##10d repeat1
data_10d_repeat1<-data.frame(X1560_3755_10d_1560=data$X1560_3755_10d_rep1_1560,
                     X1560_3755_10d_3755=data$X1560_3755_10d_rep1_3755,
                     X3755_1560_10d_1560=data$X3755_1560_10d_rep1_1560,
                     X3755_1560_10d_3755=data$X3755_1560_10d_rep1_3755)
rownames(data_10d_repeat1)<-data$genes_id
data_10d_repeat1<-data_10d_repeat1[rowSums(data_10d_repeat1[,c(1,2)])>10&rowSums(data_10d_repeat1[,c(3,4)])>10,] #6792
data_10d_repeat1$p1<-apply(data_10d_repeat1[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_10d_repeat1$p2<-apply(data_10d_repeat1[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_10d_repeat1$adp1<-p.adjust(data_10d_repeat1$p1,method = "fdr")
data_10d_repeat1$adp2<-p.adjust(data_10d_repeat1$p2,method = "fdr")
data_10d_repeat1<-subset(data_10d_repeat1,data_10d_repeat1$adp1<0.01&data_10d_repeat1$adp2<0.01)
data_10d_repeat1$state<-"NA"
data_10d_repeat1$state[data_10d_repeat1$X1560_3755_10d_1560/(data_10d_repeat1$X1560_3755_10d_1560+data_10d_repeat1$X1560_3755_10d_3755)>=0.8&
                       data_10d_repeat1$X3755_1560_10d_3755/(data_10d_repeat1$X3755_1560_10d_1560+data_10d_repeat1$X3755_1560_10d_3755)>=0.8]<-"MEG" #583
data_10d_repeat1$state[data_10d_repeat1$X1560_3755_10d_3755/(data_10d_repeat1$X1560_3755_10d_3755+data_10d_repeat1$X1560_3755_10d_1560)>=0.7&
                       data_10d_repeat1$X3755_1560_10d_1560/(data_10d_repeat1$X3755_1560_10d_1560+data_10d_repeat1$X3755_1560_10d_3755)>=0.7]<-"PEG"


## 10d repeat2
data_10d_repeat2<-data.frame(X1560_3755_10d_1560=data$X1560_3755_10d_rep2_1560,
                             X1560_3755_10d_3755=data$X1560_3755_10d_rep2_3755,
                             X3755_1560_10d_1560=data$X3755_1560_10d_rep2_1560,
                             X3755_1560_10d_3755=data$X3755_1560_10d_rep2_3755)
rownames(data_10d_repeat2)<-data$genes_id
data_10d_repeat2<-data_10d_repeat2[rowSums(data_10d_repeat2[,c(1,2)])>10&rowSums(data_10d_repeat2[,c(3,4)])>10,] #6510
data_10d_repeat2$p1<-apply(data_10d_repeat2[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_10d_repeat2$p2<-apply(data_10d_repeat2[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_10d_repeat2$adp1<-p.adjust(data_10d_repeat2$p1,method = "fdr")
data_10d_repeat2$adp2<-p.adjust(data_10d_repeat2$p2,method = "fdr")
data_10d_repeat2<-subset(data_10d_repeat2,data_10d_repeat2$adp1<0.01&data_10d_repeat2$adp2<0.01)
data_10d_repeat2$state<-"NA"
data_10d_repeat2$state[data_10d_repeat2$X1560_3755_10d_1560/(data_10d_repeat2$X1560_3755_10d_1560+data_10d_repeat2$X1560_3755_10d_3755)>=0.8&
                       data_10d_repeat2$X3755_1560_10d_3755/(data_10d_repeat2$X3755_1560_10d_1560+data_10d_repeat2$X3755_1560_10d_3755)>=0.8]<-"MEG" #583
data_10d_repeat2$state[data_10d_repeat2$X1560_3755_10d_3755/(data_10d_repeat2$X1560_3755_10d_3755+data_10d_repeat2$X1560_3755_10d_1560)>=0.7&
                       data_10d_repeat2$X3755_1560_10d_1560/(data_10d_repeat2$X3755_1560_10d_1560+data_10d_repeat2$X3755_1560_10d_3755)>=0.7]<-"PEG"

##10d repeat3
data_10d_repeat3<-data.frame(X1560_3755_10d_1560=data$X1560_3755_10d_rep3_1560,
                             X1560_3755_10d_3755=data$X1560_3755_10d_rep3_3755,
                             X3755_1560_10d_1560=data$X3755_1560_10d_rep3_1560,
                             X3755_1560_10d_3755=data$X3755_1560_10d_rep3_3755)
rownames(data_10d_repeat3)<-data$genes_id
data_10d_repeat3<-data_10d_repeat3[rowSums(data_10d_repeat3[,c(1,2)])>10&rowSums(data_10d_repeat3[,c(3,4)])>10,] #6400
data_10d_repeat3$p1<-apply(data_10d_repeat3[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_10d_repeat3$p2<-apply(data_10d_repeat3[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_10d_repeat3$adp1<-p.adjust(data_10d_repeat3$p1,method = "fdr")
data_10d_repeat3$adp2<-p.adjust(data_10d_repeat3$p2,method = "fdr")
data_10d_repeat3<-subset(data_10d_repeat3,data_10d_repeat3$adp1<0.01&data_10d_repeat3$adp2<0.01)
data_10d_repeat3$state<-"NA"
data_10d_repeat3$state[data_10d_repeat3$X1560_3755_10d_1560/(data_10d_repeat3$X1560_3755_10d_1560+data_10d_repeat3$X1560_3755_10d_3755)>=0.8&
                       data_10d_repeat3$X3755_1560_10d_3755/(data_10d_repeat3$X3755_1560_10d_1560+data_10d_repeat3$X3755_1560_10d_3755)>=0.8]<-"MEG" #583
data_10d_repeat3$state[data_10d_repeat3$X1560_3755_10d_3755/(data_10d_repeat3$X1560_3755_10d_3755+data_10d_repeat3$X1560_3755_10d_1560)>=0.7&
                       data_10d_repeat3$X3755_1560_10d_1560/(data_10d_repeat3$X3755_1560_10d_1560+data_10d_repeat3$X3755_1560_10d_3755)>=0.7]<-"PEG"


# 13days #################################################################################################
##13d repeat1
data_13d_repeat1<-data.frame(X1560_3755_13d_1560=data$X1560_3755_13d_rep1_1560,
                             X1560_3755_13d_3755=data$X1560_3755_13d_rep1_3755,
                             X3755_1560_13d_1560=data$X3755_1560_13d_rep1_1560,
                             X3755_1560_13d_3755=data$X3755_1560_13d_rep1_3755)
rownames(data_13d_repeat1)<-data$genes_id
data_13d_repeat1<-data_13d_repeat1[rowSums(data_13d_repeat1[,c(1,2)])>10&rowSums(data_13d_repeat1[,c(3,4)])>10,] #6732
data_13d_repeat1$p1<-apply(data_13d_repeat1[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_13d_repeat1$p2<-apply(data_13d_repeat1[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_13d_repeat1$adp1<-p.adjust(data_13d_repeat1$p1,method = "fdr")
data_13d_repeat1$adp2<-p.adjust(data_13d_repeat1$p2,method = "fdr")
data_13d_repeat1<-subset(data_13d_repeat1,data_13d_repeat1$adp1<0.01&data_13d_repeat1$adp2<0.01)
data_13d_repeat1$state<-"NA"
data_13d_repeat1$state[data_13d_repeat1$X1560_3755_13d_1560/(data_13d_repeat1$X1560_3755_13d_1560+data_13d_repeat1$X1560_3755_13d_3755)>=0.8&
                         data_13d_repeat1$X3755_1560_13d_3755/(data_13d_repeat1$X3755_1560_13d_1560+data_13d_repeat1$X3755_1560_13d_3755)>=0.8]<-"MEG" #583
data_13d_repeat1$state[data_13d_repeat1$X1560_3755_13d_3755/(data_13d_repeat1$X1560_3755_13d_3755+data_13d_repeat1$X1560_3755_13d_1560)>=0.7&
                         data_13d_repeat1$X3755_1560_13d_1560/(data_13d_repeat1$X3755_1560_13d_1560+data_13d_repeat1$X3755_1560_13d_3755)>=0.7]<-"PEG"


## 13d repeat2
data_13d_repeat2<-data.frame(X1560_3755_13d_1560=data$X1560_3755_13d_rep2_1560,
                             X1560_3755_13d_3755=data$X1560_3755_13d_rep2_3755,
                             X3755_1560_13d_1560=data$X3755_1560_13d_rep2_1560,
                             X3755_1560_13d_3755=data$X3755_1560_13d_rep2_3755)
rownames(data_13d_repeat2)<-data$genes_id
data_13d_repeat2<-data_13d_repeat2[rowSums(data_13d_repeat2[,c(1,2)])>10&rowSums(data_13d_repeat2[,c(3,4)])>10,] #6833
data_13d_repeat2$p1<-apply(data_13d_repeat2[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_13d_repeat2$p2<-apply(data_13d_repeat2[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_13d_repeat2$adp1<-p.adjust(data_13d_repeat2$p1,method = "fdr")
data_13d_repeat2$adp2<-p.adjust(data_13d_repeat2$p2,method = "fdr")
data_13d_repeat2<-subset(data_13d_repeat2,data_13d_repeat2$adp1<0.01&data_13d_repeat2$adp2<0.01)
data_13d_repeat2$state<-"NA"
data_13d_repeat2$state[data_13d_repeat2$X1560_3755_13d_1560/(data_13d_repeat2$X1560_3755_13d_1560+data_13d_repeat2$X1560_3755_13d_3755)>=0.8&
                         data_13d_repeat2$X3755_1560_13d_3755/(data_13d_repeat2$X3755_1560_13d_1560+data_13d_repeat2$X3755_1560_13d_3755)>=0.8]<-"MEG" #583
data_13d_repeat2$state[data_13d_repeat2$X1560_3755_13d_3755/(data_13d_repeat2$X1560_3755_13d_3755+data_13d_repeat2$X1560_3755_13d_1560)>=0.7&
                         data_13d_repeat2$X3755_1560_13d_1560/(data_13d_repeat2$X3755_1560_13d_1560+data_13d_repeat2$X3755_1560_13d_3755)>=0.7]<-"PEG"

##13d repeat3
data_13d_repeat3<-data.frame(X1560_3755_13d_1560=data$X1560_3755_13d_rep3_1560,
                             X1560_3755_13d_3755=data$X1560_3755_13d_rep3_3755,
                             X3755_1560_13d_1560=data$X3755_1560_13d.rep3_1560,
                             X3755_1560_13d_3755=data$X3755_1560_13d.rep3_3755)
rownames(data_13d_repeat3)<-data$genes_id
data_13d_repeat3<-data_13d_repeat3[rowSums(data_13d_repeat3[,c(1,2)])>10&rowSums(data_13d_repeat3[,c(3,4)])>10,] #6792
data_13d_repeat3$p1<-apply(data_13d_repeat3[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_13d_repeat3$p2<-apply(data_13d_repeat3[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_13d_repeat3$adp1<-p.adjust(data_13d_repeat3$p1,method = "fdr")
data_13d_repeat3$adp2<-p.adjust(data_13d_repeat3$p2,method = "fdr")
data_13d_repeat3<-subset(data_13d_repeat3,data_13d_repeat3$adp1<0.01&data_13d_repeat3$adp2<0.01)
data_13d_repeat3$state<-"NA"
data_13d_repeat3$state[data_13d_repeat3$X1560_3755_13d_1560/(data_13d_repeat3$X1560_3755_13d_1560+data_13d_repeat3$X1560_3755_13d_3755)>=0.8&
                         data_13d_repeat3$X3755_1560_13d_3755/(data_13d_repeat3$X3755_1560_13d_1560+data_13d_repeat3$X3755_1560_13d_3755)>=0.8]<-"MEG" #583
data_13d_repeat3$state[data_13d_repeat3$X1560_3755_13d_3755/(data_13d_repeat3$X1560_3755_13d_3755+data_13d_repeat3$X1560_3755_13d_1560)>=0.7&
                         data_13d_repeat3$X3755_1560_13d_1560/(data_13d_repeat3$X3755_1560_13d_1560+data_13d_repeat3$X3755_1560_13d_3755)>=0.7]<-"PEG"


# 15days #################################################################################################
# 15 Days repeat1
data_15d_repeat1<-data.frame(X1560_3755_15d_1560=data$X1560_3755_15d_rep1_1560,
                             X1560_3755_15d_3755=data$X1560_3755_15d_rep1_3755,
                             X3755_1560_15d_1560=data$X3755_1560_15d_rep1_1560,
                             X3755_1560_15d_3755=data$X3755_1560_15d_rep1_3755)
rownames(data_15d_repeat1)<-data$genes_id
data_15d_repeat1<-data_15d_repeat1[rowSums(data_15d_repeat1[,c(1,2)])>10&rowSums(data_15d_repeat1[,c(3,4)])>10,] #6867
data_15d_repeat1$p1<-apply(data_15d_repeat1[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_15d_repeat1$p2<-apply(data_15d_repeat1[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_15d_repeat1$adp1<-p.adjust(data_15d_repeat1$p1,method = "fdr")
data_15d_repeat1$adp2<-p.adjust(data_15d_repeat1$p2,method = "fdr")
data_15d_repeat1<-subset(data_15d_repeat1,data_15d_repeat1$adp1<0.01&data_15d_repeat1$adp2<0.01)
data_15d_repeat1$state<-"NA"
data_15d_repeat1$state[data_15d_repeat1$X1560_3755_15d_1560/(data_15d_repeat1$X1560_3755_15d_1560+data_15d_repeat1$X1560_3755_15d_3755)>=0.8&
                         data_15d_repeat1$X3755_1560_15d_3755/(data_15d_repeat1$X3755_1560_15d_1560+data_15d_repeat1$X3755_1560_15d_3755)>=0.8]<-"MEG" #583
data_15d_repeat1$state[data_15d_repeat1$X1560_3755_15d_3755/(data_15d_repeat1$X1560_3755_15d_3755+data_15d_repeat1$X1560_3755_15d_1560)>=0.7&
                         data_15d_repeat1$X3755_1560_15d_1560/(data_15d_repeat1$X3755_1560_15d_1560+data_15d_repeat1$X3755_1560_15d_3755)>=0.7]<-"PEG"


## 15d repeat2
data_15d_repeat2<-data.frame(X1560_3755_15d_1560=data$X1560_3755_15d_rep2_1560,
                             X1560_3755_15d_3755=data$X1560_3755_15d_rep2_3755,
                             X3755_1560_15d_1560=data$X3755_1560_15d_rep2_1560,
                             X3755_1560_15d_3755=data$X3755_1560_15d_rep2_3755)
rownames(data_15d_repeat2)<-data$genes_id
data_15d_repeat2<-data_15d_repeat2[rowSums(data_15d_repeat2[,c(1,2)])>10&rowSums(data_15d_repeat2[,c(3,4)])>10,] #6792
data_15d_repeat2$p1<-apply(data_15d_repeat2[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_15d_repeat2$p2<-apply(data_15d_repeat2[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_15d_repeat2$adp1<-p.adjust(data_15d_repeat2$p1,method = "fdr")
data_15d_repeat2$adp2<-p.adjust(data_15d_repeat2$p2,method = "fdr")
data_15d_repeat2<-subset(data_15d_repeat2,data_15d_repeat2$adp1<0.01&data_15d_repeat2$adp2<0.01)
data_15d_repeat2$state<-"NA"
data_15d_repeat2$state[data_15d_repeat2$X1560_3755_15d_1560/(data_15d_repeat2$X1560_3755_15d_1560+data_15d_repeat2$X1560_3755_15d_3755)>=0.8&
                         data_15d_repeat2$X3755_1560_15d_3755/(data_15d_repeat2$X3755_1560_15d_1560+data_15d_repeat2$X3755_1560_15d_3755)>=0.8]<-"MEG" #583
data_15d_repeat2$state[data_15d_repeat2$X1560_3755_15d_3755/(data_15d_repeat2$X1560_3755_15d_3755+data_15d_repeat2$X1560_3755_15d_1560)>=0.7&
                         data_15d_repeat2$X3755_1560_15d_1560/(data_15d_repeat2$X3755_1560_15d_1560+data_15d_repeat2$X3755_1560_15d_3755)>=0.7]<-"PEG"

##15d repeat3
data_15d_repeat3<-data.frame(X1560_3755_15d_1560=data$X1560_3755_15d_rep3_1560,
                             X1560_3755_15d_3755=data$X1560_3755_15d_rep3_3755,
                             X3755_1560_15d_1560=data$X3755_1560_15d_rep3_1560,
                             X3755_1560_15d_3755=data$X3755_1560_15d_rep3_3755)
rownames(data_15d_repeat3)<-data$genes_id
data_15d_repeat3<-data_15d_repeat3[rowSums(data_15d_repeat3[,c(1,2)])>10&rowSums(data_15d_repeat3[,c(3,4)])>10,] #6875
data_15d_repeat3$p1<-apply(data_15d_repeat3[,1:4],1,function(x)chisq.test(c(x[1],x[2]),p = c(2/3,1/3))$p.value)
data_15d_repeat3$p2<-apply(data_15d_repeat3[,1:4],1,function(x)chisq.test(c(x[4],x[3]),p = c(2/3,1/3))$p.value)
data_15d_repeat3$adp1<-p.adjust(data_15d_repeat3$p1,method = "fdr")
data_15d_repeat3$adp2<-p.adjust(data_15d_repeat3$p2,method = "fdr")
data_15d_repeat3<-subset(data_15d_repeat3,data_15d_repeat3$adp1<0.01&data_15d_repeat3$adp2<0.01)
data_15d_repeat3$state<-"NA"
data_15d_repeat3$state[data_15d_repeat3$X1560_3755_15d_1560/(data_15d_repeat3$X1560_3755_15d_1560+data_15d_repeat3$X1560_3755_15d_3755)>=0.8&
                         data_15d_repeat3$X3755_1560_15d_3755/(data_15d_repeat3$X3755_1560_15d_1560+data_15d_repeat3$X3755_1560_15d_3755)>=0.8]<-"MEG" #583
data_15d_repeat3$state[data_15d_repeat3$X1560_3755_15d_3755/(data_15d_repeat3$X1560_3755_15d_3755+data_15d_repeat3$X1560_3755_15d_1560)>=0.7&
                         data_15d_repeat3$X3755_1560_15d_1560/(data_15d_repeat3$X3755_1560_15d_1560+data_15d_repeat3$X3755_1560_15d_3755)>=0.7]<-"PEG"



##### Bias dsitribution plot ###########################################################################################
library(ggpubr)

data_10d_repeat1$date<-"10DAP"
data_10d_repeat1$replicate<-"replicate1"
data_10d_repeat2$date<-"10DAP"
data_10d_repeat2$replicate<-"replicate2"
data_10d_repeat3$date<-"10DAP"
data_10d_repeat3$replicate<-"replicate3"

data_13d_repeat1$date<-"13DAP"
data_13d_repeat1$replicate<-"replicate1"
data_13d_repeat2$date<-"13DAP"
data_13d_repeat2$replicate<-"replicate2"
data_13d_repeat3$date<-"13DAP"
data_13d_repeat3$replicate<-"replicate3"

data_15d_repeat1$date<-"15DAP"
data_15d_repeat1$replicate<-"replicate1"
data_15d_repeat2$date<-"15DAP"
data_15d_repeat2$replicate<-"replicate2"
data_15d_repeat3$date<-"15DAP"
data_15d_repeat3$replicate<-"replicate3"

library(data.table) 
df_plot<-rbindlist(list(data_10d_repeat1,data_10d_repeat2,data_10d_repeat3,
               data_13d_repeat1,data_13d_repeat2,data_13d_repeat3,
               data_15d_repeat1,data_15d_repeat2,data_15d_repeat3),
               use.names = FALSE)
bias_plot(as.data.frame(df_plot))


##### Imprinted genes common################################################# 
common_MEGs_10d<-Reduce(intersect, list(rownames(data_10d_repeat1)[data_10d_repeat1$state=="MEG"],
                                        rownames(data_10d_repeat2)[data_10d_repeat2$state=="MEG"],
                                        rownames(data_10d_repeat3)[data_10d_repeat3$state=="MEG"]))
common_PEGs_10d<-Reduce(intersect, list(rownames(data_10d_repeat1)[data_10d_repeat1$state=="PEG"],
                                        rownames(data_10d_repeat2)[data_10d_repeat2$state=="PEG"],
                                        rownames(data_10d_repeat3)[data_10d_repeat3$state=="PEG"]))

common_MEGs_13d<-Reduce(intersect, list(rownames(data_13d_repeat1)[data_13d_repeat1$state=="MEG"],
                                        rownames(data_13d_repeat2)[data_13d_repeat2$state=="MEG"],
                                        rownames(data_13d_repeat3)[data_13d_repeat3$state=="MEG"]))
common_PEGs_13d<-Reduce(intersect, list(rownames(data_13d_repeat1)[data_13d_repeat1$state=="PEG"],
                                        rownames(data_13d_repeat2)[data_13d_repeat2$state=="PEG"],
                                        rownames(data_13d_repeat3)[data_13d_repeat3$state=="PEG"]))

common_MEGs_15d<-Reduce(intersect, list(rownames(data_15d_repeat1)[data_15d_repeat1$state=="MEG"],
                                        rownames(data_15d_repeat2)[data_15d_repeat2$state=="MEG"],
                                        rownames(data_15d_repeat3)[data_15d_repeat3$state=="MEG"]))
common_PEGs_15d<-Reduce(intersect, list(rownames(data_15d_repeat1)[data_15d_repeat1$state=="PEG"],
                                        rownames(data_15d_repeat2)[data_15d_repeat2$state=="PEG"],
                                        rownames(data_15d_repeat3)[data_15d_repeat3$state=="PEG"]))


## PLOT veen
library(ggVennDiagram)
MEGs <- list(A=common_MEGs_10d,B=common_MEGs_13d,C=common_MEGs_15d)
MEGs_plot<-ggVennDiagram(x = MEGs,
                         category.names = c("10DAP","13DAP","15DAP"),
                         lty="dashed",
                         color="black",
                         size=2) +
  scale_fill_gradient(low="white",high = "red")+
  ggtitle("Maternally expressed genes")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))


PEGs <- list(A=common_PEGs_10d,B=common_PEGs_13d,C=common_PEGs_15d)
PEGs_plot<-ggVennDiagram(x = PEGs,
                         category.names = c("10DAP","13DAP","15DAP"),
                         lty="dashed",
                         color="black",
                         size=2) +
  scale_fill_gradient(low="white",high = "red")+
  ggtitle("Paternally expressed genes")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
library(ggpubr)
ggarrange(MEGs_plot,PEGs_plot,labels = 'AUTO')
