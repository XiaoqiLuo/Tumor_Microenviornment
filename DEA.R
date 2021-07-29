rm(list = ls())
args<-commandArgs(T)
#control<-read.csv('./test/control.txt',header = FALSE)[,1]
control<-read.csv(args[1],header = FALSE)[,1]
#treat<-read.csv('./test/treat.txt',header = FALSE)[,1]
treat<-read.csv(args[2],header = FALSE)[,1]
p<-args[3]
needed_packages<-c('TCGAbiolinks','SummarizedExperiment',
                   'edgeR','pheatmap','ggplot2','ggrepel')
if (length(which(!(needed_packages %in% rownames(installed.packages()))))>0){
  print('Please install package(s) listed below:')
  (needed_packages[!(needed_packages %in% rownames(installed.packages()))])
}

library(TCGAbiolinks)
library(SummarizedExperiment)
query.exp <- GDCquery(project = p,
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq",
                      file.type = "normalized_results",
                      experimental.strategy = "RNA-Seq")
GDCdownload(query.exp)
exp<-GDCprepare(query = query.exp, save = TRUE, 
                save.filename = "Exp.rda",remove.files.prepared = TRUE)
query.exp <- GDCquery(project = p,
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq",
                      file.type = "normalized_results",
                      experimental.strategy = "RNA-Seq")

print('Data downloaded!')
##extract partial sample
subexp<-exp[,which(exp$sample %in% c(control,treat))]
dim(subexp)
#grouping
group<-rep("control",dim(subexp)[2])
group[which(subexp$sample %in% treat)]<-'treat'
table(group)
subexp@colData@listData$analyType<-group

# #preparation of DEG analysis
idx <- subexp@colData@listData$analyType %in% c("control")
idx2 <- subexp@colData@listData$analyType %in% c("treat")
dataPrep <- TCGAanalyze_Preprocessing(object = subexp, cor.cut = 0.6)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  qnt.cut = 0.25,
                                  method='quantile')
DEA <- TCGAanalyze_DEA(mat1 = dataFilt[,idx],
                                mat2 = dataFilt[,idx2],
                                Cond1type = "control",
                                Cond2type = "treat",method = "glmLRT")
DEG<-DEA[which((DEA$PValue<0.05) & abs(DEA$logFC)>1),]

colnames(dataFilt)<-substring(colnames(dataFilt),1,16)
dataFilt<-log2(dataFilt+1)
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  write.table(dataFilt,"./result/RNAExp.txt",sep = '\t',quote = FALSE)
  }else{
  dir.create("result")
  write.table(dataFilt,"./result/RNAExp.txt",sep = '\t',quote = FALSE)
  }
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "DEA"))!=0 ){
    write.table(DEG,"./result/DEA/DEG.txt",sep = '\t',quote = FALSE)
  }else{
    dir.create("./result/DEA")
    write.table(DEG,"./result/DEA/DEG.txt",sep = '\t',quote = FALSE)
  }
}else{
  dir.create("result")
  dir.create("./result/DEA")
  write.table(DEG,"./result/DEA/DEG.txt",sep = '\t',quote = FALSE)
}
print('DEA step1 done.')
print('DEG file saved.')
library(pheatmap)
library(ggplot2)
col_anno<-data.frame(type=rep('control',dim(dataFilt)[2]))
rownames(col_anno)<-colnames(dataFilt)
col_anno[which(rownames(col_anno) %in% treat),1]<-'treat'
hm<-pheatmap(log2(dataFilt[which(rownames(dataFilt) %in% rownames(DEG)),]+1),
             show_colnames = FALSE,show_rownames = FALSE,
             scale = 'row',
             color = colorRampPalette(colors = c("blue","white","red"))(100),
             annotation_col = col_anno,cluster_cols = TRUE)
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "DEA"))!=0 ){
    ggsave(hm,filename = "./result/DEA/heatmap.png")
  }else{
    dir.create("DEA")
    ggsave(hm,filename = "./result/DEA/heatmap.png")
  }
}else{
  dir.create("result")
  dir.create("./result/DEA")
  ggsave(hm,filename = "./result/DEA/heatmap.png")
}
print('DEA step2 done.')
print('heatmap saved.')
library(ggplot2)
library(ggrepel)

DEG$threshold = factor(ifelse(DEG$FDR < 0.05 & abs(DEG$logFC) >= 1, ifelse(DEG$logFC>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
volcano<-ggplot(DEG,aes(x=logFC,y=-log10(FDR),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank()
  )+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "DEA"))!=0 ){
    ggsave(volcano,filename = "./result/DEA/volcano.png")
  }else{
    dir.create("DEA")
    ggsave(volcano,filename = "./result/DEA/volcano.png")
  }
}else{
  dir.create("result")
  dir.create("./result/DEA")
  ggsave(volcano,filename = "./result/DEA/volcano.png")
}
print('DEA step3 done.')
print('volcano saved.')

####stemness##############################################################################
st<-TCGAanalyze_Stemness(stemSig = PCBC_stemSig,
                     dataGE = dataFilt, 
                     annotation = "sampleType")

if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "Stemness"))!=0 ){
    write.table(st,"./result/Stemness/stemness.txt",sep = '\t',quote = FALSE,row.names = FALSE)
  }else{
    dir.create("./result/Stemness")
    write.table(st,"./result/Stemness/stemness.txt",sep = '\t',quote = FALSE,row.names = FALSE)
  }
}else{
  dir.create("result")
  dir.create("./result/Stemness")
  write.table(st,"./result/Stemness/stemness.txt",sep = '\t',quote = FALSE,row.names = FALSE)
}
