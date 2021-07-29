####input##################################################################
args<-commandArgs(T)
i<-args[1]
control<-read.csv(args[2],header = FALSE)[,1]
treat<-read.csv(args[3],header = FALSE)[,1]

####follow up#############################################################
download.file('https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp','followup.txt')
followup<-read.csv('followup.txt',sep = '\t')
####packages##############################################################
needed_packages<-c('ggsignif','ggplot2','estimate','xCell','patchwork','survival','survminer')
if (length(which(!(needed_packages %in% rownames(installed.packages()))))>0){
  print('Please install package(s) listed below:')
  (needed_packages[!(needed_packages %in% rownames(installed.packages()))])
}

####estimate################################################################
library(estimate)
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "Immune"))!=0 ){
    filterCommonGenes(input.f=i, 
                      output.f="./result/Immune/genes.gct", 
                      id="GeneSymbol")
  }else{
    dir.create("./result/Immune")
    filterCommonGenes(input.f=i, 
                      output.f="./result/Immune/genes.gct", 
                      id="GeneSymbol")
  }
}else{
  dir.create("result")
  dir.create("./result/Immune")
  filterCommonGenes(input.f=i, 
                    output.f="./result/Immune/genes.gct", 
                    id="GeneSymbol")
}

estimateScore(input.ds = "./result/Immune/genes.gct",
              output.ds="./result/Immune/estimate_score.gct", 
              platform="illumina")
unlink('./result/Immune/genes.gct')
scores=read.table("./result/Immune/estimate_score.gct",skip = 2,header = T)
rownames(scores)<-scores$NAME
scores<-scores[,-c(1,2)]
colnames(scores)<-gsub('[.]','-',colnames(scores))
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "Immune"))!=0 ){
    write.table(scores,'./result/Immune/estimate_score.txt',sep = '\t',quote = FALSE)
    unlink('./result/Immune/estimate_score.gct')
  }else{
    dir.create("./result/Immune")
    write.table(scores,'./result/Immune/estimate_score.txt',sep = '\t',quote = FALSE)
    unlink('./result/Immune/estimate_score.gct')
  }
}else{
  dir.create("result")
  dir.create("./result/Immune")
  write.table(scores,'./result/Immune/estimate_score.txt',sep = '\t',quote = FALSE)
  unlink('./result/Immune/estimate_score.gct')
}

####estimate plot############################################################
library(ggplot2)
library(ggsignif)
scores<-as.data.frame(t(scores))
scores$group<-rep('control',dim(scores)[1])
scores$group[which(rownames(scores) %in% treat)]<-'treat'
compaired <- list(c('control','treat'))

vioGroup<-function(df,xl,yl,fl){
  ggplot(df,aes_string(x=xl,y=yl,fill=fl))+
    geom_violin(position = position_dodge(), 
                alpha=0.8,width=1) + 
    guides(fill=F)+xlab('')+
    #facet_wrap(~ ind, scales="free")+
    theme_bw()+
    theme(
      legend.title = element_blank(),
      panel.grid.major=element_blank(),panel.grid.minor=element_blank()
    )+#scale_fill_brewer(palette = "PRGn")+
    theme(legend.position = "bottom")+
    #geom_boxplot(position = position_dodge(width = 10))+
    geom_signif(comparisons = compaired,
                map_signif_level = F,test = t.test)
}
stromal<-vioGroup(scores,'group','StromalScore','group')
immune<-vioGroup(scores,'group','ImmuneScore','group')
escore<-vioGroup(scores,'group','ESTIMATEScore','group')
library(patchwork)
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "Immune"))!=0 ){
    ggsave(stromal+immune+escore,filename = "./result/Immune/estimate-ttest.png")
  }else{
    dir.create("Immune")
    ggsave(stromal+immune+escore,filename = "./result/Immune/estimate-ttest.png")
  }
}else{
  dir.create("result")
  dir.create("./result/Immune")
  ggsave(stromal+immune+escore,filename = "./result/Immune/estimate-ttest.png")
}
####estimate survival#################################################################
library(survival)
require(survminer)
surv<-followup[,c(1,26:33)]
scores$bcr<-substring(rownames(scores),1,15)
scores<-merge(scores,surv,by.x='bcr',by.y='sample')
#args[4]='mean'
if(args[4]=='mean'){
  scores$stromalScore<-ifelse(scores$StromalScore<=mean(scores$StromalScore),'L','H')
  scores$immuneScore<-ifelse(scores$ImmuneScore<=mean(scores$ImmuneScore),'L','H')
  scores$estimateScore<-ifelse(scores$ESTIMATEScore<=mean(scores$ESTIMATEScore),'L','H')
}else if(args[4]=='median'){
  scores$stromalScore<-ifelse(scores$StromalScore<=median(scores$StromalScore),'L','H')
  scores$immuneScore<-ifelse(scores$ImmuneScore<=median(scores$ImmuneScore),'L','H')
  scores$estimateScore<-ifelse(scores$ESTIMATEScore<=median(scores$ESTIMATEScore),'L','H')
}
formulae <- list()
element<-c("stromalScore","immuneScore","estimateScore")
formulae <- lapply(element, function(x) as.formula(paste0('Surv(',args[5],'.time',',',args[5],')', " ~",x)))
fits <- surv_fit(formulae, data = scores)
p <- ggsurvplot_list(fits,title=c('      StromalScore','      ImmuneScore','      EstimateScore'),
                     data = scores,
                     risk.table = FALSE,
                     pval = TRUE,ylab = "Recurrence Free Probability",
                     legend.title = ""
)
estSur<-arrange_ggsurvplots(p[1:length(p)], print = TRUE, ncol = 2, nrow = 2)
estSur
if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "Immune"))!=0 ){
    ggsave(estSur,filename = "./result/Immune/estimate-survival.png")
  }else{
    dir.create("Immune")
    ggsave(estSur,filename = "./result/Immune/estimate-survival.png")
  }
}else{
  dir.create("result")
  dir.create("./result/Immune")
  ggsave(estSur,filename = "./result/Immune/estimate-survival.png")
}
####xCell#######################################################################
#devtools::install_github('dviraran/xCell')
library(xCell)
library(survival)
library(survminer)
xCell<-xCellAnalysis(read.csv(i,sep = '\t',check.names = FALSE), rnaseq = TRUE)
xCell<-as.data.frame(t(xCell))
plist<-vector()
cellttest<-as.data.frame(matrix(ncol = 4,nrow = 64))
colnames(cellttest)<-c('CellName','ControlMean','TreatMean','PValue')
for (j in 1:(ncol(xCell)-3)) {
  cellttest[j,1]<-colnames(xCell)[j]
  cellttest[j,2:3]<-t.test(xCell[which(rownames(xCell) %in% control),j],xCell[which(rownames(xCell) %in% treat),j])$estimate
  cellttest[j,4]<-t.test(xCell[which(rownames(xCell) %in% control),j],xCell[which(rownames(xCell) %in% treat),j])$p.value
}

if(length(which(list.files(full.names = TRUE) %in% "./result"))!=0){
  if(length(which(list.files(paste(getwd(),"./result",sep = '')) %in% "Immune"))!=0 ){
    write.table(cellttest,'./result/Immune/xCell_ttest.txt',sep = '\t',quote = FALSE,row.names = FALSE)
  }else{
    dir.create("./result/Immune")
    write.table(cellttest,'./result/Immune/xCell_ttest.txt',sep = '\t',quote = FALSE,row.names = FALSE)
  }
}else{
  dir.create("result")
  dir.create("./result/Immune")
  write.table(cellttest,'./result/Immune/xCell_ttest.txt',sep = '\t',quote = FALSE,row.names = FALSE)
}
vCell<-stack(xCell)
vCell$sample<-rep(rownames(xCell),67)
vCell$group<-'control'
vCell$group[which(vCell$sample %in% treat)]<-'treat'

if(length(which(list.files(paste(getwd(),"./result/Immune",sep = '')) %in% 'Immune'))==0){
  dir.create("./result/Immune/xCell-ttest")
  print("created file.")
}
for(k in which(cellttest$PValue<0.05)){
  v<-vioGroup(vCell[which(vCell$ind==cellttest$CellName[k]),],'group','values','group')
  ggsave(v,filename = paste0("./result/Immune/xCell-ttest/",cellttest$CellName[k],"-ttest.png"))
}

#survival
xCell$bcr<-substring(rownames(xCell),1,15)
xCell<-merge(xCell,surv,by.x='bcr',by.y='sample')
rownames(xCell)<-xCell$bcr
xCell<-xCell[,-1]
formulae <- list()
cellname<-colnames(xCell)
colnames(xCell)<-gsub(' ','',colnames(xCell))
colnames(xCell)<-gsub('-','',colnames(xCell))
colnames(xCell)<-gsub('[+]','',colnames(xCell))
element<-colnames(xCell)[1:64]
a<-xCell
for (m in 1:(ncol(xCell)-11)) {
  n1<-mean(xCell[,m])
  n2<-median(xCell[,m])
  for(h in 1:nrow(xCell)){
    if(args[4]=='mean'){
      xCell[h,m]<-ifelse(xCell[h,m]<=n1,'L','H')
    }else{
      n<-median(xCell[,m])
      xCell[h,m]<-ifelse(xCell[h,m]<=n2,'L','H')
    }
  }
}

formulae <- lapply(element, function(x) as.formula(paste0('Surv(',args[5],'.time',',',args[5],')', " ~",x)))
fits <- surv_fit(formulae, data = xCell)
p <- ggsurvplot_list(fits,title=cellname[1:64],
                     data = xCell,
                     risk.table = FALSE,
                     pval = TRUE,ylab = "Recurrence Free Probability",
                     #break.time.by = 500,#legend = c(0.85,0.25),
                     #ggtheme = theme_minimal(),
                     legend.title = ""
)
estSur<-arrange_ggsurvplots(p[1:length(p)], print = TRUE, ncol = 2, nrow = 2)
pdf('./result/Immune/cell-survival.pdf')
estSur
dev.off()

####expr##############################################################################
deg<-strsplit(args[6], ",")[[1]]
expr<-read.csv(i,sep = '\t',check.names = FALSE)
if(length(deg[!(deg %in% rownames(expr))])!=0){
  print(paste(deg[!(deg %in% rownames(expr))],'is not included in expression matrix!'))
}
deg<-deg[which(deg %in% rownames(expr))]
degcor<-expr[which(rownames(expr) %in% deg),]
degcor<-as.data.frame(t(degcor))
xCell<-xCellAnalysis(read.csv(i,sep = '\t',check.names = FALSE), rnaseq = TRUE)
xCell<-as.data.frame(t(xCell))
degcor<-merge(degcor,xCell,by='row.names')
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
for (k in 1:length(deg)) {
  if(length(which(list.files(paste(getwd(),"./result/Immune",sep = ''))==deg[k]))==0){
    dir.create(paste0('./result/Immune/',deg[k]))
  }
  for (j in (length(deg)+2):(ncol(degcor)-3)) {
    ci<-ggplot(degcor, aes(x = degcor[,which(colnames(degcor)==deg[k])], y = degcor[,j])) + 
      geom_point()+
      geom_smooth(method = "lm", color = "black", fill = "lightgray")+
      labs(x = colnames(degcor)[k], y = colnames(degcor)[j])+
           stat_cor(method = "pearson", 
           label.x = 3, label.y = 30)
    ggsave(ci,filename = paste0('./result/Immune/',deg[k],'/',deg[k],'-',colnames(degcor)[j],'.png'))
  }
}






