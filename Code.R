### Author : Qingyuan Sun, Shandong University. 
### Email : sunqingyuansdu@163.com.


##### 1 Nomalization of bulk RNA sequceing data #####


library(limma)
setwd("")

rt1=read.table("id trans gtex brain.txt",sep="\t",header=T,check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)

rt2=read.table("symbol.txt",sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)

sameGene=intersect( row.names(data1),row.names(data2) )
data=cbind(data1[sameGene,],data2[sameGene,])

outTab=normalizeBetweenArrays(data)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="merge.txt",sep="\t",quote=F,col.names=F)





##### 2 Differential expression analysis #####


library("limma")

setwd("")                 
inputFile="input.txt"                                           
fdrFilter=0.05                                                  
logFCfilter=1                                                    
conNum=1157                                                      
treatNum=697                                                  

outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)


##### 3 Uni-COX #####

library(survival)
pFilter=0.05                                                     
setwd("")       
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1,)  

outTab=data.frame()
sigGenes=c("futime","fustat")
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)




##### 4 Multi-COX #####


library(survival)   
library(survminer)

setwd("")
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)  
rt$futime=rt$futime/365

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)

pdf("forest.pdf",width = 5,height = 5)
plot<-ggforest(multiCox)
plot
dev.off()


##### 5 K-M Curve #####

setwd("")   
library(survival)
library(survminer)
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

pdf(file="survival.pdf",onefile = FALSE,
    width = 4.5,           
    height =4.5)          
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p<0.0001"),
           pval.size=4,
           risk.table=F,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()

##### 6 Processing of single cell sequceing data #####

##### 1 read
data<-Read10X_h5("Glioma_GSE131928_10X_expression.h5",use.names = T)
meta<-read_tsv("Glioma_GSE131928_10X_CellMetainfo_table.tsv")

data<-as.data.frame(data)
pbmc=CreateSeuratObject(counts = data,project = "seurat")
pbmc$TISCH_cluster<-meta$Cluster
pbmc$malignancy<-meta$`Celltype (malignancy)`
pbmc$cell_type_major<-meta$`Celltype (major-lineage)`
pbmc$cell_type_minor<-meta$`Celltype (minor-lineage)`
pbmc$patient<-meta$Patient
pbmc$sample<-meta$Sample
pbmc$age<-meta$Age
pbmc$gender<-meta$Gender
saveRDS(pbmc,"GSE131928.rds")
rm(data,meta)
sqy<-pbmc
rm(pbmc)
gc()


#### 2 dimension reduction
library(devtools)
library(harmony)
gc()

sqy[["percent.mt"]] <- PercentageFeatureSet(object = sqy, pattern = "^MT-")

VlnPlot(object = sqy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident") 
sqy=subset(x = sqy, subset = nFeature_RNA > 50 & percent.mt < 10) 

plot1 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
rm(plot1,plot2)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sqy <- CellCycleScoring(sqy, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sqy@meta.data[1:5,]
sqy<-NormalizeData(sqy,verbose = T)   

sqy<-FindVariableFeatures(sqy,selection.method = "vst", nfeatures = 2000)   
sqy<-ScaleData(sqy,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = T)
sqy<-RunPCA(sqy,verbose = T,npcs = 70)
ElbowPlot(sqy,ndims = 70) 

sqy <- FindNeighbors(object = sqy, dims = 1:30) 
sqy <- FindClusters(object = sqy, resolution = 1) 
sqy <- RunUMAP(object = sqy, dims = 1:30) 

pdf(file = "UMAP_cluster.pdf",width=5,height = 4)
UMAPPlot(object = sqy, label = TRUE, cols=cors)
dev.off()

save.image("24.1.14.RData")

### 3 pathway score
library("KEGGREST") 
listDatabases()  
gs<-keggGet('hsa04140')
gs[[1]]$GENE 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
pathways <- genes[1:length(genes)%%3 ==2] 
pathways <- data.frame(pathways)  
rm(gs,genes)
gene<-as.list(pathways)
sqy<-AddModuleScore(object = sqy,features = gene)
pathways_name<-"Autophagy"
colnames(sqy@meta.data)[21]<-pathways_name
rm(pathways,i,gene,pathways_name)

### 4 annotation
sqy@active.ident<-sqy$seurat_clusters
sqy<-RenameIdents(sqy,
                  "0"="Malignant",
                  "1"="Mono/Macro",
                  "2"="Malignant",
                  "3"="Malignant",
                  "4"="Malignant",
                  "5"="Malignant",
                  "6"="Mono/Macro",
                  "7"="Malignant",
                  "8"="Malignant",
                  "9"="Mono/Macro",
                  "10"="Malignant",
                  "11"="Mono/Macro",
                  "12"="Malignant",
                  "13"="Malignant",
                  "14"="Mono/Macro",
                  "15"="Mono/Macro",
                  "16"="Malignant",
                  "17"="Malignant",
                  "18"="Mono/Macro",
                  "19"="Oligodendrocyte",
                  "20"="Malignant",
                  "21"="Mono/Macro",
                  "22"="T/NK",
                  "23"="Malignant",
                  "24"="Malignant",
                  "25"="Malignant",
                  "26"="Malignant",
                  "27"="Mono/Macro",
                  "28"="Malignant",
                  "29"="Oligodendrocyte",
                  "30"="Malignant",
                  "31"="Malignant",
                  "32"="Mono/Macro"
)
sqy$cell_type<-sqy@active.ident
sqy@active.ident<-sqy$cell_type
pdf(file = "UMAP_cluster_celltype.pdf",width=5,height = 4)
UMAPPlot(object = sqy, label = TRUE, cols=cors)
dev.off()

Mono_Macro<-c("CD68","CD14")
Glioma<-c("SOX9","SOX2")
NK<-c("NKG7")
T<-c("CD3D","CD3E")
oligodendrocyte<-c("MOG","MAG","MBP")
DotPlot(sqy,features = c(Glioma,Mono_Macro,oligodendrocyte,T,NK),group.by = "cell_type")+
  RotatedAxis()

#### 5 CANX sub-group
table(sqy$cell_type)
malignant<-sqy[,sqy$cell_type%in%"Malignant"]
malignant@active.ident<-malignant$seurat_clusters
pdf(file = "UMAP_cluster_malignant.pdf",width=5,height = 4)
UMAPPlot(object = malignant, label = TRUE, cols=cors)
dev.off()
DotPlot(malignant,features = "CANX")
VlnPlot(malignant,features = "CANX",cols = cors,sort = TRUE)+
  geom_boxplot()+
  stat_compare_means(label.x = 3)+
  NoLegend()
FeaturePlot(malignant,features = "CANX",label = TRUE)
malignant<-RenameIdents(malignant,
                        "23"="High_CANX",
                        "8"="High_CANX",
                        "5"="High_CANX",
                        "0"="High_CANX",
                        "20"="High_CANX",
                        "13"="High_CANX",
                        "10"="High_CANX",
                        "16"="High_CANX",
                        "24"="Low_CANX",
                        "30"="Low_CANX",
                        "12"="Low_CANX",
                        "28"="Low_CANX",
                        "26"="Low_CANX",
                        "31"="Low_CANX",
                        "4"="Low_CANX",
                        "17"="Low_CANX",
                        "2"="Low_CANX",
                        "7"="Low_CANX",
                        "3"="Low_CANX",
                        "25"="Low_CANX"
)
malignant$CANX_type<-malignant@active.ident
table(malignant@active.ident)
malignant@active.ident<-malignant$CANX_type
pdf(file = "UMAP_cluster_malignant_canx_exp.pdf",width=5,height = 4)
UMAPPlot(object = malignant, label = TRUE, cols=c("#FF7A2D","#4262FF"))
dev.off()



##### 7 Enrichment analysis ######

library("clusterProfiler")
library("org.Hs.eg.db")


setwd("")             
rt=read.table("id.txt",sep="\t",header=T,check.names=F)         
rt=rt[is.na(rt[,"entrezID"])==F,]                       
gene=rt$entrezID

#GO富集分析
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)     



##### 8 Calculation of TMZ IC50 #####

library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="symbol.txt"  
setwd("")   

rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = data,
              batchCorrect = 'eb',      #"eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')

### Author : Qingyuan Sun, Shandong University. 
### Email : sunqingyuansdu@163.com.
