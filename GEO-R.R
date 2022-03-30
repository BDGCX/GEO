#GEO数据挖掘
#R包安装
rm(list = ls())
options(stringsAsFactors = F)
options()$repos
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 
BiocManager::install("GEOquery",ask = F,update = F)
#1.GEO数据集下载,读取表达矩阵
library(GEOquery)
#无需每次都下载，只需要load即可
f<- "GSE30529_eSet.Rdata"
if(!file.exists(f)){
  gset <- getGEO('GSE30529', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE30529_eSet.Rdata')  ## 载入数据
#需每次下载的方法
gset <- getGEO("GSE42872",destdir = ".",
               AnnotGPL=F, getGPL=F)#后边两个代表不下载注释和平台文件soft
gset


#读取表达矩阵
exprSet <- read.table("GSE30529_series_matrix.txt.gz",header = T,sep = "\t",
                quote = "",fill = T,comment.char = "!")

#2.常用另一种读取表达矩阵的方法
a <- gset[[1]]
exprSet <- exprs(a)
class(exprSet)
# 3.生信技能书下载GEO数据集表达矩阵的封装函数
downGSE <- function(studyID = "GSE1009", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  
  exprSet <- exprs(eSet[[1]])
  pdata <- pData(eSet[[1]])
  
  write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
  
}
downGSE("GSE42872")

##探针id转换
#首先找到对应平台所需的R包并安装加载http://www.bio-info-trainee.com/1399.html
options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror
BiocManager::install("hgu133a2.db",ask = F,update = F)
library(hgu133a2.db)
ids <- toTable(hgu133a2SYMBOL)
View(ids)
#可以查看一些信息，可省略
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))
#查看对应关系，确定是否有探针不对应基因名的情况
exprSet1 <- exprSet
##GEO的探针名带有引号时
for (i in 1:nrow(exprSet1) ){
      x=exprSet1[i,1]  # 赋值
      x=as.character(x) #化作字符串
      a=gsub('["]', '', x)  #去双引号
      exprSet1[i,1]=a  #给矩阵重新赋值
}
rownames(exprSet1) <- exprSet1[,1]
exprSet1 <- exprSet1[,-1]
table(rownames(exprSet1)%in%ids$probe_id)
#对没有对应基因名的探针进行过滤，如都对应则省略
dim(exprSet1)#过滤前探针数
exprSet1 <- exprSet1[rownames(exprSet1)%in%ids$probe_id,]
dim(exprSet1)#过滤后探针数
ids=ids[match(rownames(exprSet1),ids$probe_id),]
head(ids)
exprSet1[1:5,1:5]
#根据基因symbol对表达矩阵进行分类后取平均值最大的探针，过滤有多个探针的基因
tmp <- by(exprSet1,ids$symbol,
      function(x) rownames(x)[which.max(rowMeans(x))] )
probes <- as.character(tmp)
dim(exprSet1)
exprSet1 <- exprSet1[rownames(exprSet1) %in% probes ,]
dim(exprSet1)#又过滤了有多个探针的基因
colnames(ids) <- c("X.ID_REF.","symbol")
tmp1 <- merge(ids,exprSet1,by='X.ID_REF.')
tmp1 <- tmp1[,-1]
rownames(tmp1) <- tmp1[,1]
tmp1 <- tmp1[,-1]
##DEG 差异基因分析
rm(list = ls()) 
options()$repos
library(GEOquery)
gset <- getGEO("GSE42872",destdir = ".",
               AnnotGPL=F, getGPL=F)

#定义group_list
a <- gset[[1]]
pd <- pData(a)
library(stringr)
group_list=str_split(pd$title,',',simplify = T)[,1]#取哪一列取决于title
table(group_list)

##其它取法
x <- pd$title[1]
x
strsplit(x,' ')[[1]][4]
unlist(lapply(pd$title,function(x){
  strsplit(x,' ')[[1]][4]
}))
group_list=unlist(lapply(pd$title,function(x){
  strsplit(x,' ')[[1]][4]
}))


# DEG by limma(差异分析)
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(tmp1)
design #分组矩阵
#contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
#paste0(unique(group_list),collapse = "-") #查看比较次序是否正确，应对照组在后
contrast.matrix<-makeContrasts("D-C",levels = design)
contrast.matrix
##这个矩阵声明，我们要把Vemurafenib与Control组进行差异分析比较
##step1
fit <- lmFit(tmp1,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput <-  topTable(fit2, coef=1, n=Inf)
nrDEG <-  na.omit(tempOutput)
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
View(nrDEG)
deg <- nrDEG[nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC)>2,]
#包装函数step1,2,3
deg <- function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  View(nrDEG)
  return(nrDEG)
}
deg(exprSet,design,contrast.matrix)

## heatmap 
library(pheatmap)
choose_gene=head(rownames(nrDEG),25)
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)

##volcano plot
#plot(nrDEG$logFC,-log10nrDEG$padj) 最基础的火山图
DEG <- nrDEG
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change <- as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
head(DEG)
library(ggplot2)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) 
## corresponding to the levels(res$change)
print(g)

##结构注释(富集分析)-clusterProfiler
library(clusterProfiler)#此包更新比较频繁
#gene和geneList是两个数据集
# get the universal genes and sDEG 

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

### Using MSigDB gene set collections
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)[,1:6]
gseaplot(egmt2, geneSetID = "EXTRACELLULAR_REGION")

## GO analysis 

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)[,1:6]
ego2 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(ego2)[,1:6]

## KEGG pathway analysis
#基因id转换
gene <- head(rownames(nrDEG),1000)#取差异基因的前1000个为例
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,#基因id必须是ENS_id
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)[,1:6]

#GSE分析
data(geneList, package="DOSE")
boxplot(geneList)
boxplot(nrDEG$logFC)#从箱线图可看出nrDEG$logFC中的数值
#与geneList相同，但没有名字，geneList名字是ENTREZID
#以下为如何从nrDEG$logFC得到DOSE中的geneList
geneList <- nrDEG$logFC
names(geneList) <- rownames(nrDEG)
geneList <- sort(geneList,decreasing = T)

gene.df <- bitr(names(geneList), 
                fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
tmp <- data.frame(SYMBOL=names(geneList),
                  logfc=as.numeric(geneList))
tmp <- merge(tmp,gene.df,by="SYMBOL")
geneList <- tmp$logfc
names(geneList) <- gene.df$ENTREZID
geneList#与DOSE包中的geneList相同

#进行GSE分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)[,1:6]
gseaplot(kk2, geneSetID = "hsa04145")

### Reactome pathway analysis
library(ReactomePA)
pp <- enrichPathway(gene         = gene,
                    organism     = 'human',
                    pvalueCutoff = 0.05)
head(pp)[,1:6]
pp2 <- gsePathway(geneList     = geneList,
                  organism     = 'human',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
head(pp2)[,1:6]

### Disease analysis
library(DOSE)
dd <- enrichDO(gene         = gene)
head(dd)[,1:6]
dd2 <- gseDO(geneList     = geneList,
             nPerm        = 1000,
             minGSSize    = 120,
             pvalueCutoff = 0.05,
             verbose      = FALSE)
head(dd2)[,1:6]

b=matrix(as.numeric(a),nrow=nrow(a))#转化为数值矩阵，提前保存行列名
