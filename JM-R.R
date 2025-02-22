##作业1 ensembl_id转换
a=read.table('e1.txt')
View(a)
library(org.Hs.eg.db)
?org.Hs.eg.db
ls("package:org.Hs.eg.db")
g2s=toTable(org.Hs.egSYMBOL)
View(g2s)
head(a)
g2e=toTable(org.Hs.egENSEMBL)
View(g2e)
View(g2s)
library(stringr)
x=a$V1[1]
x
strsplit(x,'[.]')
strsplit(as.character(x),'[.]')
strsplit(as.character(x),'[.]')[[1]]
strsplit(as.character(x),'[.]')[[1]][1]
lapply(a$V1,function(x){
  strsplit(as.character(x),'[.]')[[1]][1]
})
unlist(lapply(a$V1,function(x){
  strsplit(as.character(x),'[.]')[[1]][1]
})
)
head(g2e)
a$ensembl_id=unlist(lapply(a$V1,function(x){
  strsplit(as.character(x),'[.]')[[1]][1]
})
)
View(a)
tmp=merge(a,g2e,by='ensembl_id')
View(tmp)
g2s=toTable(org.Hs.egSYMBOL);head(g2s)
tmp=merge(tmp,g2s,by='gene_id')
View(tmp)

##作业2 probe_id转换，首先找到对应平台GPL所需的R包并安装加载http://www.bio-info-trainee.com/1399.html
a=read.table('e2.txt')
a=read.table('e2.txt')
View(a)
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
rm(list = ls())
options(stringsAsFactors = F)
a=read.table('e2.txt')
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
View(ids)
View(a)
View(a)
View(ids)
head(ids)
tmp=merge(ids,a,by='probe_id')
colnames(a)='probe_id'
View(a)
tmp=merge(ids,a,by='probe_id')
View(tmp)
match(a$probe_id,ids$probe_id)
ids[match(a$probe_id,ids$probe_id),]
head(ids)
tmp1=merge(ids,a,by='probe_id')
tmp2=ids[match(a$probe_id,ids$probe_id),]
tmp1==tmp2

##作业3 找R包内置数据集TP53表达量并绘制boxplot
rm(list = ls())
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(CLL))
data(sCLLex) #内置数据集使用data()
sCLLex
exprSet=exprs(sCLLex) #使用exprs()获得表达矩阵
library(hgu95av2.db)
View(exprSet)
ids=toTable(hgu95av2SYMBOL)
head(ids)
View(ids)
exprSet['1939_at',]
pd=pData(sCLLex) #分组信息位于pData()的disease中
View(pd)
boxplot(exprSet['1939_at',] ~ pd$Disease)
boxplot(exprSet['1974_s_at',] ~ pd$Disease)
boxplot(exprSet['31618_at',] ~ pd$Disease)
boxplot(exprSet['1939_at',] ~ pd$Disease)

##作业4 BRCA1基因在TCGA数据库的乳腺癌数据集
#需要使用cBioportal网站下载文件a
rm(list = ls())
options(stringsAsFactors = F)
a=read.table('e4-plot.txt')
a=read.table('e4-plot.txt',sep = '\t')
a=read.table('e4-plot.txt',sep = '\t',fill = T)
View(a)
View(a)
View(a)
a=read.table('e4-plot.txt',sep = '\t',fill = T,header = T)
View(a)
View(a)
colnames(a)=c('id','subtype','expression','mut')
dat=a
library(ggstatsplot)
ggbetweenstats(data =dat, x = subtype,  y = expression)
ggsave('plot-again-BRCA1-TCGA-BRCA-cbioportal.png')
library(ggplot2)
ggsave('plot-again-BRCA1-TCGA-BRCA-cbioportal.png')

##作业5  TP53基因在TCGA数据库的乳腺癌数据集的表达量分组看其是否影响生存
#需要使用Oncolnc网站下载文件a
rm(list = ls())
options(stringsAsFactors = F)
a=read.table('BRCA_7157_50_50.csv',sep = ',',fill = T,header = T)
View(a)
dat=a
library(ggplot2)  #生存分析
library(survival)
library(survminer)
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=dat)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsave('survival_TP53_in_BRCA_TCGA.png') #生存分析
#subtype分组信息位于作业4的文件中
b=read.table('e4-plot.txt',sep = '\t',fill = T,header = T)
View(b)
View(a)
#使a和b列名一致
colnames(b)=c('Patient','subtype','expression','mut')
head(b)
substring(b$Patient,1,12)
b$Patient=substring(b$Patient,1,12)
tmp=merge(a,b,by='Patient')
View(tmp)
#分亚型做生存分析
table(tmp$subtype)
dat=tmp[tmp$subtype=='BRCA_Basal',]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=dat)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsave('survival_TP53_in_BRCA_Basal_TCGA.png')
table(tmp$subtype)
dat=tmp[tmp$subtype=='BRCA_Her2',]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=dat)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsave('survival_TP53_in_BRCA_Her2_TCGA.png')
table(tmp$subtype)
s_b <- function(dat){
  library(ggplot2)
  library(survival)
  library(survminer)
  table(dat$Status)
  dat$Status=ifelse(dat$Status=='Dead',1,0)
  sfit <- survfit(Surv(Days, Status)~Group, data=dat)
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  # ggsave('survival_TP53_in_BRCA_Her2_TCGA.png')
}
dat=tmp[tmp$subtype=='BRCA_Basal',]
s_b(dat)
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=dat)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
dat=tmp[tmp$subtype=='BRCA_Basal',]
s_b(dat)
tmp1=tmp[tmp$subtype=='BRCA_Basal',]
s_b(tmp1)
rm(dat)
tmp1=tmp[tmp$subtype=='BRCA_Basal',]
s_b(tmp1)
s_b <- function(dat){
  library(ggplot2)
  library(survival)
  library(survminer)
  table(dat$Status)
  dat$Status=ifelse(dat$Status=='Dead',1,0)
  sfit <- survfit(Surv(Days, Status)~Group, data=dat)
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
}
tmp1=tmp[tmp$subtype=='BRCA_Basal',]
s_b(tmp1)
s_b <- function(x){
  library(ggplot2)
  library(survival)
  library(survminer)
  table(x$Status)
  x$Status=ifelse(x$Status=='Dead',1,0)
  sfit <- survfit(Surv(Days, Status)~Group, data=x)
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
}
tmp1=tmp[tmp$subtype=='BRCA_Basal',]
s_b(tmp1)
s_b <- function(x){
  library(ggplot2)
  library(survival)
  library(survminer)
  #table(x$Status)
  x$Status=ifelse(x$Status=='Dead',1,0)
  sfit <- survfit(Surv(Days, Status)~Group, data=x)
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
}
s_b(tmp1)
table(tmp$subtype)
lapply(tmp$subtype,function(type){
  x=tmp[tmp$subtype==type,]
  library(ggplot2)
  library(survival)
  library(survminer)
  #table(x$Status)
  x$Status=ifelse(x$Status=='Dead',1,0)
  sfit <- survfit(Surv(Days, Status)~Group, data=x)
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
})
table(tmp$subtype)
type='BRCA_LumB'
x=tmp[tmp$subtype==type,]
library(ggplot2)
library(survival)
library(survminer)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
table(tmp$subtype)
type='BRCA_Normal'
x=tmp[tmp$subtype==type,]
library(ggplot2)
library(survival)
library(survminer)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
table(tmp$subtype)
type='BRCA_Basal'
x=tmp[tmp$subtype==type,]
library(ggplot2)
library(survival)
library(survminer)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
table(tmp$subtype)
type='BRCA_Her2'
x=tmp[tmp$subtype==type,]
library(ggplot2)
library(survival)
library(survminer)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
table(tmp$subtype)
type='BRCA_LumA'
x=tmp[tmp$subtype==type,]
library(ggplot2)
library(survival)
library(survminer)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)

##作业6  下载GSE17215的表达矩阵并提取以下基因绘制热图
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 包的下载代码，注意查看下载文件的大小，检查数据
f='GSE17215_eSet.Rdata'
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE17215', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE17215_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
View(dat)
#基因探针id转换
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
dat=dat[ids$probe_id,]
dat[1:4,1:4]
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]
dim(dat)
View(dat)

ng='ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T'
strsplit(ng,' ') #将ng按照空格分割
strsplit(ng,' ')[[1]]
ng=strsplit(ng,' ')[[1]]
dat[ng,]#会报错，因为有些是很多年前做的，有些基因检测不到
ng %in%  rownames(dat) #判断基因是否在dat中
table(ng %in%  rownames(dat))
ng=ng[ng %in%  rownames(dat)]
dat[ng,]
dat=dat[ng,]
#画热图，值比较大，先log2(一下)
dat=log2(dat)
pheatmap::pheatmap(dat)
pheatmap::pheatmap(dat,scale = 'row')

##作业7  下载GSE24673的表达矩阵计算样本的相关性、绘制热图并标记样本分组信息 
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据
f='GSE24673_eSet.Rdata'
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE24673', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE24673_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
View(dat)
pd=pData(a)
View(pd)
group_list=c('rbc','rbc','rbc',
             'rbn','rbn','rbn',
             'rbc','rbc','rbc',
             'normal','normal')
dat[1:4,1:4]
M=cor(dat)
View(M)
pheatmap::pheatmap(M)
tmp=data.frame(g=group_list)
rownames(tmp)=colnames(M)
pheatmap::pheatmap(M,annotation_col = tmp)
View(pd)

##作业8  找到GPL6244 platform of Affymetrix Human Gene 1.0 ST Array
##对应的bioconductor包并安装(去生信菜鸟团查找对应R包下载即可)
options()$repos
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("hugene10sttranscriptcluster.db",ask = F,update = F)#注释包后边要加".db"

##作业9  下载数据集GSE42872的表达矩阵，并分别挑选出所有样本的（平均表达量/sd/mad/)
##最大的探针并找到它们对应的基因
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据
f='GSE42872_eSet.Rdata'
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE42872_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)

sort(apply(dat,1,mean),decreasing = T)[1]
boxplot(dat)
sort(apply(dat,1,sd),decreasing = T)[1]
sort(apply(dat,1,mad),decreasing = T)[1]
#作业8重复
options()$repos
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("hugene10sttranscriptcluster.db",ask = F,update = F)

##作业10  下载数据集GSE42872的表达矩阵，并根据分组使用limma包做差异分析，得到差异结果矩阵
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据
f='GSE42872_eSet.Rdata'
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE42872_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
# (平均表达量/sd/mad/)最大的探针
boxplot(dat)
View(pd)

x=pd$title[1]
x
strsplit(x,' ')[[1]][4]
unlist(lapply(pd$title,function(x){
  strsplit(x,' ')[[1]][4]
}))
group_list=unlist(lapply(pd$title,function(x){
  strsplit(x,' ')[[1]][4]
}))
exprSet=dat
exprSet[1:4,1:4]
exprSet=dat
exprSet[1:4,1:4]
# DEG by limma
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix<-makeContrasts("progres.-stable",levels = design)
contrast.matrix
##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
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
