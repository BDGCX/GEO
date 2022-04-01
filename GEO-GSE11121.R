rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
# 注意查看下载文件的大小，检查数据 
f='GSE11121_eSet.Rdata'
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE11121', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE11121_eSet.Rdata')  ## 载入数据
class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])
gset


# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
# GPL6244
dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
dat=log2(dat+1)
boxplot(dat,las=2)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
## 挑选一些感兴趣的临床表型。
phe=pd[,c(43:46,48)]
library(stringr)
head(phe)
group_list=phe$`grade:ch1`
table(group_list)

dat[1:4,1:4] 

# GPL96
# 如果这个 GPL96 芯片平台，没有对应的数据R包，就使用下面的代码
if(F){
  library(GEOquery)
  #Download GPL file, put it in the current directory, and load it:
  gpl <- getGEO('GPL96', destdir=".")
  colnames(Table(gpl))  
  head(Table(gpl)[,c(1,15)]) ## you need to check this , which column do you need
  probe2gene=Table(gpl)[,c(1,15)]
  head(probe2gene)
  library(stringr)  
  save(probe2gene,file='probe2gene.Rdata')
}
# 
# load(file='probe2gene.Rdata')
# ids=probe2gene 
# 因为有R包，所以就直接载入R包即可。
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL) #toTable这个函数：通过看hgu133plus2.db这个包的说明书知道提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
head(ids) #head为查看前六行

head(ids)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat),]

dat[1:4,1:4]   
dat=dat[ids$probe_id,] 

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息


save(dat,phe,group_list,file = 'step1-output.Rdata')

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
table(group_list)
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
dat=as.data.frame(dat)#将matrix转换为data.frame
dat=cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
# The variable group_list (index = 54676) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')
# 因为我们的 group_list 分组，指的是 grade，从PCA图可以看到，PC2就是跟grade相关的那些基因。
# 这里大家可以尝试通过WGCNA把基因分成不同模块，然后看具体哪一个模块跟这个grade特别相关



rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'step1-output.Rdata') #此步为一个小插曲，即计算一下从第一行开是计算每一行的sd值，知道最后一行所需要的时间
dat[1:4,1:4] 

cg=names(tail(sort(apply(dat,1,sd)),1000))#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
n=t(scale(t(dat[cg,]))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group_list=group_list)
rownames(ac)=colnames(n) #把ac的行名给到n的列名，即对每一个探针标记上分组信息

pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')
# 前面的PCA图看的，只有PC2是跟我们的分组变量grade相关的
# 所以这个时候的top1000的sd基因，是不可能非常明显的区分不同grade的病人。


rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
# 每次都要检测数据
dat[1:4,1:4] 
group_list=paste0('grade',group_list)
table(group_list) #table函数，查看group_list中的分组个数
group_list=as.factor(group_list)
#通过为每个数据集绘制箱形图，比较数据集中的数据分布
boxplot(dat[1,]~group_list) #按照group_list分组画箱线图

bp=function(g){         #定义一个函数g，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat[1,]) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])
dim(dat)

library(limma)
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
## 上面是limma包用法的一种方式 
options(digits = 4) #设置全局的数字有效位数为4
#topTable(fit,coef=2,adjust='BH') 
topTable(fit,coef=2,adjust='BH') 
topTable(fit,coef=3,adjust='BH') 
## 但是上面的用法做不到随心所欲的指定任意两组进行比较

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=dat
rownames(design)=colnames(exprSet)
design
colnames(design)
contrast.matrix<-makeContrasts(contrasts=c("grade3-grade1",'grade2-grade1'),
                               levels = design)
contrast.matrix ##这个矩阵声明，我们要把 Tumor 组跟 Normal 进行差异分析比较

##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
##这一步很重要，大家可以自行看看效果

fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
nrDEG1 = topTable(fit2, coef="grade3-grade1", n=Inf)
nrDEG2 = topTable(fit2, coef='grade2-grade1', n=Inf)

# 后面的代码，我并没有调试，因为不属于本项目讲解内容
# 我们关注的应该是生存分析

# 整个第4步骤代码，我并没有调试，因为不属于本项目讲解内容
# 请大家不要运行。

rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'deg.Rdata')
head(deg)
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=1.5
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
deg$symbol=rownames(deg)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)


## KEGG pathway analysis
### 做KEGG数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
if(T){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  dotplot(kk.up );ggsave('kk.up.dotplot.png')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  dotplot(kk.down );ggsave('kk.down.dotplot.png')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  dotplot(kk.diff );ggsave('kk.diff.dotplot.png')
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  source('functions.R')
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = 'kegg_up_down.png')
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')
  
  
}

### GO database analysis 
### 做GO数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
{
  
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  
  if(F){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file = 'go_enrich_results.Rdata')
    
  }
  
  
  load(file = 'go_enrich_results.Rdata')
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
}


##生存分析 

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
table(group_list)
# 每次都要检测数据
dat[1:4,1:4]
head(phe)
# R里面实现生存分析非常简单！

# 用my.surv <- surv(OS_MONTHS,OS_STATUS=='DECEASED')构建生存曲线。
# 用kmfit2 <- survfit(my.surv~TUMOR_STAGE_2009)来做某一个因子的KM生存曲线。
# 用 survdiff(my.surv~type, data=dat)来看看这个因子的不同水平是否有显著差异，其中默认用是的logrank test 方法。
# 用coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung) 来检测自己感兴趣的因子是否受其它因子(age,gender等等)的影响。


exprSet=dat
dim(exprSet)
colnames(phe)
colnames(phe)=c('event','grade','node','size','time')
phe = as.data.frame(apply(phe,2,as.numeric))
boxplot(phe$size)
phe$size=ifelse(phe$size>median(phe$size),'big','small')

library(survival)
library(survminer)
# 利用ggsurvplot快速绘制漂亮的生存曲线图
sfit <- survfit(Surv(time, event)~size, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
## 多个 ggsurvplots作图生存曲线代码合并代码公布
sfit1=survfit(Surv(time, event)~size, data=phe)
sfit2=survfit(Surv(time, event)~grade, data=phe)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1,pval =TRUE, data = phe, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2,pval =TRUE, data = phe, risk.table = TRUE)
# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
# 可以看到grade跟生存显著相关，而size跟病人生存的关系并不显著。

## 挑选感兴趣的基因做差异分析
phe$CBX4=ifelse(exprSet['CBX4',]>median(exprSet['CBX4',]),'high','low')
table(phe$CBX4)
ggsurvplot(survfit(Surv(time, event)~CBX4, data=phe), conf.int=F, pval=TRUE)
phe$CBX6=ifelse(exprSet['CBX6',]>median(exprSet['CBX6',]),'high','low')
table(phe$CBX6)
ggsurvplot(survfit(Surv(time, event)~CBX6, data=phe), conf.int=F, pval=TRUE)
## 批量生存分析 使用  logrank test 方法
mySurv=with(phe,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  #gene=exprSet[1,]
  phe$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p=sort(log_rank_p)
head(log_rank_p)
boxplot(log_rank_p) 
phe$H6PD=ifelse(exprSet['H6PD',]>median(exprSet['H6PD',]),'high','low')
table(phe$H6PD)
ggsurvplot(survfit(Surv(time, event)~H6PD, data=phe), conf.int=F, pval=TRUE)

# 前面如果我们使用了WGCNA找到了跟grade相关的基因模块，然后在这里就可以跟生存分析的显著性基因做交集
# 这样就可以得到既能跟grade相关，又有临床预后意义的基因啦。

## 批量生存分析 使用 coxph 回归方法
colnames(phe)
mySurv=with(phe,Surv(time, event))
cox_results <-apply(exprSet , 1 , function(gene){
  group=ifelse(gene>median(gene),'high','low')
  survival_dat <- data.frame(group=group,grade=phe$grade,size=phe$size,stringsAsFactors = F)
  m=coxph(mySurv ~ grade + size + group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results=t(cox_results)
table(cox_results[,4]<0.05)
cox_results[cox_results[,4]<0.05,]

length(setdiff(rownames(cox_results[cox_results[,4]<0.05,]),
               names(log_rank_p[log_rank_p<0.05])
))
length(setdiff( names(log_rank_p[log_rank_p<0.05]),
                rownames(cox_results[cox_results[,4]<0.05,])
))
length(unique( names(log_rank_p[log_rank_p<0.05]),
               rownames(cox_results[cox_results[,4]<0.05,])
))
save(log_rank_p,cox_results,exprSet,phe,file = 'surviva.Rdata')
