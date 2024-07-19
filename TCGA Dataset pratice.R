
### TCGA data Pre-processing ###

# set working directory 
setwd("xena")
library(tidyverse)
counts1 = read.table(file = 'TCGA-LIHC.htseq_counts.tsv', sep = '\t', header = TRUE) 
rownames(counts1) <- counts1[,1] 
counts1 = counts1[,-1]
counts1 <- counts1[,substr(colnames(counts1),14,16)%in% c("01A","11A")]
table(substr(colnames(counts1),14,16))

# Retain the first 15 digits of the line name
rownames(counts1) <- substr(rownames(counts1),1,15)
counts <- ceiling(2^(counts1)-1)
# Output 
write.table(counts,"counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.csv(counts, file = "counts.csv")

### TCGA Differential analysis ###

setwd("xena")
counts <- read.table("counts.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] 

# Get the intersection of row names
comgene <- intersect(rownames(counts),rownames(Ginfo))
counts <- counts[comgene,]
class(counts)
class(comgene)
Ginfo <- Ginfo[comgene,]
a <- rownames(counts)
b <- rownames(Ginfo)
identical(a,b)

counts$Gene <- as.character(Ginfo$genename)   # Gene Symbol 

counts <- counts[!duplicated(counts$Gene),]   # Remove duplicates
rownames(counts) <- counts$Gene   # Change the row name to Gene Symbol
ncol(Ginfo)
nrow
counts <- counts[,-ncol(counts)]   # Remove the last column
write.table(counts, file = "LIHC_counts_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

tumor <- colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A <- counts[,tumor]
write.table(counts_01A, file = "LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

### TCGA Differential analysis ###
library(tidyverse)
#Download BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)

counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == "01A","T","N"),levels = c("N","T"))) %>% 
  column_to_rownames("sample")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "LIHC_DEG.rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)

####COX regression analysis####

setwd("cox")

install.packages("survival")
install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)

surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 

surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]

expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]

res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)

deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
surv.expr <- cbind(surv,deg_expr)

#Cox analysis 
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}


write.table(Coxoutput, file = "cox results.txt",sep = "\t",row.names = F,col.names = T,quote = F)

#select top gene 
pcutoff <- 0.001
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] # Remove genes whose p-value is less than the threshold
topgene <- topgene[1:10,]

#3. 绘制森林图
##3.1 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林图
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "12" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()


####ESTIMATE immune score#####
setwd("TCGA ESTIMATE")  

library(utils) 
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#计算免疫评分
filterCommonGenes(input.f = "LIHC_fpkm_mRNA_01A.txt",   
                  output.f = "LIHC_fpkm_mRNA_01A.gct",  
                  id = "GeneSymbol")   
estimateScore("LIHC_fpkm_mRNA_01A.gct",   
              "LIHC_fpkm_mRNA_01A_estimate_score.txt",   
              platform="affymetrix")  #Default Platform

#3. Output the score for each sample
result <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
write.table(result, file = "LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分



#### ROC ####

setwd("ROC")
library(tidyverse)
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 

surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]

write.table(surv, file = "survival.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]

##Extract the 10 genes from the last cox plot

gene <- c("SPP1","PAGE1","G6PD","MAGEA4",'CDCA8',
          'TRIM54','KIF2C','KIF20A','ANLN',"SLC7A11")
exp10 <- expr[gene,] %>% t() %>% as.data.frame()

#Integrating
exp_sur <- cbind(exp10,surv)
write.table(exp_sur, file = "exp_sur.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


install.packages("ROCR")
install.packages("rms")
library(ROCR)
library(rms)

##Build ROC prediction model
ROC1 <- prediction(exp_sur$SPP1,exp_sur$OS)
ROC2 <- performance(ROC1,"tpr","fpr")   
AUC <- performance(ROC1,"auc")

#Assign AUC based on the results
AUC<- 0.5604839 

# Plotting the ROC Curve
plot(ROC2,
     col="red",   
     xlab="False positive rate", ylab="True positive rate",  
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=3)  
dev.off()

####timeROC####
setwd("timeROC")

install.packages("timeROC")
install.packages("survival")
library(timeROC)
library(survival)
library(tidyverse)


exp_sur <- read.table("exp_sur.txt", header=T,sep="\t", check.names=F, row.names=1)
exp_sur$OS.time <- exp_sur$OS.time/365
exp_sur_01A <- exp_sur[substr(rownames(exp_sur),14,16) == "01A",]
write.table(exp_sur_01A, file = "exp_sur_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# Constructing ROC curve function
ROC3 <- timeROC(T=exp_sur_01A$OS.time,   # Outcome time
                delta=exp_sur_01A$OS,   # Outcome indicator
                marker=exp_sur_01A$SPP1,   #Predictor variable
                cause=1,   #Positive outcome indicator value
                weighting="marginal",   #Calculation method, default is marginal
                times=c(1, 3, 5),   #Time point, select 1-year, 3-year and 5-year survival rate
                iid=TRUE)
ROC3

# Plotting the ROC Curve
plot(ROC3,
     time=1, col="red")  
plot(ROC3,
     time=3, col="green", add=TRUE)   #add refers to whether to add to the previous picture
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)   

dev.off()



#### TCGA differential analysis heat map ####
setwd("xena")
library(tidyverse)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
library(pheatmap)
cg = rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[cg,]
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

#### TCGA differential analysis volcano map ####
setwd("xena")
library(tidyverse)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)


install.packages("ggpubr")
install.packages("ggthemes")
library(ggpubr)
library(ggthemes)

DEG$logP <- -log10(DEG$padj)
ggscatter(DEG,
          x = "log2FoldChange", y = "logP") +
  theme_base()


#Add gene up- and down-regulation information
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base()

#Add dividing line
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
dev.off()

#Add gene tag
DEG$Label = ""   #Add a new label column
DEG <- DEG[order(DEG$padj), ]   #Sort the p-values from small to large
DEG$Gene <- rownames(DEG)
#Select the 5 genes with the smallest fdr value among the highly expressed genes
up.genes <- head(DEG$Gene[which(DEG$change == "UP")], 5)
#Select the 5 genes with the smallest fdr value among the lowly expressed genes
down.genes <- head(DEG$Gene[which(DEG$change == "DOWN")], 5)
#Merge up.genes and down.genes and add them to Label
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
