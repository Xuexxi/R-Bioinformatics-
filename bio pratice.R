# 设置工作目录
setwd("C:/Users/lizhu/Desktop/Bioinformatics 1/xena") 
library('tidyverse')
# 读取文件到R中，以excel表格的形式,后面两个参数不用管，以后读取tsv文件可以就用这个格式
counts_1 = read.table(file='TCGA-LIHC.htseq_counts.tsv',sep = '\t',header=TRUE)

 # “[]”举例
#如何提取（提取前三行），使用[行，列]
#counts_3rows = counts_1[1:3,]
#如何提取前三列
#counts_3col = counts_1[,1:3]

# 实操
#将第一列的123456框..去掉，改为基因名框
rownames(counts_1) = counts_1[,1]
#将第一例基因名去掉
counts_1 = counts_1[,-1]

 # substr 举例
#substr("wanglihong",1,4) #提取第一位到第四位

 # table 举例
table(substr(colnames(counts_1),14,16))

# "c"一个集合
# in% 符号用于判断是否属于
counts_1 = counts_1[,substr(colnames(counts_1),14,16) %in% c("01A","11A")]

#保留基因名的前十五位数
rownames(counts_1) = substr(rownames(counts_1),1,15)

# ceiling() 带小数的数字进一位

counts <- ceiling(2^(counts_1)-1)

#数据框输出为文本格式
write.table(counts,"counts.txt",sep = "\t",row.names=T,col.names=NA,quote=F)

                  # 6.9 - 6.24 代码全部能正常运行，我真牛逼！ L3
#提取行
GeneInfo <-GeneInfo[(GeneInfo$genetype == "protein_coding"),]

#提取共有的行名
comgene <- intersect(rownames(counts),rownames(GeneInfo))
counts <- counts[comgene,]
GeneInfo <- GeneInfo[comgene,]

#验证两个数据框的行名是否完全一致                          
a<- rownames(counts)
b<- rownames(GeneInfo)
identical(a,b)

#Genefo数据框中的genename列添加到counts数据框中
counts $ Gene <- as.character(GeneInfo$genename)
#去除重复的基因名称行
counts <- counts [!duplicated (counts$Gene),]
# 将行名变为Gene symbol
rownames(counts) <- counts $ Gene

# ncol()
# nrow()

# 去除最后一列 
counts <- counts[,-ncol(counts)]
write.table(counts,file="LIHC_counts_mRNA_all.txt",sep = "\t",row.names=T,col.names=NA,quote=F)

# 提取列名的第14到第16位字符， 确定是否等于01A
tumor <- colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A <-counts[,tumor]
write.table(counts_01A,file = "LIHC_counts_mRNA_all_01A.txt",sep = "\t",row.names=T,col.names=NA,quote=F)

# 差异分析
library(tidyverse)
#安装BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
counts = counts[apply(counts,1,function(x) sum(x>1)>32),]
conditions = data.frame(sample = colnames(counts),
                        group = factor(ifelse(substr(colnames(counts),14,16)=='01A','T','N'),
                                       levels = c('N','T'))) %>% column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(countData = counts,colData = conditions,design = ~ group)

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "LIHC_DEG.rda")
res_deseq2 <- as.data.frame(res)%>% arrange(padj) %>% dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)

# 6.26 数据筛选 p值>0.05具有显著差异，
setwd("C:/Users/lizhu/Desktop/Bioinformatics 1/xena") 
library(tidyverse)
fpkm_01A <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm_11A <- read.table("LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
load("C:/Users/lizhu/Desktop/Bioinformatics 1/xena/LIHC_DEG.rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)

gene <- c("LIN28B","CTAG2","REG3A")
a <- fpkm_01A[gene,] #提取与gene变量对应的行数据
b <- fpkm_11A[gene,] #提取与gene变量对应的行数据
a <- t(a) # 行列位置转换
b <- t(b) # 行列位置转换
class(a)
a <- as.data.frame(a) # 转换为csv文件
b <- as.data.frame(b) # 转换为csv文件
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "01A.csv")
write.csv(b, file = "11A.csv")

# 加载ggplot2包
library(ggplot2)
library(tidyr)
# 创建箱型图

ggplot(a, aes(x = gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression Boxplot", x = "Gene", y = "Expression Level")

# 7.4-7.5 GEO 数据整理
setwd("C:/Users/lizhu/Desktop/Bioinformatics 1/GSE84402") 
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)

gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)

#提取子集 ## [[ ]] 用于提取列表中的单个元素，而不是返回一个子列表
gset[[1]]

#通过pData函数获取分组信息
pdata <- pData(gset[[1]]) # pData 是一个函数，用于提取表达式集对象 # gset[[1]] 提取列表 gset 中的第一个表达式集对象
table(pdata$source_name_ch1) # source_name_ch1 是 pdata 数据框中的一个列名，通常表示样本来源或类别 
# table 函数用于对 source_name_ch1 列中的不同值进行计数，生成一个频数表
# 这行代码将统计 source_name_ch1 列中每个不同值的出现次数。
library(stringr)

#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor","normal")
# str_detect()这是一个字符串检测函数，用于检查 pdata$source_name_ch1 列中的每个元素是否包含字符串 "hepatocellular carcinoma"。
# ifelse(条件, 值1, 值2):这是一个条件语句函数。如果条件为真，则返回值1；否则返回值2。
# 在这段代码中，如果 str_detect 的结果为 TRUE，则返回 "tumor"；如果为 FALSE，则返回 "normal"。

# 向量型转化为分子型
group_list = factor(group_list, levels = c("normal","tumor"))
# 设置分类标签的顺序，normal为第一个，tumor为第二个exp <- exprs(gset[[1]])
# 假设 group_list 中有无效的颜色名称
group_list <- c("red", "blue", "tumor")
valid_colors <- c("red", "blue", "green") 
group_list <- ifelse(group_list %in% colors(), group_list, valid_colors[1])
boxplot(exp, outline=FALSE, notch=TRUE, col=group_list, las=2)
dev.off()

# 数据校正
library(limma)
exp=normalizeBetweenArrays(exp) # 消除系统差异
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp) # 查找最大值和最小值
exp <- log2(exp+1) # 最大最小值差距太多，做log转换缩小差距，+1来消除数值为0的数据
range(exp)
dev.off()

#使用R包转换第一列id 
index = gset[[1]]@annotation
if(!require("hgu133plus2.db"))
  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)

# id转换
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id = rownames(exp))
exp <- exp %>% inner_join(ids,by= "probe_id")
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
exp <- exp[,-(29:30)]
write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 7.8 

