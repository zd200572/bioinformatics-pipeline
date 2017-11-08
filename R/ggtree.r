# 运行环境R3.4.1，ggtree版本1.8.1
# 安装ggtree包，末安装者将FALSE改为TRUE即可
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggtree"))
}

# 设置工作目录：选择 Session - Set working directory - To source file location,
# 我们的脚本ggtree.r位于测试数据文件夹根本目录，执行上面的操作可调置工作目录为脚本所在文件夹
rm(list=ls()) 
# 设置工作文件夹进入result，我们使用的大部分文件在此目录
setwd("result")
# 加载ggtree包
library("ggtree")

# 读入分析相关文件
# 读取树文件
tree <- read.tree("rep_seqs_k5.tree")
# 读取树物种注释信息
tax <- read.table("rep_seqs_k5.tax", row.names=1)
# 物种注释等级标签，共七级，但细菌末分类物种太多，一般只能在门、纲、目水平比较确定
colnames(tax) = c("kingdom","phylum","class","order")

# 按门水平建树并上色
## 给每个OTU按门分类分组，此处可以更改为其它分类级别，如纲、目等，即phylum替换为order或class即可
groupInfo <- split(row.names(tax), tax$phylum) # OTU and phylum for group
## 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
# 画树，按组上色
ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(size=3)



# 画圈图
## tiplab2保证标签自动角度，默认无图例，要显示需要+theme
pdf(file="ggtree_circle_color.pdf", width=9, height=5)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+geom_tiplab2(size=3)+ theme(legend.position = "right")
dev.off()



# 树+丰度热图
# 思路：矩形树右端添加每个样品的表达丰度。
## 读取OTU表
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")
## 读取实验设计
design = read.table("design.txt", header=T, row.names= 1, sep="\t")
## 取实验设计和OTU表中的交集:样本可能由于实验或测序量不足而舍弃掉，每次分析都要筛选数据
idx=intersect(rownames(design),colnames(otu_table))
sub_design=design[idx,]
## 按实验设计的样品顺序重排列
otu_table=otu_table[,idx]
## 将OTU表count转换为百分比
norm = t(t(otu_table)/colSums(otu_table,na=T)) * 100 # normalization to total 100
## 筛选树中OTU对应的数据
tax_per = norm[rownames(tax),]

## 保存树图于变量，align调置树OTU文字对齐，linesize设置虚线精细
p = ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(size=3, align=TRUE, linesize=.5)
p
pdf(file="ggtree_heat_sample.pdf", width=9, height=5)
## 添加数字矩阵
## offset设置两者间距，用于解决图重叠问题；width设置热图相对树图的宽度，解决热图和树图大小关系；font.size设置热图文字大小，解决文字过大重叠；colnames_angle调整热图标签角度，解决文字重叠问题；hjust调整热图标签位置，解决文字与热图重叠问题。
gheatmap(p, tax_per, offset = .15, width=3, font.size=3, colnames_angle=-45, hjust=-.1)
dev.off()



# 树+ 组均值热图
## 有时样本过多也无法展示和阅读，需要求各组均值展示：需要将分组信息添加至样品相对丰度表，再分类汇总
## 提取实验设计中的分组信息
sampFile = as.data.frame(sub_design$genotype,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
## OTU表转置，让样品名为行
mat_t = t(tax_per)
## 合并分组信息至丰度矩阵，并去除样品名列
mat_t2 = merge(sampFile, mat_t, by="row.names")[,-1]
## 按组求均值
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
## 去除非数据列并转置
mat_mean_final = do.call(rbind, mat_mean)[-1,]
## 重命名列名为组名
colnames(mat_mean_final) = mat_mean$group

## 按组均值热图
pdf(file="ggtree_heat_group.pdf", width=7, height=5)
gheatmap(p, mat_mean_final, offset = .05, width=1, font.size=3, hjust=-.1)
dev.off()




# Reference 原始代码测试过程
# 运行环境R3.4.1，ggtree版本1.8.1
# Work well in R3.3.3
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggtree","colorspace"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result")
library("ggtree")
library("colorspace")

# 读取树物种注释
tree <- read.tree("rep_seqs_k5.tree")
tax <- read.table("rep_seqs_k5.tax", row.names=1)
# 物种注释等级标签，共七级，但末知物种太多，一般只能在门、纲、目水平比较确定
colnames(tax) = c("kingdom","phylum","class","order")

# 按门水平建树上色
## 给每个OTU按门分类分组，此处可以更改为其它分类级别，如纲、目等，即phylum替换为order或class即可
groupInfo <- split(row.names(tax), tax$phylum) # OTU and phylum for group
## 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
# 直接画树
#ggtree(tree,aes(color=tax$phylum)) # 要求labels必须与树的结果+边数一致
# 画树，按组上色，
ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(size=3)

# 怎么图例没有了，去掉scale_color_manual
# 画圈图，tiplab2保证标签自动角度，默认无图例，要显示需要+theme
pdf(file="ggtree_circle_color.pdf", width=9, height=5)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+geom_tiplab2(size=3)+ theme(legend.position = "right")
dev.off()
# 用如下代码可以配色与之前保持一致
#+  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$phylum))+1)))
# 旧版代码，lables与theme冲突 
#  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$phylum))+1)), breaks=1:length(unique(tax$phylum)), labels=levels(tax$phylum))
#+geom_tiplab2(size=3)+  theme(legend.position = "right")
# 使用如下命令或在Rstudio中直接保存图片为指定大小或格式
# pdf(file="ggtree_phylum.pdf", width=8, height=8)
# dev.off()



# 树结果丰度热图
# 思路：矩形树右端添加每个样品的表达丰度。
## 读取OTU表
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")
## 读取实验设计
design = read.table("design.txt", header=T, row.names= 1, sep="\t")
## 取实验设计和OTU表中的交集:样本可能由于实验或测序量不足而舍弃掉，每次分析都要筛选数据
idx=intersect(rownames(design),colnames(otu_table))
sub_design=design[idx,]
## 按实验设计的样品顺序重排列
otu_table=otu_table[,idx]
## 将OTU表count转换为百分比
norm = t(t(otu_table)/colSums(otu_table,na=T)) * 100 # normalization to total 100
## 筛选树中OTU对应的数据
tax_per = norm[rownames(tax),]

# p = ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(size=3)
# gheatmap(p, tax_per, offset = .2, width=3, font.size=3, colnames_angle=-45, hjust=-.1)

## 保存树图于变量，align调置树OTU文字对齐，linesize设置虚线精细
p = ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(size=3, align=TRUE, linesize=.5)
p
pdf(file="ggtree_heat_sample.pdf", width=9, height=5)
## 添加数字矩阵
## offset设置两者间距，用于解决图重叠问题；width设置热图相对树图的宽度，解决热图和树图大小关系；font.size设置热图文字大小，解决文字过大重叠；colnames_angle调整热图标签角度，解决文字重叠问题；hjust调整热图标签位置，解决文字与热图重叠问题。
gheatmap(p, tax_per, offset = .15, width=3, font.size=3, colnames_angle=-45, hjust=-.1)
dev.off()



# 树+ 柱状图
## 只有热图适合展示大量样本，显示样品间重复性和组间差异关系；其它柱、饼图只适合展示各组均值
## 求各种均值：需要将分组信息添加至样品相对丰度表，再分类汇总
## 提取实验设计中的分组信息
sampFile = as.data.frame(sub_design$genotype,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
## OTU表转置，让样品名为行
mat_t = t(tax_per)
## 合并分组信息至丰度矩阵，并去除样品名列
mat_t2 = merge(sampFile, mat_t, by="row.names")[,-1]
## 按组求均值
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
## 去除非数据列并转置
mat_mean_final = do.call(rbind, mat_mean)[-1,]
## 重命名列名为组名
colnames(mat_mean_final) = mat_mean$group

## 按组均值热图
pdf(file="ggtree_heat_group.pdf", width=7, height=5)
gheatmap(p, mat_mean_final, offset = .05, width=1, font.size=3, hjust=-.1)
dev.off()



# 树+ 堆叠柱状图
biocLite(c("ggstance"))
library("ggstance")

# ？OTUID与stackplot重叠，stackplot有边框；
library("reshape2")
data_all = as.data.frame(melt(mat_mean_final, id.vars=row.names))
colnames(data_all)=c("OTU","Groups","value") # 分组名不要用group，容易冲突报错
p2 <- facet_plot(p, panel = 'Stacked Barplot', data = data_all, 
                 geom = geom_barh, 
                 mapping = aes(x = value, fill = Groups), 
                 stat='identity' )
pdf(file="ggtree_stack_group.pdf", width=7, height=5)
p2
dev.off()



# 树+箱线图
# ？一直画出图不规则，或无法找到分组OTUs，只按样品分组成功但不对应
data_all2 = as.data.frame(melt(norm, id.vars=row.names))
colnames(data_all2)=c("otus","Samples","value") # 分组名不要用group，容易冲突报错
p3 <- facet_plot(p2, panel="Boxplot", data=data_all2, geom_boxploth, 
                 mapping = aes(x=value, color= Samples))
p3
