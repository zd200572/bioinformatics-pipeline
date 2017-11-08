options(stringsAsFactors=F)
library(ggplot2)
library(grid)
working_dir <-"d:\picrust"
setwd(working_dir)

if(dir.exists("vioplot")){
}else{
	dir.create("vioplot")
	dir.create("vioplot/quantile")
	dir.create("vioplot/p_value")
}

###读取数据###
###输出格式的三行的格式
group_file <- read.csv("group.tab",sep="\t",header=F)
violin_data <- read.csv("Galaxy56-[Categorize_by_function_on_data_31].tab",sep="\t",header=T,check.names=F,skip=1)
rownames(violin_data) <- violin_data[,dim(violin_data)[2]]
violin_data <- violin_data[,-c(1,dim(violin_data)[2])]

violin_sum <- unlist(apply(violin_data,2,sum))
violin_per <- violin_data
for (i in 1:dim(violin_per)[2]){
	violin_per[,i] <- violin_data[,i]/violin_sum[i]
}

###得到大分类的名字
get_violin <- function(class,group_name,group_sample,class_data){
	class <- paste("^",class,sep="")
	group_sample_split <- unlist(strsplit(group_sample,","))
	group_class_data <- class_data[,match(group_sample_split,colnames(class_data))]
	group_class_data <- group_class_data[grep(class,rownames(group_class_data)),]

	return_data <- NULL
	for(i in 1:dim(group_class_data)[1]){
		each_data <- unlist(group_class_data[i,])
		each_data <- cbind(each_data,group_name,rownames(group_class_data)[i])
		return_data <- rbind(return_data,each_data)
	}
	colnames(return_data) <-c("per","group","class")
	return(return_data)
}


big_class_data_frame <- as.data.frame(strsplit(rownames(violin_data),";"))
big_class <- unique(unlist(big_class_data_frame[1,]))
pdf_out <- paste("vioplot/","kegg_L2_",big_class,".pdf",sep="")
pdf_out <- gsub(" ","_",pdf_out)
for (i in 1:length(big_class)){
	plot_data <- NULL
	quantile_matrix <- matrix(nrow=5,ncol=dim(group_file[1]))
	rownames(quantile_matrix) <- c("min","1/4","2/4","3/4","max")
	colnames(quantile_matrix) <- group_file[,1]

	for (m in 1:dim(group_file)[1]){
		each_data <- get_violin(big_class[i],group_file[m,1],group_file[m,2],violin_per)
		quantile_matrix[,m] <- quantile(as.numeric(each_data[,1]))
		plot_data <- rbind(plot_data,each_data)
		plot_data <- as.data.frame(plot_data)
		plot_data[,1] <- as.numeric(plot_data[,1])
		rownames(plot_data) <- 1:dim(plot_data)[1]
	}

	###计算相关性###
	max_result <-max(plot_data[,1])*0.8
	small_class <- unique(plot_data[,3])
	small_class_matrix <- cbind(small_class,"p_value","sig")
	colnames(small_class_matrix) <- c("type","p_value","sig")
	for (tt in 1:length(small_class)){
		small_data <- plot_data[grep(small_class[tt],plot_data[,3], fixed=TRUE),]
		test1 <- small_data[grep(group_file[1,1],small_data[,2]),1]
		test2 <- small_data[grep(group_file[2,1],small_data[,2]),1]
		test_result <- t.test(test1,test2)
		p_value <- round(test_result$p.value,6)
		if (p_value <= 0.01){
			sig = "**"
		}
		if (p_value <= 0.05 && p_value > 0.01){
			sig = "*"
		}
		if (p_value > 0.05){
			sig = "NA"
		}
		small_class_matrix[tt,2] <- paste("p_value: ",p_value,sep="")
		small_class_matrix[tt,3] <- paste("sig: ",sig,sep="")
	}

	pdf_width=15
	pdf_height=10

	if (length(small_class) > 10 ){
		pdf_width = 20
	}
	if (length(small_class) > 20 ){
		pdf_width = 25
	}

	y_title = "The Relative Abundance"
	figure_title = paste("KEGG_level2: ",big_class[i],sep="")
	quantile_out = paste("vioplot/quantile/",big_class[i],"_quantile.xls",sep="");
	p_value_out = paste("vioplot/p_value/",big_class[i],"_p_value.xls",sep="");
	#print(quantile_out)


	pdf(pdf_out[i],width=pdf_width,height=pdf_height)
	plot <- ggplot(plot_data, aes(x=class, y=per, fill=group)) + geom_violin(position=position_dodge(1))
	plot <- plot + labs(title=figure_title,y=y_title,x="") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust=1)) + theme(plot.margin = unit(c(5,5,5,5),"cm"))
	plot <- plot + geom_boxplot(width=0.05, outlier.shape = NA, position=position_dodge(1))
	for(nn in 1:dim(small_class_matrix)[1]){
		plot <- plot + annotate("text",label=small_class_matrix[nn,3],x=small_class_matrix[nn,1],y=max_result, size= 3)
	}
	plot(plot)
	dev.off()

	###输出
	write.table(quantile_matrix,quantile_out,sep="\t",col.names=NA)
	write.table(small_class_matrix,p_value_out,sep="\t",row.names=F)

}