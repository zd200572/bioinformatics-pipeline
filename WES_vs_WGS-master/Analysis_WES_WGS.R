# Aziz Belkadi july 2014
# Includes 5 functions :
#   1 Make_all_variants_graphs: Needs ggplot2 and scale packages. Make statistics and draw graphs for variant coverage, genotype quality and Minor allele ratio
#   2 Make_graph_discordant : Needs reshape and ggplot2 packages. Make statistics and draw graphs for concordant and discordant variant between WES and WGS.
#   3 Merge_and_write_table : Needs reshape and ggplot2 packages. Make statistics on concordant and discordant variant genotypes
#   4 Make_graph_coverage : Needs reshape and ggplot2 packages. Make statitics and graphs on gene coverage. Input files need to be arranged as columns containing percetnage of base paires covered by more than 8X tab separated and one gene by raw.
#   5 Exons_statistics : Needs biomaRt package. takes two parameters : a string type that could be : {'protein_coding','lincRNA','miRNA','snoRNA'} and a boolean plus50 to take into acount or not the 50Bp flanking regions. Needs exome kit interval list in /resources/ repertory. The interval list format should be : chr:beg-end. Make statistics on exons and RNA included in the WES kit.
# WES vcf files should be in separated files in the WES/ repertory
# WGS vcf files should be in separated files in the WGS/ repertory
# Works on both UG and HC vcf files
# aziz.belkadi@inserm.fr



Make_all_variants_graphs<-function(){ # Analyse vcf files and make Depth, GQ and MRR graphes for all WES and WGS variants
	
	library(ggplot2)
	library(scales)



	WES_files <- Sys.glob("WES/*.vcf")
	WGS_files <- Sys.glob("WGS/*.vcf")

	if(length(WES_files) != length(WGS_files)){
	
		stop("Different number of WES and WGS files")
	
	} else {
	
	for (i in 1 : length(WES_files)){
		
		exo=read.table(WES_files[i], sep="\t")
		gen=read.table(WGS_files[i], sep="\t")
		
		exo = exo[,c(1,2,3,4,5,10)]
		gen = gen[,c(1,2,3,4,5,10)]
		
		colnames(exo)=c("chr","pos","rs","ref","alt","geno1")
		colnames(gen)=c("chr","pos","rs","ref","alt","geno2")
		
		v_exo = strsplit(as.character(unlist(exo$geno1)),":")
		v1_exo = suppressWarnings(do.call(rbind, v_exo))
		
		
		if(dim(v1_exo)[2]==4){
			exo = exo[(v1_exo[,1] != v1_exo[,2])&(v1_exo[,1] != v1_exo[,3])&(v1_exo[,1] != v1_exo[,4]) ,]
			v1_exo = v1_exo[(v1_exo[,1] != v1_exo[,2])&(v1_exo[,1] != v1_exo[,3])&(v1_exo[,1] != v1_exo[,4]),]
			p_exo = strsplit(as.character(unlist(v1_exo[,2])),",")
			p1_exo = suppressWarnings(do.call(rbind, p_exo))
			exo[,7:8]=v1_exo[,1:2]
			exo[,9] = as.numeric(p1_exo[,1]) + as.numeric(p1_exo[,2])
			exo[,10]=v1_exo[,3]
		}else{	
			exo = exo[(v1_exo[,1] != v1_exo[,2])&(v1_exo[,1] != v1_exo[,3])&(v1_exo[,1] != v1_exo[,4])& (v1_exo[,1] != v1_exo[,5]),]
			v1_exo = v1_exo[(v1_exo[,1] != v1_exo[,2])&(v1_exo[,1] != v1_exo[,3])&(v1_exo[,1] != v1_exo[,4])&(v1_exo[,1] != v1_exo[,5]),]
			exo[,7:10]=v1_exo[,1:4]
		}
		exo=exo[exo[,7]!="1/2" & exo[,7]!="2/1" & exo[,7]!="2/2",]
		
		
		
			
		
		vcf_name = do.call(rbind,strsplit(do.call(rbind, strsplit(WES_files[i],"/"))[2],"[.]"))[1]
		exo[,11] = vcf_name
		
		
		ex = exo[exo[,7]=="0/1",]

		
		v_exo=strsplit(as.character(unlist(ex[,8])),",")
		v1_exo = do.call(rbind, v_exo)
		ex[,12:13]=as.numeric(v1_exo)

		ex[,14]=ifelse( ex[,12] < ex[,13], ex[,12] / (ex[,12]+ex[,13]) , ex[,13]/(ex[,12]+ex[,13]))
		ex = ex[!is.na(ex[,14]),]
		
		
		
		v_gen = strsplit(as.character(unlist(gen$geno2)),":")
		v1_gen = suppressWarnings(do.call(rbind, v_gen))
		
		
		if(dim(v1_gen)[2]==4){
			gen = gen[(v1_gen[,1] != v1_gen[,2])&(v1_gen[,1] != v1_gen[,3])&(v1_gen[,1] != v1_gen[,4]),]
			v1_gen = v1_gen[(v1_gen[,1] != v1_gen[,2])&(v1_gen[,1] != v1_gen[,3])&(v1_gen[,1] != v1_gen[,4]),]
			p_gen = strsplit(as.character(unlist(v1_gen[,2])),",")
			p1_gen = suppressWarnings(do.call(rbind, p_gen))
			gen[,7:8]=v1_gen[,1:2]
			gen[,9] = as.numeric(p1_gen[,1]) + as.numeric(p1_gen[,2])
			gen[,10]=v1_gen[,3]
		}else{	
			gen = gen[(v1_gen[,1] != v1_gen[,2])&(v1_gen[,1] != v1_gen[,3])&(v1_gen[,1] != v1_gen[,4])&(v1_gen[,1] != v1_gen[,5]),]
			v1_gen = v1_gen[(v1_gen[,1] != v1_gen[,2])&(v1_gen[,1] != v1_gen[,3])&(v1_gen[,1] != v1_gen[,4])&(v1_gen[,1] != v1_gen[,5]),]
			gen[,7:10]=v1_gen[,1:4]
		}
		gen=gen[gen[,7]!="1/2" & gen[,7]!="2/1" & gen[,7]!="2/2",]
		
		
		
		
		vcf_name = do.call(rbind,strsplit(do.call(rbind, strsplit(WGS_files[i],"/"))[2],"[.]"))[1]
		gen[,11] = vcf_name
		
		ge = gen[gen[,7]=="0/1",]

		v_gen=strsplit(as.character(unlist(ge[,8])),",")
		v1_gen = do.call(rbind, v_gen)
		ge[,12:13]=as.numeric(v1_gen)
		
		ge[,14]=ifelse(ge[,12] < ge[,13], ge[,12] / (ge[,12]+ge[,13]) , ge[,13]/(ge[,12]+ge[,13]))
		ge = ge[!is.na(ge[,14]),]

		
		
		
		if(i==1){
			exo1=exo
			gen1=gen
			ex1<-ex
			ge1<-ge
		}else{
			exo1=rbind(exo1,exo)
			gen1=rbind(gen1,gen)
			ex1<-rbind(ex1,ex)
			ge1<-rbind(ge1,ge)
		}

		
		if(i==length(WES_files)){
		
			colnames(exo1)[9]="Depth"
			colnames(exo1)[10]="GQ"
			colnames(exo1)[11]="Ind"
			
			colnames(ex1)[11]="Ind"
			colnames(ex1)[14]="MRR"
			
			
			colnames(gen1)[9]="Depth"
			colnames(gen1)[10]="GQ"
			colnames(gen1)[11]="Ind"
			
			colnames(ge1)[11]="Ind"
			colnames(ge1)[14]="MRR"
			
			
			exo1[,12] = "WES"
			gen1[,12] = "WGS"
			
			colnames(exo1)[6] = "geno"
			colnames(gen1)[6] = "geno"
			
			bot = rbind(exo1,gen1)
			colnames(bot)[12]="Method"
			
			
			ex1[,15] = "WES"
			ge1[,15] = "WGS"
			
			colnames(ex1)[6] = "geno"
			colnames(ge1)[6] = "geno"
			
			both = rbind(ex1,ge1)
			colnames(both)[15]="Method"
			
			
			
			
			
			g=ggplot(exo1[as.numeric(exo1$Depth)<=100,], aes(as.numeric(Depth), colour=Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m1 = max(ymax[[1]]$y)
			g=ggplot(gen1[as.numeric(gen1$Depth)<=100,], aes(as.numeric(Depth), colour=Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m2 = max(ymax[[1]]$y)
			g=ggplot(bot[as.numeric(bot$Depth)<=100,], aes(as.numeric(Depth), group=Method, colour=Method)) + stat_density(aes(y=..count..), position="identity",geom="line") 
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m3 = max(ymax[[1]]$y)/6
			max_depth=max(c(m1,m2,m3))
			
			for(i in seq(4000,120000,4000)){
				if((max_depth %% i) == max_depth) {
					val = i
					break
				}
			}


			g=ggplot(exo1[as.numeric(exo1$Depth)<=100,], aes(as.numeric(Depth), colour=Ind)) + stat_density(aes(y=..count..), position="identity",geom="line")  + scale_x_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) + scale_y_continuous(limits=c(0,val), breaks=seq(0,val,val/5), labels=scientific_format()(seq(0,val,val/5))) + theme_bw()		
			suppressWarnings(ggsave("all_WES_Depth.pdf", g, useDingbats=FALSE))
			
			g=ggplot(gen1[as.numeric(gen1$Depth)<=100,], aes(as.numeric(Depth), colour=Ind)) + stat_density(aes(y=..count..), position="identity",geom="line")  + scale_x_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) + scale_y_continuous(limits=c(0,val), breaks=seq(0,val,val/5), labels=scientific_format()(seq(0,val,val/5))) + theme_bw()
			suppressWarnings(ggsave("all_WGS_Depth.pdf", g, useDingbats=FALSE))
						
			g=ggplot(bot[as.numeric(bot$Depth)<=100,], aes(as.numeric(Depth), group=Method, colour=Method)) + stat_density(aes(y=..count..), position="identity",geom="line")  + scale_x_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) + scale_y_continuous(limits=c(0,val*6), breaks=seq(0,val*6,val * 6 / 5), labels=scientific_format()(seq(0,val,val/5))) + theme_bw()
			suppressWarnings(ggsave("Merged_Depth.pdf", g, useDingbats=FALSE))
			
			
			
			
			
			g=ggplot(exo1, aes(as.numeric(GQ), colour=Ind)) +  geom_freqpoly(binwidth = diff(range(as.numeric(exo1$GQ)))/30)
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m1 = max(ymax[[1]]$y)
			g=ggplot(gen1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(gen1$GQ)))/30)
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m2 = max(ymax[[1]]$y)
			g=ggplot(bot, aes(as.numeric(GQ), group=Method, colour=Method)) + geom_freqpoly(size=1.5, binwidth = diff(range(as.numeric(bot$GQ)))/30) 
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m3 = max(ymax[[1]]$y)/6
			max_GQ=max(c(m1,m2,m3))
			
			
			for(i in c(1000,10000,100000,1000000)){
				if((max_GQ %% i) == max_GQ) {
					val = i
					break
				}
			}

			
			
			
			
			
			
			if(val == 1000){
			
			g=ggplot(exo1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(exo1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100)) + scale_y_log10(limits = c(10^1,10^3), breaks = c(10,100,1000), labels=c(scientific_format()(10),scientific_format()(100),scientific_format()(1000)))+ theme_bw()
			suppressWarnings(ggsave("all_WES_GQ.pdf", g, useDingbats=FALSE))
						
			g=ggplot(gen1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(gen1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^3), breaks = c(10,100,1000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000))) + theme_bw()
			suppressWarnings(ggsave("all_WGS_GQ.pdf", g, useDingbats=FALSE))
			
			g=ggplot(bot, aes(as.numeric(GQ), group=Method, colour=Method)) + geom_freqpoly(size=1.5, binwidth = diff(range(as.numeric(bot$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(60^1, 6 * 10^3), breaks = c(60,600,6000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000)))+ theme_bw()
			suppressWarnings(ggsave("Merged_GQ.pdf", g, useDingbats=FALSE)) } else if(val == 10000){
			
			g=ggplot(exo1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(exo1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^4), breaks = c(10,100,1000,10000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000)))+ theme_bw()
			suppressWarnings(ggsave("all_WES_GQ.pdf", g, useDingbats=FALSE))
						
			g=ggplot(gen1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(gen1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^4), breaks = c(10,100,1000,10000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000))) + theme_bw()
			suppressWarnings(ggsave("all_WGS_GQ.pdf", g, useDingbats=FALSE))
			
			g=ggplot(bot, aes(as.numeric(GQ), group=Method, colour=Method)) + geom_freqpoly(size=1.5, binwidth = diff(range(as.numeric(bot$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(60^1, 6 * 10^4), breaks = c(60,600,6000,60000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000)))+ theme_bw()
			suppressWarnings(ggsave("Merged_GQ.pdf", g, useDingbats=FALSE)) } else if(val == 100000){
			
			g=ggplot(exo1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(exo1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^5), breaks = c(10,100,1000,10000,100000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000),scientific_format()(100000)))+ theme_bw()
			suppressWarnings(ggsave("all_WES_GQ.pdf", g, useDingbats=FALSE))
						
			g=ggplot(gen1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(gen1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^5), breaks = c(10,100,1000,10000,100000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000),scientific_format()(100000))) + theme_bw()
			suppressWarnings(ggsave("all_WGS_GQ.pdf", g, useDingbats=FALSE))
			
			g=ggplot(bot, aes(as.numeric(GQ), group=Method, colour=Method)) + geom_freqpoly(size=1.5, binwidth = diff(range(as.numeric(bot$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(60^1, 6 * 10^5), breaks = c(60,600,6000,60000,600000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000),scientific_format()(100000)))+ theme_bw()
			suppressWarnings(ggsave("Merged_GQ.pdf", g, useDingbats=FALSE)) } else{
			
			g=ggplot(exo1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(exo1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^6), breaks = c(10,100,1000,10000,100000,1000000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000),scientific_format()(100000),scientific_format()(1000000)))+ theme_bw()
			suppressWarnings(ggsave("all_WES_GQ.pdf", g, useDingbats=FALSE))
						
			g=ggplot(gen1, aes(as.numeric(GQ), colour=Ind)) + geom_freqpoly(binwidth = diff(range(as.numeric(gen1$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100))+ scale_y_log10(limits = c(10^1,10^6), breaks = c(10,100,1000,10000,100000,1000000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000),scientific_format()(100000),scientific_format()(1000000))) + theme_bw()
			suppressWarnings(ggsave("all_WGS_GQ.pdf", g, useDingbats=FALSE))
			
			g=ggplot(bot, aes(as.numeric(GQ), group=Method, colour=Method)) + geom_freqpoly(size=1.5, binwidth = diff(range(as.numeric(bot$GQ)))/30) + scale_x_continuous(limits=c(0,105), breaks=c(0,20,40,60,80,100)) + scale_y_log10(limits = c(60^1, 6 * 10^6), breaks = c(60,600,6000,60000,600000,6000000), labels=c(scientific_format()(10), scientific_format()(100),scientific_format()(1000),scientific_format()(10000),scientific_format()(100000),scientific_format()(1000000)))+ theme_bw()
			suppressWarnings(ggsave("Merged_GQ.pdf", g, useDingbats=FALSE)) }
			
			
			
			
			
			
			g=ggplot(ex1, aes(as.numeric(MRR), group = Ind, color = Ind)) + stat_density(aes(y = ..count..),position="identity",geom="line") 
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m1 = max(ymax[[1]]$y)/10
			
			g=ggplot(ge1, aes(as.numeric(MRR), group = Ind, color = Ind)) + stat_density(aes(y = ..count..),position="identity",geom="line")
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m2 = max(ymax[[1]]$y)/10
			
			g=ggplot(both, aes(as.numeric(MRR), group = Method, color = Method)) + stat_density(aes(y = ..count..),size= 1.5 , position="identity",geom="line")  
			ggb<-ggplot_build(g)
			ymax<-(ggb$data)
			m3 = max(ymax[[1]]$y)/60
			
			max_MRR=max(c(m1,m2,m3))
			
			
			
			for(i in seq(5000,150000,5000)){
				if((max_MRR %% i) == max_MRR) {
					val = i
					break
				}
			}

			
			
			
			g=ggplot(ex1, aes(as.numeric(MRR), group = Ind, color = Ind)) + stat_density(aes(y = ..count..),position="identity",geom="line") + scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + scale_y_continuous(limits=c(0,val*10), breaks=seq(0,val*10,val*10/5), labels=scientific_format()(seq(0,val,val/5))) + theme_bw()
			suppressWarnings(ggsave("all_het_WES_MRR.pdf", g, useDingbats=FALSE))
			
			g=ggplot(ge1, aes(as.numeric(MRR), group = Ind, color = Ind)) + stat_density(aes(y = ..count..),position="identity",geom="line") + scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + scale_y_continuous(limits=c(0,val*10), breaks=seq(0,val*10,val*10/5), labels=scientific_format()(seq(0,val,val/5)))  + theme_bw()
			suppressWarnings(ggsave("all_het_WGS_MRR.pdf", g, useDingbats=FALSE))
			
			g=ggplot(both, aes(as.numeric(MRR), group = Method, color = Method)) + stat_density(aes(y = ..count..),size= 1.5 , position="identity",geom="line")  + scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + scale_y_continuous(limits=c(0,val*60), breaks=seq(0,val*60,val*60/5), labels=scientific_format()(seq(0,val,val/5)))  + theme_bw()
			suppressWarnings(ggsave("Merged_MRR.pdf", g, useDingbats=FALSE))
			
			}


		}


	}

}












Make_graph_discordant<-function(){ # Analyse discordant variants in WES and WGS and make Depth, GQ and MRR graphes for WES and WGS 

	library(reshape)
	library(ggplot2)
	

	WES_files <- Sys.glob("WES/*.vcf")
	WGS_files <- Sys.glob("WGS/*.vcf")
	

	if(length(WES_files) != length(WGS_files)){
	
		stop("Different number of WES and WGS files")
	
	} else {
	
		for (i in 1 : length(WES_files)){
			
			exo=read.table(WES_files[i], sep="\t")
			gen=read.table(WGS_files[i], sep="\t")
		
			exo = exo[,c(1,2,3,4,5,10)]
			gen = gen[,c(1,2,3,4,5,10)]
		
			colnames(exo)=c("chr","pos","rs","ref","alt","geno1")
			colnames(gen)=c("chr","pos","rs","ref","alt","geno2")
						
			
			
			vcf_name = do.call(rbind,strsplit(do.call(rbind, strsplit(WGS_files[i],"/"))[2],"[.]"))[1]
			exo[,7]=vcf_name
			gen[,7]=vcf_name
			
			colnames(exo)[7]="WES_Ind"
			colnames(gen)[7]="WGS_Ind"

			


			v=strsplit(as.character(unlist(exo$geno1)),":")
			v1 = suppressWarnings(do.call(rbind, v))
			
			if(dim(v1)[2]==4){
				exo = exo[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				p = strsplit(as.character(unlist(v1[,2])),",")
				p1 = suppressWarnings(do.call(rbind, p))
				exo[,8:9]=v1[,1:2]
				exo[,10] = as.numeric(p1[,1]) + as.numeric(p1[,2])
				exo[,11]=v1[,3]
			}else{	
				exo = exo[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				exo[,8:11]=v1[,1:4]
			}
			
			colnames(exo)[8]="WES_genotype"
			colnames(exo)[9]="WES_Allele"
			colnames(exo)[10]="WES_Depth"
			colnames(exo)[11]="WES_GQ"


			
			

			v=strsplit(as.character(unlist(gen$geno2)),":")
			v1 = suppressWarnings(do.call(rbind, v))
			
			if(dim(v1)[2]==4){
				gen = gen[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				p = strsplit(as.character(unlist(v1[,2])),",")
				p1 = suppressWarnings(do.call(rbind, p))
				gen[,8:9]=v1[,1:2]
				gen[,10] = as.numeric(p1[,1]) + as.numeric(p1[,2])
				gen[,11]=v1[,3]
			}else{	
				gen = gen[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				gen[,8:11]=v1[,1:4]
			}
			

			colnames(gen)[8]="WGS_genotype"
			colnames(gen)[9]="WGS_Allele"
			colnames(gen)[10]="WGS_Depth"
			colnames(gen)[11]="WGS_GQ"

			
			m=merge(gen, exo)

			if(i==1){
				m1=m
				}else{
				m1=rbind(m1,m)
			}
		

			if(i==length(WES_files)){
			
			
			

			
			
				
			
				mdata = m1[((m1$WES_genotype == "1/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "0/1" & m1$WGS_genotype == "1/1")) & as.numeric(m1$WES_Depth) < 100,]
				g = ggplot(mdata, aes(x=as.numeric(WES_Depth) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WES_Depth) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WES_WGS_discordant_depth.pdf", g, useDingbats=FALSE)
				
				mdata = m1[(m1$WES_genotype == "0/1" & m1$WGS_genotype == "1/1") & as.numeric(m1$WES_Depth) < 100,]

				v=strsplit(as.character(unlist(mdata$WES_Allele)),",")
				v1 = do.call(rbind, v)
				mdata[,17:18] = v1
				mdata[,17] = as.numeric(mdata[,17])
				mdata[,18] = as.numeric(mdata[,18])


				mdata[,19]=ifelse( mdata[,17] < mdata[,18], mdata[,17] / (mdata[,17]+ mdata[,18]) , mdata[,18]/(mdata[,17] + mdata[,18]))
				colnames(mdata)[19]="WES_MRR"
		
				g = ggplot(mdata, aes(x=as.numeric(WES_MRR) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WES_MRR) , colour=WES_Ind, group=WES_Ind)) + scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WES_WGS_discordant_MRR.pdf", g, useDingbats=FALSE)
				
				
				
				
				
				
				
				
				
				
				mdata = m1[((m1$WES_genotype == "0/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "1/1" & m1$WGS_genotype == "1/1")) & as.numeric(m1$WES_Depth) < 100,]
				g = ggplot(mdata, aes(x=as.numeric(WES_Depth) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WES_Depth) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WES_WGS_concordant_depth.pdf", g, useDingbats=FALSE)
				
				mdata = m1[(m1$WES_genotype == "0/1" & m1$WGS_genotype == "0/1") & as.numeric(m1$WES_Depth) < 100,]

				v=strsplit(as.character(unlist(mdata$WES_Allele)),",")
				v1 = do.call(rbind, v)
				mdata[,17:18] = v1
				mdata[,17] = as.numeric(mdata[,17])
				mdata[,18] = as.numeric(mdata[,18])

				mdata[,19]=ifelse( mdata[,17] < mdata[,18], mdata[,17] / (mdata[,17]+ mdata[,18]) , mdata[,18]/(mdata[,17] + mdata[,18]))
				colnames(mdata)[19]="WES_MRR"
		
				g = ggplot(mdata, aes(x=as.numeric(WES_MRR) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WES_MRR) , colour=WES_Ind, group=WES_Ind)) + scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WES_WGS_concordant_MRR.pdf", g, useDingbats=FALSE)
				
				
				
				
				
				
				
				
				
				
				
				mdata = m1[((m1$WES_genotype == "1/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "0/1" & m1$WGS_genotype == "1/1")) & as.numeric(m1$WGS_Depth) < 100,]
				g = ggplot(mdata, aes(x=as.numeric(WGS_Depth) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WGS_Depth) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WGS_WES_discordant_depth.pdf", g, useDingbats=FALSE)
				
				mdata = m1[(m1$WGS_genotype == "0/1" & m1$WES_genotype == "1/1") & as.numeric(m1$WES_Depth) < 100,]
				
				v=strsplit(as.character(unlist(mdata$WGS_Allele)),",")
				v1 = do.call(rbind, v)
				mdata[,17:18] = v1
				mdata[,17] = as.numeric(mdata[,17])
				mdata[,18] = as.numeric(mdata[,18])

				mdata[,19]=ifelse( mdata[,17] < mdata[,18], mdata[,17] / (mdata[,17]+ mdata[,18]) , mdata[,18]/(mdata[,17] + mdata[,18]))
				colnames(mdata)[19]="WGS_MRR"
		
				g = ggplot(mdata, aes(x=as.numeric(WGS_MRR) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WGS_MRR) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WGS_WES_discordant_MRR.pdf", g, useDingbats=FALSE)
				
				
				
				
				
				
				
				
				
				mdata = m1[((m1$WES_genotype == "0/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "1/1" & m1$WGS_genotype == "1/1")) & as.numeric(m1$WGS_Depth) < 100,]
				g = ggplot(mdata, aes(x=as.numeric(WGS_Depth) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WGS_Depth) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WGS_WES_concordant_depth.pdf", g, useDingbats=FALSE)
				
				
				mdata = m1[(m1$WGS_genotype == "0/1" & m1$WES_genotype == "0/1") & as.numeric(m1$WES_Depth) < 100,]

				v=strsplit(as.character(unlist(mdata$WGS_Allele)),",")
				v1 = do.call(rbind, v)
				mdata[,17:18] = v1
				mdata[,17] = as.numeric(mdata[,17])
				mdata[,18] = as.numeric(mdata[,18])

				mdata[,19]=ifelse( mdata[,17] < mdata[,18], mdata[,17] / (mdata[,17]+ mdata[,18]) , mdata[,18]/(mdata[,17] + mdata[,18]))
				colnames(mdata)[19]="WGS_MRR"
		
				g = ggplot(mdata, aes(x=as.numeric(WGS_MRR) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WGS_MRR) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") +  scale_x_continuous(limits=c(0,0.5),breaks=c(0,1/7,1/6,1/5,1/4,1/3,1/2), labels=c("0","1/7","1/6","1/5","1/4","1/3","1/2")) + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WGS_WES_concordant_MRR.pdf", g, useDingbats=FALSE)
				
				
				
				
				
				
				
				
				mdata = m1[((m1$WES_genotype == "1/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "0/1" & m1$WGS_genotype == "1/1")) ,]
				g = ggplot(mdata, aes(x=as.numeric(WES_GQ) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WES_GQ) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WES_WGS_discordant_GQ.pdf", g, useDingbats=FALSE)
				
				
				mdata = m1[((m1$WES_genotype == "0/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "1/1" & m1$WGS_genotype == "1/1")) ,]
				g = ggplot(mdata, aes(x=as.numeric(WES_GQ) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WES_GQ) , colour=WES_Ind, group=WES_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WES_WGS_concordant_GQ.pdf", g, useDingbats=FALSE)
				
				
				
				
				mdata = m1[((m1$WES_genotype == "1/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "0/1" & m1$WGS_genotype == "1/1")) ,]
				g = ggplot(mdata, aes(x=as.numeric(WGS_GQ) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WGS_GQ) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WGS_WES_discordant_GQ.pdf", g, useDingbats=FALSE)
				
				
				mdata = m1[((m1$WES_genotype == "0/1" & m1$WGS_genotype == "0/1") | (m1$WES_genotype == "1/1" & m1$WGS_genotype == "1/1")) ,]
				g = ggplot(mdata, aes(x=as.numeric(WGS_GQ) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") 
				ggb<-ggplot_build(g)
				ymax<-(ggb$data)
				l = max(ymax[[1]]$y)
				g = ggplot(mdata, aes(x=as.numeric(WGS_GQ) , colour=WGS_Ind, group=WGS_Ind)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_y_continuous(limits=c(0,l), breaks=seq(0,l,l/5), labels=seq(0,1,0.2)) + theme_bw()
				ggsave("WGS_WES_concordant_GQ.pdf", g, useDingbats=FALSE)

				
		
			}
		
		}
	}
	
}












Merge_and_write_table<-function(){ # Make statistics on WES and WGS variants

	library(reshape)
	library(ggplot2)
	

	WES_files <- Sys.glob("WES/*.vcf")
	WGS_files <- Sys.glob("WGS/*.vcf")
	
	all_var_table = matrix(0, ncol=9, nrow=length(WES_files))
	merge_table = matrix(0, ncol=6, nrow=length(WGS_files))

	colnames(all_var_table)=c("ind","WES_all_var","WES_complex","WES_homo","WES_hetero","WGS_all_var","WGS_complex","WGS_homo","WGS_hetero")	
	colnames(merge_table)=c("Present_in_both","homo_in_both","hetero_in_both","homo_in_WES_hetero_in_WGS","homo_in_WGS_hetero_in_WES","Ind")
	
	if(length(WES_files) != length(WGS_files)){
	
		stop("Different number of WES and WGS files")
	
	} else {
	
		for (i in 1 : length(WES_files)){
			
			exo=read.table(WES_files[i], sep="\t")
			gen=read.table(WGS_files[i], sep="\t")
		
			exo = exo[,c(1,2,3,4,5,10)]
			gen = gen[,c(1,2,3,4,5,10)]
		
			colnames(exo)=c("chr","pos","rs","ref","alt","geno1")
			colnames(gen)=c("chr","pos","rs","ref","alt","geno2")
			exo1 = exo
			gen1 = gen
			
			
			vcf_name = do.call(rbind,strsplit(do.call(rbind, strsplit(WGS_files[i],"/"))[2],"[.]"))[1]
			merge_table[i,6] = vcf_name
			
			
			


			v=strsplit(as.character(unlist(exo1$geno1)),":")
			v1 = suppressWarnings(do.call(rbind, v))
			
			if(dim(v1)[2]==4){
				exo1 = exo1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				p = strsplit(as.character(unlist(v1[,2])),",")
				p1 = suppressWarnings(do.call(rbind, p))
				exo1[,7:8]=v1[,1:2]
				exo1[,9] = as.numeric(p1[,1]) + as.numeric(p1[,2])
				exo1[,10]=v1[,3]
			}else{	
				exo1 = exo1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				exo1[,7:10]=v1[,1:4]
			}
			

			v=strsplit(as.character(unlist(gen1$geno2)),":")
			v1 = suppressWarnings(do.call(rbind, v))
			if(dim(v1)[2]==4){
				gen1 = gen1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4]),]
				p = strsplit(as.character(unlist(v1[,2])),",")
				p1 = suppressWarnings(do.call(rbind, p))
				gen1[,7:8]=v1[,1:2]
				gen1[,9] = as.numeric(p1[,1]) + as.numeric(p1[,2])
				gen1[,10]=v1[,3]
			}else{	
				gen1 = gen1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				v1 = v1[(v1[,1] != v1[,2])&(v1[,1] != v1[,3])&(v1[,1] != v1[,4])& (v1[,1] != v1[,5]),]
				gen1[,7:10]=v1[,1:4]
			}

			
			all_var_table[i,1] =  vcf_name
			all_var_table[i,2] = dim(exo1)[1]
			all_var_table[i,3] = dim(exo1[exo1[,7]!="0/1" & exo1[,7]!="1/1",])[1]
			all_var_table[i,4] = dim(exo1[exo1[,7]=="1/1",])[1]
			all_var_table[i,5] = dim(exo1[exo1[,7]=="0/1",])[1]
			all_var_table[i,6] = dim(gen1)[1]
			all_var_table[i,7] = dim(gen1[gen1[,7]!="0/1" & gen1[,7]!="1/1",])[1]
			all_var_table[i,8] = dim(gen1[gen1[,7]=="1/1",])[1]
			all_var_table[i,9] = dim(gen1[gen1[,7]=="0/1",])[1]

			
			
			
			exo1=exo1[exo1[,7]!="1/2" & exo1[,7]!="2/1" & exo1[,7]!="2/2",]
			
			gen1=gen1[gen1[,7]!="1/2" & gen1[,7]!="2/1" & gen1[,7]!="2/2",]
			
			
			m=merge(gen,exo, all.x=F, all.y=F)
			merge_table[i,1]=dim(m)[1]
		
		
		
		
			m_exo = strsplit(as.character(unlist(m$geno1)),":")
			m1_exo = suppressWarnings(do.call(rbind, m_exo))
			m[,8:11]=m1_exo[,1:4]
		
				
			m_gen = strsplit(as.character(unlist(m$geno2)),":")
			m1_gen = suppressWarnings(do.call(rbind, m_gen))
			m[,12:15]=m1_gen[,1:4]
			
			
			m = m[(m[,8] != m[,11])&(m[,12] != m[,15]) ,]
			m = m[(m[,8] == "0/1")|(m[,8] == "1/1") & (m[,12] == "0/1")|(m[,12] == "1/1"),]


			m1 = m[m[,8]=="1/1" & m[,12]=="1/1",]
			merge_table[i,2]=dim(m1)[1]
			m1 = m[m[,8]=="0/1" & m[,12]=="0/1",]
			merge_table[i,3]=dim(m1)[1]
			m1 = m[m[,8]=="1/1" & m[,12]=="0/1",]
			merge_table[i,4]=dim(m1)[1]
			m1 = m[m[,8]=="0/1" & m[,12]=="1/1",]
			merge_table[i,5]=dim(m1)[1]
		
		}
		
	write.table(merge_table, "Merged_variants.txt", sep="\t", row.names=F)
	write.table(all_var_table, "all_variant.txt", sep="\t", row.names=F)
	merge_table=data.frame(merge_table)
	mdata=melt(merge_table, id="Ind")
	colnames(mdata)=c("Ind","Genotypes","Number_of_variants")
	mdata$Number_of_variants=as.numeric(as.character(unlist(mdata$Number_of_variants)))
	g = ggplot(mdata, aes(x=Ind, y=Number_of_variants, fill=Genotypes, colour=Genotypes)) + geom_bar(stat="identity", position=position_dodge())  + theme_bw()
	ggsave("All_shared_variants_genotypes.pdf", g, useDingbats=FALSE,  width = 15, height =5)
	
	}

}





Make_graph_coverage<-function(){ # Analyse gene coverage in WES and WGS and makes graphes. GATK WES_coverage.sample_gene_summary and WGS_coverage.sample_gene_summary should be in Coverage/WES_coverage and Coverage/WGS_coverage respectively

	library(reshape)
	library(ggplot2)

	exo = read.table("Coverage/WES_coverage/WES_coverage.sample_gene_summary", sep="\t", h=T)
	gen = read.table("Coverage/WGS_coverage/WGS_coverage.sample_gene_summary", sep="\t", h=T)
	exo = exo[-1,]
	gen = gen[-1,]
	
	tmp = strsplit(colnames(exo),"_")
	tmp1 = suppressWarnings(do.call(rbind, tmp))
	exo1 = exo[,tmp1[,1]=="Gene" | (tmp1[,3]=="above" & tmp1[,4]=="8")]
	
	for( i in 2:dim(exo1)[2]){
		j = (6 * i) - 8
		colnames(exo1)[i] = tmp1[j,1]
		}
	
	tmp = strsplit(colnames(gen),"_")
	tmp1 = suppressWarnings(do.call(rbind, tmp))
	gen1 = gen[,tmp1[,1]=="Gene" | (tmp1[,4]=="above" & tmp1[,5]=="8")]
	
	for( i in 2:dim(gen1)[2]){
		j = (6 * i) - 8
		colnames(gen1)[i] = paste(tmp1[j,1], tmp1[j,2],sep="_")
		}
	
	exo1 = melt(exo1,id="Gene")
	gen1 = melt(gen1,id="Gene")
	
	colnames(exo1)=c("Gene","Ind","Percentage_coverage_above_8")
	colnames(gen1)=c("Gene","Ind","Percentage_coverage_above_8")
	
	

	g=ggplot(exo1, aes(as.numeric(Percentage_coverage_above_8), colour=Ind, group=Ind)) + stat_density(aes(y = ..count..), position="identity",geom="line") + scale_x_continuous(limits=c(90,100),breaks=seq(90,100,2), labels=c("90%","92%","94%","96%","98%","100%")) + scale_y_continuous(limits=c(0,10000),breaks=seq(2000,10000,2000)) + theme_bw()
	suppressWarnings(ggsave("WES_gene_coverage.pdf", g, useDingbats=FALSE))
			

	g=ggplot(gen1, aes(as.numeric(Percentage_coverage_above_8), colour=Ind, group=Ind)) + stat_density(aes(y = ..count..), position="identity",geom="line") + scale_x_continuous(limits=c(90,100),breaks=seq(90,100,2), labels=c("90%","92%","94%","96%","98%","100%")) + scale_y_continuous(limits=c(0,10000),breaks=seq(2000,10000,2000)) + theme_bw()
	suppressWarnings(ggsave("WGS_gene_coverage.pdf", g, useDingbats=FALSE))
	
	exo1[,4]="WES"
	gen1[,4]="WGS"

	merged = rbind(exo1, gen1)
	colnames(merged)[4]="Method"
	
	
	g=ggplot(merged, aes(as.numeric(Percentage_coverage_above_8), colour=Method, group=Method)) + stat_density(aes(y=..count..), position="identity",geom="line") + scale_x_continuous(limits=c(90,100),breaks=seq(90,100,2), labels=c("90%","92%","94%","96%","98%","100%")) + scale_y_continuous(limits=c(0,60000),breaks=seq(12000,60000,12000), labels=seq(2000,10000,2000)) + theme_bw()
	suppressWarnings(ggsave("merged_gene_coverage.pdf", g, useDingbats=FALSE))
	

}


Exons_statistics<-function(type='protein_coding', plus50 = False){ # Make statistics on protein coding exons and RNAs coveraed by the WES kit type={'protein_coding','lincRNA','miRNA','snoRNA'}. plus50 is a logical variable to check or not plus 50 bp flanking each side of each WES captured-region. Before runing the function, declare a new class "resul" : setClass(Class="result", representation(totaly="numeric", partialy="numeric",excluded="numeric",somme="numeric"))

library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
goids = getBM(attributes=c('chromosome_name','exon_chrom_start','exon_chrom_end','ensembl_gene_id','gene_biotype','external_gene_id'), filters='biotype',  values=type, mart=ensembl)
goids = goids[goids$chromosome_name=="1" | goids$chromosome_name=="2" | goids$chromosome_name=="3" | goids$chromosome_name=="4" | goids$chromosome_name=="5" | goids$chromosome_name=="6" | goids$chromosome_name=="7" | goids$chromosome_name=="8" | goids$chromosome_name=="9" | goids$chromosome_name=="10" | goids$chromosome_name=="11" | goids$chromosome_name=="12" | goids$chromosome_name=="13" | goids$chromosome_name=="14" | goids$chromosome_name=="15" | goids$chromosome_name=="16" | goids$chromosome_name=="17" | goids$chromosome_name=="18" | goids$chromosome_name=="19" | goids$chromosome_name=="20" | goids$chromosome_name=="21" | goids$chromosome_name=="22" | goids$chromosome_name=="X" | goids$chromosome_name=="Y" ,]
goids = goids[,1:3]
co = goids[,1]
co[goids[,1]=="X"] = "23"
co[goids[,1]=="Y"] = "24"
goids[,1]=co

colnames(goids) = c("Chr","Begin","End")




goids1 = matrix(0, ncol = 3, nrow=dim(goids)[1])
colnames(goids1) = colnames(goids)
goids = goids[order(goids[,1],goids[,2],goids[,3]),]

i = 1
l=1
while ( i < dim(goids)[1]){
     chr = goids$Chr[i]
     beg = goids$Begin[i]
     end = goids$End[i]
     for( j in (i+1):dim(goids)[1]){
	     if(goids$Chr[j] == chr & goids$Begin[j] == beg){
            i = i + 1
			break} else if((goids$Chr[j] == chr & goids$Begin[j] > beg)| goids$Chr[j] != chr){
                 goids1[l,] = c(chr,beg,end)
				 l = l + 1
				 if((goids$Chr[j] == chr & goids[j,3]> end)| goids$Chr[j] != chr) {
				 i = j 
				 break} 
}
}
}


goids1 = goids1[1:(which(goids1[,1]==0)[1]-1),]






interval = read.table("/resources/SureSelect-v4plusUTR_baits.interval_list")

v=strsplit(as.character(unlist(interval[,1])),":")
v1 = suppressWarnings(do.call(rbind, v))
v=strsplit(as.character(unlist(v1[,2])),"-")
v2 = suppressWarnings(do.call(rbind, v))

v1[,2] = v1[,1]
v1[v1[,1]=="X",2] = "23"
v1[v1[,1]=="Y",2] = "24"


interval[,1] = v1[,2]
interval[,2] = v2[,1]
interval[,3] = v2[,2]



interval[,1] = as.numeric(interval[,1])
interval[,2] = as.numeric(interval[,2])
interval[,3] = as.numeric(interval[,3])

if(plus50){
interval[,2] = interval[,2] - 50
interval[,3] = interval[,3] + 50
} 

colnames(interval)=c("Chr","Begin","End")


interval1 = matrix(0, ncol = 3, nrow=dim(interval)[1])

j=1
sw=0
i = 1
while ( i < dim(interval)[1]){
if(sw == 0){
     chr = interval$Chr[i]
     beg = interval$Begin[i]
     end = interval$End[i]
     sw = 1}else{
         if(interval$Chr[i+1] == chr & interval$Begin[i+1] <=end){
             end = max(interval$End[i+1],end)
			i = i + 1} else {
                interval1[j,] = c(chr,beg,end)
	          j = j + 1
                 sw = 0
		 i = i + 1}
     }
 }
interval1 = interval1[1:(which(interval1[,1]==0)[1]-1),]







exo = matrix(0,ncol=4, nrow=((2 * dim(interval1)[1]) + (2 * dim(goids1)[1])))
exo[1:(dim(goids1)[1]),1] = goids1[,2]
exo[(dim(goids1)[1] + 1) : (2 * dim(goids1)[1]),1] = goids1[,3]
exo[((2 * dim(goids1)[1]) + 1) : ((2 * dim(goids1)[1]) + dim(interval1)[1]),1] = interval1[,2]
exo[((2 * dim(goids1)[1]) + dim(interval1)[1] + 1) : ((2 * dim(interval1)[1]) + (2 * dim(goids1)[1])),1] = interval1[,3]



exo[1:(dim(goids1)[1]),2] = 1
exo[(dim(goids1)[1] + 1) : (2 * dim(goids1)[1]),2] = 2
exo[((2 * dim(goids1)[1]) + 1) : ((2 * dim(goids1)[1]) + dim(interval1)[1]),2] = 3
exo[((2 * dim(goids1)[1]) + dim(interval1)[1] + 1) : ((2 * dim(interval1)[1]) + (2 * dim(goids1)[1])),2] = 4



exo[1:(dim(goids1)[1]),3] = goids1[,1]
exo[(dim(goids1)[1] + 1) : (2 * dim(goids1)[1]),3] = goids1[,1]
exo[((2 * dim(goids1)[1]) + 1) : ((2 * dim(goids1)[1]) + dim(interval1)[1]),3] = interval1[,1]
exo[((2 * dim(goids1)[1]) + dim(interval1)[1] + 1) : ((2 * dim(interval1)[1]) + (2 * dim(goids1)[1])),3] = interval1[,1]


exo[1:(dim(goids1)[1]),4] = c(1:(dim(goids1)[1]))
exo[(dim(goids1)[1] + 1) : (2 * dim(goids1)[1]),4] = c(1:(dim(goids1)[1]))

exo = data.frame(exo)
for (i in 1:4){
exo[,i]=as.numeric(as.character((exo[,i])))
}

exo1 = exo[order(exo[,3],exo[,1],exo[,2]),]


sw = 0
excluded = 0
partialy = 0
totaly=0


for (i in 1:24){
  
   if(i==1){r=0}else{r=j+1}   
   for (j in (r+1):((dim(exo1[exo1[,3]==i,])[1]-1)+r)){
   if(exo1[j,3] == i & exo1[j,2]==3){sw=1}
   if(exo1[j,3] == i & exo1[j,2]==4){sw=0}
   if(exo1[j,3] == i & exo1[j,2]==1){
     code = exo1[j,4]  
     for(l in j+1:((dim(exo1[exo1[,3]==i,])[1])+r)){
       if(exo1[l,3] == i & (exo1[l,2]==3 | exo1[l,2]==4)){partialy=partialy + 1
       break}
       if(exo1[l,3] == i & exo1[l,2]==2 & code == exo1[l,4] & sw == 0){
         excluded = excluded + 1
         break}
       if(exo1[l,3] == i & exo1[l,2]==2 & code == exo1[l,4] & sw == 1){
         totaly = totaly + 1
         break}
    }
   }
   }
}
     
	 
	 totaly
	 partialy
	 excluded
	 return(new("result",
          totaly=totaly,
          partialy=partialy,
		  excluded=excluded,
somme=sum(totaly,partialy,excluded)))
	 

}
