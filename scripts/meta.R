#################################################################
# Function: WGS meta Figure
# Call: Rscript meta.R -i abund_file -o outfile
# R packages used: optparse, reshape2, ggplot2,RColorBrewer, grDevices
# Last update: 2021-05-10, Zhang Wen
###Usage：Rscript D:\E盘\Project\中国健康人基线数据\2020WGS\bin\R画图\meta.R -i D:\E盘\Project\中国健康人基线数据\2020WGS\bin\R画图\Genus.txt -m D:\E盘\Project\中国健康人基线数据\2020WGS\bin\R画图\example.meta.txt
# ###Bug1:画alpha图时显著性未展示
#################################################################
# install necessary libraries
p <- c("optparse","reshape2","ggplot2","RColorBrewer","grid","scales","vegan","agricolae","gridExtra","dplyr","ggrepel","gggenes","ggsignif","pheatmap","reshape2","picante","corrplot","ape","multcomp","patchwork","factoextra")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
## clean R environment
rm(list = ls())
setwd('./')
## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--abund_file"), type="character", help="Input feature table with relative abundance (*.Abd) [Required]"),
  make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Optional]"),
  make_option(c("-o", "--out_dir"), type="character", default='MetaFigure', help="Output directory [default %default]"),
  make_option(c("-p", "--prefix"), type="character", default='Out', help="Output file prefix [default %default]"),
  make_option(c("-t", "--threshold"), type="double", default=0.01, help="Average value threshold [Optional, default %default]"),
  make_option(c("--cutoff_positivesample", "-c"), type="numeric", default=0, help="cutoff for positive sample number  [Optional, default %default]"),
  make_option(c("--cutoff_readpercentage", "-r"), type="numeric", default=0, help=" cutoff for read number/percentage  [Optional, default %default]"),
  make_option(c("--overturn", "-f"), type="character", default=F, help=" Each row for a sample (Default); -turn T Each column for a sample  [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$abund_file)) stop('Please input a feature table (*.Abd)')
# load data
matrixfile <- opts$abund_file
mapfile <- opts$meta_data
ave_t <- opts$threshold
outpath <- opts$out_dir
cut_i<-opts$cutoff_positivesample
cut_j<-opts$cutoff_readpercentage
turn<-opts$overturn
dir.create(outpath)
			library("ggplot2") # load related packages
			library("grid")
			library("scales")
			library("vegan")
			library("agricolae")
			library("gridExtra")
			library("dplyr")
			library("ggrepel")
			library("gggenes")
			library("ggsignif")
			library("pheatmap")
			library(reshape2)
			library(picante)
			library(corrplot)
			library(ape)
			library(multcomp)
			library(patchwork)
			library(factoextra)
			library(RColorBrewer)

#------------------------------------------------------------------------------------
#####画分布图Distribution
otu <- read.table(matrixfile,header = T, row.names = 1,sep="\t")

disbar <-as.data.frame(otu/rowSums(otu))
disbar <- t(disbar)

disbar <- disbar[names(sort(rowSums(disbar),decreasing = T)),]

 Unclassified_other <- disbar[which(apply(disbar,1,mean) <= ave_t),]

invisible(if (sum( Unclassified_other) ==0 ) ( Unclassified_big <- disbar))
invisible(if (sum( Unclassified_other) !=0 ) ( Unclassified_big <- disbar[-(which(apply(disbar,1,mean) <= ave_t)),]))

widforpdf <- ncol(disbar)
 Unclassified_other <- as.matrix( Unclassified_other)
if (dim( Unclassified_other)[2] ==1 )  Unclassified_other <- t( Unclassified_other)

if (is.null( Unclassified_other)==F) {
  disbar <- rbind( Unclassified_big,"Other"=c(colSums( Unclassified_other)),deparse.level = 2)
}

if (mean(colSums(disbar))<0.9999) {                         #Complete to 100%
  if (rownames(disbar)[nrow(disbar)]=="Other") {
    disbar[nrow(disbar),] <- disbar[nrow(disbar),]+(1-colSums(disbar))
  }
  else {
    disbar <- rbind(disbar,"other"=(sapply((1-colSums(disbar)),function(x)max(x,0))),deparse.level = 2)
  }
}
colours <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),sample(rainbow(length(colnames(t(disbar)))),length(colnames(t(disbar)))))
meltdata <- melt(abs(data.matrix(t(disbar))),varnames=c("Samples","Cutline"),value.name="Relative_Abundance")

if(turn == T){disbar=t(disbar);}###翻转矩阵
#-----------------------------------------------------------------------------------------

data_map <- read.table(mapfile,header = T, row.names= 1,sep="\t",na.strings = "NA")

pp<-ggplot(meltdata,aes(x=Samples,y=Relative_Abundance,fill=Cutline))+
    geom_bar(stat='identity')+ ylab("Relative Abundance")+ xlab("Samples")+
    scale_x_discrete(limits=c(colnames(disbar)))+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","25%","50%","75%","100%"))+
    guides(fill = guide_legend(reverse=TRUE ,ncol = (ceiling(nrow(disbar)/35))))+
    scale_fill_manual (values=colours) +
    theme(legend.position="right",axis.text.x=element_text(size=12,colour="black",angle=90,vjust=0.5),
          axis.text.y=element_text(size=12,colour="black"), axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),panel.grid.major=element_line(colour=NA))
  suppressWarnings(ggsave(paste(outpath, "/", opts$prefix, ".distribution.pdf",sep=""),plot=pp,width=ceiling(16+widforpdf/8),height=10, limitsize=FALSE))


for (i in 1:ncol(data_map)) {
  
    tempframe <- data.frame(meltdata,dv=data_map[,i])
    pp<-ggplot(tempframe,aes(x=Samples,y=Relative_Abundance,fill=Cutline))+
      geom_bar(stat = "identity")+ ylab("Relative Abundance")+ xlab("Samples")+
      #  scale_x_discrete(limits=c(colnames(disbar)))+
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","25%","50%","75%","100%"))+
      guides(fill = guide_legend(reverse=TRUE ,ncol = (ceiling(nrow(disbar)/35))))+
      scale_fill_manual (values=colours) +
      facet_grid(~ dv,scales="free_x",space="free_x")+ 
      theme(legend.position="right",axis.text.x=element_text(size=12,colour="black",angle=90,vjust=0.5),
            axis.text.y=element_text(size=12,colour="black"), axis.title.x=element_text(size=16),
            axis.title.y=element_text(size=16),panel.grid.major=element_line(colour=NA))
    suppressWarnings(ggsave(paste(outpath,"/", opts$prefix, ".", colnames(data_map)[i],".distribution.pdf",sep=""),plot=pp,height=10,width=ceiling(16+widforpdf/8),limitsize=FALSE))
  }


 
 ##构建函数
alpha <- function(x, tree = NULL, base = exp(1)) {
       est <- estimateR(x)
       Richness <- est[1, ]
       Chao1 <- est[2, ]
       ACE <- est[4, ]
       Shannon <- diversity(x, index = 'shannon', base = base)
       Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
       Pielou <- Shannon / log(Richness, base)
       goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
       
       result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
       if (!is.null(tree)) {
              PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
              names(PD_whole_tree) <- 'PD_whole_tree'
              result <- cbind(result, PD_whole_tree)
       }
       result
}


####
#data_map <- read.table(mapfile,header = T, row.names= 1,sep="\t")
otu <- read.table(matrixfile,header = T, row.names = 1,sep="\t")
			#去除一些低丰度OTU，在至少i个样本中具有至少j个序列
			#otu_its_t <- otu[which(rowSums(otu >0) >= 1),]
			#otu <- t(otu)
			#otu <- otu[,which((otu > cut_j) > cut_i)]
			otu_its_t=as.data.frame(otu/rowSums(otu)*100)
			otu_its_t <- otu_its_t[,which(colSums(otu_its_t > cut_j) > cut_i)] 
			otuf<-paste(outpath,"/","otu_filter.csv", sep="")
			write.csv(otu_its_t, otuf, quote = FALSE)
			otu_final<-floor(otu_its_t*1000000)
			
			alpha_all <- alpha(otu_final)
			p<-paste(outpath,"/alpha.csv", sep="")
			write.csv(alpha_all, p, quote = FALSE)
for (i in 1:ncol(data_map)) {
	i
	group=factor(data_map[,i])
	group<-group[!duplicated(group)]
	group<-na.omit(group)
	c<-length(group)
	#compaired <- list( c(group[1],group[2]))
	info <- data_map

			cmp<-merge(info,alpha_all,  by = "row.names", all = TRUE)
			
			#cmp<-na.omit(cmp)
	if(c==2){
	compaired <- list(c(group[1],group[2]))
		cmp1<-subset(cmp,cmp[,i+1]==group[1])
			cmp2<-subset(cmp,cmp[,i+1]==group[2])
			cmp=rbind(cmp1,cmp2)
	}
	if(c==3){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]))
		cmp1<-subset(cmp,cmp[,i+1]==group[1])
			cmp2<-subset(cmp,cmp[,i+1]==group[2])
			cmp3<-subset(cmp,cmp[,i+1]==group[3])
			cmp=rbind(cmp1,cmp2,cmp3)
	}
	if(c==4){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]))
			cmp1<-subset(cmp,cmp[,i+1]==group[1])
			cmp2<-subset(cmp,cmp[,i+1]==group[2])
			cmp3<-subset(cmp,cmp[,i+1]==group[3])
			cmp4<-subset(cmp,cmp[,i+1]==group[4])
			cmp=rbind(cmp1,cmp2,cmp3,cmp4)
	}
	if(c==5){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]),c(group[1],group[5]),c(group[2],group[5]),c(group[3],group[5]),c(group[4],group[5]))
		cmp1<-subset(cmp,cmp[,i+1]==group[1])
			cmp2<-subset(cmp,cmp[,i+1]==group[2])
			cmp3<-subset(cmp,cmp[,i+1]==group[3])
			cmp4<-subset(cmp,cmp[,i+1]==group[4])
			cmp5<-subset(cmp,cmp[,i+1]==group[5])
			cmp=rbind(cmp1,cmp2,cmp3,cmp4,cmp5)
	}
		
	
			###计算
			
			
			
			p<-paste(outpath,"/",colnames(data_map)[i],"_otu_alpha.csv", sep="")
			write.csv(cmp, p, quote = FALSE)
			t<-factor(cmp[,i+1])
			t<-na.omit(t)
			

			richness = ggplot(cmp, aes(x=factor(t), y=Richness,color=t))+
			  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
			  labs(x="Groups", y="Richness")+theme(legend.position='none') +
			      geom_signif(comparisons = compaired,
					step_increase = 0.3,
					map_signif_level = F,
					test = wilcox.test)


			  ace = ggplot(cmp, aes(x=factor(t), y=ACE,color=t))+
			  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
			  labs(x="Groups", y="ACE1")+theme(legend.position='none') +
			      geom_signif(comparisons = compaired,
					step_increase = 0.3,
					map_signif_level = F,
					test = wilcox.test)

			shannon = ggplot(cmp, aes(x=factor(t), y=Shannon,color=t))+
			  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
			  labs(x="Groups", y="Shannon")+theme(legend.position='none') +
			      geom_signif(comparisons = compaired,
					step_increase = 0.3,
					map_signif_level = F,
					test = wilcox.test)

			pielou = ggplot(cmp, aes(x=factor(t), y=Pielou,color=t))+
			  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
			  labs(x="Groups", y="Pielou")+theme(legend.position='none') +
			      geom_signif(comparisons = compaired,
					step_increase = 0.3,
					map_signif_level = F,
					test = wilcox.test)

			simpson = ggplot(cmp, aes(x=factor(t), y=Simpson,color=t))+
			  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
			  labs(x="Groups", y="Simpson")+theme(legend.position='none') +
			      geom_signif(comparisons = compaired,
					step_increase = 0.3,
					map_signif_level = F,
					test = wilcox.test)

			chao = ggplot(cmp, aes(x=factor(t), y=Chao1,color=t))+
			  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
			  labs(x="Groups", y="Chao1")+theme(legend.position='none') +
			      geom_signif(comparisons = compaired,
					step_increase = 0.3,
					map_signif_level = F,
					test = wilcox.test)

			com= grid.arrange(richness,ace,shannon,simpson,pielou,chao,ncol=2)
			ggsave(paste(outpath,"/Group_",colnames(data_map)[i],".OTU_sum.pdf", sep=""), com, width = 10, height = 18)
			#ggsave(paste(outpath,"/",group1,"_vs_",group2,".OTU_richness.pdf", sep=""), richness, width = 5, height = 6)
			#ggsave(paste(outpath,"/",group1,"_vs_",group2,".OTU_ace.pdf", sep=""), ace, width = 5, height = 6)
			#ggsave(paste(outpath,"/",group1,"_vs_",group2,".OTU_shannon.pdf", sep=""), shannon, width = 5, height = 6)
			#ggsave(paste(outpath,"/",group1,"_vs_",group2,".OTU_simpson.pdf", sep=""), simpson, width = 5, height = 6)
			#ggsave(paste(outpath,"/",group1,"_vs_",group2,".OTU_pielou.pdf", sep=""), pielou, width = 5, height = 6)
			#ggsave(paste(outpath,"/",group1,"_vs_",group2,".OTU_chao.pdf", sep=""), chao, width = 5, height = 6)
		
	
}
######画bray distance热力图

		otubray=as.data.frame(otu_its_t/rowSums(otu_its_t)*100)

		bray_curtis = vegdist(otubray, method = "bray")

		bray_curtis= as.matrix(bray_curtis)
		dis=bray_curtis

		#bray_curtis
		p<-paste(outpath,"/","bray_curtis.csv", sep="")
		write.csv(bray_curtis, p, quote = FALSE)
		# 准备行/列注释
		for (i in 1:ncol(data_map)) {
			i
		anno_col = data.frame(Group =data_map[,i], row.names = rownames(data_map))
		anno_row = data.frame(Group =data_map[,i], row.names = rownames(data_map))


		# 绘图dist Heatmap
		pheatmap(dis,
		  treeheight_col=15,treeheight_row=15,
		  annotation_names_row= T,annotation_names_col=T,
		  annotation_col = anno_col,annotation_row = anno_row,
		  filename = paste(outpath,"/Group_",colnames(data_map)[i],".heatmap.Bray-Curtis.jpg", sep=""),width=30, height=40,
		  fontsize=7,display_numbers=T)
		pheatmap(dis,
		  treeheight_col=15,treeheight_row=15,
		  annotation_names_row= T,annotation_names_col=T,
		  annotation_col = anno_col,annotation_row = anno_row,
		  filename = paste(outpath,"/Group_",colnames(data_map)[i],".heatmap.Bray-Curtis.pdf", sep=""),width=30, height=40,
		  fontsize=7,display_numbers=F)
		 }


#####绘Heatmap

	for (i in 1:ncol(data_map)) {
		anno_row = data.frame(Group = data_map[,i], row.names = rownames(data_map))
		scale_test <- apply(otu_its_t, 2, function(x) log2(x+1))



		pheatmap(mat=otu_its_t,
		  treeheight_col=15,treeheight_row=15,
		  annotation_names_row= T,
		 annotation_row = anno_row,
		  filename = paste(outpath,"/Group_",colnames(data_map)[i],".heatmap.sample.all.pdf", sep=""),width=20, height=40)

		pheatmap(mat=scale_test , treeheight_col=15,treeheight_row=15,
		  annotation_names_row= T,
		 annotation_row = anno_row,
			 filename=paste(outpath,"/Group_",colnames(data_map)[i],".heatmap.sample.log2.pdf", sep=""),width=20, height=40)
	}
####组内差异和组间差异比较
#生成group矩阵 
a=paste(outpath,"/bray_curtis.csv", sep="")
b=paste(outpath,"/","group.csv", sep="")
		write.csv(data_map, b, quote = FALSE)		
		a
		b
#data_map <- read.table(mapfile,header = T, row.names= 1,sep="\t")
c<-ncol(data_map)
c
for (i in 1:c) {
	i
	command<-paste("matrixtodis.pl ",a," ",b," ",i)

	system2('perl', command)
	dis <- read.table('dis.csv',header = T, row.names = 1,sep=",")
	p<-paste(outpath,"/",colnames(data_map)[i],"_bray_curtis.csv", sep="")
	write.csv(dis, p, quote = FALSE)
	t<-factor(dis[,5])
	tt<-t[!duplicated(t)]
	group<-na.omit(tt)
	c<-length(group)
	#compaired <- list( c(group[1],group[2]))
	if(c==2){
	compaired <- list(c(group[1],group[2]))
	}
	if(c==3){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]))
	}
	if(c==4){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]))
	}
	if(c==5){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]),c(group[1],group[5]),c(group[2],group[5]),c(group[3],group[5]),c(group[4],group[5]))
	}
	
				compaired 
	group = ggplot(dis, aes(x=factor(t), y=Bray,color=t))+ geom_violin(aes(fill = t),trim=FALSE) +
				  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.2, fill="transparent") + 
				  labs(x="Groups", y="Bray-Curtis distance")+theme(legend.position='none') +
				      geom_signif(comparisons = compaired,
						step_increase = 0.3,
						map_signif_level = T,
						test = wilcox.test)

	ggsave(paste(outpath,"/GroupType_",colnames(data_map)[i],"_Bray_dis.pdf", sep=""), group, width = 5, height = 6)
	group = ggplot(dis, aes(x=factor(t), y=Bray,color=t))+ geom_violin(aes(fill = t),trim=FALSE) +
				  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.2, fill="transparent") + 
				  labs(x="Groups", y="Bray-Curtis distance")+theme(legend.position='none')

	ggsave(paste(outpath,"/GroupType_",colnames(data_map)[i],"_Bray_dis_v2.pdf", sep=""), group, width = 5, height = 6)


	dis <- read.table('dis_v2.csv',header = T, row.names = 1,sep=",")
	p<-paste(outpath,"/",colnames(data_map)[i],"_bray_curtis_v2.csv", sep="")
	write.csv(dis, p, quote = FALSE)
	t<-factor(dis[,5])
	tt<-t[!duplicated(t)]
	group<-na.omit(tt)
	c<-length(group)
	#compaired <- list( c(group[1],group[2]))
	if(c==2){
	compaired <- list(c(group[1],group[2]))
	}
	if(c==3){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]))
	}
	if(c==4){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]))
	}
	if(c==5){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]),c(group[1],group[5]),c(group[2],group[5]),c(group[3],group[5]),c(group[4],group[5]))
	}
	
				compaired 
	group = ggplot(dis, aes(x=factor(t), y=Bray,color=t))+ geom_violin(aes(fill = t),trim=FALSE) +
				  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.2, fill="transparent") + 
				  labs(x="Groups", y="Bray-Curtis distance")+theme(legend.position='none') +
				      geom_signif(comparisons = compaired,
						step_increase = 0.3,
						map_signif_level = T,
						test = wilcox.test)

	ggsave(paste(outpath,"/GroupType_",colnames(data_map)[i],"_Bray_dis_v3.pdf", sep=""), group, width = 5, height = 6)
	group = ggplot(dis, aes(x=factor(t), y=Bray,color=t))+ geom_violin(aes(fill = t),trim=FALSE) +
				  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.2, fill="transparent") + 
				  labs(x="Groups", y="Bray-Curtis distance")+theme(legend.position='none')

	ggsave(paste(outpath,"/GroupType_",colnames(data_map)[i],"_Bray_dis_v4.pdf", sep=""), group, width = 5, height = 6)
}

###不同因素的影响 报错 未运行成功 主要问题存在于格式转化
mt<-function(otu,env){
  library(vegan)
  library(dplyr)
  vars <- colnames(env)
  models<-list()
  for (i in seq_along(vars)){
    otu_bray<-vegdist(otu,method = "bray")
    otu_bray= as.matrix(otu_bray)
    env_dis<-vegdist(env[vars[i]],method = "euclidean")
    env_dis= as.matrix(env_dis)
    model <- mantel(otu_bray,env_dis, permutations=999)
    name <- vars[i]
    statistic <- model$statistic
    signif <- model$signif
    models[[i]] <- data.frame(name = name, statistic = statistic, signif = signif, row.names = NULL)
  }
  models %>%  bind_rows()
}


command<-paste("groupformat.pl ",mapfile)

	system2('perl', command)
group <- read.table('group_format.txt',header = T, row.names = 1,sep=",")
group<-group[-(ncol(group))]
mantRpTotal<-mt(otu_its_t,group)
b=paste(outpath,"/","mantRpTotal.csv", sep="")
write.csv(mantRpTotal, b, quote = FALSE) ###不同生活因素的影响率


#


###PCA图

data <- read.table(matrixfile, header=T, row.names=NULL,sep="\t")



rownames_data <- make.names(data[,1],unique=T)
data <- data[,-1,drop=F]
rownames(data) <- rownames_data
data <- data[,which(colSums(data > cut_j) > cut_i)] 
#data <- data[rowSums(data)>0,]
# 去掉方差为0 的行，这些本身没有意义，也妨碍后续运算
data <- data[apply(data, 1, var)!=0,]
mads <- apply(data, 1, mad)
data <- data[rev(order(mads)),]
# 此处未选
# data <- data[1:500,]
dim(data)

data_t <- t(data)
variableL <- ncol(data_t)


  sample <- data_map
  for (i in 1:ncol(data_map)) {
  data_t_m <- cbind(data, sample[,i])
  rownames(data_t_m) <- data_t_m$Row.names
  data_t <- data_t_m[,-1]

# By default, prcomp will centralized the data using mean.
# Normalize data for PCA by dividing each data by column standard deviation.
# Often, we would normalize data.
# Only when we care about the real number changes other than the trends,
# `scale` can be set to TRUE.
# We will show the differences of scaling and un-scaling effects.
pca <- prcomp(data_t[,1:variableL])
p1=fviz_eig(pca, addlabels = TRUE)

p2=fviz_pca_ind(pca, col.ind=sample[,i], mean.point=F, addEllipses = T, legend.title="Groups")

ggsave(paste(outpath,"/GroupType_",colnames(data_map)[i],"_PCA.pdf", sep=""), p2, width = 5, height = 6)
}

###画PcoA图
for (i in 1:ncol(data_map)) {
	i
	group=factor(data_map[,i])
	group<-group[!duplicated(group)]
	group<-na.omit(group)
	c<-length(group)
	#compaired <- list( c(group[1],group[2]))
		groups<-data_map
		cmp<-cbind(otu_its_t,groups[,i])	
			
			#cmp<-na.omit(cmp)
	if(c==2){
	compaired <- c(group[1],group[2])
		cmp1<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[1])
			cmp2<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[2])
			cmp=rbind(cmp1,cmp2)
	}
	if(c==3){
		compaired <- c(group[1],group[2],group[3])
		cmp1<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[1])
			cmp2<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[2])
			cmp3<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[3])
			cmp=rbind(cmp1,cmp2,cmp3)
		
	}
	if(c==4){
		compaired <- c(group[1],group[2],group[3],group[4])
			cmp1<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[1])
			cmp2<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[2])
			cmp3<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[3])
			cmp4<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[4])
			cmp=rbind(cmp1,cmp2,cmp3,cmp4)
	}
	if(c>=5){
		compaired <- c(group[1],group[2],group[3],group[4],group[5])
		cmp1<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[1])
			cmp2<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[2])
			cmp3<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[3])
			cmp4<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[4])
			cmp5<-subset(cmp,cmp[,(1+ncol(otu_its_t))]==group[5])
			cmp=rbind(cmp1,cmp2,cmp3,cmp4,cmp5)
		
	}
	
			###计算
			
			t<-ncol(cmp)-1

data<-data.frame(cmp[,-1])
data<-data.frame(data[,1:t-1])
rownames(data) <- c(1:nrow(data))
groups <- data.frame(rownames(cmp),cmp[,t+1])
colnames(groups) <- c("V1","V2")

length=nrow(groups)
times1=length%/%8
res1=length%%8
times2=length%/%5
res2=length%%5
col1=rep(1:8,times1)
col=c(col1,1:res1)
pich1=rep(c(21:24),times2)
pich=c(pich1,15:(15+res2))
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",                "#ADD1E5")

p<-vegdist(data, na.rm=TRUE)
pcoa<- pcoa(p, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$V2)
colnames(plotdata) <-c("sample","PC1","PC2","Group")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)

plotdata$Group <- factor(plotdata$Group,levels =group[!duplicated(group)])


yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~Group,data = plotdata)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
test$Group <- factor(test$Group,levels =group[!duplicated(group)])

p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 3,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=10,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")

p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 3,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=10,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")
p2<-ggplot(plotdata, aes(PC1, PC2)) +
xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  guides(fill = guide_legend(ncol = 1))+

   stat_ellipse(geom = "polygon",
               aes(fill=Group),
               alpha=0.3)+
  scale_fill_manual(values=cbbPalette,name = "Group") +geom_point(aes(fill=Group),size=3,pch = 21)

  otu.adonis=adonis(data~V2,data = groups,distance = "bray")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],  "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),  "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
            size = 3) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)

  p5
  ggsave(paste(outpath,"/GroupType_",colnames(data_map)[i],"_PCoA.pdf", sep=""), p5, width = 10, height = 12)
  		
	
}


####各种属依次比较
outpath2<-paste(outpath,"/EachGenus", sep="")
dir.create(outpath2)
data<-otu_its_t
 for (i in 1:ncol(data_map)) {
	t<-factor(data_map[,i])
	tt<-t[!duplicated(t)]
	group<-na.omit(tt)
	c<-length(group)
	#compaired <- list( c(group[1],group[2]))
	if(c==2){
	compaired <- list(c(group[1],group[2]))
	}
	if(c==3){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]))
	}
	if(c==4){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]))
	}
	if(c==5){
		compaired <- list(c(group[1],group[2]),c(group[1],group[3]),c(group[2],group[3]),c(group[1],group[4]),c(group[2],group[4]),c(group[3],group[4]),c(group[1],group[5]),c(group[2],group[5]),c(group[3],group[5]),c(group[4],group[5]))
	}
	for (ii in 1:ncol(data)) {
	group = ggplot(data, aes(x=factor(t), y=data[,ii],color=t))+ geom_violin(aes(fill = t),trim=FALSE) +
				  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.2, fill="transparent") + 
				  labs(x="Groups", y=colnames(data)[ii])+theme(legend.position='none') +
				      geom_signif(comparisons = compaired,
						step_increase = 0.1,
						map_signif_level = T,
						test = wilcox.test)
						
		ggsave(paste(outpath2,"/",colnames(data_map)[i],"_",colnames(data)[ii],".pdf", sep=""), group, width = 5, height = 6)
		ggsave(paste(outpath2,"/",colnames(data_map)[i],"_",colnames(data)[ii],".jpg", sep=""), group, width = 5, height = 6)
	}
 }