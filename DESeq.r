########################################################################################
##Usage: Rscript + DESeq_script_yiduclound + inputdata+output_dir+Case_num+Control_num+Pvalue+logFC
##注意修改conditions<-factor(rep(c("normal","case"),c(2,2)))相关参数
########################################################################################

library(stringr)
args<-commandArgs(TRUE)
input_data<-as.character(args[1])
output_dir<-as.character(args[2])
Case_num<-as.numeric(args[3])
Control_num<-as.numeric(args[4])
Pvalue<-as.numeric(args[5])
logFC<-as.numeric(args[6])

setwd(output_dir)

rawdata <-read.table(file=input_data,header=T,sep="\t",row.names=1)
##rawdata <- rawdata[-((dim(rawdata)[1]-5):dim(rawdata)[1]),]
rawdata <-rawdata[which(rowSums(rawdata) > 0),]

rawdata <-t(rawdata)

library(gmodels)
library(ggplot2)
data.pca<-fast.prcomp(rawdata)
a<-summary(data.pca)
tmp<-a$importance
pro1<-as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2<-as.numeric(sprintf("%.3f",tmp[2,2]))*100

pc=as.data.frame(a$x)
pc$names=rownames(pc)

##########################################

##########################################
level<-c(rep("Control",Control_num),rep("Case",Case_num))
pc$group=level
write.csv(pc,file="PCA_data.csv",row.names=F)

pdf(file="PCA_graph.pdf",width=9,height=9,onefile=FALSE)
pca=ggplot(pc,aes(PC1,PC2))+geom_point(size=3,aes(shape=group,color=group))+geom_text(size=2,aes(PC1,PC2,label=names))
xlab=paste("PC1(",pro1,"%)",sep="")
ylab=paste("PC2(",pro2,"%)",sep="")
pca2=pca+labs(x=xlab,y=ylab,title="PCA")+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10))
pca2+geom_hline(yintercept=0,linetype=4,color="grey")+geom_vline(xintercept=0,linetype=4,color="grey")
dev.off()

rawdata <-read.table(file=input_data,header=T,sep="\t",row.names=1)
##rawdata <- rawdata[-((dim(rawdata)[1]-5):dim(rawdata)[1]),]
rawdata <-rawdata <-rawdata[which(rowSums(rawdata) > 0),]

library("DESeq2")
conditions<-factor(rep(c("Control","Case"),c(Control_num,Case_num)),level=c("Control","Case"))
colData<-as.data.frame(conditions)
rownames(colData)<-colnames(rawdata)
dds0 <- DESeqDataSetFromMatrix(rawdata,colData,~ conditions)
dds0<-DESeq(dds0,betaPrior=FALSE)
res<-results(dds0,contrast=c("conditions","Case","Control"),alpha=0.9999)
res <-cbind(ID=rownames(res),as.data.frame(res))

res<-subset(res,padj<0.9999)
res<-res[order(res$padj),]
res<-na.omit(res)
write.csv(res,file="all_gene_DESeq.csv",row.names=F)

pdf(file="grap_MA.pdf",width=12,height=10,onefile=FALSE)
data=data.frame(results(dds0,contrast=c("conditions","Case","Control"),alpha=0.9999))
data=cbind(ID=rownames(data),as.data.frame(data))
data$colour <- ifelse(data$pvalue < Pvalue & abs(data$log2FoldChange)>= logFC & data$baseMean >= 1,ifelse(data$log2FoldChange> logFC,'red','blue'),'gray')
LogPvalue<--(log(Pvalue,10))
color <- c(red = "red",gray = "gray",blue = "blue")
P<-ggplot (data,aes(y=log2FoldChange,x=-log10(pvalue),colour=colour))+ylim(-10,10)+xlim(0,20)+scale_color_manual(values=color)+geom_hline(yintercept=logFC,linetype=3)+geom_hline(yintercept=-logFC,linetype=3)+geom_vline(xintercept=LogPvalue,linetype=4)+xlab("log2fold change")+xlab("-log10pvalue")+geom_point()
P<-P+theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
print(P)
dev.off()


data_DEG<-subset(data,(log2FoldChange >=logFC | log2FoldChange <=-logFC) & ( pvalue < Pvalue ) & (data$baseMean >= 1) )
write.csv(data_DEG,file="DEG_gene.csv",row.names=F,quote=F)

library("vsn")
library("pheatmap")

res_de <- subset(data,( data$pvalue < Pvalue ) & (data$baseMean >= 1) )
res_de_up <- subset(res_de, res_de$log2FoldChange>=logFC)
res_de_up_id <- as.vector(res_de_up$ID)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=-logFC)
res_de_dw_id <- as.vector(res_de_dw$ID)
red_de_all <- c(res_de_up_id, res_de_dw_id)

normalized_counts <- counts(dds0, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
red_de_expr <- normalized_counts[rownames(normalized_counts) %in% red_de_all,]

pdf(file="grap_pheatmap.pdf",width=12,height=10,onefile=FALSE)
pheatmap(red_de_expr, cluster_rows=TRUE, show_rownames=FALSE,scale="row",fontsize_row=7,fontsize_col=12,cluster_cols=FALSE)
dev.off()

write.csv(red_de_expr,file="pheatmap_data.csv")
