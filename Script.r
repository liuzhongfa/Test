
library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(limma)
###1.TCGA data deal####
rm(list = ls())
exprSet <- read.table('1.TCGA data deal/mRNA.symbol.txt',
                      sep = '\t',check.names = F,header = T)

exprSet0 <- exprSet[,c(1,which(substr(colnames(exprSet),14,16)=='11A'),
                       which(substr(colnames(exprSet),14,16)=='01A') ,
                       which(substr(colnames(exprSet),14,16)=='02A'))]

exprSet1 <- aggregate(.~id,exprSet0,max)
rownames(exprSet1) <- exprSet1$id
exprSet1 <- exprSet1[,-1]
exprSet2 <- exprSet1[rowMeans(exprSet1)>0,]
write.table(cbind(id = row.names(exprSet2), exprSet2),file="1.TCGA data deal/LIHC_mRNA.symbolExprUniMax.txt",sep="\t",quote = F,row.names=F)
save(exprSet2, file = '1.TCGA data deal/LIHC_mRNA.symbolExprUniMax.Rdata')


load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMax.Rdata')
Expr <- exprSet2
exprSet1 <- Expr

tcgaReplicateFilter = function(tsb, analyte_target=c("DNA","RNA"), 
                               decreasing=TRUE, analyte_position=20, 
                               plate=c(22,25), portion=c(18,19), 
                               filter_FFPE=FALSE, full_barcode=FALSE){
  # basically, user provide tsb and analyte_target is fine. If you
  # want to filter FFPE samples, please set filter_FFPE and full_barcode
  # all to TRUE, and tsb must have nchar of 28
  
  analyte_target = match.arg(analyte_target)
  # Strings in R are largely lexicographic
  # see ??base::Comparison
  
  # filter FFPE samples
  # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html> 
  if(full_barcode & filter_FFPE){
    ffpe = c("TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01", 
             "TCGA-44-2656-01B-06D-A276-05", "TCGA-44-2656-01B-06D-A27C-26", 
             "TCGA-44-2656-01B-06R-A277-07", "TCGA-44-2662-01B-02D-A271-08", 
             "TCGA-44-2662-01B-02D-A273-01", "TCGA-44-2662-01B-02R-A277-07", 
             "TCGA-44-2665-01B-06D-A271-08", "TCGA-44-2665-01B-06D-A273-01", 
             "TCGA-44-2665-01B-06D-A276-05", "TCGA-44-2665-01B-06R-A277-07", 
             "TCGA-44-2666-01B-02D-A271-08", "TCGA-44-2666-01B-02D-A273-01", 
             "TCGA-44-2666-01B-02D-A276-05", "TCGA-44-2666-01B-02D-A27C-26", 
             "TCGA-44-2666-01B-02R-A277-07", "TCGA-44-2668-01B-02D-A271-08", 
             "TCGA-44-2668-01B-02D-A273-01", "TCGA-44-2668-01B-02D-A276-05", 
             "TCGA-44-2668-01B-02D-A27C-26", "TCGA-44-2668-01B-02R-A277-07", 
             "TCGA-44-3917-01B-02D-A271-08", "TCGA-44-3917-01B-02D-A273-01", 
             "TCGA-44-3917-01B-02D-A276-05", "TCGA-44-3917-01B-02D-A27C-26", 
             "TCGA-44-3917-01B-02R-A277-07", "TCGA-44-3918-01B-02D-A271-08", 
             "TCGA-44-3918-01B-02D-A273-01", "TCGA-44-3918-01B-02D-A276-05", 
             "TCGA-44-3918-01B-02D-A27C-26", "TCGA-44-3918-01B-02R-A277-07", 
             "TCGA-44-4112-01B-06D-A271-08", "TCGA-44-4112-01B-06D-A273-01", 
             "TCGA-44-4112-01B-06D-A276-05", "TCGA-44-4112-01B-06D-A27C-26", 
             "TCGA-44-4112-01B-06R-A277-07", "TCGA-44-5645-01B-04D-A271-08", 
             "TCGA-44-5645-01B-04D-A273-01", "TCGA-44-5645-01B-04D-A276-05", 
             "TCGA-44-5645-01B-04D-A27C-26", "TCGA-44-5645-01B-04R-A277-07", 
             "TCGA-44-6146-01B-04D-A271-08", "TCGA-44-6146-01B-04D-A273-01", 
             "TCGA-44-6146-01B-04D-A276-05", "TCGA-44-6146-01B-04D-A27C-26", 
             "TCGA-44-6146-01B-04R-A277-07", "TCGA-44-6146-01B-04R-A27D-13", 
             "TCGA-44-6147-01B-06D-A271-08", "TCGA-44-6147-01B-06D-A273-01", 
             "TCGA-44-6147-01B-06D-A276-05", "TCGA-44-6147-01B-06D-A27C-26", 
             "TCGA-44-6147-01B-06R-A277-07", "TCGA-44-6147-01B-06R-A27D-13", 
             "TCGA-44-6775-01C-02D-A271-08", "TCGA-44-6775-01C-02D-A273-01", 
             "TCGA-44-6775-01C-02D-A276-05", "TCGA-44-6775-01C-02D-A27C-26", 
             "TCGA-44-6775-01C-02R-A277-07", "TCGA-44-6775-01C-02R-A27D-13", 
             "TCGA-A6-2674-01B-04D-A270-10", "TCGA-A6-2674-01B-04R-A277-07", 
             "TCGA-A6-2677-01B-02D-A270-10", "TCGA-A6-2677-01B-02D-A274-01", 
             "TCGA-A6-2677-01B-02D-A27A-05", "TCGA-A6-2677-01B-02D-A27E-26", 
             "TCGA-A6-2677-01B-02R-A277-07", "TCGA-A6-2684-01C-08D-A270-10", 
             "TCGA-A6-2684-01C-08D-A274-01", "TCGA-A6-2684-01C-08D-A27A-05", 
             "TCGA-A6-2684-01C-08D-A27E-26", "TCGA-A6-2684-01C-08R-A277-07", 
             "TCGA-A6-3809-01B-04D-A270-10", "TCGA-A6-3809-01B-04D-A274-01", 
             "TCGA-A6-3809-01B-04D-A27A-05", "TCGA-A6-3809-01B-04D-A27E-26", 
             "TCGA-A6-3809-01B-04R-A277-07", "TCGA-A6-3810-01B-04D-A270-10", 
             "TCGA-A6-3810-01B-04D-A274-01", "TCGA-A6-3810-01B-04D-A27A-05", 
             "TCGA-A6-3810-01B-04D-A27E-26", "TCGA-A6-3810-01B-04R-A277-07", 
             "TCGA-A6-5656-01B-02D-A270-10", "TCGA-A6-5656-01B-02D-A274-01", 
             "TCGA-A6-5656-01B-02D-A27A-05", "TCGA-A6-5656-01B-02D-A27E-26", 
             "TCGA-A6-5656-01B-02R-A277-07", "TCGA-A6-5656-01B-02R-A27D-13", 
             "TCGA-A6-5659-01B-04D-A270-10", "TCGA-A6-5659-01B-04D-A274-01", 
             "TCGA-A6-5659-01B-04D-A27A-05", "TCGA-A6-5659-01B-04D-A27E-26", 
             "TCGA-A6-5659-01B-04R-A277-07", "TCGA-A6-6650-01B-02D-A270-10", 
             "TCGA-A6-6650-01B-02D-A274-01", "TCGA-A6-6650-01B-02D-A27A-05", 
             "TCGA-A6-6650-01B-02D-A27E-26", "TCGA-A6-6650-01B-02R-A277-07", 
             "TCGA-A6-6650-01B-02R-A27D-13", "TCGA-A6-6780-01B-04D-A270-10", 
             "TCGA-A6-6780-01B-04D-A274-01", "TCGA-A6-6780-01B-04D-A27A-05", 
             "TCGA-A6-6780-01B-04D-A27E-26", "TCGA-A6-6780-01B-04R-A277-07", 
             "TCGA-A6-6780-01B-04R-A27D-13", "TCGA-A6-6781-01B-06D-A270-10", 
             "TCGA-A6-6781-01B-06D-A274-01", "TCGA-A6-6781-01B-06D-A27A-05", 
             "TCGA-A6-6781-01B-06R-A277-07", "TCGA-A6-6781-01B-06R-A27D-13", 
             "TCGA-A7-A0DB-01C-02D-A272-09", "TCGA-A7-A0DB-01C-02R-A277-07", 
             "TCGA-A7-A0DB-01C-02R-A27D-13", "TCGA-A7-A13D-01B-04D-A272-09", 
             "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A13D-01B-04R-A27D-13", 
             "TCGA-A7-A13E-01B-06D-A272-09", "TCGA-A7-A13E-01B-06R-A277-07", 
             "TCGA-A7-A13E-01B-06R-A27D-13", "TCGA-A7-A26E-01B-06D-A272-09", 
             "TCGA-A7-A26E-01B-06D-A275-01", "TCGA-A7-A26E-01B-06D-A27B-05", 
             "TCGA-A7-A26E-01B-06R-A277-07", "TCGA-A7-A26E-01B-06R-A27D-13", 
             "TCGA-A7-A26J-01B-02D-A272-09", "TCGA-A7-A26J-01B-02D-A275-01", 
             "TCGA-A7-A26J-01B-02D-A27B-05", "TCGA-A7-A26J-01B-02D-A27F-26", 
             "TCGA-A7-A26J-01B-02R-A277-07", "TCGA-A7-A26J-01B-02R-A27D-13", 
             "TCGA-B2-3923-01B-10D-A270-10", "TCGA-B2-3923-01B-10R-A277-07", 
             "TCGA-B2-3923-01B-10R-A27D-13", "TCGA-B2-3924-01B-03D-A270-10", 
             "TCGA-B2-3924-01B-03D-A274-01", "TCGA-B2-3924-01B-03D-A27A-05", 
             "TCGA-B2-3924-01B-03D-A27E-26", "TCGA-B2-3924-01B-03R-A277-07", 
             "TCGA-B2-3924-01B-03R-A27D-13", "TCGA-B2-5633-01B-04D-A270-10", 
             "TCGA-B2-5633-01B-04D-A274-01", "TCGA-B2-5633-01B-04D-A27A-05", 
             "TCGA-B2-5633-01B-04D-A27E-26", "TCGA-B2-5633-01B-04R-A277-07", 
             "TCGA-B2-5633-01B-04R-A27D-13", "TCGA-B2-5635-01B-04D-A270-10", 
             "TCGA-B2-5635-01B-04D-A274-01", "TCGA-B2-5635-01B-04D-A27A-05", 
             "TCGA-B2-5635-01B-04D-A27E-26", "TCGA-B2-5635-01B-04R-A277-07", 
             "TCGA-B2-5635-01B-04R-A27D-13", "TCGA-BK-A0CA-01B-02D-A272-09", 
             "TCGA-BK-A0CA-01B-02D-A275-01", "TCGA-BK-A0CA-01B-02D-A27B-05", 
             "TCGA-BK-A0CA-01B-02D-A27F-26", "TCGA-BK-A0CA-01B-02R-A277-07", 
             "TCGA-BK-A0CA-01B-02R-A27D-13", "TCGA-BK-A0CC-01B-04D-A272-09", 
             "TCGA-BK-A0CC-01B-04D-A275-01", "TCGA-BK-A0CC-01B-04D-A27B-05", 
             "TCGA-BK-A0CC-01B-04R-A277-07", "TCGA-BK-A0CC-01B-04R-A27D-13", 
             "TCGA-BK-A139-01C-08D-A272-09", "TCGA-BK-A139-01C-08D-A275-01", 
             "TCGA-BK-A139-01C-08D-A27B-05", "TCGA-BK-A139-01C-08D-A27F-26", 
             "TCGA-BK-A139-01C-08R-A277-07", "TCGA-BK-A139-01C-08R-A27D-13", 
             "TCGA-BK-A26L-01C-04D-A272-09", "TCGA-BK-A26L-01C-04D-A275-01", 
             "TCGA-BK-A26L-01C-04D-A27B-05", "TCGA-BK-A26L-01C-04D-A27F-26", 
             "TCGA-BK-A26L-01C-04R-A277-07", "TCGA-BK-A26L-01C-04R-A27D-13", 
             "TCGA-BL-A0C8-01B-04D-A271-08", "TCGA-BL-A0C8-01B-04D-A273-01", 
             "TCGA-BL-A0C8-01B-04D-A276-05", "TCGA-BL-A0C8-01B-04D-A27C-26", 
             "TCGA-BL-A0C8-01B-04R-A277-07", "TCGA-BL-A0C8-01B-04R-A27D-13", 
             "TCGA-BL-A13I-01B-04D-A271-08", "TCGA-BL-A13I-01B-04D-A276-05", 
             "TCGA-BL-A13I-01B-04R-A277-07", "TCGA-BL-A13I-01B-04R-A27D-13", 
             "TCGA-BL-A13J-01B-04D-A271-08", "TCGA-BL-A13J-01B-04D-A273-01", 
             "TCGA-BL-A13J-01B-04D-A276-05", "TCGA-BL-A13J-01B-04D-A27C-26", 
             "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13")
    
    tsb = setdiff(tsb, tsb[which(tsb %in% ffpe)])
  }
  
  # find repeated samples
  sampleID = substr(tsb, start = 1, stop = 15)
  dp_samples = unique(sampleID[duplicated(sampleID)])
  
  if(length(dp_samples)==0){
    message("ooo Not find any duplicated barcodes, return original input..")
    tsb
  }else{
    uniq_tsb = tsb[! sampleID %in% dp_samples]
    dp_tsb = setdiff(tsb, uniq_tsb)
    
    add_tsb = c()
    
    # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
    # if analyte_target = "DNA"
    # analyte:  D > G,W,X
    if(analyte_target == "DNA"){
      for(x in dp_samples){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        if(any(analytes == "D") & !(all(analytes == "D"))){
          aliquot = mulaliquots[which(analytes == "D")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }
        
      }
    }else{
      # if analyte_target = "RNA"
      # analyte: H > R > T 
      for(x in dp_samples){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        if(any(analytes == "H") & !(all(analytes == "H"))){
          aliquot = mulaliquots[which(analytes == "H")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
          
        }else if(any(analytes == "R") & !(all(analytes == "R"))){
          aliquot = mulaliquots[which(analytes == "R")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }else if(any(analytes == "T") & !(all(analytes == "T"))){
          aliquot = mulaliquots[which(analytes == "T")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }
        
      }
    }
    
    
    if(length(dp_tsb) == 0){
      message("ooo Filter barcodes successfully!")
      c(uniq_tsb, add_tsb)
    }else{
      # filter according to portion number
      sampleID_res = substr(dp_tsb, start=1, stop=15)
      dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
      
      for(x in dp_samples_res){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        portion_codes = substr(mulaliquots,
                               start = portion[1],
                               stop = portion[2])
        portion_keep = sort(portion_codes, decreasing = decreasing)[1]
        if(!all(portion_codes == portion_keep)){
          if(length(which(portion_codes == portion_keep)) == 1){
            add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }else{
            dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)])
          }
          
        }
      }
      
      if(length(dp_tsb)==0){
        message("ooo Filter barcodes successfully!")
        c(uniq_tsb, add_tsb)
      }else{
        # filter according to plate number
        sampleID_res = substr(dp_tsb, start=1, stop=15)
        dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
        for(x in dp_samples_res){
          mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
          plate_codes = substr(mulaliquots,
                               start = plate[1],
                               stop = plate[2])
          plate_keep = sort(plate_codes, decreasing = decreasing)[1]
          add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)])
          dp_tsb = setdiff(dp_tsb, mulaliquots)
        }
        
        if(length(dp_tsb)==0){
          message("ooo Filter barcodes successfully!")
          c(uniq_tsb, add_tsb)
        }else{
          message("ooo barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
          c(uniq_tsb, add_tsb)
        }
      }
    }
  }
}
dup.filt <- tcgaReplicateFilter(colnames(exprSet1),'RNA')
exprSet2 <- exprSet1[,dup.filt]

tumor <- colnames(exprSet2)[which(substr(colnames(exprSet2),14,16)=='01A' |
                                    substr(colnames(exprSet2),14,16)=='02A')]
tumordata <- exprSet2[,tumor]
range(tumordata)
tumordata2 <- log2(tumordata+1)
range(tumordata2)
write.table(cbind(id = row.names(tumordata2), tumordata2),file="1.TCGA data deal/LIHC_mRNA.symbolExprUniMax_onlyTumor.txt",
            sep="\t",quote=F,col.names=T, row.names = F)
save(tumordata2, file = '1.TCGA data deal/LIHC_mRNA.symbolExprUniMax_onlyTumor.RData')

tumordata3 <- as.data.frame(t(tumordata2))
sampleID <- substr(row.names(tumordata3),1,12)
tumordata3$id <- sampleID
tumordata3$fullid <- row.names(tumordata3)

preclin <- read.table('1.TCGA data deal/clinicaldelete0.txt',
                      header=T,sep="\t",check.names=F)
merge <- merge(preclin, tumordata3, by = 'id')
merge <- na.omit(merge)
merge <- merge %>% dplyr::select(fullid,everything()) 

merge2 <- merge[order(merge$fullid),]
rownames(merge2) <- merge2$fullid
tmp <- merge2[,c('id','fullid')]
tmp1 <- aggregate(.~id,tmp,min)
merge2 <- merge2[tmp1$fullid,]
merge2 <- na.omit(merge2)
write.table(merge2, file = '1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.txt', sep="\t", row.names=F, col.names = T, quote=F)
save(merge2, file = '1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')



###2.Diffexpr####
library(limma)
library(pheatmap)
rm(list = ls())

load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMax.Rdata')
pvalue <-0.05
logFoldChange <- 1

design <- read.table("2.Diffexpr/LIHC_sampleTypelimma.txt",header=T,sep="\t",check.names=F,row.names = 1)
dat <- exprSet2[,row.names(design)]
range(dat)
dat <- log2(dat+1)
dat <- na.omit(dat)
range(dat)

###boxplot
pdf(file="2.Diffexpr/LIHC_mRNA.symbolExprUniMax_boxplot.pdf",width=50,height=5)
par(cex = 0.7)
n.sample = ncol(dat)
if(n.sample>40)par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(dat, col = cols, main="expression value")
dev.off()

design
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.table(cbind(Symbol=rownames(allDiff),allDiff),file="2.Diffexpr/TNlimmaOut.txt",sep="\t",row.names = F,quote = F)

diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffSig),diffSig), file="2.Diffexpr/TNlimmadiffSig.txt",sep="\t",row.names = F,quote = F)

diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
write.table(cbind(Symbol=rownames(diffUp),diffUp), file="2.Diffexpr/TNlimmadiffUP.txt",sep="\t",row.names = F,quote = F)

diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffDown),diffDown), file="2.Diffexpr/TNlimmadiffDown.txt",sep="\t",row.names = F,quote = F)

##volcano plot
pvalue <-0.05
logFC <- 1
allDiff <- read.table('2.Diffexpr/TNlimmaOut.txt',header = T,check.names = F,row.names = 1,sep = '\t')
allDiff$Significant <- ifelse(allDiff$P.Value<pvalue & abs(allDiff$logFC)>= logFC,
                              ifelse(allDiff$logFC> logFC,'up','down'),'no')
mycol <- c("#3CB371","#3D3D3D","#FF4500")
library(ggplot2)
pdf(file="2.Diffexpr/TNdiffexp.volcano.pdf",width=6,height=5.5)
p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
  geom_point(size=1.2)+theme_bw()+
  scale_color_manual(values = mycol,name='Significant')+
  labs(title="DEGs",x="Log2FC",y="-log10 (P.value)")+
  geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
  geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
print(p)
dev.off()

###heatmap
group_list=c(rep('Normal',50),rep('Tumor',371))
diff <- rownames(diffSig)
diffexp <- dat[diff,]
write.table(cbind(id=row.names(diffexp),diffexp), file="2.Diffexpr/TNlimmadiffSigExpr.txt",sep="\t",row.names = F, col.names = T, quote = F)

difftop100 <- head(diff,100)
diffexp <- dat[difftop100,]

annotation_col <- data.frame(Type = factor(group_list,levels = c("Normal","Tumor")))
rownames(annotation_col) <- colnames(diffexp)
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = '2.Diffexpr/TNdiffexp.heatmap.pdf',width = 22,height = 12)  
pheatmap(diffexp,cellwidth = 3,cellheight = 7,
         method="pearson",
         scale="row", 
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(color.key)(50),
         show_colnames=F,show_rownames =F,
         annotation_col = annotation_col,
         #treeheight_row = "0",treeheight_col = "0",
         border_color = "NA")
dev.off()


###Diff AREGene heatmap#####
##venn plot use TBtools
diffAREG <- read.table('3.DiffAREG/diffAREGene.diffsigGene.and.AREGene.common.189.txt', check.names = F, header = F)
diffAREGexpr <- diffexp[diffAREG$V1,]
write.table(cbind(id=row.names(diffAREGexpr),diffAREGexpr), file="diffAREGeneexpr.txt",sep="\t",row.names = F, col.names = T, quote = F)

group_list=c(rep('Normal',50),rep('Tumor',371))
annotation_col <- data.frame(Type = factor(group_list,levels = c("Normal","Tumor")))
rownames(annotation_col) <- colnames(diffAREGexpr)
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = 'diffAREGene.heatmap.pdf',width = 22,height = 12)  
pheatmap(diffAREGexpr,cellwidth = 3,cellheight = 4,
         method="pearson", 
         scale="row", 
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(color.key)(50),
         show_colnames=F,show_rownames =F,
         annotation_col = annotation_col,
         #treeheight_row = "0",treeheight_col = "0",
         border_color = "NA")
dev.off()


###Diff AREGene enrich#####
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(GOplot))
suppressMessages(library(GOSemSim))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ReactomePA))

diffAREG <- read.table('3.DiffAREG/diffAREGene.diffsigGene.and.AREGene.common.189.txt')
diffAREG.id <- diffAREG$V1
gene <- mapIds(x = org.Hs.eg.db, keys = diffAREG.id, column = 'ENTREZID', keytype = 'SYMBOL', multiVals='filter')
gene2 <- as.data.frame(gene)
gene2 <- cbind(row.names(gene2), gene2)
colnames(gene2) <- c('SYMBOLID', 'ENTREZID')
write.table(gene2, file = '3.DiffAREG/geneSymbol2EntreID.txt', sep = '\t', quote = F, row.names = F)

GO <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               ont = "all",
               readable = T)

KEGG <- enrichKEGG(gene = gene, 
                   organism = "hsa", 
                   keyType = "kegg", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2)

write.table(GO, file = "3.DiffAREG/diffAREGeneEnrichGO.txt", quote = F, row.names = F, sep = '\t')               
write.table(KEGG, file = "3.DiffAREG/diffAREGeneEnrichKEGG.txt", quote = F, row.names = F, sep = '\t')

pdf(file="3.DiffAREG/diffAREGeneEnrichGO.bubble.pdf",width = 9,height = 10)
dotplot(GO, showCategory = 10, split="ONTOLOGY") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()

pdf(file="3.DiffAREG/diffAREGeneEnrichKEGG.bar.pdf",width = 10,height = 6)
barplot(KEGG,showCategory = 15,font.size = 15)
dev.off()

## Diff AREGene Down Up enrich  
diffAREG <- read.table('3.DiffAREG/diffAREG_down.txt')
diffAREG.id <- diffAREG$V1
gene <- mapIds(x = org.Hs.eg.db, keys = diffAREG.id, column = 'ENTREZID', keytype = 'SYMBOL', multiVals='filter')

GO <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               ont = "all",
               readable = T)

KEGG <- enrichKEGG(gene = gene, 
                   organism = "hsa", 
                   keyType = "kegg", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2)

write.table(GO, file = "3.DiffAREG/diffAREGeneEnrichGO_down.txt", quote = F, row.names = F, sep = '\t')               
write.table(KEGG, file = "3.DiffAREG/diffAREGeneEnrichKEGG_down.txt", quote = F, row.names = F, sep = '\t')

pdf(file="3.DiffAREG/diffAREGeneEnrichGO_down.bubble.pdf",width = 9,height = 10)
dotplot(GO, showCategory = 10, split="ONTOLOGY") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()

pdf(file="3.DiffAREG/diffAREGeneEnrichKEGG_down.bar.pdf",width = 10,height = 6)
barplot(KEGG,showCategory = 15,font.size = 15)
dev.off()


###Prepare cox input####
rm(list = ls())
options(stringsAsFactors = F)

##TCGA data
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')

diffAREG <- read.table('3.DiffAREG//diffAREGene.diffsigGene.and.AREGene.common.189.txt')
diffAREGexpr <- cbind(merge2[, c(1:11)], merge2[, diffAREG$V1])
write.table(diffAREGexpr, file = '4.Prepare cox input/diffAREGeneexprClin_onlyTumor.txt', sep="\t", row.names=F, col.names = T, quote=F)

##DEO data
probeExpr <- read.table('4.Prepare cox input/GSE14520_probeExpr.txt', header = T, sep = '\t', check.names = F)
GPL <- read.table('4.Prepare cox input/GPL3921-25447probe-symble.txt', header = T, sep = '\t', check.names = F)
symbolExpr <- merge(GPL, probeExpr, by = 'ID_REF')
symbolExpr <- na.omit(symbolExpr)
symbolExpr <- symbolExpr[, -1]
range(symbolExpr[,-1])
write.table(symbolExpr, file = '4.Prepare cox input/GSE14520_symbolExpr.txt', sep="\t", row.names = F, quote = F)

gene<-read.table("4.Prepare cox input/GSE14520_symbolExpr.txt",sep="\t",header=T)
data3 <- aggregate( . ~ id, data=gene, max)  
write.table(data3,file="4.Prepare cox input/GSE14520_symbolExprUniMax.txt",sep="\t",quote = F,row.names=F)

Expr <- read.table('4.Prepare cox input/GSE14520_symbolExprUniMax.txt', header = T, sep = '\t', check.names = F, row.names = 1)
clin <- read.table('4.Prepare cox input/GSE14520_clinical.txt', header = T , sep = '\t', check.names = F)
range(Expr)
Expr <- as.data.frame(t(Expr))
Expr <- cbind(id = row.names(Expr), Expr)
Exprclin <- merge(clin, Expr, by = 'id')
Exprclin  <- na.omit(Exprclin)
write.table(Exprclin, file = '4.Prepare cox input/GSE14520_Exprclin.txt', sep="\t", row.names = F, quote = F)
row.names(Exprclin) <- Exprclin$id
Exprclin <- Exprclin[, -1]
save(Exprclin, file = '4.Prepare cox input/GSE14520_Exprclin.RData')


###COX model#####

library(survivalROC)
library(survival)
library(survminer)
library(caret)
library(tibble)
library(dplyr)
library(randomForestSRC)
library(forestplot)
library(pheatmap)
library(glmnet)

rm(list = ls())

rtt0 <- read.table('4.Prepare cox input/diffAREGeneexprClin_onlyTumor.txt',header = T,sep = '\t',check.names = F, row.names = 1)
rtt <- rtt0[, -c(1,4:10)]
range(rtt[,3:ncol(rtt)])

GSE14520 <- read.table('4.Prepare cox input/GSE14520_Exprclin.txt',header = T,sep = '\t',check.names = F, row.names = 1)
range(GSE14520[,3:ncol(GSE14520)])     

overlap <- intersect(colnames(rtt), colnames(GSE14520))

rtt <- rtt[, overlap]
GSE14520 <- GSE14520[, overlap]

set.seed(18)

coxinput <- rtt
outTab=data.frame()

######uni cox####
print('start unicox')
for(i in colnames(coxinput[,3:ncol(coxinput)])){
  cox <- coxph(Surv(futime, fustat) ~ coxinput[,i], data = coxinput)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="5.CoxModule/uniCox.txt",sep="\t",row.names=F,quote=F)

outTab[,c(2:ncol(outTab))] <- apply(outTab[,c(2:ncol(outTab))],2,function(x){as.numeric(x)})
outTab <- outTab[order(outTab$pvalue),]
uniCoxSig =  outTab[(outTab$pvalue < 0.05),]
write.table(uniCoxSig,file="5.CoxModule/uniCoxSig0.05.txt", sep="\t", row.names=F, quote=F)

## uni cox forest plot
uniFor <- read.table("5.CoxModule/uniCoxSig0.05.txt", header=T, sep="\t",row.names=1,check.names=F)
gene <- rownames(uniFor)
hr <- sprintf("%.3f",uniFor$"HR")
hrLow  <- sprintf("%.3f",uniFor$"HR.95L")
hrHigh <- sprintf("%.3f",uniFor$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(uniFor$pvalue<0.001, "<0.001", sprintf("%.3f", uniFor$pvalue))

pdf(file="5.CoxModule/uniCoxSigForest.pdf", width = 8,height = 16)
n <- nrow(uniFor)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#####lasso####
unicox <- outTab[outTab$pvalue <= 0.05,]$id

lassoinput <- coxinput[,c('futime','fustat',unicox)]

x=as.matrix(lassoinput[,c(3:ncol(lassoinput))])
y=data.matrix(Surv(lassoinput$futime,lassoinput$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 5000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 5000)

pdf(file = "5.CoxModule/lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

pdf(file = "5.CoxModule/cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
text(log(cvfit$lambda.min),12,cex=0.8,
     labels = paste0('lambda.min = \n',round(cvfit$lambda.min,3)))
text(log(cvfit$lambda.1se),12.2,cex=0.8,
     labels = paste0('lambda.lse = \n',round(cvfit$lambda.1se,3)))
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(as.matrix(coef)!= 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]

#######multi cox#####
multi <- coxinput[,c('futime','fustat',lassoGene)]
multi[,"futime"]=multi[,"futime"]/365

multiCox=coxph(Surv(futime, fustat) ~ ., data = multi)
multiCox=step(multiCox,direction = "both")#根据AIC值过滤，向前向后方法进一步优化
multiCoxSum=summary(multiCox)
outTab1=data.frame()
outTab1=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

outTab1=cbind(id=gsub("`","",row.names(outTab1)),outTab1)
outTab1 <- data.frame(outTab1)
outTab1[,c(2:ncol(outTab1))] <- apply(outTab1[,c(2:ncol(outTab1))],2,function(x){as.numeric(x)})
modelgene <- outTab1$id
modelgene=gsub("`","",row.names(outTab1))
outTab1 <- outTab1[order(outTab1$pvalue),]
write.table(outTab1,file="5.CoxModule/multiCox.txt",sep="\t",row.names=F,quote=F)

###multi cox forest plot
multiFor <- read.table("5.CoxModule/multiCox.txt", header=T, sep="\t",row.names=1,check.names=F)
multiFor <- multiFor[, -1]
gene <- rownames(multiFor)
gene=gsub("`","",gene)
hr <- sprintf("%.3f",multiFor$"HR")
hrLow  <- sprintf("%.3f",multiFor$"HR.95L")
hrHigh <- sprintf("%.3f",multiFor$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(multiFor$pvalue<0.001, "<0.001", sprintf("%.3f", multiFor$pvalue))

pdf(file="5.CoxModule/multiCoxForest.pdf", width = 8,height = 5)
n <- nrow(multiFor)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

###Train set riskscore####
riskScore=predict(multiCox,type="risk")
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
write.table(cbind(id=rownames(multi), multi[,outCol],riskScore),
            file="5.CoxModule/trainRisk.txt",
            sep="\t",
            quote=F,
            row.names=F)
rt <- read.table('5.CoxModule/trainRisk.txt',header = T,sep = '\t',check.names = F)

rt$risk <- ifelse(rt$riskScore>=median(rt$riskScore),'High','Low')

rt <- rt[order(rt$risk,decreasing = T),]
write.table(rt, file="5.CoxModule/trainRisk.txt",sep="\t",quote=F,row.names=F)

###Validation set riskscore####
print('start GSE14520 validation')

GSE14520gene <- c("futime","fustat",modelgene)
GSE14520.exp <- GSE14520[,GSE14520gene]
GSE14520.exp <- na.omit(GSE14520.exp)
GSE14520.exp[,"futime"]=GSE14520.exp[,"futime"]/365
write.table(cbind(id = rownames(GSE14520.exp), GSE14520.exp),'5.CoxModule/GSE14520.coxinput.txt',sep = '\t',quote = F,row.names = F)

riskScore=predict(multiCox,type="risk",newdata = GSE14520.exp)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
write.table(cbind(id = rownames(GSE14520.exp), GSE14520.exp[,outCol],riskScore),
            file="5.CoxModule/GSE14520.risk.txt",
            sep="\t",
            quote=F,
            row.names=F)
rt <- read.table('5.CoxModule/GSE14520.risk.txt',header = T,sep = '\t',check.names = F)

rt$risk <- ifelse(rt$riskScore>=median(rt$riskScore),'High','Low')
rt <- rt[order(rt$riskScore,decreasing = T),]
write.table(rt, file="5.CoxModule/GSE14520.risk.txt",sep="\t",quote=F,row.names=F)

###Plot######
library(survivalROC)
library(survival)
library(survminer)
library(caret)
library(tibble)
library(dplyr)
library(pheatmap)
#####Train set####
######Risk curve、survival######
datalast <- read.table('5.CoxModule/trainRisk.txt',header = T,check.names = T,sep = '\t')
risk=datalast[order(datalast$riskScore),]
median <- median(datalast$riskScore)
riskClass <- datalast[, "risk"]
lowLength <- length(riskClass[riskClass == "Low"])
highLength <- length(riskClass[riskClass == "High"])

pdf(file = "5.CoxModule/train.riskscore+state.pdf", width = 16, height = 13)
par(mfrow = c(2, 1), cex.axis = 2, mar = c(2, 5, 2, 2))
col <- c()
col[sort(datalast$riskScore) <= median] <- "blue"
col[sort(datalast$riskScore) > median] <- "red"
plot(sort(datalast$riskScore), col = col, type = "p", pch = 20, cex = 2, 
     mgp = c(2.5, 1, 0), cex.lab = 1.8, cex.main = 2,
     xlab = 'Patients (increasing risk socre)',
     ylab = "Risk Score", 
     main = 'train set')
legend('topleft', c('Low_risk','High_risk'), cex = 1.5, col=c('blue','red'), pch = 16:16, bty='n')
box(lwd = 2)
abline(h = median, v = lowLength, lty = 2)

datalastSORT <- datalast[order(datalast[, "riskScore"]), ]
col <- c()
col[datalastSORT[, "fustat"] == 0] <- "blue"
col[datalastSORT[, "fustat"] == 1] <- "red"
par(mar = c(2, 5, 1, 2))
plot(datalastSORT[, "futime"], col = col, pch = 20, cex = 2, axes = T, 
     mgp = c(2.5, 1, 0), cex.lab = 1.8, 
     xlab = "Patients (increasing risk socre)",
     ylab = "Survival time (years)")
legend("topleft", c("Alive", "Dead"), cex = 1.5, col = c("blue", "red"),  pch = 16:16, bty='n')
box(lwd = 2)
abline(v = lowLength, lty = 2)
dev.off()

######Heatmap######
data <- read.table('5.CoxModule/trainRisk.txt', header=T,sep="\t",check.names=F,row.names = 1)
data <- data[order(data$risk, decreasing = T),]
data_p <- data[,-c(1:2, ncol(data)-1, ncol(data))]
data_p <- t(data_p)
annotation_col = data.frame(Risk = factor(data$risk,levels = c('Low','High')))
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
rownames(annotation_col) <- colnames(data_p)

pdf('5.CoxModule/train.riskheatmap.pdf',width =10,height = 6)
pheatmap(data_p, cellwidth = 1.4, cellheight = 50,
         method="pearson", 
         scale="row", 
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(color.key)(100),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         border_color = "NA",main = 'train Set')
dev.off()

####### KM survival curve######
survival <- read.table('5.CoxModule/trainRisk.txt',header = T,sep = '\t',check.names = F)

fit <- survfit(Surv(futime, fustat) ~risk,data = survival)
pdf(file = "5.CoxModule/train.survival.pdf",width = 7,height = 7)
ggsurvplot(fit, data = survival,
           conf.int = TRUE, 
           conf.int.style = "ribbon",
           risk.table = T, 
           pval=TRUE,
           tables.height = 0.25,
           censor = T, 
           palette = c("red", "blue"), 
           xlab = "Time (Years)",
           legend.title = "Risk",
           title="Kaplan-Meier Curve for Survival",
           font.legend = 15)
dev.off()

######## ROC curve#####
riskdata <- read.table('5.CoxModule/trainRisk.txt',header = T,sep = '\t',check.names = F)

pdf(file='5.CoxModule/trainROC 1 3 5年.pdf', width = 7,height = 7)
aucText=c()
rt <- riskdata
rocCol = rainbow(4)
roc <- survivalROC(Stime = riskdata$futime,
                   status = riskdata$fustat,
                   marker = riskdata$riskScore,
                   predict.time = 1,method='KM')

par(mar=c(4,4,1,1),mgp=c(2,0.5,0))
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 1.5, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.5)
aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1,lty=2)
j=1

for(i in c(3,5)){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1*i, method = 'KM')
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 1.5)
}

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


####Validation set#####
######Risk curve、survival######
datalast <- read.table('5.CoxModule/GSE14520.risk.txt',header = T,check.names = T,sep = '\t')
risk=datalast[order(datalast$riskScore),]
median <- median(datalast$riskScore)
riskClass <- datalast[, "risk"]
lowLength <- length(riskClass[riskClass == "Low"])
highLength <- length(riskClass[riskClass == "High"])

pdf(file = "5.CoxModule/GSE14520.riskscore+state.pdf", width = 16, height = 13)
par(mfrow = c(2, 1), cex.axis = 2, mar = c(2, 5, 2, 2))
col <- c()
col[sort(datalast$riskScore) <= median] <- "blue"
col[sort(datalast$riskScore) > median] <- "red"
plot(sort(datalast$riskScore), col = col, type = "p", pch = 20, cex = 2, 
     mgp = c(2.5, 1, 0), cex.lab = 1.8, cex.main = 2,
     xlab = 'Patients (increasing risk socre)',
     ylab = "Risk Score", 
     main = 'GSE14520 set')
legend('topleft', c('Low_risk','High_risk'), cex = 1.5, col=c('blue','red'), pch = 16:16, bty='n')
box(lwd = 2)
abline(h = median, v = lowLength, lty = 2)

datalastSORT <- datalast[order(datalast[, "riskScore"]), ]
col <- c()
col[datalastSORT[, "fustat"] == 0] <- "blue"
col[datalastSORT[, "fustat"] == 1] <- "red"
par(mar = c(2, 5, 1, 2))
plot(datalastSORT[, "futime"], col = col, pch = 20, cex = 2, axes = T, 
     mgp = c(2.5, 1, 0), cex.lab = 1.8, 
     xlab = "Patients (increasing risk socre)",
     ylab = "Survival time (years)")
legend("topleft", c("Alive", "Dead"), cex = 1.5, col = c("blue", "red"),  pch = 16:16, bty='n')
box(lwd = 2)
abline(v = lowLength, lty = 2)
dev.off()

#########Heatmap######
data <- read.table('5.CoxModule/GSE14520.risk.txt', header=T,sep="\t",check.names=F,row.names = 1)
data <- data[order(data$risk, decreasing = T),]
data_p <- data[,-c(1:2, ncol(data)-1, ncol(data))]
data_p <- t(data_p)
annotation_col = data.frame(Risk = factor(data$risk,levels = c('Low','High')))  
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
rownames(annotation_col) <- colnames(data_p)

pdf('5.CoxModule/GSE14520.riskheatmap.pdf',width =10,height = 6)
pheatmap(data_p, cellwidth = 1.5, cellheight = 50,
         method="pearson", 
         scale="row", 
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(color.key)(100),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         border_color = "NA",main = 'GSE14520 Set')
dev.off()

####### KM survival curve######
survival <- read.table('5.CoxModule/GSE14520.risk.txt',header = T,sep = '\t',check.names = F)
fit <- survfit(Surv(futime, fustat) ~risk,data = survival)

pdf(file = "5.CoxModule/GSE14520.survival.pdf",width = 7,height = 7)
ggsurvplot(fit, data = survival,
           #ggtheme = theme_bw(), #想要网格就运行这行
           conf.int = TRUE, #不画置信区间，想画置信区间就把F改成T
           conf.int.style = "ribbon",#置信区间的类型，还可改为ribbon
           risk.table = T, # 绘制累计风险曲线
           pval=TRUE,
           #tables.theme = theme_void(),
           tables.height = 0.25,
           censor = T, #不显示观察值所在的位置
           palette = c("red", "blue"), #线的颜色对应高、低(自定义调色板)
           xlab = "Time (Years)",
           legend.title = "Risk",#基因名写在图例题目的位置
           title="Kaplan-Meier Curve for Survival",
           font.legend = 15)#图例的字体大小
dev.off()

########  ROC curve #####
riskdata <- read.table('5.CoxModule/GSE14520.risk.txt',header = T,sep = '\t',check.names = F)

pdf(file='5.CoxModule/GSE14520ROC 1 3 5年.pdf', width = 7,height = 7)
aucText=c()
rt <- riskdata
rocCol = rainbow(4)
roc <- survivalROC(Stime = riskdata$futime,
                   status = riskdata$fustat,
                   marker = riskdata$riskScore,
                   predict.time = 1,method='KM')

par(mar=c(4,4,1,1),mgp=c(2,0.5,0))
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 1.5, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.5)
aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1,lty=2)
j=1

for(i in c(3,5)){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1*i, method = 'KM')
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 1.5)
}

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


###Clinic and Risk#########
rm(list = ls())
library(ggpubr)
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')
preclin <-merge2
clin <- preclin[, c(1, 5:11)]
trainriskScore <- read.table('5.CoxModule/trainRisk.txt',header = T,sep = '\t',check.names = F)
 
rt <- merge(clin, trainriskScore, by.x = 'fullid', by.y = 'id')
rt <- na.omit(rt)
head(rt)
needclin <- rt[,c('fullid', 'futime', 'fustat', 'age', 'gender', 'grade', 'stage', 'T', 'N', 'M' ,'riskScore', 'risk')]
write.table(needclin,file="6.Clinic and Risk/needclin.txt", sep="\t", row.names=F, quote=F) 

needclin <- read.table("6.Clinic and Risk/needclin.txt",header=T,sep="\t",row.names=1,check.names=F)
grade <- subset(needclin, needclin$grade != 'unknow')
stage <- subset(needclin, needclin$stage != 'unknow')
ST <- subset(needclin, needclin$T != 'TX' & needclin$T != 'unknow')
SN <- subset(needclin, c(needclin$N != 'unknow' & needclin$N != 'NX'))
SM <- subset(needclin, c(needclin$M != 'unknow' & needclin$M != 'MX'))
data <- needclin
median(data$age)
data$age[which(data$age <= 61)] <- '<= 61'
data$age[which(data$age > 61)] <- '> 61'
# my_comparisons <- list(c('<= 61','> 61'))
# my_comparisons <- list(c('G1','G2'), c('G1','G3'), c('G1','G4'), c('G2','G3'), c('G2','G4'), c('G3', 'G4'))
# my_comparisons <- list(c('MALE','FEMALE'))
# my_comparisons <- list(c('Stage I-II','Stage III-IV'))
# my_comparisons <- list(c('EarlyStage','LateStage'))
# my_comparisons <- list(c('T12','T34'))
# my_comparisons <- list(c('N0','N1'))
my_comparisons <- list(c('M0','M1'))
# my_comparisons <- list(c("T1","T2"), c("T1","T3"), c("T1","T4"), c("T2","T3"), c("T2","T4"), c("T3","T4"))
# my_comparisons <- list(c("N0","N1"), c("N0","N2"), c("N0","N3"), c("N1","N2"), c("N1","N3"), c("N2","N3"))
# my_comparisons <- list(c('Stage I','Stage II'), c('Stage I','Stage III'), c('Stage I','Stage IV'), c('Stage II','Stage III'), c('Stage II','Stage IV'), c('Stage III','Stage IV'))

##不填充箱线图
pdf(file="6.Clinic and Risk/M.riskscore.pdf",width=6,height=6)
p=ggboxplot(SM, x="M", y="riskScore", merge = "flip",fill = 'M', #color = "gender", add="jitter",
            ylab="riskScore",
            xlab="",
            # order = c('<= 61','> 61'),
            # order = c('G1','G2', 'G3', 'G4'),
            # order = c('Low Grade','High Grade'),
            # order = c('EarlyStage','LateStage'),
            # order = c('T12','T34'),
            # order = c('N0','N12'),
            order = c('M0','M1'),
            # order = c("T1","T2","T3","T4"),
            # order = c("N0","N1"),
            # order = c('Stage I','Stage II','Stage III','Stage IV'),
            # palette = c("green", "blue", "red","slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2","LightSkyBlue"))  ##palette = "jco" 调整颜色较暗淡  c("green","blue","red")
            palette = c('#008B45','#FF0000','#FF00FF', '#FFD700','#9999FF', '#D1BBFF', '#EE7942','#FF88C2',  
                        '#FF6347','#33FFFF','#33FF33', '#FFA500', '#ADFF2F',
                        '#EE30A7', '#CCEEFF', '#0000AA','#770077','#9370DB'))
p=p+rotate_x_text(0)
p+stat_compare_means(aes(group=T),
                     method = 'wilcox.test',
                     label = "p.signif",
                     size = 6,
                     comparisons = my_comparisons,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  # stat_compare_means(size = 5, method = 'kruskal.test', label.y = 15, label.x = 0.8) +
  font('xlab',face = 'bold',size=14)+font('ylab',face = 'bold',size=14)+
  font('x.text',face = 'bold',size=14)+font('y.text',face = 'bold',size=14)+
  font('legend.title',face = 'bold',size=14)+font('legend.text',face = 'bold')
dev.off()


###Prognosis#######
needclin <- read.table('7.Prognosis/needclin3.txt',header = T,sep = '\t',check.names = T)
needclin <- needclin[, -12]
###uni cox
pFilter = 0.05
sigGenes = c("futime", "fustat")
outTab2=data.frame()
for(i in colnames(needclin[,4:ncol(needclin)])){
  cox <- coxph(Surv(futime, fustat) ~ needclin[,i], data = needclin)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab2=rbind(outTab2,
                cbind(id=i,
                      HR=coxSummary$conf.int[,"exp(coef)"],
                      HR.95L=coxSummary$conf.int[,"lower .95"],
                      HR.95H=coxSummary$conf.int[,"upper .95"],
                      pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  if(coxP < pFilter){
    sigGenes = c(sigGenes,i)
  }
}
write.table(outTab2,file="7.Prognosis/indep_uniCoxresult.txt", sep="\t", row.names=F, quote=F)
uniSigExp = needclin[,sigGenes]
uniSigExp = cbind(id = needclin$fullid, uniSigExp)
write.table(uniSigExp, file = "7.Prognosis/indep_uniSigExp.txt", sep="\t", row.names=F, quote=F)

###forest plot
rt <- read.table("7.Prognosis/indep_uniCoxresult.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file="7.Prognosis/indep_uniCoxForest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

####multi Cox
needclin2 <- read.table('7.Prognosis/indep_uniSigExp.txt',header = T,sep = '\t',check.names = T,row.names = 1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = needclin2)
multiCox =step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab3 = data.frame()
outTab3 = cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab3 = cbind(id=row.names(outTab3),outTab3)
write.table(outTab3,file="7.Prognosis/indep_multiCoxresult.txt",sep="\t",row.names=F,quote=F)

###forest plot
rt <- read.table("7.Prognosis/indep_multiCoxresult.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file="7.Prognosis/indep_multiCoxForest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


######nomogram #########
rm(list = ls())
library(rms)
rt <- read.table("6.Clinic and Risk/needclin.txt",header=T,sep="\t",row.names=1,check.names=F)
rt <- rt[, c('futime', 'fustat', 'T', 'riskScore')]
rt <- subset(rt, rt$T != 'TX' & rt$T != 'unknow')
rt[,"futime"]=rt[,"futime"]*365
pbc <-rt
dd <- datadist(pbc)
options(datadist="dd")
options(na.action="na.delete")
summary(pbc$futime)
coxpbc <- cph(formula = Surv(futime,fustat) ~ T+riskScore, data=pbc, x=T, y=T, surv = T, na.action=na.delete, time.inc = 365)
surv <- Survival(coxpbc) 
surv1 <- function(x) surv(365,x)
surv3 <- function(x) surv(1095,x)
surv5 <- function(x) surv(1825,x)

x <- nomogram(coxpbc,fun = list(surv1,surv3,surv5),lp=T,
              funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
              maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("7.Prognosis/nomogram_classical.pdf",width = 16,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

######calibration curve############
f1 <- cph(formula = Surv(futime,fustat) ~ T+riskScore, data=pbc, x=T, y=T, surv = T, na.action=na.delete, time.inc = 365) 
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=95,B=1000) #参数m=50表示每组50个样本进行重复计算

f3 <- cph(formula = Surv(futime,fustat) ~ T+riskScore, data=pbc, x=T, y=T, surv = T, na.action=na.delete, time.inc = 1095)  
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=95,B=1000)

f5 <- cph(formula = Surv(futime,fustat) ~ T+riskScore, data=pbc, x=T, y=T, surv = T, na.action=na.delete, time.inc = 1825) 
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=95,B=1000)

###slope
library(stringr)
caldat <- data.frame(summary(cal1))
cal1rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -1) ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -1))[["coefficients"]][["(Intercept)"]]
caldat <- data.frame(summary(cal3))
cal3rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -3)[1:6])[["coefficients"]][["(Intercept)"]]
caldat <- data.frame(summary(cal5))
cal5rate <- lm( str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -10, -3)[1:6])[["coefficients"]][["(Intercept)"]]

###c-index
set.seed(123)
v <- validate(coxpbc, dxy=TRUE, B=1000)
Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.corrected']
orig_Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.orig']
bias_corrected_c_index  <- abs(Dxy)/2+0.5  
orig_c_index <- abs(orig_Dxy)/2+0.5  
bias_corrected_c_index
orig_c_index


pdf("7.Prognosis/calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 1,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal5,lwd = 2,lty = 1,errbar.col = c("darkgreen"),
     xlim = c(0,1),ylim= c(0,1),col = c("darkgreen"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("darkgreen"), pch = 16)
abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", 
       legend = c("1-year","3-year","5-year"),
       col =c("#2166AC","#B2182B","darkgreen"), 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()

pdf("7.Prognosis/calibration_5y.pdf",width = 6,height = 6)
plot(cal5,
     lwd = 2,
     lty = 1,
     errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%) of 5 year",ylab = "Observed OS (%) of 5 year",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) 
lines(cal5[,c('mean.predicted',"KM")], 
      type = 'b', 
      lwd = 2, 
      pch = 16, 
      col = c("#FF3030"))
mtext("")
box(lwd = 1) 
abline(0,1,lty = 3,
       lwd = 2, 
       col = c("#224444")
) 
dev.off()

###GSEA input ######
##GSEA use local version GSEA_4.2.1
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')
data1 <- merge2[,-c(1:11)]
ss <- read.table('5.CoxModule/trainRisk.txt', header = T, sep = '\t', check.names = F)
data2 <- as.data.frame(t(data1))
data3 <- data2[, ss$id]
write.table(cbind(id = row.names(data3), data3), file = 'LIHC_mRNA.symbolExprUniMaxClin_onlyTumor_GSEAinput.txt', sep = '\t', quote = F, row.names = T)

###Inflammatory ####
rm(list = ls())
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')
allgeneExpr <- merge2
allgeneExpr1 <- allgeneExpr[, -c(1:11)]
range(allgeneExpr1)

gmtFile="9.Inflammatory/APM.gmt"
rt=as.data.frame(t(allgeneExpr1))
rt=as.matrix(rt)

geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

ssgseaScore=gsva(rt, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
range(ssgseaScore)
ssgseaOut <- as.data.frame(ssgseaScore)
write.table(cbind(id = row.names(ssgseaOut), ssgseaOut),file="9.Inflammatory/IFN-γssGseaOut.txt",sep="\t",quote=F,col.names=T, row.names = F)
riskScore <- read.table('5.CoxModule/trainRisk.txt', 
                        header = T, sep = '\t', check.names = F)
INF <- as.data.frame(t(ssgseaOut))
INF <- cbind(id = row.names(INF), INF)
mergeINF <- merge(INF, riskScore, by = 'id')

my_comparisons <- list(c('Low','High'))
pdf(file="9.Inflammatory/APM.risk.pdf",width=6,height=6)
p=ggboxplot(mergeINF, x="risk", y="APM", merge = "flip",fill = 'risk', #color = "gender", add="jitter",
            ylab="APM Score",
            xlab="",
            order = c('Low','High'),
            palette = c('#008B45','#FF0000','#FF00FF', '#FFD700','#9999FF', '#D1BBFF', '#EE7942','#FF88C2',  
                        '#FF6347','#33FFFF','#33FF33', '#FFA500', '#ADFF2F',
                        '#EE30A7', '#CCEEFF', '#0000AA','#770077','#9370DB'))
p=p+rotate_x_text(0)
p+stat_compare_means(aes(group=T),
                     method = 'wilcox.test',
                     label = "p.signif",
                     size = 6,
                     comparisons = my_comparisons,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  # stat_compare_means(size = 5, method = 'kruskal.test', label.y = 15, label.x = 0.8) +
  font('xlab',face = 'bold',size=14)+font('ylab',face = 'bold',size=14)+
  font('x.text',face = 'bold',size=14)+font('y.text',face = 'bold',size=14)+
  font('legend.title',face = 'bold',size=14)+font('legend.text',face = 'bold')
dev.off()

###inflammatory factor references DOI:https://doi.org/10.1016/j.immuni.2018.03.023
immudata <- read.table('9.Inflammatory/TCR-BCR_mmc2.txt', sep = '\t', header = T, check.names = F)
LIHC <- subset(immudata, immudata$`TCGA Study` == 'LIHC')
riskScore <- read.table('5.CoxModule/trainRisk.txt', 
                        header = T, sep = '\t', check.names = F, row.names = 1)

sampleID <- substr(row.names(riskScore),1,12)
riskScore$`TCGA Participant Barcode` <- sampleID
riskScore$fullid <- row.names(riskScore)
riskmerge <- merge(riskScore, LIHC, by = 'TCGA Participant Barcode')
plotdata <- riskmerge[, c('fullid', 'risk', 'Leukocyte Fraction', 'Lymphocyte Infiltration Signature Score', 'IFN-gamma Response',
                          'BCR Richness', 'TCR Richness')]
my_comparisons <- list(c('Low','High'))
pdf(file="9.Inflammatory/IFN-γ.risk.pdf",width=6,height=6)
p=ggboxplot(plotdata, x="risk", y="IFN-gamma Response", merge = "flip",fill = 'risk', #color = "gender", add="jitter",
            ylab="IFN-gamma Response",
            xlab="",
            order = c('Low','High'),
            palette = c('#008B45','#FF0000','#FF00FF', '#FFD700','#9999FF', '#D1BBFF', '#EE7942','#FF88C2',  
                        '#FF6347','#33FFFF','#33FF33', '#FFA500', '#ADFF2F',
                        '#EE30A7', '#CCEEFF', '#0000AA','#770077','#9370DB'))
p=p+rotate_x_text(0)
p+stat_compare_means(aes(group=T),
                     method = 'wilcox.test',
                     label = "p.signif",
                     size = 6,
                     comparisons = my_comparisons,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  # stat_compare_means(size = 5, method = 'kruskal.test', label.y = 15, label.x = 0.8) +
  font('xlab',face = 'bold',size=14)+font('ylab',face = 'bold',size=14)+
  font('x.text',face = 'bold',size=14)+font('y.text',face = 'bold',size=14)+
  font('legend.title',face = 'bold',size=14)+font('legend.text',face = 'bold')
dev.off()

### CYT
CYT <- c('GZMA', 'PRF1')
CYTexp <- allgeneExpr[, CYT]
CYTexp$CYTscore <- apply(CYTexp, 1, mean)
CYTexp <- cbind(id = row.names(CYTexp), CYTexp)
riskScore <- read.table('5.CoxModule/trainRisk.txt', 
                        header = T, sep = '\t', check.names = F)
mergeCYT <- merge(CYTexp, riskScore, by = 'id')
my_comparisons <- list(c('Low','High'))
pdf(file="9.Inflammatory/CYT.risk.pdf",width=6,height=6)
p=ggboxplot(mergeCYT, x="risk", y="CYTscore", merge = "flip",fill = 'risk', #color = "gender", add="jitter",
            ylab="CYT Score",
            xlab="",
            order = c('Low','High'),
            palette = c('#008B45','#FF0000','#FF00FF', '#FFD700','#9999FF', '#D1BBFF', '#EE7942','#FF88C2',  
                        '#FF6347','#33FFFF','#33FF33', '#FFA500', '#ADFF2F',
                        '#EE30A7', '#CCEEFF', '#0000AA','#770077','#9370DB'))
p=p+rotate_x_text(0)
p+stat_compare_means(aes(group=T),
                     method = 'wilcox.test',
                     label = "p.signif",
                     size = 6,
                     comparisons = my_comparisons,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  # stat_compare_means(size = 5, method = 'kruskal.test', label.y = 15, label.x = 0.8) +
  font('xlab',face = 'bold',size=14)+font('ylab',face = 'bold',size=14)+
  font('x.text',face = 'bold',size=14)+font('y.text',face = 'bold',size=14)+
  font('legend.title',face = 'bold',size=14)+font('legend.text',face = 'bold')
dev.off()

###xCELL Immune #####

library(xCell)
library(ggpubr)
rm(list = ls())
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')
allgeneExpr <- merge2
allgeneExpr1 <- allgeneExpr[, -c(1:11)]
range(allgeneExpr1)
exprMatrix <- as.data.frame(t(allgeneExpr1))
xcellScore <- xCellAnalysis(exprMatrix, rnaseq = TRUE)
xcellout <- as.data.frame(xcellScore)
range(xcellout)
write.table(cbind(id = row.names(xcellout), xcellout), file = '9.Inflammatory/xCellout.txt', sep = '\t', row.names=F, quote=F)

riskScore <- read.table('5.CoxModule/trainRisk.txt', 
                        header = T, sep = '\t', check.names = F)
cell <- as.data.frame(t(xcellout))
cell <- cbind(id = row.names(cell), cell)
mergecell <- merge(cell, riskScore, by = 'id')
row.names(mergecell) <- mergecell$id
my_comparisons <- list(c('Low','High'))
pdf(file="9.Inflammatory/Fibroblasts.risk.pdf",width=6,height=6)
p=ggboxplot(mergecell, x="risk", y="Fibroblasts", merge = "flip",fill = 'risk', #color = "gender", add="jitter",
            ylab="Fibroblasts",
            xlab="",
            order = c('Low','High'),
            palette = c('#008B45','#FF0000','#FF00FF', '#FFD700','#9999FF', '#D1BBFF', '#EE7942','#FF88C2',  
                        '#FF6347','#33FFFF','#33FF33', '#FFA500', '#ADFF2F',
                        '#EE30A7', '#CCEEFF', '#0000AA','#770077','#9370DB'))
p=p+rotate_x_text(0)
p+stat_compare_means(aes(group=T),
                     method = 'wilcox.test',
                     label = "p.signif",
                     size = 6,
                     comparisons = my_comparisons,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  # stat_compare_means(size = 5, method = 'kruskal.test', label.y = 15, label.x = 0.8) +
  font('xlab',face = 'bold',size=14)+font('ylab',face = 'bold',size=14)+
  font('x.text',face = 'bold',size=14)+font('y.text',face = 'bold',size=14)+
  font('legend.title',face = 'bold',size=14)+font('legend.text',face = 'bold')
dev.off()

library(ggpubr)
rt0 <- read.table("9.Inflammatory/xCellout.txt",sep="\t",header=T,row.names=1,check.names=F) 
groups <- read.table("5.CoxModule/trainRisk.txt",
                     header = TRUE,sep="\t",check.names = F, row.names = 1)
rt=t(rt0[,row.names(groups)])
data1=cbind(rt, groups)

outTab=data.frame()
for(i in colnames(data1[,1:(ncol(data1)-13)])){
  rt1=data1[,c(i,"risk")]
  colnames(rt1)=c("Abundance","risk")
  outTab=rbind(outTab,cbind(rt1,cell=i))
}
write.table(outTab,file="9.Inflammatory/xCellout_forplotdata.txt",sep="\t",row.names=F,quote=F)

data2=read.table("9.Inflammatory/xCellout_forplotdata.txt",sep="\t",header=T,check.names=F)   
data2$risk=factor(data2$risk, levels=c("Low","High"))
p=ggboxplot(data2, x="cell", y="Abundance", fill = 'risk',
            ylab="Abundance",
            xlab="",
            palette = c('#00DAE0','#FF9289',"green","blue","red")) 
p=p+rotate_x_text(45)
pdf(file="9.Inflammatory/xCellout_boxplot.pdf",width=20,height=8)                        
p+stat_compare_means(aes(group=risk),
                     method = 'wilcox.test',
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")))
dev.off()

###Immune checkpoint#####
rm(list = ls())
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')
allgeneExpr <- merge2
allgeneExpr1 <- allgeneExpr[, -c(1:11)]
range(allgeneExpr1)

immugene <- read.table('10.Immune checkpoint/immuneCheckpoint.txt')
immugeneexpr <- allgeneExpr1[,immugene$V1]
write.table(cbind(id = row.names(immugeneexpr), immugeneexpr),file="10.Immune checkpoint/immuneCheckpointexpr.txt",sep="\t",row.names=F,quote=F)

group <- read.table('5.CoxModule/trainRisk.txt', 
                    header = T, check.names = F, sep = '\t', row.names = 1)
Type=group[row.names(immugeneexpr),]
data1=cbind(immugeneexpr,Type)
data1 <- data1[row.names(group),]
outTab=data.frame()
for(i in colnames(data1[,1:7])){
  rt1=data1[,c(i,"risk")]
  colnames(rt1)=c("Abundance","Risk")
  outTab=rbind(outTab,cbind(rt1,cell=i))
}
write.table(outTab,file="10.Immune checkpoint/immuneCheckpointAbundanceforPlot.txt",sep="\t",row.names=F,quote=F)

###pheatmap
library(pheatmap)
data2 <- data1[, c(1:7)]
data2 <- as.data.frame(t(data2))
group_list=c(rep('Low',181),rep('High',182))
annotation_col <- data.frame(Risk = factor(group_list,levels = c("Low","High")))
rownames(annotation_col) <- colnames(data2)
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = '10.Immune checkpoint/immuneCheckpoint.heatmap.pdf',width = 10,height = 4) 
pheatmap(data2,cellwidth = 1.3,cellheight = 20,
         method="pearson", 
         scale="row", 
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(color.key)(50),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         border_color = "NA")
dev.off()

##box plot
data2=read.table("10.Immune checkpoint/immuneCheckpointAbundanceforPlot.txt",sep="\t",header=T,check.names=F)   
data2$Risk=factor(data2$Risk, levels=c("Low","High"))
pdf(file="10.Immune checkpoint/immuneCheckpointBoxplot.pdf",width=12,height=8)  
p=ggboxplot(data2, x="cell", y="Abundance", color = "Risk", add="jitter",
            ylab="Abundance",
            xlab="",
            palette = c('#008B45','#FF0000','#FF00FF', '#FFD700','#9999FF', '#D1BBFF', '#EE7942','#FF88C2',  
                        '#FF6347','#33FFFF','#33FF33', '#FFA500', '#ADFF2F',
                        '#EE30A7', '#CCEEFF', '#0000AA','#770077','#9370DB'))
p=p+rotate_x_text(45)
p+stat_compare_means(aes(group=Risk),
                     method = 'wilcox.test',
                     label = "p.signif",
                     size = 6,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")))
dev.off()


###ssGSEA####
library(limma)
library(e1071)
library(preprocessCore)
library(tidyr)
library(ggplot2)
library(ggstatsplot)
library(ggExtra)
library(ggpubr)

rm(list = ls())
library(GSVA)
library(GSEABase)
load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMaxClin_onlyTumor.RData')
allgeneExpr <- merge2
allgeneExpr1 <- allgeneExpr[, -c(1:11)]
range(allgeneExpr1)

gmtFile="11.ssGSEA/mmc3.gmt"
rt=as.data.frame(t(allgeneExpr1))
rt=as.matrix(rt)

geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

ssgseaScore=gsva(rt, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
range(ssgseaScore)
ssgseaOut <- as.data.frame(ssgseaScore)
write.table(cbind(id = row.names(ssgseaOut), ssgseaOut),file="11.ssGSEA/ssgseaOut.txt",sep="\t",quote=F,col.names=T, row.names = F)

###### heatmap plot
library(pheatmap)    

outdata <- read.table('11.ssGSEA/ssgseaOut.txt',header=T,sep="\t",check.names=F, row.names = 1)
group <- read.table('5.CoxModule/trainRisk.txt',
                    header = T, check.names = F, sep = '\t', row.names = 1)
outdata1 <- as.data.frame(t(outdata))
plotdata <- outdata1[row.names(group),]
plotdata1 <- as.data.frame(t(plotdata))
library(pheatmap)
group_list=c(rep('Low',181),rep('High',182))
annotation_col <- data.frame(risk = factor(group_list,levels = c("Low","High")))
rownames(annotation_col) <- colnames(plotdata1)
pdf("11.ssGSEA/ssGSEAscoreheatmap.pdf",height=6,width=12)
pheatmap(plotdata1,
         annotation=annotation_col,
         cellwidth = 1.5,cellheight = 12,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = T,
         cluster_cols = F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         # gaps_col = c(247,383),
         show_colnames=F,show_rownames =T,
         treeheight_row = "20",treeheight_col = "0",
         fontsize_col=3)
dev.off()

###boxplot
library(ggpubr)
rt0 <- read.table("11.ssGSEA/ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F) 
groups <- read.table("5.CoxModule/trainRisk.txt",
                     header = TRUE,sep="\t",check.names = F, row.names = 1)
rt=t(rt0[,row.names(groups)])
data1=cbind(rt, groups)

outTab=data.frame()
for(i in colnames(data1[,1:28])){
  rt1=data1[,c(i,"risk")]
  colnames(rt1)=c("Abundance","risk")
  outTab=rbind(outTab,cbind(rt1,cell=i))
}
write.table(outTab,file="11.ssGSEA/ssgseaOut_forplotdata.txt",sep="\t",row.names=F,quote=F)

data2=read.table("11.ssGSEA/ssgseaOut_forplotdata.txt",sep="\t",header=T,check.names=F)   
data2$risk=factor(data2$risk, levels=c("Low","High"))
p=ggboxplot(data2, x="cell", y="Abundance", fill = 'risk',
            ylab="Abundance",
            xlab="",
            palette = c('#00DAE0','#FF9289',"green","blue","red")) 
p=p+rotate_x_text(45)
pdf(file="11.ssGSEA/ssgseaOut_boxplot.pdf",width=12,height=8)                        
p+stat_compare_means(aes(group=risk),
                     method = 'wilcox.test',
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")))
dev.off()


##riskscore and immu cell
library("psych")
data3 <- data1[, c(1:28,36)]
Norm.interest.corr <- corr.test(data3, method="spearman", ci=F)
Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))
Correlation <- as.data.frame(as.table(Norm.interest.corr$r))
Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]
colnames(Cor.table) <- c("Cell1","Cell2","cor","p.value")
write.table(Cor.table,file="11.ssGSEA/ssGSEAimmucellRiskCor.txt",sep="\t",quote=F,row.names=F)

immurisk <- subset(Cor.table, Cor.table$Cell2 == 'riskScore')
immurisk <- subset(immurisk, immurisk$Cell1 != 'riskScore')

immurisk$pstar <- ifelse(immurisk$p.value < 0.05,
                     ifelse(immurisk$p.value < 0.01,"**","*"),
                     "")
pdf(file="11.ssGSEA/ssGSEAimmucellRiskCorriskScore.pdf",width=8,height=3)
ggplot(immurisk, aes(Cell1, Cell2)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 6)+
  theme_minimal()+ theme(panel.grid = element_blank(), 
  )+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
dev.off()


###DEAREGene boxpolt in different set####
####plot input data deal 
probeExpr <- read.table('12.DEAREGene boxpolt/GSE76427_probeExpr.txt', header = T, sep = '\t', check.names = F)
GPL <- read.table('12.DEAREGene boxpolt/GPL10558-50081proble-symbol.txt', header = T, sep = '\t', check.names = F)
symbolExpr <- merge(GPL, probeExpr, by = 'ID_REF')
symbolExpr <- na.omit(symbolExpr)
symbolExpr <- symbolExpr[, -1]
range(symbolExpr[,-1])
write.table(symbolExpr, file = '12.DEAREGene boxpolt/GSE76427_symbolExpr.txt', sep="\t", row.names = F, quote = F)

gene<-read.table("12.DEAREGene boxpolt/GSE76427_symbolExpr.txt",sep="\t",header=T, check.names = F)
data3 <- aggregate( . ~ id, data=gene, max)  
range(data3[,-1])
row.names(data3) <- data3$id
data3 <- data3[,-1]
data4 <- log2(data3)
range(data4)
write.table(cbind(id=row.names(data4), data4),file="12.DEAREGene boxpolt/GSE76427_symbolExprUniMax.txt",sep="\t",quote = F,row.names=F)

load('1.TCGA data deal/LIHC_mRNA.symbolExprUniMax.Rdata')
Expr <- exprSet2
Expr <- read.table('1.TCGA data deal/LIHC_mRNA.symbolExprUniMax.txt',header=T,sep="\t",check.names=F,row.names = 1)
aimGene <- read.table('5.CoxModule/multiCox.txt', header = T, sep = '\t', check.names = F,row.names = 1)
aimGeneexpr <- Expr[row.names(aimGene),]
aimGeneexpr <- na.omit(aimGeneexpr)
range(aimGeneexpr)
aimGeneexpr <- log2(aimGeneexpr +1)
aimGeneexpr2 <- as.data.frame(t(aimGeneexpr))
write.table(cbind(id = row.names(aimGeneexpr),aimGeneexpr), file = '12.DEAREGene boxpolt/TCGAaimGeneexpr.txt', sep = '\t', quote = F, row.names = F)

##boxplot
library(ggpubr)
rt0 <- read.table("12.DEAREGene boxpolt/TCGAaimGeneexpr.txt",
                  sep="\t",header=T,row.names=1,check.names=F) 
groups <- read.table("1.TCGA data deal/LIHC_sampleType.txt",
                     header = TRUE,sep="\t",check.names = F,row.names = 1)
range(rt0)
rt=t(rt0[,row.names(groups)])
data1=cbind(rt, groups)

outTab=data.frame()
for(i in colnames(data1[,1:(ncol(data1)-1)])){
  rt2=data1[,c(i,"type")]
  colnames(rt2)=c("Abundance","cluster")
  outTab=rbind(outTab,cbind(rt2,gene=i))
}
write.table(outTab,file="12.DEAREGene boxpolt/TCGAaimGeneexpr_forplotdata.txt",sep="\t",row.names=F,quote=F)
data2=read.table("12.DEAREGene boxpolt/TCGAaimGeneexpr_forplotdata.txt",sep="\t",header=T,check.names=F) 
pdf(file="12.DEAREGene boxpolt/TCGAaimGeneexpr_boxplot.pdf",width=10,height=6)
p=ggboxplot(data2, x="gene", y="Abundance", fill = 'cluster',
            ylab="Abundance",
            xlab="",#ylim=c(0,12),
            palette = c("#00DAE0","#FF9289","red")) 
p=p+rotate_x_text(45)
p+stat_compare_means(aes(group=cluster),
                     method = 'wilcox.test',
                     label = "p.signif", #label.y = 12,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")))
dev.off()

###DiffAREG venn Hypergeometric Distribution#####
AREGene <- read.table('3.DiffAREG/AREGene.txt')
DEG <- read.table('2.Diffexpr/TNlimmaOut.txt', header = T, check.names = F, sep = '\t')
biotype <- read.csv('13.Venn Hypergeometric Distribution/mart_export.csv', header = T, check.names = F)
mergetype <- merge(DEG, biotype, by.x = 'Symbol', by.y = 'Gene name')
mRNA <- subset(mergetype, mergetype$`Gene type`== "protein_coding")
length(unique(mRNA$Symbol))
intersectARG <- intersect(AREGene$V1, unique(mRNA$Symbol))
DEGsig <- read.table('2.Diffexpr/TNlimmadiffSig.txt', header = T, check.names = F, sep = '\t')
DEGsig1 <- merge(DEGsig, biotype, by.x = 'Symbol', by.y = 'Gene name')
DEGsig1 <- subset(DEGsig1, DEGsig1$`Gene type`== "protein_coding")
intersectARGsig <- intersect(AREGene$V1, unique(DEGsig1$Symbol))

n=length(unique(DEGsig1$Symbol))
k=length(intersectARGsig)
N=length(unique(mRNA$Symbol))
M=length(intersectARG)

pvalue <- phyper(k, M, N-M, n, lower.tail=FALSE)
pvalue

### Immune cell venn Hypergeometric Distribution####
n=11
k=3
N=28
M=3

pvalue <- phyper(k, M, N-M, n, lower.tail=FALSE)
pvalue
