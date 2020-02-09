# BLCA signature
# Data proprocessing
# input RNA file 
rna <- read.table('HiSeq.txt', header=T, row.names=1, sep='\t')
# and take off first row cause we don't need it
rna <- rna[-1,]

# first, remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]

library(dplyr)
rna <- rna %>% 
  group_by(ID) %>% 
  summarise_all(mean)

# ID2symbol
idmap <- read.delim("gencode.v23.annotation.gene.probemap",as.is=T)
head(idmap)
rownames(idmap) <- idmap$id
rna$gsym <- idmap[rownames(rna),]$gene
write.csv(rna,"easy_input_expr.csv", quote = F)

# protein-coding genes, long non-coding RNAs
library(biomaRt)
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl", 
                   host = "www.ensembl.org")
listAttributes(ensembl)[62,]
feature_info <- getBM(attributes = c("gene_biotype",
                      "hgnc_symbol"), 
                      filters = "hgnc_symbol", 
                      values = rownames(exprdf_uniq), mart = ensembl)

if (nrow(exprdf_uniq) != nrow(feature_info)){
  library(dplyr)
  diffName<-setdiff(rownames(exprdf_uniq),feature_info[,2])
  length(diffName)
  head(diffName)
}
write.csv(feature_info[,c(2,1)],"gene_biotype.csv",quote = F,row.names = F)

# gene biotype
unique(feature_info$gene_biotype)

mRNA <-"protein_coding"
lncRNA <- paste("non_coding","3prime_overlapping_ncRNA","antisense","lincRNA","sense_intronic","sense_overlapping","macro_lncRNA","bidirectional_promoter_lncRNA",sep = "|")

mRNA.list<-feature_info[grepl(mRNA, feature_info$gene_biotype),]
write.table(mRNA.list,"mRNA.list.txt",quote = F,row.names = F, col.names = F)

lncRNA.list<-feature_info[grepl(lncRNA, feature_info$gene_biotype),]
write.table(lncRNA.list,"lncRNA.list.txt",quote = F,row.names = F, col.names = F)

exp_sym <- read.csv('easy_input_expr.csv', header=T, row.names=1)

mRNA_expr <- exprdf_uniq[as.character(mRNA.list[,2]),]
write.csv(mRNA_expr,"mRNA_expr.csv",quote = F,row.names = T)

lncRNA_expr <- exprdf_uniq[as.character(lncRNA.list[,2]),]
write.csv(lncRNA_expr,"lncRNA_expr.csv",quote = F,row.names = T)

set.seed(1234567)
sub<-sample(1:400, 200)
mRNA_train<-mRNA_expr[,sub]
mRNA_test<-mRNA_expr[,-sub]

lncRNA_train<-lncRNA_expr[,sub]
lncRNA_test<-lncRNA_expr[,-sub]

miRNA_train<-miRNA_expr[,sub]
miRNA_test<-miRNA_expr[,-sub]

# Selection of prognostic RNAs
# Univariate Cox analysis
outTab=data.frame()

library(survival)
rt<-read.table("clinicalExp.txt", header = T, sep = "\t", row.names = 1, check.names = F)
rt1<-log2(rt[,3:ncol(rt)]+1)
rt<-cbind(rt[,1:2],rt1)
rt[,"futime"]=rt[,"futime"]/365

for(i in colnames(rt[,3:ncol(rt)])){
  cox<-coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  outTab = rbind(outTab,cbind(gene=i, 
                              z=coxSummary$coefficients[,"z"],
                              coef=coxSummary$coefficients[,"coef"],
                              pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                              HR=coxSummary$coefficients[,"exp(coef)"],
                              HRCILL=coxSummary$conf.int[,"lower .95"],
                              HRCIUL=coxSummary$conf.int[,"upper .95"]))
}

write.table(outTab,file = "univariateCox.txt",sep = "\t", row.names = F, quote = F)

# Multivariate analysis
library(survminer)

rt<-read.table("multiInput.txt", header = T, sep = "\t", row.names = 1, check.names = F)
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
rt[,"futime"]=rt[,"futime"]/365

cox<-coxph(Surv(futime, fustat) ~ ., data = rt)
cox<-step(cox,direction = "both")
riskScore<-predict(cox,type="risk", newdata=rt)
summary<-summary(cox)

coef=summary[[7]][,1] 
p=summary[[7]][,5]
coef2=coef[which(p<0.05)]
marker<-names(coef2)

rt2<-rt[,marker]
rt<-cbind(rt[,1:2],rt2)

coxGene=rownames(summary$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore), "high", "low"))
write.table(cbind(id=rownames(cbind(rt,riskScore,risk)), cbind(rt,riskScore,risk)),
            file = "risk.txt",
            sep = "\t",
            quote = F,
            row.names = F
)
write.table(cbind(id=coxGene,summary$coefficients),
            file = "coxResult.xls",
            sep = "\t",
            quote = F,
            row.names = F
)
coxGene

## ggsurvplot
library(survival)
library(survminer)

dat<-read.table("risk.txt", header = T, sep = "\t", row.names = 1, check.names = F)
time<-as.numeric(dat$futime)
status<-as.numeric(dat$fustat)

values <- dat[,13]
group<-ifelse(values>=median(values), 'high', 'low')
fit2 <- survfit(Surv(time, status) ~ group, data = dat)
ggsurvplot(fit2,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.style = "ribbon",
           xlab = "Time in years",
           ggtheme = theme_light(),
           risk.table = "abs_pct",
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           ncensor.plot = TRUE,
           surv.median.line = "hv",
           legend.labs = c("High_risk", "Low_risk"),
           palette = c("red", "blue")
)

fit <- survfit(Surv(time, event) ~ cluster,
               data = dat)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,15),        # present narrower X axis, but not affect
  # survival estimates. 
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

## risk score
classi_risk = as.numeric(ifelse(values>=median(values), 1, 0))
Cox_risk = coxph(Surv(time, status)~classi_risk)
summary(Cox)

library(ggplot2)
survdata <- read.table("exp_5.txt",sep = "\t",header=T,stringsAsFactors = F,na.strings = "")
survdata$state <- factor(survdata$status_type)
sg <- ggplot(survdata,aes(x=seq,y=os_years,color=state))+geom_point(shape=21,sort(aes(survdata$seq)))+ylab("Survival time in years")    
final <- sg + theme(axis.line = element_line(color = "black"),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank())
p = final + labs(color="state")+scale_color_manual(values = c("red","blue"))
p + geom_vline(aes(xintercept=100), colour="#BB0000", linetype="dashed")

# heatmap for the five RNAs
df<-read.table("exp_5_200_seq.txt",row.names = 1,header = T,as.is = T)
library(pheatmap)
annotation_col = data.frame(
  Risk = rep(c("Low_risk","High_risk"),c(100,100)),
  row.names = colnames(df)
)

ann_colors = list(
  Risk = c(Low_risk = "blue", High_risk = "red")
)

# Combined the clinical data
clinical <- read.csv("pro_simple_clinical.csv", as.is = T)
## Age
age <- clinical[,4]
classi_age = as.numeric(ifelse(age > median(age), 1, 0))
Cox = coxph(Surv(time, status)~classi_age)
cox_age = c(summary(Cox)$coefficients[c(5,2)], summary(Cox)$conf.int[c(3,4)])
cox_age
summary(Cox)

## Gender
classi_gender = as.numeric(ifelse(clinical$gender == "FEMALE", 1, 0))
Cox = coxph(Surv(time, status)~classi_gender)
cox_gender = c(summary(Cox)$coefficients[c(5,2)], summary(Cox)$conf.int[c(3,4)])
cox_gender
summary(Cox)

## Clinical Stage
classi_stage = as.numeric(ifelse(clinical$stage %in% unlist(lapply(c("Stage III", "Stage IV", "Stage X"), function(s) paste(s, c("", "A", "B", "C", "D"), sep = ""))), 1, 
                                 ifelse(clinical$stage %in% unlist(lapply(c("Stage I", "Stage II"), function(s) paste(s, c("", "A", "B", "C"), sep = ""))), 0, NA)))
Cox = coxph(Surv(time, status)~classi_stage)
cox_stage = c(summary(Cox)$coefficients[c(5,2)], summary(Cox)$conf.int[c(3,4)])
cox_stage
summary(Cox)

## pM
pM <- clinical[,5]
classi_M = as.numeric(ifelse(pM == "1", 1, 0))
Cox = coxph(Surv(time, status)~classi_M)
cox_M = c(summary(Cox)$coefficients[c(5,2)], summary(Cox)$conf.int[c(3,4)])
cox_M
summary(Cox)

## multi_cox
Cox = coxph(Surv(time, status)~classi_age + classi_M + classi_N + classi_stage + classi_risk)
summary(Cox)

# Robust likelihood-based survival modeling
library(rbsurv)
fit <- rbsurv(time=futime, status=fustat, x=rt, method="efron", max.n.genes=20, n.iter=100, n.fold=3)
fit$model

# Random Forests for Survival, Regression, and Classification (RF-SRC)
obj <- rfsrc(Surv(futime, fustat) ~ ., data = rt, ntree = 100, tree.err=TRUE, importance = TRUE)
plot(obj)

# permutation and combination method 
dat = read.csv('gene_nine.csv')
rownames(dat) = as.vector(dat[, 1])
dat = dat[, -1]
os1 = read.csv('os.csv')
os = os1[, 2] 
group = apply(dat, 1, function(values1){
  group=ifelse(values1>median(values1),'high','low')})

C = vector(9, mode = 'list')
for (i in 1:9) {
  C[[i]] = combn(9, i)
}


library(survival)
score_res = matrix(NA, 200, 511) 
score_res = as.data.frame(score_res)
l = 1
for (i in 1:9) {
  for (j in 1:dim(C[[i]])[2]) {
    my.surv <- Surv(os,rep(1,length(os)))
    name = colnames(group)[c(C[[i]][, j])]
    m = coxph(my.surv ~ ., data =  data.frame(my.surv, group[, c(C[[i]][, j])]))
    coe = coef(m)
    score = rep(0, dim(dat)[2])
    for (k in 1:dim(C[[i]])[1]) {
      score = dat[C[[i]][k, j], ] * coe[k] + score
    }
    score_res[, l] = t(score)
    colnames(score_res)[l] = paste0(name, collapse = "_")
    l = l + 1
  }
}

rownames(score_res) = colnames(dat)

write.csv(score_res, file = 'res.csv')

require(pROC)
outcome = os1[, 3]
auc = c()
for (i in 1:511) {
  score = score_res[, i]
  rocobj = roc(outcome, score)
  auc[i] = auc(rocobj)
}

write.csv(auc, file = 'auc.csv')

# NMF
data = read.table("tumor_train_mRNA_8969_200.txt",sep='\t',header=T,row.names=1)
datExpr = t(data[order(apply(data,1,mad), decreasing = T)[1:500],])
mRNA_eset<-t(datExpr)
probeid<-rownames(mRNA_eset) 
tumor_exprs<-cbind(probeid, mRNA_eset) 
write.table(tumor_exprs,file="mRNA_eset.txt",sep='\t',quote=F,row.names=F)
# x: Rows are genes, colums are samples. There is no rownames and colnames.
x <- read.table ("0_mRNA_eset.txt",header=FALSE,sep="\t")
res <- nmf(x, 2:10, method='brunet', nrun = 30, seed = 123456)
plot(res)
# Estimation of the rank: Consensus matrices computed from 30 runs for each value of r.
consensusmap(res, annCol = x, labCol = NA, labRow = NA)
#select the k value (k = 3) of the max cophenetic correlation coefficient
res_3 <- nmf(x, 3)
# Compute the featureScore and extractFeatures
s <- featureScore(res_3)
summary(s)
s1 <- extractFeatures(res_3, method="max")
str(s1)
mRNA = read.table("mRNA.txt",sep='\t',header=T,row.names=1)
q<-rownames(mRNA)
meta2<-mRNA[y,]
probeid2<-rownames(meta2) 
mRNA<-cbind(probeid2, meta2) 
write.table(mRNA,file="metagene_eset_max.txt",sep='\t',quote=F,row.names=F)

# CancerSubtypes
library(CancerSubtypes)
n <- read.table ("tumor_train_mRNA_8969_200.txt", header=TRUE,sep="\t")
n <- as.matrix(n)
b <- read.table ("0_tumor_train_mRNA_8969_200.txt",header=FALSE,sep="\t")
b <- as.matrix(b)

rownames(b) <- n[,1]
col<-read.table("colnames.txt",sep='\t',header=T,row.names=1)
colnames(b)<- rownames(col)

GeneExp<-read.table("pro_metagene_eset_mRNA_172_200.txt",header=T,sep="\t",row.names=1)
lncRNAExp<-read.table("pro_metagene_eset_lncRNA_15_200.txt",header=T,sep="\t",row.names=1)
miRNAExp<-read.table("pro_metagene_eset_miRNA_42_200.txt",header=T,sep="\t",row.names=1)
a<-as.matrix(GeneExp)
b<-as.matrix(lncRNAExp)
c<-as.matrix(miRNAExp)
clin<-read.table("clindata.txt",header=T,sep="\t",row.names=1)

Tumor=list(GeneExp=GeneExp,lncRNAExp=lncRNAExp,miRNAExp=miRNAExp)
result=ExecuteSNF(Tumor, clusterNum=3, K=20, alpha=0.5, t=20)
group=result$group
distanceMatrix=result$distanceMatrix
silhouette=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
silhouette
plot(silhouette, col = c("red", "orange", "purple"))

drawHeatmap(a,group,silhouette=silhouette,color = colorRampPalette(c("red", "white", "purple"))(300), scale="max_min",Title="Bladder Cancer Gene Expression")
drawHeatmap(a,group,silhouette=silhouette,scale="max_min",
            color="-RdYlBu",Title="Bladder Cancer Gene Expression")

# A statistical method for testing the significance of clustering results.
sigclust1=sigclustTest(miRNAExp,group, nsim=500, nrep=1, icovest=3)
sigclust2=sigclustTest(miRNAExp,group, nsim=1000, nrep=1, icovest=1)
p_value=survAnalysis(mainTitle="Bladder_cancer_SNF",time,status,group,distanceMatrix=distanceMatrix,similarity=TRUE)

# heatmap
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(rtracklayer)
library(SummarizedExperiment)
library(clusterProfiler)
library(RColorBrewer)
library(maftools)
library(circlize)
library(matrixStats)
library(GetoptLong)
library(GenomicRanges)
library(TCGAbiolinks)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 


targetMut <- c("TP53", "KRAS", "BRAF", "EGFR") 
nonsilentmutation <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation") 

blcamut <- GDCquery_Maf(tumor = "BLCA", pipelines = "varscan2")

crcmut <- data.frame(sample = str_sub(blcamut$Tumor_Sample_Barcode, 1, 12), blcamut)
all_mutSample <- unique(crcmut$sample)

table(crcmut$Variant_Classification)
crcmut <- crcmut[crcmut$Variant_Classification %in% nonsilentmutation, ]

crc_mut_gene <- crcmut[crcmut$Hugo_Symbol %in% targetMut,]

TP53mutSample <- unique(crc_mut_gene[crc_mut_gene$Hugo_Symbol == "TP53", "sample"])  
TP53wildSample <- setdiff(all_mutSample, TP53mutSample) 

KRASmutSample <- unique(crc_mut_gene[crc_mut_gene$Hugo_Symbol == "KRAS", "sample"])
KRASwildSample <- setdiff(all_mutSample, KRASmutSample)

BRAFmutSample <- unique(crc_mut_gene[crc_mut_gene$Hugo_Symbol == "BRAF", "sample"])
BRAFwildSample <- setdiff(all_mutSample, BRAFmutSample)

EGFRmutSample <- unique(crc_mut_gene[crc_mut_gene$Hugo_Symbol == "EGFR", "sample"])
EGFRwildSample <- setdiff(all_mutSample, EGFRmutSample)

targetAnno<-read.table("targetAnno.txt", header=TRUE, sep="\t")

targetAnno[targetAnno$ID %in% TP53mutSample, "TP53mut"] <- "mutant"
targetAnno[targetAnno$ID %in% TP53wildSample, "TP53mut"] <- "wildtype" 

targetAnno[targetAnno$ID %in% KRASmutSample, "KRASmut"] <- "mutant"
targetAnno[targetAnno$ID %in% KRASwildSample, "KRASmut"] <- "wildtype" 

targetAnno[targetAnno$ID %in% BRAFmutSample, "BRAFmut"] <- "mutant"
targetAnno[targetAnno$ID %in% BRAFwildSample, "BRAFmut"] <- "wildtype" 

targetAnno[targetAnno$ID %in% EGFRmutSample, "EGFRmut"] <- "mutant"
targetAnno[targetAnno$ID %in% EGFRwildSample, "EGFRmut"] <- "wildtype"

summary(targetAnno)
write.table(targetAnno,file="targetAnno_1.txt",sep='\t',quote=F,row.names=F)

## 'heatmapinput'
setwd("C:/Users/liudanqi/Desktop/bind/heatmap")
heatmapinput <- read.csv("heatmapinput.csv", row.names = 1)
ml <- heatmapinput[, c(2:230)] 
ml <- as.data.frame(t(apply(ml, 2, scale))) 
colnames(ml) <- rownames(heatmapinput)

col_fun <- colorRamp2(c(-5, 0, 5), c("#377EB8", "white", "#E41A1C"))

h1 <- Heatmap(ml, cluster_rows = TRUE, cluster_columns = TRUE, clustering_method_columns = "ward.D2",show_row_names = FALSE, show_column_names = FALSE,
              clustering_distance_columns = "euclidean", 
              clustering_distance_rows = "euclidean",
              clustering_method_rows  = "ward.D2")
tree <- column_dend(h1)
ind <- cutree(as.hclust(tree), k = 3)[order.dendrogram(tree)] 
table(ind)

heatmapinput$Cluster <- ind[heatmapinput$barcode]

heatmapinput$Cluster <- str_replace(heatmapinput$Cluster, "1", "Cluster 1")
heatmapinput$Cluster <- str_replace(heatmapinput$Cluster, "2", "Cluster 2")
heatmapinput$Cluster <- str_replace(heatmapinput$Cluster, "3", "Cluster 3")

# annotation of clinical information 
Cluster <- heatmapinput[, "Cluster"]
mRNA.cluster <- heatmapinput[, "mRNA.cluster"]
lncRNA.cluster <- heatmapinput[, "lncRNA.cluster"]
microRNA.cluster <- heatmapinput[, "microRNA.cluster"]
RPPA.cluster <- heatmapinput[, "RPPA.cluster"]
Hypermethylation.cluster <- heatmapinput[, "Hypermethylation.cluster"]
Hypomethylation.cluster <- heatmapinput[, "Hypomethylation.cluster"]
Mutation.process.cluster <- heatmapinput[, "Mutation.process.cluster"]
SMG.SCNA.cluster <- heatmapinput[, "SMG.SCNA.cluster"]
Histological.subtype <- heatmapinput[, "Histological.subtype"]
TP53_mutation <- heatmapinput[, "TP53mut"]
KRAS_mutation <- heatmapinput[, "KRASmut"]
BRAF_mutation <- heatmapinput[, "BRAFmut"]
EGFR_mutation <- heatmapinput[, "EGFRmut"]
Survival <- heatmapinput[, "Survival"] 
Gender <- heatmapinput[, "Gender"] 
pstage <- heatmapinput[, "pstage"]
Histological_grade <- heatmapinput[, "Histological_grade"]
pT <- heatmapinput[, "pT"]
pN <- heatmapinput[, "pN"]
pM <- heatmapinput[, "pM"]

ha = HeatmapAnnotation(Cluster = Cluster, mRNA.cluster = mRNA.cluster, lncRNA.cluster = lncRNA.cluster, microRNA.cluster = microRNA.cluster, RPPA.cluster = RPPA.cluster, Hypermethylation.cluster = Hypermethylation.cluster, Hypomethylation.cluster = Hypomethylation.cluster, Mutation.process.cluster = Mutation.process.cluster, SMG.SCNA.cluster = SMG.SCNA.cluster, Histological.subtype = Histological.subtype, TP53_mutation = TP53_mutation, KRAS_mutation = KRAS_mutation, BRAF_mutation = BRAF_mutation, EGFR_mutation = EGFR_mutation,
                       Survival = Survival, Gender = Gender, pstage = pstage, Histological_grade = Histological_grade, pT = pT, pN = pN, pM = pM, 
                       col = list(Cluster = c("Cluster 1" = "#3FA538", "Cluster 2" = "#9FD29BFF", "Cluster 3" = "#C4868E"), mRNA.cluster = c("Basal_squamous" = "#E5554D", "Luminal" = "#C4868E", "Luminal_infiltrated" = "#AEB6CE", "Luminal_papillary" = "#2B3D44", "Neuronal" = "#FCF732FF"), lncRNA.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44", "NA" = "#FCF732FF"), microRNA.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44"), RPPA.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44", "5" = "#FCF732FF", "ND" = "#3FA538"), Hypermethylation.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44", "5" = "#FCF732FF"), Hypomethylation.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44", "5" = "#FCF732FF"), Mutation.process.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44", "NA" = "#FCF732FF"), SMG.SCNA.cluster = c("1" = "#E5554D", "2" = "#C4868E", "3" = "#AEB6CE", "4" = "#2B3D44", "NA" = "#FCF732FF"), 
                                  Histological.subtype = c("Papillary" = "#C4868E", "Non-Papillary" = "#97A8C7", "ND" = "#B0B0FFFF"), 
                                  TP53_mutation = c("mutant" = "black", "wildtype" = "grey"),KRAS_mutation = c("mutant" = "black", "wildtype" = "grey"),BRAF_mutation = c("mutant" = "black", "wildtype" = "grey"),EGFR_mutation = c("mutant" = "black", "wildtype" = "grey"),
                                  Survival = c("Alive" = "#3FA538", "Dead" = "#E00115"), 
                                  Gender = c("MALE" = "#E5554D", "FEMALE" = "#C4868E"), 
                                  pstage =  c("Stage I" = "#B0B0FFFF", "Stage II" = "#B0FFB0FF", 
                                              "Stage III" = "#F7E897FF",
                                              "Stage IV" = "#FF6060FF"),
                                  Histological_grade =  c("High Grade" = "#B0B0FFFF", "Low Grade" = "#6060FFFF"),
                                  pT = c("T0" = "#F7E897FF", "T1" = "#E5554D", "T2" = "#C4868E", "T3" = "#AEB6CE", "T4" = "#2B3D44", "TX" = "#FCF732FF"),
                                  pN = c("N0" = "#E5554D", "N1" = "#C4868E", "N2" = "#AEB6CE", "N3" = "#2B3D44", "NX" = "#FCF732FF"),
                                  pM = c("M0" = "#E5554D", "M1" = "#C4868E", "MX" = "#AEB6CE")),na_col = "white",
                       show_legend = rep(TRUE, 21),
                       annotation_height = unit(rep(5, 21), "mm"),
                       annotation_legend_param = list(
                         Cluster = list(title = "Cluster"),
                         mRNA.cluster = list(title = "mRNA.cluster"),
                         lncRNA.cluster = list(title = "lncRNA.cluster"),
                         microRNA.cluster = list(title = "microRNA.cluster"),
                         RPPA.cluster = list(title = "RPPA.cluster"),
                         Hypermethylation.cluster = list(title = "Hypermethylation.cluster"),
                         Hypomethylation.cluster = list(title = "Hypomethylation.cluster"),
                         Mutation.process.cluster = list(title = "Mutation.process.cluster"),
                         SMG.SCNA.cluster = list(title = "SMG.SCNA.cluster"),
                         Histological.subtype = list(title = "Histological.subtype"),
                         TP53_mutation = list(title = "TP53 mutation"),
                         KRAS_mutation = list(title = "KRAS mutation"),
                         BRAF_mutation = list(title = "BRAF mutation"),
                         EGFR_mutation = list(title = "EGFR mutation"),
                         Survival = list(title = "Survival"),
                         Gender = list(title = "Gender"),
                         pstage = list(title = "pstage"),
                         Histological_grade = list(title = "Histological_grade"),
                         pT = list(title = "pT"),
                         pN = list(title = "pN"),
                         pM = list(title = "pM")
                       ))

ht <- Heatmap(ml, col = col_fun, 
              name = "TCGA_BLCA 400 samples",
              cluster_rows = TRUE, cluster_columns = TRUE,          
              show_row_names = FALSE, show_column_names = FALSE,
              bottom_annotation = ha, column_title = qq("TCGA_BLCA 400 samples (n = @{ncol(ml)})"),
              clustering_method_columns = "ward.D2",
              clustering_distance_columns = "euclidean", 
              clustering_distance_rows = "euclidean",
              clustering_method_rows  = "ward.D2", column_dend_height = unit(30, "mm")
)

pdf("heatmap.pdf", 30, 26)
draw(ht, annotation_legend_side = "left", heatmap_legend_side = "right")

annotation_titles <- c(
  Cluster = "Cluster",
  mRNA.cluster = "mRNA.cluster",
  lncRNA.cluster = "lncRNA.cluster",
  microRNA.cluster = "microRNA.cluster",
  RPPA.cluster = "RPPA.cluster",
  Hypermethylation.cluster = "Hypermethylation.cluster",
  Hypomethylation.cluster = "Hypomethylation.cluster",
  Mutation.process.cluster = "Mutation.process.cluster",
  SMG.SCNA.cluster = "SMG.SCNA.cluster",
  Histological.subtype = "Histological.subtype",
  TP53_mutation = "TP53 mutation",
  KRAS_mutation = "KRAS mutation",
  BRAF_mutation = "BRAF mutation",
  EGFR_mutation = "EGFR mutation",
  Survival = "Survival",
  Gender = "Gender",
  pstage = "Stage", 
  Histological_grade = "Histological_grade",
  pT = "pT",
  pN = "pN",
  pM = "pM"
)
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

# Functional annotation 
# GO and KEGG analyses
library(org.Hs.eg.db)
library(clusterProfiler)
data<-read.table(file = "cor_p.txt",row.names = 1,header = T, sep = "\t")
genelist<-rownames(data)
s.EntrezID<- bitr(genelist, fromType="SYMBOL", toType="ENTREZID",OrgDb="org.Hs.eg.db")
head(s.EntrezID)
dim(s.EntrezID)

s.ego <- enrichGO(gene = s.EntrezID$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = 'ALL', pAdjustMethod = 'fdr',pvalueCutoff = 0.05, qvalueCutoff = 0.05, keyType = 'ENTREZID')
head(s.ego)
go.gsym <- setReadable(s.ego, 'org.Hs.eg.db', 'ENTREZID')
x<-go.gsym@result
id<-rownames(x) 
tumor_exprs<-cbind(id, x) 
write.table(tumor_exprs,file="s.ego.txt",sep='\t',quote=F,row.names=F)
save(s.ego, go.gsym, file = "ego.RData")
class(go.gsym)
str(go.gsym)

## Integration of GO analysis
go_dot<-dotplot(go.gsym,split="ONTOLOGY",showCategory=5)+facet_grid(ONTOLOGY~.)
go_dot

# KEGG analysis
kk <- enrichKEGG(gene         = s.EntrezID$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05
)
head(kk)
kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
# save the file
write.csv(kk.gsym,"kegg.csv", quote = F, row.names = F)
save(kk, kk.gsym, file = "kk_KEGG.RData")
## dotplot
dotplot(kk)

# GSVA 
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

# input the file
msigdbr_show_species()

h <- msigdbr(species = "Homo sapiens",
             category = "H") # hallmark gene set

h <- select(h, gs_name, gene_symbol) %>% 
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol)) 

# Remove duplicate genes in each geneset
gs <- lapply(h, unique)

# Remove genes that have appeared in two or more pathways
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))

gs <- gs[lapply(gs, length) > 0]
head(gs)
# save the file
save(gs, file = "hallmark.gs.RData")

## GSVA分析
(load("hallmark.gs.RData")) 

# input the gene expression matrix
gsym.expr <- read.csv("easy_input_expr.csv", row.names = 1)
head(gsym.expr)

# GSVA analysis
gsva_es <- gsva(as.matrix(gsym.expr), gs)
pheatmap::pheatmap(gsva_es)
head(gsva_es)

# save the file
write.csv(gsva_es, "gsva_output.csv", quote = F)

# Differential expression analysis of the pathways
# group
group_list <- data.frame(sample = colnames(gsva_es), group = c(rep("a", 100), rep("b", 100)))
head(group_list)

# contrast
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

# Building a difference comparison matrix
contrast.matrix <- makeContrasts(b-a, levels = design)

# Difference analysis, b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

# Save the results of the limma analysis of the pathway to a file
write.csv(x, "gsva_limma.csv", quote = F)

# output the t value
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, "easy_input2_for39bar.csv", quote = F, row.names = F)

# plot
df <- read.csv("easy_input2_for39bar.csv")
head(df)

# grouping according to the score values
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

# sorting according to the score values
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, 
             size = 0.3) + 
  
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),
            size = 3, 
            hjust = "inward" ) +  
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, tumor \n versus non-malignant")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) 

ggsave("gsva.pdf", width = 6, height = 8)

# forestplot
library(tidyverse)
library(survival)
library(forestplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
data <- read.csv("subgroup.csv")
head(data)
dim(data)[1] + 5 

data$mean<-(data$CI_low_0.95 + data$CI_up_0.95)/2
tabletext <- cbind(c("Training.set..n.200.",NA,NA, data$Training.set..n.200., NA),
                   c("Univariable.analysis", NA, NA, data$Univariable.analysis, NA),
                   c("Hazard Ratio\n(95% CI)", NA, NA, ifelse(!is.na(data$Univariable.analysis), paste(round(data$HR, 2), " (", round(data$CI_low_0.95, 2), " to ", round(data$CI_up_0.95, 2),")", sep=""), NA), NA),
                   c("P-value", NA, NA, data$P, NA))
tabletext
forestplot(labeltext=tabletext,
           mean=c(NA,NA,1,data$mean,NA),
           lower=c(NA,NA,1,data$CI_low_0.95,NA),
           upper=c(NA,NA,1,data$CI_up_0.95,NA))

pdf("forestplot.pdf",width=12,height = 7)
forestplot(labeltext=tabletext, 
           mean=c(NA,NA,1,data$mean,NA),#HR
           lower=c(NA,NA,1,data$CI_low_0.95,NA), 
           upper=c(NA,NA,1,data$CI_up_0.95,NA),
           graph.pos=3,
           graphwidth = unit(.4,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="steelblue", lines="black", zero = "black"),
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,
           zero=1,
           lwd.zero=2,
           grid = structure(c(data[1,]$mean), gp = gpar(col = "black", lty=2,lwd=2)),
           lwd.xaxis=2, 
           xlab="Hazard Ratio",
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),
                           "14" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           lineheight = unit(.75,"cm"),
           colgap = unit(0,"cm"),
           mar=unit(rep(1.25, times = 4), "cm"),
           new_page = F
)
dev.off()

# histogram
library(ggplot2)
ggplot(data=dat)+
  geom_bar(
    mapping = aes(x=cluster,fill=neoplasm_histologic_grade),
    position = 'fill'
  )

