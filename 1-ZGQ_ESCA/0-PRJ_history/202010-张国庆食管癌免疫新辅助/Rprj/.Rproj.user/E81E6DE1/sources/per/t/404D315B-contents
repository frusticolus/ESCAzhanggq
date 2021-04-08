library(ggplot2)
require(cowplot)
#require(tidyverse)
require(ggsci)
require(ggpubr)
options(stringsAsFactors = F)

clin = read.table("data/sample_clinical_info.xls", header = T, sep = "\t")
head(clin)

myboxplot <- function(dt, biomarker="GEP"){
  require(ggplot2)
  require(ggpubr)
  comp1 = combn(unique(dt$Response_imaging),2,simplify = FALSE)
  dt$Response_imaging = factor(dt$Response_imaging,levels = c("PR","SD","PD"))
  # print(dt$Response_imaging)
  p1<-ggboxplot(dt, x = "Response_imaging", y = biomarker,
                color = "Response_imaging", palette = "nejm",
                add = "point")+
    stat_compare_means(comparisons = comp1,method="t.test")
  
  comp2 = combn(unique(dt$Response1),2,simplify = FALSE)
  dt$Response1 = factor(dt$Response1,levels = c(1,2,3))
  p2<-ggboxplot(dt, x = "Response1", y = biomarker,
                color = "Response1", palette = "nejm",
                add = "point")+
    stat_compare_means(comparisons = comp2,method="t.test")
  
  comp3 = combn(unique(dt$Response2),2,simplify = FALSE)
  dt$Response2 = factor(dt$Response2,levels = c(1,2,3,4))
  p3<-ggboxplot(dt, x = "Response2", y = biomarker,
                color = "Response2", palette = "nejm",
                add = "point")+
    stat_compare_means(comparisons = comp3,method="t.test")
  
  dt_sub = dt[which(dt$Response_irPRC != ""),]
  dt_sub$irPRC = ifelse(dt_sub$Response_irPRC == "noMPR","noMPR","MPR+pCR")
  comp4 = combn(unique(dt_sub$irPRC),2,simplify = FALSE)
  dt_sub$irPRC = factor(dt_sub$irPRC,levels = c("MPR+pCR","noMPR"))
  p4<-ggboxplot(dt_sub, x = "irPRC", y = biomarker,
                color = "irPRC", palette = "nejm",
                add = "point")+
    stat_compare_means(comparisons = comp4,method="t.test")
  
  fig <- ggarrange(p1,p4,p2,p3,ncol=2,nrow=2,labels=c("A","B","C","D"))
  print(fig)
}

corplot <- function(dt, biomarker="GEP"){
  p1<-ggplot(dt, aes_string(x="RTR", y=biomarker)) +
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson")+
    labs(x="Tumor reduction ratio")
  p2<-ggplot(dt, aes_string(x="Residential_tumor", y=biomarker)) +
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson")+
    labs(x="Residential tumor")
  fig <- ggarrange(p1,p2,ncol=2,nrow=1,labels=c("A","B"))
  print(fig)
}

corplot1 <- function(dt, biomarker="GEP"){
  require(ggplot2)
  require(ggpubr)
  p1<-ggplot(dt, aes_string(x="RTR", y=biomarker)) +
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson")+
    labs(x="Tumor reduction ratio")
  print(p1)
}


###  RNA biomarkers
gep = read.table("data/GEP.tcell_inflamed.GEPscore.xls", header = T, sep = "\t")
gep_clin = merge(gep,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "Sample",by.y = "Alias")
corplot1(gep_clin, biomarker = "GEP")

pdf("output/GEP.pdf",width = 9,height = 9)
myboxplot(gep_clin,biomarker = "GEP")
dev.off()

sig = read.table("data/all_immune_signature_ssgsea.txt", header = T, sep = "\t")
sig = as.data.frame(t(sig))
sig$Sample = rownames(sig)
#colnames(sig) = gsub(" ","_",fixed = TRUE, colnames(sig))

sig_clin = merge(sig,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "Sample",by.y = "Alias")
pdf("output/immune_signature_TRR.pdf",width = 9,height = 9)
for (i in 2:9){
  corplot1(sig_clin, biomarker = colnames(sig_clin)[i])
}
dev.off()
sig_clin = sig_clin[which(sig_clin$Response_imaging != ""),]
pdf("output/immune_signature.pdf",width = 9,height = 9)
for (i in 2:9){
  myboxplot(sig_clin, biomarker = colnames(sig_clin)[i])
}
dev.off()

cell = read.table("data/all_immune_cell_ssgsea.txt", header = T, sep = "\t")
cell = as.data.frame(t(cell))
cell$Sample = rownames(cell)
colnames(cell) = gsub(" ","_",fixed = TRUE, colnames(cell))
cell_clin = merge(cell,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "Sample",by.y = "Alias")
pdf("output/immune_cell_TRR.pdf",width = 9,height = 9)
for (i in 2:ncol(cell)){
  corplot1(cell_clin, biomarker = colnames(cell_clin)[i])
}
dev.off()
cell_clin = cell_clin[which(cell_clin$Response_imaging != ""),]
pdf("output/immune_cell.pdf",width = 9,height = 9)
for (i in 1:ncol(cell)-1){
  myboxplot(cell_clin, biomarker = colnames(cell)[i])
}
dev.off()

cell22 = read.table("data/CIBERSORTx_ESPD_absolute_Results.txt", header = T, sep = "\t")
colnames(cell22) = gsub(".","_",fixed = TRUE, colnames(cell22))
cell22_clin = merge(cell22[,c(1:23,27)],clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                  by.x = "Patient",by.y = "Alias")
pdf("output/immune_cell22_TRR.pdf",width = 9,height = 9)
for (i in 2:24){
  corplot1(cell22_clin, biomarker = colnames(cell22_clin)[i])
}
dev.off()
cell22_clin = cell22_clin[which(cell22_clin$Response_imaging != ""),]
pdf("output/immune_cell22.pdf",width = 9,height = 9)
for (i in 2:24){
  myboxplot(cell22_clin, biomarker = colnames(cell22_clin)[i])
}
dev.off()

ips = read.table("data/IPS.txt", header = T, sep = "\t")
ips_clin = merge(ips,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                    by.x = "SAMPLE",by.y = "Alias")
pdf("output/IPS_TRR.pdf",width = 9,height = 9)
for (i in 2:ncol(ips)){
  corplot1(ips_clin, biomarker = colnames(ips_clin)[i])
}
dev.off()
ips_clin = ips_clin[which(ips_clin$Response_imaging != ""),]
pdf("output/IPS.pdf",width = 9,height = 9)
for (i in 2:ncol(ips)){
  myboxplot(ips_clin, biomarker = colnames(ips_clin)[i])
}
dev.off()


tide = read.table("data/TIDE.result.txt", header = T, sep = "\t")
tide_clin = merge(tide,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "Sample",by.y = "Alias")
pdf("output/TIDE_TRR.pdf",width = 9,height = 9)
for (i in 2:ncol(tide)){
  corplot1(tide_clin, biomarker = colnames(tide_clin)[i])
}
dev.off()
tide_clin = tide_clin[which(tide_clin$Response_imaging != ""),]
pdf("output/TIDE.pdf",width = 9,height = 9)
for (i in 2:ncol(tide)){
  myboxplot(tide_clin, biomarker = colnames(tide_clin)[i])
}
dev.off()

### DNA biomarkers
tmnb = read.table("data/oncoprint_io.txt", header = T, sep = "\t")
tmnb_clin = merge(tmnb,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                  by.x = "Sample",by.y = "Alias")
pdf("output/TMB_TNB_TRR.pdf",width = 12,height = 6)
for (i in 2:3){
  corplot(tmnb_clin, biomarker = colnames(tmnb_clin)[i])
}
dev.off()
tmnb_clin = tmnb_clin[which(tmnb_clin$Response_imaging != ""),]
pdf("output/TMB_TNB.pdf",width = 9,height = 9)
for (i in 2:3){
  myboxplot(tmnb_clin, biomarker = colnames(tmnb_clin)[i])
}
dev.off()
myboxplot(tmnb_clin, biomarker = "MSI")


cin = read.table("data/baseline_CIN.xls", header = T, sep = "\t")
cin_clin = merge(cin,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                  by.x = "Sample",by.y = "Alias")
corplot(cin_clin, biomarker = "CIN")
cin_clin = cin_clin[which(cin_clin$Response_imaging != ""),]
pdf("output/CIN.pdf",width = 9,height = 9)
myboxplot(cin_clin, biomarker = "CIN")
dev.off()

hed = read.table("data/Patient_HED.txt", header = T, sep = "\t")
hed_clin = merge(hed,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "SampleID",by.y = "Alias")
corplot(hed_clin, biomarker = "HED")

hed_clin = hed_clin[which(hed_clin$Response_imaging != ""),]
pdf("output/HED.pdf",width = 9,height = 9)
myboxplot(hed_clin, biomarker = "HED")
dev.off()

neo = read.table("data/baseline.neo.nonsyn.xls", header = T, sep = "\t")
neo_clin = merge(neo,clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "Alias",by.y = "Alias")

corplot(neo_clin, biomarker = "NNI")

neo_clin = neo_clin[which(neo_clin$Residential_tumor != "NA"),]
pdf("output/Neoantigen_frequency.pdf",width = 9,height = 9)
myboxplot(neo_clin, biomarker = "NNI")
myboxplot(neo_clin, biomarker = "NMI")
dev.off()

ith = read.table("data/ITH.txt", header = T, sep = "\t")
ith_clin = merge(ith[,c("Sample","ITH")],clin[,c("Alias","RTR","Residential_tumor","Response_imaging","Response_irPRC","Response1","Response2")],
                 by.x = "Sample",by.y = "Alias")
corplot(ith_clin, biomarker = "ITH")

ith_clin = ith_clin[which(ith_clin$Response_imaging != ""),]
myboxplot(ith_clin, biomarker = "ITH")


